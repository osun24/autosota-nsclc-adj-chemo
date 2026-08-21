"""Fit an RSF to pinned MSigDB Reactome genes plus clinical covariates.

Only training and validation CSVs are accepted. The repository's sealed test
set is deliberately outside this program's data contract.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import platform
from pathlib import Path
import pickle
import time

import numpy as np
import pandas as pd
import sklearn
import sksurv
from lifelines import KaplanMeierFitter
from lifelines.utils import restricted_mean_survival_time
from sklearn.impute import SimpleImputer
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv

from .config import (
    CLINICAL_COVARIATES,
    DEFAULT_DATA_DIR,
    DEFAULT_RSF_PARAMS,
    DEFAULT_RUNS_DIR,
    DEFAULT_TRAIN_CSV,
    DEFAULT_VALID_CSV,
    MSIGDB_FILENAME,
    MSIGDB_VERSION,
    OUTCOME_COLUMNS,
)
from .msigdb import parse_gmt, retrieve_collection, sha256_file


def _assert_training_or_validation_path(path: str | Path) -> None:
    if "test" in str(path).lower():
        raise ValueError(f"Sealed test paths are not accepted: {path}")


def _headers(path: Path) -> list[str]:
    return list(pd.read_csv(path, nrows=0).columns)


def _prepare_frame(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    mapped = out["Adjuvant Chemo"].map({"OBS": 0, "ACT": 1})
    out["Adjuvant Chemo"] = mapped.where(mapped.notna(), out["Adjuvant Chemo"])
    for column in out.columns:
        out[column] = pd.to_numeric(out[column], errors="coerce")
    if out[OUTCOME_COLUMNS].isna().any().any():
        raise ValueError("Survival outcome columns contain missing or non-numeric values")
    if not set(out["OS_STATUS"].astype(int).unique()).issubset({0, 1}):
        raise ValueError("OS_STATUS must be binary")
    if (out["OS_MONTHS"] <= 0).any():
        raise ValueError("OS_MONTHS must be positive")
    return out


def _outcome(df: pd.DataFrame) -> np.ndarray:
    return Surv.from_arrays(
        event=df["OS_STATUS"].astype(bool).to_numpy(),
        time=df["OS_MONTHS"].astype(float).to_numpy(),
    )


def select_available_reactome_genes(
    gmt_path: str | Path, train_columns: list[str], valid_columns: list[str]
) -> tuple[list[str], int]:
    gene_sets = parse_gmt(gmt_path)
    all_genes = {gene for members in gene_sets.values() for gene in members}
    shared = set(train_columns).intersection(valid_columns)
    return sorted(all_genes.intersection(shared)), len(all_genes)


def _dataset_summary(df: pd.DataFrame) -> dict:
    return {
        "n": int(len(df)),
        "events": int(df["OS_STATUS"].sum()),
        "censored": int(len(df) - df["OS_STATUS"].sum()),
        "treatment_counts": {
            str(int(key)): int(value)
            for key, value in df["Adjuvant Chemo"].value_counts().sort_index().items()
        },
    }


def _bootstrap_cindex_interval(
    event: np.ndarray,
    time_values: np.ndarray,
    risk: np.ndarray,
    draws: int = 1000,
    seed: int = 20260816,
) -> dict:
    """Patient-level percentile interval for an already-fitted prediction vector."""
    rng = np.random.default_rng(seed)
    estimates: list[float] = []
    n = len(risk)
    for _ in range(draws):
        indices = rng.integers(0, n, size=n)
        try:
            estimate = concordance_index_censored(
                event[indices], time_values[indices], risk[indices]
            )[0]
        except ValueError:
            continue
        if np.isfinite(estimate):
            estimates.append(float(estimate))
    if not estimates:
        raise RuntimeError("No valid validation C-index bootstrap draws")
    lower, upper = np.percentile(estimates, [2.5, 97.5])
    return {
        "method": "patient-level percentile bootstrap",
        "requested_draws": int(draws),
        "valid_draws": len(estimates),
        "seed": int(seed),
        "lower_95": float(lower),
        "upper_95": float(upper),
    }


def _alignment_rmst_difference(
    time_values: np.ndarray,
    event: np.ndarray,
    aligned: np.ndarray,
    tau: float = 60.0,
) -> dict:
    """Compute RMST(aligned) minus RMST(not aligned) on validation rows."""
    aligned = np.asarray(aligned, dtype=bool)
    if int(aligned.sum()) == 0 or int((~aligned).sum()) == 0:
        raise ValueError("Delta RMST requires both aligned and non-aligned patients")
    km_aligned = KaplanMeierFitter().fit(
        time_values[aligned], event_observed=event[aligned]
    )
    km_not_aligned = KaplanMeierFitter().fit(
        time_values[~aligned], event_observed=event[~aligned]
    )
    aligned_rmst = float(restricted_mean_survival_time(km_aligned, t=tau))
    not_aligned_rmst = float(restricted_mean_survival_time(km_not_aligned, t=tau))
    return {
        "tau_months": float(tau),
        "aligned_n": int(aligned.sum()),
        "not_aligned_n": int((~aligned).sum()),
        "aligned_rmst_months": aligned_rmst,
        "not_aligned_rmst_months": not_aligned_rmst,
        "delta_rmst_months": aligned_rmst - not_aligned_rmst,
        "definition": (
            "RMST(observed treatment matches RSF recommendation) - "
            "RMST(observed treatment does not match RSF recommendation)"
        ),
    }


def _bootstrap_delta_rmst_interval(
    time_values: np.ndarray,
    event: np.ndarray,
    aligned: np.ndarray,
    tau: float = 60.0,
    draws: int = 1000,
    seed: int = 20260816,
) -> dict:
    """Validation-patient bootstrap interval for the fixed recommendations."""
    rng = np.random.default_rng(seed)
    estimates: list[float] = []
    n = len(aligned)
    for _ in range(draws):
        indices = rng.integers(0, n, size=n)
        sampled_alignment = aligned[indices]
        if int(sampled_alignment.sum()) == 0 or int((~sampled_alignment).sum()) == 0:
            continue
        estimate = _alignment_rmst_difference(
            time_values[indices], event[indices], sampled_alignment, tau=tau
        )["delta_rmst_months"]
        if np.isfinite(estimate):
            estimates.append(float(estimate))
    if not estimates:
        raise RuntimeError("No valid validation delta-RMST bootstrap draws")
    lower, upper = np.percentile(estimates, [2.5, 97.5])
    return {
        "method": "patient-level percentile bootstrap",
        "requested_draws": int(draws),
        "valid_draws": len(estimates),
        "seed": int(seed),
        "lower_95": float(lower),
        "upper_95": float(upper),
    }


def run(
    train_csv: str | Path = DEFAULT_TRAIN_CSV,
    valid_csv: str | Path = DEFAULT_VALID_CSV,
    data_dir: str | Path = DEFAULT_DATA_DIR,
    runs_dir: str | Path = DEFAULT_RUNS_DIR,
    n_estimators: int = int(DEFAULT_RSF_PARAMS["n_estimators"]),
    random_state: int = int(DEFAULT_RSF_PARAMS["random_state"]),
    save_model: bool = False,
) -> dict:
    train_csv, valid_csv = Path(train_csv), Path(valid_csv)
    _assert_training_or_validation_path(train_csv)
    _assert_training_or_validation_path(valid_csv)
    collection = retrieve_collection(data_dir)
    gmt_path = Path(data_dir) / MSIGDB_FILENAME

    train_header, valid_header = _headers(train_csv), _headers(valid_csv)
    missing_clinical = [
        col for col in CLINICAL_COVARIATES
        if col not in train_header or col not in valid_header
    ]
    if missing_clinical:
        raise ValueError(f"Missing clinical covariates: {missing_clinical}")
    genes, collection_gene_count = select_available_reactome_genes(
        gmt_path, train_header, valid_header
    )
    if not genes:
        raise ValueError("No pinned Reactome gene symbols match the expression columns")

    requested = OUTCOME_COLUMNS + CLINICAL_COVARIATES + genes
    train_df = _prepare_frame(pd.read_csv(train_csv, usecols=requested))
    valid_df = _prepare_frame(pd.read_csv(valid_csv, usecols=requested))

    # All preprocessing decisions are fitted on training rows only. Clinical
    # variables remain prespecified; only training-constant gene columns drop.
    gene_values = train_df[genes].to_numpy(dtype=np.float32)
    with np.errstate(all="ignore"):
        gene_ranges = np.nanmax(gene_values, axis=0) - np.nanmin(gene_values, axis=0)
    kept_genes = [gene for gene, span in zip(genes, gene_ranges) if np.isfinite(span) and span > 0]
    dropped_constant_genes = sorted(set(genes) - set(kept_genes))
    feature_names = CLINICAL_COVARIATES + kept_genes

    imputer = SimpleImputer(strategy="median")
    x_train = imputer.fit_transform(
        train_df[feature_names].to_numpy(dtype=np.float32)
    ).astype(np.float32, copy=False)
    x_valid = imputer.transform(
        valid_df[feature_names].to_numpy(dtype=np.float32)
    ).astype(np.float32, copy=False)

    params = dict(DEFAULT_RSF_PARAMS)
    params["n_estimators"] = int(n_estimators)
    params["random_state"] = int(random_state)
    model = RandomSurvivalForest(**params)
    started = time.perf_counter()
    model.fit(x_train, _outcome(train_df))
    fit_seconds = time.perf_counter() - started
    valid_risk = model.predict(x_valid)
    valid_cindex = float(
        concordance_index_censored(
            valid_df["OS_STATUS"].astype(bool).to_numpy(),
            valid_df["OS_MONTHS"].astype(float).to_numpy(),
            valid_risk,
        )[0]
    )
    valid_event = valid_df["OS_STATUS"].astype(bool).to_numpy()
    valid_time = valid_df["OS_MONTHS"].astype(float).to_numpy()
    valid_cindex_interval = _bootstrap_cindex_interval(
        valid_event, valid_time, valid_risk
    )

    # For each validation patient, recommend the treatment producing lower
    # predicted risk when only the treatment indicator is counterfactually set.
    valid_act = valid_df.copy()
    valid_act["Adjuvant Chemo"] = 1
    valid_obs = valid_df.copy()
    valid_obs["Adjuvant Chemo"] = 0
    x_valid_act = imputer.transform(
        valid_act[feature_names].to_numpy(dtype=np.float32)
    ).astype(np.float32, copy=False)
    x_valid_obs = imputer.transform(
        valid_obs[feature_names].to_numpy(dtype=np.float32)
    ).astype(np.float32, copy=False)
    valid_risk_act = model.predict(x_valid_act)
    valid_risk_obs = model.predict(x_valid_obs)
    recommendation = (valid_risk_act < valid_risk_obs).astype(int)
    observed_treatment = valid_df["Adjuvant Chemo"].astype(int).to_numpy()
    aligned = observed_treatment == recommendation
    delta_rmst = _alignment_rmst_difference(
        valid_time, valid_event, aligned, tau=60.0
    )
    delta_rmst["confidence_interval"] = _bootstrap_delta_rmst_interval(
        valid_time, valid_event, aligned, tau=60.0
    )

    run_id = datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%SZ")
    run_dir = Path(runs_dir) / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    (run_dir / "feature_names.txt").write_text(
        "\n".join(feature_names) + "\n", encoding="utf-8"
    )
    (run_dir / "reactome_genes_used.txt").write_text(
        "\n".join(kept_genes) + "\n", encoding="utf-8"
    )
    medians = pd.DataFrame(
        {"feature": feature_names, "training_median": imputer.statistics_}
    )
    medians.to_csv(run_dir / "training_imputation_medians.csv", index=False)
    pd.DataFrame({
        "validation_row": np.arange(len(valid_df)),
        "OS_STATUS": valid_event.astype(int),
        "OS_MONTHS": valid_time,
        "observed_treatment": observed_treatment,
        "observed_treatment_risk": valid_risk,
        "counterfactual_risk_act": valid_risk_act,
        "counterfactual_risk_obs": valid_risk_obs,
        "recommended_treatment": recommendation,
        "treatment_aligned": aligned.astype(int),
    }).to_csv(run_dir / "validation_predictions.csv", index=False)
    if save_model:
        with (run_dir / "rsf_model.pkl").open("wb") as handle:
            pickle.dump({"model": model, "imputer": imputer}, handle)

    result = {
        "schema_version": 1,
        "run_id": run_id,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "analysis": "Random Survival Forest: MSigDB Reactome genes + clinical covariates",
        "collection": collection,
        "inputs": {
            "train_csv": str(train_csv.resolve()),
            "train_sha256": sha256_file(train_csv),
            "validation_csv": str(valid_csv.resolve()),
            "validation_sha256": sha256_file(valid_csv),
            "test_data_used": False,
        },
        "cohorts": {
            "training": _dataset_summary(train_df),
            "validation": _dataset_summary(valid_df),
        },
        "features": {
            "clinical_covariates": CLINICAL_COVARIATES,
            "n_clinical_covariates": len(CLINICAL_COVARIATES),
            "n_unique_genes_in_collection": collection_gene_count,
            "n_genes_shared_with_expression_data": len(genes),
            "n_training_constant_genes_dropped": len(dropped_constant_genes),
            "training_constant_genes_dropped": dropped_constant_genes,
            "n_reactome_genes_used": len(kept_genes),
            "n_total_features": len(feature_names),
            "selection_rule": (
                "Union of pinned GMT symbols present in both train and validation; "
                "training-constant genes removed without using outcomes"
            ),
            "imputation": "training-set median; applied unchanged to validation",
        },
        "model": {
            "class": "sksurv.ensemble.RandomSurvivalForest",
            "parameters": params,
            "fit_seconds": fit_seconds,
            "saved": bool(save_model),
        },
        "validation": {
            "harrell_cindex": valid_cindex,
            "metric": "sksurv.metrics.concordance_index_censored",
            "confidence_interval": valid_cindex_interval,
            "delta_rmst": delta_rmst,
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
            "scikit_survival": sksurv.__version__,
        },
    }
    (run_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    result["run_dir"] = str(run_dir.resolve())
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-csv", type=Path, default=DEFAULT_TRAIN_CSV)
    parser.add_argument("--valid-csv", type=Path, default=DEFAULT_VALID_CSV)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--runs-dir", type=Path, default=DEFAULT_RUNS_DIR)
    parser.add_argument("--n-estimators", type=int, default=DEFAULT_RSF_PARAMS["n_estimators"])
    parser.add_argument("--random-state", type=int, default=DEFAULT_RSF_PARAMS["random_state"])
    parser.add_argument("--save-model", action="store_true")
    args = parser.parse_args()
    result = run(
        train_csv=args.train_csv,
        valid_csv=args.valid_csv,
        data_dir=args.data_dir,
        runs_dir=args.runs_dir,
        n_estimators=args.n_estimators,
        random_state=args.random_state,
        save_model=args.save_model,
    )
    print(json.dumps({
        "run_dir": result["run_dir"],
        "validation": result["validation"],
        "n_total_features": result["features"]["n_total_features"],
        "fit_seconds": result["model"]["fit_seconds"],
    }, indent=2))


if __name__ == "__main__":
    main()
