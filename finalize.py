"""Human-only sealed test evaluation for arena artifacts.

Run this manually after `affyfRMATest.csv` has been added locally. The
autonomous loop must never import or execute this file.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib")

import numpy as np
import pandas as pd
import xgboost as xgb
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test
from lifelines.utils import restricted_mean_survival_time
from sklearn.model_selection import train_test_split

from xgb_arena import prepare as xgb_prepare
from xgb_arena import train as xgb_train

REPO_ROOT = Path(__file__).resolve().parent
TEST_CSV = REPO_ROOT / "affyfRMATest.csv"
XGB_RUNS_DIR = REPO_ROOT / "xgb_arena" / "runs"


def _latest_run_dir(runs_dir: Path) -> Path:
    candidates = [p for p in runs_dir.iterdir() if p.is_dir() and (p / "metadata.json").exists()]
    if not candidates:
        raise FileNotFoundError(f"No run directories with metadata.json found in {runs_dir}")
    return sorted(candidates, key=lambda p: p.stat().st_mtime, reverse=True)[0]


def _read_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def _needed_gene_names(feature_names: list[str]) -> list[str]:
    genes = []
    for feat in feature_names:
        base = feat.split("*ACT", 1)[0] if "*ACT" in feat else feat
        if base not in xgb_prepare.CLINICAL_VARS and base not in {"OS_STATUS", "OS_MONTHS"}:
            genes.append(base)
    return sorted(set(genes))


def _load_test_for_xgb(feature_names: list[str]) -> pd.DataFrame:
    if not TEST_CSV.exists():
        raise FileNotFoundError(
            f"{TEST_CSV} does not exist. Add the sealed test CSV locally, then rerun finalize.py manually."
        )
    genes = _needed_gene_names(feature_names)
    test_raw = pd.read_csv(TEST_CSV)
    test_df = xgb_prepare.preprocess_split(test_raw, xgb_prepare.CLINICAL_VARS, genes)
    return test_df.sort_values(by=["OS_MONTHS", "OS_STATUS"], ascending=[False, False]).reset_index(drop=True)


def _clin_cols_from_metadata(metadata: dict) -> list[str]:
    if "clin_cols" in metadata:
        return list(metadata["clin_cols"])
    if "clin_pretx" in metadata:
        return ["Adjuvant Chemo", *list(metadata["clin_pretx"])]
    return list(xgb_prepare.CLINICAL_VARS)


def _xgb_gene_stability(run_dir: Path, genes_main: list[str]) -> dict:
    current = set(genes_main)
    overlaps = []
    for other_meta in XGB_RUNS_DIR.glob("*/metadata.json"):
        other_dir = other_meta.parent
        if other_dir == run_dir:
            continue
        try:
            other = json.loads(other_meta.read_text())
            other_genes = set(other.get("genes_main", []))
        except Exception:
            continue
        if not other_genes:
            continue
        union = current | other_genes
        overlaps.append(len(current & other_genes) / len(union) if union else np.nan)
    return {
        "gene_set_jaccard_vs_other_runs_mean": float(np.nanmean(overlaps)) if overlaps else np.nan,
        "gene_set_jaccard_n_comparisons": int(len(overlaps)),
    }


def _best_ntree(metadata: dict, key: str, fallback: int = 0) -> int:
    result = metadata.get("result", {})
    return int(metadata.get(key, result.get(key, fallback)))


def _format_pvalue(pvalue: float) -> str:
    if np.isnan(pvalue):
        return "NA"
    if pvalue < 0.0001:
        return f"{pvalue:.2e}"
    return f"{pvalue:.4f}"


def _km_alignment_summary(
    df: pd.DataFrame,
    alignment: np.ndarray,
    figure_path: Path,
    cindex: float,
    tau: float = 60,
) -> dict:
    mask_aligned = pd.Series(alignment, index=df.index).astype(bool)
    mask_not_aligned = ~mask_aligned
    if int(mask_aligned.sum()) == 0 or int(mask_not_aligned.sum()) == 0:
        return {
            "test_rmst_diff": 0.0,
            "km_logrank_pvalue": np.nan,
            "km_figure": "",
        }

    kmf_aligned = KaplanMeierFitter()
    kmf_not_aligned = KaplanMeierFitter()
    kmf_aligned.fit(
        durations=df.loc[mask_aligned, "OS_MONTHS"],
        event_observed=df.loc[mask_aligned, "OS_STATUS"],
        label="Aligned with XGBoost recommendation",
    )
    kmf_not_aligned.fit(
        durations=df.loc[mask_not_aligned, "OS_MONTHS"],
        event_observed=df.loc[mask_not_aligned, "OS_STATUS"],
        label="Not aligned with XGBoost recommendation",
    )

    logrank = logrank_test(
        df.loc[mask_aligned, "OS_MONTHS"],
        df.loc[mask_not_aligned, "OS_MONTHS"],
        event_observed_A=df.loc[mask_aligned, "OS_STATUS"],
        event_observed_B=df.loc[mask_not_aligned, "OS_STATUS"],
    )
    pvalue = float(logrank.p_value)
    rmst_aligned = float(restricted_mean_survival_time(kmf_aligned, t=tau))
    rmst_not_aligned = float(restricted_mean_survival_time(kmf_not_aligned, t=tau))
    rmst_diff = rmst_aligned - rmst_not_aligned

    figure_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(10, 6))
    ax = kmf_aligned.plot(ci_show=True)
    kmf_not_aligned.plot(ax=ax, ci_show=True)
    ax.set_title("Kaplan-Meier Survival Curves by Treatment Alignment")
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival Probability")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.text(0.08, 0.20, f"Log-rank p-value: {_format_pvalue(pvalue)}", transform=ax.transAxes)
    ax.text(0.08, 0.14, f"C-index: {float(cindex):.4f}", transform=ax.transAxes)
    ax.text(0.08, 0.08, f"5-year RMST difference: {rmst_diff:.2f} months", transform=ax.transAxes)
    add_at_risk_counts(kmf_aligned, kmf_not_aligned, ax=ax)
    plt.tight_layout()
    plt.savefig(figure_path, dpi=200, bbox_inches="tight")
    plt.close()

    return {
        "test_rmst_diff": float(rmst_diff),
        "km_logrank_pvalue": pvalue,
        "km_figure": str(figure_path.relative_to(REPO_ROOT) if figure_path.is_relative_to(REPO_ROOT) else figure_path),
    }


def _fit_tlearner_from_metadata(
    metadata: dict,
    train_df: pd.DataFrame,
    valid_df: pd.DataFrame,
    feature_names: list[str],
) -> tuple[xgb.Booster, xgb.Booster, int, int]:
    clin_pretx = list(
        metadata.get("clin_pretx") or [c for c in _clin_cols_from_metadata(metadata) if c != "Adjuvant Chemo"]
    )
    params0 = dict(metadata["xgb_params"])
    params1 = dict(metadata["xgb_params"])
    params0["seed"] = int(params0.get("seed", 7))
    params1["seed"] = int(params1.get("seed", 7)) + 1
    num_boost_round = int(metadata["num_boost_round"])
    early_stopping_rounds = int(metadata["early_stopping_rounds"])

    trainval_df = pd.concat([train_df, valid_df], axis=0, ignore_index=True)
    trainval_df = trainval_df.sort_values(
        by=["OS_MONTHS", "OS_STATUS"], ascending=[False, False]
    ).reset_index(drop=True)
    weights, _, _ = xgb_prepare.compute_iptw(trainval_df, covariate_cols=clin_pretx, act_col="Adjuvant Chemo")

    boosters: list[xgb.Booster] = []
    best_ntrees: list[int] = []
    for arm, params in [(0, params0), (1, params1)]:
        arm_mask = trainval_df["Adjuvant Chemo"].astype(int).to_numpy() == arm
        arm_df = trainval_df.loc[arm_mask].reset_index(drop=True)
        arm_weights = weights[arm_mask]
        idx = np.arange(len(arm_df))
        events = arm_df["OS_STATUS"].astype(int).to_numpy()
        stratify = events if np.bincount(events, minlength=2).min() >= 2 else None
        tr_idx, va_idx = train_test_split(
            idx,
            test_size=0.25,
            random_state=42 + arm,
            stratify=stratify,
        )
        X_arm = arm_df[feature_names].to_numpy(dtype=np.float32)
        dtr = xgb_train.make_dmatrix(
            X_arm[tr_idx],
            arm_df.loc[tr_idx, "OS_MONTHS"].values,
            arm_df.loc[tr_idx, "OS_STATUS"].values,
            arm_weights[tr_idx],
            feature_names,
        )
        dva = xgb_train.make_dmatrix(
            X_arm[va_idx],
            arm_df.loc[va_idx, "OS_MONTHS"].values,
            arm_df.loc[va_idx, "OS_STATUS"].values,
            None,
            feature_names,
        )
        booster, _ = xgb_train.train_xgb_cox(dtr, dva, params, num_boost_round, early_stopping_rounds)
        best_ntree = booster.best_iteration + 1 if booster.best_iteration is not None else num_boost_round
        boosters.append(xgb_prepare.slice_booster_to_best_iteration(booster, best_ntree))
        best_ntrees.append(int(best_ntree))

    return boosters[0], boosters[1], best_ntrees[0], best_ntrees[1]


def evaluate_xgb_tlearner(run_dir: Path, refit_train_valid: bool, km_figure_path: Path) -> dict:
    metadata = json.loads((run_dir / "metadata.json").read_text())
    feature_names = (
        metadata.get("feature_names") or metadata.get("feat_names_t") or _read_lines(run_dir / "feature_names.txt")
    )
    genes_main = metadata.get("genes_main") or _read_lines(run_dir / "genes_main.txt")
    genes_inter = metadata.get("genes_inter") or []
    test_df = _load_test_for_xgb(feature_names)

    if refit_train_valid:
        train_df, valid_df = xgb_prepare.load_train_valid()
        booster0, booster1, best0, best1 = _fit_tlearner_from_metadata(metadata, train_df, valid_df, feature_names)
        model_source = "refit_train_valid"
    else:
        booster0 = xgb.Booster()
        booster1 = xgb.Booster()
        booster0.load_model(run_dir / "xgb_model_arm0.json")
        booster1.load_model(run_dir / "xgb_model_arm1.json")
        best0 = _best_ntree(metadata, "best_ntree_arm0")
        best1 = _best_ntree(metadata, "best_ntree_arm1")
        model_source = "saved_validation_model"

    X_test = test_df[feature_names].to_numpy(dtype=np.float32)
    dtest = xgb.DMatrix(X_test, feature_names=feature_names)
    risk0 = booster0.predict(dtest, iteration_range=(0, int(best0)), output_margin=True)
    risk1 = booster1.predict(dtest, iteration_range=(0, int(best1)), output_margin=True)
    test_ci = xgb_prepare.cindex(
        risk0,
        test_df["OS_MONTHS"].to_numpy(dtype=float),
        test_df["OS_STATUS"].to_numpy(dtype=int),
    )

    model_rec = (risk1 < risk0).astype(int)
    actual = test_df["Adjuvant Chemo"].astype(int).to_numpy()
    alignment = actual == model_rec
    km = _km_alignment_summary(test_df, alignment, km_figure_path, cindex=test_ci)

    stability = _xgb_gene_stability(run_dir, genes_main)
    return {
        "model": "xgb_tlearner",
        "run_dir": str(run_dir.relative_to(REPO_ROOT)),
        "model_source": model_source,
        "test_ci": float(test_ci),
        "test_rmst_diff": float(km["test_rmst_diff"]),
        "km_logrank_pvalue": float(km["km_logrank_pvalue"]),
        "km_figure": km["km_figure"],
        "c_for_benefit": np.nan,
        "c_for_benefit_note": "Not computed in Phase 1 scaffold; add a versioned definition before using for claims.",
        "n_features": int(len(feature_names)),
        "k_main": int(len(genes_main)),
        "k_int": int(len(genes_inter)),
        **stability,
    }


def evaluate_xgb(run_dir: Path, refit_train_valid: bool, km_figure_path: Path) -> dict:
    metadata = json.loads((run_dir / "metadata.json").read_text())
    if metadata.get("arena") == "xgb_tlearner":
        return evaluate_xgb_tlearner(run_dir, refit_train_valid, km_figure_path)
    if metadata.get("arena") not in {None, "xgb"}:
        raise NotImplementedError(
            f"finalize.py does not yet support arena={metadata.get('arena')!r} in {run_dir}"
        )

    feature_names = metadata.get("feature_names") or _read_lines(run_dir / "feature_names.txt")
    genes_main = metadata.get("genes_main") or _read_lines(run_dir / "genes_main.txt")
    genes_inter = metadata.get("genes_inter") or _read_lines(run_dir / "genes_inter.txt")
    clin_cols = _clin_cols_from_metadata(metadata)

    test_df = _load_test_for_xgb(feature_names)

    if refit_train_valid:
        train_df, valid_df = xgb_prepare.load_train_valid()
        booster, best_ntree, feature_names, genes_main, genes_inter, clin_cols = xgb_train.fit_model_from_metadata(
            metadata,
            train_df,
            valid_df,
            refit_train_valid=True,
        )
        model_source = "refit_train_valid"
    else:
        booster = xgb.Booster()
        booster.load_model(run_dir / "xgb_model.json")
        best_ntree = _best_ntree(metadata, "best_ntree")
        model_source = "saved_validation_model"

    X_test = xgb_prepare.build_matrix_from_feature_names(test_df, feature_names)
    pred = xgb_prepare.predict_xgb_risk(booster, X_test, feature_names, best_ntree)
    test_ci = xgb_prepare.cindex(
        pred,
        test_df["OS_MONTHS"].to_numpy(dtype=float),
        test_df["OS_STATUS"].to_numpy(dtype=int),
    )
    del genes_main, genes_inter, clin_cols
    test_counterfactual_treated = test_df.copy()
    test_counterfactual_treated["Adjuvant Chemo"] = 1
    test_counterfactual_untreated = test_df.copy()
    test_counterfactual_untreated["Adjuvant Chemo"] = 0
    X_treated = xgb_prepare.build_matrix_from_feature_names(test_counterfactual_treated, feature_names)
    X_untreated = xgb_prepare.build_matrix_from_feature_names(test_counterfactual_untreated, feature_names)
    risk_treated = xgb_prepare.predict_xgb_risk(booster, X_treated, feature_names, best_ntree)
    risk_untreated = xgb_prepare.predict_xgb_risk(booster, X_untreated, feature_names, best_ntree)
    model_rec = np.where(risk_treated < risk_untreated, 1, 0)
    alignment = test_df["Adjuvant Chemo"].to_numpy(int) == model_rec
    km = _km_alignment_summary(test_df, alignment, km_figure_path, cindex=test_ci)

    genes_main = metadata.get("genes_main") or _read_lines(run_dir / "genes_main.txt")
    genes_inter = metadata.get("genes_inter") or _read_lines(run_dir / "genes_inter.txt")
    stability = _xgb_gene_stability(run_dir, genes_main)
    return {
        "model": "xgb",
        "run_dir": str(run_dir.relative_to(REPO_ROOT)),
        "model_source": model_source,
        "test_ci": float(test_ci),
        "test_rmst_diff": float(km["test_rmst_diff"]),
        "km_logrank_pvalue": float(km["km_logrank_pvalue"]),
        "km_figure": km["km_figure"],
        "c_for_benefit": np.nan,
        "c_for_benefit_note": "Not computed in Phase 1 scaffold; add a versioned definition before using for claims.",
        "n_features": int(len(feature_names)),
        "k_main": int(len(genes_main)),
        "k_int": int(len(genes_inter)),
        **stability,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Human-only sealed test evaluation.")
    parser.add_argument("--xgb-run-dir", type=Path, default=None, help="Run directory under xgb_arena/runs.")
    parser.add_argument(
        "--no-refit-train-valid",
        action="store_true",
        help="Evaluate the saved validation model instead of refitting on Train+Val before test.",
    )
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "final_comparison.csv")
    parser.add_argument(
        "--km-fig",
        type=Path,
        default=None,
        help="Path for the Kaplan-Meier treatment-alignment figure.",
    )
    args = parser.parse_args()

    xgb_run_dir = args.xgb_run_dir or _latest_run_dir(XGB_RUNS_DIR)
    km_figure_path = args.km_fig or args.out.with_name(f"{args.out.stem}_km.png")
    row = evaluate_xgb(
        xgb_run_dir.resolve(),
        refit_train_valid=not args.no_refit_train_valid,
        km_figure_path=km_figure_path.resolve(),
    )
    table = pd.DataFrame([row])
    print(table.to_string(index=False))
    table.to_csv(args.out, index=False)
    print(f"\nWrote comparison table to {args.out}")
    if row.get("km_figure"):
        print(f"Wrote Kaplan-Meier figure to {row['km_figure']}")


if __name__ == "__main__":
    main()
