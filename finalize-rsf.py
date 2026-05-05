"""Human-only sealed test evaluation for RSF arena artifacts.

Run this manually after affyfRMATest.csv has been added locally. The
autonomous loop must never import or execute this file.

Supports both single-model runs and 10-seed ensemble runs (iter_011+).
Ensemble runs are detected automatically from metadata["result"]["seed_panel_cis"].
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test
from lifelines.utils import restricted_mean_survival_time

from rsf_arena import prepare as rsf_prepare
from rsf_arena import train as rsf_train

REPO_ROOT = Path(__file__).resolve().parent
TEST_CSV = REPO_ROOT / "affyfRMATest.csv"
RSF_RUNS_DIR = REPO_ROOT / "rsf_arena" / "runs"

# Must match train.py FINAL_SEEDS used in iter_011+
FINAL_SEEDS = [7, 13, 21, 37, 53, 71, 89, 97, 113, 127]


def latest_run_dir(runs_dir: Path = RSF_RUNS_DIR) -> Path:
    candidates = [p for p in runs_dir.iterdir() if p.is_dir() and (p / "metadata.json").exists()]
    if not candidates:
        raise FileNotFoundError(f"No RSF run directories with metadata.json found in {runs_dir}")
    return sorted(candidates, key=lambda p: p.stat().st_mtime, reverse=True)[0]


def _needed_gene_names(feature_names: list[str]) -> list[str]:
    genes = []
    for feat in feature_names:
        base = feat.split("*ACT", 1)[0] if "*ACT" in feat else feat
        if base not in rsf_prepare.CLINICAL_VARS and base not in {"OS_STATUS", "OS_MONTHS"}:
            genes.append(base)
    return sorted(set(genes))


def _load_test(feature_names: list[str]) -> pd.DataFrame:
    if not TEST_CSV.exists():
        raise FileNotFoundError(
            f"{TEST_CSV} does not exist. Add the sealed test CSV locally, then rerun manually."
        )
    test_raw = pd.read_csv(TEST_CSV)
    test_df = rsf_prepare.preprocess_split(test_raw, rsf_prepare.CLINICAL_VARS, _needed_gene_names(feature_names))
    return test_df.sort_values(["OS_MONTHS", "OS_STATUS"], ascending=[False, False]).reset_index(drop=True)


def _gene_stability(run_dir: Path, genes_main: list[str]) -> dict:
    current = set(genes_main)
    overlaps = []
    for meta_path in RSF_RUNS_DIR.glob("*/metadata.json"):
        if meta_path.parent == run_dir:
            continue
        try:
            other = json.loads(meta_path.read_text())
            other_genes = set(other.get("genes_main", []))
        except Exception:
            continue
        if other_genes:
            overlaps.append(len(current & other_genes) / len(current | other_genes))
    return {
        "gene_set_jaccard_vs_other_runs_mean": float(np.nanmean(overlaps)) if overlaps else np.nan,
        "gene_set_jaccard_n_comparisons": int(len(overlaps)),
    }


def _format_pvalue(pvalue: float) -> str:
    if np.isnan(pvalue):
        return "NA"
    if pvalue < 0.0001:
        return f"{pvalue:.2e}"
    return f"{pvalue:.4f}"


def _km_alignment_plot(
    df: pd.DataFrame,
    alignment: np.ndarray,
    figure_path: Path,
    cindex: float,
    rmst_diff: float,
    logrank_p: float,
) -> None:
    mask_a = pd.Series(alignment, index=df.index).astype(bool)
    mask_n = ~mask_a
    kmf_a = KaplanMeierFitter().fit(
        df.loc[mask_a, "OS_MONTHS"],
        event_observed=df.loc[mask_a, "OS_STATUS"],
        label="Aligned with RSF recommendation",
    )
    kmf_n = KaplanMeierFitter().fit(
        df.loc[mask_n, "OS_MONTHS"],
        event_observed=df.loc[mask_n, "OS_STATUS"],
        label="Not aligned with RSF recommendation",
    )
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(10, 6))
    ax = kmf_a.plot(ci_show=True)
    kmf_n.plot(ax=ax, ci_show=True)
    ax.set_title("Kaplan-Meier Survival Curves by Treatment Alignment (RSF Ensemble)")
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival Probability")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.text(0.08, 0.20, f"Log-rank p-value: {_format_pvalue(logrank_p)}", transform=ax.transAxes)
    ax.text(0.08, 0.14, f"C-index: {cindex:.4f}", transform=ax.transAxes)
    ax.text(0.08, 0.08, f"5-year RMST difference: {rmst_diff:.2f} months", transform=ax.transAxes)
    add_at_risk_counts(kmf_a, kmf_n, ax=ax)
    plt.tight_layout()
    plt.savefig(figure_path, dpi=200, bbox_inches="tight")
    plt.close()


def _ensemble_predict(
    metadata: dict,
    train_df: pd.DataFrame,
    valid_df: pd.DataFrame,
    eval_df: pd.DataFrame,
    seeds: list[int] = FINAL_SEEDS,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Fit one RSF per seed on train+valid; return averaged risk scores on eval_df."""
    eval_tr = eval_df.copy()
    eval_tr["Adjuvant Chemo"] = 1
    eval_co = eval_df.copy()
    eval_co["Adjuvant Chemo"] = 0
    risks, risks_tr, risks_co = [], [], []
    feat_names_out: list[str] = []
    for i, seed in enumerate(seeds, 1):
        print(f"  [Ensemble] Fitting seed {seed} ({i}/{len(seeds)})...")
        meta = dict(metadata)
        meta["rsf_params"] = dict(metadata["rsf_params"])
        meta["rsf_params"]["random_state"] = int(seed)
        m, fn, _, _, _ = rsf_train.fit_model_from_metadata(meta, train_df, valid_df, refit_train_valid=True)
        X = rsf_prepare.build_matrix_from_feature_names(eval_df, fn)
        X_tr = rsf_prepare.build_matrix_from_feature_names(eval_tr, fn)
        X_co = rsf_prepare.build_matrix_from_feature_names(eval_co, fn)
        risks.append(rsf_prepare.predict_rsf_risk(m, X))
        risks_tr.append(rsf_prepare.predict_rsf_risk(m, X_tr))
        risks_co.append(rsf_prepare.predict_rsf_risk(m, X_co))
        feat_names_out = fn
    return np.mean(risks, axis=0), np.mean(risks_tr, axis=0), np.mean(risks_co, axis=0), feat_names_out


def evaluate_rsf(run_dir: Path, refit_train_valid: bool = True, km_figure_path: Path | None = None) -> dict:
    run_dir = run_dir.resolve()
    metadata = json.loads((run_dir / "metadata.json").read_text())
    feature_names = metadata.get("feature_names") or [
        line.strip() for line in (run_dir / "feature_names.txt").read_text().splitlines() if line.strip()
    ]
    genes_main = list(metadata.get("genes_main", []))
    genes_inter = list(metadata.get("genes_inter", []))
    clin_cols = list(metadata["clin_cols"])
    test_df = _load_test(feature_names)

    # Detect whether this run used a seed ensemble (iter_011+).
    is_ensemble = len(metadata.get("result", {}).get("seed_panel_cis", [])) > 1

    if is_ensemble and refit_train_valid:
        print(f"[Finalize] Ensemble run detected ({len(FINAL_SEEDS)} seeds). Refitting on train+valid...")
        train_df, valid_df = rsf_prepare.load_train_valid()
        ens_r, ens_r_tr, ens_r_co, feature_names = _ensemble_predict(
            metadata, train_df, valid_df, test_df, seeds=FINAL_SEEDS
        )
        test_ci = rsf_prepare.cindex(
            ens_r, test_df["OS_MONTHS"].to_numpy(float), test_df["OS_STATUS"].to_numpy(int)
        )
        model_rec = np.where(ens_r_tr < ens_r_co, 1, 0)
        model_source = "ensemble_refit_train_valid"
        model_label = "rsf_ensemble"
    elif refit_train_valid:
        print("[Finalize] Single-model run. Refitting on train+valid...")
        train_df, valid_df = rsf_prepare.load_train_valid()
        model, feature_names, genes_main, genes_inter, clin_cols = rsf_train.fit_model_from_metadata(
            metadata, train_df, valid_df, refit_train_valid=True
        )
        test_tr = test_df.copy(); test_tr["Adjuvant Chemo"] = 1
        test_co = test_df.copy(); test_co["Adjuvant Chemo"] = 0
        X_test = rsf_prepare.build_matrix_from_feature_names(test_df, feature_names)
        X_tr = rsf_prepare.build_matrix_from_feature_names(test_tr, feature_names)
        X_co = rsf_prepare.build_matrix_from_feature_names(test_co, feature_names)
        test_ci = rsf_prepare.cindex(
            rsf_prepare.predict_rsf_risk(model, X_test),
            test_df["OS_MONTHS"].to_numpy(float), test_df["OS_STATUS"].to_numpy(int),
        )
        model_rec = np.where(
            rsf_prepare.predict_rsf_risk(model, X_tr) < rsf_prepare.predict_rsf_risk(model, X_co), 1, 0
        )
        model_source = "refit_train_valid"
        model_label = "rsf"
    else:
        print("[Finalize] Loading saved validation model (no refit)...")
        with open(run_dir / "rsf_model.pkl", "rb") as f:
            model = pickle.load(f)
        test_tr = test_df.copy(); test_tr["Adjuvant Chemo"] = 1
        test_co = test_df.copy(); test_co["Adjuvant Chemo"] = 0
        X_test = rsf_prepare.build_matrix_from_feature_names(test_df, feature_names)
        X_tr = rsf_prepare.build_matrix_from_feature_names(test_tr, feature_names)
        X_co = rsf_prepare.build_matrix_from_feature_names(test_co, feature_names)
        test_ci = rsf_prepare.cindex(
            rsf_prepare.predict_rsf_risk(model, X_test),
            test_df["OS_MONTHS"].to_numpy(float), test_df["OS_STATUS"].to_numpy(int),
        )
        model_rec = np.where(
            rsf_prepare.predict_rsf_risk(model, X_tr) < rsf_prepare.predict_rsf_risk(model, X_co), 1, 0
        )
        model_source = "saved_validation_model"
        model_label = "rsf"

    actual = test_df["Adjuvant Chemo"].to_numpy(int)
    alignment = actual == model_rec

    if int(alignment.sum()) > 0 and int((~alignment).sum()) > 0:
        km_a = KaplanMeierFitter().fit(
            test_df.loc[alignment, "OS_MONTHS"], event_observed=test_df.loc[alignment, "OS_STATUS"]
        )
        km_n = KaplanMeierFitter().fit(
            test_df.loc[~alignment, "OS_MONTHS"], event_observed=test_df.loc[~alignment, "OS_STATUS"]
        )
        test_rmst_diff = float(restricted_mean_survival_time(km_a, t=60) - restricted_mean_survival_time(km_n, t=60))
        logrank = logrank_test(
            test_df.loc[alignment, "OS_MONTHS"],
            test_df.loc[~alignment, "OS_MONTHS"],
            event_observed_A=test_df.loc[alignment, "OS_STATUS"],
            event_observed_B=test_df.loc[~alignment, "OS_STATUS"],
        )
        logrank_p = float(logrank.p_value)
    else:
        test_rmst_diff = 0.0
        logrank_p = float("nan")

    if km_figure_path is not None and int(alignment.sum()) > 0 and int((~alignment).sum()) > 0:
        _km_alignment_plot(test_df, alignment, km_figure_path, test_ci, test_rmst_diff, logrank_p)
        print(f"[Finalize] Wrote KM figure to {km_figure_path}")

    return {
        "model": model_label,
        "run_dir": str(run_dir.relative_to(REPO_ROOT)),
        "model_source": model_source,
        "test_ci": float(test_ci),
        "test_rmst_diff": float(test_rmst_diff),
        "km_logrank_pvalue": float(logrank_p),
        "km_figure": str(km_figure_path) if km_figure_path else "",
        "c_for_benefit": float("nan"),
        "c_for_benefit_note": "Placeholder pending a versioned Van Klaveren-style implementation.",
        "n_features": int(len(feature_names)),
        "k_main": int(len(genes_main)),
        "k_int": int(len(genes_inter)),
        **_gene_stability(run_dir, genes_main),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Human-only sealed test evaluation for RSF.")
    parser.add_argument(
        "--rsf-run-dir",
        type=Path,
        default=None,
        help="Run directory under rsf_arena/runs (default: most recently modified).",
    )
    parser.add_argument(
        "--no-refit-train-valid",
        action="store_true",
        help="Evaluate the saved validation model instead of refitting on train+valid before test.",
    )
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "final_rsf_comparison.csv")
    parser.add_argument(
        "--km-fig",
        type=Path,
        default=None,
        help="Path for the Kaplan-Meier treatment-alignment figure (PNG).",
    )
    args = parser.parse_args()
    run_dir = args.rsf_run_dir or latest_run_dir()
    km_figure_path = args.km_fig or args.out.with_name(f"{args.out.stem}_km.png")
    row = evaluate_rsf(
        run_dir.resolve(),
        refit_train_valid=not args.no_refit_train_valid,
        km_figure_path=km_figure_path.resolve(),
    )
    table = pd.DataFrame([row])
    print(table.to_string(index=False))
    table.to_csv(args.out, index=False)
    print(f"\nWrote RSF comparison table to {args.out}")
    if row.get("km_figure"):
        print(f"Wrote Kaplan-Meier figure to {row['km_figure']}")


if __name__ == "__main__":
    main()
