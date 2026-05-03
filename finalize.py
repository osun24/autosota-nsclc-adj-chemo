"""Human-only sealed test evaluation for arena artifacts.

Run this manually after `affyfRMATest.csv` has been added locally. The
autonomous loop must never import or execute this file.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb

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


def evaluate_xgb(run_dir: Path, refit_train_valid: bool) -> dict:
    metadata = json.loads((run_dir / "metadata.json").read_text())
    feature_names = metadata.get("feature_names") or _read_lines(run_dir / "feature_names.txt")
    genes_main = metadata.get("genes_main") or _read_lines(run_dir / "genes_main.txt")
    genes_inter = metadata.get("genes_inter") or _read_lines(run_dir / "genes_inter.txt")
    clin_cols = metadata["clin_cols"]

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
        best_ntree = int(metadata.get("best_ntree", metadata["result"].get("best_ntree", 0)))
        model_source = "saved_validation_model"

    X_test = xgb_prepare.build_matrix_from_feature_names(test_df, feature_names)
    pred = xgb_prepare.predict_xgb_risk(booster, X_test, feature_names, best_ntree)
    test_ci = xgb_prepare.cindex(
        pred,
        test_df["OS_MONTHS"].to_numpy(dtype=float),
        test_df["OS_STATUS"].to_numpy(dtype=int),
    )
    test_rmst_diff = xgb_prepare.compute_alignment_rmst_diff_xgb(
        booster,
        test_df,
        genes_main=genes_main,
        genes_inter=genes_inter,
        dup_inter=int(metadata.get("dup_inter", 1)),
        feature_names=feature_names,
        best_ntree=best_ntree,
        clin_cols=clin_cols,
        tau=60,
    )
    stability = _xgb_gene_stability(run_dir, genes_main)
    return {
        "model": "xgb",
        "run_dir": str(run_dir.relative_to(REPO_ROOT)),
        "model_source": model_source,
        "test_ci": float(test_ci),
        "test_rmst_diff": float(test_rmst_diff),
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
    args = parser.parse_args()

    xgb_run_dir = args.xgb_run_dir or _latest_run_dir(XGB_RUNS_DIR)
    row = evaluate_xgb(xgb_run_dir.resolve(), refit_train_valid=not args.no_refit_train_valid)
    table = pd.DataFrame([row])
    print(table.to_string(index=False))
    table.to_csv(args.out, index=False)
    print(f"\nWrote comparison table to {args.out}")


if __name__ == "__main__":
    main()
