"""Human-only sealed test evaluation for RSF arena artifacts."""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

from rsf_arena import prepare as rsf_prepare
from rsf_arena import train as rsf_train

REPO_ROOT = Path(__file__).resolve().parent
TEST_CSV = REPO_ROOT / "affyfRMATest.csv"
RSF_RUNS_DIR = REPO_ROOT / "rsf_arena" / "runs"


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
        raise FileNotFoundError(f"{TEST_CSV} does not exist. Add the sealed test CSV locally, then rerun manually.")
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


def evaluate_rsf(run_dir: Path, refit_train_valid: bool = True) -> dict:
    run_dir = run_dir.resolve()
    metadata = json.loads((run_dir / "metadata.json").read_text())
    feature_names = metadata.get("feature_names") or [
        line.strip() for line in (run_dir / "feature_names.txt").read_text().splitlines() if line.strip()
    ]
    genes_main = list(metadata.get("genes_main", []))
    genes_inter = list(metadata.get("genes_inter", []))
    clin_cols = list(metadata["clin_cols"])
    test_df = _load_test(feature_names)

    if refit_train_valid:
        train_df, valid_df = rsf_prepare.load_train_valid()
        model, feature_names, genes_main, genes_inter, clin_cols = rsf_train.fit_model_from_metadata(
            metadata, train_df, valid_df, refit_train_valid=True
        )
        model_source = "refit_train_valid"
    else:
        with open(run_dir / "rsf_model.pkl", "rb") as f:
            model = pickle.load(f)
        model_source = "saved_validation_model"

    X_test = rsf_prepare.build_matrix_from_feature_names(test_df, feature_names)
    pred = rsf_prepare.predict_rsf_risk(model, X_test)
    test_ci = rsf_prepare.cindex(pred, test_df["OS_MONTHS"].to_numpy(float), test_df["OS_STATUS"].to_numpy(int))
    test_rmst_diff = rsf_prepare.compute_alignment_rmst_diff_rsf(
        model,
        test_df,
        genes_main=genes_main,
        genes_inter=genes_inter,
        dup_inter=int(metadata.get("dup_inter", 1)),
        feature_names=feature_names,
        clin_cols=clin_cols,
        tau=60,
    )
    return {
        "model": "rsf",
        "run_dir": str(run_dir.relative_to(REPO_ROOT)),
        "model_source": model_source,
        "test_ci": float(test_ci),
        "test_rmst_diff": float(test_rmst_diff),
        "c_for_benefit": np.nan,
        "c_for_benefit_note": "Placeholder pending a versioned Van Klaveren-style implementation.",
        "n_features": int(len(feature_names)),
        "k_main": int(len(genes_main)),
        "k_int": int(len(genes_inter)),
        **_gene_stability(run_dir, genes_main),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Human-only sealed test evaluation for RSF.")
    parser.add_argument("--rsf-run-dir", type=Path, default=None)
    parser.add_argument("--no-refit-train-valid", action="store_true")
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "final_rsf_comparison.csv")
    args = parser.parse_args()
    run_dir = args.rsf_run_dir or latest_run_dir()
    row = evaluate_rsf(run_dir, refit_train_valid=not args.no_refit_train_valid)
    table = pd.DataFrame([row])
    print(table.to_string(index=False))
    table.to_csv(args.out, index=False)
    print(f"\nWrote RSF comparison table to {args.out}")


if __name__ == "__main__":
    main()
