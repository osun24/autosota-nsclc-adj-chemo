"""Human-only sealed test evaluation for DeepSurv arena artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from deepsurv_arena import prepare as ds_prepare
from deepsurv_arena import train as ds_train

REPO_ROOT = Path(__file__).resolve().parent
TEST_CSV = REPO_ROOT / "affyfRMATest.csv"
DEEPSURV_RUNS_DIR = REPO_ROOT / "deepsurv_arena" / "runs"


def latest_run_dir(runs_dir: Path = DEEPSURV_RUNS_DIR) -> Path:
    candidates = [p for p in runs_dir.iterdir() if p.is_dir() and (p / "metadata.json").exists()]
    if not candidates:
        raise FileNotFoundError(f"No DeepSurv run directories with metadata.json found in {runs_dir}")
    return sorted(candidates, key=lambda p: p.stat().st_mtime, reverse=True)[0]


def _needed_gene_names(feature_names: list[str]) -> list[str]:
    genes = []
    for feat in feature_names:
        base = feat.split("*ACT", 1)[0] if "*ACT" in feat else feat
        if base not in ds_prepare.CLINICAL_VARS and base not in {"OS_STATUS", "OS_MONTHS"}:
            genes.append(base)
    return sorted(set(genes))


def _load_test(feature_names: list[str]) -> pd.DataFrame:
    if not TEST_CSV.exists():
        raise FileNotFoundError(f"{TEST_CSV} does not exist. Add the sealed test CSV locally, then rerun manually.")
    test_raw = pd.read_csv(TEST_CSV)
    test_df = ds_prepare.preprocess_split(test_raw, ds_prepare.CLINICAL_VARS, _needed_gene_names(feature_names))
    return test_df.sort_values(["OS_MONTHS", "OS_STATUS"], ascending=[False, False]).reset_index(drop=True)


def _gene_stability(run_dir: Path, genes_main: list[str]) -> dict:
    current = set(genes_main)
    overlaps = []
    for meta_path in DEEPSURV_RUNS_DIR.glob("*/metadata.json"):
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


def evaluate_deepsurv(run_dir: Path, refit_train_valid: bool = True) -> dict:
    run_dir = run_dir.resolve()
    metadata = json.loads((run_dir / "metadata.json").read_text())
    feature_names = metadata.get("feature_names") or [
        line.strip() for line in (run_dir / "feature_names.txt").read_text().splitlines() if line.strip()
    ]
    genes_main = list(metadata.get("genes_main", []))
    genes_inter = list(metadata.get("genes_inter", []))
    clin_cols = list(metadata["clin_cols"])
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    test_df = _load_test(feature_names)

    if refit_train_valid:
        train_df, valid_df = ds_prepare.load_train_valid()
        model, transform, feature_names, genes_main, genes_inter, clin_cols = ds_train.fit_model_from_metadata(
            metadata, train_df, valid_df, refit_train_valid=True, device=device
        )
        model_source = "refit_train_valid"
    else:
        transform_npz = np.load(run_dir / "tabular_transform.npz")
        transform = {k: transform_npz[k] for k in transform_npz.files}
        params = metadata["deepsurv_params"]
        model = ds_train.DeepSurvMLP(len(feature_names), params["hidden_layers"], params["dropout"]).to(device)
        model.load_state_dict(torch.load(run_dir / "deepsurv_model.pt", map_location=device))
        model.eval()
        model_source = "saved_validation_model"

    X_test = ds_prepare.apply_tabular_transform(ds_prepare.build_matrix_from_feature_names(test_df, feature_names), transform)
    pred = ds_prepare.predict_deepsurv_risk(model, X_test, device)
    test_ci = ds_prepare.cindex(pred, test_df["OS_MONTHS"].to_numpy(float), test_df["OS_STATUS"].to_numpy(int))
    test_rmst_diff = ds_prepare.compute_alignment_rmst_diff_deepsurv(model, test_df, feature_names, transform, device, tau=60)
    return {
        "model": "deepsurv",
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
    parser = argparse.ArgumentParser(description="Human-only sealed test evaluation for DeepSurv.")
    parser.add_argument("--deepsurv-run-dir", type=Path, default=None)
    parser.add_argument("--no-refit-train-valid", action="store_true")
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "final_deepsurv_comparison.csv")
    args = parser.parse_args()
    run_dir = args.deepsurv_run_dir or latest_run_dir()
    row = evaluate_deepsurv(run_dir, refit_train_valid=not args.no_refit_train_valid)
    table = pd.DataFrame([row])
    print(table.to_string(index=False))
    table.to_csv(args.out, index=False)
    print(f"\nWrote DeepSurv comparison table to {args.out}")


if __name__ == "__main__":
    main()
