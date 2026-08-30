"""Human-only sealed-test evaluation for the frozen Reactome RSF T-learner.

This script is deliberately separate from the autonomous arena. It loads the
sealed ``affyfRMATest.csv`` only after the human operator passes
``--confirm-sealed-test``. Do not import or execute it from an autonomous
search loop, and do not rerun it after inspecting test results.

The frozen run-20 clinical-plus-genomic T-learner is refit on the pooled
development cohort once for every base random seed from 0 through 50. For
each seed, the ACT-arm forest uses the arena's original ``seed + 100_003``
convention. The reported observed alignment contrast is

    KM RMST(aligned with recommendation, 60 months)
    - KM RMST(not aligned with recommendation, 60 months).

This observed contrast is descriptive and unadjusted; it is not the arena's
cross-fitted IPCW-AIPW estimand and remains vulnerable to confounding.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gc
import hashlib
from importlib import metadata as importlib_metadata
import json
import os
from pathlib import Path
import platform
import sys

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
from sksurv.metrics import concordance_index_censored

from reactome_tlearner_arena import integrity, prepare


REPO_ROOT = Path(__file__).resolve().parent
TEST_CSV = REPO_ROOT / "affyfRMATest.csv"
ARENA_DIR = REPO_ROOT / "reactome_tlearner_arena"
FROZEN_RUN_ID = "run_020_20260824T021215Z"
FROZEN_RUN_DIR = ARENA_DIR / "runs" / FROZEN_RUN_ID
FROZEN_RESULT_SHA256 = "9f81e80f2530890740142df05c841d1ead74695352e4ea114872a1aa36d15087"
FROZEN_SNAPSHOT_SHA256 = "5e47a33629d2d77d3fe81f6b9e49b6ebed2694a28d22c8930738c8a29df65996"
FINAL_SEEDS = tuple(range(51))
TAU_MONTHS = 60.0

FROZEN_GENES = (
    "OR2S2",
    "TAAR5",
    "ANO2",
    "OR10H2",
    "HTR6",
    "OR2J2",
    "OR3A2",
    "OR52A1",
    "FGF22",
    "OR10C1",
    "FGF4",
    "OR2F1",
    "TAS2R13",
    "ADRA2B",
    "FGF16",
    "CASR",
)

FROZEN_CANDIDATE = {
    "name": "tlearner_low_floor_module16",
    "selector": "dr_gene",
    "n_genes": 16,
    "representation": "module",
    "module_count": 1,
    "benefit_threshold_months": 0.0,
    "tlearner": {
        "observation": {
            "n_estimators": 1000,
            "max_depth": 9,
            "min_samples_leaf": 8,
            "min_samples_split": 16,
            "max_features": 1.0,
        },
        "act": {
            "n_estimators": 1000,
            "max_depth": 7,
            "min_samples_leaf": 12,
            "min_samples_split": 24,
            "max_features": 1.0,
        },
    },
}


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_frozen_contract(run_dir: Path = FROZEN_RUN_DIR) -> tuple[dict, prepare.CandidateSpec]:
    """Verify that the requested artifact is exactly the frozen run-20 winner."""
    integrity.verify_lock()
    run_dir = run_dir.resolve()
    if run_dir != FROZEN_RUN_DIR.resolve():
        raise RuntimeError(f"Only frozen {FROZEN_RUN_ID} may be finalized: {run_dir}")

    best_run = (ARENA_DIR / "best_run.txt").read_text(encoding="utf-8").strip()
    if best_run != FROZEN_RUN_ID:
        raise RuntimeError(f"best_run.txt names {best_run!r}, expected {FROZEN_RUN_ID!r}")

    result_path = run_dir / "result.json"
    snapshot_path = run_dir / "train_snapshot.py"
    if _sha256_file(result_path) != FROZEN_RESULT_SHA256:
        raise RuntimeError("Frozen run-20 result.json hash does not match the finalization contract")
    if _sha256_file(snapshot_path) != FROZEN_SNAPSHOT_SHA256:
        raise RuntimeError("Frozen run-20 train snapshot hash does not match the finalization contract")

    result = json.loads(result_path.read_text(encoding="utf-8"))
    if result.get("run_id") != FROZEN_RUN_ID:
        raise RuntimeError("Frozen result has the wrong run_id")
    if not bool(result.get("eligible")) or not all(result.get("gates", {}).values()):
        raise RuntimeError("Frozen run-20 result is no longer eligible on every gate")
    if result.get("candidate") != FROZEN_CANDIDATE:
        raise RuntimeError("Frozen run-20 candidate specification changed")
    if tuple(result.get("gene_selection", {}).get("full_development_selection", [])) != FROZEN_GENES:
        raise RuntimeError("Frozen run-20 full-development gene panel changed")
    if float(result.get("budget", {}).get("tau_months", np.nan)) != TAU_MONTHS:
        raise RuntimeError("Frozen run-20 RMST horizon is not 60 months")
    if bool(result.get("data", {}).get("test_used")):
        raise RuntimeError("Frozen search artifact unexpectedly reports prior test use")
    return result, prepare.validate_candidate(result["candidate"])


def _required_columns() -> list[str]:
    return [*prepare.OUTCOME_COLUMNS, *prepare.CLINICAL_COLUMNS, *FROZEN_GENES]


def _read_analysis_frame(path: Path) -> pd.DataFrame:
    required = _required_columns()
    header = set(pd.read_csv(path, nrows=0).columns)
    missing = sorted(set(required) - header)
    if missing:
        raise ValueError(f"{path.name} is missing frozen-pipeline columns: {missing}")
    raw = pd.read_csv(path, usecols=required)
    # The locked search preparer performs the same ACT mapping and numeric
    # validation. Test genes may be missing and are imputed from development.
    return prepare.base._prepare_frame(raw)[required]


def _load_development(result: dict) -> pd.DataFrame:
    expected = result["data"]
    if _sha256_file(prepare.TRAIN_CSV) != expected["train_sha256"]:
        raise RuntimeError("Development train CSV no longer matches frozen run 20")
    if _sha256_file(prepare.VALID_CSV) != expected["validation_sha256"]:
        raise RuntimeError("Development validation CSV no longer matches frozen run 20")
    train = _read_analysis_frame(prepare.TRAIN_CSV)
    valid = _read_analysis_frame(prepare.VALID_CSV)
    development = pd.concat([train, valid], axis=0, ignore_index=True)
    if len(development) != int(expected["development_n"]):
        raise RuntimeError("Pooled development row count changed")
    if int(development["OS_STATUS"].sum()) != int(expected["development_events"]):
        raise RuntimeError("Pooled development event count changed")
    if int(development[prepare.TREATMENT].sum()) != int(expected["development_act"]):
        raise RuntimeError("Pooled development ACT count changed")
    return development


def _load_sealed_test(test_csv: Path) -> pd.DataFrame:
    test_csv = test_csv.resolve()
    if test_csv in {prepare.TRAIN_CSV.resolve(), prepare.VALID_CSV.resolve()}:
        raise ValueError("The sealed-test input cannot be a development CSV")
    if not test_csv.exists():
        raise FileNotFoundError(
            f"{test_csv} does not exist. Add the sealed test CSV locally, then rerun manually."
        )
    test = _read_analysis_frame(test_csv)
    if len(test) == 0:
        raise ValueError("The sealed test cohort is empty")
    if set(test[prepare.TREATMENT].unique()) != {0, 1}:
        raise ValueError("The sealed test cohort must contain both observation and ACT patients")
    if set(test["OS_STATUS"].unique()) != {0, 1}:
        raise ValueError("The sealed test cohort must contain events and censoring")
    return test


def _format_pvalue(pvalue: float) -> str:
    if not np.isfinite(pvalue):
        return "NA"
    if pvalue < 0.0001:
        return f"{pvalue:.2e}"
    return f"{pvalue:.4f}"


def _alignment_summary(frame: pd.DataFrame, alignment: np.ndarray, tau: float = TAU_MONTHS) -> dict:
    aligned = np.asarray(alignment, dtype=bool)
    if aligned.shape != (len(frame),):
        raise ValueError("Alignment vector has the wrong shape")
    n_aligned = int(aligned.sum())
    n_not_aligned = int((~aligned).sum())
    if n_aligned == 0 or n_not_aligned == 0:
        return {
            "aligned_n": n_aligned,
            "not_aligned_n": n_not_aligned,
            "aligned_rmst_60_months": np.nan,
            "not_aligned_rmst_60_months": np.nan,
            "alignment_rmst_difference_60_months": np.nan,
            "km_logrank_pvalue": np.nan,
        }

    km_aligned = KaplanMeierFitter().fit(
        frame.loc[aligned, "OS_MONTHS"],
        event_observed=frame.loc[aligned, "OS_STATUS"],
    )
    km_not_aligned = KaplanMeierFitter().fit(
        frame.loc[~aligned, "OS_MONTHS"],
        event_observed=frame.loc[~aligned, "OS_STATUS"],
    )
    rmst_aligned = float(restricted_mean_survival_time(km_aligned, t=tau))
    rmst_not_aligned = float(restricted_mean_survival_time(km_not_aligned, t=tau))
    logrank = logrank_test(
        frame.loc[aligned, "OS_MONTHS"],
        frame.loc[~aligned, "OS_MONTHS"],
        event_observed_A=frame.loc[aligned, "OS_STATUS"],
        event_observed_B=frame.loc[~aligned, "OS_STATUS"],
    )
    return {
        "aligned_n": n_aligned,
        "not_aligned_n": n_not_aligned,
        "aligned_rmst_60_months": rmst_aligned,
        "not_aligned_rmst_60_months": rmst_not_aligned,
        "alignment_rmst_difference_60_months": rmst_aligned - rmst_not_aligned,
        "km_logrank_pvalue": float(logrank.p_value),
    }


def _fit_predict_seed(
    development: pd.DataFrame,
    test: pd.DataFrame,
    spec: prepare.CandidateSpec,
    seed: int,
) -> tuple[dict, np.ndarray]:
    model0, model1, transformer, support = prepare._fit_tlearner(
        development,
        list(FROZEN_GENES),
        spec,
        int(seed),
        clinical_only=False,
        smoke=False,
    )
    risk, rmst0, rmst1, _ = prepare._tlearner_predictions(
        model0, model1, transformer, test, support
    )
    benefit = rmst1 - rmst0
    recommendation = prepare.recommendations(benefit, spec.benefit_threshold_months)
    actual = test[prepare.TREATMENT].to_numpy(dtype=int)
    alignment = actual == recommendation
    cindex = float(
        concordance_index_censored(
            test["OS_STATUS"].astype(bool).to_numpy(),
            test["OS_MONTHS"].to_numpy(dtype=float),
            risk,
        )[0]
    )
    summary = {
        "seed": int(seed),
        "observation_forest_random_state": int(seed),
        "act_forest_random_state": int(seed) + 100_003,
        "test_cindex": cindex,
        "act_recommended_fraction": float(np.mean(recommendation)),
        "predicted_benefit_mean_months": float(np.mean(benefit)),
        "predicted_benefit_iqr_months": float(np.percentile(benefit, 75) - np.percentile(benefit, 25)),
        **_alignment_summary(test, alignment),
    }
    del model0, model1, transformer, risk, rmst0, rmst1, benefit, recommendation
    gc.collect()
    return summary, alignment


def _plot_histogram(seed_results: pd.DataFrame, path: Path) -> float:
    values = seed_results["alignment_rmst_difference_60_months"].to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise RuntimeError("No seed produced a finite observed alignment RMST difference")
    median = float(np.median(finite))
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(10, 6))
    plt.hist(finite, bins="auto", color="#4472C4", edgecolor="white", linewidth=1.0)
    plt.axvline(median, color="#C00000", linestyle="--", linewidth=2, label=f"Median: {median:.2f} months")
    plt.title("Frozen Run-20 T-learner: Seed Robustness on the Sealed Test Set")
    plt.xlabel("Observed alignment RMST difference at 60 months (months)")
    plt.ylabel("Number of random seeds")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    return median


def _plot_median_km(
    frame: pd.DataFrame,
    alignment: np.ndarray,
    seed_result: dict,
    path: Path,
) -> None:
    aligned = np.asarray(alignment, dtype=bool)
    if int(aligned.sum()) == 0 or int((~aligned).sum()) == 0:
        raise RuntimeError("The median seed does not define both KM alignment groups")
    km_aligned = KaplanMeierFitter().fit(
        frame.loc[aligned, "OS_MONTHS"],
        event_observed=frame.loc[aligned, "OS_STATUS"],
        label="Aligned with run-20 recommendation",
    )
    km_not_aligned = KaplanMeierFitter().fit(
        frame.loc[~aligned, "OS_MONTHS"],
        event_observed=frame.loc[~aligned, "OS_STATUS"],
        label="Not aligned with run-20 recommendation",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(10, 6))
    ax = km_aligned.plot(ci_show=True)
    km_not_aligned.plot(ax=ax, ci_show=True)
    seed = int(seed_result["seed"])
    difference = float(seed_result["alignment_rmst_difference_60_months"])
    pvalue = float(seed_result["km_logrank_pvalue"])
    ax.set_title(f"Treatment Alignment for Median Seed Result (base seed {seed})")
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.text(0.06, 0.18, f"Log-rank p-value: {_format_pvalue(pvalue)}", transform=ax.transAxes)
    ax.text(0.06, 0.12, f"Observed 60-month RMST difference: {difference:.2f} months", transform=ax.transAxes)
    ax.text(
        0.06,
        0.06,
        f"Aligned n={int(seed_result['aligned_n'])}; not aligned n={int(seed_result['not_aligned_n'])}",
        transform=ax.transAxes,
    )
    add_at_risk_counts(km_aligned, km_not_aligned, ax=ax)
    plt.tight_layout()
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()


def _package_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for distribution in (
        "numpy",
        "pandas",
        "matplotlib",
        "lifelines",
        "scikit-learn",
        "scikit-survival",
    ):
        try:
            versions[distribution] = importlib_metadata.version(distribution)
        except importlib_metadata.PackageNotFoundError:
            versions[distribution] = "not-installed"
    return versions


def _display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path.resolve())


def evaluate_frozen_run(
    *,
    test_csv: Path,
    run_dir: Path,
    out_dir: Path,
) -> dict:
    """Run the complete prespecified 51-seed sealed-test evaluation once."""
    result, spec = _load_frozen_contract(run_dir)
    out_dir = out_dir.resolve()
    seed_csv = out_dir / "seed_results.csv"
    histogram_path = out_dir / "alignment_rmst_seed_histogram.png"
    km_path = out_dir / "median_alignment_kaplan_meier.png"
    manifest_path = out_dir / "finalization_manifest.json"
    occupied = [path for path in (seed_csv, histogram_path, km_path, manifest_path) if path.exists()]
    if occupied:
        raise FileExistsError(
            "Refusing to overwrite prior sealed-test outputs: " + ", ".join(map(str, occupied))
        )

    # Check for an earlier completed evaluation before the sealed CSV is opened.
    development = _load_development(result)
    test = _load_sealed_test(test_csv)

    rows: list[dict] = []
    alignments: dict[int, np.ndarray] = {}
    for index, seed in enumerate(FINAL_SEEDS, start=1):
        print(f"[Finalize] Fitting frozen run 20 for base seed {seed} ({index}/{len(FINAL_SEEDS)})...", flush=True)
        row, alignment = _fit_predict_seed(development, test, spec, seed)
        rows.append(row)
        alignments[int(seed)] = alignment

    table = pd.DataFrame(rows).sort_values("seed").reset_index(drop=True)
    if table["seed"].tolist() != list(FINAL_SEEDS):
        raise RuntimeError("Seed evaluation did not cover every integer from 0 through 50")
    finite = table[np.isfinite(table["alignment_rmst_difference_60_months"])].copy()
    if finite.empty:
        raise RuntimeError("Every seed produced a degenerate alignment comparison")
    distribution_median = float(np.median(finite["alignment_rmst_difference_60_months"]))
    finite["distance_from_distribution_median"] = np.abs(
        finite["alignment_rmst_difference_60_months"] - distribution_median
    )
    median_row = finite.sort_values(
        ["distance_from_distribution_median", "seed"], kind="mergesort"
    ).iloc[0].drop(labels="distance_from_distribution_median").to_dict()
    median_row = {
        key: value.item() if isinstance(value, np.generic) else value
        for key, value in median_row.items()
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    histogram_median = _plot_histogram(table, histogram_path)
    if not np.isclose(histogram_median, distribution_median, rtol=0.0, atol=1e-12):
        raise RuntimeError("Histogram median does not match the seed-result distribution")
    _plot_median_km(test, alignments[int(median_row["seed"])], median_row, km_path)
    table.to_csv(seed_csv, index=False)

    manifest = {
        "schema_version": 1,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "warning": (
            "Human-only one-shot sealed-test evaluation. Observed KM alignment RMST is descriptive, "
            "unadjusted, and not the development IPCW-AIPW estimand."
        ),
        "frozen_run_id": FROZEN_RUN_ID,
        "frozen_result_sha256": FROZEN_RESULT_SHA256,
        "frozen_train_snapshot_sha256": FROZEN_SNAPSHOT_SHA256,
        "candidate": FROZEN_CANDIDATE,
        "genes": list(FROZEN_GENES),
        "seeds": list(FINAL_SEEDS),
        "tau_months": TAU_MONTHS,
        "test": {
            "path": _display_path(test_csv),
            "sha256": _sha256_file(test_csv.resolve()),
            "n": int(len(test)),
            "events": int(test["OS_STATUS"].sum()),
            "act": int(test[prepare.TREATMENT].sum()),
        },
        "median_alignment_rmst_difference_60_months": distribution_median,
        "median_seed_result": median_row,
        "median_seed_selection": "Seed whose observed alignment RMST difference is closest to the distribution median; ties use the smaller seed.",
        "finite_seed_results": int(len(finite)),
        "artifacts": {
            "seed_results_csv": _display_path(seed_csv),
            "histogram_png": _display_path(histogram_path),
            "median_kaplan_meier_png": _display_path(km_path),
        },
        "environment": {
            "python": sys.version,
            "executable": sys.executable,
            "platform": platform.platform(),
            "packages": _package_versions(),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Human-only 51-seed sealed-test evaluation for frozen Reactome RSF T-learner run 20."
    )
    parser.add_argument("--test-csv", type=Path, default=TEST_CSV)
    parser.add_argument("--run-dir", type=Path, default=FROZEN_RUN_DIR)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=REPO_ROOT / "final_rsf_tlearner_test",
        help="New directory for the seed table, histogram, KM plot, and manifest.",
    )
    parser.add_argument(
        "--confirm-sealed-test",
        action="store_true",
        help="Required acknowledgement that this is the human-authorized one-shot sealed-test evaluation.",
    )
    args = parser.parse_args()
    if not args.confirm_sealed_test:
        parser.error(
            "Refusing to open sealed data without --confirm-sealed-test. "
            "Review the frozen contract and output paths first."
        )

    manifest = evaluate_frozen_run(
        test_csv=args.test_csv,
        run_dir=args.run_dir,
        out_dir=args.out_dir,
    )
    median = manifest["median_seed_result"]
    print("\n[Finalize] Completed the one-shot 51-seed evaluation.")
    print(
        "[Finalize] Median observed 60-month alignment RMST difference: "
        f"{float(manifest['median_alignment_rmst_difference_60_months']):.3f} months."
    )
    print(
        "[Finalize] Kaplan-Meier uses the closest seed result: "
        f"base seed {int(median['seed'])}, "
        f"{float(median['alignment_rmst_difference_60_months']):.3f} months."
    )
    for label, path in manifest["artifacts"].items():
        print(f"[Finalize] {label}: {path}")


if __name__ == "__main__":
    main()
