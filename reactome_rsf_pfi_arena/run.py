"""Fixed-budget launcher for Reactome-wide RSF PFI experiments."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

try:
    from . import integrity, prepare
except ImportError:
    import integrity
    import prepare


ARENA_DIR = Path(__file__).resolve().parent
TRAIN_PATH = ARENA_DIR / "train.py"
RUNS_DIR = ARENA_DIR / "runs"
LEDGER_PATH = ARENA_DIR / "search_ledger.jsonl"
BEST_PATH = ARENA_DIR / "best_run.txt"


def _experiments(rows: list[dict]) -> set[int]:
    return {int(row["experiment"]) for row in rows if "experiment" in row}


def _validate_result(result: dict, *, smoke: bool, artifact_dir: Path) -> None:
    required = {
        "schema_version", "phase", "candidate", "objective", "screening",
        "permutation_importance", "reduced_search", "final", "reward",
        "secondary_score", "eligible", "data", "budget", "artifacts",
    }
    missing = sorted(required - set(result))
    if missing:
        raise RuntimeError(f"Candidate result is missing keys: {missing}")
    if int(result["schema_version"]) != 1:
        raise RuntimeError("Candidate result uses the wrong schema")
    if bool(result["data"].get("test_used")):
        raise RuntimeError("Candidate reports test access")
    if bool(result["data"].get("confirmatory_validation_used")):
        raise RuntimeError("The former validation cohort must be adaptive development data")
    if bool(result["budget"].get("smoke_reductions_applied")) != bool(smoke):
        raise RuntimeError("Candidate result smoke marker does not match launcher mode")
    for key in (
        "max_experiments", "wall_time_seconds", "cv_folds", "screening_trials",
        "reduced_trials", "permutation_repeats", "max_top_genes",
        "max_trees_per_forest", "tau_months", "forest_seed", "fold_seed",
        "permutation_seed",
    ):
        if result["budget"].get(key) != prepare.BUDGET.get(key):
            raise RuntimeError(f"Candidate changed fixed budget field: {key}")
    expected_n = list(range(3)) if smoke else list(range(int(prepare.BUDGET["max_top_genes"]) + 1))
    if result["reduced_search"].get("evaluated_top_n") != expected_n:
        raise RuntimeError("Candidate did not evaluate every required top-N panel")
    if not smoke and int(result["screening"]["n_reactome_genes"]) != int(
        result["data"]["reactome_genes_before_smoke_reduction"]
    ):
        raise RuntimeError("Full screening did not use every eligible Reactome gene")
    for name, metadata in result["artifacts"].items():
        path = artifact_dir / name
        if not path.exists():
            raise RuntimeError(f"Candidate did not create artifact: {name}")
        if integrity.sha256_file(path) != metadata.get("sha256"):
            raise RuntimeError(f"Candidate artifact hash mismatch: {name}")


def _stream_process(command: list[str], timeout: int) -> tuple[int, str]:
    process = subprocess.Popen(
        command,
        cwd=prepare.REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    try:
        output, _ = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        process.terminate()
        try:
            output, _ = process.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            output, _ = process.communicate()
        output += "\n[launcher] wall-time budget exceeded\n"
        print(output, end="", flush=True)
        return 124, output
    print(output, end="", flush=True)
    return int(process.returncode), output


def _is_better(current: dict, prior: dict | None) -> bool:
    if prior is None:
        return True
    primary = float(current["reward"])
    prior_primary = float(prior["reward"])
    if primary > prior_primary + prepare.NUMERICAL_TOLERANCE:
        return True
    return (
        abs(primary - prior_primary) <= prepare.NUMERICAL_TOLERANCE
        and float(current["secondary_score"]) > float(prior["secondary_score"])
    )


def run(*, smoke: bool = False) -> dict:
    integrity.verify_lock()
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    train_sha = integrity.sha256_file(TRAIN_PATH)
    rows = integrity.read_ledger(LEDGER_PATH)
    if smoke:
        experiment = 0
    else:
        experiments = _experiments(rows)
        if len(experiments) >= int(prepare.BUDGET["max_experiments"]):
            raise RuntimeError("The fixed experiment budget is exhausted")
        completed_shas = {
            row.get("train_sha256") for row in rows if row.get("event") == "started"
        }
        if train_sha in completed_shas:
            raise RuntimeError("This exact train.py has already consumed an experiment")
        experiment = max(experiments, default=0) + 1
        integrity.append_ledger(LEDGER_PATH, {
            "event": "started",
            "experiment": experiment,
            "started_at_utc": datetime.now(timezone.utc).isoformat(),
            "train_sha256": train_sha,
        })

    work_dir = Path(tempfile.mkdtemp(prefix="candidate_", dir=RUNS_DIR))
    result_path = work_dir / "candidate_result.json"
    command = [
        sys.executable,
        "-m",
        "reactome_rsf_pfi_arena.train",
        "--artifact-dir",
        str(work_dir),
        "--result-path",
        str(result_path),
    ]
    if smoke:
        command.append("--smoke")
    started = time.monotonic()
    return_code, output = _stream_process(
        command, int(prepare.BUDGET["wall_time_seconds"])
    )
    elapsed = time.monotonic() - started
    if return_code != 0:
        if not smoke:
            integrity.append_ledger(LEDGER_PATH, {
                "event": "finished",
                "experiment": experiment,
                "status": "failed",
                "return_code": return_code,
                "wall_clock_seconds": elapsed,
                "train_sha256": train_sha,
            })
        raise RuntimeError(
            f"Candidate failed with return code {return_code}; work retained at {work_dir}"
        )
    result = json.loads(result_path.read_text(encoding="utf-8"))
    _validate_result(result, smoke=smoke, artifact_dir=work_dir)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    if smoke:
        destination = RUNS_DIR / f"smoke_{timestamp}"
        work_dir.replace(destination)
        print(f"\n[launcher] saved {destination}", flush=True)
        return result

    run_id = f"run_{experiment:03d}_{timestamp}"
    destination = RUNS_DIR / run_id
    result["run_id"] = run_id
    result["wall_clock_seconds"] = elapsed
    result_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (work_dir / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (work_dir / "stdout.log").write_text(output, encoding="utf-8")
    shutil.copy2(TRAIN_PATH, work_dir / "train_snapshot.py")
    work_dir.replace(destination)
    integrity.append_ledger(LEDGER_PATH, {
        "event": "finished",
        "experiment": experiment,
        "status": "completed",
        "run_id": run_id,
        "reward": float(result["reward"]),
        "secondary_score": float(result["secondary_score"]),
        "wall_clock_seconds": elapsed,
        "train_sha256": train_sha,
    })

    prior = None
    if BEST_PATH.exists() and BEST_PATH.read_text(encoding="utf-8").strip():
        prior_name = BEST_PATH.read_text(encoding="utf-8").strip()
        prior_path = RUNS_DIR / prior_name / "result.json"
        if prior_path.exists():
            prior = json.loads(prior_path.read_text(encoding="utf-8"))
    if _is_better(result, prior):
        BEST_PATH.write_text(run_id + "\n", encoding="utf-8")
    print(f"\n[launcher] saved {destination}", flush=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--smoke", action="store_true",
        help="Reduced infrastructure check; does not consume budget",
    )
    args = parser.parse_args()
    run(smoke=args.smoke)


if __name__ == "__main__":
    main()

