"""Fixed-budget launcher for agent-edited Reactome RSF v2 experiments."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
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
DIAGNOSTIC_PATH = ARENA_DIR / "diagnostic_leader.txt"


def _experiments(rows: list[dict]) -> set[int]:
    return {int(row["experiment"]) for row in rows if "experiment" in row}


def _validate_result(result: dict, *, smoke: bool) -> None:
    required = {
        "schema_version", "phase", "candidate", "reward", "eligible",
        "gates", "generalization", "development_oof", "repeat_results",
        "gene_selection", "data", "budget",
    }
    missing = sorted(required - set(result))
    if missing:
        raise RuntimeError(f"Candidate result is missing keys: {missing}")
    if int(result["schema_version"]) != 2:
        raise RuntimeError("Candidate result does not use the v2 schema")
    if bool(result["data"].get("test_used")):
        raise RuntimeError("Candidate reports test access")
    if bool(result["data"].get("confirmatory_validation_used")):
        raise RuntimeError("V2 has no confirmatory validation cohort")
    if bool(result["budget"].get("smoke_reductions_applied")) != bool(smoke):
        raise RuntimeError("Candidate result smoke marker does not match launcher mode")
    expected_repeats = 1 if smoke else int(prepare.BUDGET["outer_repeats"])
    if len(result["repeat_results"]) != expected_repeats:
        raise RuntimeError("Candidate returned the wrong number of outer repeats")
    for key in (
        "max_experiments", "wall_time_seconds", "outer_folds",
        "outer_repeats", "forest_seeds", "clinical_forest",
    ):
        if result["budget"].get(key) != prepare.BUDGET.get(key):
            raise RuntimeError(f"Candidate changed fixed budget field: {key}")
    panels = [("development_oof", result["development_oof"])] + [
        (f"repeat_{item.get('repeat')}", item) for item in result["repeat_results"]
    ]
    for cohort, panel in panels:
        models = panel.get("models", {})
        for model_name in ("clinical", "clinical_plus_genomic"):
            usage = models.get(model_name, {}).get("act_usage", {})
            required_usage = {
                "tree_contains_act_split_fraction",
                "patient_tree_paths_traversing_act_fraction",
                "per_patient_different_terminal_tree_fraction",
            }
            if required_usage - set(usage):
                raise RuntimeError(
                    f"{cohort}/{model_name} is missing locked ACT-use diagnostics"
                )
            terminal = usage["per_patient_different_terminal_tree_fraction"]
            if {"mean", "median", "p10", "p90", "nonzero_patient_fraction"} - set(terminal):
                raise RuntimeError(
                    f"{cohort}/{model_name} has an incomplete ACT terminal-node distribution"
                )


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
            row.get("train_sha256") for row in rows
            if row.get("event") == "started"
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

    handle, temp_name = tempfile.mkstemp(prefix="candidate_", suffix=".json", dir=RUNS_DIR)
    os.close(handle)
    result_path = Path(temp_name)
    command = [
        sys.executable, "-m", "reactome_rsf_arena_v2.train",
        "--result-path", str(result_path),
    ]
    if smoke:
        command.append("--smoke")
    started = time.monotonic()
    return_code, output = _stream_process(command, int(prepare.BUDGET["wall_time_seconds"]))
    elapsed = time.monotonic() - started
    if return_code != 0:
        if not smoke:
            integrity.append_ledger(LEDGER_PATH, {
                "event": "finished", "experiment": experiment,
                "status": "failed", "return_code": return_code,
                "wall_clock_seconds": elapsed, "train_sha256": train_sha,
            })
        raise RuntimeError(f"Candidate failed with return code {return_code}")
    result = json.loads(result_path.read_text(encoding="utf-8"))
    _validate_result(result, smoke=smoke)
    if smoke:
        smoke_path = RUNS_DIR / f"smoke_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}.json"
        result_path.replace(smoke_path)
        return result

    run_id = f"run_{experiment:03d}_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
    run_dir = RUNS_DIR / run_id
    run_dir.mkdir(parents=False, exist_ok=False)
    result["run_id"] = run_id
    result["wall_clock_seconds"] = elapsed
    (run_dir / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (run_dir / "stdout.log").write_text(output, encoding="utf-8")
    (run_dir / "train_snapshot.py").write_bytes(TRAIN_PATH.read_bytes())
    result_path.replace(run_dir / "candidate_result.json")
    integrity.append_ledger(LEDGER_PATH, {
        "event": "finished", "experiment": experiment, "status": "completed",
        "run_id": run_id, "reward": float(result["reward"]),
        "eligible": bool(result["eligible"]), "wall_clock_seconds": elapsed,
        "train_sha256": train_sha,
    })

    prior_best = None
    if BEST_PATH.exists() and BEST_PATH.read_text(encoding="utf-8").strip():
        prior_name = BEST_PATH.read_text(encoding="utf-8").strip()
        prior_file = RUNS_DIR / prior_name / "result.json"
        if prior_file.exists():
            prior_best = json.loads(prior_file.read_text(encoding="utf-8"))
    if bool(result["eligible"]) and (
        prior_best is None or float(result["reward"]) > float(prior_best["reward"])
    ):
        BEST_PATH.write_text(run_id + "\n", encoding="utf-8")

    prior_diagnostic = None
    if DIAGNOSTIC_PATH.exists() and DIAGNOSTIC_PATH.read_text(encoding="utf-8").strip():
        prior_name = DIAGNOSTIC_PATH.read_text(encoding="utf-8").strip()
        prior_file = RUNS_DIR / prior_name / "result.json"
        if prior_file.exists():
            prior_diagnostic = json.loads(prior_file.read_text(encoding="utf-8"))
    current_score = float(result["generalization"]["score_before_eligibility"])
    prior_score = (
        float(prior_diagnostic["generalization"]["score_before_eligibility"])
        if prior_diagnostic is not None else float("-inf")
    )
    if current_score > prior_score:
        DIAGNOSTIC_PATH.write_text(run_id + "\n", encoding="utf-8")
    print(f"\n[launcher] saved {run_dir}", flush=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", action="store_true", help="Reduced infrastructure check; does not consume budget")
    args = parser.parse_args()
    run(smoke=args.smoke)


if __name__ == "__main__":
    main()
