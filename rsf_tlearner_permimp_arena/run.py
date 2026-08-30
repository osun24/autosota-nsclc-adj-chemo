"""Human-operated, cumulative-budget launcher for RSF T-learner PFI experiments."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

from . import integrity, prepare


ARENA_DIR = Path(__file__).resolve().parent
TRAIN_PATH = ARENA_DIR / "train.py"
RUNS_DIR = ARENA_DIR / "runs"
LEDGER_PATH = ARENA_DIR / "search_ledger.jsonl"
BEST_PATH = ARENA_DIR / "best_run.txt"
DIAGNOSTIC_PATH = ARENA_DIR / "diagnostic_leader.txt"
NOMINEE_PATH = ARENA_DIR / "test_nominee.txt"


def _candidate() -> dict:
    spec = importlib.util.spec_from_file_location("rsf_tlearner_permimp_arena._candidate_runtime", TRAIN_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load train.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return dict(module.CANDIDATE)


def cumulative_elapsed(rows: list[dict]) -> float:
    return float(sum(float(row.get("wall_clock_seconds", 0.0)) for row in rows if row.get("event") == "finished"))


def next_timeout(rows: list[dict]) -> int:
    remaining = int(prepare.BUDGET["cumulative_wall_time_seconds"] - cumulative_elapsed(rows))
    return max(0, min(int(prepare.BUDGET["max_run_wall_time_seconds"]), remaining))


def _experiments(rows: list[dict]) -> set[int]:
    return {int(row["experiment"]) for row in rows if "experiment" in row}


def _validate_result(result: dict, artifact_dir: Path, smoke: bool) -> None:
    required = {"schema_version", "candidate", "reward", "eligible", "gates", "generalization", "development_oof", "repeat_results", "gene_selection", "permutation_importance", "diagnostics", "data", "environment", "budget", "artifacts"}
    if required - set(result):
        raise RuntimeError(f"Candidate result missing {sorted(required - set(result))}")
    if result["data"].get("test_used") or result["data"].get("confirmatory_validation_used"):
        raise RuntimeError("Candidate reports forbidden validation/test access")
    if bool(result["budget"].get("smoke_reductions_applied")) != smoke:
        raise RuntimeError("Smoke marker mismatch")
    for key in ("max_experiments", "cumulative_wall_time_seconds", "max_run_wall_time_seconds", "outer_folds", "outer_repeats", "inner_pfi_folds", "forest_seeds", "patient_bootstraps", "run20_reward_to_beat", "clinical_tlearner"):
        if result["budget"].get(key) != prepare.BUDGET.get(key):
            raise RuntimeError(f"Candidate changed fixed budget field {key}")
    if result["environment"] != integrity.current_environment():
        raise RuntimeError("Candidate environment record is not launch environment")
    for relative in result["artifacts"].values():
        if relative.endswith("sha256"):
            continue
    required_files = {"candidate.json", "observation_pfi.csv", "act_pfi.csv", "joint_pfi.csv", "observation_pfi.png", "act_pfi.png", "joint_pfi.png", "permutation_importance.png", "selected_panels.json", "diagnostics.json", "pool_lineage.json"}
    missing = sorted(name for name in required_files if not (artifact_dir / name).exists())
    if missing:
        raise RuntimeError(f"Candidate artifacts missing {missing}")


def _stream(command: list[str], timeout: int) -> tuple[int, str]:
    process = subprocess.Popen(command, cwd=prepare.REPO_ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    try:
        output, _ = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        process.terminate()
        try:
            output, _ = process.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill(); output, _ = process.communicate()
        output += "\n[launcher] per-run/cumulative wall-time limit exceeded\n"
        print(output, end="", flush=True)
        return 124, output
    print(output, end="", flush=True)
    return int(process.returncode), output


def _load_result_pointer(path: Path) -> dict | None:
    if not path.exists() or not path.read_text(encoding="utf-8").strip():
        return None
    result_path = RUNS_DIR / path.read_text(encoding="utf-8").strip() / "result.json"
    return json.loads(result_path.read_text(encoding="utf-8")) if result_path.exists() else None


def qualifies_for_nomination(result: dict) -> bool:
    return bool(result.get("eligible")) and float(result.get("reward", float("-inf"))) > float(prepare.BUDGET["run20_reward_to_beat"])


def validate_round_contract(experiment: int, candidate: prepare.CandidateSpec, parent_dir: Path | None) -> None:
    if experiment == 1 and candidate.parent_run_id is not None:
        raise RuntimeError("The unique root run must have parent_run_id null")
    if experiment > 1 and candidate.parent_run_id is None:
        raise RuntimeError("Every run after the root must name a completed parent")
    if parent_dir is not None and not (parent_dir / "result.json").exists():
        raise RuntimeError("Named parent run is not completed")


def run(*, smoke: bool = False) -> dict:
    integrity.verify_lock()
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    rows = integrity.read_ledger(LEDGER_PATH)
    raw_candidate = _candidate()
    candidate = prepare.validate_candidate(raw_candidate)
    parent_dir = RUNS_DIR / candidate.parent_run_id if candidate.parent_run_id else None
    train_sha = integrity.sha256_file(TRAIN_PATH)
    if smoke:
        if candidate.parent_run_id is not None:
            raise RuntimeError("Smoke requires a root candidate")
        experiment = 0
        timeout = int(prepare.BUDGET["max_run_wall_time_seconds"])
    else:
        experiments = _experiments(rows)
        if len(experiments) >= int(prepare.BUDGET["max_experiments"]):
            raise RuntimeError("The 20-attempt budget is exhausted")
        timeout = next_timeout(rows)
        if timeout <= 0:
            raise RuntimeError("The cumulative six-hour budget is exhausted")
        if train_sha in {row.get("train_sha256") for row in rows if row.get("event") == "started"}:
            raise RuntimeError("This exact train.py already consumed an attempt")
        experiment = max(experiments, default=0) + 1
        validate_round_contract(experiment, candidate, parent_dir)
        integrity.append_ledger(LEDGER_PATH, {"event": "started", "experiment": experiment, "started_at_utc": datetime.now(timezone.utc).isoformat(), "train_sha256": train_sha, "parent_run_id": candidate.parent_run_id, "environment": integrity.current_environment()})
    temp_dir = Path(tempfile.mkdtemp(prefix="candidate_", dir=RUNS_DIR))
    result_path = temp_dir / "candidate_result.json"
    command = [sys.executable, "-m", "rsf_tlearner_permimp_arena.train", "--artifact-dir", str(temp_dir), "--result-path", str(result_path)]
    if parent_dir is not None: command.extend(["--parent-run-dir", str(parent_dir)])
    if smoke: command.append("--smoke")
    started = time.monotonic(); return_code, output = _stream(command, timeout); elapsed = time.monotonic() - started
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    if return_code != 0:
        if not smoke:
            failure_dir = RUNS_DIR / f"failed_{experiment:03d}_{stamp}"
            temp_dir.replace(failure_dir); (failure_dir / "stdout.log").write_text(output, encoding="utf-8"); (failure_dir / "train_snapshot.py").write_bytes(TRAIN_PATH.read_bytes())
            integrity.append_ledger(LEDGER_PATH, {"event": "finished", "experiment": experiment, "status": "failed", "return_code": return_code, "wall_clock_seconds": elapsed, "train_sha256": train_sha, "artifact_dir": failure_dir.name})
        else:
            shutil.rmtree(temp_dir)
        raise RuntimeError(f"Candidate failed with return code {return_code}")
    result = json.loads(result_path.read_text(encoding="utf-8")); _validate_result(result, temp_dir, smoke)
    if smoke:
        smoke_dir = RUNS_DIR / f"smoke_{stamp}"; temp_dir.replace(smoke_dir)
        (smoke_dir / "stdout.log").write_text(output, encoding="utf-8"); (smoke_dir / "train_snapshot.py").write_bytes(TRAIN_PATH.read_bytes())
        print(f"\n[launcher] smoke saved {smoke_dir}", flush=True)
        return result
    run_id = f"run_{experiment:03d}_{stamp}"
    lineage_path = temp_dir / "pool_lineage.json"
    lineage = json.loads(lineage_path.read_text(encoding="utf-8")); lineage["run_id"] = run_id
    lineage_path.write_text(json.dumps(lineage, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    result["run_id"] = run_id; result["wall_clock_seconds"] = elapsed
    result["artifacts"]["pool_lineage_sha256"] = integrity.sha256_file(lineage_path)
    (temp_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (temp_dir / "stdout.log").write_text(output, encoding="utf-8"); (temp_dir / "train_snapshot.py").write_bytes(TRAIN_PATH.read_bytes())
    run_dir = RUNS_DIR / run_id; temp_dir.replace(run_dir)
    integrity.append_ledger(LEDGER_PATH, {"event": "finished", "experiment": experiment, "status": "completed", "run_id": run_id, "reward": float(result["reward"]), "score_before_eligibility": float(result["generalization"]["score_before_eligibility"]), "eligible": bool(result["eligible"]), "beats_run20_reward": bool(result["beats_run20_reward"]), "wall_clock_seconds": elapsed, "train_sha256": train_sha})
    prior_best = _load_result_pointer(BEST_PATH)
    if result["eligible"] and (prior_best is None or float(result["reward"]) > float(prior_best["reward"])):
        BEST_PATH.write_text(run_id + "\n", encoding="utf-8")
    prior_diagnostic = _load_result_pointer(DIAGNOSTIC_PATH)
    if prior_diagnostic is None or float(result["generalization"]["score_before_eligibility"]) > float(prior_diagnostic["generalization"]["score_before_eligibility"]):
        DIAGNOSTIC_PATH.write_text(run_id + "\n", encoding="utf-8")
    if qualifies_for_nomination(result):
        prior_nominee = _load_result_pointer(NOMINEE_PATH)
        if prior_nominee is None or float(result["reward"]) > float(prior_nominee["reward"]):
            NOMINEE_PATH.write_text(run_id + "\n", encoding="utf-8")
    if not NOMINEE_PATH.read_text(encoding="utf-8").strip():
        print("NO_NEW_TEST_NOMINEE", flush=True)
    print(f"[launcher] saved {run_dir}", flush=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--smoke", action="store_true")
    run(smoke=parser.parse_args().smoke)


if __name__ == "__main__":
    main()
