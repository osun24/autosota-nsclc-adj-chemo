"""Editable candidate for the progressive RSF T-learner PFI arena."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import prepare


CANDIDATE = {
    "name": "root_universe_screen_t100d4_pool1024",
    "parent_run_id": None,
    "max_panel_genes": {"observation": 32, "act": 32},
    "next_pool_size": {"observation": 1024, "act": 480},
    "benefit_threshold_months": 0.0,
    "screening_tlearner": {
        "observation": {
            "n_estimators": 100,
            "max_depth": 4,
            "min_samples_leaf": 16,
            "min_samples_split": 32,
            "max_features": "sqrt",
        },
        "act": {
            "n_estimators": 100,
            "max_depth": 4,
            "min_samples_leaf": 8,
            "min_samples_split": 16,
            "max_features": "sqrt",
        },
    },
    "tlearner": {
        "observation": {
            "n_estimators": 1000,
            "max_depth": 9,
            "min_samples_leaf": 8,
            "min_samples_split": 16,
            "max_features": 0.35,
        },
        "act": {
            "n_estimators": 1000,
            "max_depth": 7,
            "min_samples_leaf": 12,
            "min_samples_split": 24,
            "max_features": 0.35,
        },
    },
}


def run(*, artifact_dir: Path, parent_run_dir: Path | None = None, smoke: bool = False) -> dict:
    return prepare.evaluate_candidate(
        CANDIDATE, artifact_dir=artifact_dir,
        parent_run_dir=parent_run_dir, smoke=smoke,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--result-path", type=Path, required=True)
    parser.add_argument("--parent-run-dir", type=Path)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    result = run(artifact_dir=args.artifact_dir, parent_run_dir=args.parent_run_dir, smoke=args.smoke)
    args.result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
