"""Agent-editable candidate for the Reactome RSF T-learner arena."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

try:
    from . import prepare
except ImportError:
    import prepare


CANDIDATE = {
    "name": "tlearner_dr8_two_modules",
    "selector": "dr_pathway",
    "n_genes": 8,
    "representation": "module",
    "module_count": 2,
    "benefit_threshold_months": 0.0,
    "tlearner": {
        "observation": {
            "n_estimators": 1000,
            "max_depth": 9,
            "min_samples_leaf": 8,
            "min_samples_split": 16,
            "max_features": 0.5,
        },
        "act": {
            "n_estimators": 1000,
            "max_depth": 7,
            "min_samples_leaf": 12,
            "min_samples_split": 24,
            "max_features": 0.5,
        },
    },
}


def run(*, smoke: bool = False) -> dict:
    return prepare.evaluate_candidate(CANDIDATE, smoke=smoke)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--result-path", type=Path)
    args = parser.parse_args()
    result = run(smoke=args.smoke)
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.result_path:
        args.result_path.write_text(payload, encoding="utf-8")
    print(payload)


if __name__ == "__main__":
    main()

