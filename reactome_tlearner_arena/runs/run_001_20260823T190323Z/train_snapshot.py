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
    "name": "tlearner_raw4_drgene_locked_geometry",
    "selector": "dr_gene",
    "n_genes": 4,
    "representation": "raw",
    "module_count": 4,
    "benefit_threshold_months": 0.0,
    "tlearner": {
        # Observation arm is geometrically identical to the locked clinical
        # learner, so the four raw genes are the only difference in this arm.
        "observation": {
            "n_estimators": 1000,
            "max_depth": 9,
            "min_samples_leaf": 8,
            "min_samples_split": 16,
            "max_features": 1.0,
        },
        # ACT arm carries ~114 patients per fitting partition and is the
        # noisiest term in RMST1 - RMST0, so it is regularized harder than
        # the locked clinical ACT arm (depth 7, leaf 12, split 24).
        "act": {
            "n_estimators": 1000,
            "max_depth": 6,
            "min_samples_leaf": 16,
            "min_samples_split": 32,
            "max_features": 1.0,
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

