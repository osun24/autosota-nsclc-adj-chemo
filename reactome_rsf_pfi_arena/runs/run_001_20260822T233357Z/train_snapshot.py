"""Agent-editable search specification for the Reactome-wide RSF PFI arena."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

try:
    from . import prepare
except ImportError:
    import prepare


# Gene symbols and patient identifiers are intentionally absent. The locked
# evaluator supplies every eligible Reactome gene to the screening model and
# derives the reduced panel solely from held-out permutation importance.
# Held-out permutation importance costs one panel per
# (fold x used feature x repeat), and the number of features a forest actually
# splits on grows with trees x splits-per-tree, so PFI wall time scales roughly
# with trees^2 x splits. The screening space is therefore capped well below the
# tree budget; capacity is spent in the reduced stage, where the feature count
# is at most nineteen clinical columns plus the selected panel.
#
# The reduced stage enqueues the controlled N = 0..32 curve using the first
# element of every reduced list, so those lists are ordered with the intended
# reference geometry first rather than smallest first.
CANDIDATE = {
    "name": "cost_bounded_screen_v1",
    "screening_space": {
        "n_estimators": [50, 100],
        "max_depth": [3, 4],
        "min_samples_leaf": [24],
        "split_leaf_multiplier": [2],
        "max_features": ["sqrt", 0.025, 0.05],
    },
    "reduced_space": {
        "n_estimators": [400, 200, 600],
        "max_depth": [6, 4, 8],
        "min_samples_leaf": [16, 8, 32],
        "split_leaf_multiplier": [2, 3],
        "max_features": ["sqrt", 0.5, 1.0],
    },
    "plot_top_k": 32,
}


def run(*, artifact_dir: str | Path, smoke: bool = False) -> dict:
    return prepare.evaluate_candidate(CANDIDATE, artifact_dir=artifact_dir, smoke=smoke)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--result-path", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    result = run(artifact_dir=args.artifact_dir, smoke=args.smoke)
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.result_path.write_text(payload, encoding="utf-8")
    print(payload)


if __name__ == "__main__":
    main()

