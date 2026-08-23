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
# with trees^2 x splits. Measured constants: a panel costs about
# 0.02 + 0.00145 * n_estimators seconds, and the used-feature count follows
# 8666 * (1 - exp(-trees * splits_per_tree / 8666)) with roughly seven splits
# per tree at depth four. The screening space is sized from those constants.
#
# The screening forest's own policy metric is degenerate: with treatment as one
# candidate among thousands of genes it is never chosen for a split, so every
# screening geometry returns the same all-observation policy value. Screening
# capacity is spent on the C-index permutation ranking, which is what actually
# reaches the reduced stage.
#
# Near-greedy feature selection crowds treatment out of the trees and collapses
# the policy to all-observation even at panel size, so the reduced stage is
# restricted to the sqrt-scale family.
#
# The reduced stage enqueues the controlled N = 0..32 curve using the first
# element of every reduced list, so those lists are ordered with the intended
# reference geometry first rather than smallest first.
CANDIDATE = {
    "name": "reference_600_trees_v6_relaunch",
    "screening_space": {
        "n_estimators": [150, 250],
        "max_depth": [4, 6],
        "min_samples_leaf": [16],
        "split_leaf_multiplier": [2],
        "max_features": [0.25, 0.5],
    },
    "reduced_space": {
        "n_estimators": [600, 400, 200],
        "max_depth": [6, 4, 8],
        "min_samples_leaf": [16, 8, 32],
        "split_leaf_multiplier": [2, 3],
        "max_features": ["sqrt", "log2", 0.25],
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

