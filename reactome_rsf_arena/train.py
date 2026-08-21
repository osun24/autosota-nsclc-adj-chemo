"""Agent-editable candidate for the Reactome causal-RMST RSF arena.

Only this Python file is edited during an autonomous search.  Frozen loading,
folds, metrics, budgets, and artifact handling live in ``prepare.py`` and
``run.py``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

try:
    from . import prepare
except ImportError:
    import prepare


CANDIDATE = {
    "name": "baseline_dr_pathway_12",
    "selector": "dr_pathway",
    "n_genes": 12,
    "benefit_threshold_months": 0.25,
    "rsf": {
        "n_estimators": 300,
        "max_depth": 6,
        "min_samples_leaf": 16,
        "min_samples_split": 32,
        "max_features": 0.50,
    },
}


def run(*, smoke: bool = False) -> dict:
    """Evaluate exactly one candidate under the frozen train-only protocol."""
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
