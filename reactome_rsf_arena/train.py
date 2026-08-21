"""Agent-editable candidate for the Reactome causal-RMST RSF arena.

Only this Python file is edited during an autonomous search.  Frozen loading,
folds, metrics, budgets, and artifact handling live in ``prepare.py`` and
``run.py``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from . import prepare
except ImportError:
    import prepare


# Stability-selection controls.  These govern only how the train-only DR gene
# ranking is aggregated; they never touch the estimand, gates, or budget.
STABILITY_SUBSAMPLES = 40
STABILITY_TOP_K = 300
STABILITY_SEED = 20260823
# The DR pseudo-outcome is heavy tailed: the 0.05 propensity clip and the 0.05
# censoring-survival floor each admit weights up to 20, so a few patients would
# otherwise decide every gene score.  Winsorizing bounds each patient's
# influence without dropping anyone.
WINSOR_PERCENT = 5.0


CANDIDATE = {
    "name": "comppairs_winsor5_dr_gene_4_zerothresh",
    "selector": "dr_gene",
    "n_genes": 4,
    "benefit_threshold_months": 0.00,
    "rsf": {
        "n_estimators": 600,
        "max_depth": 9,
        "min_samples_leaf": 8,
        "min_samples_split": 16,
        "max_features": 1.00,
    },
}


def _stratified_halves(
    fit: pd.DataFrame, rng: np.random.Generator
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split the fitting partition into two disjoint (OS_STATUS, ACT) halves.

    Stratifying keeps both halves evaluable by the locked propensity, censoring,
    and outcome nuisance models, and keeps the two DR rankings comparable.
    """
    keys = np.asarray(list(zip(
        fit["OS_STATUS"].to_numpy(int), fit[prepare.TREATMENT].to_numpy(int)
    )))
    left: list[int] = []
    right: list[int] = []
    for stratum in np.unique(keys, axis=0):
        members = np.flatnonzero((keys == stratum).all(axis=1))
        shuffled = rng.permutation(members)
        cut = len(shuffled) // 2
        left.extend(shuffled[:cut].tolist())
        right.extend(shuffled[cut:].tolist())
    return (
        fit.iloc[sorted(left)].reset_index(drop=True),
        fit.iloc[sorted(right)].reset_index(drop=True),
    )


def _robust_gamma(gamma: np.ndarray) -> np.ndarray:
    """Winsorize the DR benefit pseudo-outcome inside the scoring partition."""
    low, high = np.percentile(gamma, [WINSOR_PERCENT, 100.0 - WINSOR_PERCENT])
    return np.clip(gamma, low, high)


def _dr_ranking(
    part: pd.DataFrame, available: list[str], inner_folds: int
) -> np.ndarray:
    """Locked train-only DR benefit scores for ``available``, as an array."""
    gamma = _robust_gamma(prepare.cross_fitted_benefit_pseudo_outcome(part, inner_folds))
    scores = prepare._gene_effect_scores(part, available, gamma)
    return np.asarray([scores[gene] for gene in available], dtype=float)


def stability_select_genes(
    fit: pd.DataFrame,
    pathways: dict[str, tuple[str, ...]],
    genes: list[str],
    spec: "prepare.CandidateSpec",
    *,
    smoke: bool = False,
) -> list[str]:
    """Rank genes by simultaneous selection across complementary half-samples.

    Each draw splits the fitting partition into two disjoint stratified halves
    and recomputes the locked cross-fitted treatment-benefit pseudo-outcome and
    clinical-adjusted gene scores independently in each.  A gene scores only
    when it reaches the top of *both* rankings, so evidence that rests on a few
    influential patients cannot promote it.  Everything is computed inside the
    fitting partition.
    """
    available = genes[:160] if smoke else genes
    draws = 4 if smoke else STABILITY_SUBSAMPLES
    inner_folds = 2 if smoke else int(prepare.BUDGET["inner_folds"])
    top_k = max(1, min(STABILITY_TOP_K, len(available) // 2))
    rng = np.random.default_rng(STABILITY_SEED)

    index = {gene: position for position, gene in enumerate(available)}
    both_counts = np.zeros(len(available), dtype=float)
    rank_total = np.zeros(len(available), dtype=float)

    for _ in range(draws):
        left, right = _stratified_halves(fit, rng)
        chosen = []
        for part in (left, right):
            values = _dr_ranking(part, available, inner_folds)
            order = np.argsort(-values, kind="stable")
            rank_total[order] += np.arange(len(available), dtype=float)
            flags = np.zeros(len(available), dtype=bool)
            flags[order[:top_k]] = True
            chosen.append(flags)
        both_counts += (chosen[0] & chosen[1]).astype(float)

    mean_rank = rank_total / float(2 * draws)
    ordered = sorted(
        available,
        key=lambda gene: (-both_counts[index[gene]], mean_rank[index[gene]], gene),
    )
    return ordered[: spec.n_genes]


def run(*, smoke: bool = False) -> dict:
    """Evaluate exactly one candidate under the frozen train-only protocol."""
    return prepare.evaluate_candidate(
        CANDIDATE,
        selector=lambda fit, pathway_map, available, candidate: stability_select_genes(
            fit, pathway_map, available, candidate, smoke=smoke
        ),
        smoke=smoke,
    )


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
