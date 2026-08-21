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
STABILITY_FRACTION = 0.70
STABILITY_TOP_K = 100
STABILITY_SEED = 20260823


CANDIDATE = {
    "name": "stability40_dr_gene_4",
    "selector": "dr_gene",
    "n_genes": 4,
    "benefit_threshold_months": 0.25,
    "rsf": {
        "n_estimators": 600,
        "max_depth": 9,
        "min_samples_leaf": 8,
        "min_samples_split": 16,
        "max_features": 0.25,
    },
}


def _stratified_subsample(
    fit: pd.DataFrame, fraction: float, rng: np.random.Generator
) -> pd.DataFrame:
    """Draw a subsample of the fitting partition preserving outcome/ACT strata.

    Stratifying on (OS_STATUS, ACT) keeps every subsample evaluable by the
    locked propensity, censoring, and outcome nuisance models.
    """
    keys = list(zip(fit["OS_STATUS"].to_numpy(int), fit[prepare.TREATMENT].to_numpy(int)))
    picked: list[int] = []
    for stratum in sorted(set(keys)):
        members = np.flatnonzero(np.asarray([key == stratum for key in keys]))
        take = max(int(round(fraction * len(members))), min(len(members), 8))
        picked.extend(rng.choice(members, size=take, replace=False).tolist())
    return fit.iloc[sorted(picked)].reset_index(drop=True)


def stability_select_genes(
    fit: pd.DataFrame,
    pathways: dict[str, tuple[str, ...]],
    genes: list[str],
    spec: "prepare.CandidateSpec",
    *,
    smoke: bool = False,
) -> list[str]:
    """Rank genes by how reproducibly they top a train-only DR benefit ranking.

    The whole computation lives inside ``fit``.  Each subsample recomputes the
    locked cross-fitted treatment-benefit pseudo-outcome and the locked
    clinical-adjusted gene scores, so a gene is promoted only when its
    association with treatment benefit survives resampling, not when a single
    split happens to favour it.
    """
    available = genes[:160] if smoke else genes
    draws = 4 if smoke else STABILITY_SUBSAMPLES
    inner_folds = 2 if smoke else int(prepare.BUDGET["inner_folds"])
    top_k = min(STABILITY_TOP_K, len(available))
    rng = np.random.default_rng(STABILITY_SEED)

    counts = np.zeros(len(available), dtype=float)
    rank_total = np.zeros(len(available), dtype=float)
    index = {gene: position for position, gene in enumerate(available)}

    for _ in range(draws):
        subsample = _stratified_subsample(fit, STABILITY_FRACTION, rng)
        gamma = prepare.cross_fitted_benefit_pseudo_outcome(subsample, inner_folds)
        scores = prepare._gene_effect_scores(subsample, available, gamma)
        order = sorted(available, key=lambda gene: (-scores[gene], gene))
        for rank, gene in enumerate(order):
            rank_total[index[gene]] += rank
        for gene in order[:top_k]:
            counts[index[gene]] += 1.0

    mean_rank = rank_total / float(draws)
    ordered = sorted(
        available,
        key=lambda gene: (-counts[index[gene]], mean_rank[index[gene]], gene),
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
