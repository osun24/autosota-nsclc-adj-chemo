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
STABILITY_SUBSAMPLES = 60
STABILITY_FRACTION = 0.70
STABILITY_TOP_K = 60
STABILITY_TOP_PATHWAYS = 50
STABILITY_SEED = 20260823


CANDIDATE = {
    "name": "reactome_stability60_dr_gene_4_greedy",
    "selector": "dr_gene",
    "n_genes": 4,
    "benefit_threshold_months": 0.25,
    "rsf": {
        "n_estimators": 600,
        "max_depth": 9,
        "min_samples_leaf": 8,
        "min_samples_split": 16,
        "max_features": 1.00,
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

    index = {gene: position for position, gene in enumerate(available)}
    # Reactome membership restricted to the genes actually offered to us.
    members = {
        name: np.asarray([index[gene] for gene in genes_in if gene in index], dtype=int)
        for name, genes_in in pathways.items()
    }
    members = {name: idx for name, idx in members.items() if idx.size}
    pathway_names = sorted(members)
    top_pathways = min(STABILITY_TOP_PATHWAYS, len(pathway_names))

    gene_counts = np.zeros(len(available), dtype=float)
    pathway_counts = np.zeros(len(pathway_names), dtype=float)
    rank_total = np.zeros(len(available), dtype=float)

    for _ in range(draws):
        subsample = _stratified_subsample(fit, STABILITY_FRACTION, rng)
        gamma = prepare.cross_fitted_benefit_pseudo_outcome(subsample, inner_folds)
        scores = prepare._gene_effect_scores(subsample, available, gamma)
        values = np.asarray([scores[gene] for gene in available], dtype=float)

        order = np.argsort(-values, kind="stable")
        rank_total[order] += np.arange(len(available), dtype=float)
        gene_counts[order[:top_k]] += 1.0

        # Pathway evidence: mean of the pathway's top-3 member scores, matching
        # the locked dr_pathway aggregation but recomputed per subsample.
        pathway_scores = np.asarray([
            float(np.mean(np.sort(values[members[name]])[-3:]))
            for name in pathway_names
        ], dtype=float)
        for position in np.argsort(-pathway_scores, kind="stable")[:top_pathways]:
            pathway_counts[position] += 1.0

    gene_frequency = gene_counts / float(draws)
    pathway_frequency = np.zeros(len(available), dtype=float)
    for position, name in enumerate(pathway_names):
        share = pathway_counts[position] / float(draws)
        if share > 0:
            idx = members[name]
            pathway_frequency[idx] = np.maximum(pathway_frequency[idx], share)

    combined = gene_frequency + pathway_frequency
    mean_rank = rank_total / float(draws)
    ordered = sorted(
        available,
        key=lambda gene: (-combined[index[gene]], mean_rank[index[gene]], gene),
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
