"""Agent-editable candidate for the Reactome RSF T-learner arena."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

try:
    from . import prepare
except ImportError:
    import prepare


CANDIDATE = {
    "name": "tlearner_stability_selected_module16",
    # Recorded enum value.  The effective selection is the custom
    # pathway-anchored stability selector below, injected through the locked
    # evaluator's own ``selector`` hook.  It is a global gene ranking
    # restricted to a Reactome-pathway pool, so ``dr_gene`` is the closest of
    # the three allowed labels.
    "selector": "dr_gene",
    "n_genes": 16,
    "representation": "module",
    "module_count": 1,
    "benefit_threshold_months": 0.0,
    # Both arms match the locked clinical T-learner exactly, so the single
    # genomic module column is the only difference between the two panels.
    # Experiments 7 and 8 established that locked geometry is at or near the
    # optimum for the C+G panel in both directions.
    "tlearner": {
        "observation": {
            "n_estimators": 1000,
            "max_depth": 9,
            "min_samples_leaf": 8,
            "min_samples_split": 16,
            "max_features": 1.0,
        },
        "act": {
            "n_estimators": 1000,
            "max_depth": 7,
            "min_samples_leaf": 12,
            "min_samples_split": 24,
            "max_features": 1.0,
        },
    },
}

MIN_PATHWAY_MEMBERS = 25
POOL_PATHWAYS = 5
PER_PATHWAY_TOP = 8
# Stability selection.  Ten experiments agree that the binding defect is that
# the top of a noisy 8,647-gene ranking does not reproduce: experiment 1
# measured Jaccard 0.070 at n_genes=4, and experiment 10 showed that changing
# 2-5 genes in 3 of 8 folds swings the repeat range by 3.1 months.  Rather
# than re-rank once on the whole fitting partition, re-run the entire
# pathway-anchored pick on many subsamples of it and keep the genes chosen
# most often.  Selection frequency is a far lower-variance statistic than any
# single ranking, which is the point of the method.
STABILITY_DRAWS = 40
STABILITY_SUBSAMPLE = 0.80
STABILITY_SEED = 20260824


def _prepare_scoring_inputs(fit, available: list[str], *, smoke: bool):
    """Fit-only matrices reused across every stability subsample."""
    inner_folds = 2 if smoke else int(prepare.BUDGET["inner_folds"])
    gamma = prepare.cross_fitted_benefit_pseudo_outcome(fit, inner_folds)
    clinical = fit[prepare.NUISANCE_COLUMNS].to_numpy(float)
    values = fit[available].to_numpy(dtype=float)
    medians = np.nanmedian(values, axis=0)
    medians = np.where(np.isfinite(medians), medians, 0.0)
    bad = ~np.isfinite(values)
    if bad.any():
        values[bad] = np.take(medians, np.where(bad)[1])
    return gamma, clinical, values


def _signed_scores_on(gamma, clinical, values) -> np.ndarray:
    """Signed partial correlation of each gene column with DR benefit.

    Mirrors the locked ``_gene_effect_scores`` residualization exactly but
    keeps the sign that the locked version discards with ``np.abs``.
    """
    clinical = StandardScaler().fit_transform(clinical)
    design = np.column_stack([np.ones(len(clinical)), clinical])
    gamma_residual = gamma - Ridge(alpha=1.0).fit(clinical, gamma).predict(clinical)
    gamma_norm = max(float(np.linalg.norm(gamma_residual)), 1e-12)
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    residual = values - design @ coefficients
    denominators = np.linalg.norm(residual, axis=0) * gamma_norm
    numerators = gamma_residual @ residual
    return np.divide(
        numerators, denominators,
        out=np.zeros_like(numerators), where=denominators > 0,
    )


def _pathway_anchored_pick(
    scores: np.ndarray,
    index_of: dict[str, int],
    pathway_members: list[tuple[str, np.ndarray]],
    available: list[str],
    n_genes: int,
) -> list[str]:
    """One pathway-anchored, sign-coherent panel for a single score vector."""
    ranked = sorted(
        ((float(scores[members].mean()), name, members) for name, members in pathway_members),
        key=lambda item: (-item[0], item[1]),
    )
    pool: list[int] = []
    seen: set[int] = set()
    for _, _, members in ranked[:POOL_PATHWAYS]:
        order = members[np.argsort(-scores[members], kind="stable")][:PER_PATHWAY_TOP]
        for position in order:
            if int(position) not in seen:
                seen.add(int(position))
                pool.append(int(position))
    pool.sort(key=lambda position: (-scores[position], available[position]))
    picked = [position for position in pool if scores[position] > 0][:n_genes]
    if len(picked) < n_genes:
        chosen = set(picked)
        for position in np.argsort(-scores, kind="stable"):
            if len(picked) >= n_genes:
                break
            if int(position) not in chosen:
                chosen.add(int(position))
                picked.append(int(position))
    return [available[position] for position in picked[:n_genes]]


def select_stability_benefit_module(
    fit,
    pathways: dict[str, tuple[str, ...]],
    genes: list[str],
    spec,
    *,
    smoke: bool = False,
) -> list[str]:
    """Pathway-anchored stability selector, fit-only.

    The cross-fitted DR benefit pseudo-outcome is computed once on the
    fitting partition.  The whole pathway-anchored, sign-coherent pick is
    then repeated on ``STABILITY_DRAWS`` random subsamples of that partition,
    and the genes chosen most often are returned.  Nothing outside the
    fitting partition is touched and no gene is hard-coded.
    """
    available = genes[:160] if smoke else list(genes)
    index_of = {gene: position for position, gene in enumerate(available)}
    available_set = set(available)
    gamma, clinical, values = _prepare_scoring_inputs(fit, available, smoke=smoke)

    pathway_members: list[tuple[str, np.ndarray]] = []
    for name, members in pathways.items():
        present = [index_of[gene] for gene in members if gene in available_set]
        if len(present) >= MIN_PATHWAY_MEMBERS:
            pathway_members.append((name, np.asarray(sorted(present), dtype=int)))

    n = len(fit)
    size = max(int(round(STABILITY_SUBSAMPLE * n)), 50)
    draws = 5 if smoke else STABILITY_DRAWS
    rng = np.random.default_rng(STABILITY_SEED)
    counts = np.zeros(len(available), dtype=float)
    score_total = np.zeros(len(available), dtype=float)
    for _ in range(draws):
        rows = rng.choice(n, size=size, replace=False)
        scores = _signed_scores_on(gamma[rows], clinical[rows], values[rows])
        score_total += scores
        for gene in _pathway_anchored_pick(
            scores, index_of, pathway_members, available, spec.n_genes
        ):
            counts[index_of[gene]] += 1.0

    mean_scores = score_total / float(draws)
    order = sorted(
        range(len(available)),
        key=lambda position: (-counts[position], -mean_scores[position], available[position]),
    )
    return [available[position] for position in order[: spec.n_genes]]


def run(*, smoke: bool = False) -> dict:
    def selector(fit, pathways, genes, spec):
        return select_stability_benefit_module(fit, pathways, genes, spec, smoke=smoke)

    return prepare.evaluate_candidate(CANDIDATE, selector=selector, smoke=smoke)


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
