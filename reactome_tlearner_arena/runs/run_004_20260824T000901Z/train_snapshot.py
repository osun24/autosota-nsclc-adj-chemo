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
    "name": "tlearner_bidirectional_pool_module16",
    # Recorded enum value.  The effective selection is the custom
    # pathway-anchored, sign-stratified selector below, injected through the
    # locked evaluator's own ``selector`` hook.  It is a global gene ranking
    # restricted to a Reactome-pathway pool, so ``dr_gene`` is the closest of
    # the three allowed labels.
    "selector": "dr_gene",
    "n_genes": 16,
    "representation": "module",
    # Two modules: one averaging the benefit-increasing genes, one averaging
    # the benefit-decreasing genes.  The locked transformer groups by rank
    # position via ``array_split``, so returning 8 positives followed by 8
    # negatives makes each module sign-pure.
    "module_count": 2,
    "benefit_threshold_months": 0.0,
    # Both arms match the locked clinical T-learner exactly, so the genomic
    # module columns are the only difference between the two panels.
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

# Minimum number of development-present members a Reactome pathway must have
# before its aggregate benefit score is trusted.  Small sets win the ranking
# on noise; averaging over >= 25 members shrinks the aggregate's variance.
MIN_PATHWAY_MEMBERS = 25
# The pool for each direction is the union of the top ``PER_PATHWAY_TOP``
# members of each of ``POOL_PATHWAYS`` anchor pathways.  Experiment 3 showed
# that spreading a fold's panel across several anchors is what suppresses
# repeat-to-repeat heterogeneity, so the quota structure is kept intact.
POOL_PATHWAYS = 5
PER_PATHWAY_TOP = 8


def _signed_benefit_scores(fit, available: list[str], *, smoke: bool) -> dict[str, float]:
    """Fit-only signed partial correlation of each gene with DR benefit.

    Mirrors the locked ``_gene_effect_scores`` residualization exactly, but
    keeps the sign that the locked version discards with ``np.abs``.
    """
    inner_folds = 2 if smoke else int(prepare.BUDGET["inner_folds"])
    gamma = prepare.cross_fitted_benefit_pseudo_outcome(fit, inner_folds)
    clinical = StandardScaler().fit_transform(fit[prepare.NUISANCE_COLUMNS].to_numpy(float))
    design = np.column_stack([np.ones(len(fit)), clinical])
    gamma_residual = gamma - Ridge(alpha=1.0).fit(clinical, gamma).predict(clinical)
    gamma_norm = max(float(np.linalg.norm(gamma_residual)), 1e-12)
    scores: dict[str, float] = {}
    chunk_size = 512
    for start in range(0, len(available), chunk_size):
        chunk = available[start:start + chunk_size]
        values = fit[chunk].to_numpy(dtype=float)
        medians = np.nanmedian(values, axis=0)
        medians = np.where(np.isfinite(medians), medians, 0.0)
        bad = ~np.isfinite(values)
        if bad.any():
            values[bad] = np.take(medians, np.where(bad)[1])
        coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
        residual = values - design @ coefficients
        denominators = np.linalg.norm(residual, axis=0) * gamma_norm
        numerators = gamma_residual @ residual
        signed = np.divide(
            numerators, denominators,
            out=np.zeros_like(numerators), where=denominators > 0,
        )
        scores.update({gene: float(value) for gene, value in zip(chunk, signed)})
    return scores


def _directional_pool(
    ranked_pathways: list[tuple[float, str, list[str]]],
    scores: dict[str, float],
    sign: int,
) -> list[str]:
    """Bounded gene pool for one benefit direction.

    ``sign`` is +1 for benefit-increasing genes (anchored on the pathways with
    the most positive mean signed score) and -1 for benefit-decreasing genes
    (anchored on the most negative).  Within each anchor only the strongest
    ``PER_PATHWAY_TOP`` members in that direction are contributed.
    """
    anchors = ranked_pathways[:POOL_PATHWAYS] if sign > 0 else ranked_pathways[-POOL_PATHWAYS:]
    pool: list[str] = []
    seen: set[str] = set()
    for _, _, present in anchors:
        best = sorted(present, key=lambda gene: (-sign * scores[gene], gene))[:PER_PATHWAY_TOP]
        for gene in best:
            if gene not in seen:
                seen.add(gene)
                pool.append(gene)
    return sorted(pool, key=lambda gene: (-sign * scores[gene], gene))


def select_bidirectional_benefit_modules(
    fit,
    pathways: dict[str, tuple[str, ...]],
    genes: list[str],
    spec,
    *,
    smoke: bool = False,
) -> list[str]:
    """Pathway-anchored, sign-stratified fit-only gene selector.

    1. Score every available gene by its signed partial correlation with the
       cross-fitted DR benefit pseudo-outcome (fitting partition only).
    2. Score every Reactome pathway with at least ``MIN_PATHWAY_MEMBERS``
       present genes by the *mean signed* score of its members.
    3. Build a bounded pool from the top ``POOL_PATHWAYS`` anchors for the
       benefit-increasing direction and from the bottom ``POOL_PATHWAYS``
       anchors for the benefit-decreasing direction.
    4. Return half the panel from each direction, positives first, so the
       locked ``array_split`` grouping yields two sign-pure modules.
    """
    available = genes[:160] if smoke else list(genes)
    available_set = set(available)
    scores = _signed_benefit_scores(fit, available, smoke=smoke)

    ranked_pathways: list[tuple[float, str, list[str]]] = []
    for name, members in pathways.items():
        present = [gene for gene in members if gene in available_set]
        if len(present) < MIN_PATHWAY_MEMBERS:
            continue
        ranked_pathways.append(
            (float(np.mean([scores[gene] for gene in present])), name, present)
        )
    ranked_pathways.sort(key=lambda item: (-item[0], item[1]))

    half = spec.n_genes // 2
    quotas = ((1, half), (-1, spec.n_genes - half))
    selected: list[str] = []
    chosen: set[str] = set()
    for sign, quota in quotas:
        taken = 0
        for gene in _directional_pool(ranked_pathways, scores, sign):
            if taken >= quota:
                break
            if gene not in chosen and sign * scores[gene] > 0:
                chosen.add(gene)
                selected.append(gene)
                taken += 1
        if taken < quota:
            # Direction's pool exhausted: fall back to the global ranking,
            # still strongest-in-direction first.
            ordered = sorted(available, key=lambda item: (-sign * scores[item], item))
            for gene in ordered:
                if taken >= quota:
                    break
                if gene not in chosen:
                    chosen.add(gene)
                    selected.append(gene)
                    taken += 1
    return selected[: spec.n_genes]


def run(*, smoke: bool = False) -> dict:
    def selector(fit, pathways, genes, spec):
        return select_bidirectional_benefit_modules(fit, pathways, genes, spec, smoke=smoke)

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
