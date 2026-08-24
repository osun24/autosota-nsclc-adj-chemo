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
    "name": "tlearner_wide_anchor_module16",
    # Recorded enum value.  The effective selection is the custom
    # pathway-anchored, sign-coherent selector below, injected through the
    # locked evaluator's own ``selector`` hook.  It is a global gene ranking
    # restricted to a Reactome-pathway pool, so ``dr_gene`` is the closest
    # of the three allowed labels.
    "selector": "dr_gene",
    "n_genes": 16,
    "representation": "module",
    "module_count": 1,
    "benefit_threshold_months": 0.0,
    # Both arms match the locked clinical T-learner exactly, so the single
    # genomic module column is the only difference between the two panels.
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
# The candidate pool is the union of the top ``PER_PATHWAY_TOP`` members of
# each of the top ``POOL_PATHWAYS`` pathways.  Experiment 2 pooled whole
# pathways, so one pathway flip replaced the entire pool and two of eight
# folds jumped to an unrelated family.  Bounding each pathway's contribution
# keeps the pool small while spreading anchor risk across many pathways.
#
# Experiment 3 revealed that the mechanism is anchor *spreading*, not gene
# identity: Jaccard stayed flat at 0.14 while the repeat range fell 61%,
# because a fold that switches its top anchor still shares the remaining
# anchors with the other folds.  Experiment 5 pushes that mechanism further -
# eight anchors contributing three genes each, so one anchor flip perturbs
# 3 of 16 panel slots instead of 8 of 16.
POOL_PATHWAYS = 8
PER_PATHWAY_TOP = 3


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


def select_pathway_benefit_module(
    fit,
    pathways: dict[str, tuple[str, ...]],
    genes: list[str],
    spec,
    *,
    smoke: bool = False,
) -> list[str]:
    """Pathway-anchored, sign-coherent fit-only gene selector.

    1. Score every available gene by its signed partial correlation with the
       cross-fitted DR benefit pseudo-outcome (fitting partition only).
    2. Score every sufficiently large Reactome pathway by the *mean signed*
       score of its members.  A directionally coherent pathway is a real
       aggregate signal, and a mean over >= 25 members is far lower variance
       than any individual gene score.
    3. Pool the top ``PER_PATHWAY_TOP`` members of each of the top
       ``POOL_PATHWAYS`` pathways, bounding the pool at 40 genes.
    4. Return the ``n_genes`` pool members with the most positive scores, so
       the single averaged module is sign-coherent rather than self-cancelling.
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

    pool: list[str] = []
    seen: set[str] = set()
    for _, _, present in ranked_pathways[:POOL_PATHWAYS]:
        best = sorted(present, key=lambda gene: (-scores[gene], gene))[:PER_PATHWAY_TOP]
        for gene in best:
            if gene not in seen:
                seen.add(gene)
                pool.append(gene)

    ordered = sorted(pool, key=lambda gene: (-scores[gene], gene))
    selected = [gene for gene in ordered if scores[gene] > 0][: spec.n_genes]

    if len(selected) < spec.n_genes:
        # Pool exhausted of positive-direction genes: fall back to the global
        # ranking, still taking the most positive scores first.
        chosen = set(selected)
        for gene in sorted(available, key=lambda item: (-scores[item], item)):
            if len(selected) >= spec.n_genes:
                break
            if gene not in chosen:
                chosen.add(gene)
                selected.append(gene)
    return selected[: spec.n_genes]


def run(*, smoke: bool = False) -> dict:
    def selector(fit, pathways, genes, spec):
        return select_pathway_benefit_module(fit, pathways, genes, spec, smoke=smoke)

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
