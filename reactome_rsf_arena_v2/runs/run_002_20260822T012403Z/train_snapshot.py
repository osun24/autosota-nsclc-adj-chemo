"""Agent-editable candidate for the Reactome causal-RMST RSF v2 arena.

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
from scipy import sparse

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
# Near-floor probes have a tiny residual norm, so a few samples that recur in
# every subsample can win the score and look "stable" while carrying no signal.
# Only genes the array actually measures are eligible for selection.
GENE_IQR_PERCENTILE = 50.0
# Pathway means are the stable unit: averaging over members shrinks the
# sampling variance that makes any single gene's DR score irreproducible.
MIN_PATHWAY_MEMBERS = 12
TOP_PATHWAYS = 20
# Spread the block over several pathways: folds then need only share one
# pathway anywhere in a short list rather than agree on a single top choice.
MAX_GENES_PER_PATHWAY = 4


CANDIDATE = {
    "name": "v2_pathway8_signed_modules",
    "selector": "dr_gene",
    "n_genes": 8,
    "representation": "module",
    "module_count": 2,
    "benefit_threshold_months": 0.00,
    "rsf": {
        "n_estimators": 1000,
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


def _detectable_genes(
    fit: pd.DataFrame, available: list[str]
) -> tuple[list[str], dict[str, float]]:
    """Keep genes in the upper half of within-partition expression spread.

    Also returns each kept gene's spread.  Expression spread is a covariate
    summary over hundreds of patients, so it is near-identical across fitting
    partitions -- which is what makes it a reproducible way to order the members
    of a pathway the DR signal has already selected.
    """
    values = fit[available].to_numpy(dtype=float)
    spread = np.nanpercentile(values, 75, axis=0) - np.nanpercentile(values, 25, axis=0)
    cutoff = float(np.nanpercentile(spread, GENE_IQR_PERCENTILE))
    keep = [gene for gene, width in zip(available, spread) if width >= cutoff]
    widths = {gene: float(width) for gene, width in zip(available, spread)}
    return (keep or list(available)), widths


def _robust_gamma(gamma: np.ndarray) -> np.ndarray:
    """Winsorize the DR benefit pseudo-outcome inside the scoring partition."""
    low, high = np.percentile(gamma, [WINSOR_PERCENT, 100.0 - WINSOR_PERCENT])
    return np.clip(gamma, low, high)


def _signed_gene_scores(
    frame: pd.DataFrame, genes: list[str], gamma: np.ndarray
) -> np.ndarray:
    """Sign-preserving replica of ``prepare._gene_effect_scores``.

    The locked helper returns the *absolute* clinical-adjusted partial
    correlation between each gene and the benefit pseudo-outcome, which is the
    right quantity for ranking but discards the direction of the association.
    A module feature is a mean of standardized members, so members whose
    benefit-association runs in opposite directions cancel inside it.  This
    returns the same statistic with its sign, computed with the identical
    design, ridge residualization, chunking, and NaN rule, so ``np.abs`` of the
    result reproduces the locked ranking exactly.  Everything is computed
    inside the fitting partition.
    """
    clinical = frame[prepare.NUISANCE_COLUMNS].to_numpy(dtype=float)
    clinical = prepare.StandardScaler().fit_transform(clinical)
    design = np.column_stack([np.ones(len(frame)), clinical])
    gamma_residual = gamma - prepare.Ridge(alpha=1.0).fit(clinical, gamma).predict(clinical)
    gamma_norm = max(float(np.linalg.norm(gamma_residual)), 1e-12)
    signed = np.empty(len(genes), dtype=float)
    chunk_size = 512
    for start in range(0, len(genes), chunk_size):
        chunk = genes[start:start + chunk_size]
        values = frame[chunk].to_numpy(dtype=float)
        medians = np.nanmedian(values, axis=0)
        bad = ~np.isfinite(values)
        if bad.any():
            values[bad] = np.take(medians, np.where(bad)[1])
        coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
        residual = values - design @ coefficients
        denominators = np.linalg.norm(residual, axis=0) * gamma_norm
        numerators = gamma_residual @ residual
        signed[start:start + len(chunk)] = np.divide(
            numerators, denominators,
            out=np.zeros_like(numerators), where=denominators > 0,
        )
    return signed


def _dr_ranking(
    part: pd.DataFrame, available: list[str], inner_folds: int
) -> np.ndarray:
    """Locked train-only DR benefit scores for ``available``, signed.

    ``np.abs`` of this array is exactly the locked ranking statistic; the sign
    is carried alongside so module composition can be made coherent.
    """
    gamma = _robust_gamma(prepare.cross_fitted_benefit_pseudo_outcome(part, inner_folds))
    return _signed_gene_scores(part, available, gamma)


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
    offered = genes[:160] if smoke else genes
    available, gene_spread = _detectable_genes(fit, offered)
    draws = 4 if smoke else STABILITY_SUBSAMPLES
    inner_folds = 2 if smoke else int(prepare.BUDGET["inner_folds"])
    top_pathways = min(TOP_PATHWAYS, max(1, len(pathways) // 4))
    rng = np.random.default_rng(STABILITY_SEED)

    index = {gene: position for position, gene in enumerate(available)}
    # Sparse pathway-membership matrix over the detectable gene pool, so a
    # pathway's mean member score is one matrix-vector product per half-sample.
    names, rows, cols = [], [], []
    for name in sorted(pathways):
        positions = [index[gene] for gene in pathways[name] if gene in index]
        if len(positions) < MIN_PATHWAY_MEMBERS:
            continue
        rows.extend([len(names)] * len(positions))
        cols.extend(positions)
        names.append(name)
    if not names:
        return sorted(available, key=lambda gene: gene)[: spec.n_genes]
    membership = sparse.csr_matrix(
        (np.ones(len(rows), dtype=float), (rows, cols)),
        shape=(len(names), len(available)),
    )
    sizes = np.asarray(membership.sum(axis=1)).ravel()

    pathway_counts = np.zeros(len(names), dtype=float)
    score_total = np.zeros(len(available), dtype=float)
    signed_total = np.zeros(len(available), dtype=float)

    for _ in range(draws):
        left, right = _stratified_halves(fit, rng)
        flags = []
        for part in (left, right):
            signed = _dr_ranking(part, available, inner_folds)
            signed_total += signed
            values = np.abs(signed)
            score_total += values
            pathway_scores = (membership @ values) / sizes
            chosen = np.zeros(len(names), dtype=bool)
            chosen[np.argsort(-pathway_scores, kind="stable")[:top_pathways]] = True
            flags.append(chosen)
        pathway_counts += (flags[0] & flags[1]).astype(float)

    mean_score = score_total / float(2 * draws)
    mean_signed = signed_total / float(2 * draws)
    pathway_order = sorted(
        range(len(names)),
        key=lambda position: (
            -pathway_counts[position],
            -float(np.mean(mean_score[membership[position].indices])),
            names[position],
        ),
    )
    def _sign_ordered(chosen: list[str]) -> list[str]:
        """Group benefit-positive genes ahead of benefit-negative ones.

        ``FeatureTransformer`` builds modules by ``np.array_split`` over this
        list in order, so sorting by descending signed association makes each
        module an internally coherent score instead of a mean over members that
        cancel.
        """
        return sorted(chosen, key=lambda gene: (-mean_signed[index[gene]], gene))

    selected: list[str] = []
    for position in pathway_order:
        members = membership[position].indices
        taken = 0
        ordered_members = sorted(
            members, key=lambda item: (-gene_spread[available[item]], available[item])
        )
        for gene_position in ordered_members:
            if taken >= MAX_GENES_PER_PATHWAY:
                break
            gene = available[gene_position]
            if gene not in selected:
                selected.append(gene)
                taken += 1
            if len(selected) >= spec.n_genes:
                return _sign_ordered(selected[: spec.n_genes])
    for gene in sorted(available, key=lambda item: (-mean_score[index[item]], item)):
        if len(selected) >= spec.n_genes:
            break
        if gene not in selected:
            selected.append(gene)
    return _sign_ordered(selected[: spec.n_genes])


def run(*, smoke: bool = False) -> dict:
    """Evaluate one candidate under pooled repeated development cross-fitting."""
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
