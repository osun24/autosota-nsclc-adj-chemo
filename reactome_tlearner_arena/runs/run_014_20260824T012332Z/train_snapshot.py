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
    "name": "tlearner_modal_anchor_scored_module16",
    # Recorded enum value.  The effective selection is the custom
    # pathway-anchored stability selector below, injected through the locked
    # evaluator's own ``selector`` hook.  It is a global gene ranking
    # restricted to a Reactome-pathway pool, so ``dr_gene`` is the closest of
    # the three allowed labels.
    "selector": "dr_gene",
    "n_genes": 16,
    "representation": "module",
    # Back to one column.  Experiment 12 established that a second module
    # amplifies partition sensitivity - the repeat range went 3.69 to 5.40 -
    # and the reward takes the MINIMUM repeat LCB minus that range, so
    # strengthening the better repeat is worthless.
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
# Two modal anchors contributing eight genes each.
POOL_PATHWAYS = 2
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


def _subsample_statistics(
    gamma, clinical, values,
    pathway_members: list[tuple[str, np.ndarray]],
    n_genes: int,
    draws: int,
    seed: int,
):
    """Frequency statistics over subsamples of the fitting partition.

    Returns, per subsample-repeat: how often each pathway lands in the top
    ``POOL_PATHWAYS`` of the anchor ranking, how often each gene's signed
    score is positive, and each gene's mean signed score.  Counts are bounded
    statistics with far lower variance than the ranks this search has been
    using at every level.
    """
    n = len(gamma)
    size = max(int(round(STABILITY_SUBSAMPLE * n)), 50)
    rng = np.random.default_rng(seed)
    anchor_counts: dict[str, int] = {}
    positive_counts = np.zeros(values.shape[1], dtype=float)
    score_total = np.zeros(values.shape[1], dtype=float)
    for _ in range(draws):
        rows = rng.choice(n, size=size, replace=False)
        scores = _signed_scores_on(gamma[rows], clinical[rows], values[rows])
        score_total += scores
        positive_counts += (scores > 0).astype(float)
        ranked = sorted(
            ((float(scores[members].mean()), name) for name, members in pathway_members),
            key=lambda item: (-item[0], item[1]),
        )
        for _, name in ranked[:POOL_PATHWAYS]:
            anchor_counts[name] = anchor_counts.get(name, 0) + 1
    return anchor_counts, positive_counts, score_total / float(draws)


def select_modal_anchor_module(
    fit,
    pathways: dict[str, tuple[str, ...]],
    genes: list[str],
    spec,
    *,
    smoke: bool = False,
) -> list[str]:
    """Frequency-selected, sign-coherent, fit-only Reactome module.

    Rank is the statistic that has been unstable at every level of this
    search: unstable across folds (experiment 1), across biologically
    structured pools (experiment 6), and under 2-5 gene perturbations
    (experiment 10).  This selector removes rank from both decisions.

    1. The anchor pathways are the ones that appear most *often* in the top
       ``POOL_PATHWAYS`` across subsamples, not the ones that rank highest in
       any single fit.
    2. Within an anchor, members are ordered by **mean signed score** across
       subsamples.  Experiment 13 ordered by sign-consistency instead, which
       maximised constancy (repeat range 0.1806) but cost most of the
       increment (min +0.568).  The anchor choice is what produces the large
       fold-level jumps; the within-anchor ordering is a far smaller
       perturbation, so restoring signal strength here should recover
       magnitude while keeping the modal anchor's constancy.
    """
    available = genes[:160] if smoke else list(genes)
    index_of = {gene: position for position, gene in enumerate(available)}
    available_set = set(available)
    gamma, clinical, values = _prepare_scoring_inputs(fit, available, smoke=smoke)

    pathway_members: list[tuple[str, np.ndarray]] = []
    members_by_name: dict[str, np.ndarray] = {}
    for name, members in pathways.items():
        present = [index_of[gene] for gene in members if gene in available_set]
        if len(present) >= MIN_PATHWAY_MEMBERS:
            array = np.asarray(sorted(present), dtype=int)
            pathway_members.append((name, array))
            members_by_name[name] = array

    draws = 5 if smoke else STABILITY_DRAWS
    anchor_counts, positive_counts, mean_scores = _subsample_statistics(
        gamma, clinical, values, pathway_members, spec.n_genes, draws, STABILITY_SEED
    )

    anchors = sorted(anchor_counts, key=lambda name: (-anchor_counts[name], name))[:POOL_PATHWAYS]
    quota = spec.n_genes // max(len(anchors), 1) if anchors else spec.n_genes
    selected: list[int] = []
    chosen: set[int] = set()
    for name in anchors:
        members = members_by_name[name]
        order = sorted(
            (int(position) for position in members if mean_scores[position] > 0),
            key=lambda position: (-mean_scores[position], available[position]),
        )
        taken = 0
        for position in order:
            if taken >= quota:
                break
            if position not in chosen:
                chosen.add(position)
                selected.append(position)
                taken += 1

    if len(selected) < spec.n_genes:
        # Anchors exhausted of positive-direction genes: fall back to the
        # global frequency ranking, still positive-direction first.
        order = sorted(
            range(len(available)),
            key=lambda position: (
                -positive_counts[position], -mean_scores[position], available[position]
            ),
        )
        for position in order:
            if len(selected) >= spec.n_genes:
                break
            if position not in chosen:
                chosen.add(position)
                selected.append(position)
    return [available[position] for position in selected[: spec.n_genes]]


def run(*, smoke: bool = False) -> dict:
    def selector(fit, pathways, genes, spec):
        return select_modal_anchor_module(fit, pathways, genes, spec, smoke=smoke)

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
