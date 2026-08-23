# Reactome causal survival forest arena

This arena is a clean methodological successor to `../reactome_rsf_arena_v2`.
V2 remains an immutable completed search. This arena replaces its survival
forest S-learner with R `grf::causal_survival_forest`, while retaining its
pooled-development, repeated-cross-fit, locked-policy-evaluation design.

## Model and estimand

For the 60-month horizon the causal forest estimates

```text
tau60(x) = E[min(T(1), 60) - min(T(0), 60) | X=x]
d(x)     = 1{tau60(x) > 0}
```

`X` contains pretreatment clinical covariates and, for C+G, fold-selected
Reactome genes or fit-only modules. ACT is passed separately as `W`; it is
never a feature and no gene-by-ACT products are constructed.

The clinical and C+G causal forests receive identical candidate geometry,
all-feature split sampling, honest estimation, and seeds 42, 43, and 44. Thus
the agent cannot weaken only the clinical comparator. A separate matched GRF
survival forest at locked seed 42 supplies the secondary held-out prognostic
risk ranking; it is not needlessly repeated across causal-forest seeds.

## Independent objective

The causal forests learn policies only on each outer fitting partition. A
locked Python nuisance pipeline independently produces assessment-fold
IPCW-AIPW arm scores shared by the two policies. For each repeated four-fold
outer cross-fit:

```text
A60(d)      = V60(d) - V60(1-d)
increment_r = A60(C+G) - A60(clinical)
reward_if_eligible = min_r(selection-adjusted LCB_r)
                     - range_r(increment_r)
```

An ineligible candidate receives the `-1000000` sentinel, but a continuous
diagnostic leader remains available for scientific iteration. Eligibility
requires repeat-consistent genomic value, C-index noninferiority, seed
agreement and CATE correlation, nonzero genomic variable importance, gene
stability, nontrivial effects, and treatment overlap.

The former train and validation CSVs are pooled as adaptive development data.
Their source label is used to stratify folds and to report source-specific OOF
increments, but neither source is confirmatory and source metrics do not enter
the reward. Test data is outside the arena contract and is run exactly once by
a human after a winner is frozen, with no experimentation afterward.

## Commands

```bash
/Users/owensun/miniconda3/bin/python -m reactome_csf_arena.run --smoke
/Users/owensun/miniconda3/bin/python -m reactome_csf_arena.run
```

Smoke uses one seed, one repeat, two folds, four selector subsamples, 40 trees
per top-level forest, and 100 patient bootstraps. It verifies infrastructure
without consuming the 20-experiment ledger. Full runs permit at most 1,000
trees in each causal or prognostic forest and have a 30-minute wall limit.
Smoke is SHA-locked to the shipped baseline `train.py`; after the agent edits
the candidate it cannot use smoke for unbudgeted scientific feedback.

The budget explicitly counts GRF's internal work. At the maximum setting a
full run fits 48 causal forests, their 96 internal 250-tree survival/censoring
nuisance forests, and 16 seed-42 prognostic forests: at most 88,000 constructed
trees total, with no individual forest above 1,000. Supplying the locked
cross-fitted `W.hat` prevents GRF from adding an internal propensity forest.
The survival nuisances use a fixed 121-point grid through month 60.

## Interpretation

This remains an observational analysis. Causal interpretation requires
adequate measured confounding control, positivity, consistency, and
conditionally independent censoring. Repeated row-level cross-fitting measures
internal stability, not transport across studies or expression platforms.
