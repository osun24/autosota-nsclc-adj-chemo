# Reactome causal-RMST RSF arena

This is an autoresearch-style arena for developing a Reactome-informed random
survival forest treatment policy. It keeps the small, reviewable contract used
by `rsf_arena`:

- `prepare.py`: frozen train/validation loading, nested feature selection,
  nuisance estimation, RSF fitting, counterfactual RMST prediction, and scoring;
- `train.py`: the only agent-editable Python file;
- `program.md`: the autonomous-loop protocol;
- `red_lines.md`: scientific and data-access restrictions;
- `run.py`: fixed wall-time/budget launcher and append-only ledger;
- `log.md`: human-readable experiment record;
- `runs/`: immutable per-experiment results and `train.py` snapshots.

The agent can use `affyfRMATrain.csv` and `affyfRMAValidation.csv`. Test data is
outside the arena contract and is evaluated once by a human after the search.
The pinned MSigDB Reactome collection and checksum are reused from
`../reactome_rsf/data/`.

## Objective

For policy `d`, define its 60-month counterfactual alignment value as

```text
A60(d) = V60(d) - V60(1-d)
```

where `V60(d)` is the IPCW-AIPW estimate of restricted survival time if the
target population followed `d`. Each candidate compares the clinical-plus-gene
policy with a matched clinical-only policy:

```text
increment = A60(C+G) - A60(C)
reward = min(train_OOF_selection_LCB, validation_selection_LCB)
         - abs(train_OOF_increment - validation_increment)
```

The LCB is multiplicity-adjusted for the fixed maximum of 20 experiments. A
candidate receives `-1000000` unless all policy-value, C-index noninferiority,
overlap, gene-stability, effect-size, and seed-stability gates pass on both
train OOF and validation.

RSF recommendations are based on integrated counterfactual survival curves:

```text
recommend ACT if predicted_RMST60(ACT) - predicted_RMST60(OBS) > threshold
```

They never use the time-independent RSF mortality score for treatment choice.
Harrell C-index remains secondary.

## Commands

Use the repository's analysis environment:

```bash
/Users/owensun/miniconda3/bin/python -m reactome_rsf_arena.run --smoke
/Users/owensun/miniconda3/bin/python -m reactome_rsf_arena.run
```

The smoke command uses reduced folds, genes, trees, seeds, and bootstraps and
does not consume the 20-experiment budget. Its metrics are infrastructure-only.
Every full run consumes one ledger slot, including failures and timeouts. An
identical `train.py` cannot consume another full slot.

The 25-minute wall-clock limit includes all feature selection, fitting, and
evaluation. A candidate may use 100 to 1,000 trees per forest. The fixed
scientific budget is in `budget.json`; changing it invalidates the arena lock.

The genomic RSF receives gene main effects only. It contains no explicit
gene-by-ACT product columns; nonlinear treatment-effect modification is learned
through tree splits involving the ACT indicator and genomic features.

Every result includes locked ACT-use diagnostics for both the clinical and
clinical-plus-genomic forests on train OOF and validation: the fraction of
trees containing an ACT split; the fraction of counterfactual patient-tree
paths traversing ACT (averaged over the ACT=0 and ACT=1 paths); and the
distribution across patients of the fraction of trees whose ACT=0 and ACT=1
copies reach different terminal nodes. The autonomous agent must copy the C+G
values into every completed `log.md` iteration. These are transparency
diagnostics, not eligibility gates.

## Interpretation

This is observational treatment-policy development. IPCW-AIPW estimates still
require consistency, adequate measured-confounder adjustment, positivity, and
conditionally independent censoring. A test result supports transport to the
held-out random split, not necessarily to a different study or expression
platform. Study-level external validation remains necessary for a clinical
treatment recommendation claim.
