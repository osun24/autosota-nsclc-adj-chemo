# Reactome causal-RMST RSF arena v2

V2 is a clean continuation of the completed v1 search. V1 remains unchanged in
`../reactome_rsf_arena`; its negative outcome and reusable findings are
summarized in `v1_findings.md`.

## What changed

- The former train and validation CSVs are pooled into one 1,034-patient
  adaptive development cohort. Validation is no longer represented as an
  untouched confirmation set.
- Each candidate is evaluated by two repeated four-fold outer cross-fits.
- Selected genes can enter the genomic RSF as one to four fold-fitted module
  scores. Each gene is imputed and standardized on the fitting partition only;
  ordered groups are averaged into module features on assessment rows.
- The clinical comparator is independently locked at 1,000 trees, depth 9,
  leaf size 8, split size 16, and all-feature split sampling. The agent cannot
  weaken it through `train.py`.
- `best_run.txt` is updated only by an eligible improvement. The continuous
  pre-eligibility leader is tracked separately in `diagnostic_leader.txt` and
  is never a test nominee merely because it leads that diagnostic.

The genomic model still contains no explicit gene-by-ACT products. A raw gene
representation is allowed for controlled comparisons, while the starting
candidate uses eight selected genes compressed into two ordered modules.

## Objective

For policy `d`, let `A60(d)=V60(d)-V60(1-d)`, where `V60` is the cross-fitted
IPCW-AIPW 60-month restricted-time value. For each outer repeat:

```text
increment_r = A60(clinical+genomic) - A60(clinical)
reward_if_eligible = min_r(selection-adjusted LCB_r)
                     - (max_r increment_r - min_r increment_r)
```

An ineligible candidate receives `-1000000`. Eligibility requires positive and
directionally consistent repeat increments, positive genomic alignment,
genomic value at least as high as the locked clinical and best constant
policies, C-index noninferiority, seed agreement, gene stability, nontrivial
benefits, and treatment-overlap diagnostics in every repeat.

The decision threshold in `train.py` applies only to the genomic policy. The
locked clinical comparator uses zero months. Harrell C-index is secondary.

## Commands

```bash
/Users/owensun/miniconda3/bin/python -m reactome_rsf_arena_v2.run --smoke
/Users/owensun/miniconda3/bin/python -m reactome_rsf_arena_v2.run
```

Smoke uses one seed, one repeat, two folds, reduced selector work, 40 trees,
and 100 bootstraps. It does not consume the 20-run ledger. Full runs have a
25-minute wall limit and at most 1,000 trees per forest.

The autonomous search may use both development CSVs. Test data is outside the
arena contract and is evaluated once by a human only after an eligible run has
been frozen.

## Interpretation

Repeated cross-fitting measures internal stability; the repeats reuse patients
and are not independent external validation. IPCW-AIPW conclusions require
adequate measured-confounder adjustment, positivity, consistency, and
conditionally independent censoring. Row-level performance does not establish
transport across source studies or expression platforms.
