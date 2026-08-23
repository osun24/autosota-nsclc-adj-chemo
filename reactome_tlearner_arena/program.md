# Autonomous program: Reactome RSF T-learner

## Goal

Maximize the locked eligible reward from `python -m reactome_tlearner_arena.run`.
Find a small stable Reactome panel whose two-arm RSF policy improves on the
locked clinical-only T-learner in both repeated development cross-fits.

Before experiment 1, read `README.md`, `red_lines.md`, all of `log.md`,
`budget.json`, and `train.py`.

## Editable surface

- Edit only `train.py`; append prespecified hypotheses/results to `log.md`.
- Never edit locked evaluator, launcher, budget, manifest, ledger, artifacts,
  the v2 dependency, or data.
- Selectors receive only the current fitting partition, Reactome membership,
  available genes, and validated candidate specification.
- Do not hard-code genes.

## Fixed resources

- 20 full experiments; 35 minutes each.
- Two repeated four-fold outer cross-fits and three inner nuisance folds.
- Seeds 42, 43, 44; no forest above 1,000 trees.
- At most 16 genes and 34 baseline features.
- 4,000 paired patient bootstraps per repeat.
- A separately locked clinical T-learner; only C+G arm geometry is editable.

After three parameter-only experiments, the next attempt must change selector
or representation logic. Smaller panels/modules and stronger ACT-arm
regularization are the preferred first axes.

## Required result record

Record eligibility, gates, reward, diagnostic score, both repeat increments
and LCBs, repeat range, C-index, alignment, gene Jaccard, ACT recommendation
rate, arm patient/event support, seed agreement, seed benefit correlation, and
observation/ACT genomic split fractions.

Freeze only an eligible run named by `best_run.txt`. `diagnostic_leader.txt`
guides hypotheses but is not a test nominee. If no eligible run exists, report
`NO_ELIGIBLE_CANDIDATE` and preserve the test.

