# Project: NSCLC adjuvant chemo HTE — Random Survival Forest arena

## Goal
Improve the `train.py` pipeline along TWO objectives jointly:
- `val_ci` (Val Harrell's C-index, higher better)
- `val_rmst_diff` (Val 5-yr RMST diff under counterfactual recommendation, higher better)
Both must improve beyond bootstrap uncertainty for a candidate to count as BETTER.

## Hard Rules
See `red_lines.md`. The test set remains sealed; `prepare.py` and `train.py` must only load Train, Validation, and `LOOCV_Genes2.csv`.

## Per-Iteration Workflow
1. Read `log.md`. Identify the last 3 iterations and their `type` (PARAM | CODE | ALGO).
2. If all 3 most recent COMPLETED iterations were PARAM, propose a CODE or ALGO idea this iteration.
3. Pick one CLEARED idea from the idea library below, or generate a new one.
4. Write an idea entry to `log.md` with type, hypothesis, expected effect, and a red-line audit BEFORE editing `train.py`.
5. `git add . && git commit -m "iter_NNN_pre: <idea_id>"` to snapshot pre-state.
6. Edit `train.py`. Run `python rsf_arena/train.py`. Wall clock budget: 25 min.
7. On completion, append the result to `log.md` using the schema below.
8. WORSE on both objectives should be reverted. MIXED may stay for at most 3 follow-up iterations.

## log.md Entry Schema
```markdown
### iter_NNN — <short_title>
- type: PARAM | CODE | ALGO
- idea_id: <slug>
- hypothesis: <one sentence>
- changed_files: <list>
- val_ci: 0.XXXX ± 0.XXXX
- val_rmst_diff: X.XX ± X.XX (months)
- n_features: <int>
- verdict: BETTER | WORSE | MIXED
- one_line_lesson: <text>
```

## Idea Library

### ALGO
- `interaction_enabled_rsf`: Reintroduce gene x ACT interactions under a fixed total feature budget. Risk: medium. Admissibility: CLEARED.
- `honest_rsf_split`: Fit feature/ranking decisions on one train split and RSF trees on the other to reduce adaptive overfit. Risk: medium. Admissibility: CLEARED.
- `stability_selection_genes`: Replace univariate ranking with train-only stability selection before RSF. Risk: low. Admissibility: CLEARED.
- `s_t_ensemble`: Average counterfactual risk deltas from RSF S-learner and a treatment-arm T-learner. Risk: medium. Admissibility: CLEARED.
- `permutation_null_rmst`: Add a validation-only treatment-label permutation null for RMST alignment. Risk: low. Admissibility: CLEARED.

### CODE
- `uno_ipcw_cindex`: Add Uno's IPCW C-index alongside Harrell's. Risk: low. Admissibility: CLEARED.
- `bootstrap_finalize_ci`: At finalization time, compute 95% percentile bootstrap CIs on sealed-test metrics. Risk: low. Admissibility: CLEARED.
- `artifact_manifest`: Add a manifest hash for feature names, genes, and model params. Risk: low. Admissibility: CLEARED.
- `seed_panel_report`: Keep per-seed RSF validation scores in run metadata. Risk: low. Admissibility: CLEARED.

### PARAM
- `smaller_leaf`: Explore `min_samples_leaf` below 10 for more flexible HTE. Risk: medium. Admissibility: CLEARED.
- `more_trees`: Increase `n_estimators` to 1000-2000 for final candidate reruns. Risk: low but slow. Admissibility: CLEARED.
- `aggressive_mtry`: Try `max_features` fractions from 0.1 to 0.4. Risk: low. Admissibility: CLEARED.

## Finalization
Do not run `../finalize-rsf.py` or `../finalize-all.py` during the autonomous loop. They are HUMAN ONLY and may load `affyfRMATest.csv` after the user adds that file locally.
