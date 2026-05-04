# Project: NSCLC adjuvant chemo HTE — DeepSurv arena

## Goal
Improve the `train.py` pipeline along TWO objectives jointly:
- `val_ci` (Val Harrell's C-index, higher better)
- `val_rmst_diff` (Val 5-yr RMST diff under counterfactual recommendation, higher better)
Both must improve beyond bootstrap uncertainty for a candidate to count as BETTER.

## Hard Rules
See `red_lines.md`. The test set remains sealed; `prepare.py` and `train.py` must only load Train, Validation, and `LOOCV_Genes2.csv`.

## Runtime Note
This arena is GPU-friendly but CPU-safe. On this machine, base Python reports CUDA unavailable, so the default budget is conservative. If a CUDA GPU is available in a later session, the human may raise `DEEPSURV_ARENA_EPOCHS`, `DEEPSURV_ARENA_N_TRIALS`, and `DEEPSURV_ARENA_BOOTSTRAPS` before launch.

## Per-Iteration Workflow
1. Read `log.md`. Identify the last 3 iterations and their `type` (PARAM | CODE | ALGO).
2. If all 3 most recent COMPLETED iterations were PARAM, propose a CODE or ALGO idea this iteration.
3. Pick one CLEARED idea from the idea library below, or generate a new one.
4. Write an idea entry to `log.md` with type, hypothesis, expected effect, and a red-line audit BEFORE editing `train.py`.
5. `git add . && git commit -m "iter_NNN_pre: <idea_id>"` to snapshot pre-state.
6. Edit `train.py`. Run `python deepsurv_arena/train.py`. Wall clock budget: 25 min on CPU unless the human approves a GPU budget.
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
- `t_learner_deepsurv`: Fit separate DeepSurv models for ACT and OBS, then compare counterfactual risks. Risk: high (small ACT subgroup). Admissibility: CLEARED after red-line audit.
- `s_t_ensemble`: Average ITE predictions from S-learner and T-learner DeepSurv. Risk: medium once T-learner exists. Admissibility: CLEARED.
- `pathway_first_layer`: Replace gene inputs with train-only GSVA/Hallmark pathway features. Risk: medium dependency risk. Admissibility: CLEARED pending dependency check.
- `group_l1_interactions`: Use separate first-layer L1 penalties for main genes and gene x ACT interactions. Risk: low. Admissibility: CLEARED.
- `permutation_null_rmst`: Add a validation-only treatment-label permutation null for RMST alignment. Risk: low. Admissibility: CLEARED.

### CODE
- `uno_ipcw_cindex`: Add Uno's IPCW C-index alongside Harrell's. Risk: low. Admissibility: CLEARED.
- `artifact_manifest`: Add a manifest hash for feature names, genes, network params, and transform arrays. Risk: low. Admissibility: CLEARED.
- `seed_panel_report`: Report multiple random seeds for the final candidate without changing the objective definitions. Risk: low. Admissibility: CLEARED.
- `batch_riskset_training`: Move from full-risk-set CPU training to event-balanced mini-batches for larger features. Risk: medium. Admissibility: CLEARED.

### PARAM
- `wider_arch`: Try `128`, `128-64`, and `128-64-32` hidden layers. Risk: medium on CPU. Admissibility: CLEARED.
- `stronger_dropout`: Raise dropout range to 0.2-0.6. Risk: low. Admissibility: CLEARED.
- `longer_epochs_gpu`: Increase epochs to 200-500 only when CUDA is available. Risk: low with GPU, high on CPU. Admissibility: CLEARED pending hardware.

## Finalization
Do not run `../finalize-deepsurv.py` or `../finalize-all.py` during the autonomous loop. They are HUMAN ONLY and may load `affyfRMATest.csv` after the user adds that file locally.
