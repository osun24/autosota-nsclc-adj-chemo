# Project: NSCLC adjuvant chemo HTE — XGBoost arena

## Goal
Improve the `train.py` pipeline along TWO objectives jointly:
- `val_ci` (Val Harrell's C-index, higher better)
- `val_rmst_diff` (Val 5-yr RMST diff under counterfactual recommendation, higher better)
Both must improve beyond bootstrap uncertainty for a candidate to count as BETTER.

## Hard Rules
1. **Test set is sealed.** `affyfRMATest.csv` is read exactly once, at the very end, by the human. The agent's `prepare.py` must not load, reference, or expose test data in any form. Any code path that opens the test file during the loop is a violation.
2. **No leakage in any preprocessing step.** Gene ranking, stability selection, IPTW propensity model, scalers, encoders, feature selection — fit on the bootstrap's training rows only. Validation rows pass through transforms fitted on training rows.
3. **Never drop censored patients.** All n in the analytic cohort are used.
4. **Metric definitions are versioned.** Harrell's C from `sksurv.metrics.concordance_index_censored`, RMST at tau=60 months from `lifelines.utils.restricted_mean_survival_time`. If the agent wants to *add* a metric (e.g., Uno's IPCW C, C-for-benefit from Van Klaveren 2018), it adds it as an additional column in `log.md`, never replacing the existing two.
5. **Treatment recommendation must be counterfactual.** Predict risk under ACT=1 and ACT=0 for every patient; recommend the lower-risk arm. Do not regress on observed outcomes within treatment arms as a shortcut.
6. **No regimen-level claims.** ACT is binary. Don't infer or assume cisplatin vs. carboplatin, doublet vs. single-agent, dose intensity, or schedule effects.
7. **Improvement claims require both objectives to hold.** A new candidate is "BETTER" only if val_ci AND val_rmst_diff each beat the previous best by more than the bootstrap SE / IQR of the metric. One-objective wins are MIXED, not BETTER.
8. **No hardcoded predictions, no oracle features, no test-set peeking.** If the agent finds itself "fixing" the metric to climb the score, that's a violation.

## Per-Iteration Workflow
1. Read `log.md`. Identify the last 3 iterations and their `type` (PARAM | CODE | ALGO).
2. If all 3 most recent COMPLETED iterations were PARAM, you MUST propose a CODE or ALGO idea this iteration. This is the Leap Path rule from AutoSOTA.
3. Pick one CLEARED idea from the idea library below, or generate a new one.
4. Write an idea entry to `log.md` with type, hypothesis, expected effect, and a red-line audit BEFORE editing `train.py`.
5. `git add . && git commit -m "iter_NNN_pre: <idea_id>"` to snapshot pre-state.
6. Edit `train.py`. Run `python xgb_arena/train.py`. Wall clock budget: 25 min.
7. On completion, append the result to `log.md` in the schema below.
8. If verdict == WORSE on both objectives, `git revert HEAD` and choose a different idea. If MIXED, leave the change in place but flag it; pursue at most 3 follow-up iterations to debug before rolling back (the AutoSOTA honeymoon period).

## log.md Entry Schema
Use this exactly.

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
- `t_learner_split`: Fit separate XGB-Cox on ACT=1 and ACT=0 subsets; predict ITE as risk_treated - risk_untreated. Risk: medium (small ACT subgroup, n approximately 223 in full cohort, smaller per bootstrap). Admissibility: CLEARED.
- `c_for_benefit_metric`: Add Van Klaveren's C-for-benefit as a third reported metric. Risk: low (additive, no protocol change). Admissibility: CLEARED.
- `stability_selection_genes`: Replace univariate ranking with stability selection (LASSO-Cox on 100 subsamples at frac 0.5; keep genes selected in at least 60% of runs). Risk: low. Admissibility: CLEARED.
- `gsva_pathway_features`: Reduce 13K genes to approximately 50 Hallmark pathway scores via GSVA before ranking. Risk: medium (requires extra dependency). Admissibility: CLEARED pending dependency check.
- `s_t_ensemble`: Average ITE predictions from S-learner and T-learner. Risk: low once `t_learner_split` is implemented. Admissibility: CLEARED.
- `causal_survival_forest_stack`: Train a causal-survival-forest-style learner outside XGB, then stack its estimated treatment benefit as an extra training-only-derived feature. Risk: high (method complexity, leakage audit required). Admissibility: CLEARED only after a written red-line audit.
- `permutation_null_rmst`: Add a validation-only treatment-label permutation null for RMST alignment so large RMST gains are interpreted against a null. Risk: low (additive). Admissibility: CLEARED.

### CODE
- `enable_dup_inter`: The `dup_inter` parameter is currently dead code (forced to 1). Either re-enable as a tunable, or remove the dead branch. Risk: low. Admissibility: CLEARED.
- `uno_ipcw_cindex`: Add Uno's IPCW C-index alongside Harrell's. Risk: low. Admissibility: CLEARED.
- `bootstrap_test_ci`: At `finalize.py` time, compute 95% percentile bootstrap CIs on test metrics. Risk: low. Admissibility: CLEARED.
- `stratified_es_split`: Stratify the early-stopping split by event x ACT rather than event alone. Risk: low. Admissibility: CLEARED.
- `artifact_manifest`: Add a manifest hash for feature names, genes, and model params in every run directory. Risk: low. Admissibility: CLEARED.

### PARAM
- `wider_eta`: Allow eta in `[0.005, 0.15]`. Risk: low. Admissibility: CLEARED.
- `stronger_subsample`: Set `colsample_bytree` in `[0.3, 0.7]` for more aggressive regularization when p >> n. Risk: low. Admissibility: CLEARED.
- `more_regularized_depth`: Restrict `max_depth` to 2-4 and increase `min_child_weight` lower bound. Risk: low. Admissibility: CLEARED.
- `interaction_budget_probe`: Shift feature budget from main effects to gene x ACT interactions while keeping total feature budget fixed. Risk: medium. Admissibility: CLEARED.

## Anti-Stagnation
If you have gone 5 iterations without a BETTER verdict on both objectives, explicitly try an ALGO move from the library. If the library is exhausted, read the latest survival-HTE literature (causal survival forest, BCF survival, deep counterfactual survival models) and propose a new ALGO entry, with red-line audit, before implementing.

## Finalization
Do not run `../finalize.py` during the autonomous loop. It is HUMAN ONLY and is allowed to load `affyfRMATest.csv` once the user has added that file locally.
