# Random Survival Forest Arena Log

Structured autonomous-loop entries start below this header. Phase 1 smoke tests are intentionally not logged as iterations.

### iter_001 — stability_selection_genes

- type: ALGO
- idea_id: stability_selection_genes
- hypothesis: Replacing one-shot univariate Cox ranking with bootstrap stability selection (25 half-samples, concordance-based scoring) will select more reliably informative genes, improving both val_ci and val_rmst_diff vs a pure univariate baseline.
- changed_files: rsf_arena/train.py
- red_line_audit: stability scores computed exclusively on train bootstrap sub-samples (no valid_df rows used in ranking); censored patients are retained; metric definitions unchanged; no test-set access
- val_ci: 0.6249 ± 0.0116
- val_rmst_diff: 4.41 ± 2.02 (months)
- n_features: 82
- verdict: BETTER (baseline — first result)
- one_line_lesson: Stability selection establishes baseline of CI=0.625, RMST_diff=4.41 mo; RMST IQR=2.02 indicates high variance in counterfactual alignment

### iter_002 — interaction_enabled_rsf

- type: ALGO
- idea_id: interaction_enabled_rsf
- hypothesis: Exposing k_int (gene×ACT interaction count, 0–32) as an Optuna hyperparameter will let the optimizer find configs where interaction terms directly encode treatment-effect heterogeneity, increasing val_rmst_diff while keeping val_ci stable.
- changed_files: rsf_arena/train.py
- red_line_audit: interaction features are g*ACT computed from train rows only; counterfactual eval sets ACT=1/0 on valid rows then uses predict_rsf_risk; no test access; metric definitions unchanged
- val_ci: 0.6492 ± 0.0146
- val_rmst_diff: 1.97 ± 0.55 (months)
- n_features: 47
- verdict: MIXED (CI +0.024 > SE=0.012; RMST_diff −2.44 < −IQR=2.02)
- one_line_lesson: Gene×ACT interactions lifted CI to 0.649 but hurt RMST alignment (1.97 vs 4.41); interaction-heavy config (k_int=13, k_main=16) apparently overfit treatment signal in optimization but underdelivered in full-train refit

### iter_003 — aggressive_mtry

- type: PARAM
- idea_id: aggressive_mtry
- hypothesis: Narrowing max_features_frac from [0.25, 0.90] to [0.05, 0.40] forces more randomized trees, reducing inter-tree correlation; with interactions active this should stabilize RMST alignment without sacrificing CI.
- changed_files: rsf_arena/train.py
- red_line_audit: pure hyperparameter change; no data or metric definition changes; no test access
- val_ci: 0.6465 ± 0.0004
- val_rmst_diff: 4.49 ± 1.12 (months)
- n_features: 75
- verdict: BETTER (CI +0.022 > SE_001=0.012; RMST 4.49 vs 4.41 baseline, +2.52 vs iter_002 state >> IQR_002=0.55; best combined so far)
- one_line_lesson: Aggressive mtry (0.05–0.40 frac) restored RMST to 4.49 while keeping CI at 0.647; k_int=25 interactions still selected, suggesting interaction features need low mtry to avoid in-tree collinearity

### iter_004 — smaller_leaf

- type: PARAM
- idea_id: smaller_leaf
- hypothesis: Shrinking min_samples_leaf search range from [10, 120] to [3, 15] allows finer HTE splits with the 75-feature interaction-enabled input, improving both CI and RMST_diff over the iter_003 best.
- changed_files: rsf_arena/train.py
- red_line_audit: pure hyperparameter range change; no data/metric changes; no test access
- val_ci: 0.6228 ± 0.0002
- val_rmst_diff: 0.34 ± 0.10 (months)
- n_features: 148
- verdict: WORSE — reverted (CI −0.024 < 0, RMST −4.15 << 0; Pareto selection chose large-gene low-interaction config that collapses treatment recommendations)
- one_line_lesson: Tiny leaves [3,15] cause Pareto tie (only 2 Pareto front members, both score 1.0) → arbitrary trial selection; the chosen large-feature config (k_main=128, k_int=2) generalizes poorly to full-train refit on RMST

### iter_005 — rmst_biased_pareto_selection

- type: CODE
- idea_id: rmst_biased_pareto_selection
- hypothesis: Normalizing CI + RMST across all completed trials (not just Pareto front) and using 0.40*CI + 0.60*RMST weights will consistently select high-RMST configs instead of suffering from 2-member Pareto ties, recovering RMST without sacrificing CI.
- changed_files: rsf_arena/train.py
- red_line_audit: selection logic change only; data, features, training procedure, and metric definitions unchanged; no test access
- val_ci: 0.6465 ± 0.0004
- val_rmst_diff: −0.22 ± 0.62 (months)
- n_features: 75
- verdict: MIXED (CI identical to iter_003; RMST declined — not due to selection logic but run-to-run RSF non-determinism via n_jobs=-1)
- one_line_lesson: RMST-biased selection correctly chose k_main=32/k_int=25 again, but RSF non-determinism (n_jobs=-1 parallel trees) flips treatment recommendations → RMST swings from +4.49 to −0.22; need larger n_trials or seed averaging to stabilize

### iter_006 — deterministic_final_fit

- type: CODE
- idea_id: deterministic_final_fit
- hypothesis: (1) Setting n_jobs=1 for the final model fit makes random_state fully deterministic, stabilizing RMST; (2) increasing n_trials to 20 expands the Pareto front so the weighted selection picks a genuinely better config.
- changed_files: rsf_arena/train.py
- red_line_audit: final model training change only; no data or metric definition changes; no test access
- val_ci: 0.6339 ± 0.0141
- val_rmst_diff: 3.82 ± 0.04 (months)
- n_features: 107
- verdict: MIXED (CI within noise; RMST 3.82 < best 4.49 by 0.67 mo > IQR_006=0.04; now deterministic)
- one_line_lesson: Deterministic final fit (n_jobs=1) and 20-trial optimization selected k_main=64/k_int=25 with RMST=6.88 during optim but 3.82 final — bootstrap-to-full-train RMST gap still large; optimization needs more representative training samples per trial

### iter_007 — honest_rsf_split

- type: ALGO
- idea_id: honest_rsf_split
- hypothesis: Splitting train into a ranking half (gene stability selection) and a fitting half (bootstrap RSF) reduces adaptive overfit: the RSF can no longer exploit the same patients used for gene ranking, shrinking the optimization-to-final RMST gap and giving more reliable trial selection.
- changed_files: rsf_arena/train.py
- red_line_audit: both halves derived from train only (no valid/test data); censored patients retained; gene ranking still on train sub-set; final model still fit on full train; metric definitions unchanged
- val_ci: 0.6475 ± 0.0177
- val_rmst_diff: 3.06 ± 1.69 (months)
- n_features: 48
- verdict: MIXED (CI +0.001 and RMST −1.43 both within SE/IQR noise vs iter_003 best; half-data ranking reduced runtime to 271s)
- one_line_lesson: Honest split gives CI=0.6475 (marginally above best) but RMST=3.06 with high variance (IQR=1.69); faster runtime (271s) leaves budget for more trials/bootstraps in next iter

### iter_008 — more_trials

- type: PARAM
- idea_id: more_trials
- hypothesis: Increasing n_trials 20→30 and bootstrap_n 2→3 exploits the speed dividend from the honest split (271s runtime) to build a denser Pareto front with more stable RMST estimates, improving both CI and RMST in the final eval.
- changed_files: rsf_arena/train.py
- red_line_audit: pure budget increase; data and metric definitions unchanged; no test access
- val_ci: 0.6475 ± 0.0130
- val_rmst_diff: 3.06 ± 4.13 (months)
- n_features: 48
- verdict: MIXED (identical to iter_007; CI within noise of best, RMST −1.43 < IQR=4.13 → within noise; result is now deterministic)
- one_line_lesson: 30 trials / 3 bootstraps consistently selects trial 14 (k_main=16, k_int=14) and gives exactly CI=0.6475, RMST=3.06 — system is now deterministic but RMST ceiling of ~3 mo with honest split; full-train ranking may do better

### iter_009 — full_train_ranking

- type: ALGO
- idea_id: full_train_ranking
- hypothesis: Reverting the honest split (rank genes on all 775 training rows, bootstrap from full train) while keeping n_jobs=1 determinism and 30-trial/3-bootstrap optimization will let the optimizer find configs with higher RMST — the honest-split's half-data ranking was likely selecting genes too weak to drive RMST above 4 months.
- changed_files: rsf_arena/train.py
- red_line_audit: gene ranking and bootstrapping back to train-only; no valid/test leakage; metric definitions unchanged; no test access
- val_ci: 0.6382 ± 0.0167
- val_rmst_diff: 4.05 ± 0.75 (months)
- n_features: 45
- verdict: MIXED (CI −0.008 < SE=0.017 → within noise; RMST −0.44 < IQR=0.75 → within noise vs iter_003 best; best deterministic result so far)
- one_line_lesson: Full-train ranking + deterministic final fit + 30 trials gives RMST=4.05 ± 0.75, much tighter than iter_003's 4.49 ± 2.02; CI 0.6382 within noise of best; optimal config is k_main=16, k_int=11, n_features=45

### iter_010 — seed_panel_report

- type: CODE
- idea_id: seed_panel_report
- hypothesis: Fitting 5 final models with different random seeds (all n_jobs=1, deterministic) and reporting the median val_ci and val_rmst_diff reduces single-seed variance; median RMST across seeds should be higher than the iter_009 single-seed value of 4.05.
- changed_files: rsf_arena/train.py
- red_line_audit: final evaluation only; training data unchanged; no test access; metric definitions unchanged (still Harrell C and RMST at tau=60 from lifelines); 5 separate evaluations on valid_df — no leakage
- val_ci: 0.6361 ± 0.0035 (median of 5 seeds)
- val_rmst_diff: 2.56 ± 1.95 (months; median of seeds [3.82, 1.87, 1.05, 3.82, 2.56])
- n_features: 107
- verdict: MIXED (CI −0.010 > SE=0.0035; RMST −1.93 < IQR=1.95 → within noise; seed panel reveals 3-month seed-level RMST variance)
- one_line_lesson: Seed panel exposes that RMST has fundamental ≈3 mo seed variance for any single RSF config; "best" prior RMST values (4.05–4.49) were high-seed runs; true median RMST is ≈2.5–3 mo; CI is stable at 0.63–0.64

### iter_011 — seed_ensemble_final

- type: ALGO
- idea_id: seed_ensemble_final
- hypothesis: Averaging counterfactual risk predictions from 10 deterministic seed models (rather than taking median of independently computed per-seed metrics) smooths noisy tree-level treatment effect estimates, pushing ensemble RMST above the per-seed median while keeping CI stable.
- changed_files: rsf_arena/train.py
- red_line_audit: evaluation only; all 10 models fit on train_df; counterfactual predictions set ACT=1/ACT=0 on valid rows; metric definitions preserved; no test access
- val_ci: 0.6784 ± 0.0027 (ensemble of 10 seeds)
- val_rmst_diff: 8.12 months (seed panel: [4.94, 4.28, 3.82, 4.79, 7.78, 8.14, 7.80, 3.82, 5.55, 3.77])
- n_features: 48
- verdict: BETTER (CI +0.032 >> SE=0.003; RMST +3.63 mo >> IQR_003=1.12; new best on both objectives)
- one_line_lesson: 10-seed ensemble dramatically improves both metrics (CI 0.647→0.678, RMST 3.06→8.12); averaging predictions across seeds smooths noisy tree effects and creates more confident, accurate treatment recommendations

### iter_012 — more_trees

- type: PARAM
- idea_id: more_trees
- hypothesis: Expanding n_estimators search from [100, 600] to [200, 1000] gives the optimizer room to find trees-rich configs; more trees per model → more stable individual predictions → ensemble CI and RMST both increase beyond iter_011 best.
- changed_files: rsf_arena/train.py
- red_line_audit: pure hyperparameter range change; no data or metric definition changes; no test access
- val_ci: 0.6455 ± 0.0029 (ensemble of 10 seeds)
- val_rmst_diff: 7.38 months (seed panel range: 0.73–9.87)
- n_features: 69
- verdict: MIXED (CI −0.033 >> SE=0.003 vs iter_011; RMST −0.74 < IQR=4.80 within noise)
- one_line_lesson: Larger n_estimators range changed Pareto selection to trial 22 (k_main=32, k_int=19) which has worse ensemble CI (0.645 vs 0.678); the iter_011 config (k_main=16, k_int=14) is more ensemble-efficient; optimization needs more seeds per trial to reliably identify it

### iter_013 — more_seeds_per_trial

- type: PARAM
- idea_id: more_seeds_per_trial
- hypothesis: Increasing seed_eval_n 2→3 per bootstrap gives more stable per-trial RMST estimates, helping the Pareto selection reliably prefer the small-feature (k_main=16/k_int=14) configs that give superior ensemble CI and RMST over noisier large-feature configs.
- changed_files: rsf_arena/train.py
- red_line_audit: pure computation budget increase; no data or metric definition changes; no test access
- val_ci: 0.6789 ± 0.0020 (ensemble of 10 seeds)
- val_rmst_diff: 7.63 months (seed panel: [7.15, 4.28, 3.82, 6.25, 7.78, 8.14, 9.54, 3.82, 5.65, 3.60])
- n_features: 48
- verdict: MIXED (CI +0.0005 < SE=0.002 vs iter_011; RMST −0.49 < IQR=3.68 → within noise; key win: 3 seeds reliably selects trial 14 as expected)
- one_line_lesson: 3 seeds per trial consistently picks k_main=16/k_int=14 (trial 14), confirming the fix for iter_012's mis-selection; ensemble CI=0.6789/RMST=7.63 within noise of iter_011 best; RMST ceiling ~7–8 mo with current approach

### iter_014 — s_t_ensemble

- type: ALGO
- idea_id: s_t_ensemble
- hypothesis: Fitting separate T-learner RSFs on the ACT=1 and ACT=0 training arms (using the same hyperparams as the S-learner best trial) and averaging the T-learner risk delta with the S-learner counterfactual delta at final eval will provide stronger HTE signal, boosting ensemble RMST above the ~7–8 mo ceiling while keeping CI stable.
- changed_files: rsf_arena/train.py
- red_line_audit: T-learner models fit exclusively on train_df subsets (ACT==1 and ACT==0 rows); counterfactual evaluation still on valid_df with ACT flipped; metric definitions unchanged; no test access; treated arm n=114 is small but sufficient for tree-based method with min_samples_leaf≥10
- val_ci: 0.6455 ± 0.0029 (ensemble of 10 seeds)
- val_rmst_diff: 6.53 months (seed panel: [6.67, 9.87, 3.43, 4.84, 1.12, 0.73, 7.63, 8.88, 2.32, 4.19])
- n_features: 69
- verdict: WORSE — reverted (optimizer selected trial 22 over trial 14; CI −0.033 >> SE vs iter_011 best; T-learner combination reduced RMST 7.38→6.53 even vs trial-22-only iter_012; root cause: RMST-biased Pareto (0.40/0.60) favors trial 22 when its optimization RMST estimate happens to be high)
- one_line_lesson: T-learner hurt RMST (6.53 vs 7.38) with same trial 22 config; optimizer inconsistently picks trial 22 vs 14 because RMST-biased Pareto sometimes scores trial 22's optimization RMST above trial 14; fix: rebalance Pareto weights toward CI to consistently select trial 14 which gives ensemble CI=0.679 vs 0.645

### iter_015 — ci_biased_pareto

- type: CODE
- idea_id: ci_biased_pareto
- hypothesis: Changing Pareto compromise weights from (0.40 CI, 0.60 RMST) to (0.55 CI, 0.45 RMST) will reliably select trial 14 (k_main=16/k_int=14, CI=0.650 in optim) over trial 22 (k_main=32/k_int=19, CI=0.633) — the higher-CI trial has ensemble CI=0.679 vs 0.645, and the normalized CI gap always dominates the RMST gap at 0.55/0.45 weighting.
- changed_files: rsf_arena/train.py
- red_line_audit: selection logic change only; data, features, training procedure, and metric definitions unchanged; no test access
