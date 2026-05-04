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
- val_ci: PENDING
- val_rmst_diff: PENDING
- n_features: PENDING
- verdict: PENDING
- one_line_lesson: PENDING
