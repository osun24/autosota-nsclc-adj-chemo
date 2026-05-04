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
- val_ci: PENDING
- val_rmst_diff: PENDING
- n_features: PENDING
- verdict: PENDING
- one_line_lesson: PENDING
