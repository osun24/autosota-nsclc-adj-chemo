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
- val_ci: PENDING
- val_rmst_diff: PENDING
- n_features: PENDING
- verdict: PENDING
- one_line_lesson: PENDING
