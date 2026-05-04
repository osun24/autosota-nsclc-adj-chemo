# Random Survival Forest Arena Log

Structured autonomous-loop entries start below this header. Phase 1 smoke tests are intentionally not logged as iterations.

### iter_001 — stability_selection_genes

- type: ALGO
- idea_id: stability_selection_genes
- hypothesis: Replacing one-shot univariate Cox ranking with bootstrap stability selection (25 half-samples, concordance-based scoring) will select more reliably informative genes, improving both val_ci and val_rmst_diff vs a pure univariate baseline.
- changed_files: rsf_arena/train.py
- red_line_audit: stability scores computed exclusively on train bootstrap sub-samples (no valid_df rows used in ranking); censored patients are retained; metric definitions unchanged; no test-set access
- val_ci: PENDING
- val_rmst_diff: PENDING
- n_features: PENDING
- verdict: PENDING
- one_line_lesson: PENDING
