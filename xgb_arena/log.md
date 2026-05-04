# XGBoost Arena Log

Structured autonomous-loop entries start below this header. The Phase 1 smoke test is intentionally not logged as an iteration.

---

### iter_001 — stability_selection_genes
- type: ALGO
- idea_id: stability_selection_genes
- hypothesis: Replacing univariate Cox ranking with LASSO-Cox stability selection (100 half-subsamples, ≥60% frequency threshold) will produce a more robust gene feature pool, improving both val_ci and val_rmst_diff.
- changed_files: xgb_arena/train.py

**Red-line audit (pre-edit):**
- Gene selection uses `train_df` only; `valid_df` is not consulted during stability selection. ✓
- `rank_genes_univariate` is called inside `run()` using `train_df` (train CSV). The new `stability_selection_genes` replaces this call with identical scoping. ✓
- No path that touches `affyfRMATest.csv` is added or implied. ✓
- `valid_df` continues to be used only for final evaluation (`evaluate_on_valid`). ✓
- IPTW, scalers, and encoders remain fitted on boot training rows only (unchanged). ✓
- Censored patients are not dropped (unchanged). ✓
- Treatment recommendation logic (counterfactual risk comparison) is unchanged. ✓

- val_ci: 0.6792 ± 0.0093
- val_rmst_diff: 3.82 ± 0.64 (months)
- n_features: 42 (k_main=16, k_int=8)
- verdict: BASELINE (no prior to compare)
- one_line_lesson: Stability selection found only 5/2261 genes at ≥60% threshold (alpha too aggressive); it still acts as a soft ranker, and the final model stopped at best_ntree=8 — very early stopping suggests the XGB is under-training or the ES criterion is misaligned.

*Notes:* best_ntree=8 on final refit is anomalous; Optuna bootstrap median was 214 trees for the chosen trial. RMST diff locked at 3.82 across most trials, suggesting the model mostly recommends OBS for all patients.
