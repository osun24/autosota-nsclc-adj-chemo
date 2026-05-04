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

---

### iter_002 — wider_eta
- type: PARAM
- idea_id: wider_eta
- hypothesis: Expanding eta from [0.01, 0.12] to [0.005, 0.15] allows slower-learning models with more trees, addressing the anomalous best_ntree=8 and potentially improving generalisation.
- changed_files: xgb_arena/train.py

**Red-line audit (pre-edit):**
- Parameter range change only; no new data paths or features introduced. ✓
- Gene selection, IPTW, val split all unchanged. ✓
- No test-set paths introduced. ✓

- val_ci: 0.6599 ± 0.0355
- val_rmst_diff: 3.82 ± 1.46 (months)
- n_features: 58 (k_main=32, k_int=8)
- verdict: MIXED (val_ci dropped 0.019 vs baseline; rmst same)
- one_line_lesson: Wider eta made things slightly worse on CI (best_ntree 8→12, still anomalously low); the root cause is not the eta range but likely the ES-vs-validation overlap — val_df is used for both early stopping and evaluation.

---

### iter_003 — stratified_es_split
- type: CODE
- idea_id: stratified_es_split
- hypothesis: Carving a 20% ES monitor from train_df (stratified by event×ACT), using it for all early stopping, and reserving val_df purely for evaluation will eliminate the ES/eval overlap that causes best_ntree=8 anomaly.
- changed_files: xgb_arena/train.py

- val_ci: 0.6498 ± 0.0114
- val_rmst_diff: 2.56 ± 0.0 (months)
- n_features: 114 (k_main=96, k_int=0)
- verdict: WORSE (CI −0.029 vs baseline; RMST −1.26 months) → REVERTED
- one_line_lesson: Separating ES from val_df fixed best_ntree (8→342) but caused RMST collapse because Optuna found k_int=0 (no interactions) when optimizing against internal ES — the interaction hyperparameter needs a lower bound or RMST should be computed during bootstrap to guide k_int selection.

---

### iter_004 — t_learner_split
- type: ALGO
- idea_id: t_learner_split
- hypothesis: Fitting separate XGB-Cox models on ACT=0 and ACT=1 subsets produces better counterfactual risk estimates than the S-learner, breaking the RMST=3.82 lock where the model recommends OBS for everyone.
- changed_files: xgb_arena/train.py

**Red-line audit (pre-edit):**
- model_0 fitted on ACT=0 rows of bootstrap(train_df) only; model_1 on ACT=1 rows of bootstrap(train_df) only. ✓
- Gene ranking unchanged; uses same stability-selection train_df gene rank. ✓
- Features: clin_pretx + genes_main (ACT column excluded from T-learner feature matrix). ✓
- val_df used only for evaluation (CI and RMST); no val_df rows enter the fit or ES phase. ✓
- No test-set paths introduced. ✓
- Censored patients not dropped in either arm. ✓
- Counterfactual: predict_0(x) and predict_1(x) for all val patients; recommend ACT if risk_1 < risk_0. ✓
- CI: computed from average risk (predict_0 + predict_1)/2 across val patients. ✓

*Result to be filled in after run.*
