# Reactome-wide RSF permutation-importance arena

This arena runs the requested two-stage Random Survival Forest workflow on the
pooled 1,034-patient development cohort while keeping the test set sealed.

1. Tune an RSF containing the locked clinical covariates and **every eligible
   Reactome gene** shared by train and validation (currently 8,647 before the
   locked nonconstant check).
2. Refit the winning screening configuration in fixed cross-validation folds
   and compute held-out paired permutation feature importance. The primary
   importance is the loss in five-year IPCW-AIPW policy-versus-anti-policy RMST
   difference; loss in Harrell C-index is secondary.
3. Rank genes by those held-out importances, evaluate every top-N panel from
   clinical-only (`N=0`) through `N=32`, and continue hyperparameter optimization
   over N and the reduced RSF geometry.
4. Refit the selected clinical+gene panel on all development rows and emit the
   final covariates, parameters, CV metrics, importance table, and plot.

The clinical columns are always retained. `Adjuvant Chemo` is the treatment
feature used to form paired counterfactual predictions. Its primary PFI is
reported as zero because the policy calculation deliberately overwrites it
with ACT=0 and ACT=1; its observed-row C-index PFI remains defined.

## Commands

```bash
/Users/owensun/miniconda3/bin/python -m reactome_rsf_pfi_arena.run --smoke
/Users/owensun/miniconda3/bin/python -m reactome_rsf_pfi_arena.run
```

Smoke uses 64 genes, two folds, 20-tree forests, one permutation repeat, and
N=0..2. It verifies infrastructure and does not consume the experiment ledger.
A full run can take hours because genuine PFI perturbs thousands of features.

## Run artifacts

- `result.json`: complete machine-readable result and limitations
- `final_spec.json`: chosen top N, genes, covariates, and RSF parameters
- `final_covariates.txt`: clinical variables followed by selected genes
- `permutation_importance.csv`: all clinical and Reactome features, including
  exact zeros for features unused by every tree in their assessment fold
- `permutation_importance.png`: RMST and C-index PFI panels
- `top_n_results.csv`: every reduced-stage CV trial, including all N=0..32
- `final_rsf.pkl`: full-development fitted model
- `training_imputation_medians.csv`: fit transform required for later use

The CV scores are adaptive development estimates, not external confirmation.
Permutation importance measures dependence of this fitted pipeline, not a
causal effect or a uniquely identifiable biological mechanism.

