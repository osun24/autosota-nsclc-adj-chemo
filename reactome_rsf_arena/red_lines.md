# Red lines

1. **Test is human-only.** The autonomous loop may use train and validation.
   It must not open, enumerate, hash, summarize, infer, or reference any test
   dataset or test-derived artifact. Test evaluation occurs once after search.
2. **Only `train.py` is agent-editable Python.** The agent may append `log.md`.
   It must not change or bypass locked preparation, metrics, budgets, launcher,
   integrity manifest, ledger, data, or prior run artifacts.
3. **All adaptive preprocessing is fit-only.** Imputation, propensity and
   censoring models, gene ranking, pathway scoring, feature selection, and RSFs
   are fitted without assessment-fold or validation outcomes. Validation is
   prediction/evaluation only.
4. **Never drop censored patients.** The IPCW-AIPW objective uses every patient.
5. **The estimand is immutable.** `tau=60`, `A60(d)=V60(d)-V60(1-d)`, the
   genomic-minus-clinical paired increment, train/validation gap penalty,
   multiplicity LCB, bootstrap draws, and gates cannot change during search.
6. **Recommendations use counterfactual RMST.** Integrate RSF survival curves
   under ACT=1 and ACT=0. Do not use `RandomSurvivalForest.predict()` mortality
   scores or floating-point signs as treatment benefit.
7. **No explicit gene-by-ACT products.** The RSF must learn treatment-feature
   interactions through its tree structure. Candidate schemas, custom selectors,
   and feature builders must not add duplicated or explicit `gene*ACT` columns.
8. **ACT-use reporting is mandatory.** Every completed iteration must document
   the locked clinical-plus-genomic ACT split-tree fraction, counterfactual
   patient-tree path-traversal fraction, and per-patient different-terminal-node
   tree-fraction distribution on both train OOF and validation. These are
   diagnostics and must not be presented as proof of treatment-effect validity.
9. **Genes must be predictive, not merely prognostic.** Selection must target a
   train-only treatment-benefit pseudo-outcome. No hard-coded gene lists,
   validation-ranked genes, patient indices, predictions, or oracle features.
10. **Matched comparator.** Clinical and C+G models use the same folds, patients,
   seeds, weights, RSF parameters, resamples, nuisance procedure, and clinical
   variables. Only the gene feature block differs.
11. **No metric shopping.** `reward` is the leaderboard. C-index, raw increments,
   ACT rate, and cohort-specific values are gates/diagnostics, not substitute
   objectives.
12. **Budget is final.** At most 20 unique full experiments and 25 minutes per
   experiment. Failures/timeouts consume a slot. Do not run unlogged sweeps,
   duplicate candidates, direct `train.py` scientific runs, or rewrite history.
13. **Claims remain observational.** No regimen-level, causal certainty, or
   clinical deployment claim is allowed without adequate external/randomized
    evidence. Similar validation and test results do not prove transportability
    to a new study or platform.
