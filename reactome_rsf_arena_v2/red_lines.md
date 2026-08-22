# Red lines

1. **Test is human-only.** The loop may load only the former train and
   validation CSVs as pooled development data. It must not open, enumerate,
   hash, summarize, infer, or reference any test dataset or test-derived file.
2. **Only `train.py` is agent-editable Python.** The agent may append `log.md`.
   It must not change locked preparation, scoring, budgets, launcher, manifest,
   ledger, prior artifacts, or data.
3. **All adaptive transforms are fit-only.** Imputation, scaling, module
   construction, nuisances, gene/pathway selection, and RSFs use no assessment
   outcomes or covariate summaries.
4. **Never drop censored patients.** IPCW-AIPW evaluation includes everyone.
5. **The estimand is immutable.** Tau, policy alignment, paired genomic
   increment, worst-repeat LCB, repeat-range penalty, bootstrap, and gates do
   not change during search.
6. **Recommendations use counterfactual RMST.** Mortality predictions and
   numerical floating-point signs are not treatment benefits.
7. **No explicit gene-by-ACT products.** Native forest splits or locked
   fit-only module features must express effect modification.
8. **Genes must be predictive candidates, not merely prognostic.** Selection
   targets train-only treatment-benefit pseudo-outcomes. No hard-coded symbols,
   patient indices, assessment rankings, or oracle features.
9. **The comparator cannot be weakened.** Clinical geometry and its zero-month
   threshold are locked. Only genomic forest settings and selector logic in
   `train.py` may change.
10. **ACT-use reporting is mandatory.** Every iteration records tree-split,
    patient-path, and different-terminal-node diagnostics for C+G in the pooled
    OOF summary and each repeat. These are diagnostics, not causal proof.
11. **No unbudgeted adaptive probes.** Any structural, expression, outcome,
    treatment, policy, or model probe intended to guide a candidate must run
    through the full launcher and consume a slot. Locked smoke is infrastructure
    verification only.
12. **No metric shopping.** Eligible `reward` is the leaderboard. The
    diagnostic leader guides hypotheses but is not a frozen candidate.
13. **Budget is final.** At most 20 unique full attempts and 25 minutes each;
    failures and timeouts consume a slot. Do not rerun identical `train.py`.
14. **Claims remain observational.** No treatment-effect certainty, regimen
    claim, or deployment recommendation is allowed without appropriate
    external or randomized evidence.
