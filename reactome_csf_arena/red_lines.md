# Red lines

1. **Test is human-only and terminal.** The agent may load only the former
   train and validation CSVs. It must not open, enumerate, hash, summarize,
   infer, or reference test data. After the one human test run, no experiment,
   model revision, threshold change, or candidate substitution is permitted.
2. **Only `train.py` is agent-editable code.** The agent may append `log.md`.
   Locked Python/R code, budgets, manifests, ledgers, artifacts, and data do not
   change.
3. **All adaptation is fit-only.** Imputation, scaling, modules, nuisances,
   gene/pathway selection, and forests use no outer-assessment outcome or
   covariate summary.
4. **ACT is treatment, not X.** It is passed only as `W` to the causal forest.
   No ACT feature or explicit gene-by-ACT product is allowed.
5. **The estimand is immutable.** Horizon 60, zero threshold, policy alignment,
   paired genomic increment, selection-adjusted LCB, repeat-range penalty,
   bootstrap, and gates cannot change during search.
6. **Policy learning and grading remain separate.** GRF predictions create the
   recommendation; common locked assessment-fold IPCW-AIPW scores grade both
   policies. GRF training scores cannot become the reward.
7. **Never drop censored patients.** Evaluation includes every assessment row.
8. **Genes must target effect modification.** Selection uses only fit-partition
   treatment-benefit signals and Reactome structure. No hard-coded symbols,
   patient indices, assessment rankings, or oracle features.
9. **Matched comparison cannot be weakened.** Every tunable forest parameter
   applies identically to clinical and C+G. Honesty, all-feature mtry, seeds,
   and zero threshold are locked.
10. **Diagnostics are mandatory.** Every iteration explicitly documents ACT's
    W-only role, the inapplicability of RSF split/path/terminal diagnostics,
    seed policy agreement, seed CATE correlation, CATE scale, recommendation
    rate, and genomic variable-importance fraction for pooled OOF and repeats.
11. **C-index is secondary.** It is reported from a separate matched prognostic
    survival forest and cannot enter reward; only its locked noninferiority gate
    remains.
12. **No unbudgeted probes or internal tuning.** Scientific feedback runs only
    through the full launcher and consumes a slot. GRF `tune.parameters` stays
    `none`; smoke is infrastructure-only.
13. **No metric shopping.** Eligible reward is the leaderboard. Source-specific
    metrics and the diagnostic leader guide hypotheses but cannot nominate a
    test candidate.
14. **Fixed compute is final.** At most 20 attempts, 30 minutes each, 1,000
    trees per forest, and the locked forest accounting in `budget.json`.
    Failures/timeouts consume a slot; identical `train.py` cannot rerun.
15. **Claims remain observational.** No treatment-effect certainty, regimen
    claim, or deployment recommendation is supported without appropriate
    external or randomized evidence.
