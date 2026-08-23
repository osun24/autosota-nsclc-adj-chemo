# Red lines

1. Test data is human-only and may not be opened, enumerated, hashed, or used.
2. Only `train.py` is agent-editable code; only `log.md` may be appended.
3. Selection, imputation, scaling, modules, weights, and forests are fit-only.
4. ACT selects the training arm and is never an input feature.
5. Both arm learners predict every assessment patient; recommendations use the
   difference in their counterfactual 60-month RMST predictions.
6. The IPCW-AIPW estimand, paired bootstrap, repeat penalty, and gates are fixed.
7. Censored patients remain in policy evaluation.
8. Gene selection targets fit-only treatment-benefit signals; no hard-coded
   genes, assessment rankings, patient indices, or oracle features.
9. The locked clinical T-learner cannot be weakened or threshold-tuned.
10. Arm support, seed agreement/correlation, benefit scale, ACT exclusion, and
    per-arm genomic split fractions are mandatory diagnostics.
11. Smoke is infrastructure-only. Every adaptive scientific probe consumes a
    full experiment; failures and timeouts consume a slot.
12. Eligible reward is the leaderboard; a diagnostic leader is not a winner.
13. At most 20 attempts, 35 minutes each, and 1,000 trees per forest.
14. Claims remain observational and require the usual causal assumptions.

