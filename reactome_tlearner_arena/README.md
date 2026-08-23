# Reactome RSF T-learner arena

This arena follows `../reactome_rsf_arena_v2` but replaces each survival-forest
S-learner with a T-learner: an observation-arm RSF and an ACT-arm RSF are fit
separately in every outer training fold. ACT is absent from `X`; the policy is
`1{RMST1(x) - RMST0(x) > threshold}` at 60 months.

The former train and validation cohorts remain pooled adaptive development
data. Two repeated four-fold cross-fits produce held-out policies, and the
same locked IPCW-AIPW arm scores grade the clinical and clinical-plus-genomic
policies. The objective remains the worst selection-adjusted repeat LCB minus
the repeat increment range. The clinical-only T-learner is independently
locked, including stronger regularization for its smaller ACT arm.

The model records arm-level patient/event support, seed policy agreement,
seed benefit correlation, benefit scale, and genomic split fractions. The old
S-learner ACT split/path/terminal diagnostics are inapplicable because ACT is
never a feature.

```bash
/Users/owensun/miniconda3/bin/python -m reactome_tlearner_arena.run --smoke
/Users/owensun/miniconda3/bin/python -m reactome_tlearner_arena.run
```

Smoke uses one repeat, two folds, one seed, 40 trees per forest, and 100
bootstraps. It does not consume the 20-run ledger and is SHA-locked to the
shipped baseline. Full runs fit at most 96 forests / 96,000 trees and have a
35-minute wall limit. Test data remains human-only after an eligible winner is
frozen.

