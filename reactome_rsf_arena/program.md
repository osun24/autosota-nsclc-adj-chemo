# Autonomous program: Reactome causal-RMST RSF

## Goal

Maximize the single locked `reward` emitted by `python -m
reactome_rsf_arena.run`. The scientific goal is a small, stable Reactome gene
set whose ACT policy increment is positive and similar on training OOF and
validation, so that held-out test performance is plausibly comparable.

Read `README.md`, `red_lines.md`, all of `log.md`, and the current `train.py`
before the first experiment.

## Editable surface

During the loop:

- edit only `train.py`;
- append experiment hypotheses and results to `log.md`;
- never edit `prepare.py`, `run.py`, `integrity.py`, `budget.json`,
  `lock_manifest.json`, the ledger, prior runs, or data.

The simplest experiments change `CANDIDATE`. Algorithm experiments may define
a train-only gene-selector function in `train.py` and pass it to
`prepare.evaluate_candidate`. The callback receives only the fitting partition,
Reactome membership, available training genes, and the validated candidate
specification. It must not contain hard-coded gene symbols.

## Fixed resources

- 20 unique full experiments total;
- 25 minutes maximum per full experiment;
- 3 outer training folds and 3 inner nuisance folds;
- forest seeds 42, 43, and 44;
- at most 1,000 trees per forest;
- at most 32 genes;
- at most 51 total model features (19 clinical + 32 genes);
- 4,000 paired patient bootstraps per cohort.

Do not run `train.py` directly for scientific feedback. Use the launcher so the
budget, timeout, lock, snapshot, and ledger are enforced.

## Iteration workflow

1. Review the best eligible reward and the last three completed attempts.
2. Choose one interpretable hypothesis. After three consecutive PARAM attempts,
   the next attempt must be ALGO or CODE.
3. Append the prespecified hypothesis and red-line audit to `log.md` before
   editing.
4. Edit `train.py` only.
5. Run `python -m reactome_rsf_arena.run` once.
6. Append the run ID, reward, eligibility, failed gates, train increment,
   validation increment, absolute gap, both C-indices, gene Jaccard, and lesson.
   Also report all three locked ACT-use diagnostics for the clinical-plus-genomic
   forest on both train OOF and validation: tree split fraction, patient-tree
   path traversal fraction, and the per-patient terminal-difference distribution.
7. Keep a candidate only if its reward improves. An ineligible result is not a
   scientific win regardless of point estimates.

Do not optimize a secondary diagnostic after the primary reward disappoints.
Use diagnostics to form the next prespecified hypothesis.

## Candidate priorities

Prefer changes in this order:

1. robust gene selection (`dr_gene`, `dr_pathway`, `dr_hybrid`);
2. smaller stable gene counts;
3. clinically prespecified RMST benefit threshold;
4. leaf size, depth, and mtry regularization;
5. tree count only when seed instability suggests Monte Carlo noise.

Do not force an arbitrary ACT recommendation rate. All-OBS or all-ACT behavior
is rejected only when it fails value, effect-size, or evaluability gates.

## Result schema for `log.md`

```markdown
### iter_NNN — title
- type: PARAM | CODE | ALGO
- hypothesis: one sentence
- changed: exact `train.py` change
- red_line_audit: concise audit
- run_id: run_NNN_...
- eligible: true | false
- failed_gates: [...]
- reward: X.XXX
- train_increment: X.XXX months
- validation_increment: X.XXX months
- absolute_gap: X.XXX months
- train_cindex_cg: 0.XXX
- validation_cindex_cg: 0.XXX
- gene_jaccard: 0.XXX
- train_act_usage_cg: tree_split=0.XXX; path_traversal=0.XXX; terminal_difference_mean=0.XXX, median=0.XXX, p10=0.XXX, p90=0.XXX, nonzero_patients=0.XXX
- validation_act_usage_cg: tree_split=0.XXX; path_traversal=0.XXX; terminal_difference_mean=0.XXX, median=0.XXX, p10=0.XXX, p90=0.XXX, nonzero_patients=0.XXX
- verdict: SEARCH_LEADER | NOT_LEADER | INELIGIBLE
- lesson: one sentence
```

The two `act_usage_cg` lines are mandatory in every completed iteration,
including ineligible runs. Copy them from
`models.clinical_plus_genomic.act_usage` in the immutable run result. The
locked evaluator also records the same diagnostics for the clinical-only
reference model; do not substitute those values for the required C+G lines.

## Completion

After 20 experiments—or earlier if no credible improvement remains—freeze the
best eligible run named in `best_run.txt`. A human performs the one-shot test
evaluation. The agent must not request, locate, infer, or access test data.
