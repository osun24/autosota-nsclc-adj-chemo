# Autonomous program: Reactome causal-RMST RSF v2

## Goal

Maximize the locked eligible `reward` from `python -m
reactome_rsf_arena_v2.run`. Seek a small, stable Reactome gene module whose
ACT-policy increment is positive in both repeated development cross-fits.

Before the first experiment, read `README.md`, `red_lines.md`,
`v1_findings.md`, all of `log.md`, and `train.py`.

## Editable surface

- Edit only `train.py`.
- Append prespecified hypotheses and completed results to `log.md`.
- Never edit locked files, ledgers, prior run artifacts, or data.
- Custom selectors receive only a fitting partition, Reactome membership,
  available development genes, and the validated candidate specification.
  They must not contain hard-coded gene symbols.

## Fixed resources

- 20 unique full experiments; 25 minutes per experiment.
- Two repeated four-fold outer cross-fits and three inner nuisance folds.
- Forest seeds 42, 43, and 44; no more than 1,000 trees per forest.
- At most 16 selected genes and 35 actual model features.
- 4,000 paired patient bootstraps per repeat.
- Independently locked clinical comparator; agent-tunable genomic forest.

Use only the launcher. Do not run `train.py` or helper probes directly for
scientific feedback.

## Iteration workflow

1. Review `best_run.txt`, `diagnostic_leader.txt`, and the last three attempts.
2. Choose one interpretable hypothesis that does not repeat a refuted v1 idea.
3. Append the hypothesis, exact planned change, and red-line audit before editing.
4. Edit `train.py` only and launch exactly one full run.
5. Complete the entry with eligibility, gates, reward, continuous diagnostic
   score, both repeat increments/LCBs, repeat range, pooled C-index and policy
   alignment, gene Jaccard, and all required ACT-use diagnostics.
6. Treat only an eligible reward improvement as a search leader. Use an
   ineligible diagnostic improvement solely to form the next hypothesis.

After three consecutive parameter-only attempts, the next must change selector
or representation logic. Prefer module representation, stable detectable
pathways, smaller gene sets, genomic threshold calibration, and genomic forest
regularization—in that order. Do not weaken the comparator or force an ACT rate.

## Result schema for `log.md`

```markdown
### iter_NNN — title
- type: PARAM | CODE | ALGO
- hypothesis: one sentence
- changed: exact train.py change
- red_line_audit: concise audit
- run_id: run_NNN_...
- eligible: true | false
- failed_gates: [...]
- reward: X.XXX
- diagnostic_score: X.XXX
- repeat_1_increment: X.XXX months
- repeat_1_lcb: X.XXX months
- repeat_2_increment: X.XXX months
- repeat_2_lcb: X.XXX months
- repeat_range: X.XXX months
- development_alignment_cg: X.XXX months
- development_cindex_cg: 0.XXX
- gene_jaccard: 0.XXX
- development_act_usage_cg: tree_split=0.XXX; path_traversal=0.XXX; terminal_difference_mean=0.XXX, median=0.XXX, p10=0.XXX, p90=0.XXX, nonzero_patients=0.XXX
- repeat_1_act_usage_cg: tree_split=...; path_traversal=...; terminal_difference_mean=..., median=..., p10=..., p90=..., nonzero_patients=...
- repeat_2_act_usage_cg: tree_split=...; path_traversal=...; terminal_difference_mean=..., median=..., p10=..., p90=..., nonzero_patients=...
- verdict: SEARCH_LEADER | DIAGNOSTIC_LEADER | NOT_LEADER | INELIGIBLE
- lesson: one sentence
```

All three ACT-use lines are mandatory even for ineligible runs.

## Completion

After 20 slots—or earlier if no credible hypothesis remains—freeze only the
eligible run named in `best_run.txt`. If that file is empty, report
`NO_ELIGIBLE_CANDIDATE` and preserve the test. `diagnostic_leader.txt` is never
authorization for test evaluation.
