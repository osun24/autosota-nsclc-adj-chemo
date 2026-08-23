# Autonomous program: Reactome-wide RSF permutation importance

## Goal

Optimize an RSF over all development-available Reactome genes plus the locked
clinical covariates. The primary objective is five-year IPCW-AIPW
policy-versus-anti-policy RMST difference; Harrell C-index is the strict
secondary tie-breaker. Then compute held-out permutation importance and search
every top-N gene panel from N=0 through N=32 to deliver the final covariates and
RSF hyperparameters.

Before the first experiment, read `README.md`, `red_lines.md`, all of `log.md`,
and `train.py`.

## Editable surface

- Edit only `train.py` during autonomous experiments.
- Append hypotheses and completed results to `log.md`.
- Never edit locked files, ledgers, run artifacts, or data.
- Do not place gene symbols or patient identifiers in `train.py`.

## Fixed resources

- 20 unique full experiments; four hours per experiment.
- Three fixed stratified CV folds.
- All eligible Reactome genes in the screening forest.
- Three held-out permutation repeats.
- Every N in 0..32 receives at least one reduced-stage evaluation.
- At most 1,000 trees per forest.

## Iteration workflow

1. Review `best_run.txt` and the last three completed attempts.
2. State one testable search-space or regularization hypothesis in `log.md`.
3. Audit the red lines before changing `train.py`.
4. Change only the candidate search specification in `train.py`.
5. Launch exactly one full run with `python -m reactome_rsf_pfi_arena.run`.
6. Record the screening optimum, final top N, final RMST difference, final
   C-index, selected genes, and runtime.
7. Treat RMST as primary. Compare C-index only when RMST scores are tied within
   the locked numerical tolerance.

## Result schema for `log.md`

```markdown
### iter_NNN — title
- type: PARAM | SEARCH
- hypothesis: one sentence
- changed: exact train.py change
- red_line_audit: concise audit
- run_id: run_NNN_...
- screening_rmst_difference: X.XXX months
- screening_cindex: 0.XXX
- final_top_n: N
- final_rmst_difference: X.XXX months
- final_cindex: 0.XXX
- final_genes: [...]
- verdict: SEARCH_LEADER | NOT_LEADER
- lesson: one sentence
```

## Completion

Freeze the lexicographic leader named in `best_run.txt`. The run directory's
`final_spec.json`, `final_covariates.txt`, `permutation_importance.csv`,
`permutation_importance.png`, and `top_n_results.csv` are the handoff. The
sealed test remains unavailable until a human performs a one-shot evaluation.

