# Red lines

1. **Test is human-only.** Search code may load only `affyfRMATrain.csv` and
   `affyfRMAValidation.csv`. It must never open, enumerate, hash, summarize, or
   infer from a test dataset or a test-derived artifact.
2. **Only `train.py` is agent-editable Python during a search.** Append-only
   notes may be added to `log.md`; locked preparation, metrics, budgets,
   launcher, manifest, ledger, and prior artifacts do not change.
3. **All adaptive transforms are fit-only.** Imputation, RSF fitting, and
   nuisance models are fit on the fitting rows. Permutation importance is
   computed only on held-out fold rows.
4. **The screening forest contains every eligible Reactome gene.** A smoke run
   may use its declared reduced gene pool, but a full run may not prefilter,
   rank, or select genes before the all-gene forest.
5. **Never drop censored patients.** IPCW-AIPW scoring includes every row.
6. **The estimand is fixed.** The primary metric is the cross-fitted 60-month
   IPCW-AIPW policy-versus-anti-policy RMST difference. Harrell C-index is the
   secondary, lexicographic objective; it is not combined with RMST using an
   arbitrary weighted sum.
7. **Recommendations use paired counterfactual RMST.** For every patient,
   predict with ACT set to zero and one while holding all pretreatment
   covariates fixed. Numerical ties recommend observation.
8. **Permutation preserves counterfactual pairs.** A feature receives one
   patient permutation per repeat and that same permutation is used in its ACT
   zero and ACT one copies. Outcomes, treatment, and nuisance scores are never
   permuted.
9. **Top-N is honest about zero.** Clinical-only (`N=0`) and every integer gene
   count through the configured maximum (default 32) must receive at least one
   reduced-stage CV evaluation before a final N is chosen.
10. **No metric shopping or hard-coded biology.** Genes are ordered only by
    held-out permutation RMST importance, with C-index importance and symbol as
    deterministic tie-breakers. No symbols, row IDs, predictions, or outcomes
    may be hard-coded in `train.py`.
11. **Budget is final.** At most 20 unique full attempts, four hours per attempt,
    and 1,000 trees per forest. Failed and timed-out full attempts consume a
    slot; identical `train.py` snapshots cannot be rerun.
12. **Claims remain observational.** Internal cross-validation and permutation
    importance do not establish causal validity, external transportability, or
    a deployable treatment rule.

