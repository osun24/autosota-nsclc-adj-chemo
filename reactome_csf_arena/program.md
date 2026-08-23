# Autonomous program: Reactome GRF causal survival forest

## Goal

Maximize the locked eligible reward from `python -m reactome_csf_arena.run`.
Find a small, stable Reactome gene panel that improves a 60-month ACT policy
over the matched clinical-only causal survival forest and remains consistent
across repeated outer cross-fits.

Before experiment 1, read `README.md`, `red_lines.md`, `v2_findings.md`, all of
`log.md`, `budget.json`, and `train.py`.

## Editable surface

- Edit only `train.py`.
- Append prespecified hypotheses and completed results to `log.md`.
- Never edit locked R/Python evaluation code, budgets, manifests, ledgers,
  artifacts, or data.
- Selectors receive only the current fitting partition, Reactome membership,
  available development genes, and the validated candidate specification.
- Do not hard-code gene symbols.

## Fixed resources

- 20 unique full experiments; 30 minutes per experiment.
- Two repeated four-fold outer cross-fits and three inner nuisance folds.
- Seeds 42, 43, and 44; at most 1,000 trees in every forest.
- At most 16 selected genes and 34 treatment-effect X features.
- 4,000 paired patient bootstraps per repeat.
- Honest forests, all-feature `mtry`, zero-month policy threshold, and identical
  geometry for clinical and C+G are locked.
- Internal GRF tuning is disabled and cannot create unmetered searches.
- Budget accounting includes GRF's two internal nuisance survival forests per
  causal forest; the supplied cross-fitted `W.hat` suppresses a third nuisance
  forest for propensity.

Use only the launcher. Scientific probes must consume a full slot. Smoke is
SHA-locked to the shipped baseline and cannot be run after candidate edits.

## Iteration workflow

1. Review `best_run.txt`, `diagnostic_leader.txt`, and the latest attempts.
2. Prespecify one interpretable hypothesis and exact change in `log.md`.
3. Audit every red line before editing.
4. Edit `train.py` and launch exactly one full run.
5. Record reward, failed gates, both repeat increments and LCBs, repeat range,
   source gap, alignment, C-index, gene Jaccard, and all mandatory CSF
   diagnostics below.
6. Promote only an eligible reward improvement. An ineligible diagnostic
   leader can guide the next hypothesis but is never a test nominee.

After three parameter-only attempts, change selector or representation logic.
Start from the inherited pathway-stability selector with eight raw genes.
Compare raw genes against fit-only modules before adding selector complexity.
Do not add policy trees or SHAP selection unless a prespecified result identifies
a specific failure they can address.

## Mandatory per-iteration diagnostics

```text
act_mechanism: W supplied separately; ACT absent from X
rsf_act_split/path/terminal: NA by design
development_csf_cg: seed_agreement=...; seed_tau_correlation=...;
  genomic_vimp_fraction=...; benefit_iqr=...; median_abs_benefit=...;
  nontrivial_fraction=...; act_recommended=...
repeat_1_csf_cg: same fields
repeat_2_csf_cg: same fields
```

Also record the clinical versions. These replace—not silently omit—the prior
RSF ACT-use diagnostics, which are undefined because ACT is `W`, not an X
feature.

## Completion

Freeze only the eligible run named in `best_run.txt`. If none exists, report
`NO_ELIGIBLE_CANDIDATE` and preserve the test. Before test access, verify the
human finalizer on synthetic data. The held-out test is evaluated once at the
very end; no model, threshold, selector, or code changes may follow.
