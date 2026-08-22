# Findings inherited from v1

V1 exhausted 20 slots: 19 completed candidates and one infrastructure failure.
No candidate was eligible. This is a negative result for the v1 S-learner RSF
and search space, not proof that genomic treatment heterogeneity is absent.

## Findings to build on

- Raw gene-block width crowded ACT out of survival-tree splits. ACT tree use
  fell from roughly 0.34 with four genes to 0.09 with twenty genes.
- A reproducible eight-gene signal was obtainable at the Reactome pathway
  level by selecting pathways with train-only doubly robust scores and ordering
  detectable members with a stable covariate summary.
- Winsorized per-gene stability alone selected near-floor probes. Detectability
  filtering is mandatory; do not reintroduce those artifacts.
- The genomic counterfactual benefit scale was narrower than the clinical
  scale. A shared 0.25-month threshold often silenced the genomic policy; zero
  threshold restored evaluability but did not make v1 eligible.
- Every v1 candidate had a negative train-OOF genomic policy increment and a
  negative multiplicity-adjusted lower bound. V1 run 020 was diagnostically
  interesting, not validated and not eligible.

## V2 design response

- The former train and validation cohorts are pooled as one adaptive
  development cohort. Neither retains confirmatory status.
- Two repeated four-fold outer cross-fits replace the train-versus-validation
  comparison. The untouched test remains human-only.
- Selected genes may enter as fold-fitted module scores, so eight identified
  genes need not compete with ACT as eight separate split variables.
- The clinical comparator geometry is locked independently. The agent may tune
  only the genomic forest and selector, preventing comparator weakening.
- Eligibility and promotion require a positive, repeat-consistent genomic
  increment. Sentinel ties never nominate an ineligible run.

Do not repeat v1's depth-only, raw-width, shared-threshold, unfiltered
winsorization, or pathway-smoothing experiments without a materially new
mechanism and a prespecified reason.
