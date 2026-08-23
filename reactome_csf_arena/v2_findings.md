# Findings inherited from the completed RSF arenas

The predecessor searches are immutable evidence, not experiment slots in this
arena. V1 found no eligible candidate. V2 completed 20 slots and froze its only
eligible candidate, but its worst-repeat selection-adjusted LCB remained
negative. This is motivation for a new treatment-effect architecture, not
evidence that the V2 diagnostic leader is test-eligible.

## Reuse

- Detectability filtering is mandatory. Unfiltered stability selection favored
  near-floor probes whose apparent reproducibility came from a few recurring
  samples.
- The most reproducible selection unit was a Reactome pathway scored with a
  fit-only doubly robust 60-month benefit pseudo-outcome.
- Eight genes was a useful starting width. Exact gene identity stability and
  policy-value stability were distinct, so retain both diagnostics.
- Winsorization controlled the heavy-tailed pseudo-outcome without dropping
  censored patients. Retain the non-finite guard added after the V2 Cox failure.

## Do not transplant mechanically

- RSF ACT crowd-out was architectural: ACT competed with continuous genes for
  ordinary survival splits. In this arena ACT is the separate treatment `W`,
  so split/path/terminal ACT diagnostics and `max_features` remedies are
  inapplicable.
- One-module compression helped ACT remain visible to the S-learner but may
  cancel gene-specific effect modification. Start with eight raw genes; compare
  modules as a controlled representation experiment.
- V2's genomic-only threshold tuning is retired. The CSF policy threshold is
  locked at zero because the estimand contains no treatment-cost term.
- V2 allowed different genomic forest geometry from the locked comparator.
  Here one candidate geometry is applied to both clinical and C+G, isolating
  genomic information more cleanly.

## Starting scientific question

Does changing only the learner—from a prognostic S-learner to an honest RMST
treatment-effect forest—turn the inherited stable Reactome selector into a
repeat-consistent positive genomic policy increment? Establish that baseline
before inventing a new selector.
