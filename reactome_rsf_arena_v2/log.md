# Reactome causal-RMST RSF v2 experiment log

Fixed budget: 20 full experiments. The starting candidate carries forward the
legitimate pathway-level stability result from v1; it has not consumed a v2
experiment. Smoke metrics are infrastructure-only.

Read `v1_findings.md` before proposing the first hypothesis. Every adaptive
probe, including train-only structural or covariate probes, must consume a full
experiment unless it is part of the locked smoke check.

Append each prespecified hypothesis before editing `train.py`, then complete
the same entry after the launcher returns using the schema in `program.md`.

## Infrastructure notes (non-budgeted, synthetic-data timings only)

Wall-time feasibility was established with `--smoke` (22 s) plus library
benchmarks on **synthetic** matrices of the same shape as the arena workload.
No development row, outcome, treatment, or expression value was inspected, and
no gene identity or metric was produced, so these are launcher-feasibility
facts rather than adaptive probes.

- A 1,000-tree `RandomSurvivalForest` fit is dominated by the number of
  *continuous* split candidates, not the raw feature count. On a 775-row
  fitting partition: 19 clinical features (18 binary dummies + Age) ~6 s;
  +2 continuous modules ~15 s; +8 raw continuous genes ~46 s.
- A full run fits 48 forests (2 repeats x 4 folds x 2 models x 3 seeds), of
  which 24 are the locked 1,000-tree clinical comparator (~2.4 min total).
- The stability selector runs 9 x 40 x 2 = 720 `_dr_ranking` calls at ~0.34 s
  each (~4-5 min).
- Projected totals: `module_count=2` ~15 min; `module_count=1` ~13 min;
  `module_count=4` ~19 min; **`representation="raw"` with 8 genes at 1,000
  trees ~26 min, i.e. over the 25-minute wall.** Raw-representation controls
  must reduce `n_estimators` to roughly 600 or fewer.

### iter_001 — locked starting candidate, full-budget baseline
- type: PARAM
- hypothesis: The v1-derived eight-gene / two-module starting candidate, evaluated once at full budget (2 repeats x 4 folds, 3 seeds, 1,000 trees, 40 stability draws, 4,000 bootstraps), establishes the reference increment, bootstrap dispersion, and ACT-use profile that every later v2 hypothesis is measured against; it is expected to be ineligible on `all_repeat_genomic_increment_positive`.
- changed: none — `train.py` is run exactly as shipped. This is the reference point, not a tuning attempt, and it consumes slot 1 of 20.
- red_line_audit: No file edited. Launcher-only execution via `python -m reactome_rsf_arena_v2.run`. Development CSVs only; no test path is opened, named, or hashed. Estimand, gates, bootstrap, and comparator geometry untouched. Timings above came from synthetic matrices, not from development data.
- run_id: run_001_20260822T010814Z
- eligible: false
- failed_gates: [all_repeat_genomic_increment_positive, all_repeat_genomic_value_at_least_clinical]
- reward: -1000000.000
- diagnostic_score: -8.278
- repeat_1_increment: -1.581 months
- repeat_1_lcb: -5.266 months
- repeat_2_increment: -1.858 months
- repeat_2_lcb: -8.000 months
- repeat_range: 0.277 months
- development_alignment_cg: 2.212 months
- development_cindex_cg: 0.677
- gene_jaccard: 0.131
- development_act_usage_cg: tree_split=0.502; path_traversal=0.116; terminal_difference_mean=0.116, median=0.095, p10=0.024, p90=0.236, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.489; path_traversal=0.114; terminal_difference_mean=0.114, median=0.095, p10=0.022, p90=0.230, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.514; path_traversal=0.118; terminal_difference_mean=0.118, median=0.095, p10=0.027, p90=0.242, nonzero_patients=1.000
- verdict: INELIGIBLE
- lesson: Module compression did not solve v1's crowding-out -- adding only two continuous module features halves the forest's operational use of ACT (tree_split 0.823 -> 0.502, terminal_difference 0.239 -> 0.116) and with it the counterfactual benefit scale (median |b| 1.427 -> 0.599), which fully accounts for the -1.72-month increment.

#### Reference numbers from iter_001 (baseline for all later runs)
- Locked clinical comparator: alignment 3.931, value 48.139, C-index 0.681,
  ACT recommended 0.373, median |benefit| 1.427, ACT tree_split 0.823,
  terminal_difference_mean 0.239.
- Clinical+genomic: alignment 2.212, value 47.602, C-index 0.677,
  ACT recommended 0.338, median |benefit| 0.599.
- Constant policies: all-ACT 47.389, all-observation 45.604 (ATE ~= +1.79).
- Wall clock 649.6 s of the 1,500 s limit.

#### Two mechanisms this baseline isolates
1. **Continuous-feature split bias.** ACT is binary and offers the log-rank
   splitter one candidate threshold per node; each standardized module score
   offers ~n. With `max_features=1.0` both modules are evaluated at every node,
   so they outcompete ACT nearly everywhere. The clinical forest carries one
   continuous feature (Age); the genomic forest carries three. ACT usage and
   benefit scale both fall by ~55%, and the genomic C-index does *not* rise
   (0.677 vs 0.681), so the modules are paying the full crowding cost while
   adding almost no information.
2. **Possible sign cancellation inside a module.** `prepare._gene_effect_scores`
   returns the *absolute* partial correlation with the DR benefit pseudo-outcome,
   and `FeatureTransformer` averages standardized members. Genes whose
   benefit-association runs in opposite directions therefore cancel inside a
   module. Module composition is agent-controllable, because
   `np.array_split` slices the selector's returned list in order.

#### Height of the eligibility bar (derived from run 001, not a probe)
`selection_lcb` is the 0.25th bootstrap percentile, i.e. roughly
`mean - 2.81 sd`. Observed sd of the paired increment was 1.33 and 2.04 months,
so a positive worst-repeat LCB needs a mean increment of roughly 3.7-5.7
months against a clinical alignment of only 3.93. Writing the increment as
`(2/n) * sum over policy-disagreeing patients of +/- gamma`, its bootstrap sd
scales as `sqrt(f) * sd(gamma)` while its mean scales as `f * m`, so the
signal-to-noise ratio grows as `sqrt(f) * m`: broad-but-correct deviation from
the clinical policy is favoured over a narrow one. With sd(gamma) implied to be
~40 months by the observed dispersion, eligibility needs `sqrt(f) * m` of order
7 months. This is a demanding bar and may well not be reachable; that outcome
would be a legitimate negative result, not a reason to weaken any gate.

### iter_002 — sign-coherent module composition
- type: ALGO
- hypothesis: The eight selected genes are currently grouped into modules by pathway-and-spread order while their association with the DR benefit pseudo-outcome is scored by *absolute* partial correlation, so oppositely-associated genes cancel inside a standardized mean; ordering the same selected set by *signed* benefit association before `np.array_split` makes module 1 a coherent benefit-positive score and module 2 a coherent benefit-negative score, which should raise the genomic counterfactual benefit scale and the increment without changing which genes are selected.
- changed: Added `_signed_gene_scores()` to `train.py`, an exact sign-preserving replica of `prepare._gene_effect_scores` (same standardized-clinical design, same Ridge(alpha=1.0) gamma residualization, same chunked lstsq, same NaN-to-column-median rule). `_dr_ranking` now returns `np.abs()` of it, so pathway stability ranking, the tie-break and the fallback fill are bit-identical to iter_001 and the selected gene *set* is unchanged. A new `signed_total` accumulator averages the signed score over the same 40 x 2 half-samples, and `stability_select_genes` returns the selected genes sorted by descending mean signed association. `CANDIDATE["name"]` -> `v2_pathway8_signed_modules`. No other field changed; `module_count` stays 2.
- red_line_audit: Only `train.py` edited. Signed scores are computed exclusively inside the fitting partition passed to the selector, from the locked cross-fitted train-only benefit pseudo-outcome — no assessment rows, outcomes, or covariate summaries. No hard-coded gene symbols, patient indices, or assessment rankings. The detectability filter is retained per v1's mandate. No gene-by-ACT product is formed; only the order of a gene list changes. Estimand, gates, bootstrap, thresholds, and the locked clinical comparator are untouched. Development CSVs only.
- run_id: run_002_20260822T012403Z
- eligible: false
- failed_gates: [all_repeat_genomic_increment_positive, all_repeat_genomic_value_at_least_clinical, all_repeat_genomic_value_at_least_best_constant]
- reward: -1000000.000
- diagnostic_score: -7.815
- repeat_1_increment: -1.731 months
- repeat_1_lcb: -7.096 months
- repeat_2_increment: -2.224 months
- repeat_2_lcb: -7.322 months
- repeat_range: 0.493 months
- development_alignment_cg: 1.954 months
- development_cindex_cg: 0.670
- gene_jaccard: 0.131
- development_act_usage_cg: tree_split=0.517; path_traversal=0.119; terminal_difference_mean=0.119, median=0.094, p10=0.027, p90=0.246, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.500; path_traversal=0.117; terminal_difference_mean=0.117, median=0.091, p10=0.024, p90=0.234, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.534; path_traversal=0.120; terminal_difference_mean=0.120, median=0.098, p10=0.030, p90=0.258, nonzero_patients=1.000
- verdict: INELIGIBLE
- lesson: Sign-coherent grouping is not the binding defect -- it swapped exactly one gene between modules at the full-development selection, moved the increment by -0.26 months (to -1.977), cost 0.007 C-index, widened the repeat range to 0.493, and left ACT usage untouched at tree_split 0.517, so pathway coherence is the more valuable grouping principle and composition tweaks at this granularity are inside the noise.

#### What iter_002 rules out
- The signed statistic is exact (`abs(signed)` reproduced the locked ranking to
  0.0 error on a synthetic check), so the gene *set* and pathway ranking were
  identical to iter_001 by construction. Only grouping changed.
- Grouping by sign moved just one gene across the module boundary:
  iter_001 gave `[CHRDL1, BMP2, TGFBR3, SMURF2]` + `[TBL1XR1, NRIP1, CHD9, RORA]`
  (a TGF-beta/BMP block and a nuclear-receptor-coregulator block); iter_002 gave
  `[SMURF2, CHRDL1, BMP2, RORA]` + `[NRIP1, TBL1XR1, CHD9, TGFBR3]`, which mixes
  the two pathway blocks. The result was *worse* on increment, C-index and
  repeat range. This reinforces v1's finding that the **pathway is the stable
  unit** and shows sign cancellation was not what was limiting the candidate.
- A one-gene composition change producing a 0.26-month increment swing, against
  a repeat range of 0.28-0.49 and a bootstrap sd of 1.7-1.9, means module
  composition at this granularity is within noise. Only structural changes
  (feature count, forest geometry, policy threshold) can plausibly move the
  increment by the several months eligibility requires.

#### Standing obstacle after two slots
ACT operational use in the genomic forest is pinned near tree_split 0.50-0.52
and terminal_difference 0.116-0.119, against the locked comparator's 0.82 and
0.239, in both runs and both repeats. Neither run has yet touched the feature
geometry that causes it. That is the next target.
