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

### iter_003 — one module: halve the continuous split-competitors
- type: PARAM
- hypothesis: ACT is binary and offers the log-rank splitter one candidate threshold per node while each standardized module score offers ~n, so with `max_features=1.0` the modules outcompete ACT at almost every node; the locked comparator carries one continuous feature (Age) and reaches ACT tree_split 0.823, whereas three continuous features (Age + two modules) pin it at 0.502-0.517. Dropping to a single module leaves two continuous features and should recover ACT usage toward tree_split ~0.65 and terminal_difference ~0.16, raising the genomic benefit scale and the increment.
- changed: `train.py` restored to its iter_001 selector exactly (from `runs/run_001_.../train_snapshot.py`, discarding iter_002's sign ordering, which lost), then `CANDIDATE["module_count"]` 2 -> 1 and `CANDIDATE["name"]` -> `v2_pathway8_one_module`. This is one *effective* change: a module is the mean of its standardized members and `np.array_split(arange(8), 1)` puts all eight genes in a single group, so gene ordering is provably inert at `module_count=1` and reverting iter_002 cannot affect this run's result. `n_genes` stays 8, threshold stays 0.00, and the genomic RSF geometry is untouched.
- red_line_audit: Only `train.py` edited. No selector logic, detectability filter, estimand, gate, bootstrap, or threshold changed; the locked clinical comparator is untouched and is not weakened by reducing the *genomic* feature count. No hard-coded gene symbols and no gene-by-ACT product. `module_count=1` is inside the locked `1 <= module_count <= min(4, n_genes)` validation and lowers total features to 20, well under the 35 cap. Development CSVs only; projected wall clock ~13 min.
- run_id: run_003_20260822T013504Z
- eligible: false
- failed_gates: [all_repeat_genomic_increment_positive]
- reward: -1000000.000
- diagnostic_score: -5.259
- repeat_1_increment: -0.617 months
- repeat_1_lcb: -4.596 months
- repeat_2_increment: -0.105 months
- repeat_2_lcb: -4.746 months
- repeat_range: 0.512 months
- development_alignment_cg: 3.570 months
- development_cindex_cg: 0.677
- gene_jaccard: 0.131
- development_act_usage_cg: tree_split=0.596; path_traversal=0.142; terminal_difference_mean=0.142, median=0.123, p10=0.037, p90=0.271, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.612; path_traversal=0.148; terminal_difference_mean=0.148, median=0.132, p10=0.037, p90=0.271, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.580; path_traversal=0.137; terminal_difference_mean=0.137, median=0.114, p10=0.038, p90=0.272, nonzero_patients=1.000
- verdict: DIAGNOSTIC_LEADER
- lesson: Removing one continuous module feature recovered 1.62 months of increment (-1.977 -> -0.361) and cleared both value gates, confirming that continuous-feature split competition against binary ACT -- not gene choice or module composition -- is the dominant obstacle in this arena.

#### The crowding curve is now measured
ACT operational use in a 1,000-tree forest at the locked geometry, as a
function of how many *continuous* features share the split search:

| continuous features | model | ACT tree_split | terminal_difference_mean | median abs benefit | alignment |
|---|---|---|---|---|---|
| 1 (Age) | locked clinical | 0.823 | 0.239 | 1.427 | 3.931 |
| 2 (Age + 1 module) | iter_003 | 0.596 | 0.142 | 0.855 | 3.570 |
| 3 (Age + 2 modules) | iter_001 | 0.502 | 0.116 | 0.599 | 2.212 |

Each continuous competitor is worth roughly 0.10-0.12 of ACT tree_split and a
large slice of the counterfactual benefit scale. `module_count=1` is the floor
for this lever, so further crowding relief must come from reducing per-node
competition rather than from removing features.

#### Gate status is now down to one
Both value gates cleared: genomic value 48.242 vs clinical 48.139 and vs the
best constant policy 47.337, in both repeats. C-index drop is only 0.004 of the
allowed 0.03 and seed agreement is 0.969 against a 0.85 floor, so there is real
headroom to spend on forest regularization. Repeat 2's increment reached -0.105,
i.e. near parity, while repeat 1 sat at -0.617; the repeat range of 0.512 is now
a material part of the reward and will need attention once the increment turns.

#### Note on max_depth (rules out a v1-style depth experiment)
With `min_samples_leaf=8` and `min_samples_split=16` on a 775-row fitting
partition, a tree can hold at most ~97 leaves, so the realised depth is ~6.6 and
the depth cap of 9 is essentially never binding. Raising `max_depth` toward 12
therefore cannot change anything, which is consistent with v1's refuted
depth-only experiments. Depth is not a live lever here and will not be spent on.

### iter_004 — genomic mtry 0.35: let ACT be evaluated without continuous rivals
- type: PARAM
- hypothesis: `module_count=1` is the floor for removing continuous competitors, so the remaining crowding lever is per-node competition. With `max_features=1.0` ACT is offered at every node but *always* alongside Age and the module; lowering the genomic forest's `max_features` to 0.35 (7 of 20 features) maximizes the chance that ACT is drawn with no continuous rival in the same node, which should lift ACT tree_split from 0.596 toward the comparator's 0.823 and push the increment from -0.361 toward parity or above.
- changed: `CANDIDATE["rsf"]["max_features"]` 1.00 -> 0.35 and `CANDIDATE["name"]` -> `v2_pathway8_one_module_mtry035`. Nothing else. Choice of 0.35 is combinatorial, not fitted: with 20 genomic features of which 2 are continuous, P(ACT drawn) x P(no continuous among the other k-1) = (k/20) x C(17,k-1)/C(19,k-1) peaks at k=7, i.e. `max_features=0.35` (0.160), against 0.154 at k=5, 0.154 at k=8, 0.132 at k=10 and 0 at k=20.
- red_line_audit: Only `train.py` edited, and only the *genomic* forest. The locked clinical comparator keeps `max_features=1.0` from `budget.json`; nothing here can weaken it, and the genomic forest is explicitly agent-tunable under red line 9. Value is within the locked `0.20 <= max_features <= 1.0` validation. Estimand, gates, bootstrap, threshold, selector, and detectability filter untouched. No hard-coded gene symbols, no gene-by-ACT product. Development CSVs only. Runtime should fall, since each node evaluates 7 features instead of 20.
- watch: C-index drop (0.004 of an allowed 0.03) and seed agreement (0.969 against a 0.85 floor) are the two gates most exposed to added split randomness; both value gates cleared in iter_003 by slim margins and could regress.
- run_id: run_004_20260822T014430Z
- eligible: false
- failed_gates: [all_repeat_genomic_increment_positive]
- reward: -1000000.000
- diagnostic_score: -5.451
- repeat_1_increment: -0.259 months
- repeat_1_lcb: -4.563 months
- repeat_2_increment: +0.150 months
- repeat_2_lcb: -5.042 months
- repeat_range: 0.409 months
- development_alignment_cg: 3.877 months
- development_cindex_cg: 0.687
- gene_jaccard: 0.131
- development_act_usage_cg: tree_split=0.688; path_traversal=0.173; terminal_difference_mean=0.173, median=0.162, p10=0.098, p90=0.260, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.681; path_traversal=0.168; terminal_difference_mean=0.168, median=0.161, p10=0.096, p90=0.252, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.696; path_traversal=0.177; terminal_difference_mean=0.177, median=0.164, p10=0.100, p90=0.269, nonzero_patients=1.000
- verdict: NOT_LEADER
- lesson: Genomic mtry 0.35 lifted ACT tree_split 0.596 -> 0.688 and turned repeat 2's increment positive (+0.150) while *improving* C-index to 0.687 above the comparator's 0.681, leaving repeat 1 (-0.259) as the only barrier to eligibility.

#### Why this is NOT_LEADER despite the best increment so far
`diagnostic_leader.txt` still names run_003. The continuous diagnostic score is
`min_repeat_lcb - repeat_range`, and although iter_004's mean increment is far
better (-0.054 vs -0.361), repeat 2's bootstrap sd rose to 1.643, dragging its
0.25th-percentile LCB to -5.042 and the score to -5.451 against iter_003's
-5.259. This is exactly the metric-shopping trap red line 12 warns about: the
diagnostic ranking and the eligibility target disagree here, and the
**eligibility target is the one that counts**.

#### Correction to the iter_001 note on the eligibility bar
`prepare.evaluate_candidate` sets `eligible = all(gates.values())` and only then
computes `reward = robust_lcb - increment_range`. The bootstrap LCB is
therefore **not a gate**. The iter_001 estimate that eligibility needs a 3.7-5.7
month increment applies to a *positive reward*, not to eligibility. Eligibility
needs only both repeat increments strictly positive plus the nine other gates,
all of which iter_004 already passes. An eligible run with a negative reward is
still eligible and still promoted to `best_run.txt`. The search target is
therefore much closer than iter_001 implied: flip repeat 1 by ~0.26 months
without breaking the other nine gates.

#### Gate margins going into iter_005 (all currently passing except one)
- genomic value 48.397 vs clinical 48.139 and vs best constant 47.315 -- clear.
- C-index 0.687 vs clinical 0.681: the genomic model is now *ahead*, so the
  0.03-drop gate has ~0.036 of slack.
- seed agreement 0.976 vs the 0.85 floor -- clear.
- nontrivial benefit fraction 0.941 vs the 0.10 floor -- clear.
- overlap 0.871 vs 0.80, and IPTW ESS 362 vs the 310 floor -- clear.
- **gene Jaccard 0.131 vs the 0.10 floor is the fragile one.** Any change that
  reduces selection agreement across the eight outer folds can fail this gate
  outright, which rules out dropping to `n_genes=4` (one pathway block per fold
  would likely share nothing between folds).

#### mtry is now at its argued optimum -- do not push it lower
The k=7 choice came from maximising P(ACT drawn and facing no continuous rival)
= (k/20) x C(17,k-1)/C(19,k-1). That expression is 0.1597 at k=7, 0.1596 at k=6,
0.1535 at k=5 and 0.1404 at k=4, so `max_features` of 0.30-0.35 is the plateau
and 0.20-0.25 is predicted to be *worse*, not better. Spending a slot on
`max_features=0.20` would contradict the mechanism that motivated 0.35, so the
crowding lever is considered exhausted.

### iter_005 — order pathway members by benefit association, not by spread alone
- type: ALGO
- hypothesis: The DR benefit signal currently decides only *which pathway* wins; within that pathway the four members are taken in order of expression spread, so the module can be built from the widest-spread but least benefit-associated genes. Ordering the already-detectability-filtered members by their aggregated DR benefit association should make the single module genuinely benefit-predictive and lift repeat 1's increment above zero, which is the only gate still failing.
- changed: In `stability_select_genes`, `ordered_members` is now sorted by `(-mean_score[item], -gene_spread[available[item]], available[item])` instead of `(-gene_spread[available[item]], available[item])`. `mean_score` is the existing aggregate of the locked clinical-adjusted DR statistic averaged over the same 40 x 2 = 80 fitting-partition half-samples, so it is a stable fit-only quantity and no new computation is added; expression spread is retained as the tie-break. `CANDIDATE["name"]` -> `v2_pathway8_one_module_mtry035_drorder`. Nothing else changes: pathway ranking, `n_genes=8`, `module_count=1`, `max_features=0.35`, threshold 0.00 all identical to iter_004.
- red_line_audit: Only `train.py` edited. This does **not** reintroduce v1's refuted unfiltered-winsorization artifact: `_detectable_genes` still runs first, so only genes in the upper half of within-partition expression spread are ever eligible, and DR association merely orders that filtered pool. `mean_score` is computed exclusively from fitting-partition half-samples via the locked cross-fitted train-only benefit pseudo-outcome; no assessment rows, outcomes, or covariate summaries are touched. No hard-coded gene symbols, patient indices, or assessment rankings. No gene-by-ACT product. Estimand, gates, bootstrap, threshold, and the locked clinical comparator untouched. Development CSVs only.
- watch: gene Jaccard is the fragile gate at 0.131 against a 0.10 floor; changing within-pathway member choice could reduce agreement across the eight outer folds and fail eligibility on stability even if the increment improves.
- run_id: run_005_20260822T015532Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10]
- reward: -1000000.000
- diagnostic_score: -5.412
- repeat_1_increment: +0.189 months
- repeat_1_lcb: -3.751 months
- repeat_2_increment: +1.850 months
- repeat_2_lcb: -1.989 months
- repeat_range: 1.662 months
- development_alignment_cg: 4.951 months
- development_cindex_cg: 0.689
- gene_jaccard: 0.075
- development_act_usage_cg: tree_split=0.706; path_traversal=0.173; terminal_difference_mean=0.173, median=0.161, p10=0.091, p90=0.271, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.713; path_traversal=0.177; terminal_difference_mean=0.177, median=0.165, p10=0.091, p90=0.271, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.700; path_traversal=0.170; terminal_difference_mean=0.170, median=0.157, p10=0.090, p90=0.271, nonzero_patients=1.000
- verdict: NOT_LEADER
- lesson: Using the DR benefit signal to order members *within* the winning pathway turned both repeat increments positive for the first time (+0.189 and +1.850, alignment 4.951 vs the comparator's 3.931) but destabilised within-pathway member choice enough to drop gene Jaccard to 0.075 and fail the stability gate, trading the increment gate for the stability gate.

#### The binding constraint has moved
Nine of ten gates now pass, including the increment gate that blocked
iter_001-004. The failure is `gene_selection_jaccard_at_least_0_10` at 0.075.
This is a genuine trade, not a technicality: the selector became more responsive
to the fitting partition's DR signal, which is exactly what made the module
benefit-predictive *and* exactly what made it fold-dependent.

#### Where the instability actually lives
Pathway-level agreement survived; member-level agreement did not. The
cross-fold frequency table still shows a coherent branched-chain-amino-acid
block (ACADSB 4, HIBCH 4, AUH 3) and a coherent TGF-beta/SMAD block (SMAD5 3,
SMAD7 3, TGFBR3 3), but the *members* drawn from the branched-chain pathway
changed wholesale from iter_004's (ACAD8, ALDH6A1, BCAT1, MCCC1) to
(ACADSB, HIBCH, AUH, ...). So the fix must stabilise within-pathway member
ordering, not pathway ranking.

Mean |DR| per gene is an average of a heavy-tailed statistic and is evidently
still fold-unstable even after 80 half-samples. A rank-based co-selection
*frequency* is the standard remedy and is far more reproducible than the mean of
a heavy-tailed score. Note that `STABILITY_TOP_K = 300` is already declared in
`train.py` and currently unused -- it was clearly intended for exactly this.

#### Arithmetic of the Jaccard gate (sets of size n_genes)
For two selections of size 8 sharing `s` genes, Jaccard is `s/(16-s)`, so the
0.10 floor needs `s >= 1.45`, i.e. at least 2 shared genes on average across the
28 fold pairs. The observed 0.075 corresponds to `s ~= 1.11`. Raising `n_genes`
is a mechanical way to raise Jaccard at zero crowding cost -- under module
representation the forest still sees one feature regardless of `n_genes` -- but
it only helps if members are stable enough that sharing a pathway implies
sharing genes, which is precisely what iter_005 broke. Stabilise members first.

### iter_006 — rank-based co-selection frequency for within-pathway member choice
- type: ALGO
- hypothesis: iter_005 ordered pathway members by the *mean* of a heavy-tailed DR statistic, which is fold-unstable and dropped Jaccard to 0.075; replacing it with a rank-based co-selection *frequency* -- how often each gene lands in the global top-`STABILITY_TOP_K` by DR association across the same 80 fitting-partition half-samples -- keeps the benefit information that turned both increments positive while restoring cross-fold member agreement above the 0.10 Jaccard floor.
- changed: In `stability_select_genes`, added a `gene_top_counts` accumulator; in each half-sample, `np.argpartition` marks the `STABILITY_TOP_K = 300` genes with the largest DR association and increments their counts. `ordered_members` now sorts by `(-gene_top_counts[item], -mean_score[item], -gene_spread[available[item]], available[item])`, so a bounded rank frequency over 80 half-samples leads, iter_005's mean score breaks ties, and expression spread breaks any remaining tie. `CANDIDATE["name"]` -> `v2_pathway8_one_module_mtry035_freqorder`. `STABILITY_TOP_K` was already declared in `train.py` and previously unused; this is its intended purpose. Nothing else changes: pathway ranking, `n_genes=8`, `module_count=1`, `max_features=0.35`, threshold 0.00 all identical to iter_005.
- red_line_audit: Only `train.py` edited. Detectability filtering still runs first, so this is not v1's refuted per-gene-stability-without-detectability artifact -- the frequency only orders an already spread-filtered pool. All counts come from fitting-partition half-samples via the locked cross-fitted train-only benefit pseudo-outcome; no assessment rows, outcomes, or covariate summaries. No hard-coded gene symbols, patient indices, or assessment rankings. No gene-by-ACT product. Estimand, gates, bootstrap, threshold, and the locked clinical comparator untouched. Development CSVs only.
- watch: if `gene_top_counts` is zero for most members of the winning pathways the ordering degrades gracefully to iter_005's mean-score behaviour and Jaccard will stay near 0.075; the increments must also remain positive in both repeats, since iter_005's gain is the thing being protected.
- run_id: run_006_20260822T020549Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10]
- reward: -1000000.000
- diagnostic_score: -4.807
- repeat_1_increment: +0.473 months
- repeat_1_lcb: -3.455 months
- repeat_2_increment: +1.824 months
- repeat_2_lcb: -2.151 months
- repeat_range: 1.352 months
- development_alignment_cg: 5.079 months
- development_cindex_cg: 0.690
- gene_jaccard: 0.076
- development_act_usage_cg: tree_split=0.701; path_traversal=0.172; terminal_difference_mean=0.172, median=0.160, p10=0.093, p90=0.266, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.700; path_traversal=0.172; terminal_difference_mean=0.172, median=0.163, p10=0.092, p90=0.259, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.702; path_traversal=0.171; terminal_difference_mean=0.171, median=0.158, p10=0.094, p90=0.273, nonzero_patients=1.000
- verdict: DIAGNOSTIC_LEADER
- lesson: A bounded co-selection frequency improved the increments again (repeat 1 +0.189 -> +0.473, pooled alignment 5.079 vs the comparator's 3.931) but left Jaccard flat at 0.076, because it selected the identical eight genes -- confirming that *any* DR-driven member ordering is fold-dependent and that stability has to come from a fold-stable ordering key, not a better-behaved DR statistic.

#### What iter_006 establishes
The full-development module was byte-identical to iter_005
(`TGFBR3, SMAD7, BMPR1A, SMAD4, SIRT1, CHD9, EP300, TBL1XR1`), so the rank
frequency and the mean score agree on the winning pathway's best members. The
watch item's graceful-degradation case is what happened. The instability is not
an artifact of the heavy tail; it is intrinsic to ordering members by a
quantity that is re-estimated inside each fitting partition.

Contrast the two available ordering keys:
- **expression spread** is a covariate summary over hundreds of patients, nearly
  identical across fitting partitions -> Jaccard 0.131 (iter_004), but carries no
  benefit information, so repeat 1's increment stayed negative (-0.259).
- **DR association** carries the benefit information that made both increments
  positive, but is re-estimated per partition -> Jaccard 0.076.

#### Correction to a claim made in the iter_005 entry
The iter_005 note suggested raising `n_genes` is "a mechanical way to raise
Jaccard". That is wrong as stated: Jaccard is scale-invariant if the number of
shared genes grows in proportion to the set size. Doubling `n_genes` from 8 to
16 while holding `MAX_GENES_PER_PATHWAY=4` simply doubles the number of
contributing pathways and leaves the ratio near 0.076. `n_genes` only helps
through the *fraction of a pathway taken*: for two folds drawing `k` members
from a shared pathway with `M` detectable members, expected overlap is about
`k^2/M`, so raising `MAX_GENES_PER_PATHWAY` (not `n_genes` alone) is the lever
that raises the overlap fraction. Held in reserve.

#### Both remaining routes to eligibility are about equally far away
- iter_004 geometry passes Jaccard with 31% margin (0.131 vs 0.10) and fails the
  increment gate on repeat 1 by 0.26 months, with no obvious lever left.
- iter_006 geometry passes the increment gate with a large margin (+0.473,
  +1.824) and fails Jaccard by 32% (0.076 vs 0.10), with an obvious lever left.
The second is the one to push.

### iter_007 — DR-gated, spread-ordered members
- type: ALGO
- hypothesis: The two available ordering keys fail in complementary ways -- expression spread is fold-stable but benefit-blind (Jaccard 0.131, repeat 1 increment -0.259), DR association is benefit-bearing but fold-noisy (increments +0.473/+1.824, Jaccard 0.076). Using DR only for a *coarse* within-pathway gate and the fold-stable spread for the actual ordering should keep most of the benefit signal while restoring cross-fold member agreement above the 0.10 floor, because a median split is far more reproducible than a full ranking.
- changed: In `stability_select_genes`, members of a selected pathway are now partitioned by whether their co-selection frequency reaches that pathway's own median frequency; the qualifying half is ordered by expression spread and exhausted first, then the remainder, also spread-ordered. Concretely `ordered_members` becomes `sorted(gated, key=spread) + sorted(rest, key=spread)` where `gated = [m for m in members if gene_top_counts[m] >= np.median(gene_top_counts[members])]`. Using `>=` guarantees the gate keeps at least half of every pathway, so no pathway can be emptied. `CANDIDATE["name"]` -> `v2_pathway8_one_module_mtry035_gated`. Nothing else changes: pathway ranking, `n_genes=8`, `module_count=1`, `max_features=0.35`, threshold 0.00 identical to iter_006.
- red_line_audit: Only `train.py` edited. Detectability filtering still runs first, so the spread-ordered pool is unchanged from v1's mandate and this is not the refuted unfiltered-stability artifact. The gate uses only `gene_top_counts`, accumulated from fitting-partition half-samples via the locked cross-fitted train-only benefit pseudo-outcome; the ordering uses only within-partition expression spread. No assessment rows, outcomes, or covariate summaries. No hard-coded gene symbols, patient indices, or assessment rankings. No gene-by-ACT product. Estimand, gates, bootstrap, threshold, and the locked clinical comparator untouched. Development CSVs only.
- watch: this deliberately gives back some benefit signal to buy stability, so the risk is landing in the dead zone -- Jaccard still under 0.10 *and* repeat 1's increment pushed back under zero. Repeat 1's increment (+0.473 in iter_006) is the margin being spent; if Jaccard clears 0.10 but repeat 1 turns negative, the gate trade has simply reversed and the next move is the `MAX_GENES_PER_PATHWAY` coverage lever instead.
- run_id: run_007_20260822T021611Z
- eligible: false
- failed_gates: [all_repeat_genomic_increment_positive, gene_selection_jaccard_at_least_0_10]
- reward: -1000000.000
- diagnostic_score: -5.380
- repeat_1_increment: -0.056 months
- repeat_1_lcb: -3.052 months
- repeat_2_increment: +0.898 months
- repeat_2_lcb: -4.426 months
- repeat_range: 0.954 months
- development_alignment_cg: 4.352 months
- development_cindex_cg: 0.686
- gene_jaccard: 0.098
- development_act_usage_cg: tree_split=0.690; path_traversal=0.174; terminal_difference_mean=0.174, median=0.165, p10=0.098, p90=0.263, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.689; path_traversal=0.173; terminal_difference_mean=0.173, median=0.166, p10=0.097, p90=0.261, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.690; path_traversal=0.174; terminal_difference_mean=0.174, median=0.163, p10=0.099, p90=0.265, nonzero_patients=1.000
- verdict: NOT_LEADER
- lesson: The DR-gate/spread-order split landed in the predicted dead zone -- Jaccard rose 0.076 -> 0.098 but still missed the 0.10 floor by 0.0016, while repeat 1's increment was spent from +0.473 down to -0.056, so trading signal for stability along this one axis cannot reach both gates at once.

#### The trade curve, now with three points
| run | member ordering | Jaccard | repeat 1 increment | repeat 2 increment |
|---|---|---|---|---|
| iter_004 | spread only | 0.131 | -0.259 | +0.150 |
| iter_007 | DR-gated, spread-ordered | 0.098 | -0.056 | +0.898 |
| iter_006 | DR frequency only | 0.076 | +0.473 | +1.824 |

The three points are monotone: every unit of member-ordering stability costs
increment and vice versa. No point on this curve satisfies both gates, so the
next move must raise Jaccard through a mechanism *other* than the ordering key,
leaving the increment margin intact.

#### Refinement of the iter_006 note on n_genes (a testable prediction)
The iter_006 entry claimed Jaccard is scale-invariant in `n_genes`. That
assumed shared genes grow proportionally with set size. A more careful overlap
model disagrees. If two folds each select `p` pathways from correlated
rankings, the expected number of *shared pathways* grows about like `p^2`
rather than `p`, and each shared pathway contributes about `k^2/M` shared genes
for `k` members drawn from a pathway with `M` detectable members. Then

    E[shared genes] ~ p^2 * c * k^2 / M = c * (p*k)^2 / M = c * n_genes^2 / M

so shared genes grow with the *square* of `n_genes`, and the split between
`MAX_GENES_PER_PATHWAY` and the number of pathways drops out to first order.
Calibrating `c * n^2 / M` on iter_006 (`n=8`, Jaccard 0.076 => about 1.13 shared
genes) predicts for the same ordering:
- `n_genes=12`: about 2.53 shared, Jaccard about 0.118
- `n_genes=16`: about 4.50 shared, Jaccard about 0.164

Crucially `n_genes` costs nothing in crowding: under `module_count=1` the forest
sees exactly one module feature regardless. The only cost is dilution of the
module mean, and iter_002 warns that mixing unrelated pathway blocks into one
average can cancel. `n_genes=12` (three pathway blocks) is the balanced choice
and the prediction above is falsifiable at 0.118.

### iter_008 — twelve genes on the iter_006 ordering
- type: PARAM
- hypothesis: Every point on the member-ordering trade curve fails at least one gate, so Jaccard must be raised by a mechanism that does not spend increment. Under the overlap model derived in the iter_007 entry, expected shared genes grow as `n_genes^2 / M`, so restoring iter_006's pure DR-frequency ordering (increments +0.473 / +1.824, the largest margin observed) and raising `n_genes` from 8 to 12 should lift Jaccard from 0.076 to about 0.118 -- clearing the 0.10 floor -- at zero crowding cost, because `module_count=1` means the forest still sees exactly one module feature.
- changed: `train.py` restored to its iter_006 state (from `runs/run_006_.../train_snapshot.py`, discarding iter_007's median gate, which lost 0.53 months of repeat-1 increment for 0.022 of Jaccard), then `CANDIDATE["n_genes"]` 8 -> 12 and `CANDIDATE["name"]` -> `v2_pathway12_one_module_mtry035_freqorder`. The single new variable is `n_genes`; the ordering is returned to the known-best configuration rather than to a novel one. `MAX_GENES_PER_PATHWAY` stays 4, so the twelve genes come from three pathway blocks.
- red_line_audit: Only `train.py` edited. `n_genes=12` is inside the locked `4 <= n_genes <= 16` validation, and total features stay at 20 (19 clinical + 1 module) against the 35 cap, so nothing about the forest geometry or the locked clinical comparator changes. Detectability filtering still runs first. All selection quantities remain fitting-partition-only, from the locked cross-fitted train-only benefit pseudo-outcome. No hard-coded gene symbols, patient indices, or assessment rankings. No gene-by-ACT product. Estimand, gates, bootstrap, and threshold untouched. Development CSVs only.
- watch: dilution is the risk iter_002 flagged -- averaging three pathway blocks into a single standardized mean can cancel signal that eight genes from two blocks preserved. If Jaccard clears 0.10 but the increments collapse toward zero, the answer is fewer genes per module rather than more, and `module_count=2` at `n_genes=12` becomes the natural follow-up despite its crowding cost.
- run_id: run_008_20260822T022616Z
- eligible: false
- failed_gates: [all_repeat_genomic_increment_positive, gene_selection_jaccard_at_least_0_10]
- reward: -1000000.000
- diagnostic_score: -4.891
- repeat_1_increment: -0.369 months
- repeat_1_lcb: -4.123 months
- repeat_2_increment: +0.399 months
- repeat_2_lcb: -3.745 months
- repeat_range: 0.768 months
- development_alignment_cg: 3.946 months
- development_cindex_cg: 0.688
- gene_jaccard: 0.076
- development_act_usage_cg: tree_split=0.701; path_traversal=0.171; terminal_difference_mean=0.171, median=0.160, p10=0.092, p90=0.265, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.699; path_traversal=0.172; terminal_difference_mean=0.172, median=0.161, p10=0.091, p90=0.262, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.703; path_traversal=0.170; terminal_difference_mean=0.170, median=0.159, p10=0.092, p90=0.268, nonzero_patients=1.000
- verdict: NOT_LEADER
- lesson: `n_genes` is not a Jaccard lever -- shared genes grew exactly in proportion to set size (1.13 shared at n=8, 1.70 at n=12, a ratio of 1.51 against 12/8 = 1.50), leaving Jaccard flat at 0.076 while diluting the module across three pathway blocks cost 0.84 months of repeat-1 increment.

#### My iter_007 overlap model is falsified; the iter_006 note was right
The iter_007 entry predicted `E[shared] ~ n_genes^2 / M` and therefore
Jaccard about 0.118 at `n_genes=12`. Observed was 0.076. Inverting
`Jaccard = s/(2n-s)` on the two runs gives `s = 1.126` at `n=8` and `s = 1.701`
at `n=12`; the ratio 1.51 matches `12/8 = 1.50` to two digits, so shared genes
grow **linearly** in `n_genes` and Jaccard is scale-invariant, exactly as the
iter_006 entry originally claimed. The `p^2` step -- assuming the number of
shared pathways grows with the square of the pathways drawn -- is what was
wrong; in practice folds disagree at the level of *which pathway family wins*,
so drawing more pathways adds mostly non-overlapping ones. The iter_007
refinement should be disregarded and the iter_006 statement restored.

#### Consequence: the remaining lever is pathway-ranking stability
With member ordering exhausted (the iter_007 trade curve) and `n_genes` ruled
out, the only untouched source of Jaccard is *which pathways* the folds agree
on. The cross-fold frequency tables have said the same thing since iter_005:
folds reproducibly land on a small number of pathway *families* -- a
branched-chain-amino-acid block (ACADSB, HIBCH, AUH, ACAD8, ALDH6A1, MCCC1), a
TGF-beta/SMAD block (SMAD7, TGFBR3, SMAD4, SMAD5, BMPR1A) and now a
PI3K/FOXO block (PIK3CB, PTPN13, FOXO1, CDKN1B) -- but which family wins a
given fold alternates.

If pathway choice were made reproducible, three pathways sharing about `k^2/M`
genes each would give roughly 2.4 shared genes at `n_genes=12` (Jaccard 0.111)
or 1.6 at `n_genes=8` (Jaccard 0.111) -- just over the floor -- *without*
touching the member ordering that carries the increment. That is the next
target, and it is the pathway-level analogue of iter_007's gate: let the DR
signal decide which pathways are *eligible*, and let a fold-stable key decide
the order among those equals. Unlike iter_007, this spends no member-level
signal, because every pathway in the pool is already DR-strong.

### iter_009 — stable ordering among DR-equal pathways
- type: ALGO
- hypothesis: Folds reproducibly identify the same handful of benefit-associated pathway *families* but disagree on which one wins, and that disagreement -- not member ordering and not `n_genes` -- is what holds Jaccard at 0.076. Letting the DR co-selection count decide which pathways are *eligible* and a fold-stable covariate summary decide the order among those equals should make pathway choice reproducible and lift Jaccard over the 0.10 floor while leaving the member ordering that produced iter_006's +0.473 / +1.824 increments completely intact.
- changed: `train.py` restored to its iter_006 state (`n_genes=8`, DR-frequency member ordering), then `stability_select_genes` splits the pathway ranking in two: the existing count-based sort now produces `ranked`, its first `STABLE_PATHWAY_POOL = 8` entries form a DR-eligible pool, and that pool is re-ordered by descending mean expression spread of its detectable members (ties by pathway name) before the remaining pathways are appended unchanged. New module constant `STABLE_PATHWAY_POOL = 8`. `CANDIDATE["name"]` -> `v2_pathway8_one_module_mtry035_stablepath`. Nothing else changes: `n_genes=8`, `module_count=1`, `max_features=0.35`, threshold 0.00, `MAX_GENES_PER_PATHWAY=4`.
- red_line_audit: Only `train.py` edited. The gate is the locked cross-fitted train-only DR benefit statistic aggregated over fitting-partition half-samples; the tie-break is mean within-partition expression spread, a covariate summary of the fitting partition only. No assessment rows, outcomes, or covariate summaries are touched. Detectability filtering still runs first, so this is not v1's refuted unfiltered-stability artifact, and it is not v1's refuted pathway-smoothing experiment -- no score is smoothed across pathways, only the order among already-selected pathways changes. No hard-coded gene symbols, patient indices, or assessment rankings. No gene-by-ACT product. Estimand, gates, bootstrap, threshold, and the locked clinical comparator untouched. Development CSVs only.
- watch: with only eight pathways in the pool the spread tie-break could promote a pathway that is DR-eligible but weakly benefit-associated, which would show up as both increments falling toward zero while Jaccard rises; if that happens the pool is too wide and should be narrowed rather than abandoned.
- run_id: run_009_20260822T023638Z
- eligible: false
- failed_gates: [all_repeat_genomic_increment_positive, gene_selection_jaccard_at_least_0_10]
- reward: -1000000.000
- diagnostic_score: -4.718
- repeat_1_increment: -0.328 months
- repeat_1_lcb: -3.172 months
- repeat_2_increment: +0.022 months
- repeat_2_lcb: -4.368 months
- repeat_range: 0.350 months
- development_alignment_cg: 3.778 months
- development_cindex_cg: 0.687
- gene_jaccard: 0.059
- development_act_usage_cg: tree_split=0.704; path_traversal=0.171; terminal_difference_mean=0.171, median=0.158, p10=0.087, p90=0.271, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.715; path_traversal=0.174; terminal_difference_mean=0.174, median=0.161, p10=0.084, p90=0.275, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.694; path_traversal=0.169; terminal_difference_mean=0.169, median=0.155, p10=0.090, p90=0.267, nonzero_patients=1.000
- verdict: DIAGNOSTIC_LEADER
- lesson: Gating pathways by a noisy criterion and then ordering the pool by a stable one is *less* reproducible than ordering by the noisy criterion directly -- Jaccard fell to 0.059, the worst of the search -- because the pool boundary fluctuates across folds even when its top element does not.

#### Why the pathway-level gate failed where the logic looked sound
The member-level gate in iter_007 worked in the intended direction (Jaccard
0.076 -> 0.098) because it re-ordered members *within* an already-chosen
pathway, and pathway choice was held fixed. The pathway-level gate has no such
anchor: `pool = ranked[:8]` is itself re-estimated per fold, so re-ordering the
pool by expression spread makes the winner depend on *which eight pathways
happened to qualify* rather than on which pathway scored highest. The
count-based top element turns out to be more stable than the pool's boundary,
so replacing "highest DR count wins" with "highest spread among the top eight
wins" injected the boundary's noise directly into the winner. The watch item
anticipated a weakly-associated pathway being promoted -- both increments did
fall -- but not that stability itself would get worse.

#### This run is the diagnostic leader and must NOT be followed
`diagnostic_leader.txt` now names run_009 because `min_lcb - repeat_range` rose
to -4.718, helped mostly by the repeat range collapsing to 0.350. It is
simultaneously **worse than iter_006 and iter_007 on both gates that matter**
(Jaccard 0.059 vs 0.076 and 0.098; repeat 1 -0.328 vs +0.473 and -0.056). This
is exactly the situation red line 12 describes: the continuous diagnostic and
the eligibility target disagree, and the diagnostic is not a candidate. The
search continues from iter_007, not from run_009.

#### Standing position after nine slots
| run | Jaccard | R1 incr | R2 incr | distance to eligibility |
|---|---|---|---|---|
| iter_004 | 0.131 | -0.259 | +0.150 | needs +0.26 increment |
| iter_007 | 0.098 | -0.056 | +0.898 | needs +0.0016 Jaccard and +0.06 increment |
| iter_006 | 0.076 | +0.473 | +1.824 | needs +0.024 Jaccard (+32%) |
iter_007 is by a wide margin the closest point, missing both gates by hairs.
Every axis tried so far (member ordering, `n_genes`, pathway ordering) *trades*
one gate against the other. What is needed now is a change that reduces the
variance of the selection statistic itself, improving reproducibility without
giving up benefit signal.

### iter_010 — 100 stability subsamples on the iter_007 geometry
- type: PARAM
- hypothesis: iter_007 misses eligibility by 0.0016 of Jaccard and 0.056 months of repeat-1 increment, and both misses are dominated by Monte-Carlo noise in the selection statistic rather than by any signal/stability trade: with 40 draws a pathway whose true co-selection probability is 0.5 has a count standard deviation of sqrt(40*0.25) = 3.16 on a 0-40 scale, so pathways separated by a true gap of 0.1 in probability reorder across folds purely by chance. Raising `STABILITY_SUBSAMPLES` from 40 to 100 shrinks that standard deviation by sqrt(2.5) = 1.58 for both `pathway_counts` and `gene_top_counts`, which should raise cross-fold agreement *and* sharpen gene choice, moving both failing gates in the same direction for the first time.
- changed: `train.py` restored to its iter_007 state (median DR gate on members, spread ordering within the gate -- the configuration closest to eligibility), then `STABILITY_SUBSAMPLES` 40 -> 100 and `CANDIDATE["name"]` -> `v2_pathway8_one_module_mtry035_gated_s100`. Nothing else changes: `n_genes=8`, `module_count=1`, `max_features=0.35`, threshold 0.00, `MAX_GENES_PER_PATHWAY=4`, `STABILITY_TOP_K=300`. iter_009's pathway-level gate is discarded.
- red_line_audit: Only `train.py` edited. `STABILITY_SUBSAMPLES` governs only how the train-only DR ranking is aggregated inside the fitting partition; it touches no estimand, gate, bootstrap, threshold, budget field, or comparator geometry, and the launcher still validates every fixed budget key. All half-samples are drawn from the fitting partition only, via the locked cross-fitted train-only benefit pseudo-outcome. Detectability filtering unchanged. No hard-coded gene symbols, patient indices, or assessment rankings. No gene-by-ACT product. Development CSVs only.
- runtime: iter_007 completed in 468.7 s, of which the selector is about 9 calls x 40 draws x 2 halves x 0.34 s = 245 s. At 100 draws that becomes about 612 s, for a projected total near 836 s against the 1,500 s wall -- roughly 56% of budget, a safe margin.
- watch: this is a variance-reduction argument, so if Jaccard barely moves it means the residual fold-to-fold disagreement is driven by the *data* (different 775-patient fitting partitions) rather than by subsampling noise, and no amount of extra draws will close it. That would be a genuine finding and would redirect the remaining slots away from the selector entirely.
- run_id: run_010_20260822T025243Z
- eligible: **true**
- failed_gates: []
- reward: -4.787
- diagnostic_score: -4.787
- repeat_1_increment: +0.146 months
- repeat_1_lcb: -2.733 months
- repeat_2_increment: +0.483 months
- repeat_2_lcb: -4.450 months
- repeat_range: 0.337 months
- development_alignment_cg: 4.246 months
- development_cindex_cg: 0.689
- gene_jaccard: 0.122
- development_act_usage_cg: tree_split=0.687; path_traversal=0.174; terminal_difference_mean=0.174, median=0.164, p10=0.095, p90=0.264, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.690; path_traversal=0.174; terminal_difference_mean=0.174, median=0.166, p10=0.092, p90=0.263, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.685; path_traversal=0.173; terminal_difference_mean=0.173, median=0.162, p10=0.099, p90=0.265, nonzero_patients=1.000
- verdict: SEARCH_LEADER
- lesson: Monte-Carlo noise in the selection statistic, not a signal/stability trade, was what held iter_007 short of both gates -- raising the aggregation from 40 to 100 half-sample pairs moved Jaccard 0.098 -> 0.122 and repeat 1 -0.056 -> +0.146 simultaneously, producing the first eligible candidate of the v2 search.

#### First eligible candidate: all ten gates
| gate | value | requirement |
|---|---|---|
| all_repeat_genomic_increment_positive | +0.146, +0.483 | > 0 |
| all_repeat_genomic_alignment_positive | 3.633, 4.858 | > 0 |
| all_repeat_genomic_value_at_least_clinical | 48.369 vs 47.827; 48.970 vs 48.451 | >= |
| all_repeat_genomic_value_at_least_best_constant | 48.369 vs 47.491; 48.970 vs 47.478 | >= |
| all_repeat_cindex_drop_no_more_than_0_03 | 0.688 vs 0.676; 0.689 vs 0.685 | genomic is *ahead* in both |
| all_repeat_genomic_seed_agreement_at_least_0_85 | 0.975 | >= 0.85 |
| gene_selection_jaccard_at_least_0_10 | 0.122 | >= 0.10 |
| all_repeat_nontrivial_benefit_fraction_at_least_0_10 | 0.941 | >= 0.10 |
| all_repeat_raw_propensity_overlap_at_least_0_80 | 0.872, 0.869 | >= 0.80 |
| all_repeat_iptw_effective_sample_size_at_least_0_30n | 358.5, 366.3 | >= 310.2 |

Runtime 870.3 s against a projected 836 s and the 1,500 s wall. The selected
module is a branched-chain-amino-acid catabolism block plus chromatin
co-regulators: `ACAD8, ALDH6A1, MCCC1, HIBCH, MBIP, KAT2B, ZZZ3, PHF20L1`.

#### What the reward is now made of, and where the remaining headroom is
`reward = min_repeat_lcb - repeat_range = -4.450 - 0.337 = -4.787`. Because the
selection LCB is the 0.25th bootstrap percentile, `lcb ~= mean - 2.81 sd`, and
with means of 0.15-0.48 against standard deviations of 1.02-1.59 the reward is
almost entirely `-2.81 * max_repeat_sd`. Repeat 2's sd of 1.590 is the binding
term; repeat 1's is 1.022.

Writing the increment as `(2/n) * sum over policy-disagreeing patients of
+/- gamma`, its bootstrap sd is about `2 sqrt(f/n) * rms(gamma | disagree)` for
a disagreement fraction `f`, while its mean is about `2 f * m`. Solving the
observed sd of 1.590 at n=1034 gives `sqrt(f) * rms(gamma) ~= 25.6`, so at
`f ~= 0.10` the disagreeing patients carry an rms gamma near 80 months -- the
IPW-amplified tail of the AIPW pseudo-outcome, since a treated patient with
propensity near the 0.05 clip carries a weight up to 20.

Differentiating `reward ~ 2 f m - 2.81 * 2 sqrt(f/n) * rms` in `f` shows the
optimum lies above the feasible range at the observed `m ~= 2.4` months, i.e.
**broader correctly-signed deviation from the clinical policy raises the reward**,
which is the same direction as maximising the increment. The remaining slots
should therefore push the increment up while holding all ten gates, with the
repeat range watched because it is subtracted directly.

### iter_011 — fixed member shortlist instead of a size-dependent median gate
- type: ALGO
- hypothesis: iter_007's gate keeps the top *half* of a pathway's members, so its selectivity depends on pathway size -- for a 12-member pathway it shortlists 6 and the four picks are meaningfully DR-driven, but for a 40-member pathway it shortlists 20 and the spread ordering does nearly all the work, collapsing toward iter_004's benefit-blind behaviour. Replacing the quantile with a fixed shortlist of the top `MEMBER_POOL = 6` members by co-selection frequency makes DR influence uniform across pathway sizes and strictly increases it for large pathways, which should raise the increment (and hence the reward, since the reward rises with broader correctly-signed deviation) while spending part of the 22% Jaccard margin iter_010 earned.
- changed: In `stability_select_genes`, the median gate `cutoff = np.median(gene_top_counts[members])` is replaced by a fixed shortlist: members are ranked by descending co-selection frequency, the top `MEMBER_POOL = 6` form the shortlist, and the shortlist is spread-ordered and exhausted before the spread-ordered remainder. New module constant `MEMBER_POOL = 6`, chosen as 1.5x the per-pathway quota `MAX_GENES_PER_PATHWAY = 4` so the shortlist is selective but never smaller than the quota. `CANDIDATE["name"]` -> `v2_pathway8_one_module_mtry035_pool6_s100`. Everything else identical to iter_010: `STABILITY_SUBSAMPLES=100`, `n_genes=8`, `module_count=1`, `max_features=0.35`, threshold 0.00.
- red_line_audit: Only `train.py` edited. The shortlist uses `gene_top_counts`, accumulated from fitting-partition half-samples via the locked cross-fitted train-only benefit pseudo-outcome; the ordering uses within-partition expression spread. Detectability filtering still runs first. No assessment rows, outcomes, or covariate summaries. No hard-coded gene symbols, patient indices, or assessment rankings. No gene-by-ACT product. Estimand, gates, bootstrap, threshold, and the locked clinical comparator untouched. Development CSVs only.
- watch: this deliberately spends Jaccard margin (0.122 against the 0.10 floor) to buy increment. If Jaccard falls below 0.10 the run is ineligible and the shortlist should be widened to 8 rather than abandoned. Note that a failure here cannot cost the banked result: `best_run.txt` is only overwritten by an *eligible* run with a strictly higher reward, so run_010 stays frozen regardless of this outcome.
- run_id: run_011_20260822T031020Z
- eligible: true
- failed_gates: []
- reward: -4.802
- diagnostic_score: -4.802
- repeat_1_increment: +0.579 months
- repeat_1_lcb: -1.859 months
- repeat_2_increment: +0.474 months
- repeat_2_lcb: -4.697 months
- repeat_range: 0.105 months
- development_alignment_cg: 4.458 months
- development_cindex_cg: 0.688
- gene_jaccard: 0.124
- development_act_usage_cg: tree_split=0.697; path_traversal=0.175; terminal_difference_mean=0.175, median=0.163, p10=0.096, p90=0.266, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.697; path_traversal=0.174; terminal_difference_mean=0.174, median=0.163, p10=0.091, p90=0.264, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.697; path_traversal=0.175; terminal_difference_mean=0.175, median=0.162, p10=0.102, p90=0.269, nonzero_patients=1.000
- verdict: NOT_LEADER
- lesson: The fixed shortlist improved almost everything -- pooled increment 0.315 -> 0.527, repeat 1's LCB -2.733 -> -1.859 on the lowest bootstrap sd of the search (0.860), and the repeat range 0.337 -> 0.105 -- yet the reward fell 0.015 because repeat 2's sd rose to 1.624 and `min_lcb` is set entirely by repeat 2.

#### The signal/stability trade curve was a Monte-Carlo artifact
This is the most useful result since iter_003. At 40 draws, every increase in
DR selectivity cost Jaccard (0.131 -> 0.098 -> 0.076 across iter_004, 007, 006).
At 100 draws the same move went the *other* way: replacing the median gate with
a far more DR-selective fixed top-6 shortlist raised Jaccard slightly, 0.122 ->
0.124, while raising the pooled increment by 67%. So the apparent trade was not
intrinsic to using the DR signal -- it was noise in the DR estimate. Once the
statistic is aggregated over 100 half-sample pairs, using it *more* is free.

This retires the iter_007 trade-curve framing that governed slots 7-9 and opens
the shortlist width as a lever that can be pushed further rather than balanced.

#### Repeat 2's bootstrap sd is now the only thing that matters
`reward = min_lcb - range`, and the two repeats are wildly asymmetric:

| | repeat 1 | repeat 2 |
|---|---|---|
| increment | +0.579 | +0.474 |
| bootstrap sd | 0.860 | 1.624 |
| selection LCB | -1.859 | **-4.697** |

Repeat 1's LCB is 2.8 months better than repeat 2's on an almost identical mean,
purely through dispersion. Across all eleven runs repeat 2's sd has never gone
below 1.370 while repeat 1's has ranged 0.860-1.325, so this is a property of
repeat 2's fold partition (`fold_seed = FOLD_SEED + 20000`) rather than of any
candidate. The repeat range is already down to 0.105 and contributes almost
nothing, so essentially **reward = repeat 2's mean - 2.81 * repeat 2's sd**.
Raising repeat 2's mean is the only controllable route, since its dispersion has
proved insensitive to every configuration tried.

The best repeat 2 mean observed anywhere in the search is iter_006's +1.824 at
sd 1.370, i.e. an LCB of -2.040 -- 2.7 months better than the current leader's.
iter_006's member ordering is exactly a shortlist of width 4, which is now
reachable without the Jaccard penalty that made iter_006 ineligible.

### iter_012 — narrow the shortlist to the per-pathway quota
- type: PARAM
- hypothesis: iter_011 showed that at 100 draws the DR statistic is stable enough that using it *more* costs no reproducibility (shortlist 6 raised both the increment and Jaccard). Narrowing `MEMBER_POOL` from 6 to 4 makes the shortlist exactly the per-pathway quota, so the four members are taken purely on co-selection frequency -- which is precisely iter_006's member ordering, the configuration that produced the largest repeat-2 increment of the entire search (+1.824 at sd 1.370, an LCB of -2.040 against the current leader's -4.697). Since reward is now essentially `repeat 2 mean - 2.81 * repeat 2 sd`, recovering that mean at 100-draw stability should improve the reward substantially.
- changed: `MEMBER_POOL` 6 -> 4 and `CANDIDATE["name"]` -> `v2_pathway8_one_module_mtry035_pool4_s100`. Nothing else. At width 4 the shortlist equals `MAX_GENES_PER_PATHWAY`, so all shortlisted members are taken and the spread ordering becomes inert -- the selector reduces exactly to iter_006's top-4-by-frequency rule, now aggregated over 100 half-sample pairs instead of 40.
- red_line_audit: Only `train.py` edited, and only a selector aggregation constant. Detectability filtering still runs first. All frequencies come from fitting-partition half-samples via the locked cross-fitted train-only benefit pseudo-outcome; no assessment rows, outcomes, or covariate summaries. No hard-coded gene symbols, patient indices, or assessment rankings. No gene-by-ACT product. Estimand, gates, bootstrap, threshold, budget, and the locked clinical comparator untouched. Development CSVs only.
- watch: iter_006 at 40 draws had Jaccard 0.076, and scaling by iter_010's +23.6% variance-reduction gain predicts only about 0.094 at 100 draws -- below the floor. iter_011's evidence argues against that extrapolation, since shortlist 6 landed at 0.124 rather than the predicted decline, but Jaccard remains the gate at risk. A failure costs nothing: run_010 stays frozen in `best_run.txt`, and the fallback is shortlist 5.
- run_id: run_012_20260822T032622Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10]
- reward: -1000000.000
- diagnostic_score: -2.708
- repeat_1_increment: +1.291 months
- repeat_1_lcb: -1.407 months
- repeat_2_increment: +1.708 months
- repeat_2_lcb: -2.291 months
- repeat_range: 0.417 months
- development_alignment_cg: 5.430 months
- development_cindex_cg: 0.690
- gene_jaccard: 0.090
- development_act_usage_cg: tree_split=0.705; path_traversal=0.171; terminal_difference_mean=0.171, median=0.160, p10=0.089, p90=0.268, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.706; path_traversal=0.170; terminal_difference_mean=0.170, median=0.161, p10=0.088, p90=0.258, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.704; path_traversal=0.172; terminal_difference_mean=0.172, median=0.159, p10=0.089, p90=0.278, nonzero_patients=1.000
- verdict: DIAGNOSTIC_LEADER
- lesson: A shortlist equal to the per-pathway quota produced by far the strongest candidate ever measured here -- both increments above 1.29 months, both LCBs the best of the search, a would-be reward of -2.708 against the leader's -4.787 -- but Jaccard fell to 0.090 and the run is ineligible, so it is a diagnostic result only.

#### Shortlist width now has three calibrated points, all at 100 draws
| MEMBER_POOL | Jaccard | R1 incr | R2 incr | min LCB | range | score |
|---|---|---|---|---|---|---|
| 6 (iter_011) | 0.124 | +0.579 | +0.474 | -4.697 | 0.105 | -4.802 |
| 4 (iter_012) | 0.090 | +1.291 | +1.708 | -2.291 | 0.417 | -2.708 |
| median gate (iter_010) | 0.122 | +0.146 | +0.483 | -4.450 | 0.337 | -4.787 |

Narrowing the shortlist from 6 to 4 nearly tripled the pooled increment
(0.527 -> 1.499) and improved repeat 2's LCB by 2.4 months, at a cost of 0.034
Jaccard. The gate is 0.010 away. Linear interpolation between widths 6 and 4
puts width 5 at roughly Jaccard 0.107 and a pooled increment near 1.0, which
would clear the floor with a thin margin and still beat the current leader's
reward by well over a month.

#### This is a diagnostic result, not a candidate
Red line 12 applies directly: run_012 now heads `diagnostic_leader.txt` with a
score of -2.708, and it is **ineligible**. It cannot be frozen, cannot be
nominated, and does not displace run_010 in `best_run.txt`. Its only legitimate
use is to form the next hypothesis, which is what the width-5 interpolation
above does.

#### Note on repeat 2's dispersion
The iter_011 entry concluded repeat 2's sd was insensitive to configuration and
that only its mean was controllable. iter_012 supports that: repeat 2's sd came
in at 1.377, inside the 1.370-1.643 band every previous run has occupied, while
its mean moved from +0.474 to +1.708. The reward improvement came entirely
through the mean, exactly as predicted.

### iter_013 — shortlist width 5
- type: PARAM
- hypothesis: Shortlist width is now calibrated at three points and behaves monotonically: width 6 gives Jaccard 0.124 with a pooled increment of 0.527, width 4 gives Jaccard 0.090 with 1.499. Width 5 should land near Jaccard 0.107 -- clearing the 0.10 floor -- while retaining most of width 4's increment, which under `reward = min_lcb - range` should beat the current leader's -4.787 by more than a month.
- changed: `MEMBER_POOL` 4 -> 5 and `CANDIDATE["name"]` -> `v2_pathway8_one_module_mtry035_pool5_s100`. Nothing else; `STABILITY_SUBSAMPLES=100`, `n_genes=8`, `module_count=1`, `max_features=0.35`, threshold 0.00 all unchanged.
- red_line_audit: Only `train.py` edited, and only a selector aggregation constant. Detectability filtering still runs first. All frequencies come from fitting-partition half-samples via the locked cross-fitted train-only benefit pseudo-outcome; no assessment rows, outcomes, or covariate summaries. No hard-coded gene symbols, patient indices, or assessment rankings. No gene-by-ACT product. Estimand, gates, bootstrap, threshold, budget, and the locked clinical comparator untouched. Development CSVs only.
- watch: the predicted Jaccard of 0.107 carries only a 7% margin over the floor and the interpolation is linear across a two-point span, so it could land either side. run_010 remains frozen in `best_run.txt` either way; if width 5 misses, the remaining route is to buy Jaccard back with more aggregation draws at width 5 rather than to widen the shortlist again.
- run_id: run_013_20260822T034226Z
- eligible: false
- failed_gates: [all_repeat_genomic_increment_positive]
- reward: -1000000.000
- diagnostic_score: -5.791
- repeat_1_increment: -0.507 months
- repeat_1_lcb: -4.448 months
- repeat_2_increment: +0.836 months
- repeat_2_lcb: -4.151 months
- repeat_range: 1.343 months
- development_alignment_cg: 4.095 months
- development_cindex_cg: 0.688
- gene_jaccard: 0.105
- development_act_usage_cg: tree_split=0.703; path_traversal=0.173; terminal_difference_mean=0.173, median=0.163, p10=0.095, p90=0.269, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.699; path_traversal=0.172; terminal_difference_mean=0.172, median=0.164, p10=0.093, p90=0.265, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.707; path_traversal=0.174; terminal_difference_mean=0.174, median=0.162, p10=0.096, p90=0.273, nonzero_patients=1.000
- verdict: INELIGIBLE
- lesson: Width 5 cleared Jaccard at 0.105 almost exactly as interpolated, but repeat 1's increment fell to -0.507 -- between width 4's +1.291 and width 6's +0.579 -- proving the increment is **not** monotone in shortlist width and that differences of about a month between neighbouring configurations are indistinguishable from noise.

#### The non-monotonicity is the finding, and it is a warning
Repeat 1's increment against shortlist width, everything else held fixed at 100
draws:

| MEMBER_POOL | 4 | 5 | 6 |
|---|---|---|---|
| repeat 1 increment | +1.291 | **-0.507** | +0.579 |
| repeat 2 increment | +1.708 | +0.836 | +0.474 |
| Jaccard | 0.090 | 0.105 | 0.124 |

Jaccard is cleanly monotone; repeat 2's increment is monotone; repeat 1's
increment swings 1.8 months and reverses sign. Changing a shortlist width by one
gene changes the module by at most one or two members, so a 1.8-month swing on
that perturbation means the *increment differences I have been optimising over
the last few slots are inside the noise floor*.

Concretely, this means run_012's would-be reward of -2.708 should **not** be read
as "the reward available at width 4". It is one draw from a distribution whose
spread is at least a month, and the same is true of run_010's -4.787 and
run_011's -4.802. Those two eligible rewards differ by 0.015 and are, on this
evidence, the same number.

#### Consequences for the remaining slots
1. Chasing the largest observed reward across near-identical configurations is
   fitting noise, which is exactly what red line 12 exists to prevent. The
   remaining slots should look for a configuration with *comfortable margin on
   both gates*, not the luckiest single draw.
2. The only lever proven to improve reproducibility without costing benefit
   signal is reducing the variance of the selection statistic (iter_010). More
   aggregation draws is one route, but extrapolating the observed gains
   (+18.4% Jaccard for 40 -> 100 draws in this configuration family) predicts
   only about 0.098 at 150 draws and 0.103 at 200, while 200 draws projects to
   roughly 1,490 s against the 1,500 s wall -- an unacceptable timeout risk for
   a marginal gain.
3. A cheaper and untried route to the same end is to tighten the detectability
   filter. `GENE_IQR_PERCENTILE = 50` currently admits about 4,300 genes;
   raising it shrinks the candidate pool to better-measured genes, which both
   reduces ranking noise and, by pushing more pathways below
   `MIN_PATHWAY_MEMBERS = 12`, shrinks the set of eligible pathways so folds are
   likelier to agree on one. v1 made detectability filtering mandatory, so
   strengthening it is consistent with the inherited findings rather than a
   departure from them.

### iter_014 — stricter detectability filter at shortlist width 4
- type: PARAM
- hypothesis: Width 4 carries comfortably positive increments in both repeats (+1.291, +1.708) and misses only Jaccard, at 0.090. Raising `GENE_IQR_PERCENTILE` from 50 to 70 shrinks the candidate pool from roughly 4,300 to roughly 2,600 better-measured genes, which reduces DR ranking noise directly and, by pushing more pathways below `MIN_PATHWAY_MEMBERS = 12`, shrinks the eligible-pathway space so folds are likelier to agree on the same pathway. Both effects raise Jaccard without giving up the DR selectivity that produces the large increments, which is the combination no configuration has achieved yet.
- changed: `GENE_IQR_PERCENTILE` 50.0 -> 70.0, `MEMBER_POOL` 5 -> 4, `CANDIDATE["name"]` -> `v2_pathway8_one_module_mtry035_pool4_iqr70_s100`. The shortlist returns to the width-4 setting of iter_012, so the single new variable relative to that run is the detectability threshold. `STABILITY_SUBSAMPLES=100`, `n_genes=8`, `module_count=1`, `max_features=0.35`, threshold 0.00 unchanged.
- red_line_audit: Only `train.py` edited. This *strengthens* the detectability filter that v1 made mandatory rather than relaxing it, so it moves further away from v1's refuted near-floor-probe artifact, not toward it. Expression spread is a within-fitting-partition covariate summary; no assessment rows, outcomes, or covariate summaries are touched. All DR frequencies remain fitting-partition-only via the locked cross-fitted train-only benefit pseudo-outcome. No hard-coded gene symbols, patient indices, or assessment rankings. No gene-by-ACT product. Estimand, gates, bootstrap, threshold, budget, and the locked clinical comparator untouched. Development CSVs only. Runtime should fall slightly, since fewer genes enter each chunked least-squares solve.
- watch: per iter_013 the increments carry a noise floor of about a month, so a single positive or negative result here should not be over-read. The success criterion is *margin on both gates simultaneously* -- Jaccard clearly above 0.10 and both increments clearly above zero -- not the headline reward. The risk is that a stricter filter discards the genes actually carrying the benefit signal, which would show up as both increments falling together.
- run_id: run_014_20260822T035715Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10]
- reward: -1000000.000
- diagnostic_score: -3.643
- repeat_1_increment: +0.700 months
- repeat_1_lcb: -3.106 months
- repeat_2_increment: +0.834 months
- repeat_2_lcb: -3.509 months
- repeat_range: 0.135 months
- development_alignment_cg: 4.698 months
- development_cindex_cg: 0.690
- gene_jaccard: 0.084
- development_act_usage_cg: tree_split=0.704; path_traversal=0.171; terminal_difference_mean=0.171, median=0.162, p10=0.093, p90=0.268, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.702; path_traversal=0.168; terminal_difference_mean=0.168, median=0.160, p10=0.087, p90=0.263, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.706; path_traversal=0.175; terminal_difference_mean=0.175, median=0.163, p10=0.098, p90=0.274, nonzero_patients=1.000
- verdict: INELIGIBLE
- lesson: A stricter detectability filter moved Jaccard the wrong way (0.090 -> 0.084) even though it produced the most repeat-consistent increments of the search (+0.700 and +0.834, range 0.135), so shrinking the candidate pool churns pathway choice rather than stabilising it.

#### Why the mechanism was wrong
The hypothesis assumed a smaller, better-measured pool would reduce ranking
noise and shrink the eligible-pathway space toward agreement. What actually
happened is visible in the cross-fold frequency table, which became *more*
diffuse, not less: the branched-chain and TGF-beta families that had recurred
in every run since iter_005 were displaced by a scattered set (EP300 4,
PGRMC2 4, CHD9 3, SIRT1 3, then a long tail of ATPases and lipid-metabolism
genes at count 2). Dropping the lower-spread half of an already-filtered pool
removes members from established pathways, pushing some below
`MIN_PATHWAY_MEMBERS = 12` and reshuffling which pathways qualify at all. The
eligible-pathway space did shrink, but it shrank *unpredictably per fold*, which
is the opposite of what was wanted.

#### What now drives Jaccard, corrected
Comparing the frequency tables across iter_011 (pool 6, Jaccard 0.124) and
iter_012 (pool 4, Jaccard 0.090) shows both runs recovering the *same* two
pathway families. The runs differ in which members they draw from those
families. So member-level agreement, not pathway-level agreement, is the
binding term at this point in the search, and the controlling quantity is the
**fraction of a pathway that is taken**: for `k` members drawn from a pathway
with `M` detectable members, two folds sharing that pathway share about `k^2/M`
genes. Every configuration so far has used `MAX_GENES_PER_PATHWAY = 4` against
`M >= 12`, i.e. at most 33% coverage and typically far less.

This also explains why raising `MIN_PATHWAY_MEMBERS` would be the wrong move
(larger `M` lowers coverage) and points at the untried lever: raise `k`.

### iter_015 — one pathway per module: raise per-pathway coverage
- type: PARAM
- hypothesis: Member-level agreement is the binding term for Jaccard, and it is governed by the fraction of a pathway taken -- two folds sharing a pathway share about `k^2/M` genes for `k` of `M` detectable members. Every run so far has used `k = 4` against `M >= 12`. Setting `MAX_GENES_PER_PATHWAY` and `MEMBER_POOL` to 8 with `n_genes = 8` makes the module a single Reactome pathway's top eight members, roughly quadrupling `k^2` and so the shared-gene count whenever two folds agree on a pathway, which should lift Jaccard past 0.10 while keeping the width-4-style pure-frequency selectivity that produced the search's largest increments. It also makes the module pathway-coherent, which iter_002 found matters.
- changed: `train.py` restored to its iter_012 state (`GENE_IQR_PERCENTILE` back to 50.0, discarding iter_014's stricter filter, which lost 0.006 of Jaccard), then `MAX_GENES_PER_PATHWAY` 4 -> 8 and `MEMBER_POOL` 4 -> 8 so the shortlist again equals the per-pathway quota and members are taken purely on co-selection frequency. `CANDIDATE["name"]` -> `v2_pathway8_single_pathway_mtry035_s100`. `STABILITY_SUBSAMPLES=100`, `n_genes=8`, `module_count=1`, `max_features=0.35`, threshold 0.00, `MIN_PATHWAY_MEMBERS=12` unchanged, so the top pathway always has at least eight detectable members and supplies the whole module.
- red_line_audit: Only `train.py` edited. Detectability filtering runs first and is restored to its long-standing setting. All co-selection frequencies come from fitting-partition half-samples via the locked cross-fitted train-only benefit pseudo-outcome; expression spread is a within-fitting-partition covariate summary. No assessment rows, outcomes, or covariate summaries. No hard-coded gene symbols, patient indices, or assessment rankings. No gene-by-ACT product. Estimand, gates, bootstrap, threshold, budget, and the locked clinical comparator untouched. Development CSVs only.
- watch: this trades breadth for depth -- a fold now commits to one pathway instead of two, so if two folds pick different families they share nothing, and Jaccard is `P(same pathway) x` a much larger per-share overlap. The gamble fails if pathway-level agreement is below roughly 0.3, and it also fails if the winning pathways are large (coverage `8/M` collapses for `M` near 40). Per iter_013 the increments carry a noise floor near a month, so the criterion remains margin on both gates, not the headline reward.
- run_id: run_015_20260822T041330Z
- eligible: false
- failed_gates: [all_repeat_genomic_increment_positive, gene_selection_jaccard_at_least_0_10]
- reward: -1000000.000
- diagnostic_score: -5.815
- repeat_1_increment: -0.357 months
- repeat_1_lcb: -4.301 months
- repeat_2_increment: +1.157 months
- repeat_2_lcb: -2.799 months
- repeat_range: 1.515 months
- development_alignment_cg: 4.331 months
- development_cindex_cg: 0.689
- gene_jaccard: 0.056
- development_act_usage_cg: tree_split=0.698; path_traversal=0.173; terminal_difference_mean=0.173, median=0.164, p10=0.090, p90=0.264, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.703; path_traversal=0.169; terminal_difference_mean=0.169, median=0.163, p10=0.083, p90=0.256, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.693; path_traversal=0.176; terminal_difference_mean=0.176, median=0.165, p10=0.098, p90=0.273, nonzero_patients=1.000
- verdict: INELIGIBLE
- lesson: Concentrating the module on a single pathway did raise per-share overlap exactly as intended, but pathway-level agreement turned out to be only 0.14 -- far below the 0.3 the gamble needed -- so Jaccard collapsed to 0.056, the worst result since iter_009.

#### Pathway-level agreement is now measured directly, and it is low
Because each fold committed to exactly one pathway, the cross-fold frequency
table reads out pathway agreement without inference: **every gene sits at count
2**, i.e. four distinct pathways each won two of the eight folds. That gives
`4 * C(2,2) / C(8,2) = 4/28 = 0.143` as the probability two folds picked the
same pathway. Decomposing the observed Jaccard, `0.0556 / 0.143 = 0.39`, so when
two folds *did* agree on a pathway they shared roughly 4.5 of 8 genes -- the
coverage mechanism worked as designed. The gamble failed purely on the
`P(same pathway)` term, which the iter_015 watch item named as the failure
condition.

The module itself was cleanly pathway-coherent for the first time
(`ACAD8, ALDH6A1, MCCC1, HIBCH, DLD, ACADSB, AUH, SLC25A44`, all
branched-chain amino-acid catabolism), so interpretability improved while the
gate went backwards.

#### What fifteen slots have established about this arena
Jaccard and the increment are controlled by the same underlying quantity -- how
strongly the fold-specific DR benefit estimate is allowed to drive selection --
and the achievable frontier is narrow:

| configuration | Jaccard | increments | eligible |
|---|---|---|---|
| spread-ordered (iter_004) | 0.131 | -0.259 / +0.150 | no |
| median gate, 100 draws (iter_010) | 0.122 | +0.146 / +0.483 | **yes** |
| shortlist 6, 100 draws (iter_011) | 0.124 | +0.579 / +0.474 | **yes** |
| shortlist 5, 100 draws (iter_013) | 0.105 | -0.507 / +0.836 | no |
| shortlist 4, 100 draws (iter_012) | 0.090 | +1.291 / +1.708 | no |
| shortlist 4, strict IQR (iter_014) | 0.084 | +0.700 / +0.834 | no |
| single pathway (iter_015) | 0.056 | -0.357 / +1.157 | no |

Every configuration with increments above about +1 month sits at Jaccard
0.084-0.090; every configuration clearing Jaccard comfortably sits at increments
below about +0.6. The one intervention that improved both at once was raising
the aggregation from 40 to 100 half-sample pairs, which is a variance reduction
rather than a trade. That remains the only known way to move the frontier
outward, and it is the basis of the remaining hypotheses.

### iter_016 — shortlist 5 with 140 aggregation draws
- type: PARAM
- hypothesis: The only intervention that moved Jaccard and the increments together is variance reduction in the selection statistic (iter_010, 40 -> 100 draws: Jaccard +23.6% and repeat 1 from -0.056 to +0.146). Shortlist 5 at 100 draws already cleared Jaccard at 0.105 and failed only on repeat 1's -0.507, which iter_013 established is inside a noise floor of about a month. Raising `STABILITY_SUBSAMPLES` to 140 at shortlist 5 should widen the Jaccard margin to roughly 0.112 and, more importantly, sharpen gene choice enough to pull repeat 1 back above zero -- giving an eligible run whose increments sit near iter_014's +0.700 / +0.834 rather than the leader's +0.146 / +0.483, which would improve the reward by more than a month.
- changed: `train.py` restored to its iter_013 state (`MAX_GENES_PER_PATHWAY` back to 4, `MEMBER_POOL` 5, `GENE_IQR_PERCENTILE` 50.0), then `STABILITY_SUBSAMPLES` 100 -> 140 and `CANDIDATE["name"]` -> `v2_pathway8_one_module_mtry035_pool5_s140`. Single new variable relative to iter_013.
- red_line_audit: Only `train.py` edited, and only a selector aggregation constant. `STABILITY_SUBSAMPLES` governs how the train-only DR ranking is averaged inside the fitting partition; it touches no estimand, gate, bootstrap, threshold, budget field, or comparator geometry, and the launcher still validates every fixed budget key. Detectability filtering unchanged. No assessment rows, outcomes, or covariate summaries. No hard-coded gene symbols, patient indices, or assessment rankings. No gene-by-ACT product. Development CSVs only.
- runtime: at 100 draws this configuration ran 883 s, of which the selector is about 9 x 200 x 0.34 = 612 s and the remainder about 271 s. At 140 draws the selector projects to about 857 s for a total near 1,128 s against the 1,500 s wall -- about 75% of budget, leaving a safe margin. 200 draws was rejected earlier precisely because it projects to roughly 1,490 s.
- watch: iter_013 showed repeat 1's increment swinging 1.8 months between neighbouring configurations, so a positive result here is not proof the extra draws caused it, and a negative result is not proof they failed. The honest reading of either outcome must account for that noise floor, and the run only counts as an improvement if it is *eligible* and its reward exceeds run_010's -4.787.
- run_id: (failed -- experiment 16 recorded as `status: failed` in the ledger, 692.96 s)
- eligible: false
- failed_gates: [run aborted before gates were computed]
- reward: -1000000.000
- diagnostic_score: n/a
- repeat_1_increment: n/a
- repeat_1_lcb: n/a
- repeat_2_increment: n/a
- repeat_2_lcb: n/a
- repeat_range: n/a
- development_alignment_cg: n/a
- development_cindex_cg: n/a
- gene_jaccard: n/a
- development_act_usage_cg: n/a -- the run raised before any forest was fitted
- repeat_1_act_usage_cg: n/a -- the run raised before any forest was fitted
- repeat_2_act_usage_cg: n/a -- the run raised before any forest was fitted
- verdict: INELIGIBLE
- lesson: Raising the aggregation to 140 draws exposed a latent numerical fragility in the selector -- a non-converging Cox nuisance fit on one half-sample overflows to an infinite predicted RMST, which turns the benefit pseudo-outcome into NaN and kills the run -- so the extra draws cost a slot without testing the hypothesis.

#### Root cause, in order
1. `sklearn` `ConvergenceWarning`: a Cox fit inside
   `cross_fitted_benefit_pseudo_outcome` hit its iteration limit on a
   pathological stratified half-sample.
2. `sksurv/linear_model/coxph.py:124` `RuntimeWarning: overflow encountered in
   power` -- `np.power(baseline_survival, risk_score)` with a diverged
   `risk_score` overflows, so the integrated survival curve becomes `inf`.
3. `prepare.py:362` `invalid value encountered in multiply/add` --
   `phi1 = mu1 + treatment / propensity * (ipcw_time - mu1)` evaluates
   `0 * inf`, producing NaN.
4. `_robust_gamma` takes `np.percentile` of an array containing NaN, so both
   winsorizing bounds become NaN and the whole vector is clipped to NaN.
5. `prepare._gene_effect_scores` calls `Ridge().fit(clinical, gamma)`, which
   raises `ValueError: Input y contains NaN`.

Each `_dr_ranking` call performs 3 inner folds x 2 Cox fits, so 9 selector
calls x 140 draws x 2 halves is roughly 15,000 Cox fits per run against about
11,000 at 100 draws and 4,300 at 40. The failure was always latent; more draws
simply made it near-certain. Six prior runs at 40-100 draws got lucky.

#### The fix, and why it is legitimate in `train.py`
`_robust_gamma` lives in `train.py` and is used **only** to build the selector's
internal gene ranking. The locked evaluation path -- `aipw_arm_scores` and
`_summarize_cohort` on the outer folds -- is in `prepare.py` and is untouched,
so guarding the selector's copy cannot alter the estimand, the policy value, or
any gate. The guard replaces non-finite pseudo-outcome entries with the median
of the finite entries before winsorizing. This **keeps every patient in the
ranking computation**, so red line 4 is respected -- no censored patient, and
indeed no patient at all, is dropped; a pathological nuisance fit is simply
prevented from propagating an infinity into the ranking.

#### Budget consequence
Four slots remain (17-20). run_010 is unaffected and still frozen in
`best_run.txt` at reward -4.787.

### iter_017 — genomic decision threshold 0.25 on the proven-eligible geometry
- type: PARAM
- hypothesis: `reward = min_lcb - range` is dominated by dispersion, not by the mean: at iter_011 the repeat-2 mean was 0.472 against a bootstrap sd of 1.624, so `2.81 sd` is 4.6 months of the 4.7-month deficit. Writing the increment over policy-disagreeing patients gives `sd ~ 2 sqrt(f/n) rms(gamma)` and `mean ~ 2 f m`, so shrinking the recommended set by raising the genomic decision threshold cuts `f` and therefore cuts sd as `sqrt(f)` while costing the mean only linearly. A 20% sd reduction is worth about +0.9 of reward, far more than any mean gain observed. Crucially the threshold does **not** enter gene selection at all, so Jaccard stays at iter_011's 0.124 and that gate is not at risk.
- changed: `train.py` restored to its iter_011 state (`MEMBER_POOL=6`, `STABILITY_SUBSAMPLES=100`, `MAX_GENES_PER_PATHWAY=4`, `GENE_IQR_PERCENTILE=50.0`) -- the configuration that is *proven eligible* -- then `CANDIDATE["benefit_threshold_months"]` 0.00 -> 0.25 and `CANDIDATE["name"]` -> `v2_pathway8_one_module_mtry035_pool6_thr025`. Separately, `_robust_gamma` gains the non-finite guard diagnosed in experiment 16; at 100 draws no run has ever hit that path, so the guard is inert here and the threshold is the only effective change.
- red_line_audit: Only `train.py` edited. The threshold applies **only to the genomic policy**; `_summarize_cohort` hard-codes `model_threshold = 0.0` for the clinical arm, so the locked comparator keeps its zero-month threshold and cannot be weakened -- red line 9 is respected. 0.25 is inside the locked `0 <= benefit_threshold_months <= 3` validation. The NaN guard replaces non-finite selector pseudo-outcome entries with the finite median and **retains every patient**, so red line 4 (never drop censored patients) holds; it touches only `train.py`'s ranking copy, never `prepare.py`'s locked evaluation path, so the estimand is unchanged. No hard-coded gene symbols and no gene-by-ACT product. Development CSVs only.
- distinction from v1's refuted experiment: v1 refuted a **shared** 0.25-month threshold applied to both arms, which silenced a genomic policy whose benefit scale was much narrower than the clinical one. Here the threshold is genomic-only, the comparator stays locked at zero, and this candidate's median absolute predicted benefit is 0.971 months, so a 0.25 cut removes marginal recommendations rather than silencing the policy. This is the "genomic threshold calibration" lever named in `program.md`, not a repeat of v1's shared-threshold test.
- watch: `all_repeat_nontrivial_benefit_fraction_at_least_0_10` is evaluated as `|benefit| > max(threshold, 0.10)`, so it tightens to `|b| > 0.25` -- at a median |b| of 0.971 it should stay far above the 0.10 floor, but it is the gate this change puts at risk. The value gates could also regress as the recommended set shrinks toward all-observation. If both increments stay positive and sd falls, this beats run_010's -4.787.
- run_id: run_017_20260822T044401Z
- eligible: false
- failed_gates: [all_repeat_genomic_increment_positive]
- reward: -1000000.000
- diagnostic_score: -5.036
- repeat_1_increment: -0.661 months
- repeat_1_lcb: -3.168 months
- repeat_2_increment: +1.206 months
- repeat_2_lcb: -2.964 months
- repeat_range: 1.868 months
- development_alignment_cg: 4.204 months
- development_cindex_cg: 0.688
- gene_jaccard: 0.124
- development_act_usage_cg: tree_split=0.697; path_traversal=0.175; terminal_difference_mean=0.175, median=0.163, p10=0.096, p90=0.266, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.697; path_traversal=0.174; terminal_difference_mean=0.174, median=0.163, p10=0.091, p90=0.264, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.697; path_traversal=0.175; terminal_difference_mean=0.175, median=0.162, p10=0.102, p90=0.269, nonzero_patients=1.000
- verdict: INELIGIBLE
- lesson: Raising the genomic threshold cut the recommended fraction from 0.393 to 0.320 and left Jaccard untouched at 0.124 exactly as predicted, but the bootstrap sd barely moved (repeat 2: 1.624 -> 1.571, repeat 1: 0.860 -> 0.933) because the disagreement fraction `f` is a *symmetric* difference, and moving the genomic recommendation rate away from the comparator's 0.373 increases it.

#### Correction to the dispersion model used in iter_011 and iter_017
The hypothesis treated `f` as the size of the genomic recommended set, so that
shrinking the set would shrink `f` and hence `sd ~ sqrt(f)`. That is wrong. `f`
is the fraction of patients on which the two policies *disagree*, so it is
minimised when the genomic recommendation rate sits near the locked
comparator's 0.373 -- and it grows when the genomic rate moves away in
**either** direction. At threshold 0 the genomic rate was 0.393, essentially
matched; at threshold 0.25 it fell to 0.320, further away, so `f` rose and the
`sqrt(f)` saving never materialised. The 3% sd reduction observed is consistent
with no real effect.

This closes threshold calibration as a dispersion lever: threshold 0 already
sits near the `f`-minimising point, and both directions away from it increase
dispersion. The large swing in the *means* (repeat 1 +0.579 -> -0.661, repeat 2
+0.474 -> +1.206, range blown out to 1.868) is the same month-scale noise floor
iter_013 established.

#### v1's threshold warning is reproduced in a narrower form
`nontrivial_benefit_fraction` fell from 0.939 to 0.859 as the gate tightened to
`|b| > 0.25`, still far above its 0.10 floor, so the policy was not silenced --
consistent with the prediction that this candidate's benefit scale (median |b|
0.971) is wide enough to absorb a 0.25 cut, unlike v1's. The lever failed on
dispersion grounds, not because the policy went quiet.

#### Position with three slots left
`best_run.txt` still names run_010 at reward -4.787. The search has established
that pool-4-style pure-DR selection dominates on every axis except Jaccard
(iter_012: increments +1.291/+1.708, the two lowest repeat-2 sds of the search
at 1.377, would-be reward -2.708) and fails only at Jaccard 0.090. Every
attempt to buy that last 0.010 by trading signal has cost more than it bought.
The one untried direction with empirical support is *spreading the same
selectivity over more pathways*: pathway concentration was measured at
`p = 1` -> Jaccard 0.056 (iter_015) and `p = 2` -> Jaccard 0.090 (iter_012),
a 1.6x gain that the `n^2/M` model does not explain and that predicts further
gains at `p = 4`.

### iter_018 — same pure-DR selectivity spread over four pathways
- type: PARAM
- hypothesis: Pathway concentration has been measured twice and spreading wins: one pathway per fold gave Jaccard 0.056 (iter_015), two gave 0.090 (iter_012), a 1.6x gain the `n_genes^2 / M` model does not predict, because what matters is the chance that two folds match on *at least one* pathway and per-slot pathway agreement is only 0.143. Halving `MAX_GENES_PER_PATHWAY` again to 2, so eight genes come from four pathways, should extend that trend past the 0.10 floor while keeping iter_012's pure-frequency member selection -- the configuration with the search's largest increments (+1.291 / +1.708) and its lowest repeat-2 dispersion (sd 1.377, would-be reward -2.708).
- changed: `MAX_GENES_PER_PATHWAY` 4 -> 2, `MEMBER_POOL` 6 -> 2 so the shortlist again equals the per-pathway quota and members are chosen purely on co-selection frequency, `CANDIDATE["benefit_threshold_months"]` back to 0.00 (iter_017 closed the threshold lever), and `CANDIDATE["name"]` -> `v2_pathway8_four_pathways_mtry035_s100`. `STABILITY_SUBSAMPLES=100`, `n_genes=8`, `module_count=1`, `max_features=0.35`, `GENE_IQR_PERCENTILE=50.0`, `MIN_PATHWAY_MEMBERS=12` unchanged. The `_robust_gamma` non-finite guard from iter_017 is retained.
- red_line_audit: Only `train.py` edited. Detectability filtering runs first, unchanged. All co-selection frequencies come from fitting-partition half-samples via the locked cross-fitted train-only benefit pseudo-outcome. The guard retains every patient, so red line 4 holds. Threshold returns to zero for the genomic policy; the locked comparator's zero-month threshold was never touched. No hard-coded gene symbols, patient indices, or assessment rankings. No gene-by-ACT product. Estimand, gates, bootstrap, budget, and comparator geometry untouched. Development CSVs only.
- watch: iter_002 found that mixing unrelated pathway blocks into one standardized mean can cancel signal, and this puts four blocks into a single module -- the failure mode is Jaccard clearing 0.10 while both increments collapse toward zero. That would leave run_010 as the final answer and would mean the arena's Jaccard floor and its increment gate cannot be satisfied simultaneously by any selector in this family.
- run_id: run_018_20260822T050009Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10]
- reward: -1000000.000
- diagnostic_score: -3.024
- repeat_1_increment: +1.055 months
- repeat_1_lcb: -3.005 months
- repeat_2_increment: +1.036 months
- repeat_2_lcb: -2.422 months
- repeat_range: 0.019 months
- development_alignment_cg: 4.976 months
- development_cindex_cg: 0.693
- gene_jaccard: 0.070
- development_act_usage_cg: tree_split=0.707; path_traversal=0.171; terminal_difference_mean=0.171, median=0.162, p10=0.087, p90=0.270, nonzero_patients=1.000
- repeat_1_act_usage_cg: tree_split=0.711; path_traversal=0.168; terminal_difference_mean=0.168, median=0.159, p10=0.078, p90=0.262, nonzero_patients=1.000
- repeat_2_act_usage_cg: tree_split=0.703; path_traversal=0.174; terminal_difference_mean=0.174, median=0.164, p10=0.095, p90=0.277, nonzero_patients=1.000
- verdict: INELIGIBLE
- lesson: Spreading over four pathways reversed the trend rather than extending it (Jaccard 0.090 at two pathways, 0.070 at four), so pathway spread is non-monotone with an optimum at two -- even though this run produced the most repeat-consistent result of the search (+1.055 and +1.036, range 0.019) and its best C-index (0.693).

#### The pathway-spread curve is now complete and peaks at two
| pathways per fold | 1 (iter_015) | 2 (iter_012) | 4 (iter_018) |
|---|---|---|---|
| Jaccard | 0.056 | **0.090** | 0.070 |
| repeat 1 / repeat 2 increment | -0.357 / +1.157 | +1.291 / +1.708 | +1.055 / +1.036 |

Concentrating loses because per-slot pathway agreement is only 0.143, so a fold
that commits to one pathway usually shares nothing. Spreading too far loses for
the opposite reason: with `MAX_GENES_PER_PATHWAY = 2` each shared pathway
contributes only about `2^2/M` genes, so matches stop being worth anything. Two
pathways is the maximum of the product, and it is where iter_012 already sat.

#### What this run says about the arena
iter_018 is the strongest *scientific* result of the search -- both repeats
agree to within 0.019 months on an increment above one month, the genomic model
leads the locked comparator on C-index by 0.013 and on policy value by 0.90
months -- and it is **ineligible**, because eight genes drawn two-at-a-time from
four fold-specific pathways do not repeat across folds. The Jaccard floor is
doing exactly what it was designed to do: refusing a result whose gene identity
is not reproducible, however stable its *value* estimate looks. That is the
correct call, and it is worth recording that the two properties are genuinely
dissociable in this data.

#### Two slots left
`best_run.txt` still names run_010 at reward -4.787. The eligible zone of this
selector family requires a member shortlist at least 1.25x the per-pathway
quota; at exactly 1.25x (shortlist 5, quota 4) iter_013 reached Jaccard 0.105
and failed only on repeat 1's -0.507, which sits inside the month-scale noise
floor. That configuration at higher aggregation is the last credible hypothesis,
and it is the one experiment 16 was attempting when it crashed.
