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
