# Reactome causal-RMST RSF experiment log

Fixed budget: 20 full experiments. Baseline is defined in `train.py` but has not
yet consumed a full experiment. Smoke-test metrics are not scientific feedback.

Append each prespecified hypothesis before editing `train.py`, then complete the
same entry after the launcher returns.

Every completed entry must include the mandatory `train_act_usage_cg` and
`validation_act_usage_cg` lines from `program.md`, even if the run is
ineligible. Each line reports the ACT split-tree fraction, counterfactual
patient-tree path-traversal fraction, and the mean/median/p10/p90/nonzero
summary of each patient's fraction of trees reaching different ACT-versus-OBS
terminal nodes.

## Search plan (prespecified before experiment 1)

Structural reasoning available before any full run, derived from the locked code
rather than from any scored result:

- `increment = A60(C+G) - A60(C)` and `A60(d) = mean[(2d-1)*(phi1-phi0)]`, so the
  increment is driven almost entirely by patients where the C+G policy and the
  matched clinical policy disagree. Bootstrap variance therefore grows with the
  size of that disagreement set, and `reward = min(LCB_train, LCB_valid) - gap`
  penalises large, unstable disagreements twice.
- The C+G forest sees 19 clinical columns plus `n_genes` genes. With
  `max_features=0.50` and 31 columns, the ACT indicator is offered at a split
  only about half the time it is drawn, so deep gene blocks structurally crowd
  ACT out of the trees. `nontrivial_benefit_fraction >= 0.10` is a gate on the
  C+G forest actually producing non-degenerate counterfactual contrasts, so ACT
  crowd-out is the first thing that can make a candidate ineligible.
- Consequently the search prioritises (1) eligibility, (2) small stable gene
  sets with high `max_features` so ACT competes at every split, (3) shrinking
  the train/validation gap, and only then point-estimate increments.

Infrastructure calibration (timing only, no candidate scoring, no gate or reward
evaluation): `load_train_valid` takes ~13 s, gene scoring ~0.3 s per call, and a
300-tree forest fits in ~2.3 s, so a full 24-forest experiment costs roughly two
minutes against the 1500 s cap. Wall time is not a binding constraint; the
20-experiment ledger is.

### iter_001 — shipped baseline anchor
- type: PARAM
- hypothesis: The shipped `dr_pathway`/12-gene/`max_features=0.50` baseline is
  ineligible because the ACT indicator is crowded out of the C+G forest, giving a
  degenerate counterfactual contrast that fails the nontrivial-benefit gate.
- changed: none — `CANDIDATE` is evaluated exactly as shipped to anchor the
  search with full-data ACT-use diagnostics, LCB scale, and gene stability.
- red_line_audit: no test path is referenced; only `CANDIDATE` in `train.py` is
  in scope and it is unmodified; selection stays inside the locked train-only
  `dr_pathway` routine; no gene symbols, patient indices, or validation-derived
  quantities are hard-coded; the estimand, gates, and budget are untouched.
- run_id: run_001_20260821T222205Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10, train_genomic_alignment_positive, train_genomic_value_at_least_best_constant, train_genomic_value_at_least_clinical, train_nontrivial_benefit_fraction_at_least_0_10, validation_genomic_alignment_positive, validation_genomic_value_at_least_best_constant, validation_genomic_value_at_least_clinical]
- reward: -1000000.000
- train_increment: -1.810 months
- validation_increment: -5.942 months
- absolute_gap: 4.132 months
- train_cindex_cg: 0.672
- validation_cindex_cg: 0.672
- gene_jaccard: 0.000
- train_act_usage_cg: tree_split=0.025; path_traversal=0.010; terminal_difference_mean=0.010, median=0.010, p10=0.004, p90=0.014, nonzero_patients=0.999
- validation_act_usage_cg: tree_split=0.122; path_traversal=0.045; terminal_difference_mean=0.045, median=0.044, p10=0.020, p90=0.069, nonzero_patients=1.000

- verdict: INELIGIBLE
- lesson: `dr_pathway` on the full 8,647-gene universe is unusable — its
  top-one-gene-per-pathway rule gives **Jaccard 0.000**, three disjoint fold
  selections — and the 12-gene block crowds ACT out of the C+G forest
  (`tree_split=0.025` on train OOF, C+G ACT rate 0.004 versus the matched
  clinical model's 0.059), so the C+G policy collapses to all-observation and
  loses to its own comparator.

Anchor facts carried forward (from run_001, train and validation only):

- The AIPW ATE is **positive and consistent** in both cohorts (train +0.641,
  validation +0.632 months), so all-ACT is the best constant policy in both.
  `*_genomic_value_at_least_best_constant` therefore requires the
  observation-recommended group to have non-positive mean DR benefit; an
  all-observation policy can never be eligible here.
- The matched clinical comparator is strong (alignment +1.120 train, +5.310
  validation), so the C+G model must track it rather than diverge from it.
- `phi1 - phi0` has a per-patient SD near 72 months, so a disagreement set of
  size m carries bootstrap SD near `2*72*sqrt(m)/n`. Small, stable disagreement
  is the only route to an LCB near zero; large divergence is punished twice,
  once in the LCB and once in the gap.

Mechanism probe (disclosed, ledger-free): before choosing iteration 2 geometry,
forest structure was measured on **training fold 1 only**, reporting solely the
locked ACT-use transparency diagnostics that `README.md` states are "not
eligibility gates". No policy value, reward, gate, cohort comparison, or
validation row was computed. Results are in the iteration 2 hypothesis.

### iter_002 — minimal gene block with ACT-competitive geometry
- type: PARAM
- hypothesis: ACT crowd-out is driven by the size of the gene block and by
  per-node feature competition, not by tree depth, so a 4-gene `dr_gene` block
  with `min_samples_leaf=8` and `max_features=0.25` restores a non-degenerate
  counterfactual contrast and keeps the C+G policy close enough to the strong
  clinical comparator to clear the value gates.
- changed: `CANDIDATE` -> name `dr_gene_4_actcompetitive`, `selector`
  `dr_pathway`->`dr_gene`, `n_genes` 12->4, `rsf.n_estimators` 300->600,
  `max_depth` 6->9, `min_samples_leaf` 16->8, `min_samples_split` 32->16,
  `max_features` 0.50->0.25.
- probe evidence (train fold 1, ACT-use diagnostics only, no gates/values):
  tree-split fraction rose 0.025 (anchor geometry) -> 0.340 at
  4 genes/depth 9/leaf 8/mtry 0.25; the gene-block size dominates
  (4 genes 0.34, 8 genes 0.16, 12 genes 0.11, 20 genes 0.09), `leaf=8` roughly
  doubles usage over `leaf=16`, `mtry` 0.25 beats 0.50 beats 1.00 because the
  binary ACT indicator wins a log-rank split more often against fewer
  continuous rivals, and depth beyond 9 adds nothing. Tree count is raised to
  600 purely as Monte-Carlo variance control for the seed-agreement gate, since
  `leaf=8` increases per-tree variance.
- risk accepted: 4 genes give the fewest chances for cross-fold overlap, so
  `gene_selection_jaccard_at_least_0_10` is the gate most likely to fail; that
  is exactly the question a stability-based ALGO selector would answer next.
- red_line_audit: only `CANDIDATE` changed; `dr_gene` is the locked train-only
  selector fitted on the fitting partition alone; no gene symbols, patient
  indices, validation-derived quantities, or gene-by-ACT product columns are
  introduced; estimand, gates, bootstrap, and budget untouched; no test artifact
  referenced.
- run_id: run_002_20260821T222916Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10, train_genomic_value_at_least_clinical, validation_genomic_alignment_positive, validation_genomic_value_at_least_best_constant, validation_genomic_value_at_least_clinical]
- reward: -1000000.000
- train_increment: -6.417 months
- validation_increment: -3.627 months
- absolute_gap: 2.789 months
- train_cindex_cg: 0.676
- validation_cindex_cg: 0.686
- gene_jaccard: 0.000
- train_act_usage_cg: tree_split=0.334; path_traversal=0.099; terminal_difference_mean=0.099, median=0.094, p10=0.065, p90=0.141, nonzero_patients=1.000
- validation_act_usage_cg: tree_split=0.528; path_traversal=0.141; terminal_difference_mean=0.141, median=0.139, p10=0.104, p90=0.181, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: The geometry hypothesis is confirmed — ACT tree-split use rose
  0.025 -> 0.334 (train OOF) and 0.528 (validation), and the nontrivial-benefit,
  train-alignment, and train-best-constant gates all flipped to pass — but
  `dr_gene` also returns **Jaccard 0.000**, so neither locked selector produces
  a reproducible gene set and the 4 noise genes actively degrade the policy
  relative to the matched clinical comparator (train increment -6.417).
- carried forward: with the same geometry the clinical comparator itself gets
  much stronger (train alignment +1.120 -> +8.738), so a gene block only helps
  if it is real; unstable genes are penalised twice, in the increment and in the
  `value_at_least_clinical` gates.

### iter_003 — train-only subsampling stability selection
- type: ALGO
- hypothesis: A single DR ranking of 8,647 genes on ~516 patients is dominated
  by selection noise, so ranking genes by how *often* they reach the top of an
  independently recomputed DR ranking across stratified subsamples of the
  fitting partition will recover whatever weak reproducible signal exists and
  lift `gene_selection_jaccard` above 0.10 without changing the estimand.
- changed: added `stability_select_genes` to `train.py` and passed it to
  `prepare.evaluate_candidate(..., selector=...)`. It draws 40 stratified 70%
  subsamples of the fitting partition (strata = OS_STATUS x ACT), recomputes the
  locked `cross_fitted_benefit_pseudo_outcome` and `_gene_effect_scores` inside
  each subsample, counts how often each gene lands in that subsample's top 100,
  and returns the `n_genes` genes with the highest selection frequency, ties
  broken by mean within-subsample rank then symbol. `CANDIDATE` is unchanged
  from iter_002 so the selector is the only difference.
- red_line_audit: the callback receives only the fitting partition, Reactome
  membership, training gene symbols, and the validated spec; every subsample,
  pseudo-outcome, nuisance fit, and score is computed inside the fitting
  partition, so no assessment-fold or validation outcome is touched; the
  pseudo-outcome is the locked train-only treatment-benefit target, satisfying
  the predictive-not-prognostic rule; no gene symbols, patient indices, or
  oracle features are hard-coded; no gene-by-ACT product columns are created;
  the RNG seed is a fixed constant, not tuned against validation; estimand,
  gates, bootstrap draws, and budget are untouched; no test artifact referenced.
- run_id: run_003_20260821T223312Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10, train_genomic_value_at_least_clinical, validation_genomic_alignment_positive, validation_genomic_value_at_least_best_constant, validation_genomic_value_at_least_clinical]
- reward: -1000000.000
- train_increment: -6.572 months
- validation_increment: -5.331 months
- absolute_gap: 1.241 months
- train_cindex_cg: 0.673
- validation_cindex_cg: 0.672
- gene_jaccard: 0.048
- train_act_usage_cg: tree_split=0.320; path_traversal=0.094; terminal_difference_mean=0.094, median=0.093, p10=0.066, p90=0.126, nonzero_patients=1.000
- validation_act_usage_cg: tree_split=0.493; path_traversal=0.129; terminal_difference_mean=0.129, median=0.127, p10=0.096, p90=0.164, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: Subsampling stability selection is directionally right — Jaccard rose
  0.000 -> 0.048 (UQCRH survived in two of three folds) and the train/validation
  gap halved from 2.789 to 1.241 months — but 40 draws of a gene-level ranking
  is still too noisy to clear 0.10, and the four genes still cost the policy
  ~6.6 months against the matched clinical comparator.
- carried forward: the four remaining non-Jaccard gates are all the same
  failure — the C+G policy is simply worse than the matched clinical policy
  (C+G ACT rate 0.174/0.201 versus clinical 0.363/0.417). At `max_features=0.25`
  the C+G forest sees only ~5 of 23 columns per node, so it is frequently forced
  onto a noise gene when the clinical split that the comparator used is not in
  the candidate set. Low `mtry` therefore buys ACT usage at the cost of
  decorrelating the C+G forest from its own comparator.

### iter_004 — Reactome-smoothed stability + greedy comparator-matched forest
- type: ALGO
- hypothesis: (a) Pooling stability evidence across Reactome pathway members
  averages out per-gene selection noise, so a gene promoted by both its own
  selection frequency and its pathway's selection frequency reproduces across
  folds and lifts Jaccard past 0.10; (b) `max_features=1.00` makes both matched
  forests greedy over the same column set, so the C+G forest reproduces the
  clinical comparator's splits wherever a clinical variable genuinely wins and
  deviates only where a gene beats it, shrinking the increment magnitude and
  turning the four value gates from large losses into near-ties.
- changed: `stability_select_genes` now also scores every Reactome pathway per
  subsample as the mean of its top-3 member scores (the locked `dr_pathway`
  convention), counts top-50 pathway appearances, and ranks genes by
  `gene_selection_frequency + best_containing_pathway_frequency`;
  `STABILITY_SUBSAMPLES` 40->60, `STABILITY_TOP_K` 100->60, new
  `STABILITY_TOP_PATHWAYS=50`; `CANDIDATE` name and `rsf.max_features`
  0.25->1.00. Everything else is held at iter_003 values.
- attribution note: the two changes touch disjoint reported metrics — Jaccard is
  a pure function of the selector, and the ACT-use and value diagnostics are a
  pure function of the forest given the genes — so a combined run is still
  attributable.
- red_line_audit: pathway membership is the pinned MSigDB Reactome collection
  already supplied to the callback, not an external or hard-coded list; all
  scores, pseudo-outcomes, and nuisance fits stay inside the fitting partition;
  no validation row, assessment-fold outcome, gene symbol, or patient index is
  referenced; no gene-by-ACT product column is created; `max_features=1.00` is
  applied identically to the clinical and C+G forests, preserving the matched
  comparator; estimand, gates, bootstrap, and budget untouched.
- run_id: run_004_20260821T223904Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10, train_genomic_value_at_least_clinical, validation_genomic_value_at_least_clinical]
- reward: -1000000.000
- train_increment: -3.686 months
- validation_increment: -0.249 months
- absolute_gap: 3.436 months
- train_cindex_cg: 0.662
- validation_cindex_cg: 0.666
- gene_jaccard: 0.000
- train_act_usage_cg: tree_split=0.213; path_traversal=0.058; terminal_difference_mean=0.058, median=0.054, p10=0.016, p90=0.105, nonzero_patients=1.000
- validation_act_usage_cg: tree_split=0.497; path_traversal=0.150; terminal_difference_mean=0.150, median=0.145, p10=0.077, p90=0.235, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: The two halves of the hypothesis split cleanly. (b) is confirmed and
  large: `max_features=1.00` cut the failed gates from five to three — both
  alignment gates and both best-constant gates now pass, and the validation
  increment moved -5.331 -> -0.249 months with C+G alignment +2.170 against the
  comparator's +2.419. (a) is refuted: pathway smoothing drove Jaccard 0.048 ->
  0.000, because `max`-over-containing-pathway frequency is a coarse, heavily
  tied statistic that swamps the finer gene-level stability evidence.
- carried forward: the only remaining failures are Jaccard and the two
  `value_at_least_clinical` gates, and validation is now nearly a tie (45.96 vs
  46.16). Arithmetically, at `n_genes=4` a single gene shared by all three folds
  gives Jaccard 0.143 and passes, whereas at `n_genes=8` two shared genes per
  pair are needed; the smallest gene block is therefore also the easiest
  stability target, so the fix must make the top-ranked gene reproducible rather
  than make the list longer.

### iter_005 — complementary-pairs simultaneous stability selection
- type: ALGO
- hypothesis: A gene that tops the DR ranking in one 70% subsample can still be
  driven by a handful of influential patients, so requiring *simultaneous*
  selection in two disjoint halves of the fitting partition — Shah-Samworth
  complementary pairs — keeps only genes whose treatment-benefit association
  survives being measured on completely different patients, which is exactly the
  property cross-fold reproducibility needs; the top-ranked gene should then
  repeat in all three folds and give Jaccard 0.143.
- changed: `stability_select_genes` reverted from the refuted pathway smoothing
  back to gene-level evidence, and its criterion changed from "in the top 60 of
  a 70% subsample" to "in the top 300 of *both* disjoint stratified halves of
  the fitting partition"; `STABILITY_SUBSAMPLES` 60->40, `STABILITY_TOP_K`
  60->300, `STABILITY_TOP_PATHWAYS` removed, `STABILITY_FRACTION` removed in
  favour of exact halves. `CANDIDATE` keeps every iter_004 value including the
  confirmed `max_features=1.00`, so the selector is the only difference.
- red_line_audit: both halves are drawn from the fitting partition only, so no
  assessment-fold or validation outcome enters selection; the target remains the
  locked train-only cross-fitted treatment-benefit pseudo-outcome, satisfying
  predictive-not-prognostic; no gene symbols, patient indices, or oracle
  features are hard-coded; the RNG seed is a fixed constant; no gene-by-ACT
  product columns; matched comparator, estimand, gates, bootstrap, and budget
  untouched; no test artifact referenced.
- run_id: run_005_20260821T224440Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10, train_genomic_value_at_least_clinical, validation_cindex_drop_no_more_than_0_03, validation_genomic_value_at_least_best_constant, validation_genomic_value_at_least_clinical]
- reward: -1000000.000
- train_increment: -3.354 months
- validation_increment: -0.993 months
- absolute_gap: 2.361 months
- train_cindex_cg: 0.646
- validation_cindex_cg: 0.652
- gene_jaccard: 0.048
- train_act_usage_cg: tree_split=0.250; path_traversal=0.075; terminal_difference_mean=0.075, median=0.070, p10=0.019, p90=0.135, nonzero_patients=1.000
- validation_act_usage_cg: tree_split=0.436; path_traversal=0.114; terminal_difference_mean=0.114, median=0.117, p10=0.050, p90=0.170, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: Complementary pairs did not beat plain subsampling — Jaccard is 0.048
  again, the identical value it reached in iter_003, with UQCRH surviving in
  folds 1 and 3 and fold 2 selecting a disjoint set for the third consecutive
  run — and the genes now also cost discrimination, newly failing
  `validation_cindex_drop_no_more_than_0_03` (0.652 versus 0.695).
- carried forward: making the *aggregation* more stringent cannot help while
  the underlying statistic is unstable. The locked pseudo-outcome
  `phi1 - phi0` has SD ~72 months on a quantity conceptually bounded by tau=60,
  because the propensity clip at 0.05 and the censoring-survival floor at 0.05
  each admit weights up to 20. A handful of patients therefore carry most of the
  squared variation, and the gene-versus-gamma correlation is decided by which
  extreme patients happen to have high expression — which differs in every
  fold. The instability is in the target, not in the aggregation.

### iter_006 — winsorized DR pseudo-outcome for gene scoring
- type: ALGO
- hypothesis: Gene-selection instability comes from the heavy tail of the DR
  pseudo-outcome rather than from how selections are aggregated, so winsorizing
  gamma at its within-partition 5th and 95th percentiles before scoring —
  leaving the locked scoring math and the locked pseudo-outcome construction
  untouched — will make the same complementary-pairs machinery reproduce the
  same genes across folds and lift Jaccard past 0.10, and should also recover
  the c-index the noise genes just cost.
- changed: added `WINSOR_PERCENT = 5.0` and a `_robust_gamma` helper; `_dr_ranking`
  now clips the locked `cross_fitted_benefit_pseudo_outcome` output to its own
  5th/95th percentiles before handing it to the locked
  `prepare._gene_effect_scores`. Nothing else changes: same complementary-pairs
  criterion, same 40 draws, same top-300 rule, same `CANDIDATE`.
- red_line_audit: the pseudo-outcome remains the locked train-only
  treatment-benefit target, so selection stays predictive rather than
  prognostic; percentiles are computed inside the half-sample being scored, so
  the transform is strictly fit-only and no assessment-fold or validation
  outcome is used; no censored patient is dropped — winsorizing bounds each
  patient's contribution but keeps all of them; the locked scoring routine and
  its clinical adjustment are reused unmodified; no hard-coded genes, patient
  indices, or gene-by-ACT products; estimand, gates, bootstrap, and budget
  untouched; no test artifact referenced.
- run_id: run_006_20260821T224914Z
- eligible: false
- failed_gates: [train_genomic_alignment_positive, train_genomic_value_at_least_best_constant, train_genomic_value_at_least_clinical, validation_genomic_value_at_least_best_constant]
- reward: -1000000.000
- train_increment: -5.204 months
- validation_increment: -1.620 months
- absolute_gap: 3.584 months
- train_cindex_cg: 0.651
- validation_cindex_cg: 0.666
- gene_jaccard: 0.111
- train_act_usage_cg: tree_split=0.248; path_traversal=0.050; terminal_difference_mean=0.050, median=0.047, p10=0.015, p90=0.082, nonzero_patients=1.000
- validation_act_usage_cg: tree_split=0.420; path_traversal=0.070; terminal_difference_mean=0.070, median=0.048, p10=0.019, p90=0.167, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: Confirmed and decisive for stability — winsorizing gamma lifted
  Jaccard 0.048 -> **0.111**, clearing the gene-stability gate for the first
  time (IL9 and IL22 both reproduce across folds), and it also recovered the
  validation c-index gate and passed `validation_genomic_value_at_least_clinical`
  (46.26 vs 46.16). The instability was in the target, exactly as diagnosed.
- carried forward: the binding constraint has moved. The winsorized genes raise
  the genomic model's estimated ATE (train +2.272, validation +3.032), so
  all-ACT becomes a strong constant and
  `*_genomic_value_at_least_best_constant` now demands that the
  observation-recommended group carry non-positive mean DR benefit against a
  large positive average. The C+G forest's predicted benefits are also strongly
  attenuated relative to the clinical comparator (median |benefit| 0.378 vs
  1.190 on train), so the shared 0.25-month threshold silences the genomic
  policy while barely touching the clinical one — an effectively unmatched
  comparison that the next iteration tests directly.

### iter_007 — matched decision scale via a zero-benefit threshold
- type: PARAM
- hypothesis: The C+G forest's counterfactual RMST contrasts are attenuated
  roughly threefold relative to the clinical forest's (median |benefit| 0.378
  versus 1.190 on train OOF), so a common 0.25-month threshold is not a matched
  decision rule — it silences the genomic policy while leaving the clinical one
  intact. Prespecifying the clinically neutral rule "recommend ACT whenever
  predicted counterfactual RMST60 is higher under ACT" removes that asymmetry,
  raises the C+G ACT rate toward the level the positive estimated ATE rewards,
  and should recover `train_genomic_alignment_positive` and both
  `*_value_at_least_best_constant` gates.
- changed: `CANDIDATE["benefit_threshold_months"]` 0.25 -> 0.00 and the name;
  nothing else changes, so the selector and forest are held at iter_006 values.
- red_line_audit: the threshold is applied identically to the clinical and C+G
  policies by the locked evaluator, preserving the matched comparator; the rule
  is prespecified here on the basis of the train-OOF benefit scale, not tuned to
  a validation outcome; recommendations still come from integrated
  counterfactual RMST curves, never from RSF mortality scores; no ACT rate is
  forced — the locked evaluator applies the rule and the gates judge it;
  estimand, gates, bootstrap, budget, and selector untouched; no test artifact
  referenced.
- run_id: (none — experiment 7 consumed without producing a result)
- eligible: false
- failed_gates: [n/a — no result was produced]
- reward: n/a
- train_increment: n/a
- validation_increment: n/a
- absolute_gap: n/a
- train_cindex_cg: n/a
- validation_cindex_cg: n/a
- gene_jaccard: n/a
- train_act_usage_cg: n/a — the run aborted before any cohort was summarised
- validation_act_usage_cg: n/a — the run aborted before any cohort was summarised
- verdict: INELIGIBLE
- lesson: Lost to an infrastructure failure, not to science. The host revoked
  filesystem access to the repository (macOS TCC on `~/Documents` and
  `~/Downloads`; `stat` still worked while `open`/`opendir` returned EPERM)
  while the run was assembling its result, so `sha256_file(TRAIN_CSV)` raised
  `PermissionError` after all three folds had been fitted, and the launcher's
  own failure-row write raised as well.
- ledger state: verified after access was restored — the lock manifest verifies,
  the ledger hash chain is intact across 13 rows, and experiment 7 has a
  `started` row (train sha `dbe3ce85896b`) with no `finished` row. Per red line
  12 a failure consumes a slot, so **7 of 20 are used and 13 remain**, and that
  exact `train.py` is permanently barred from consuming another slot by the
  launcher's duplicate check. The zero-threshold candidate is therefore retired
  rather than resubmitted under a new name, which would be a duplicate
  candidate; its hypothesis is revisited later only as part of a genuinely
  different candidate. The orphaned `runs/candidate_n2vwr5jo.json` temporary
  file from the aborted run is left untouched as a prior-run artifact.

### iter_008 — restrict selection to reliably measured genes
- type: ALGO
- hypothesis: Winsorizing gamma fixed the outcome side of the association but
  left the gene side unprotected, so the score is now won by near-floor probes
  whose residual norm is tiny and whose apparent signal comes from a handful of
  samples that are present in *every* subsample — reproducible and spurious at
  once. Restricting the selectable pool to genes in the upper half of
  within-partition expression variability should keep the Jaccard gain while
  replacing degenerate probes with genes the array actually measures, restoring
  `train_genomic_alignment_positive` and the value gates.
- changed: added `GENE_IQR_PERCENTILE = 50.0` and `_detectable_genes`;
  `stability_select_genes` now computes each candidate gene's interquartile
  range on the fitting partition and keeps only those at or above the median
  IQR of the offered gene universe before running the unchanged
  complementary-pairs, winsorized-gamma machinery. `CANDIDATE` is restored to
  its iter_006 values (`benefit_threshold_months` back to 0.25) so this run is a
  clean one-variable comparison against run_006.
- supporting evidence (training covariates only — no outcome, treatment,
  policy, gate, or reward quantity was computed): across the 8,647-gene
  universe the median expression SD is 0.529 and median IQR 0.644; run_006's
  selections sit at median IQR percentile **10.1** (OR3A2 0.6, IL9 2.8, MYF5
  3.4, MYH8 8.6), whereas the pre-winsorization selections of iter_003/005 sat
  at median percentile **62.0**.
- red_line_audit: the filter is a data-driven variability threshold computed
  inside the fitting partition, not a hard-coded gene list, and it names no
  symbol; it uses expression only, never an outcome, treatment, or validation
  quantity, so selection remains fit-only and still targets the locked
  train-only treatment-benefit pseudo-outcome; no patient is dropped; the RSF
  still receives raw, untransformed gene values so the matched comparator is
  unaffected; no gene-by-ACT products; estimand, gates, bootstrap, and budget
  untouched; no test artifact referenced.
- run_id: run_008_20260821T230408Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10, train_genomic_value_at_least_best_constant, train_genomic_value_at_least_clinical, validation_genomic_alignment_positive, validation_genomic_value_at_least_best_constant, validation_genomic_value_at_least_clinical]
- reward: -1000000.000
- train_increment: -4.287 months
- validation_increment: -2.623 months
- absolute_gap: 1.664 months
- train_cindex_cg: 0.655
- validation_cindex_cg: 0.673
- gene_jaccard: 0.000
- train_act_usage_cg: tree_split=0.246; path_traversal=0.068; terminal_difference_mean=0.068, median=0.064, p10=0.017, p90=0.126, nonzero_patients=1.000
- validation_act_usage_cg: tree_split=0.375; path_traversal=0.069; terminal_difference_mean=0.069, median=0.063, p10=0.022, p90=0.122, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: The filter did exactly what it was designed to do and the result is
  informative rather than encouraging — with near-floor probes excluded, Jaccard
  falls 0.111 -> 0.000. Taken with iter_006 this isolates the finding: the only
  reproducible gene signal available to a per-gene DR ranking in this cohort was
  an artifact of undetectable probes, and among reliably measured genes
  fold-level selection is entirely irreproducible. Train alignment did recover
  to +0.433 (from -0.484), so the filter improved gene *quality* while removing
  the spurious source of gene *stability*.
- carried forward: a per-gene statistic on 516 patients cannot be stabilised by
  any amount of internal resampling, because the three fitting partitions share
  only one third of their patients — internal resampling removes Monte-Carlo
  noise but not fitting-partition noise. The remaining legitimate route to
  stability is variance reduction by aggregation: a Reactome pathway mean over
  dozens of members has a far smaller sampling variance than any single gene's
  score, which is what iteration 9 tests. Note iter_004's pathway attempt failed
  for a different and now-understood reason — it used the locked top-3-member
  statistic, itself an extreme-value quantity, combined by a coarse `max` that
  produced large ties.

### iter_009 — Reactome pathway-mean stability selection
- type: ALGO
- hypothesis: A single gene's DR score on ~258 patients is too noisy to
  reproduce across fitting partitions that share only a third of their
  patients, but the mean score over all detectable members of a Reactome
  pathway averages that noise down by roughly the square root of the member
  count, so selecting the pathway that is most reproducibly top-ranked and then
  taking its best-scoring members will reproduce across folds and clear Jaccard
  while keeping the gene-quality gain of the detectability filter.
- changed: `stability_select_genes` now scores pathways, not genes, as the
  primary unit. For each complementary half it computes the locked
  winsorized-gamma gene scores over detectable genes, forms every Reactome
  pathway's *mean member score* via a sparse membership matrix (mean over all
  members, deliberately not the locked top-3 extreme-value statistic), and
  counts a pathway only when it lands in the top 20 of *both* halves; pathways
  are ranked by that simultaneous count, and genes are then taken from the
  best-ranked pathways in order of their score averaged over all 80 half-sample
  rankings. New constants `MIN_PATHWAY_MEMBERS = 12`, `TOP_PATHWAYS = 20`.
  `CANDIDATE` is unchanged from iter_008.
- red_line_audit: pathway membership is the pinned MSigDB Reactome collection
  handed to the callback, not an external or hard-coded list, and no gene symbol
  is named anywhere; every score, winsorization, half-split, and pathway mean is
  computed inside the fitting partition, so selection stays fit-only and no
  assessment-fold or validation outcome is touched; the target remains the
  locked train-only treatment-benefit pseudo-outcome, so selection is predictive
  rather than prognostic; no patient is dropped; the RSF still receives raw gene
  values and the same clinical block, preserving the matched comparator; no
  gene-by-ACT product columns; estimand, gates, bootstrap, and budget untouched;
  no test artifact referenced.
- run_id: run_009_20260821T230908Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10, train_genomic_alignment_positive, train_genomic_value_at_least_best_constant, train_genomic_value_at_least_clinical, validation_genomic_alignment_positive, validation_genomic_value_at_least_best_constant, validation_genomic_value_at_least_clinical]
- reward: -1000000.000
- train_increment: -4.737 months
- validation_increment: -3.112 months
- absolute_gap: 1.626 months
- train_cindex_cg: 0.671
- validation_cindex_cg: 0.682
- gene_jaccard: 0.000
- train_act_usage_cg: tree_split=0.243; path_traversal=0.066; terminal_difference_mean=0.066, median=0.055, p10=0.016, p90=0.135, nonzero_patients=1.000
- validation_act_usage_cg: tree_split=0.399; path_traversal=0.088; terminal_difference_mean=0.088, median=0.088, p10=0.043, p90=0.132, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: Pathway aggregation produced coherent biology but no reproducibility.
  Each fold now returns a single interpretable Reactome block — fold 1 TGF-beta
  receptor/SMAD (TGFBR3, SMAD4, SMAD5, ACVR2A), fold 2 nuclear pore and RAN
  transport (NDC1, RAN, SUMO1, NUP37), fold 3 immunoglobulin (IGLC2, IGLV3-25,
  IGKV4-1) — but each fold picks a *different* block, so Jaccard stays 0.000.
  Averaging over members reduced the variance of each pathway's score without
  making the pathway *ranking* reproducible, because the ranking is a
  competition among 1,839 near-zero means.
- carried forward: with `n_genes=4` and whole-pathway blocks, the three folds
  must agree on their single top-ranked pathway, which is a knife-edge. Taking a
  few genes from each of several top pathways instead means the folds only need
  to share *one* pathway anywhere in their short list, and one shared four-gene
  block already gives Jaccard 4/(32-4) = 0.143. This is the one structural
  change that converts a knife-edge into a tolerant criterion, and it is what
  iteration 10 tests.

### iter_010 — spread the block across several top pathways
- type: PARAM
- hypothesis: Requiring three fitting partitions to agree on a single top-ranked
  pathway is a knife-edge, but taking four genes from each of the four
  best-ranked pathways means the folds only need to share one pathway anywhere
  in a short list; a single shared four-gene block already yields Jaccard
  4/(32-4) = 0.143, so this should clear the stability gate without relying on
  the undetectable-probe artifact that iter_006 exposed.
- changed: `MAX_GENES_PER_PATHWAY = 4` added and enforced in the pathway walk,
  and `CANDIDATE["n_genes"]` 4 -> 16, so the block is drawn from four distinct
  Reactome pathways instead of one. The selector's scoring, winsorization,
  detectability filter, complementary-pairs criterion, and the forest are all
  unchanged from iter_009.
- risk accepted: 16 genes make the C+G forest 35 columns wide instead of 23, and
  the iteration-2 probe implies the ACT tree-split fraction falls from about
  0.21 to roughly 0.09, so `*_nontrivial_benefit_fraction_at_least_0_10` is the
  gate to watch; it currently passes with wide margin (0.632 train, 0.710
  validation), and the ACT-use diagnostics will show directly whether the
  margin is being spent.
- red_line_audit: only `CANDIDATE["n_genes"]`, the candidate name, and the
  per-pathway cap change; no gene symbol is hard-coded and pathway membership
  remains the pinned Reactome collection supplied to the callback; selection
  stays fit-only against the locked train-only treatment-benefit pseudo-outcome;
  16 genes plus 19 clinical columns is 35 features, inside the 51-feature and
  32-gene budget; matched comparator, estimand, gates, bootstrap, and budget
  untouched; no test artifact referenced.
- run_id: run_010_20260821T231840Z
- eligible: false
- failed_gates: [train_genomic_alignment_positive, train_genomic_value_at_least_best_constant, train_genomic_value_at_least_clinical, validation_cindex_drop_no_more_than_0_03, validation_genomic_alignment_positive, validation_genomic_value_at_least_best_constant, validation_genomic_value_at_least_clinical]
- reward: -1000000.000
- train_increment: -5.325 months
- validation_increment: -4.946 months
- absolute_gap: 0.380 months
- train_cindex_cg: 0.676
- validation_cindex_cg: 0.642
- gene_jaccard: 0.122
- train_act_usage_cg: tree_split=0.078; path_traversal=0.019; terminal_difference_mean=0.019, median=0.016, p10=0.003, p90=0.040, nonzero_patients=0.997
- validation_act_usage_cg: tree_split=0.123; path_traversal=0.024; terminal_difference_mean=0.024, median=0.019, p10=0.007, p90=0.049, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: The stability hypothesis is confirmed on legitimate genes — spreading
  the block across four pathways lifted Jaccard 0.000 -> **0.122** with a shared
  immunoglobulin block (IGKV1-17, IGKV4-1, IGLV2-14, IGLV3-19, IGLV3-25 each in
  two folds) drawn from detectable genes, not from the undetectable-probe
  artifact — and the train/validation gap fell to 0.380 months, the smallest so
  far. The accepted risk then materialised exactly as predicted: 16 genes cut
  the ACT tree-split fraction to 0.078 train / 0.123 validation, the C+G ACT
  rate collapsed to 0.000 in both cohorts, and an all-observation policy cannot
  survive `alignment_positive` or `value_at_least_best_constant` while the
  estimated ATE is positive. The 0.380 gap is an artifact of both cohorts
  collapsing to the same constant policy, not evidence of generalisation.
- carried forward: Jaccard and ACT usage are both solved, but at different block
  sizes — 16 genes buys stability and loses ACT, 4 genes buys ACT and loses
  stability. The per-pathway cap decouples them: what earns Jaccard is the
  number of *shared pathways*, not the number of genes, so taking two genes from
  each of four pathways preserves the tolerant criterion (one shared pathway
  still gives 2/(16-2) = 0.143) at half the block width.

### iter_011 — same pathway spread at half the block width
- type: PARAM
- hypothesis: Cross-fold overlap is earned by the number of *shared pathways*,
  not by the number of genes, so halving the per-pathway cap keeps the same four
  distinct Reactome blocks and the same tolerant Jaccard criterion — one shared
  pathway still gives 2/(16-2) = 0.143 — while an eight-gene block restores the
  ACT tree-split fraction from 0.078 to roughly 0.135 and with it a non-degenerate
  ACT recommendation rate, recovering the alignment, best-constant, and
  validation c-index gates.
- changed: `MAX_GENES_PER_PATHWAY` 4 -> 2 and `CANDIDATE["n_genes"]` 16 -> 8;
  the block still spans four pathways. Selector scoring, winsorization,
  detectability filter, complementary-pairs criterion, and every forest
  parameter are unchanged from iter_010.
- red_line_audit: only the candidate name, `n_genes`, and the per-pathway cap
  change; no gene symbol is hard-coded and membership remains the pinned
  Reactome collection passed to the callback; selection stays fit-only against
  the locked train-only treatment-benefit pseudo-outcome; 8 genes plus 19
  clinical columns is 27 features, inside budget; matched comparator, estimand,
  gates, bootstrap, and budget untouched; no test artifact referenced.
- run_id: run_011_20260821T232353Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10, train_genomic_alignment_positive, train_genomic_value_at_least_best_constant, train_genomic_value_at_least_clinical, validation_cindex_drop_no_more_than_0_03, validation_genomic_alignment_positive, validation_genomic_value_at_least_best_constant, validation_genomic_value_at_least_clinical]
- reward: -1000000.000
- train_increment: -4.935 months
- validation_increment: -7.254 months
- absolute_gap: 2.319 months
- train_cindex_cg: 0.671
- validation_cindex_cg: 0.644
- gene_jaccard: 0.048
- train_act_usage_cg: tree_split=0.152; path_traversal=0.040; terminal_difference_mean=0.040, median=0.033, p10=0.009, p90=0.080, nonzero_patients=1.000
- validation_act_usage_cg: tree_split=0.255; path_traversal=0.054; terminal_difference_mean=0.054, median=0.051, p10=0.023, p90=0.093, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: Halving the cap more than halved the overlap — Jaccard fell 0.122 ->
  0.048 while ACT tree-split recovered 0.078 -> 0.152 — because what two folds
  share is a *pathway*, and within a shared pathway their top-two genes need not
  be the same two. The cap cannot go below the width at which within-pathway
  gene rankings agree, so overlap must be bought with genes per pathway, not
  with more pathways.
- carried forward, and now the binding constraint: `n_genes` is squeezed from
  both sides by three different gates. At 16 genes Jaccard passes (0.122) but
  ACT collapses and validation c-index drops 0.053; at 8 genes ACT recovers but
  Jaccard fails and validation c-index still drops 0.051; only at 4 genes does
  `validation_cindex_drop_no_more_than_0_03` reliably pass (0.673 in run_008,
  0.682 in run_009). A four-gene block is therefore the *only* width that can
  satisfy the c-index, ACT-use, and nontrivial-benefit gates simultaneously, so
  the remaining question is whether four genes drawn as two genes from each of
  two pathways can reach Jaccard 0.10 — arithmetically one shared pathway
  contributing its two genes to one fold pair gives 2/(8-2) = 0.333, a mean of
  0.111, which clears the gate.

### iter_012 — four genes as two pathways of two
- type: PARAM
- hypothesis: Four genes is the only block width that clears the validation
  c-index, ACT-use, and nontrivial-benefit gates together, and drawing it as two
  genes from each of the two best-ranked pathways keeps the tolerant stability
  criterion — a single shared pathway contributing its two genes to one fold
  pair gives 2/(8-2) = 0.333 and a mean of 0.111 — so this is the narrowest
  configuration that can satisfy every structural gate at once.
- changed: `CANDIDATE["n_genes"]` 8 -> 4 and the candidate name; the per-pathway
  cap stays at 2, so the block spans two pathways. Selector scoring,
  winsorization, detectability filter, complementary-pairs criterion, and every
  forest parameter are unchanged from iter_011.
- risk accepted: the immunoglobulin pathway that supplied the shared genes in
  iter_010 and iter_011 ranked second or third, so restricting each fold to its
  top two pathways may drop it in some folds; the fold-frequency table will show
  directly whether the shared block survives.
- red_line_audit: only the candidate name and `n_genes` change; no gene symbol
  is hard-coded and membership remains the pinned Reactome collection passed to
  the callback; selection stays fit-only against the locked train-only
  treatment-benefit pseudo-outcome; matched comparator, estimand, gates,
  bootstrap, and budget untouched; no test artifact referenced.
- run_id: run_012_20260821T232817Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10, train_genomic_alignment_positive, train_genomic_value_at_least_best_constant, train_genomic_value_at_least_clinical, validation_genomic_alignment_positive, validation_genomic_value_at_least_best_constant, validation_genomic_value_at_least_clinical]
- reward: -1000000.000
- train_increment: -6.327 months
- validation_increment: -2.457 months
- absolute_gap: 3.869 months
- train_cindex_cg: 0.662
- validation_cindex_cg: 0.689
- gene_jaccard: 0.000
- train_act_usage_cg: tree_split=0.256; path_traversal=0.069; terminal_difference_mean=0.069, median=0.059, p10=0.021, p90=0.126, nonzero_patients=1.000
- validation_act_usage_cg: tree_split=0.426; path_traversal=0.094; terminal_difference_mean=0.094, median=0.088, p10=0.048, p90=0.148, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: The flagged risk materialised — restricting each fold to its top two
  pathways dropped the shared immunoglobulin block, and Jaccard fell to 0.000
  even though the c-index (0.662 train, 0.689 validation), ACT use (0.256 /
  0.426), and nontrivial-benefit gates were all comfortable. What two folds
  share is a mid-ranked pathway, not their top-ranked one, so a fold must be
  allowed to reach past its single best block for any overlap to exist.

Constraint map after twelve experiments (train and validation only):

- `train_genomic_value_at_least_clinical` has **never** passed. The matched
  clinical comparator's train value is fixed by geometry alone — 46.34 at the
  anchor's regularised geometry (run_001), 48.18 at depth 9 / leaf 8 /
  `mtry=1.00` (runs 004-012), and 50.20 at `mtry=0.25` (runs 002-003) — so the
  bar the gene block must clear is something the *forest*, not the genes, sets.
  The deficit is correspondingly smallest at the most regularised geometry
  (0.49 months in run_001) and largest where the comparator can overfit the AIPW
  noise. Regularisation, not better genes, is the lever on this gate.
- `max_features=0.25` protects the validation c-index (drop 0.016 in run_002)
  where `max_features=1.00` does not (drop 0.043 in run_005), but it also makes
  the clinical comparator luckiest of all, so it trades one gate for another.
- ACT use is extremely sensitive to leaf size once the block widens: the
  iteration-2 probe measured tree-split 0.100 at 8 genes with `leaf=8` but 0.005
  at 8 genes with `leaf=16`, so heavy leaf regularisation and a wide gene block
  cannot be combined.

### iter_013 — regularised forest with a restored shared-pathway block
- type: PARAM
- hypothesis: The gate that has never passed is set by the comparator's ability
  to overfit AIPW noise, not by gene quality, so pulling the forest back toward
  the anchor's regularised geometry (depth 6, `max_features=0.50`) should shrink
  the clinical train value from 48.18 toward the 46.34 seen at that geometry and
  bring the C+G deficit back to the 0.49 months of run_001; simultaneously
  restoring a four-gene-per-pathway cap at eight genes lets each fold reach its
  second pathway, which is where the shared immunoglobulin block lives, so one
  shared pathway gives 4/(16-4) = 0.333 for that pair and a mean of 0.111.
- changed: `MAX_GENES_PER_PATHWAY` 2 -> 4 and `CANDIDATE["n_genes"]` 4 -> 8, so
  the block spans two pathways of four; `rsf.max_depth` 9 -> 6 and
  `rsf.max_features` 1.00 -> 0.50. `min_samples_leaf` stays at 8 deliberately —
  the iteration-2 probe measured ACT tree-split collapsing from 0.100 to 0.005
  at eight genes when leaf went 8 -> 16, so leaf size is the one regularisation
  knob that cannot be turned here. Trees, split, selector scoring, winsorization,
  detectability filter, and threshold are unchanged.
- attribution note: the two changes again touch disjoint reported quantities —
  Jaccard is a pure function of the selector and `n_genes`, while the comparator
  value, c-index, and ACT-use diagnostics are pure functions of the forest.
- red_line_audit: only the candidate name, `n_genes`, the per-pathway cap, and
  two RSF hyperparameters change; the geometry is applied identically to the
  clinical and C+G forests, preserving the matched comparator; no gene symbol is
  hard-coded and membership remains the pinned Reactome collection passed to the
  callback; selection stays fit-only against the locked train-only
  treatment-benefit pseudo-outcome; 8 genes plus 19 clinical columns is 27
  features, inside budget; depth 6 and `max_features` 0.50 are inside the locked
  validation ranges; estimand, gates, bootstrap, and budget untouched; no test
  artifact referenced.
- run_id: run_013_20260821T233352Z
- eligible: false
- failed_gates: [gene_selection_jaccard_at_least_0_10, train_genomic_alignment_positive, train_genomic_value_at_least_best_constant, train_genomic_value_at_least_clinical, validation_cindex_drop_no_more_than_0_03]
- reward: -1000000.000
- train_increment: -4.746 months
- validation_increment: -0.265 months
- absolute_gap: 4.481 months
- train_cindex_cg: 0.671
- validation_cindex_cg: 0.667
- gene_jaccard: 0.022
- train_act_usage_cg: tree_split=0.090; path_traversal=0.024; terminal_difference_mean=0.024, median=0.022, p10=0.012, p90=0.042, nonzero_patients=1.000
- validation_act_usage_cg: tree_split=0.153; path_traversal=0.032; terminal_difference_mean=0.032, median=0.032, p10=0.015, p90=0.046, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: Failed gates fall to five and, for the first time, **all three
  validation value gates pass together** — alignment +1.909, value 46.45 against
  a 46.36 best constant and a 46.10 clinical comparator — with the validation
  increment down to -0.265 months. Regularisation barely moved the train
  comparator though (48.18 -> 48.05), so the train value gate still fails, and
  the validation c-index now misses by 0.004 (drop 0.034 against the 0.03
  allowance).
- carried forward, and this is the sharpest diagnosis yet: the train failure is
  a *decision-scale* failure, not a modelling failure. At this geometry the C+G
  forest's median |benefit| is 0.163 months on train and 0.212 on validation
  against the clinical forest's 0.937 and 1.279 — a sixfold mismatch — so the
  shared 0.25-month threshold sits above the entire C+G benefit distribution and
  silences it, leaving a train ACT rate of 0.017 and a policy that is
  all-observation in all but name. That is exactly why train alignment is
  -0.314, close to the -0.726 an all-observation policy would score. The
  threshold, not the forest, is what pins the train side.

### iter_014 — the only jointly feasible corner
- type: PARAM
- hypothesis: Thirteen experiments have pinned each gate to a different corner of
  the same two knobs, and exactly one corner satisfies all of them at once —
  sixteen genes across four pathways for stability (0.122 in run_010),
  `max_features=0.25` to stop those genes from costing discrimination (drop 0.016
  at 0.25 in run_002 versus 0.043 at 1.00 in run_005), depth 6 to keep the
  clinical comparator from overfitting AIPW noise, `min_samples_leaf=8` because
  ACT use collapses at wider leaves, and a zero threshold because the C+G
  benefit scale is six times narrower than the comparator's so a 0.25-month rule
  silences it. This run tests that corner.
- changed: `CANDIDATE` -> name `feasible_corner_pathway16_zerothresh`,
  `n_genes` 8 -> 16, `benefit_threshold_months` 0.25 -> 0.00,
  `rsf.max_features` 0.50 -> 0.25, `rsf.n_estimators` 600 -> 800;
  `MAX_GENES_PER_PATHWAY` stays 4, `max_depth` stays 6, `min_samples_leaf`
  stays 8. The selector is untouched. Tree count rises purely as Monte-Carlo
  variance control, since `max_features=0.25` raises per-tree variance and the
  seed-agreement gate must hold at 0.85.
- expected failure modes to read from the diagnostics: if ACT tree-split falls
  near the 0.078 of run_010 *and* the zero threshold still leaves the ACT rate
  near zero, the decision-scale explanation is wrong; if Jaccard drops below
  0.10 despite sixteen genes, `max_features` is interacting with selection,
  which it should not, since the selector never sees the forest.
- red_line_audit: only the candidate name and five `CANDIDATE` fields change;
  every parameter is applied identically to the clinical and C+G forests, so the
  comparator stays matched; the threshold is prespecified from the train-OOF
  benefit scale, not tuned to validation; recommendations remain integrated
  counterfactual RMST contrasts, never RSF mortality scores; no ACT rate is
  forced; 16 genes plus 19 clinical columns is 35 features and 800 trees is
  inside the 1,000-tree cap; no gene symbol, patient index, or validation
  quantity is hard-coded; estimand, gates, bootstrap, and budget untouched; no
  test artifact referenced.
- run_id: run_014_20260821T233818Z
- eligible: false
- failed_gates: [train_genomic_value_at_least_best_constant, train_genomic_value_at_least_clinical, validation_cindex_drop_no_more_than_0_03, validation_genomic_value_at_least_best_constant, validation_genomic_value_at_least_clinical]
- reward: -1000000.000
- train_increment: -7.750 months
- validation_increment: -2.997 months
- absolute_gap: 4.753 months
- train_cindex_cg: 0.672
- validation_cindex_cg: 0.658
- gene_jaccard: 0.122
- train_act_usage_cg: tree_split=0.052; path_traversal=0.014; terminal_difference_mean=0.014, median=0.013, p10=0.007, p90=0.022, nonzero_patients=1.000
- validation_act_usage_cg: tree_split=0.079; path_traversal=0.016; terminal_difference_mean=0.016, median=0.015, p10=0.009, p90=0.026, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: The corner splits cleanly into a confirmed half and a refuted half.
  Confirmed: the zero threshold is the right decision rule — the C+G ACT rate
  went 0.000 -> 0.160 train and 0.263 validation and **both alignment gates
  passed** (+0.332, +1.326), with Jaccard holding at 0.122 exactly as expected
  since the selector never sees the forest. Refuted: `max_features=0.25` does
  not protect discrimination at sixteen genes — the validation c-index drop
  worsened to 0.046 — and it made the clinical comparator markedly luckier
  (train value 48.05 -> 49.88), so it loses on both counts here; the 0.016 drop
  it bought in run_002 was a four-gene effect that does not survive a wide
  block. Seed agreement also fell to 0.897/0.864, close to its 0.85 floor,
  because a zero threshold puts many patients near the decision boundary.
- carried forward: with alignment solved, the wall is `validation_cindex` and
  the two `value_at_least_clinical` gates. Both respond to the same thing —
  how much the forest lets sixteen noise-carrying genes shape the risk ranking —
  and depth is the one regularisation knob not yet tried downward, since the
  iteration-2 probe showed depth barely affects ACT use (0.165/0.230/0.235 at
  depths 6/9/12) while shallower trees necessarily give genes fewer chances to
  enter a path and also keep the comparator from overfitting AIPW noise.

### iter_015 — shallow trees to stop genes shaping the risk ranking
- type: PARAM
- hypothesis: The remaining three gates all measure how much sixteen
  noise-carrying genes are allowed to shape the forest, so cutting `max_depth`
  to 4 — the one regularisation knob the iteration-2 probe showed to be nearly
  neutral for ACT use (tree-split 0.165/0.230/0.235 at depths 6/9/12) — should
  shorten every root-to-leaf path, give genes fewer chances to enter it, recover
  the validation c-index, and simultaneously keep the clinical comparator from
  overfitting the AIPW noise that inflated its value to 49.88, while restoring
  `max_features=0.50` undoes the refuted quarter-mtry setting.
- changed: `rsf.max_depth` 6 -> 4 and `rsf.max_features` 0.25 -> 0.50; the
  confirmed zero threshold, the sixteen-gene four-pathway block,
  `min_samples_leaf=8`, 800 trees, and the whole selector are held fixed.
- red_line_audit: only the candidate name and two RSF fields change, applied
  identically to both forests so the comparator stays matched; depth 4 and
  `max_features` 0.50 are inside the locked validation ranges; selection is
  untouched and remains fit-only against the locked train-only treatment-benefit
  pseudo-outcome; no gene symbol, patient index, or validation quantity is
  hard-coded; no gene-by-ACT products; estimand, gates, bootstrap, and budget
  untouched; no test artifact referenced.
- run_id: run_015_20260821T234304Z
- eligible: false
- failed_gates: [train_genomic_alignment_positive, train_genomic_value_at_least_best_constant, train_genomic_value_at_least_clinical, validation_cindex_drop_no_more_than_0_03, validation_genomic_alignment_positive, validation_genomic_value_at_least_best_constant, validation_genomic_value_at_least_clinical]
- reward: -1000000.000
- train_increment: -7.477 months
- validation_increment: -4.047 months
- absolute_gap: 3.430 months
- train_cindex_cg: 0.668
- validation_cindex_cg: 0.666
- gene_jaccard: 0.122
- train_act_usage_cg: tree_split=0.022; path_traversal=0.007; terminal_difference_mean=0.007, median=0.005, p10=0.003, p90=0.012, nonzero_patients=1.000
- validation_act_usage_cg: tree_split=0.029; path_traversal=0.006; terminal_difference_mean=0.006, median=0.005, p10=0.003, p90=0.011, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: Refuted, and the reason is a scope error in my own probe evidence. The
  depth-insensitivity I relied on (tree-split 0.165/0.230/0.235 at depths
  6/9/12) was measured at *four* genes; at sixteen genes depth is decisive,
  because a depth-4 tree has too few nodes for a 1-in-35 feature to enter, and
  ACT tree-split collapsed to 0.022 train / 0.029 validation, taking both
  alignment gates down with it. The validation c-index did improve to 0.666
  (drop 0.035, the best yet at sixteen genes) but still misses.
- carried forward: the train-OOF c-index gap is essentially zero (0.668 genomic
  versus 0.669 clinical) while the validation gap is 0.035, so the genes are not
  overfitting the *fitting* partition — their contribution simply fails to
  transport to the validation cohort. That is a transportability failure of the
  gene block itself, not a tuning failure, and no forest parameter tried so far
  removes it.

### iter_016 — every confirmed element, nothing refuted
- type: PARAM
- hypothesis: Three elements are now individually confirmed — the zero threshold
  (both alignment gates in iter_014), the sixteen-gene four-pathway block
  (Jaccard 0.122 in iter_010 and iter_014), and depth 6 with
  `max_features=0.50` (the weakest clinical comparator at 48.05 and all three
  validation value gates in iter_013) — while `max_features=0.25` and depth 4
  are refuted. Combining only the confirmed elements, at the 1,000-tree cap for
  maximum Monte-Carlo variance reduction, is the best-supported single point in
  the space and should hold alignment and stability while giving the validation
  c-index its best chance.
- changed: `rsf.max_depth` 4 -> 6, `rsf.n_estimators` 800 -> 1000; the zero
  threshold, sixteen-gene block, `min_samples_leaf=8`, `max_features=0.50`, and
  the selector are unchanged.
- red_line_audit: only the candidate name and two RSF fields change, applied
  identically to both forests so the comparator stays matched; 1,000 trees is
  exactly the locked cap and depth 6 is inside the locked range; selection is
  untouched and remains fit-only against the locked train-only treatment-benefit
  pseudo-outcome; no gene symbol, patient index, or validation quantity is
  hard-coded; no gene-by-ACT products; estimand, gates, bootstrap, and budget
  untouched; no test artifact referenced.
- run_id: run_016_20260821T234935Z
- eligible: false
- failed_gates: [train_genomic_alignment_positive, train_genomic_value_at_least_best_constant, train_genomic_value_at_least_clinical, validation_cindex_drop_no_more_than_0_03]
- reward: -1000000.000
- train_increment: -9.366 months
- validation_increment: 0.653 months
- absolute_gap: 10.018 months
- train_cindex_cg: 0.673
- validation_cindex_cg: 0.659
- gene_jaccard: 0.122
- train_act_usage_cg: tree_split=0.047; path_traversal=0.012; terminal_difference_mean=0.012, median=0.010, p10=0.005, p90=0.022, nonzero_patients=1.000
- validation_act_usage_cg: tree_split=0.076; path_traversal=0.016; terminal_difference_mean=0.016, median=0.014, p10=0.007, p90=0.027, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: Best gate count yet at four, and the validation cohort is now
  genuinely good — increment **+0.653 months**, C+G value 46.97 against a 45.86
  clinical comparator and a 46.97 best constant, alignment +2.336, all three
  validation value gates passing. But the train side collapsed the other way:
  at sixteen genes and depth 6 the ACT tree-split fraction is only 0.047 on the
  smaller fold-fitting partitions, the train ACT rate falls to 0.125, alignment
  goes to -2.990, and the gap blows out to 10.018 months. The two cohorts are
  now failing in opposite directions, which is the signature of a gene block too
  wide for the 516-patient fitting partitions to support.
- carried forward, and this closes the sixteen-gene branch: the validation
  c-index drop at sixteen genes has now been measured at four geometries —
  0.053 (depth 9, mtry 1.00), 0.046 (depth 6, mtry 0.25), 0.042 (depth 6, mtry
  0.50), and 0.035 (depth 4, mtry 0.50) — and only the depth that destroys ACT
  use gets near the 0.03 allowance. Sixteen genes cannot pass that gate at any
  workable geometry, so the block must narrow, which means the *within-pathway*
  gene choice has to become reproducible: two folds that share a pathway
  currently pick different members of it, which is exactly why eight genes gave
  Jaccard 0.022 in iter_013 while sixteen gave 0.122.

### iter_017 — reproducible within-pathway members at half the block width
- type: ALGO
- hypothesis: Pathway *choice* is already DR-driven and partly reproducible —
  two folds share a pathway — but which members they then take is decided by
  per-gene DR scores that do not reproduce, which is why eight genes gave
  Jaccard 0.022 while sixteen gave 0.122. Ordering the members of an already
  chosen pathway by a covariate summary that is near-identical across fitting
  partitions makes two folds that share a pathway take the *same* genes from it,
  so one shared pathway yields 4/(16-4) = 0.333 for that pair and a mean of
  0.111 at only eight genes — narrow enough for the validation c-index gate that
  closed the sixteen-gene branch.
- changed: within-pathway member ordering switched from mean DR score to
  descending within-partition expression IQR (`_detectable_genes` now returns
  the spreads it already computes, and the pathway walk sorts members by them);
  the final outside-pathway fallback keeps the DR ordering. `CANDIDATE["n_genes"]`
  16 -> 8. Pathway ranking, the complementary-pairs simultaneous criterion,
  winsorization, the detectability filter, and every forest parameter are
  unchanged from iter_016.
- red_line_audit: this is the one change that needs care against red line 9, so
  to be explicit — *which* pathways are eligible and how they rank is still
  decided entirely by the locked train-only treatment-benefit pseudo-outcome
  through the winsorized DR scores and the simultaneous top-20 criterion, so
  selection continues to target treatment benefit rather than prognosis; the IQR
  ordering only breaks the choice *among members of a pathway the DR signal has
  already chosen*, and it uses expression alone, never an outcome, treatment,
  survival time, or validation quantity. It is computed inside the fitting
  partition, so it stays fit-only. No gene symbol, patient index, or oracle
  feature is hard-coded; no gene-by-ACT products; 8 genes plus 19 clinical
  columns is 27 features; matched comparator, estimand, gates, bootstrap, and
  budget untouched; no test artifact referenced.
- run_id: run_017_20260821T235452Z
- eligible: false
- failed_gates: [train_genomic_value_at_least_clinical, validation_genomic_alignment_positive, validation_genomic_value_at_least_best_constant, validation_genomic_value_at_least_clinical]
- reward: -1000000.000
- train_increment: -5.682 months
- validation_increment: -2.411 months
- absolute_gap: 3.271 months
- train_cindex_cg: 0.663
- validation_cindex_cg: 0.696
- gene_jaccard: 0.111
- train_act_usage_cg: tree_split=0.132; path_traversal=0.038; terminal_difference_mean=0.038, median=0.030, p10=0.012, p90=0.075, nonzero_patients=1.000
- validation_act_usage_cg: tree_split=0.153; path_traversal=0.032; terminal_difference_mean=0.032, median=0.032, p10=0.014, p90=0.049, nonzero_patients=1.000

- verdict: NOT_LEADER
- lesson: Confirmed, and it reopens the narrow-block branch. Stabilising the
  within-pathway choice gave **Jaccard 0.111 at eight genes** — folds 2 and 3 now
  take the *same* immunoglobulin members (IGKV1-17, IGKV2D-28, IGKV4-1,
  IGLV3-25) instead of different ones — and halving the block width collapsed
  the validation c-index drop from 0.042 to **0.005**, clearing the gate that
  closed the sixteen-gene branch. ACT use stayed healthy (0.132 train, 0.153
  validation), and train alignment (+0.694) and the train best-constant gate
  both pass. Eight genes is where stability, discrimination, and ACT use finally
  coexist.
- carried forward: the four remaining failures are all value gates, and they
  have swapped cohorts relative to iter_016 — validation is now the weak side
  (alignment -0.729, value 44.93) where it was the strong one, which is the
  signature of noise rather than of a systematic defect. The one confirmed lever
  not currently applied is `max_features=1.00`, which in iter_004 cut the failed
  gates from five to three by making both matched forests greedy over the same
  columns so the C+G policy tracks its comparator; it also lowers the comparator
  itself, since the train clinical value runs 49.02 at `mtry=0.50` against 48.18
  at `mtry=1.00`, moving both sides of `value_at_least_clinical` the right way.
