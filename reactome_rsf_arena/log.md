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
