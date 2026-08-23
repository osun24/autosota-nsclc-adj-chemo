# Reactome GRF-CSF experiment log

Arena created before any full experiment. The initial candidate inherits v2's
detectability-filtered Reactome pathway stability selector, selects eight genes,
and exposes them raw to an honest causal survival forest. No search slot has
been consumed.

The former train and validation cohorts are both adaptive development data.
Test data has not been accessed by this arena.

---

## Experiment 1 (prespecified before launch)

### Hypothesis

Changing only the learner -- from v2's prognostic S-learner survival forest to
an honest RMST causal survival forest with ACT supplied as `W` -- is by itself
enough to turn the inherited detectability-filtered Reactome pathway-stability
selector (eight raw genes) into a repeat-consistent positive genomic policy
increment. Formally: `increment_r = A60(C+G) - A60(clinical) > 0` in both outer
repeats, with the worst-repeat selection-adjusted LCB above the repeat range.

### Rationale

`v2_findings.md` names this as the starting scientific question and explicitly
forbids inventing a new selector before the architectural change is measured on
its own. V2's RSF ACT crowd-out was architectural, so its remedies do not
transfer; the only way to know what the CSF architecture contributes is to hold
the selector fixed and change nothing else. This run is the reference against
which every later representation, regularization, and selector hypothesis is
scored.

### Exact change

None. `train.py` is the shipped baseline, unmodified, at
sha256 d53b4e23a695b82f8955cc4012f745fe71bf49379026052d5bfa014aeb75ad22.
Candidate `csf_pathway8_raw_honest_s100`: selector `dr_gene` wrapped by the
inherited pathway stability selection (100 subsamples, disjoint stratified
halves, detectability filter at the median gene IQR, >=12-member pathways,
top-20 pathways, <=4 genes per pathway), `n_genes=8`, `representation=raw`,
`num_trees=1000`, `min_node_size=5`, `sample_fraction=0.50`,
`honesty_fraction=0.50`, `alpha=0.05`, `imbalance_penalty=0.00`.

### Red-line audit

1. Test untouched: loader is restricted to the two development CSVs; no test
   path is read, hashed, enumerated, or referenced. OK.
2. Only `train.py` editable, and it is not edited at all this run. OK.
3. All adaptation is fit-only: selector sees only the fitting partition;
   imputers/scalers fit inside `FeatureTransformer` on the fit rows. OK.
4. ACT is `W` only; `PRETREATMENT_COLUMNS` excludes it; no gene-by-ACT terms. OK.
5. Estimand untouched: horizon 60, threshold 0, locked LCB and range penalty. OK.
6. Policy learning (GRF) and grading (locked IPCW-AIPW) stay separate. OK.
7. No censored patient is dropped anywhere in selection or evaluation. OK.
8. Selection uses only fit-partition DR benefit signal and Reactome membership;
   no hard-coded symbols. OK.
9. Identical geometry for clinical and C+G forests; both get the same spec. OK.
10. All mandatory CSF diagnostics recorded below. OK.
11. C-index reported from the separate prognostic forest, gate only. OK.
12. No probes: smoke is not run; this is one full launcher slot. OK.
13. No metric shopping: eligible reward is the leaderboard. OK.
14. Fixed compute: slot 1 of 20, 30-minute wall, <=1000 trees per forest. OK.
15. Claims remain observational. OK.

### Prespecified decision rule

- If eligible and reward > current best: promote to `best_run.txt` (launcher).
- If ineligible: read the failed gate set to choose experiment 2. Priority
  order per `program.md`: raw-vs-module representation first, then causal-forest
  regularization, and only then selector complexity.
- No post-hoc redefinition of success; the sentinel is accepted as-is.

### Operator error: slot 2 spent with no result (corrected record)

Experiment 1 ran to completion in 1282.3 s and wrote
`runs/run_001_20260823T070159Z`. While it was still running, an agent-side
process check used a malformed pattern that matched nothing, and the run was
wrongly judged to have been killed. On that false premise a second launcher was
started at 2026-08-23T07:01:28Z with a docstring-only edit of the same
baseline. Experiment 1 finished 31 seconds later, and the duplicate was
terminated once the mistake was found.

This is recorded, not worked around:

- Slot 2 is spent with no result. `_experiments` counts its `started` row, so
  18 of 20 experiments remain.
- Its `train.py` sha `98c4c6eb...` is permanently retired by that row.
- No metric, gate, or diagnostic from the duplicate was produced or observed,
  so no adaptive information entered the search from it. It was a
  behaviourally identical rerun of experiment 1, whose seeds and folds are
  fixed, so it could not have produced different numbers.
- Nothing about experiment 1's own conduct or result is affected; it is the
  unmodified shipped baseline at sha `d53b4e23...`.
- One artifact caveat: `run.py` copies `train.py` into the run directory only
  after the candidate returns, so `runs/run_001_20260823T070159Z/`
  `train_snapshot.py` captured the docstring-edited file (`98c4c6eb...`)
  rather than the file that actually ran. The authoritative record of what ran
  is the ledger's `train_sha256 = d53b4e23...`, and the two differ by that one
  docstring paragraph and nothing else. The launcher-written snapshot is left
  untouched, as locked artifacts are not edited; `train.py` in the working tree
  has been restored to the shipped baseline.

### Experiment 1 result (`run_001_20260823T070159Z`, 1282.3 s)

`reward = -1000000` (sentinel). `eligible = false`.
Diagnostic leader score before eligibility `-12.657`
= robust selection LCB `-11.845` minus repeat increment range `0.813`.
Mean source increment gap `1.713` months.

Failed gates (2 of 12):

- `all_repeat_genomic_increment_positive`
- `all_repeat_genomic_value_at_least_clinical`

Passed: alignment positive, value at least best constant, C-index drop
<= 0.03, seed agreement >= 0.85, seed tau correlation >= 0.50, genomic VIMP
>= 0.01, gene Jaccard >= 0.10, nontrivial benefit fraction >= 0.10, raw
propensity overlap >= 0.80, IPTW ESS >= 0.30n.

Per-repeat increments `A60(C+G) - A60(clinical)`:

```text
repeat_1  increment=-4.008  selection_lcb=-11.845  ci95=[-9.256, 1.131]  mean=-4.000 sd=2.686
repeat_2  increment=-3.195  selection_lcb=-10.938  ci95=[-8.940, 2.187]  mean=-3.212 sd=2.868
range=0.813  worst_lcb=-11.845
```

Source-specific OOF increments (diagnostic only, not reward):

```text
repeat_1  former_train=-4.822  former_validation=-1.571  gap=3.251
repeat_2  former_train=-3.151  former_validation=-3.326  gap=0.175
```

Both sources are negative in both repeats, so the deficit is not a
source-composition artifact.

C-index (separate matched prognostic forest, secondary):

```text
repeat_1  clinical=0.6806  C+G=0.6836   (+0.003)
repeat_2  clinical=0.6869  C+G=0.6635   (-0.023)
```

Gene stability: mean pairwise Jaccard `0.136` across the 8 outer-fold
selections. Full-development selection:
`TGFBR3, SMAD7, SMAD5, BMPR1A, TBL1XR1, NRIP1, CHD9, MEF2C`. The most
frequently reselected genes across folds were a branched-chain/valine
catabolism block (`ACAD8`, `HIBCH`, `MCCC1` at 5/8; `ALDH6A1` at 4/8), while
the full-development panel is dominated by a TGF-beta/SMAD block. Fold-level
membership therefore moves between two distinct pathway blocks.

#### Mandatory CSF diagnostics

```text
act_mechanism: W supplied separately as treatment; ACT is absent from X
rsf_act_split/path/terminal: NA by design (undefined for causal survival
  forests; ACT is W, not an X feature)

development_csf_cg: seed_agreement=0.9924; seed_tau_correlation=0.9986;
  genomic_vimp_fraction=0.8067; benefit_iqr=2.788; median_abs_benefit=1.664;
  nontrivial_fraction=0.9192; act_recommended=0.7215
development_csf_clinical: seed_agreement=0.9971; seed_tau_correlation=0.9999;
  genomic_vimp_fraction=0.0000; benefit_iqr=4.382; median_abs_benefit=2.124;
  nontrivial_fraction=0.9371; act_recommended=0.5169

repeat_1_csf_cg: seed_agreement=0.9942; seed_tau_correlation=0.9989;
  genomic_vimp_fraction=0.8138; benefit_iqr=3.157; median_abs_benefit=1.832;
  nontrivial_fraction=0.9381; act_recommended=0.6702
repeat_1_csf_clinical: seed_agreement=0.9971; seed_tau_correlation=0.9999;
  genomic_vimp_fraction=0.0000; benefit_iqr=4.265; median_abs_benefit=2.067;
  nontrivial_fraction=0.9420; act_recommended=0.5261

repeat_2_csf_cg: seed_agreement=0.9906; seed_tau_correlation=0.9984;
  genomic_vimp_fraction=0.7997; benefit_iqr=2.418; median_abs_benefit=1.496;
  nontrivial_fraction=0.9004; act_recommended=0.7727
repeat_2_csf_clinical: seed_agreement=0.9971; seed_tau_correlation=0.9999;
  genomic_vimp_fraction=0.0000; benefit_iqr=4.499; median_abs_benefit=2.181;
  nontrivial_fraction=0.9323; act_recommended=0.5077
```

Overlap: raw propensity in [0.05, 0.95] for 0.860 / 0.877 of rows; IPTW
effective sample size 360.1 / 384.5 against a 310.2 floor. Constant policies:
all-observation 45.88 / 45.90, all-ACT 45.58 / 45.46 months.

#### Interpretation

The hypothesis is rejected. Swapping the S-learner for an honest RMST causal
forest does not by itself make the inherited eight-raw-gene panel a positive
increment; the genomic policy is worse than the matched clinical policy in
both repeats by about 3-4 months of alignment.

The mechanism is visible and is not instability. Seed agreement is ~0.99 and
seed tau correlation ~0.999 in every panel, and the repeat range is only 0.81
months, so the negative increment is a stable property of the candidate, not
noise. What changes when the eight raw genes enter X is the shape of the CATE:

- Genomic variable importance takes 0.80 of the total, so with `mtry` = all
  and 18 clinical columns, the eight continuous genes absorb most splits:
  they are 31 percent of the candidate columns but 81 percent of the importance.
- The benefit IQR contracts from 4.38 to 2.79 months while the mean predicted
  benefit rises from 0.49 to 1.18, so tau shifts positive and flattens.
- ACT recommendation rises from 0.517 to 0.722, pushing the genomic policy
  toward the all-ACT constant (45.5 months) and away from the clinical
  policy's discrimination (alignment 4.52 versus 0.92).

So the genes are not adding effect modification; they are crowding out
clinical effect modification, the CSF analogue of v2's RSF ACT crowd-out. The
clinical comparator is strong, and eight noisy continuous columns dilute it.
This makes the number of genomic split candidates, not the selector's gene
identity, the first thing to test.


---

## Experiment 2 (prespecified before launch; ledger slot 3)

### Hypothesis

The eight-raw-gene deficit is caused by the *number* of genomic split
candidates, not by which genes the selector picks. Compressing the same eight
selected genes into two fit-only standardized module means should shrink
genomic variable importance well below experiment 1's 0.807, restore the
clinical CATE's discrimination (benefit IQR back toward the clinical 4.38
months and ACT recommendation back toward ~0.52), and raise the increment
`A60(C+G) - A60(clinical)` in both repeats.

Prespecified directional prediction, so the result can falsify the mechanism
rather than merely be described afterwards:

- If crowd-out is the cause, genomic VIMP falls (to roughly 0.2-0.5), benefit
  IQR rises above 2.79 months, ACT recommendation falls below 0.72, and both
  repeat increments rise above experiment 1's -4.008 and -3.195.
- If instead the eight genes carry no usable effect-modification signal at
  all, compression will leave the increments near zero-or-negative while VIMP
  drops, since averaging cannot create signal that is absent.
- If increments *fall* while VIMP drops, gene-specific direction matters and
  averaging cancels it, which is exactly the v2 warning; the next experiment
  would then attack regularization at fixed raw representation instead.

### Rationale

`program.md` puts raw-versus-module representation ahead of regularization and
ahead of any selector change, and experiment 1 supplies the specific failure
this addresses: genomic variable importance 0.807 from 31 percent of the
columns, against a clinical
comparator whose alignment is 4.52 months versus C+G's 0.92. With `mtry` fixed
at all features, eight continuous gene columns compete against 18 clinical
columns at every split, so the measured crowd-out is a representation problem
before it is a gene-identity problem. Two modules is the smallest compression
that still keeps more than one genomic axis, so gene-specific direction is
only partly averaged; `v2_findings.md` warns that full one-module compression
can cancel gene-specific effect modification.

### Exact change

In `train.py` `CANDIDATE`, only:

```text
name:           csf_pathway8_raw_honest_s100 -> csf_pathway8_module2_honest_s100
representation: raw    -> module
module_count:   8      -> 2
```

Unchanged: the selector and every one of its constants, `n_genes = 8`,
`num_trees = 1000`, `min_node_size = 5`, `sample_fraction = 0.50`,
`honesty_fraction = 0.50`, `alpha = 0.05`, `imbalance_penalty = 0.00`, the
zero-month threshold, seeds, folds, and bootstraps. The modules are built by
the locked `FeatureTransformer`: fit-only median imputation, fit-only
standardization, then `np.array_split` of the selected panel into two
contiguous groups, which follows the selector's pathway-block ordering. Total
treatment-effect features fall from 26 to 20, against the locked cap of 34.

### Red-line audit

1. Test untouched; loader still restricted to the two development CSVs. OK.
2. Only `train.py` edited; three `CANDIDATE` fields. OK.
3. Fit-only adaptation: module scaler and imputer are fit inside the locked
   transformer on fitting rows only. OK.
4. ACT stays `W`; modules are built from genes only, no ACT term, no
   gene-by-ACT product. OK.
5. Estimand untouched: horizon 60, threshold 0, LCB and range penalty. OK.
6. Learning and grading stay separate. OK.
7. No censored patient dropped. OK.
8. Selection unchanged and still fit-partition-only with Reactome structure;
   no hard-coded symbols; the module grouping is positional, not symbol-based. OK.
9. Identical geometry for clinical and C+G; the clinical forest is untouched
   by representation because it receives no genes. OK.
10. All mandatory CSF diagnostics will be recorded. OK.
11. C-index secondary, gate only. OK.
12. No probes; one full launcher slot. OK.
13. No metric shopping; eligible reward remains the leaderboard. OK.
14. Ledger slot 3 of 20; 30-minute wall; <=1000 trees per forest;
    28 <= 34 features; module_count 2 within the locked 1-4 range. OK.
15. Claims remain observational. OK.

### Prespecified decision rule

- Promote only on an eligible reward improvement.
- If ineligible, compare genomic VIMP, benefit IQR, ACT recommendation rate,
  and both increments against experiment 1 to decide which of the three
  branches above holds, then move to causal-forest regularization at the
  representation that scored better.

### Experiment 2 result (`run_003_20260823T072640Z`, 1209.8 s, ledger slot 3)

`reward = -1000000` (sentinel). `eligible = false`.
Diagnostic leader score before eligibility `-15.771`
= robust selection LCB `-14.034` minus repeat increment range `1.737`.
Mean source increment gap `2.033` months. Worse than experiment 1 on every
component, so `diagnostic_leader.txt` stays on `run_001`.

Failed gates (4 of 12, two more than experiment 1):

- `all_repeat_genomic_increment_positive`
- `all_repeat_genomic_value_at_least_clinical`
- `all_repeat_genomic_alignment_positive` (newly failing)
- `all_repeat_genomic_value_at_least_best_constant` (newly failing)

```text
repeat_1  increment=-5.458  selection_lcb=-14.034  ci95=[-11.068, -0.754]  mean=-5.478 sd=2.646
repeat_2  increment=-3.720  selection_lcb=-10.751  ci95=[ -8.907,  1.041]  mean=-3.746 sd=2.526
range=1.737  worst_lcb=-14.034
```

Source-specific OOF increments (diagnostic only):

```text
repeat_1  former_train=-5.372  former_validation=-5.714  gap=0.342
repeat_2  former_train=-2.788  former_validation=-6.511  gap=3.723
```

C-index: repeat_1 clinical `0.6806` vs C+G `0.6759`; repeat_2 clinical
`0.6869` vs C+G `0.6781`. Both within the 0.03 noninferiority gate.

Gene stability is unchanged by construction: Jaccard `0.136`, identical
selections to experiment 1, since only the representation of the selected
panel changed. Full-development modules:
`[TGFBR3, SMAD7, SMAD5, BMPR1A]` and `[TBL1XR1, NRIP1, CHD9, MEF2C]`.

#### Mandatory CSF diagnostics

```text
act_mechanism: W supplied separately as treatment; ACT is absent from X
rsf_act_split/path/terminal: NA by design (undefined for causal survival
  forests; ACT is W, not an X feature)

development_csf_cg: seed_agreement=0.9945; seed_tau_correlation=0.9996;
  genomic_vimp_fraction=0.5571; benefit_iqr=3.955; median_abs_benefit=2.143;
  nontrivial_fraction=0.9289; act_recommended=0.6465
development_csf_clinical: seed_agreement=0.9971; seed_tau_correlation=0.9999;
  genomic_vimp_fraction=0.0000; benefit_iqr=4.382; median_abs_benefit=2.124;
  nontrivial_fraction=0.9371; act_recommended=0.5169

repeat_1_csf_cg: seed_agreement=0.9955; seed_tau_correlation=0.9996;
  genomic_vimp_fraction=0.5623; benefit_iqr=4.197; median_abs_benefit=2.240;
  nontrivial_fraction=0.9304; act_recommended=0.6257
repeat_1_csf_clinical: seed_agreement=0.9971; seed_tau_correlation=0.9999;
  genomic_vimp_fraction=0.0000; benefit_iqr=4.265; median_abs_benefit=2.067;
  nontrivial_fraction=0.9420; act_recommended=0.5261

repeat_2_csf_cg: seed_agreement=0.9936; seed_tau_correlation=0.9996;
  genomic_vimp_fraction=0.5519; benefit_iqr=3.714; median_abs_benefit=2.046;
  nontrivial_fraction=0.9275; act_recommended=0.6673
repeat_2_csf_clinical: seed_agreement=0.9971; seed_tau_correlation=0.9999;
  genomic_vimp_fraction=0.0000; benefit_iqr=4.499; median_abs_benefit=2.181;
  nontrivial_fraction=0.9323; act_recommended=0.5077
```

Overlap and constants are identical to experiment 1 by construction (the
locked nuisance pipeline does not depend on the candidate): overlap 0.860 /
0.877, IPTW ESS 360.1 / 384.5, all-observation 45.88 / 45.90, all-ACT
45.58 / 45.46.

#### Interpretation: third branch of the prespecified prediction

Every *mechanistic* prediction of the crowd-out hypothesis held, and the
*policy* prediction failed:

```text
                        exp 1 (raw 8)   exp 2 (module 2)   predicted direction
genomic vimp fraction        0.807           0.557         down    -> held
benefit IQR (months)         2.788           3.955         up      -> held
ACT recommended              0.722           0.646         down    -> held
repeat increments      -4.008 / -3.195  -5.458 / -3.720    up      -> FAILED
```

So reducing genomic split competition does restore the CATE's shape, but the
restored shape is not a better policy. This is the third branch written down
before the run: gene-specific direction carries the little information the
panel has, and averaging four genes into a module mean cancels it. Two module
means still absorb 0.557 of variable importance from only 2 of 20 columns,
which also confirms that continuous columns dominate GRF's split-frequency
importance largely through split opportunity rather than signal.

The deeper reading is that the panel's problem is not representation. In both
experiments the genomic policy converges to recommending ACT for roughly
two-thirds of patients while the all-ACT constant (45.5 months) is *below* the
all-observation constant (45.9), and the genomic alignment collapses from the
clinical 4.5-5.1 months to about 0. Genes are not reordering who benefits;
they are adding variance that erases the clinical ordering and biasing the
policy toward blanket treatment.

Per the prespecified rule, representation reverts to raw, which scored better
on every component, and the next experiment attacks causal-forest
regularization.

---

## Experiment 3 (prespecified before launch; ledger slot 4)

### Hypothesis

The genomic deficit is driven by leaf-level overfitting, not by representation
and not yet by gene identity. With `min_node_size = 5`, honesty at 0.50 and
`sample_fraction = 0.50`, each tree splits on roughly 194 of the 775 fitting
rows and estimates each leaf's 60-month RMST contrast from as few as five
honest observations under censoring. Continuous gene columns can carve that
space almost arbitrarily, so the genomic forest manufactures tau variation that
does not reproduce out of fold. Raising `min_node_size` to 25 should damp that
manufactured variation and raise both repeat increments toward zero.

Prespecified directional predictions:

- Genomic ACT recommendation moves from 0.72 toward the clinical 0.52, and the
  genomic policy stops converging on blanket treatment.
- Both repeat increments rise above experiment 1's -4.008 and -3.195.
- If increments become positive in both repeats with the other ten gates still
  passing, the candidate is eligible and is promoted by the launcher.
- If increments rise but stay negative, the deficit is partly overfitting and
  partly gene identity, and the remaining gap is attributable to *which* genes
  the selector picks. Selector logic is then the next target.
- If increments do not rise at all, overfitting is not the mechanism, the
  panel carries no usable effect modification as selected, and regularization
  is retired as a lever.

Note the honest limit of this lever: as `min_node_size` grows, both forests
shrink toward a constant tau and both alignments go to zero, so regularization
alone can move the increment toward zero but cannot manufacture a positive one.
This experiment is therefore diagnostic about *how much* of the 3-4 month
deficit is noise, and only incidentally a candidate for promotion.

### Rationale

`program.md` orders causal-forest regularization directly after the
raw-versus-module comparison and before selector complexity, and experiments 1
and 2 have now settled the representation question in favour of raw. The
specific failure this addresses is measured, not assumed: in both experiments
the genomic alignment collapses from the clinical 4.5-5.1 months to about 0
while seed agreement stays at 0.99, which is the signature of variance that is
stable across seeds (all seeds see the same overfit fitting partition) but does
not carry out of fold. Leaf size is the one locked-geometry knob that directly
controls it, and it applies identically to the clinical and C+G forests, so the
comparator cannot be weakened.

### Exact change

In `train.py` `CANDIDATE`, revert experiment 2's representation and change one
forest parameter:

```text
name:            csf_pathway8_module2_honest_s100 -> csf_pathway8_raw_node25
representation:  module -> raw
module_count:    2      -> 8      (raw requires module_count == n_genes)
min_node_size:   5      -> 25
```

Unchanged: selector and all its constants, `n_genes = 8`, `num_trees = 1000`,
`sample_fraction = 0.50`, `honesty_fraction = 0.50`, `alpha = 0.05`,
`imbalance_penalty = 0.00`, zero-month threshold, seeds, folds, bootstraps.
`min_node_size = 25` is inside the locked `[3, 40]` range.

### Red-line audit

1. Test untouched. OK.
2. Only `train.py` edited; four `CANDIDATE` fields, no selector logic. OK.
3. Fit-only adaptation unchanged. OK.
4. ACT remains `W`; no ACT feature, no gene-by-ACT product. OK.
5. Estimand untouched: horizon 60, threshold 0, LCB, range penalty. OK.
6. Learning and grading separate. OK.
7. No censored patient dropped. OK.
8. Selection unchanged; no hard-coded symbols. OK.
9. `min_node_size` applies identically to the clinical and C+G forests and to
   the secondary prognostic forest, so the comparator is not weakened; if
   anything the clinical policy also benefits, which makes a positive increment
   harder, not easier. OK.
10. All mandatory CSF diagnostics will be recorded. OK.
11. C-index secondary, gate only. OK.
12. No probes; one full launcher slot. OK.
13. No metric shopping. OK.
14. Ledger slot 4 of 20; 30-minute wall; 1000 trees; 26 <= 34 features. OK.
15. Claims remain observational. OK.

### Prespecified decision rule

- Promote only on an eligible reward improvement.
- Otherwise classify the outcome into one of the three branches above and, if
  the deficit is attributable to gene identity, move to selector logic: the
  first target is the member-ordering rule, which currently ranks genes inside
  a chosen pathway by expression spread and uses the DR co-selection frequency
  only as a coarse median gate, so the widest-spread rather than the most
  benefit-associated members enter the panel.

### Experiment 3 result (`run_004_20260823T075139Z`, 1177.0 s, ledger slot 4)

`reward = -1000000` (sentinel). `eligible = false`.
Diagnostic leader score before eligibility `-4.745`
= robust selection LCB `-2.885` minus repeat increment range `1.861`.
Mean source increment gap `1.854` months. This is numerically the best score so
far, so the launcher moved `diagnostic_leader.txt` to `run_004`; the section
below argues that this leader is degenerate and must not steer the search.

Failed gates (3 of 12):

- `all_repeat_genomic_increment_positive`
- `all_repeat_genomic_alignment_positive`
- `all_repeat_genomic_value_at_least_best_constant`

`all_repeat_genomic_value_at_least_clinical` now passes, but only because the
clinical policy itself collapsed.

```text
repeat_1  increment= 0.000  selection_lcb=  0.000  ci95=[0.000, 0.000]  mean= 0.000 sd=0.000
repeat_2  increment= 1.861  selection_lcb= -2.885  ci95=[-1.608, 5.447] mean= 1.864 sd=1.788
range=1.861  worst_lcb=-2.885
```

Repeat 1's exact zero is not a rounding artifact: at `min_node_size = 25` the
clinical and C+G forests produced *identical recommendation sets* (both 0.7505
ACT), so every bootstrap draw differenced to exactly zero.

Source-specific OOF increments (diagnostic only):

```text
repeat_1  former_train= 0.000  former_validation= 0.000  gap=0.000
repeat_2  former_train= 0.932  former_validation= 4.639  gap=3.707
```

C-index: repeat_1 clinical `0.6788` vs C+G `0.6777`; repeat_2 clinical
`0.6831` vs C+G `0.6584`. Within the 0.03 gate.

Gene selection is unchanged from experiments 1 and 2 (Jaccard `0.136`, same
panel), since only forest geometry moved.

#### Mandatory CSF diagnostics

```text
act_mechanism: W supplied separately as treatment; ACT is absent from X
rsf_act_split/path/terminal: NA by design (undefined for causal survival
  forests; ACT is W, not an X feature)

development_csf_cg: seed_agreement=1.0000; seed_tau_correlation=0.9928;
  genomic_vimp_fraction=0.8249; benefit_iqr=1.193; median_abs_benefit=2.021;
  nontrivial_fraction=0.8752; act_recommended=0.8752
development_csf_clinical: seed_agreement=1.0000; seed_tau_correlation=0.9990;
  genomic_vimp_fraction=0.0000; benefit_iqr=1.365; median_abs_benefit=1.281;
  nontrivial_fraction=0.8752; act_recommended=0.7500

repeat_1_csf_cg: seed_agreement=1.0000; seed_tau_correlation=1.0000;
  genomic_vimp_fraction=0.7968; benefit_iqr=0.923; median_abs_benefit=2.352;
  nontrivial_fraction=0.7505; act_recommended=0.7505
repeat_1_csf_clinical: seed_agreement=1.0000; seed_tau_correlation=1.0000;
  genomic_vimp_fraction=0.0000; benefit_iqr=0.760; median_abs_benefit=1.431;
  nontrivial_fraction=0.7505; act_recommended=0.7505

repeat_2_csf_cg: seed_agreement=1.0000; seed_tau_correlation=0.9855;
  genomic_vimp_fraction=0.8530; benefit_iqr=1.462; median_abs_benefit=1.689;
  nontrivial_fraction=1.0000; act_recommended=1.0000
repeat_2_csf_clinical: seed_agreement=1.0000; seed_tau_correlation=0.9979;
  genomic_vimp_fraction=0.0000; benefit_iqr=1.970; median_abs_benefit=1.131;
  nontrivial_fraction=1.0000; act_recommended=0.7495
```

Overlap and constants are unchanged by construction: overlap 0.860 / 0.877,
IPTW ESS 360.1 / 384.5, all-observation 45.88 / 45.90, all-ACT 45.58 / 45.46.

#### Interpretation: the lever works by destroying both policies

The increments did rise, exactly as predicted, and the prediction about ACT
recommendation was wrong in an informative way: recommendation went *up*, not
down, to 0.875 pooled and to 1.000 in repeat 2, where the genomic policy
recommends ACT for every single patient.

```text
                          exp 1 (node 5)      exp 3 (node 25)
clinical benefit IQR          4.38                1.36
clinical alignment       5.091 / 3.948      0.465 / -2.301
clinical ACT recommended      0.517               0.750
C+G alignment            1.083 / 0.753      0.465 / -0.440
increments              -4.008 / -3.195     0.000 / +1.861
```

At `min_node_size = 25` a tree has roughly eight leaves over ~194 splitting
rows, so tau collapses toward a nearly constant positive value in both arms.
Both forests then recommend ACT for three-quarters to all patients, and both
policy values fall to the all-ACT constant (45.5 months), which is *below* the
all-observation constant (45.9). The increment improves only because the
clinical comparator was flattened from a genuinely discriminating policy
(value 48.3 versus 45.9 for all-observation) to a near-blanket one.

This is why `run_004` is a degenerate diagnostic leader. Red line 13 already
forbids it from nominating a test candidate; recorded here explicitly, it must
also not steer the next hypothesis, because its score comes from removing
information rather than adding it. The scientifically informative leader
remains `run_001`.

Three conclusions carry forward:

1. Fine leaves are *necessary* for the clinical policy. At `min_node_size = 5`
   the clinical forest reaches OOF alignment 4.5-5.1 months and value 47.7-48.3
   against constants near 45.5-45.9, so its CATE is real out-of-fold signal, not
   an artifact. Regularizing it away is not progress.
2. Regularization is retired as a lever for this deficit. It cannot manufacture
   a positive increment; between node 5 and node 25 it only trades clinical
   discrimination for a smaller gap.
3. What is left is gene identity. Across three runs the panel has been held
   fixed while representation and leaf size moved through their useful ranges,
   and the genomic policy never once ordered patients better than the clinical
   one. The next experiment therefore changes selector logic, as the experiment
   3 decision rule prespecified.

---

## Experiment 4 (prespecified before launch; ledger slot 5)

### Hypothesis

The panel fails because of *which members* of the selected pathways enter it.
The inherited selector uses its two statistics in the opposite of the useful
order: within a chosen pathway it keeps every gene at or above the median DR
co-selection count and then orders by **expression spread**, so the genes that
actually enter the panel are the widest-spread members of a benefit-associated
pathway rather than its most benefit-associated members. Widest spread is also
maximal split opportunity, which is why eight genes take 0.81 of variable
importance while adding no policy value. Ordering members by DR co-selection
frequency instead, with spread demoted to a tie-break, should produce a panel
whose genes are chosen for treatment-benefit association, and should raise both
repeat increments above experiment 1's -4.008 and -3.195.

Prespecified directional predictions:

- Genomic variable importance falls below experiment 1's 0.807, because
  benefit-ranked members are not selected for having the most split points.
- The selected panel differs from experiment 1's on most folds; if it does not
  differ at all, the hypothesis is untestable as posed and the run is reported
  as such rather than reinterpreted.
- Both repeat increments rise above -4.008 and -3.195. If they become positive
  in both repeats with the other gates passing, the candidate is eligible.
- The identified risk is gene stability: `gene_top_counts` is re-estimated
  inside every fitting partition, whereas expression spread is nearly constant
  across partitions, so the Jaccard could fall from 0.136 toward the 0.10 gate.
  A Jaccard failure with improved increments would be a real finding about the
  stability-versus-relevance trade-off, and would point to raising the
  aggregation (more subsamples or a tighter `STABILITY_TOP_K`) rather than
  reverting to spread ordering.

### Rationale

Three runs have now held the panel fixed while moving everything around it.
Representation moved through raw and modules (experiment 2) and leaf size moved
through its useful range (experiment 3), and in every configuration where the
clinical CATE survived, the genomic policy ordered patients worse than the
clinical one. `program.md` allows selector logic once the representation and
regularization branches are settled, and the experiment 3 decision rule named
this exact target in advance. This is also a swap of two statistics the
selector already computes, not new selector machinery: no policy trees, no SHAP,
no extra tuning surface.

### Exact change

In `train.py` `stability_select_genes`, replace the median-gate-then-spread
member ordering with a direct benefit ordering:

```text
before: keep members with gene_top_counts >= pathway median, order that group
        by descending expression spread, then append the rest by spread
after:  order all members by descending gene_top_counts, tie-broken by
        descending expression spread, then by gene symbol
```

`CANDIDATE` returns to the experiment 1 geometry exactly:
`representation = raw`, `module_count = 8`, `min_node_size = 5`, name
`csf_pathway8_raw_drorder`. The detectability filter, `MIN_PATHWAY_MEMBERS`,
`TOP_PATHWAYS`, `MAX_GENES_PER_PATHWAY`, `STABILITY_SUBSAMPLES`,
`STABILITY_TOP_K`, `WINSOR_PERCENT`, the pathway-level co-selection rule, and
the trailing mean-score fallback are all unchanged, so experiment 4 differs
from experiment 1 in the member-ordering rule alone.

### Red-line audit

1. Test untouched; selector still sees only the fitting partition. OK.
2. Only `train.py` edited. OK.
3. Fit-only: `gene_top_counts` and `gene_spread` are both computed inside the
   fitting partition, from its own half-samples; no assessment row, outcome, or
   covariate summary is involved. OK.
4. ACT stays `W`; the DR pseudo-outcome is the locked cross-fitted quantity and
   enters selection only, never X, and no gene-by-ACT product is built. OK.
5. Estimand untouched. OK.
6. Selection scores never become the reward; grading stays with the locked
   assessment-fold AIPW scores. OK.
7. No censored patient dropped; winsorization still replaces rather than
   removes non-finite pseudo-outcomes. OK.
8. Genes are targeted at effect modification by construction: the ordering key
   is now the treatment-benefit co-selection count itself. No hard-coded
   symbols, no patient indices, no assessment ranking. OK.
9. Forest geometry identical for clinical and C+G; unchanged from experiment 1. OK.
10. All mandatory CSF diagnostics will be recorded. OK.
11. C-index secondary. OK.
12. No probes; one full launcher slot. OK.
13. No metric shopping; `run_004` is explicitly not steering this. OK.
14. Ledger slot 5 of 20; 30-minute wall; 1000 trees; 26 <= 34 features. OK.
15. Claims remain observational. OK.

### Prespecified decision rule

- Promote only on an eligible reward improvement.
- If increments improve but stay negative, the selector direction is right and
  the next experiment tightens the same axis (panel width or aggregation),
  not a new selector family.
- If increments do not improve, the DR benefit signal itself does not identify
  effect modifiers at this sample size, and the remaining slots go to reporting
  `NO_ELIGIBLE_CANDIDATE` honestly rather than to searching for a lucky panel.

### Experiment 4 result (`run_005_20260823T081451Z`, 1181.2 s, ledger slot 5)

`reward = -1000000` (sentinel). `eligible = false`.
Diagnostic leader score before eligibility `-12.376`
= robust selection LCB `-11.443` minus repeat increment range `0.933`.
Mean source increment gap `1.219` months. Better than experiment 1's `-12.657`
on every component; still behind experiment 3's degenerate `-4.745`, which the
previous section rules out as a guide.

Failed gates (2 of 12, back to experiment 1's pair and no others):

- `all_repeat_genomic_increment_positive`
- `all_repeat_genomic_value_at_least_clinical`

```text
repeat_1  increment=-2.400  selection_lcb=-11.443  ci95=[-8.693, 3.595]  mean=-2.445 sd=3.156
repeat_2  increment=-1.467  selection_lcb= -9.403  ci95=[-7.532, 4.206]  mean=-1.504 sd=2.973
range=0.933  worst_lcb=-11.443
```

Source-specific OOF increments (diagnostic only):

```text
repeat_1  former_train=-1.798  former_validation=-4.201  gap=2.403
repeat_2  former_train=-1.475  former_validation=-1.441  gap=0.034
```

C-index: repeat_1 clinical `0.6806` vs C+G `0.6836`; repeat_2 clinical
`0.6869` vs C+G `0.6661`. Within the 0.03 gate.

Gene stability: Jaccard `0.1027`, down from `0.1362` and only just above the
0.10 gate, exactly the risk written into the prespecification. Full-development
panel `TGFBR3, SMAD7, SMAD4, SMAD5, SIRT1, CHD9, EP300, TBL1XR1`. Most
frequently reselected across the eight folds: `HIBCH` 5, `MCCC1` 4, `ACAD8` 3,
`ACADSB` 3, `AUH` 3, `TGFBR3` 3. Fold panels are drawn from three recurring
blocks: branched-chain amino-acid catabolism (`HIBCH`, `MCCC1`, `ACAD8`,
`ACADSB`, `AUH`, `ALDH6A1`), TGF-beta/SMAD (`TGFBR3`, `SMAD1/4/5/7`, `BMP2`),
and PI3K/phosphoinositide (`PIK3CB`, `PIK3R1`, `PTPN13`, `RAB4A`, `MTMR2`).

#### Mandatory CSF diagnostics

```text
act_mechanism: W supplied separately as treatment; ACT is absent from X
rsf_act_split/path/terminal: NA by design (undefined for causal survival
  forests; ACT is W, not an X feature)

development_csf_cg: seed_agreement=0.9953; seed_tau_correlation=0.9991;
  genomic_vimp_fraction=0.8214; benefit_iqr=3.190; median_abs_benefit=1.828;
  nontrivial_fraction=0.9328; act_recommended=0.6958
development_csf_clinical: seed_agreement=0.9971; seed_tau_correlation=0.9999;
  genomic_vimp_fraction=0.0000; benefit_iqr=4.382; median_abs_benefit=2.124;
  nontrivial_fraction=0.9371; act_recommended=0.5169

repeat_1_csf_cg: seed_agreement=0.9971; seed_tau_correlation=0.9992;
  genomic_vimp_fraction=0.8212; benefit_iqr=3.386; median_abs_benefit=1.943;
  nontrivial_fraction=0.9352; act_recommended=0.6809
repeat_1_csf_clinical: seed_agreement=0.9971; seed_tau_correlation=0.9999;
  genomic_vimp_fraction=0.0000; benefit_iqr=4.265; median_abs_benefit=2.067;
  nontrivial_fraction=0.9420; act_recommended=0.5261

repeat_2_csf_cg: seed_agreement=0.9936; seed_tau_correlation=0.9989;
  genomic_vimp_fraction=0.8217; benefit_iqr=2.994; median_abs_benefit=1.712;
  nontrivial_fraction=0.9304; act_recommended=0.7108
repeat_2_csf_clinical: seed_agreement=0.9971; seed_tau_correlation=0.9999;
  genomic_vimp_fraction=0.0000; benefit_iqr=4.499; median_abs_benefit=2.181;
  nontrivial_fraction=0.9323; act_recommended=0.5077
```

Overlap and constants unchanged by construction: overlap 0.860 / 0.877, IPTW
ESS 360.1 / 384.5, all-observation 45.88 / 45.90, all-ACT 45.58 / 45.46.

#### Interpretation: the selector direction is right, the shift is not fixed

Against experiment 1, which this run differs from only in the member-ordering
rule:

```text
                          exp 1            exp 4
repeat increments    -4.008 / -3.195  -2.400 / -1.467   +1.6 / +1.7 months
C+G alignment         1.083 / 0.753    2.691 / 2.481
C+G value (months)   46.272 / 46.055  47.076 / 46.919
repeat range              0.813            0.933
source gap                1.713            1.219
gene Jaccard              0.136            0.103
genomic vimp              0.807            0.821
```

Ordering pathway members by treatment-benefit co-selection instead of
expression spread recovers about 1.7 months of the deficit in both repeats and
raises the genomic policy's own value, while gene identity moves onto
biologically coherent blocks. The prediction that genomic variable importance
would fall was wrong: it stayed at 0.821. Importance is therefore not the
quantity that tracks policy quality here, and the earlier crowd-out reading
should be narrowed to say that split *opportunity* explains importance while
member *relevance* explains policy value.

What remains is a systematic positive shift in the genomic CATE, unchanged
across all four runs. The C+G forest's mean predicted benefit is 1.05-1.29
months against the clinical 0.47-0.52, so it recommends ACT for 0.68-0.71 of
patients against the clinical 0.51-0.53. Because the all-ACT constant (45.5) is
below the all-observation constant (45.9), that inflation is directly costly and
is the whole of the residual deficit: the genomic policy treats too many people
rather than ranking them wrongly.

Two candidate mechanisms for the shift, of which only the first is addressable
inside `train.py`: fewer effective observations per honest leaf once eight
continuous columns fragment the space, or residual confounding, since the
supplied `W.hat` is a locked clinical-only propensity and any gene-treatment
association beyond the clinical covariates is not orthogonalized away. The
locked nuisance pipeline is not editable, so the next experiment attacks
fragility of member choice, the axis the experiment 4 decision rule named.

---

## Experiment 5 (prespecified before launch; ledger slot 6)

### Hypothesis

Member choice is now benefit-driven but fragile: the Jaccard fell to 0.1027
against a 0.10 gate. The reason is that the ordering statistic is weaker than
the selector's own stated design. `stability_select_genes` documents that "a
gene scores only when it reaches the top of *both* rankings, so evidence that
rests on a few influential patients cannot promote it", but the code implements
that rule only for pathways: `pathway_counts` increments on the intersection of
the two complementary halves, while `gene_top_counts` increments once per half,
independently. A gene can therefore accumulate a high count from single-half
flukes. Applying the same both-halves intersection to the gene counts should
make member choice reproduce across fitting partitions (Jaccard back above
0.136) and concentrate the panel on genes whose benefit association survives
sample splitting, raising both repeat increments above -2.400 and -1.467.

Prespecified directional predictions:

- Gene Jaccard rises above experiment 4's 0.1027, and preferably above
  experiment 1's 0.1362.
- Both repeat increments rise above -2.400 and -1.467.
- The genomic ACT recommendation falls from 0.68-0.71 toward the clinical
  0.51-0.53, since the residual deficit is over-treatment rather than
  misranking.
- Failure mode to watch, and the reason the diagnostic below is prespecified:
  intersecting the two halves makes the count sparser, so if most members of a
  chosen pathway land on a count of zero the ordering silently reverts to the
  expression-spread tie-break and the run regresses toward experiment 1. The
  test is the selected panel itself: if it returns to experiment 1's
  spread-ordered membership, that is what happened, and the fix would be to
  widen `STABILITY_TOP_K` rather than to abandon the intersection.

### Rationale

The experiment 4 decision rule prespecified that a genuine improvement should
be followed by tightening the same axis, either panel width or aggregation, and
not by a new selector family. Panel width was considered and rejected as the
first move: with eight genes a Jaccard of 0.1027 corresponds to about 1.5 shared
genes per fold pair, and shrinking the panel to six would make the same overlap
score near 0.09, pushing an already-marginal gate below threshold for reasons
that have nothing to do with the science. Aggregation is the safer half of the
same axis, and this particular tightening costs no extra compute: it is the same
200 half-sample rankings, combined the way the surrounding code and the
docstring already say they should be.

### Exact change

In `train.py` `stability_select_genes`, accumulate the gene counts on the
intersection of the two complementary halves of each draw, matching the
existing pathway rule:

```text
before: for each half, gene_top_counts[top_k indices] += 1        (200 events)
after:  per draw, gene_top_counts[top_k in half A AND half B] += 1 (100 events)
```

Nothing else changes: `STABILITY_SUBSAMPLES = 100`, `STABILITY_TOP_K = 300`,
`WINSOR_PERCENT`, the detectability filter, `MIN_PATHWAY_MEMBERS`,
`TOP_PATHWAYS`, `MAX_GENES_PER_PATHWAY`, the pathway co-selection rule, the
mean-score fallback, and the whole `CANDIDATE` block including
`min_node_size = 5` and raw eight-gene representation. Experiment 5 differs
from experiment 4 in that one accumulation rule.

### Red-line audit

1. Test untouched. OK.
2. Only `train.py` edited. OK.
3. Fit-only: both halves are disjoint stratified subsets of the fitting
   partition; nothing from the assessment fold is read. OK.
4. ACT remains `W`; no ACT feature, no gene-by-ACT product. OK.
5. Estimand untouched. OK.
6. Selection statistics never enter the reward. OK.
7. No censored patient dropped. OK.
8. The count is still a fit-partition treatment-benefit signal restricted to
   Reactome pathways; the intersection only makes it stricter. No hard-coded
   symbols. OK.
9. Forest geometry unchanged and identical across arms. OK.
10. All mandatory CSF diagnostics will be recorded. OK.
11. C-index secondary. OK.
12. No probes; one full launcher slot; no new tuning surface, and the compute
    is identical to experiment 4. OK.
13. No metric shopping; the degenerate `run_004` leader is still excluded. OK.
14. Ledger slot 6 of 20; 30-minute wall; 1000 trees; 26 <= 34 features. OK.
15. Claims remain observational. OK.

### Prespecified decision rule

- Promote only on an eligible reward improvement.
- If the panel regresses to spread-ordered membership, widen `STABILITY_TOP_K`
  next, keeping the intersection.
- If stability improves but the increments do not, member choice is no longer
  the binding constraint and the remaining deficit is the CATE shift, which is
  attributable to the locked clinical-only `W.hat` and therefore not fixable
  from `train.py`; the program would then move toward an honest
  `NO_ELIGIBLE_CANDIDATE` report rather than spending slots on panels.

### Experiment 5 result (`run_006_20260823T083815Z`, 1180.0 s, ledger slot 6)

`reward = -1000000` (sentinel). `eligible = false`.
Diagnostic leader score before eligibility `-13.515`
= robust selection LCB `-12.465` minus repeat increment range `1.051`.
Mean source increment gap `1.019` months, the smallest so far. Worse than
experiment 4's `-12.376`; the hypothesis is rejected.

Failed gates (2 of 12): `all_repeat_genomic_increment_positive` and
`all_repeat_genomic_value_at_least_clinical`.

```text
repeat_1  increment=-3.640  selection_lcb=-12.465  ci95=[-9.675, 2.090]  mean=-3.670 sd=3.012
repeat_2  increment=-2.589  selection_lcb=-10.308  ci95=[-8.448, 2.962]  mean=-2.600 sd=2.924
range=1.051  worst_lcb=-12.465
```

Source-specific OOF increments: repeat_1 `-3.325` / `-4.582` (gap 1.257);
repeat_2 `-2.393` / `-3.174` (gap 0.781).

C-index: repeat_1 clinical `0.6806` vs C+G `0.6835`; repeat_2 clinical
`0.6869` vs C+G `0.6634`. Within gate.

Gene stability: Jaccard `0.1053`, essentially unchanged from experiment 4's
`0.1027` and still below experiment 1's `0.1362`. The full-development panel is
byte-identical to experiment 4's, and fold panels moved only slightly, mostly
by swapping one or two members inside the same block.

#### Mandatory CSF diagnostics

```text
act_mechanism: W supplied separately as treatment; ACT is absent from X
rsf_act_split/path/terminal: NA by design (undefined for causal survival
  forests; ACT is W, not an X feature)

development_csf_cg: seed_agreement=0.9911; seed_tau_correlation=0.9989;
  genomic_vimp_fraction=0.8096; benefit_iqr=3.062; median_abs_benefit=1.745;
  nontrivial_fraction=0.9192; act_recommended=0.6963
development_csf_clinical: seed_agreement=0.9971; seed_tau_correlation=0.9999;
  genomic_vimp_fraction=0.0000; benefit_iqr=4.382; median_abs_benefit=2.124;
  nontrivial_fraction=0.9371; act_recommended=0.5169

repeat_1_csf_cg: seed_agreement=0.9919; seed_tau_correlation=0.9990;
  genomic_vimp_fraction=0.8091; benefit_iqr=3.246; median_abs_benefit=1.816;
  nontrivial_fraction=0.9197; act_recommended=0.6567
repeat_1_csf_clinical: seed_agreement=0.9971; seed_tau_correlation=0.9999;
  genomic_vimp_fraction=0.0000; benefit_iqr=4.265; median_abs_benefit=2.067;
  nontrivial_fraction=0.9420; act_recommended=0.5261

repeat_2_csf_cg: seed_agreement=0.9903; seed_tau_correlation=0.9988;
  genomic_vimp_fraction=0.8100; benefit_iqr=2.879; median_abs_benefit=1.674;
  nontrivial_fraction=0.9188; act_recommended=0.7360
repeat_2_csf_clinical: seed_agreement=0.9971; seed_tau_correlation=0.9999;
  genomic_vimp_fraction=0.0000; benefit_iqr=4.499; median_abs_benefit=2.181;
  nontrivial_fraction=0.9323; act_recommended=0.5077
```

Overlap and constants unchanged by construction.

#### Interpretation: two findings, one of them the important one

The stated failure mode partly occurred. Requiring a gene to reach the top of
both complementary halves made the count sparse enough that ties fell through
to the expression-spread tie-break more often, and the Jaccard did not recover
(0.1027 -> 0.1053, against the predicted >0.136). The intersection is rejected
and the per-half count of experiment 4 is restored.

The more important finding is how *little* had to change for the policy to move
this much. Between experiments 4 and 5 the full-development panel is identical
and the fold panels differ by one or two members within the same biological
block, yet both repeat increments moved by more than a month
(-2.400 -> -3.640, -1.467 -> -2.589). Swapping `SMAD7` for `SMAD5` in one fold
is worth about the same as everything representation and regularization
achieved in experiments 2 and 3. That is a direct measurement of how weak the
genomic effect-modification signal is at n = 1034 with 152 treated patients:
panel-level conclusions here are not stable enough to support a claim about
which genes matter, whatever the increment eventually does.

Three quantities have now stayed fixed across all five runs and describe the
real obstacle:

```text
                     clinical        C+G (exp 1/4/5)
mean predicted tau   0.47-0.52       0.93-1.37
ACT recommended      0.51-0.53       0.66-0.74
genomic vimp         0               0.81-0.82
```

The genomic CATE is shifted positive by roughly 0.6-0.9 months regardless of
representation, leaf size, or member ordering, and with the policy threshold
locked at zero that shift is spent on treating an extra 15-20 percent of
patients while the all-ACT constant sits below all-observation. Reducing
genomic columns from eight to two in experiment 2 barely moved the shift
(1.18 -> 1.13), so it is not caused by the number of genomic split candidates.
The remaining candidate mechanism is that the supplied `W.hat` is a locked
clinical-only propensity, so any gene-treatment association beyond the clinical
covariates is not orthogonalized away inside gene-defined leaves. That is a
property of the locked nuisance pipeline, not of `train.py`.

A hazard for later experiments, recorded so no slot is wasted on it: the arena
validator accepts `sample_fraction` up to 0.70, but the locked bridge builds
every forest with `ci.group.size = 2` and honesty enabled, and `grf` rejects a
sampling fraction above 0.5 in that configuration. `sample_fraction` is
therefore treated as capped at 0.50 here.
