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
