# Human experiment log

Scientific attempts are recorded automatically in `search_ledger.jsonl`. Add human interpretation here without changing locked code.

## 2026-08-30 — environment reconstruction before experiment 1

`integrity.verify_lock()` failed on every interpreter present on the machine.
`environment.json` pins python 3.10.9 / numpy 1.23.2 / pandas 2.2.3 /
scikit-learn 1.3.1 / scikit-survival 0.22.2; the `python3` on PATH was 3.12.8
with numpy 2.4.4, pandas 2.3.3, scikit-learn 1.8.0, scikit-survival 0.27.0, and
no other stack on the machine carried scikit-survival at all. The arena could
not have run a single attempt as found.

Reconstructed environment (all experiments in this arena use this interpreter,
launched explicitly by absolute path, never via PATH):

    interpreter  /Users/owen/Documents/GitHub/autosota-nsclc-adj-chemo/.venv-arena310/bin/python
    cpython      3.10.9 (uv-managed python-build-standalone, macos-aarch64)
    numpy        1.23.2
    pandas       2.2.3
    scikit-learn 1.3.1
    scikit-survival 0.22.2
    scipy        1.11.4
    matplotlib   3.7.5
    lifelines    0.27.8

`integrity.verify_lock()` passes under it.

Reproduction status of the committed smoke (`runs/smoke_20260830T075559Z`):
NOT reproduced. That smoke recorded reward -5.7841422004299625 with fold panels
OBS=12/ACT=11 and OBS=14/ACT=11. The rebuilt environment returns reward
-1000000.0 (ineligible: increment_positive, value_at_least_clinical and
cindex_drop all fail) with panels OBS=13/ACT=10 and OBS=16/ACT=7, and
score_before_eligibility -8.526934079580498.

What does reproduce exactly, to all printed digits:

    development_n 1034, events 465, act 152, eligible genes 64
    train_sha256      b93f3cee3858716bca8add447bc837c734e203f6cb37962b1081e3c118d2c65d
    validation_sha256 824689aeaa34984727f2551bd512214164d0559d05cbec5817b0f88c9adb15a9
    reactome_sha256   5d61f289a2400cddfbb3a3353829fd2284a360bbe50f2093b566c4b7bea93341
    clinical alignment_months 6.555514   clinical harrell_cindex 0.651838

So the whole deterministic nuisance layer (StratifiedKFold, logistic
propensity, both Cox fits, IPCW, AIPW) and the locked clinical T-learner agree
bit for bit; only the arm-specific screening path diverges. The five pinned
packages are therefore not sufficient to pin this arena's screening result —
something below them (BLAS/compiler build of the scikit-survival wheel is the
obvious candidate) differs from the machine that produced the committed smoke.
This is evidenced, not proven.

Two consecutive smokes under the rebuilt environment are byte-identical
(`runs/smoke_20260830T080630Z` and `runs/smoke_20260830T080856Z`), so the
environment is self-consistent and deterministic.

`search_ledger.jsonl` was empty when this session started: zero of twenty slots
consumed, no ledger rows, no runs. The experiment series therefore starts clean
inside this environment and no series is being spliced across environments. The
frozen bar in `budget.json`, run20_reward_to_beat = -3.5746550301982576, comes
from `reactome_tlearner_arena` run_020, produced on 2026-08-24 under the older
stack; it is the arena's declared bar and is used as such, but it was not
recomputed here and is not environment-matched.

Acceptance test for any future claim that this environment has been restored:

    .venv-arena310/bin/python -m rsf_tlearner_permimp_arena.run --smoke

must print fold panels `OBS=13 ACT=10` then `OBS=16 ACT=7` and record
score_before_eligibility -8.526934079580498 on the 64-gene pool.

## Experiment 001 — root universe screen (`run_001_20260830T092246Z`)

Candidate `root_universe_screen_t100d4_pool1024`. Screening 100 trees / depth 4
(OBS leaf 16 split 32, ACT leaf 8 split 16, `sqrt`) over all 8,647 development
Reactome genes; panels capped 32/32; `next_pool_size` 1,024 OBS / 480 ACT; final
T-learner at the locked clinical geometry with `max_features` 0.35. Elapsed
2,320 s of the 3,600 s per-run cap.

Result: **ineligible**, reward -1000000.0, score_before_eligibility
-15.09211083002802, robust_selection_lcb -13.359726106606182, repeat increment
range 1.7323847234218377.

Failed gates: observation_panel_jaccard, act_panel_jaccard,
all_repeat_genomic_increment_positive, all_repeat_genomic_alignment_positive,
all_repeat_genomic_value_at_least_clinical,
all_repeat_genomic_value_at_least_best_constant.

Per repeat:

    rep 1  increment -5.4471  lcb -13.3597  boot mean -5.4810 sd 2.7218
           clinical value 48.152 align 5.098 cindex 0.656
           genomic  value 44.632 align -0.349 cindex 0.635
    rep 2  increment -3.7147  lcb -12.0702  boot mean -3.7333 sd 2.5855
           clinical value 48.526 align 5.361 cindex 0.664
           genomic  value 45.839 align  1.646 cindex 0.641

Two findings.

**1. Thirty-two raw noise-selected genes actively destroy the policy.** Clinical
alignment 5.10 / 5.36 months collapses to -0.35 / 1.65 once the panel is added,
and the C-index drops 0.021 / 0.023. The increment is not merely non-positive,
it is large and negative in both repeats. Note the sibling
`reactome_tlearner_arena` reached +2.0 to +2.7 with 16 genes, but under
`representation: "module"` (genes compressed to a single standardized mean).
This arena hard-codes `FeatureTransformer(genes, "raw", len(genes))` in
`_fit_arm`, so no compression is available and every selected gene enters as its
own raw feature. Panel size is therefore the dominant risk, and it points the
other way from the stability gate (see below).

**2. The panel-stability gate is unreachable, and the prediction is exact.**
Before launching, the arena's own selection path was reproduced offline on the
same eight outer fit sets (fit rows and inner folds only; no assessment rows, no
AIPW value, no reward) and predicted obs Jaccard 0.0034 at cap 32. The run
reported 0.003401 for both arms — an exact match. The offline harness is
therefore a zero-cost, exact predictor of this gate.

Measured cross-fit reproducibility of the selection score (Spearman between two
independent fit sets), across every reachable corner of the candidate space:

    pool 8647  40t d6            obs  0.044   act -0.034
    pool 8647 100t d4            obs  0.051   act  0.087
    pool  512 200t d4            obs  0.133   act -0.079
    pool  512 500t d4            obs -0.045   act -0.027
    pool  512 200t d7            obs  0.008   act -0.011
    pool  256 200t d4 mtry sqrt  obs -0.158   act  0.057
    pool  256 200t d4 mtry 0.5   obs  0.104   act -0.252
    pool  256 200t d4 mtry 1.0   obs  0.092   act  0.120

Full 8-context panel Jaccard, root geometry over the whole universe:

    cap  4/8/16 -> 0.0000   cap 24 -> 0.0015   cap 32 -> 0.0034

and at pool 128 with 1000-tree screening (coverage 128/128, positive rate 0.28):

    cap 8 -> 0.0024   cap 16 -> 0.0046   cap 32 -> 0.0058

Structurally, the panel Jaccard is `s*cap/(2P)` when the cap binds and
`rho*s/(2-rho*s)` when it does not, with `P` the pool size, `rho` the
positive-score rate and `s` the cross-context pool commonality. Measured `s(K)`
from the root tables: 0.020 (K=64), 0.029 (128), 0.045 (256), 0.073 (512), 0.268
(1024), 0.638 (2048), 0.822 (4096), 1.000 (8647). Commonality only appears above
K ~ 1184, the point where the ranking exhausts covered genes and pads with
never-used genes that tie at exactly zero and sort by gene symbol. A pool can
only shrink through per-context ranking, which is precisely what destroys `s`;
there is no operator in the candidate space that shrinks a pool identically
across contexts. Both branches land at 0.003-0.006 at every setting, ~20x short
of the 0.10 gate. `max_features` -> 1.0, which makes splits data-driven rather
than a feature lottery, was the last plausible mechanism and does not move it:
panel Jaccard 0.051 / 0.044 / 0.032 for sqrt / 0.5 / 1.0, all at or *below* the
chance null of 0.067.

Root cause: `budget.json` locks `inner_pfi_folds: 2` and
`permutation_repeats: 2`, so each gene's importance is a mean of **four**
samples of a discrete policy-flip quantity, and the selection score then
subtracts one standard error of those four. The 0.10 stability gate is inherited
from `reactome_tlearner_arena`, where selection ran through
`select_genes_by_dr_benefit` — a *deterministic* ridge-residual screen, which
passed the gate in 16 of 20 runs. Pairing that gate with a stochastic four-sample
forest PFI is what makes it unreachable.

Deliberately not used: setting `next_pool_size` above the coverage boundary makes
the alphabetically ordered zero-block dominate the pool, driving `s` to ~0.9 and
letting a ~128-gene pool clear the gate at Jaccard ~0.11. That passes a
*stability* gate through gene-symbol tie-breaking rather than reproducible
biology, i.e. manufacturing the metric, and repo red line 8 forbids it. Both
pool sizes here (1,024 OBS against ~1,180 coverage; 480 ACT against ~493) were
chosen to stay strictly inside the screened region for exactly this reason.

Consequence: no candidate in this arena can produce an eligible run, so
`test_nominee.txt` cannot be populated and the frozen run-20 bar
(-3.5746550301982576) is unreachable behind the gate. Remaining slots go to the
arena's secondary target, `diagnostic_leader.txt` / `score_before_eligibility`,
which does not depend on the stability gate.

## Experiment 002 — panel-size dose response, first point (`run_002_20260830T094746Z`)

Candidate `child001_panel04_mf1_thr0`, parent `run_001_20260830T092246Z`,
inheriting the 1,024 OBS / 480 ACT pools. Panels capped 4/4, final T-learner at
the locked clinical geometry with `max_features` 1.0, threshold 0.0.
`next_pool_size` 512/240. Elapsed 1,115 s.

Result: **ineligible**, reward -1000000.0, score_before_eligibility
-14.567831, robust_selection_lcb -10.2479, repeat increment range 4.3200.

Panel Jaccard predicted offline before launch as 0.000000 / 0.000000; reported
0.000000 / 0.000000. Second exact match, on a different candidate and a
different pool. The offline predictor is confirmed.

    rep 1  increment -3.8188  lcb -10.2479  boot mean -3.8919 sd 2.0784
           clinical value 48.152 align 5.098 cindex 0.6561
           genomic  value 45.684 align 1.279 cindex 0.6429  seed agreement 0.993
    rep 2  increment +0.5012  lcb  -5.3079  boot mean +0.4560 sd 2.0022
           clinical value 48.526 align 5.361 cindex 0.6641
           genomic  value 47.715 align 5.862 cindex 0.6406  seed agreement 0.992

Failed gates: the two panel-Jaccard gates,
all_repeat_genomic_increment_positive, all_repeat_genomic_value_at_least_clinical.

Dropping from 32 raw genes per arm to 4 recovers most of the damage.
all_repeat_genomic_alignment_positive and
all_repeat_genomic_value_at_least_best_constant now **pass** (they failed in
experiment 001), and repeat 2's increment turns positive (+0.50 against -3.71).
Increment by panel size, both repeats:

    32 genes/arm (exp 001, mf 0.35)   -5.4471 / -3.7147
     4 genes/arm (exp 002, mf 1.00)   -3.8188 / +0.5012

Repeat 1 is the harder of the two in both runs, and it is what still blocks
increment_positive and value_at_least_clinical. Bootstrap sd also falls with
panel size (2.72/2.59 -> 2.08/2.00), as expected when fewer patients disagree
with the clinical policy.

Note on the secondary target: score_before_eligibility = min-repeat LCB minus the
repeat increment range, and both terms go to zero as the panel empties, so the
metric is maximized in the limit by a genomic model identical to the clinical
one. Whatever tops `diagnostic_leader.txt` at the small-panel end is therefore a
property of the metric, not evidence that a model is good. It is recorded here
as such.

## Experiment 003 — panel-size dose response, extreme point (`run_003_20260830T100400Z`)

Candidate `child001_panel01_mf1_thr0`, parent `run_001_20260830T092246Z`,
panels capped 1/1, otherwise identical to experiment 002. Elapsed 917 s.

Result: **ineligible**, reward -1000000.0, score_before_eligibility
-11.127635, robust_selection_lcb -9.3888, repeat increment range 1.7389.
Panel Jaccard predicted 0.000000 / 0.000000, reported 0.000000 / 0.000000 —
third exact match. New diagnostic leader.

    rep 1  increment -1.4997  lcb -6.3827  sd 1.6842
           clinical value 48.152 align 5.098 | genomic value 47.251 align 3.599 cindex 0.6596
           genomic split fraction  obs 0.3769  act 0.2932
    rep 2  increment -3.2385  lcb -9.3888  sd 1.8431
           clinical value 48.526 align 5.361 | genomic value 47.152 align 2.123 cindex 0.6584
           genomic split fraction  obs 0.3850  act 0.3260

The dose response is **not monotone**, and that is the finding. Increments by
panel size:

    32 genes/arm   -5.4471 / -3.7147
     4 genes/arm   -3.8188 / +0.5012
     1 gene /arm   -1.4997 / -3.2385

Four genes produced the only positive increment seen so far; one gene is
negative in both repeats. The scatter is consistent with noise, not with a
dose-response curve.

**Mechanism.** With a *single* raw gene per arm, that one gene absorbs
**37.7-38.5%** of OBS splits and 29.3-32.6% of ACT splits. Eighteen clinical
pretreatment columns, largely categorical or low-cardinality, compete against one
continuous expression value with ~1,000 distinct levels, and the log-rank split
search strongly prefers the high-cardinality continuous feature. So panel size is
not the operative variable: *any* raw gene admitted to the forest captures roughly
a third of the splits and displaces clinical structure. This arena hard-codes
`FeatureTransformer(genes, "raw", len(genes))` in `_fit_arm`, so there is no
module compression available to blunt it — the sibling arena's eligible runs used
`representation: "module"` precisely here.

Combining with the selector result, the picture is closed:

1. the PFI selector has ~zero cross-fit reproducibility, so selected genes are
   effectively random;
2. any raw gene captures ~1/3 of splits regardless of panel size;
3. therefore the genomic policy is a random perturbation of the clinical policy,
   and `alignment = mean((2*pi - 1) * delta)` falls under *any* random
   perturbation, so the expected increment is negative for every panel;
4. `all_repeat_genomic_increment_positive` can then only pass by luck in both
   repeats, and observed pairs are (-5.45, -3.71), (-3.82, +0.50), (-1.50, -3.24);
5. independently, the panel-stability gate is unreachable at 0.006-0.018 against
   a 0.10 threshold.

Discrimination is not the problem: genomic C-index 0.6596 / 0.6584 against
clinical 0.6561 / 0.6641, comfortably inside the 0.03 allowance. It is the
policy value that degrades.

## Experiment 004 — observation-only panel (`run_004_20260830T102405Z`)

Candidate `child001_obsonly04_mf1_thr0`, parent `run_001_20260830T092246Z`,
panels 4 OBS / 0 ACT, locked clinical geometry, `max_features` 1.0, threshold 0.0.
Elapsed 1,111 s.

Result: **ineligible**, reward -1000000.0, score_before_eligibility
-9.339350, robust_selection_lcb -7.0303, repeat increment range 2.3091.
New diagnostic leader. Panel Jaccard predicted obs 0.000000 / act 1.000000,
reported obs 0.000000 / act 1.000000 — fourth exact match.

    rep 1  increment -3.5178  lcb -7.0303  sd 1.1968
           clinical value 48.152 align 5.098 | genomic value 46.334 align 1.580 cindex 0.6437
    rep 2  increment -1.2087  lcb -5.1390  sd 1.3183
           clinical value 48.526 align 5.361 | genomic value 47.839 align 4.152 cindex 0.6462
    genomic split fraction  obs 0.644 / 0.649,  act 0.000

Only **three** gates now fail, down from six in experiment 001: the OBS panel
Jaccard, increment_positive, and value_at_least_clinical. Emptying the ACT panel
makes the ACT Jaccard gate pass at exactly 1.0, because `_nonempty_jaccard`
filters empty panels and `_pairwise_jaccard([])` returns 1.0 — the gate is
vacuous for an arm carrying no genes.

The substantive finding is the **variance decomposition**. Zeroing the ACT panel
roughly halves the paired-bootstrap sd at the same OBS panel size:

    4 OBS / 4 ACT (exp 002)   sd 2.0784 / 2.0022
    4 OBS / 0 ACT (exp 004)   sd 1.1968 / 1.3183

The ACT arm is fit on 152 patients against the OBS arm's 882, so ACT-arm genes
contribute most of the instability in the estimated benefit and hence in the
policy. Increments improve correspondingly (-3.82/+0.50 -> -3.52/-1.21 in level,
with far less spread), and score_before_eligibility improves from -14.568 to
-9.339. This is a real modelling result and not a metric artifact: it says that
under a T-learner on this cohort, genomic terms in the small treated arm cost
more in variance than they return in signal.

Running summary:

    run  panels obs/act  mf    increments          score     sf_obs sf_act  #gates failed
    001  32 / 32         0.35  -5.447, -3.715      -15.092   0.878  0.879   6
    002   4 /  4         1.00  -3.819, +0.501      -14.568   0.647  0.587   4
    003   1 /  1         1.00  -1.500, -3.239      -11.128   0.381  0.310   4
    004   4 /  0         1.00  -3.518, -1.209       -9.339   0.644  0.000   3

## Experiment 005 — minimal perturbation (`run_005_20260830T104022Z`)

Candidate `child001_obsonly01_mf1_thr0`, parent `run_001_20260830T092246Z`,
panels 1 OBS / 0 ACT, locked clinical geometry in both arms, threshold 0.0.
Elapsed 916 s.

Result: **ineligible**, reward -1000000.0, score_before_eligibility
-5.308255, robust_selection_lcb -4.5904, repeat increment range 0.7179.
New diagnostic leader. Panel Jaccard predicted obs 0.000000 / act 1.000000,
reported identically — fifth exact match.

    rep 1  increment -0.8238  lcb -4.4490  sd 1.3507
           clinical value 48.152 | genomic value 47.735 align 4.274 cindex 0.6579
    rep 2  increment -0.1060  lcb -4.5904  sd 1.4872
           clinical value 48.526 | genomic value 48.422 align 5.255 cindex 0.6621
    genomic split fraction obs 0.377 / 0.385, act 0.000
    nontrivial benefit fraction 0.993 / 0.990

This is the **minimum-perturbation configuration reachable in this arena**, and
worth stating precisely why. With `max_panel_genes["act"] = 0` and
`tlearner["act"]` set to exactly the locked clinical ACT parameters
(1000/7/12/24, mf 1.0), the genomic ACT arm in `_fit_predict` is fit on the same
rows, with the same geometry and the same seeds, and no genes — so it is
*identical* to the clinical ACT arm and contributes exactly zero difference.
The entire genomic-vs-clinical contrast is then one raw gene entering the OBS
arm. Increments -0.8238 / -0.1060 are the cost of that single gene, and repeat 2
is within noise of zero.

There is no lever that reduces this further. Shrinking `max_features` on the OBS
arm would cut the gene's share of splits, but `max_features` applies to the whole
candidate set, so it would simultaneously perturb the OBS arm away from the
mf 1.0 clinical baseline — trading a gene-shaped perturbation for a
geometry-shaped one. Any deviation from the locked clinical geometry adds
difference rather than removing it. score_before_eligibility -5.308 is therefore
close to the arena's attainable ceiling, and the residual gap to zero is the
irreducible cost of admitting one raw gene at 38% of splits.

Running summary:

    run  panels obs/act  increments          score     sd            #gates failed
    001  32 / 32         -5.447, -3.715      -15.092   2.72 / 2.59   6
    002   4 /  4         -3.819, +0.501      -14.568   2.08 / 2.00   4
    003   1 /  1         -1.500, -3.239      -11.128   1.68 / 1.84   4
    004   4 /  0         -3.518, -1.209       -9.339   1.20 / 1.32   3
    005   1 /  0         -0.824, -0.106       -5.308   1.35 / 1.49   3

## Experiment 006 — policy threshold (`run_006_20260830T105721Z`)

Candidate `child001_obsonly01_mf1_thr1`: experiment 005 exactly, with
`benefit_threshold_months` 1.0 instead of 0.0. Elapsed 919 s.

Result: **ineligible**, reward -1000000.0, score_before_eligibility
-4.363832, robust_selection_lcb -3.7478, repeat increment range 0.6160.
New diagnostic leader.

    rep 1  increment +0.7218  lcb -3.7478  sd 1.7833
           clinical value 48.152 align 5.098 act_frac 0.372
           genomic  value 48.519 align 5.820 act_frac 0.319
    rep 2  increment +0.1058  lcb -3.5971  sd 1.3494
           clinical value 48.526 align 5.361 act_frac 0.382
           genomic  value 48.508 align 5.467 act_frac 0.345

**Both increments are positive**, so `all_repeat_genomic_increment_positive`
passes for the first time. Only **two** gates now fail:
`all_repeat_genomic_value_at_least_clinical` and the OBS panel Jaccard.

I expected the threshold to hurt and it helped; recording that, because the
reasoning that produced the wrong prediction was itself wrong. The argument was
that the clinical policy at threshold 0 already scores 5.098 against 1.5 for
treating nobody, so it has real signal, and withdrawing marginal patients from
treatment should cost. What that misses is that the threshold is applied *only*
to the genomic policy while the clinical comparator stays at 0, so raising it
does not move the model — it moves the decision boundary of one arm of the
comparison. ACT-recommended fraction falls from 0.372/0.382 to 0.319/0.345, and
the patients dropped are those with small predicted benefit, where the AIPW
score delta is on average negative. The threshold is not a tie-breaker here, it
is a genuine policy lever, and it is the only candidate knob that improves
alignment without depending on gene signal.

The remaining value failure is **0.0176 months**. Genomic minus clinical value by
repeat: +0.3667 and -0.0176. Repeat 2 misses `value_at_least_clinical` by under
two hundredths of a month on a 60-month RMST scale.

Headroom for a larger threshold looks adequate: nontrivial benefit fraction
0.9429 / 0.9323 against a 0.10 gate, median |predicted benefit| 7.50 / 7.29
months, IQR ~14 months, so a threshold of 2 months still leaves the
`nontrivial_benefit_fraction` gate far from binding. Seed agreement 0.996/0.992
and seed benefit correlation 0.999 are both comfortable.

    run  panels  thr   increments          score      #gates failed
    001  32/32   0.0   -5.447, -3.715      -15.092    6
    002   4/ 4   0.0   -3.819, +0.501      -14.568    4
    003   1/ 1   0.0   -1.500, -3.239      -11.128    4
    004   4/ 0   0.0   -3.518, -1.209       -9.339    3
    005   1/ 0   0.0   -0.824, -0.106       -5.308    3
    006   1/ 0   1.0   +0.722, +0.106       -4.364    2

## Experiment 007 — threshold overshoot (`run_007_20260830T111411Z`)

Candidate `child001_obsonly01_mf1_thr2`: experiment 006 with
`benefit_threshold_months` 2.0. Elapsed 915 s.

Result: **ineligible**, score_before_eligibility -7.863992,
robust_selection_lcb -6.4043, range 1.4597. Three gates fail again
(increment_positive, value_at_least_clinical, OBS Jaccard). Not a new leader.

    rep 1  increment -0.2120  clinical value 48.1524  genomic 48.0433  gap -0.1091
    rep 2  increment -1.6717  clinical value 48.5260  genomic 47.7296  gap -0.7964
    act_recommended_fraction 0.315 / 0.315, nontrivial benefit 0.878 / 0.881

The threshold response is single-peaked near 1.0 month:

    threshold  increments          rep-2 value gap   score
      0.0      -0.824, -0.106      -0.104            -5.308
      1.0      +0.722, +0.106      -0.018            -4.364
      2.0      -0.212, -1.672      -0.796            -7.864

At 2.0 months the ACT-recommended fraction saturates at 0.315 in both repeats —
the threshold has pushed past the useful part of the predicted-benefit
distribution and is now withdrawing treatment from patients whose AIPW delta is
positive. Experiment 006 remains the diagnostic leader.
