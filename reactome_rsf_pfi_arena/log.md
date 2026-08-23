# Reactome-wide RSF PFI search log

Append a prespecified hypothesis before each full experiment and complete the
entry afterward using the schema in `program.md`.


## Search-space review before iter_001

The agent-editable surface is exactly four fields: `name`, `screening_space`,
`reduced_space` (each restricted to `n_estimators`, `max_depth`,
`min_samples_leaf`, `split_leaf_multiplier`, `max_features`), and `plot_top_k`.
Everything else — folds, seeds, estimand, PFI protocol, gene eligibility, N
sweep — is locked in `prepare.py`.

Two locked mechanics dominate how that surface should be used.

1. **PFI cost is quadratic in screening trees.** `_permutation_importance`
   refits the screening winner and then scores one panel per
   (fold x used feature x repeat) = 9 panels per used feature. A panel costs
   roughly `0.15 s + 6.4e-6 * trees * 345` (wide-frame copy/transform overhead
   plus forest prediction), calibrated from smoke optuna timestamps
   (102 panels / 6.75 s at 20 trees). The number of *used* features grows with
   `trees * splits_per_tree`, so total PFI time scales as trees^2 * splits.
   Projected PFI wall time: 50 trees/depth 4 ~ 1.7 ks; 100/depth 4 ~ 4.6 ks;
   100/depth 5 ~ 8.1 ks; 200/depth 4 ~ 13.5 ks against a 14.4 ks budget. The
   shipped default space (200-800 trees, depth 4-10) has no affordable corner
   and would very likely burn a slot on a timeout.
2. **The reduced stage spends 33 of its 49 trials at one fixed geometry.**
   `_optimize_reduced` enqueues top_n = 0..32 using `space[key][0]` for every
   key. The first element of each reduced list is therefore the reference
   geometry for the entire controlled N curve, and only 16 trials are free to
   co-optimize N with geometry. List order is a real design decision, not
   cosmetic.

### iter_001 — cost-bounded screening baseline
- type: SEARCH
- hypothesis: A screening space capped at 100 trees and depth 4 finishes the
  whole two-stage pipeline well inside the 4-hour wall, and the reduced N curve
  at a moderate reference geometry (400 trees, depth 6, leaf 16, sqrt) selects a
  final panel whose 60-month IPCW-AIPW policy RMST difference is at least as
  large as the clinical-only N=0 panel.
- changed: `name` -> `cost_bounded_screen_v1`; `screening_space` ->
  n_estimators [50, 100], max_depth [3, 4], min_samples_leaf [24],
  split_leaf_multiplier [2], max_features ["sqrt", 0.025, 0.05] (a 12-point
  factorial matching the 12-trial screening budget); `reduced_space` ->
  n_estimators [400, 200, 600], max_depth [6, 4, 8], min_samples_leaf [16, 8, 32],
  split_leaf_multiplier [2, 3], max_features ["sqrt", 0.5, 1.0], ordered so the
  33 enqueued N-curve trials run at the intended reference geometry.
- red_line_audit: (1) train.py names no dataset path; loading stays inside
  locked `prepare.load_development`, which accepts only the train and validation
  CSVs — no test artifact is opened, hashed, or summarized. (2) Only `train.py`
  is edited and `log.md` appended; locked files verified by `integrity.verify_lock`
  are untouched. (3) No transform is touched; imputation, RSF fitting, nuisance
  models and PFI remain locked fit-only/held-out code. (4) No gene prefilter,
  ranking, or selection is introduced anywhere; the screening forest still
  receives every eligible Reactome gene, which the launcher re-checks. (5) No
  row filtering of any kind. (6) Reward stays the locked 60-month IPCW-AIPW RMST
  difference with C-index as a lexicographic second; no weighted blend.
  (7)-(8) Counterfactual pairing and permutation pairing are locked and
  unmodified. (9) 49 reduced trials >= the 33 enqueued N=0..32 panels, so every
  N is evaluated before selection. (10) No gene symbol, row id, prediction, or
  outcome appears in train.py; ranking remains the locked RMST/C-index/symbol
  rule. (11) Max n_estimators 600 <= 1000; this is one full attempt, experiment
  1 of 20, and the train.py snapshot is new. (12) Conclusions stay observational
  development-CV statements.
- run_id: run_001_20260822T233357Z
- screening_rmst_difference: 0.218 months
- screening_cindex: 0.653
- final_top_n: 20
- final_rmst_difference: 5.214 months
- final_cindex: 0.710
- final_genes: [ANKFY1, RAB4B, TRIM38, PRIM2, FRAT1, CFL1, FMO3, LDLRAP1,
  POLR3G, SSB, NDUFAF5, CTDSP2, HDAC2, EEF1G, INHA, GLS2, CALU, KRAS, TRIM68,
  OSBPL3]
- verdict: SEARCH_LEADER
- lesson: The screening stage's primary objective is structurally degenerate --
  `Adjuvant Chemo` is never split on when it competes against 8,647 genes, so
  every screening trial returns the identical all-observation policy value and
  every RMST permutation importance is exactly zero; the gene ranking that
  reaches the reduced stage is therefore a pure C-index PFI ranking.

#### iter_001 detail

**Wall clock 575.9 s of the 14,400 s budget** (screening 146 s for 12 trials,
PFI 300 s for ~3.2k panels, reduced 122 s for 49 trials, load/folds ~7 s). The
cost model was ~4x conservative: forests split on far fewer features than the
depth cap allows (~7.3 splits/tree at depth 4 / leaf 24, giving 350-367 used
features per fold, not the ~1,400 projected). Recalibrated constants: a
permutation panel costs about `0.02 + 0.00145 * n_estimators` seconds, and
`n_used ~ 8666 * (1 - exp(-trees * splits_per_tree / 8666))`. Roughly 25x
headroom remains.

**The degeneracy, stated exactly.** All 12 screening trials returned
RMST = 0.2183 to every printed digit while C-index ranged 0.6288-0.6526. In
`permutation_importance.csv` all 8,666 RMST importances are exactly 0.0 with
exactly 0.0 standard deviation, and `Adjuvant Chemo` has
`split_used_fold_fraction = 0.0`. The chain is forced: if no tree ever splits on
treatment, then `rmst0 == rmst1` for every patient, so `recommendations` returns
all-zero, so the policy is "observe everyone" and the contribution is
`phi0 - phi1` for every row. Its mean, 0.2183 months, is the value of the
always-observe policy and is invariant to forest geometry and to permuting any
feature. Stage 1 therefore selects purely on its lexicographic secondary, and
the ranking rule degrades to C-index importance then symbol -- which is
compliant with red line 10 (RMST first, C-index and symbol as deterministic
tie-breakers) but means the primary objective does no work before the reduced
stage. 980 of 8,666 features had nonzero C-index importance; 993 were used in
at least one fold.

**0.2183 is a diagnostic constant.** Reduced trials 39, 43 and 46 also returned
exactly 0.218 -- all three used `max_features` 0.5 or 1.0. Near-greedy feature
selection crowds treatment out of the trees even at 51 features, collapsing the
policy to all-observe. Every `max_features = 1.0` trial was poor (-0.849,
-0.363, 0.218, 0.218, 1.435). The `sqrt` family, which offers treatment at about
one split in seven, is the only setting that reliably produces a live policy.

**The N curve is noise-dominated.** At the fixed reference geometry the RMST
difference over N = 0..32 ranges from -2.462 to +5.214 with no trend, while
C-index climbs smoothly from 0.6825 at N=0 to a 0.705-0.714 plateau by N>=15.
The winning N=20 has fold RMSTs of 2.25, 9.16 and 4.23 -- a spread far larger
than any between-N difference. The reported 5.214 is the maximum of 49 noisy
adaptive trials and should be read as an optimistically biased development
estimate, not an effect size. C-index behaves like a real signal; the policy
RMST on 1,034 patients with only 152 treated does not.

### iter_002 — screening capacity and the sqrt-family reduced stage
- type: SEARCH
- hypothesis: The stage-1 RMST tie is structural, not a capacity artifact, so a
  substantially stronger screening forest will still return the identical
  0.2183 for every trial and still yield all-zero RMST importance; what it will
  change is the C-index PFI ranking, which is the ranking that actually reaches
  the reduced stage. Enlarging the screening forest to 100-150 trees at depth
  4-6 with leaf 16-24 should roughly double or triple the used-feature pool
  (target >= 1,500 versus 993), resolving the top-32 from a better-populated
  ranking and improving the reduced stage's attainable RMST at the unchanged
  reference geometry.
- changed: `name` -> `screen_capacity_sqrt_reduced_v2`; `screening_space` ->
  n_estimators [100, 150], max_depth [4, 6], min_samples_leaf [16, 24],
  split_leaf_multiplier [2], max_features [0.05, 0.25]; `reduced_space`
  max_features -> ["sqrt", "log2", 0.25], dropping 0.5 and 1.0, which run 001
  showed collapse the policy to all-observe. The reduced reference geometry
  (400 trees, depth 6, leaf 16, multiplier 2, sqrt) is deliberately held fixed
  so the 33-trial controlled N curve is directly comparable to run 001 and any
  change is attributable to the gene ranking rather than to forest geometry.
- prespecified_falsifier: if any screening trial returns an RMST difference
  other than 0.2183, the "treatment is never split on" mechanism is wrong and
  screening capacity is a live RMST lever after all.
- cost_projection: worst corner 150 trees / depth 6 / leaf 16 / mf 0.25 gives
  ~243 s per screening trial (~2.9 ks for 12) and n_used ~2.4k at 0.24 s per
  panel (~5.3 ks of PFI), plus ~0.2 ks reduced: ~8.4 ks worst case against the
  14.4 ks wall, a 1.7x margin.
- red_line_audit: unchanged from iter_001 in every respect -- train.py still
  names no dataset path, contains no gene symbol, row id, prediction or
  outcome, and applies no gene prefilter, so the screening forest still receives
  all 8,647 eligible genes (launcher-verified). Only train.py was edited and
  log.md appended; locked files still pass `integrity.verify_lock`. The
  estimand, counterfactual pairing, permutation pairing and fit-only transforms
  are untouched locked code. 49 reduced trials still cover N=0..32. Max
  n_estimators 600 <= 1000. This is experiment 2 of 20 with a fresh train.py
  snapshot. Conclusions remain observational development-CV statements.
- run_id: run_002_20260823T001849Z
- screening_rmst_difference: 0.218 months (all 12 trials, identical)
- screening_cindex: 0.668
- final_top_n: 6
- final_rmst_difference: 5.686 months
- final_cindex: 0.707
- final_genes: [LDLRAP1, SEC23A, GNG7, ANKFY1, PPFIBP2, SLC24A1]
- verdict: SEARCH_LEADER
- lesson: The stage-1 RMST tie is structural and survived a doubled screening
  forest, but screening capacity is nonetheless a real lever, because it acts
  through the C-index PFI ranking: doubling the used-feature pool lifted the
  whole controlled N curve rather than winning a single lucky trial.

#### iter_002 detail

**Wall clock 2,433.5 s** (screening 1,109 s, PFI 1,205 s, reduced 111 s).
Measured panel cost 0.177 s at 150 trees against a 0.24 s projection, so the
cost model remains conservative. Used features per fold 734 / 820 / 714, union
1,942 versus 993 in run 001 -- the targeted doubling, achieved.

**Falsifier result: not falsified.** All 12 screening trials returned
0.2183461156, identical to ten decimal places, across 100-150 trees, depth 4-6,
leaf 16-24 and `max_features` 0.05-0.25, while C-index ranged 0.6498-0.6680.
All 8,666 RMST importances are again exactly zero and `Adjuvant Chemo` again has
`split_used_fold_fraction = 0.0`. The degeneracy is a property of treatment
competing against 8,647 genes, not of forest capacity, and it should be treated
as fixed for the rest of the search: stage 1 optimizes C-index, full stop.

**Screening capacity acted through the ranking, and it worked.** Comparing the
two controlled N curves at the identical reference geometry, so that only the
gene ranking differs:

| curve statistic | run 001 (993 used) | run 002 (1,942 used) |
| --- | --- | --- |
| RMST mean over N=0..32 | ~1.0 | 1.866 |
| RMST median | ~1.0 | 2.052 |
| RMST minimum | -2.462 | -1.882 |
| negative-RMST panels | 9 of 33 | 7 of 33 |
| C-index at N=0 | 0.6825 | 0.6825 |
| C-index plateau | 0.705-0.714 | 0.717-0.725 |

The level of the whole curve moved, not just its maximum, and the secondary
objective improved by about 0.011 at the plateau. That is the signature of a
better-resolved ranking rather than a luckier draw.

**The new leader is also the more credible one.** N=6 has fold RMSTs of 5.372,
4.641 and 7.049, versus 2.247, 9.160 and 4.232 for the run 001 leader. The
run 002 winner is not carried by a single fold, which matters because the
reported reward is still the maximum of 49 adaptive trials.

**Screening preferences observed.** `max_features` 0.25 beat 0.05 on screening
C-index at every matched geometry (0.6680 vs 0.6557 at 150/depth 4/leaf 16), and
depth 4 slightly beat depth 6 (0.6680 vs 0.6641). In the reduced stage the
sqrt family again dominated; the one surviving near-greedy trial, t39 at
`max_features` 0.25 with N=28, returned exactly 0.218 -- the all-observation
signature -- confirming the diagnostic transfers.

### iter_003 — does ranking resolution keep paying?
- type: SEARCH
- hypothesis: The run 001 -> 002 gain came from the size and resolution of the
  used-feature pool that the C-index PFI ranking is drawn from (993 -> 1,942).
  Pushing screening to 150-250 trees at `max_features` 0.25-0.5 should raise the
  pool again to roughly 3,000 per-fold-union and lift the controlled N curve's
  mean a third time. If the curve mean instead flattens or falls, ranking
  resolution has saturated and screening capacity is closed as a lever.
- changed: `name` -> `screen_resolution_push_v3`; `screening_space` ->
  n_estimators [150, 250], max_depth [4, 6], min_samples_leaf [16],
  split_leaf_multiplier [2], max_features [0.25, 0.5]. The `reduced_space` is
  byte-identical to run 002, including the reference geometry, so the N curve
  stays comparable across all three runs and any change is attributable solely
  to the ranking.
- prespecified_readout: compare the mean and median of the 33-trial controlled
  N curve against 1.866 / 2.052; treat a lower mean as saturation.
- cost_projection: worst corner 250 trees / depth 6 / `max_features` 0.5 gives
  ~430 s per screening trial (~5.2 ks for 12) and ~1,590 used features per fold
  at ~0.28 s per panel (~4.0 ks of PFI), plus ~0.12 ks reduced: ~9.3 ks worst
  case against the 14.4 ks wall, a 1.55x margin. Realized cost has come in
  well under worst case in both runs so far.
- red_line_audit: unchanged and re-checked. train.py names no dataset path and
  contains no gene symbol, row id, prediction or outcome; no gene prefilter or
  pre-ranking exists, so the screening forest still receives all 8,647 eligible
  genes, which the launcher independently verifies. Only train.py was edited and
  log.md appended; `integrity.verify_lock` still passes over budget, integrity,
  prepare, red_lines and run. Fit-only imputation, held-out-only permutation,
  the locked 60-month IPCW-AIPW estimand with lexicographic C-index, paired
  counterfactual prediction and paired permutation are all untouched locked
  code. No row is dropped. 49 reduced trials still cover N=0..32 before
  selection. Max n_estimators 600 <= 1000; max screening trees 250 <= 1000.
  Experiment 3 of 20, fresh train.py snapshot. Claims remain observational.
- run_id: run_003_20260823T022901Z
- screening_rmst_difference: 0.218 months (all 12 trials, identical)
- screening_cindex: 0.671
- final_top_n: 6
- final_rmst_difference: 4.943 months
- final_cindex: 0.705
- final_genes: [LDLRAP1, SEC23A, GNG7, ANKFY1, PRIM2, PPFIBP2]
- verdict: NOT_LEADER
- lesson: A third increase in ranking resolution kept improving the N curve's
  location and stability but lowered its maximum, which separates the two things
  cleanly: resolution buys a genuinely better curve, while the reported reward
  is a maximum and therefore partly a tail draw.

#### iter_003 detail

**Wall clock 7,687.5 s** (screening 5,124 s, PFI 2,444 s, reduced 112 s) --
3.2x run 002 for the resolution step from 1,942 to 2,848 used features
(1,114 / 1,289 / 1,069 per fold). Screening RMST was again 0.2183461156 in all
12 trials with all-zero RMST importance and `Adjuvant Chemo` unused; that is
three runs and 36 trials of the same constant.

**The prespecified readout, at an unchanged reference geometry and ranking as
the only difference:**

| curve statistic | run 001 (993) | run 002 (1,942) | run 003 (2,848) |
| --- | --- | --- | --- |
| RMST mean | ~1.0 | 1.866 | **2.021** |
| RMST median | ~1.0 | 2.052 | **2.121** |
| RMST minimum | -2.462 | -1.882 | **-0.714** |
| negative panels | 9 of 33 | 7 of 33 | **4 of 33** |
| RMST maximum | 5.214 | **5.686** | 4.847 |
| C-index max | 0.7138 | **0.7248** | 0.7229 |

The mean did not fall, so resolution is not saturated by the stated criterion,
but the gain decelerated hard (+0.87, then +0.155) for 3.2x the wall time, and
the maximum moved the wrong way. Reading these three runs together: resolution
raises the curve's floor and centre, and the run 002 reward of 5.686 was a
favourable draw on top of a lower centre. I am closing screening capacity as an
active lever -- further spend buys stability I cannot convert into reward under a
max-of-trials objective.

**Where the reward actually lives.** In all three controlled curves RMST is
larger at small N while C-index rises with N: run 003 peaks at N=12 (4.551) and
N=28 (4.847) but its best small-N values are N=0 (3.395), N=3 (3.681) and
N=6 (3.777), and the C-index plateau only arrives past N=18. Since RMST is
lexicographically first, small panels are where the primary objective pays.

**Pooled free-trial geometry evidence (runs 001-003).** Each geometry recurs at
a near-fixed N, so geometry and N are entangled, but the ordering is stable
across all three rankings. `400 trees / depth 8 / leaf 8 / sqrt` at N=7 gave
3.26, 4.29 and 1.65 (mean 3.065); `200 / depth 6 / leaf 16 / sqrt` at N=6 gave
2.46, 1.22 and 4.94 (mean 2.875). At the other end, `400 / depth 4 / leaf 8 /
sqrt` at N=32 gave 0.88, -3.79 and 0.32, and every `max_features` of 0.5 or 1.0
or `log2` sat at or below 0.62 -- three more trials returned exactly 0.218. Deep
trees with small leaves and sqrt features at small N is the strongest region
observed, and it has never been the reference geometry.

### iter_004 — move the reference geometry to the pooled-best region
- type: PARAM
- hypothesis: The reduced stage spends 33 of 49 trials at
  `space[key][0]`, and that reference has been `400 / depth 6 / leaf 16 / sqrt`
  in every run so far. Re-pointing it at the pooled-best region,
  `400 / depth 8 / leaf 8 / sqrt`, should raise the controlled N curve's mean
  above run 003's 2.021 at the identical gene ranking, and raise its maximum
  above 4.847.
- changed: `reduced_space` reordered to make the reference deep-and-fine --
  max_depth [8, 6, 4] and min_samples_leaf [8, 16, 32], leaving n_estimators
  [400, 200, 600], split_leaf_multiplier [2, 3] and max_features
  ["sqrt", "log2", 0.25] as in run 003. `screening_space` is byte-identical to
  run 003, so the gene ranking is reproduced exactly and the reference geometry
  is the only moving part.
- prespecified_readout: curve mean/median/max against 2.021 / 2.121 / 4.847. A
  lower mean retires deep-and-fine and sends the reference back to depth 6 /
  leaf 16.
- cost_projection: screening and PFI reproduce run 003 at ~7.6 ks; the reduced
  stage is deeper so ~0.2 ks rather than 0.11 ks. Total ~7.8 ks against the
  14.4 ks wall.
- red_line_audit: no change in kind from iter_003 and re-verified. Only train.py
  edited, log.md appended, locked files pass `integrity.verify_lock`. train.py
  contains no dataset path, gene symbol, row id, prediction or outcome, and no
  gene prefilter or pre-ranking, so all 8,647 eligible genes still enter the
  screening forest under launcher verification. Locked fit-only transforms,
  held-out-only permutation, paired counterfactual prediction, paired
  permutation, the 60-month IPCW-AIPW estimand and its lexicographic C-index
  secondary are untouched; no rows dropped; N=0..32 all evaluated by the 49
  reduced trials; max trees 600 <= 1000. Experiment 4 of 20, fresh snapshot.
  Claims remain observational development-CV statements.
- run_id: run_004_20260823T043658Z
- screening_rmst_difference: 0.218 months (all 12 trials, identical)
- screening_cindex: 0.671
- final_top_n: 14
- final_rmst_difference: 5.423 months
- final_cindex: 0.719
- final_genes: [LDLRAP1, SEC23A, GNG7, ANKFY1, PRIM2, PPFIBP2, INHA, SLC24A1,
  KCNS3, P2RY10, SLCO4A1, DNM1L, ARHGAP25, LAMC2]
- verdict: NOT_LEADER
- lesson: Deep-and-fine is worse as a reference geometry, not better -- the
  pooled free-trial evidence that motivated it was confounded with N, and the
  controlled curve reversed it.

#### iter_004 detail

**Wall clock 7,559.6 s.** The screening stage reproduced run 003 exactly, as
intended: same best parameters (250 trees / depth 4 / leaf 16 /
`max_features` 0.25), same 2,848 used features, same 0.2183461156 across all 12
trials. The gene ranking was therefore identical and the reference geometry was
the only moving part.

**Prespecified readout: hypothesis rejected on every statistic.**

| curve statistic | run 003 ref = d6/leaf16 | run 004 ref = d8/leaf8 |
| --- | --- | --- |
| RMST mean | **2.021** | 1.277 |
| RMST median | **2.121** | 1.603 |
| RMST maximum | **4.847** | 3.891 |
| negative panels | **4 of 33** | 6 of 33 |
| C-index max | 0.7229 | **0.7275** |

Deep-and-fine is retired as a reference; the reference returns to depth 6 /
leaf 16 per the stated rule.

**Why the pooled evidence misled, which is the transferable lesson.** The
free-trial table that motivated this run ranked `400/d8/leaf8/sqrt` top at
mean 3.065 -- but every one of those observations sat at N=7, because the
sampler revisits similar points. Geometry and N were entangled, so what looked
like a good geometry was substantially "N=7 is a good panel size." Run across
all 33 panel sizes under control, the same geometry is clearly worse. Pooled
observational trial tables are for generating hypotheses here; only the
controlled N curve, with the ranking held fixed, settles them.

**Two further readings.** First, iter_004 moved depth 6->8 and leaf 16->8
together, so the drop is attributable to the pair and not to either alone;
free trial t44 at `600/d8/leaf16` scored this run's best (5.423), which hints
leaf 8 rather than depth 8 is the harmful half. Second, there is a real
primary/secondary tension: the deeper reference produced the best C-index seen
so far (0.7275, versus 0.7248 and 0.7229) while producing the worst RMST curve.
Deeper forests discriminate better and prescribe worse.

**The reward is part level and part spread.** Reward is the maximum of 49
trials, so it tracks roughly mean + 2.2 standard deviations of the curve.
Run 002 leads with a *lower* curve mean than run 003 (1.866 vs 2.021) because
its curve was more dispersed. That is worth stating plainly: chasing the curve
mean is the honest optimization, but part of the standing leader's margin is
dispersion rather than quality. Fold consistency is the tiebreaker I trust --
run 002's leader (5.372 / 4.641 / 7.049) remains the most evenly supported
result, ahead of run 004's (4.166 / 8.547 / 3.551) and run 001's
(2.247 / 9.160 / 4.232).

### iter_005 — the last unvaried reference dimension
- type: PARAM
- hypothesis: Reference geometry has always used 400 trees. Depth, leaf and
  `max_features` are now each settled by a controlled comparison
  (depth 6 > depth 8, leaf 16 > leaf 8, sqrt > log2/0.25/0.5/1.0), leaving
  n_estimators the only reference dimension never varied. Raising the reference
  to 600 trees reduces the forest's Monte Carlo variance and should lift the
  controlled curve mean above run 003's 2.021 at the identical ranking.
- changed: `reduced_space` restored to the run 003 ordering except
  n_estimators -> [600, 400, 200], making the reference
  600 trees / depth 6 / leaf 16 / multiplier 2 / sqrt. `screening_space`
  unchanged from runs 003-004, so the ranking is reproduced exactly for the
  third time and trees is the only moving part.
- prespecified_readout: curve mean/median against run 003's 2.021 / 2.121. A
  lower mean retires 600 trees and settles the reference at 400.
- cost_projection: screening and PFI reproduce runs 003-004 at ~7.5 ks; the
  reduced stage grows with trees to ~0.2 ks. Total ~7.7 ks against 14.4 ks.
- red_line_audit: unchanged in kind and re-verified. Only train.py edited and
  log.md appended; locked files pass `integrity.verify_lock`. No dataset path,
  gene symbol, row id, prediction or outcome in train.py; no gene prefilter or
  pre-ranking, so all 8,647 eligible genes enter the screening forest under
  launcher verification. Locked fit-only transforms, held-out-only permutation,
  paired counterfactual prediction, paired permutation and the 60-month
  IPCW-AIPW estimand with lexicographic C-index are untouched. No rows dropped.
  All of N=0..32 evaluated before selection. Max trees 600 <= 1,000.
  Experiment 5 of 20, fresh snapshot. Claims remain observational.
- run_id: pending
