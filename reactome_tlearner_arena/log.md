# Reactome RSF T-learner search log

No full experiments have been run. The shipped candidate is an infrastructure
baseline, not a scientific result.


---

## Baseline review (smoke, infrastructure only — no ledger slot)

Re-ran `run --smoke` (19.9 s wall). It reproduces the shipped
`runs/smoke_20260823T181235Z.json` exactly. Development pool: n=1034,
465 events, 152 ACT (14.7%). Reactome intersect: 8,647 usable genes across
1,839 pathways.

Smoke verdict on the shipped candidate (`dr_pathway`, 8 genes, 2 modules,
max_features 0.5): **not eligible**, increment **-1.61 months** — the genomic
panel actively *hurt*. Failing gates: `genomic_increment_positive`,
`genomic_value_at_least_clinical`, `gene_selection_jaccard_at_least_0_10`.

Four things about the smoke configuration limit what it can tell us, and they
shape the prespecification below:

1. **The seed gates are untested.** Smoke runs one seed.
   `_mean_pairwise_correlation` over a single row has no pairs and returns
   1.0 by construction, and `seed_agreement` is likewise trivially 1.0. The
   observed 1.0/1.0 is an artifact, not evidence. With three seeds,
   `seed_agreement >= 0.85` and `seed_benefit_correlation >= 0.50` are live
   risks — especially in the ACT arm, which carries only ~114 patients /
   ~70 events per 4-fold training partition.
2. **The Jaccard failure is largely an artifact.** Smoke uses two *disjoint*
   half-samples; the full run uses 4-fold partitions whose fitting sets
   overlap in 2/3 of patients. Selection stability should be far higher at
   full fidelity, so this gate is not the primary design target.
3. **Genomic capacity is the visible defect.** The two modules absorbed
   38-43% of all splits while displacing clinical structure: C-index fell
   0.652 -> 0.640 and alignment fell 6.56 -> 4.95. Eighteen of the nineteen
   pretreatment columns are binary; the modules are standardized continuous
   averages, so RSF split-selection bias hands them far more splits than
   their signal warrants. The gate only asks for a mean genomic split
   fraction >= 0.01 — three orders of magnitude below what we are spending.
4. **Module averaging cancels the signal it is built from.**
   `_gene_effect_scores` ranks genes by the *absolute* partial correlation
   with the DR benefit pseudo-outcome, and `FeatureTransformer` groups the
   returned list by rank position via `array_split`. A module is therefore a
   mean of standardized genes with *mixed* benefit directions, and
   positively- and negatively-associated genes cancel. The shipped module
   representation is testing a partly self-defeating feature.

Wall-clock is not a binding constraint. Synthetic-shape benchmarks (no real
data touched) put the 96 forests at roughly 9 minutes against the 35-minute
limit; the locked clinical forests are cheap because their features are
almost entirely binary.

## Experiment 1 (prespecified before running)

**Candidate** `tlearner_raw4_drgene_locked_geometry`

- selector `dr_gene`, `n_genes` 4, representation `raw`, `module_count` 4
- `benefit_threshold_months` 0.0
- observation arm: 1000 trees, depth 9, leaf 8, split 16, max_features 1.0
- ACT arm: 1000 trees, depth 6, leaf 16, split 32, max_features 1.0

**Hypothesis.** The shipped candidate fails because it spends far too much
model capacity on genomics and spends it on a representation that cancels
its own signal. A T-learner whose observation arm is *geometrically identical
to the locked clinical learner* and which adds only the four most
reproducibly benefit-correlated Reactome genes as raw, signed features will
produce a positive incremental alignment in both repeats.

**Reasoning per change, versus the shipped candidate.**

- `module` -> `raw`: removes the sign-cancellation defect entirely. Each gene
  is its own signed feature and the forest resolves direction natively. This
  costs more genomic columns (4 vs 2) but each column carries real signal.
- `n_genes` 8 -> 4 (the arena minimum): the program names smaller panels as
  the preferred first axis, and it is the most direct lever on the observed
  capacity overspend.
- `dr_pathway` -> `dr_gene`: `dr_pathway` ranks 1,839 pathways by the mean of
  their top-3 gene scores and then takes one gene per pathway. Ranking that
  many pathways on a noisy score is a selection-on-noise procedure, and it
  deliberately spreads the panel across *unrelated* pathways. `dr_gene` takes
  the global top-4, which is the most reproducible signal available and the
  better bet for the Jaccard gate.
- `max_features` 0.5 -> 1.0 on the observation arm: the locked clinical
  comparator runs at 1.0. Running the C+G arm at 0.5 handicaps it on the
  *clinical* dimension, so part of the shipped candidate's deficit may be a
  self-inflicted comparator mismatch rather than a genomic failure. Matching
  the geometry makes the added genes the only difference in the observation
  arm and gives a clean read.
- ACT arm depth 7->6, leaf 12->16, split 24->32: the program names stronger
  ACT-arm regularization as a preferred first axis. With ~114 patients per
  fitting partition, this arm is the noisiest term in
  `benefit = RMST1 - RMST0`, and it is the most likely cause of a
  seed-correlation failure. This is the one deliberate deviation from locked
  clinical geometry, and it is a variance-reduction move.

**Prespecified predictions and what would falsify them.**

- Mean genomic split fraction lands well below the shipped 0.38-0.43 (target
  roughly 0.10-0.25) while staying clear of the 0.01 floor.
- Incremental alignment is positive in both repeats.
- `seed_benefit_correlation` >= 0.50 with three seeds. If it comes in below
  that, ACT-arm variance is the dominant failure mode and the next axis is
  much stronger ACT-arm regularization, not panel composition.
- If the increment is still negative with capacity this low and the
  representation sign-clean, the failure is not a capacity or representation
  artifact and the next axis is the selector's *target*, not its size.

**Result.** (to be appended after the run)

**Result — run_001_20260823T190323Z, 413 s wall. NOT ELIGIBLE, reward -1e6.**

| field | repeat 1 | repeat 2 |
|---|---|---|
| incremental alignment (months) | **-1.660** | **-4.571** |
| selection LCB | -8.519 | -13.116 |
| bootstrap mean / CI95 | -1.616 / [-6.001, +2.327] | -4.565 / [-10.229, +0.500] |
| alignment C / C+G | 5.098 / 3.438 | 5.361 / 0.790 |
| value C / C+G / best constant | 48.152 / 46.515 / 45.534 | 48.526 / 45.080 / 45.493 |
| Harrell C, C / C+G (drop) | 0.6561 / 0.6508 (0.0053) | 0.6641 / 0.6395 (0.0246) |
| ACT recommended fraction (C+G) | 0.3066 | 0.2921 |
| seed agreement / benefit corr | 0.9907 / 0.9986 | 0.9907 / 0.9985 |
| predicted benefit mean / IQR / nontrivial | -4.588 / 13.757 / 0.9894 | -5.202 / 14.663 / 0.9903 |
| genomic split fraction obs / ACT | **0.6346 / 0.4703** | **0.6373 / 0.5564** |
| arm support obs (pt/ev) | 661.5 / 278.25 | 661.5 / 278.25 |
| arm support ACT (pt/ev) | 114.0 / 70.5 | 114.0 / 70.5 |
| raw propensity overlap / IPTW ESS | 0.8723 / 358.5 | 0.8694 / 366.3 |

Diagnostic score -16.027 (robust LCB -13.116, repeat range 2.911).
Gene Jaccard **0.0697**. Full-development selection:
`UQCRH, EPHA1, SPTB, MYH8`. Per-fold selections were near-disjoint; only
`REG3A` recurred more than twice (4/8 folds).

Gates: 9/13. Failing — `genomic_increment_positive`,
`genomic_value_at_least_clinical`, `genomic_value_at_least_best_constant`,
`gene_selection_jaccard_at_least_0_10`.

**Two prespecified predictions were falsified, and both are informative.**

1. *"Genomic split fraction lands around 0.10-0.25."* Wrong, and in the wrong
   direction: it rose to **0.63/0.56**, above the shipped candidate's 0.40.
   Four raw continuous genes against eighteen binary clinical dummies take
   two thirds of all splits. A companion simulation of split competition
   (synthetic, no real data) reproduces this and shows why the lever I
   reached for is the wrong one: across its entire legal range
   `max_features` moves the genomic share only 0.41 -> 0.51 at four genomic
   columns, whereas the *column count* moves it from ~0.15 (one column) to
   ~0.67 (eight). **Genomic column count, not `max_features`, is the
   capacity lever.** Matching the locked clinical `max_features` was
   therefore a cost with no compensating benefit.
2. *"The Jaccard failure was a smoke artifact of disjoint halves."* Wrong.
   At full fidelity, with 4-fold fitting sets that share two thirds of their
   patients, Jaccard was **0.0697** — still under the 0.10 gate. The top-4
   of an 8,647-gene ranking by DR-benefit partial correlation is essentially
   not reproducible. This is the single most important finding of the run:
   `n_genes = 4` is not "the preferred small panel", it is below the
   resolution of the selection signal.

**One prediction was confirmed, and it retires a constraint.** Seed agreement
0.991 and seed benefit correlation **0.9985** — far above the 0.85/0.50
gates. Averaging 1,000 trees makes the three seeds nearly identical. ACT-arm
seed instability is *not* a live failure mode, so ACT-arm regularization can
be spent on bias/variance rather than on reproducibility, and the extra
regularization I applied to the ACT arm bought nothing.

**Revised diagnosis.** The failure is not capacity *or* stability; it is that
the model is spending 63% of its splits on a signal whose top-4 does not
reproduce across folds. Both must be fixed together, and the two fixes are
the same fix: put the genomic signal into a small number of *averaged* columns
built from a *larger* gene set.

That implies module representation — but the modules must be sign-coherent,
because `_gene_effect_scores` ranks on |correlation| and `array_split` groups
by rank, so an averaged module of mixed-direction genes cancels its own
signal in expectation. Locked machinery cannot produce a sign-coherent
module, so experiment 2 changes selector logic.

## Experiment 2 (prespecified before running)

**Candidate** `tlearner_pathway_benefit_module16`

- custom pathway-anchored sign-coherent selector (below), `n_genes` 16
- representation `module`, `module_count` **1** — a single genomic column
- `benefit_threshold_months` 0.0
- **both** arms exactly the locked clinical geometry (obs 1000/9/8/16/1.0,
  ACT 1000/7/12/24/1.0), so the module column is the *only* difference
  between the clinical and C+G panels in either arm

**Selector logic** (injected through the locked evaluator's own `selector`
hook; the recorded enum stays `dr_gene`, the closest of the three allowed
labels, and the effective procedure is this one):

1. Score every available gene by the **signed** partial correlation with the
   cross-fitted DR benefit pseudo-outcome, replicating the locked
   `_gene_effect_scores` residualization exactly but keeping the sign that
   `np.abs` discards.
2. Score every Reactome pathway with **>= 25** development-present members by
   the *mean signed* score of its members.
3. Pool the members of the **top 3** pathways.
4. Return the 16 pool genes with the most positive scores.

Fit-only throughout: the selector sees only the fitting partition, Reactome
membership, the available gene list, and the validated spec. No gene is
hard-coded.

**Hypothesis.** Experiment 1's deficit is caused by spending 63% of splits on
a signal whose top-4 does not reproduce (Jaccard 0.070). Compressing a larger
gene set into one sign-coherent averaged column, drawn from a pathway-
restricted pool, will cut genomic split share to roughly 0.20, lift Jaccard
over the 0.10 gate, and turn the increment positive in both repeats.

**Why each element, tied to what experiment 1 actually showed.**

- **One module column.** The split-competition simulation and experiment 1
  agree that column count is the capacity lever: 4 raw columns -> 0.63 share.
  One column is the arena minimum and should land near 0.20. `max_features`
  is not used as a lever because it moves the share by only ~0.1 across its
  whole range, and experiment 1's repeat-2 C-index drop (0.0246 against a
  0.03 gate) leaves no room to weaken the clinical side of the C+G forest.
- **16 genes, not 4.** Experiment 1 established that 4 is below the
  resolution of the selection signal. Sixteen is the arena maximum, gives the
  module 16-fold averaging, and mechanically raises the Jaccard denominator's
  tolerance for rank noise.
- **A pathway-restricted pool.** This is the part aimed squarely at the
  Jaccard gate. Ranking 8,647 genes individually is selection on noise; a
  mean over >= 25 pathway members is a far lower-variance statistic, and
  restricting the final pick to ~100-300 pooled members raises the
  probability that two folds choose overlapping genes. It is also what the
  program actually asks for — a *small stable Reactome panel* — rather than
  16 unrelated genes.
- **Sign coherence.** Without it a one-module design is self-defeating: the
  locked ranking is on |correlation| and the locked transformer averages by
  rank position, so a mixed-direction module cancels in expectation. Taking
  the most positive scores makes the module a genuine benefit axis.
- **Locked geometry in both arms.** Experiment 1 showed seed stability is a
  non-issue (0.9985), so the extra ACT-arm regularization bought nothing and
  only muddied attribution. Reverting it makes the module column the single
  difference between the panels.

**Prespecified predictions and what would falsify them.**

- Genomic split fraction lands near 0.15-0.25 in both arms. If it is again
  above 0.4 with a single column, then split share is not controllable at all
  within this feature geometry and the whole module/raw axis is closed.
- Jaccard >= 0.10. If pathway pooling still cannot reach 0.10, the DR benefit
  signal has no reproducible gene-level content at this sample size, and the
  remaining path is to stop trying to stabilize *which* genes are chosen and
  instead make the panel's *effect* small enough to be harmless while still
  clearing the 0.01 split-fraction floor.
- Increment positive in both repeats. If the increment is still negative even
  with low capacity, a sign-coherent module and a stable panel, then the DR
  benefit signal does not carry incremental policy value over the clinical
  T-learner, and the honest next axes are the decision rule
  (`benefit_threshold_months`) rather than the panel.

**Result.** (to be appended after the run)

**Result — run_002_20260823T191143Z, 226 s wall. NOT ELIGIBLE, reward -1e6.**

| field | repeat 1 | repeat 2 |
|---|---|---|
| incremental alignment (months) | **+1.470** | **-0.153** |
| selection LCB | -4.070 | -6.142 |
| bootstrap mean / CI95 | +1.404 / [-2.282, +4.791] | -0.106 / [-4.027, +3.574] |
| alignment C / C+G | 5.098 / 6.568 | 5.361 / 5.209 |
| value C / C+G / best constant | 48.152 / 48.269 / 45.570 | 48.526 / 47.943 / 45.442 |
| anti-policy value C / C+G | 43.054 / 41.701 | 43.165 / 42.734 |
| Harrell C, C / C+G (drop) | 0.6561 / 0.6603 (**-0.0041**) | 0.6641 / 0.6528 (0.0113) |
| ACT recommended fraction C / C+G | 0.3723 / 0.3317 | 0.3820 / 0.3685 |
| seed agreement / benefit corr | 0.9929 / 0.9991 | 0.9929 / 0.9992 |
| predicted benefit mean / IQR / nontrivial | -5.160 / 16.422 / 0.9961 | -4.014 / 16.687 / 0.9942 |
| genomic split fraction obs / ACT | 0.3791 / 0.2826 | 0.3881 / 0.2884 |
| arm support obs (pt/ev) | 661.5 / 278.25 | 661.5 / 278.25 |
| arm support ACT (pt/ev) | 114.0 / 70.5 | 114.0 / 70.5 |
| raw propensity overlap / IPTW ESS | 0.8723 / 358.5 | 0.8694 / 366.3 |

Diagnostic score **-7.764** (robust LCB -6.142, repeat range 1.622), up from
-16.027. Gene Jaccard **0.1462** — gate cleared. Full-development module (one
module): `KCNC3, OR2S2, KCNV1, TAAR5, ANO2, OR10H2, HTR6, OR2J2, OR3A2,
OR52A1, OR10C1, OR2F1, OR1D2, KCNA6, OR2W1, KCNQ3`.

Gates: **11/13**. Failing — `genomic_increment_positive`,
`genomic_value_at_least_clinical`. Both failures are in repeat 2 only, and
both are small: increment -0.153 and a value shortfall of 0.583.

**Scorecard against the prespecification.**

- *Jaccard >= 0.10*: **confirmed**, 0.070 -> 0.1462. Pathway pooling is the
  right instrument for selection stability.
- *Split fraction 0.15-0.25*: **partly falsified** — 0.379/0.288. One module
  column roughly halved it from experiment 1's 0.63, confirming column count
  as the lever, but a single standardized continuous column still outdraws
  eighteen binary dummies by far more than its 1-in-19 share. My simulation
  understated the continuous-feature advantage by about 1.9x.
- *Increment positive in both repeats*: **half confirmed**. Repeat 1 reached
  +1.470 with the C-index actually *improving* (-0.0041 drop), which is the
  first direct evidence in this arena that the module carries real signal
  rather than just noise. Repeat 2 came in at -0.153.

**What the two remaining failures actually are.** Not capacity, not stability,
not seeds. The module's contribution is *real but not yet reliable*: it is
worth +1.47 months of alignment in one patient partition and -0.15 in the
other. Shrinking its influence further would shrink the gain along with the
loss and park both increments near zero, which does not pass a strictly-
positive gate either. The contribution has to become more consistent, not
smaller.

**A caveat that must survive into any write-up.** Six of the eight folds
anchored on olfactory-receptor and potassium-channel families (`OR*`, `KCN*`,
`TAAR5`, `HTR6`); the other two jumped to histone/HLA sets. Olfactory
receptors are not plausibly expressed in lung tumour tissue, so this module
is more likely tracking a low-expression/array-background axis that happens
to correlate with the DR benefit pseudo-outcome than any olfactory biology.
That is consistent with the arena's own limitation list, and it does not
affect gate arithmetic, but the panel should not be described as a
biologically interpretable pathway result.

**Two candidate next axes, and why I am taking the first.**

1. *Stabilize the pathway anchor.* Two of eight folds jumped to an unrelated
   pathway family, which is fold-level heterogeneity injected straight into
   the repeat range (1.622). The current pool is the full membership of the
   top 3 pathways, so a single pathway flip replaces the entire pool.
2. *Filter to genes with real expression dynamic range* before scoring, which
   would push selection off the near-null `OR*` families.

Axis 1 first: it targets the measured defect (fold-level pathway switching)
and it improves the reward term directly through the repeat range, whereas
axis 2 replaces a panel that is already clearing the stability gate and could
lose it. Axis 2 is the fallback if axis 1 does not close repeat 2.

## Experiment 3 (prespecified before running)

**Candidate** `tlearner_bounded_pool_module16`

Identical to experiment 2 in every respect — 16 genes, one sign-coherent
module, both arms at locked clinical geometry, threshold 0.0 — except the
construction of the candidate pool:

- `POOL_PATHWAYS` 3 -> **5**
- pool is now the union of each pathway's **top 8 members** by signed score,
  rather than each pathway's entire membership (`PER_PATHWAY_TOP = 8`)

This bounds the pool at 40 genes instead of the several hundred that whole-
membership pooling produced.

**Hypothesis.** Experiment 2's residual failure is fold-level anchor
switching, not capacity or signal strength. Under whole-pathway pooling a
single pathway flip replaces the entire pool, and two of eight folds did
exactly that. Drawing a bounded quota from five pathways makes one flip cost
one fifth of the pool instead of all of it, which will raise Jaccard, cut the
repeat range, and pull repeat 2's increment up toward repeat 1's +1.470.

**Why this and not the alternatives.**

- *Not `benefit_threshold_months`.* It looked attractive because repeat 2's
  C+G policy treats 36.9% against the clinical 38.2%, but the arithmetic kills
  it: raising the threshold flips marginal patients from treat to no-treat,
  and each flip moves value by only `-tau/n`. With a cohort-mean AIPW benefit
  of about -0.2 months, trimming even 6% of patients buys on the order of
  0.01 months against a 0.583-month shortfall. It is two orders of magnitude
  too weak to be the fix, and it would make part of any "genomic" gain an
  artifact of the decision rule rather than the panel.
- *Not lower `max_features` or fewer columns.* Already at one column, the
  arena minimum. Experiment 2 showed the module's contribution is real but
  variable (+1.47 / -0.15); shrinking its influence shrinks the gain as well
  as the loss and parks both increments near zero, which fails a strictly-
  positive gate just as surely.
- *Not yet the expression-variability filter.* It would replace a panel that
  is already clearing the Jaccard gate. Held as the fallback.
- *Not larger pathways.* Raising `MIN_PATHWAY_MEMBERS` would enlarge the pool
  and push selection back toward a global ranking, undoing the mechanism that
  produced experiment 2's Jaccard win.

**Prespecified predictions and what would falsify them.**

- Jaccard rises above 0.1462, and the repeat range falls below 1.622.
- Both repeats' increments are positive. Repeat 2 has 0.153 months to make up
  and repeat 1 has 1.470 in hand, so the run passes if bounding the pool is
  worth even a fraction of the fold-heterogeneity it is aimed at.
- If Jaccard rises but the increments do not move together, then fold-level
  anchor switching was *not* what was driving the repeat gap, the remaining
  variability is in the forests rather than the selection, and the next axis
  is the expression-variability filter — changing *which* signal the module
  tracks, not how stably it is picked.
- If Jaccard falls, bounded quotas fragment the pool across families and
  whole-pathway pooling was load-bearing; revert to experiment 2's pooling
  and take the variability filter instead.

**Result.** (to be appended after the run)

**Result — run_003_20260823T191920Z, 226 s wall. ELIGIBLE. Reward -3.9105.**

**All 13 gates pass.** `best_run.txt` and `diagnostic_leader.txt` both now
name this run.

| field | repeat 1 | repeat 2 |
|---|---|---|
| incremental alignment (months) | **+2.574** | **+1.939** |
| selection LCB | -2.406 | **-3.275** |
| bootstrap mean / CI95 | +2.543 / [-0.904, +5.906] | +1.904 / [-1.624, +5.177] |
| alignment C / C+G | 5.098 / 7.672 | 5.361 / 7.300 |
| value C / C+G / best constant | 48.152 / 48.596 / 45.514 | 48.526 / 48.881 / 45.483 |
| Harrell C, C / C+G (drop) | 0.6561 / 0.6518 (0.0043) | 0.6641 / 0.6561 (0.0080) |
| ACT recommended fraction C+G | 0.3056 | 0.3540 |
| seed agreement / benefit corr | 0.9919 / 0.9991 | 0.9936 / 0.9989 |
| predicted benefit mean / IQR / nontrivial | -5.333 / 15.474 / 0.9942 | -3.925 / 15.218 / 0.9961 |
| genomic split fraction obs / ACT | 0.3794 / 0.2851 | 0.3852 / 0.2902 |
| arm support obs (pt/ev) | 661.5 / 278.25 | 661.5 / 278.25 |
| arm support ACT (pt/ev) | 114.0 / 70.5 | 114.0 / 70.5 |
| raw propensity overlap / IPTW ESS | 0.8723 / 358.5 | 0.8694 / 366.3 |

Reward **-3.9105** = robust LCB -3.2751 minus repeat range 0.6355.
Gene Jaccard 0.1417. Full-development module: `IFNA8, KCNC3, OR2S2, KCNV1,
TAAR5, ANO2, OR10H2, PTPN1, HTR6, OR2J2, OR3A2, OR52A1, OR10C1, OR2F1,
KCNJ9, IFNA21`.

**Scorecard against the prespecification.**

- *Both increments positive*: **confirmed**. +2.574 and +1.939. Repeat 2 moved
  +2.09 months, from -0.153 to +1.939 — far more than the 0.153 it needed.
- *Repeat range below 1.622*: **confirmed**, 0.6355, a 61% reduction. This is
  where most of the reward improvement came from.
- *Jaccard above 0.1462*: **falsified** — 0.1417, essentially unchanged and a
  hair lower.

That last falsification is the interesting one, because it means my stated
mechanism was wrong even though the prediction it was attached to came true.
Bounding the pool did **not** work by making gene selection more repeatable;
Jaccard is flat and the two histone/HLA folds (r1f3, r2f2) still switched
families. What it did was change the *composition* of every fold's panel:
with only eight genes drawn per pathway, each fold's 16 genes are now spread
across five anchors instead of being dominated by one, so even a fold that
switches its top anchor still shares its remaining four anchors with the
others. `KCNC3` now appears in 7 of 8 folds, against 5 of 8 for the best gene
in experiment 2. Fold-level *heterogeneity of the module* fell without
fold-level *agreement on gene identity* rising. Jaccard was simply the wrong
instrument for what mattered; the repeat range measured it directly.

**Interpretation.** The winning candidate is the locked clinical T-learner
plus a single averaged, sign-coherent Reactome module built from a bounded
five-pathway quota. It buys about 2.3 months of incremental policy alignment
in both repeated cross-fits, at a cost of under 0.008 in Harrell C. The
olfactory/potassium-channel caveat recorded under experiment 2 stands
unchanged and applies to this panel too.

**Status: an eligible candidate exists.** The reward is negative because the
objective is a 0.25%-quantile selection LCB (alpha/20 multiplicity), roughly
2.8 bootstrap SDs below the mean, and the increment's bootstrap SD is about
1.7 months. A positive reward would need a mean increment near 4.8 months
against the 2.3 achieved. The remaining experiments go to enlarging the
increment, since `best_run.txt` only advances on a strictly higher *eligible*
reward and the frozen winner cannot be lost.

## Experiment 4 (prespecified before running)

**Candidate** `tlearner_bidirectional_pool_module16`

Identical to the eligible experiment 3 — 16 genes, bounded five-pathway
quota of 8 members each, both arms at locked clinical geometry, threshold
0.0 — except that the panel is now **sign-stratified across both benefit
directions**:

- `module_count` 1 -> **2**
- the benefit-increasing half is anchored on the **top 5** pathways by mean
  signed score, the benefit-decreasing half on the **bottom 5**
- 8 genes per direction, positives returned first, so the locked
  `array_split` grouping yields one sign-pure positive module and one
  sign-pure negative module (verified on synthetic data before launch)

**Hypothesis.** Experiments 2 and 3 threw away half the available signal.
The locked ranking is on |correlation|, and every candidate so far has kept
only the positive tail, discarding the genes most strongly associated with
*reduced* benefit even though they are exactly as informative. Giving the
forest both directions supplies a genuine contrast — the benefit axis is
approximately module 1 minus module 2 — and should raise the mean increment
above experiment 3's +2.574 / +1.939.

**Why the mean increment is the right target.** Reward is a 0.25%-quantile
selection LCB minus the repeat range. The bootstrap SD of the increment is
about 1.7 months, so the LCB sits roughly 2.8 SDs below the mean and the mean
is worth nearly month-for-month in reward. The repeat range is already down
to 0.6355, capping any further gain from that term at well under a month.
Enlarging the increment is the only lever with several months of headroom.

**Costs I am knowingly accepting.**

- A second genomic column will raise the genomic split fraction from ~0.38
  toward ~0.5. Experiment 1 associated a high split share with harm, but that
  was four *raw, individually noisy* columns; these are averaged, sign-pure,
  pathway-anchored modules, and experiment 2 already showed a module column
  can raise the C-index rather than lower it.
- Each module now averages 8 genes rather than 16, so per-column noise rises
  by roughly sqrt(2). The bet is that a real second direction beats the lost
  averaging.

**Prespecified predictions and what would falsify them.**

- Mean increment exceeds +2.574 in at least one repeat and stays positive in
  both. If both increments fall below experiment 3's, the negative direction
  carries no usable signal at this sample size and the bidirectional axis is
  closed — the panel should go back to one 16-gene module and the next axis
  becomes *which* signal is tracked, not how much of it.
- Genomic split fraction lands near 0.5 and the C-index drop stays under
  0.03. If the drop approaches the gate, two columns is more capacity than
  this feature geometry supports and column count must return to one.
- Repeat range stays below 1.0. If it blows out, the negative direction is
  less stable across partitions than the positive one, and the fix is a
  bidirectional panel with an *unequal* split favouring the stabler
  direction.

Eligibility of `run_003` is not at risk: `best_run.txt` advances only on a
strictly higher eligible reward.

**Result.** (to be appended after the run)

**Result — run_004_20260824T000901Z, 282 s wall. NOT ELIGIBLE, reward -1e6.**

| field | repeat 1 | repeat 2 |
|---|---|---|
| incremental alignment (months) | **-2.003** | **+0.868** |
| selection LCB | -10.342 | -4.531 |
| bootstrap mean / CI95 | -2.039 / [-7.373, +3.232] | +0.843 / [-2.925, +4.264] |
| alignment C / C+G | 5.098 / 3.095 | 5.361 / 6.229 |
| value C / C+G / best constant | 48.152 / 45.948 / 45.448 | 48.526 / 47.916 / 45.489 |
| Harrell C, C / C+G (drop) | 0.6561 / 0.6534 (0.0027) | 0.6641 / 0.6475 (0.0166) |
| ACT recommended fraction C+G | 0.2959 | 0.3404 |
| seed agreement / benefit corr | 0.9926 / 0.9990 | 0.9900 / 0.9986 |
| predicted benefit mean / IQR / nontrivial | -6.010 / 16.961 / 0.9942 | -4.413 / 15.224 / 0.9913 |
| genomic split fraction obs / ACT | 0.5094 / 0.4798 | 0.5099 / 0.4094 |
| arm support obs (pt/ev) | 661.5 / 278.25 | 661.5 / 278.25 |
| arm support ACT (pt/ev) | 114.0 / 70.5 | 114.0 / 70.5 |
| raw propensity overlap / IPTW ESS | 0.8723 / 358.5 | 0.8694 / 366.3 |

Diagnostic score -13.214 (robust LCB -10.342, repeat range 2.871).
Jaccard **0.0784**. Gates 10/13 — `genomic_increment_positive`,
`genomic_value_at_least_clinical`, `gene_selection_jaccard_at_least_0_10`.

**Scorecard: the central prediction was falsified.**

- *Mean increment exceeds +2.574 somewhere, positive in both*: **falsified on
  both counts**. Both repeats came in below experiment 3 (-2.003 vs +2.574,
  +0.868 vs +1.939), and repeat 1 went negative.
- *Split fraction near 0.5, C-index drop under 0.03*: **confirmed**
  (0.51/0.48 and 0.0166). Capacity behaved exactly as modelled, so capacity
  is not what broke this run.
- *Repeat range below 1.0*: **falsified**, 2.871 — worse than any run since
  experiment 1.

**Why it failed, and it is not the reason I hedged for.** I accepted two
risks going in: more split share, and less averaging per module. Both
materialised and neither was decisive. The decisive fact is an **asymmetry
between the two benefit directions that I had assumed away**. Listing each
fold's negative half makes it obvious:

- r1f1 mitochondrial (`ATP5PF, MRPS15, IMMT, APOOL, LYRM4, ATP5ME`)
- r1f2 cytosolic ribosome (`RPS4Y1, RPS13, RPL30, RPS8, RPL14, RPSA`)
- r1f4 ER glycosylation (`DPM3, DDOST, CANX, MGAT2, MOGS, EDEM1`)
- r2f2 glycosylation + ATP synthase, r2f3 back to ribosome

The negative direction has **no stable anchor at all** — it lands on a
different housekeeping complex in nearly every fold. The positive direction
recurs (`KCNC3` in 6 of 8 folds, `IFNA8` in 4). So the two tails are not
mirror images: the positive tail is dominated by a coherent, highly
reproducible axis, and the negative tail is essentially fold-specific noise
drawn from whichever large housekeeping pathway happens to score lowest.
Bolting a noise module onto a signal module halved Jaccard, blew the repeat
range out to 2.871, and cost more than the second direction could ever repay.

This also sharpens the experiment 2 caveat rather than resolving it. A
positive tail that is stable across folds while the negative tail is pure
noise is what a *technical* axis looks like — low-expression probes moving
together with array background — not what a symmetric biological
effect-modification signal would look like.

**Prespecified consequence, applied.** My falsification rule for this run
said that if both increments fell below experiment 3, the bidirectional axis
is closed, the panel returns to a single 16-gene module, and the next axis
becomes *which* signal is tracked rather than how much of it. Both fell. The
axis is closed; experiment 5 reverts to the experiment 3 selector.

## Experiment 5 (prespecified before running)

**Candidate** `tlearner_wide_anchor_module16`

Reverts to the eligible experiment 3 selector exactly — single sign-coherent
16-gene module, positive direction only, both arms at locked clinical
geometry, threshold 0.0 — and changes only the quota shape:

- `POOL_PATHWAYS` 5 -> **8**
- `PER_PATHWAY_TOP` 8 -> **3** (pool bounded at 24 genes, 16 selected)

**Hypothesis.** Experiment 3's gain came from anchor spreading, not from gene
identity: Jaccard stayed flat at 0.14 while the repeat range fell 61%,
because a fold that switches its top anchor still shares the remaining
anchors with the other folds. Pushing the same mechanism further — eight
anchors contributing three slots each, so one anchor flip perturbs 3 of 16
panel slots instead of 8 of 16 — should cut the repeat range below 0.6355
and hold or raise the increments.

**Why this rather than the alternatives.**

- *Not the bidirectional axis.* Closed by experiment 4's prespecified rule.
- *Not more genomic columns.* Experiments 1 and 4 both show that added
  columns cost more than they return in this feature geometry.
- *Not `benefit_threshold_months`.* The arithmetic ruled it out before
  experiment 3 and nothing since has changed it.
- *Not yet the expression-variability filter.* That axis is now more
  interesting than it was — experiment 4's tail asymmetry is real evidence
  that the module tracks a technical low-expression axis — but it replaces
  the one signal that has ever cleared all thirteen gates. It is worth a slot
  only after the cheap, low-risk follow-through on the known mechanism is
  spent, because that follow-through is the likeliest source of a better
  reward and cannot plausibly lose eligibility.

**Prespecified predictions and what would falsify them.**

- Repeat range falls below 0.6355 and both increments stay positive, giving a
  reward above -3.9105.
- Jaccard is again roughly flat, near 0.14. If instead Jaccard moves sharply
  with the range, my reading of experiment 3's mechanism is wrong and the two
  quantities are coupled after all.
- If the repeat range *rises*, then eight anchors reaches far enough down the
  pathway ranking to admit unstable low-ranked pathways, the spreading
  mechanism has an optimum between 5 and 8 anchors, and the remaining quota
  tuning is exhausted — the next axis becomes the expression-variability
  filter.

**Result.** (to be appended after the run)

**Result — run_005_20260824T001447Z, 226 s wall. NOT ELIGIBLE, reward -1e6.**

| field | repeat 1 | repeat 2 |
|---|---|---|
| incremental alignment (months) | **+2.023** | **+1.146** |
| selection LCB | -3.895 | -4.568 |
| bootstrap mean / CI95 | +1.983 / [-2.227, +6.089] | +1.105 / [-2.882, +4.657] |
| alignment C / C+G | 5.098 / 7.121 | 5.361 / 6.507 |
| value C / C+G / best constant | 48.152 / **47.672** / 45.515 | 48.526 / 48.311 / 45.471 |
| Harrell C, C / C+G (drop) | 0.6561 / 0.6554 (0.0007) | 0.6641 / 0.6519 (0.0122) |
| ACT recommended fraction C+G | 0.2950 | 0.3385 |
| seed agreement / benefit corr | 0.9932 / 0.9991 | 0.9961 / 0.9991 |
| predicted benefit mean / IQR / nontrivial | -6.212 / 16.164 / 0.9942 | -4.062 / 15.760 / 0.9990 |
| genomic split fraction obs / ACT | 0.3890 / 0.2883 | 0.3725 / 0.2988 |
| arm support obs (pt/ev) | 661.5 / 278.25 | 661.5 / 278.25 |
| arm support ACT (pt/ev) | 114.0 / 70.5 | 114.0 / 70.5 |
| raw propensity overlap / IPTW ESS | 0.8723 / 358.5 | 0.8694 / 366.3 |

Diagnostic score -5.446 (robust LCB -4.568, repeat range **0.8775**).
Jaccard 0.1270. Gates **12/13** — the only failure is
`genomic_value_at_least_clinical`, in repeat 1 (47.672 against 48.152).

**Scorecard.**

- *Repeat range below 0.6355, reward above -3.9105*: **falsified.** The range
  rose to 0.8775 and the run is not eligible.
- *Jaccard roughly flat near 0.14*: **confirmed**, 0.1270. Jaccard and the
  repeat range again moved independently, which supports the experiment 3
  reading even though the prediction built on it failed.

**Why widening the anchor set backfired.** Eight anchors reach further down
the pathway ranking than five, and the pathways at ranks 6-8 are exactly the
ones whose mean signed score is closest to the noise floor. The frequency
table shows the contamination directly: `H2BC15` now appears in 5 of 8 folds,
the histone family that experiments 2 and 3 saw only as an occasional whole-
fold defection has become a permanent low-rank tenant of the pool. Spreading
risk across more anchors only helps while the added anchors are themselves
stable; past that it imports instability rather than diluting it. The
mechanism has an optimum, and five anchors of eight genes is nearer to it
than eight anchors of three.

**A new failure mode worth naming.** This run is the first where alignment
and value came apart: repeat 1 has C+G alignment 7.121 against clinical
5.098 — the largest alignment gap of the whole search — while its C+G *value*
falls 0.48 below clinical. The policy separates who benefits from who does
not better than the clinical rule does, yet the patients it actually selects
are not worth more. Alignment rewards discrimination, value rewards the
chosen set; a policy that recommends ACT to only 29.5% of patients against
the clinical rule's 37% can win the first and lose the second. Any future
candidate that buys alignment by shrinking the recommended set has to be
checked against the value gate, which does not move with it.

**Prespecified consequence, applied.** My rule for this run was that a rising
repeat range means the spreading mechanism has an optimum between 5 and 8
anchors, that quota tuning is exhausted, and that the next axis is the
expression-variability filter. The range rose. Experiment 6 returns the quota
to experiment 3's shape (5 anchors, 8 genes each) and changes which genes are
eligible to be scored at all.
