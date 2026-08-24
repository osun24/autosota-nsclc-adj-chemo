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

## Experiment 6 (prespecified before running)

**Candidate** `tlearner_dynamic_range_module16`

The experiment 3 configuration in every respect — single sign-coherent
16-gene module, 5 anchors of 8 genes, both arms at locked clinical geometry,
threshold 0.0 — with one addition ahead of everything else in the selector:

- **`MIN_DYNAMIC_RANGE_QUANTILE = 0.50`**: genes whose fit-partition
  interquartile range falls below the median are dropped before any gene is
  scored and before any pathway is ranked.

Fit-only, computed inside the selector on the fitting partition, no gene
named. Verified on synthetic data with deliberately floor-compressed genes:
the filter removed all of them and none reached the panel.

**Hypothesis.** Every candidate since experiment 2 has anchored on olfactory
receptor and potassium channel families. Olfactory receptors are not
plausibly expressed in lung tumour tissue, so their measured variation is
compressed against the array detection floor and is more likely tracking
background than biology. Experiment 4 supplied the corroborating evidence:
the positive tail is stable across folds while the negative tail is
fold-specific noise, which is the signature of a technical axis rather than a
symmetric biological effect. Requiring above-median dynamic range should move
selection onto genes carrying real expression signal, and if that signal is
biological rather than technical it should be at least as reproducible and
worth more in incremental alignment.

**This is a deliberate gamble on a different signal, not a refinement.** It
replaces the only panel that has ever cleared all thirteen gates. I am
spending a slot on it because the quota-tuning axis is exhausted by my own
prespecified rule, because the remaining reward headroom is in the mean
increment rather than the repeat range, and because `best_run.txt` advances
only on a strictly higher eligible reward — experiment 3 cannot be lost.

**Prespecified predictions and what would falsify them.**

- The panel no longer contains `OR*`/`TAAR*`/`HTR*` genes. If it does, the
  olfactory families survive an above-median dynamic range cut, the
  detection-floor reading is wrong, and the caveat recorded under experiment
  2 should be softened rather than acted on.
- Jaccard stays above 0.10 and both increments stay positive. If Jaccard
  collapses toward experiment 4's 0.078, then the *stability* of the
  olfactory anchor was load-bearing and the technical axis, whatever its
  provenance, is the only reproducible signal this cohort offers at n=1034.
  That would be a genuine negative result about the data rather than about
  the model, and it would close the "which signal" axis.
- Reward above -3.9105. If both increments stay positive but the reward does
  not beat experiment 3, the filter is neutral and the next axis is forest
  geometry for the C+G panel, which the search has not yet touched.

**Result.** (to be appended after the run)

**Result — run_006_20260824T002108Z, 232 s wall. NOT ELIGIBLE, reward -1e6.**

| field | repeat 1 | repeat 2 |
|---|---|---|
| incremental alignment (months) | **-5.092** | **-2.194** |
| selection LCB | -12.702 | -8.784 |
| bootstrap mean / CI95 | -5.073 / [-10.046, -0.663] | -2.212 / [-6.701, +1.945] |
| alignment C / C+G | 5.098 / **0.006** | 5.361 / 3.167 |
| value C / C+G / best constant | 48.152 / 45.145 / 45.634 | 48.526 / 47.257 / 45.766 |
| Harrell C, C / C+G (drop) | 0.6561 / 0.6576 (**-0.0015**) | 0.6641 / 0.6553 (0.0088) |
| ACT recommended fraction C+G | 0.3327 | 0.3743 |
| seed agreement / benefit corr | 0.9907 / 0.9990 | 0.9942 / 0.9991 |
| predicted benefit mean / IQR / nontrivial | -4.444 / 14.765 / 0.9942 | -3.604 / 16.849 / 0.9932 |
| genomic split fraction obs / ACT | 0.3771 / 0.3100 | 0.3730 / 0.3446 |
| arm support obs (pt/ev) | 661.5 / 278.25 | 661.5 / 278.25 |
| arm support ACT (pt/ev) | 114.0 / 70.5 | 114.0 / 70.5 |
| raw propensity overlap / IPTW ESS | 0.8723 / 358.5 | 0.8694 / 366.3 |

Diagnostic score -15.601 (robust LCB -12.702, repeat range 2.899).
Jaccard **0.0762**. Gates 9/13.

Full-development panel: `RPLP1, RPS15, CFL1, ERBB2, OASL, MET, RPL36A,
H2AC14, RPL37, RPS15A, PLXNA3, RHOB, H2AX, RPS23, IRF1, PTPN6`.

**Scorecard: prediction 1 confirmed, prediction 2 falsified, and together
they settle the question.**

- *No `OR*`/`TAAR*`/`HTR*` genes survive*: **confirmed**, and completely. Not
  one olfactory or biogenic-amine receptor appears in any of the eight folds.
  The olfactory families sit below the median fit-partition IQR, exactly as
  the detection-floor reading predicted. That reading was right.
- *Jaccard stays above 0.10, increments stay positive*: **falsified on both
  counts.** Jaccard collapsed to 0.0762 — statistically indistinguishable
  from experiment 4's 0.0784 — and both increments went sharply negative.

**The negative result this establishes.** Restricted to genes with real
expression dynamic range, the DR benefit signal in this cohort is **not
reproducible across folds** at n=1034 with 152 treated patients. Each fold
lands on a different biologically coherent story: r1f1 ribosomal proteins,
r1f2 growth-factor/ERBB signalling (`HBEGF, EREG, ERBB2, MET, IL6, CXCL2`),
r1f3 interferon response (`IRF1, IRF9, STAT1, OASL, BST2, IFI27`), r2f1
histones. These are plausible-looking panels — which is precisely the danger.
They are what selection on noise produces when the candidate pool is
biologically structured, and the Jaccard gate is what exposes them.

Per my prespecified rule, this closes the "which signal" axis: the
low-dynamic-range axis is the only *reproducible* benefit signal this cohort
offers, and the caveat recorded under experiment 2 stands as a limitation of
the winning candidate rather than as a defect to be engineered away.

**A dissociation worth recording.** In repeat 1 the real-expression panel
*improved* the Harrell C-index (0.6576 against the clinical 0.6561) while
incremental alignment collapsed to 0.006 — the clinical and C+G policies
became effectively equivalent in discrimination. Genes that sharpen
*prognosis* are not the genes that identify *who benefits from treatment*.
The arena grades the second and the search must not be seduced by the first.

**Axis status after six experiments.**

- Panel size / representation: settled (16 genes, one sign-coherent module).
- Genomic column count: settled at 1 (experiments 1 and 4 both punished more).
- Bidirectional selection: closed by experiment 4's rule.
- Quota shape: closed by experiment 5's rule (optimum at 5 anchors x 8).
- Which signal is tracked: closed by this run's rule.
- `benefit_threshold_months`: ruled out arithmetically, twice.
- **Forest geometry for the C+G panel: never touched.** This is the only
  substantive axis left, and it is the one the program names as a preferred
  first axis in the form of ACT-arm regularization.

## Experiment 7 (prespecified before running)

**Candidate** `tlearner_regularized_act_module16`

The eligible experiment 3 configuration — single sign-coherent 16-gene
module, 5 anchors of 8 genes, threshold 0.0, observation arm at locked
clinical geometry — with one change, on the only substantive axis the search
has not touched:

- ACT arm depth 7 -> **6**, leaf 12 -> **20**, split 24 -> **40**
  (`max_features` held at 1.0)

**Hypothesis.** The ACT arm carries ~114 patients and ~70 events per fitting
partition against the observation arm's ~662 and ~278. Since
`benefit = RMST1 - RMST0`, ACT-arm variance dominates the benefit estimate,
and it is the most plausible source of the repeat-to-repeat spread that the
reward penalises directly. Regularizing the small arm harder while leaving
the well-supported observation arm at locked geometry should cut the repeat
range below 0.6355 and lift the binding LCB above -3.2751, beating
experiment 3's reward of -3.9105.

**Why this axis now.** Every other axis is closed by a prespecified rule:
panel size and representation are settled, column count is settled at one,
the bidirectional axis fell to experiment 4, quota shape to experiment 5,
which-signal to experiment 6, and `benefit_threshold_months` to arithmetic.
ACT-arm regularization is also the axis `program.md` names as preferred, and
experiment 1's confirmation that seed stability is a non-issue means the
regularization budget can be spent entirely on bias/variance.

**Why this is a variance move and not a capacity cut.** The genomic signal is
not being removed, it is being relocated. The observation arm keeps full
locked geometry and its ~0.38 genomic split fraction, so the module still
drives RMST0 and therefore still drives benefit. What changes is that the
114-patient arm stops trying to resolve genomic structure it does not have
the support to estimate.

**A calibration note on what a "better" result now means.** Experiments 3 and
5 are close relatives — one module, pathway-anchored, differing only in quota
shape — and their increments were +2.574/+1.939 and +2.023/+1.146. That
spread of roughly 0.8 months between near-identical designs is the search's
own noise floor. Any reward improvement here smaller than that should be read
as partition noise, not as evidence the change worked, and I will say so
rather than bank it.

**Prespecified predictions and what would falsify them.**

- Repeat range falls below 0.6355 and the run stays eligible, with a reward
  above -3.9105 by more than the ~0.8-month noise floor.
- Both increments stay positive. If they fall materially, the ACT arm's
  genomic resolution was load-bearing after all, coarse RMST1 collapses
  benefit toward a purely prognostic `-RMST0` rule, and ACT geometry should
  return to locked clinical.
- Seed agreement and benefit correlation stay far above their gates; if
  either moves at all, something other than variance reduction happened.

**Result.** (to be appended after the run)

**Result — run_007_20260824T002722Z, 229 s wall. NOT ELIGIBLE, reward -1e6.**

| field | repeat 1 | repeat 2 |
|---|---|---|
| incremental alignment (months) | **+0.452** | **-0.082** |
| selection LCB | -4.904 | -5.914 |
| bootstrap mean / CI95 | +0.423 / [-3.111, +3.888] | -0.096 / [-3.924, +3.282] |
| alignment C / C+G | 5.098 / 5.550 | 5.361 / 5.279 |
| value C / C+G / best constant | 48.152 / 47.771 / 45.514 | 48.526 / 48.188 / 45.615 |
| Harrell C, C / C+G (drop) | 0.6561 / 0.6495 (0.0066) | 0.6641 / 0.6534 (0.0107) |
| ACT recommended fraction C+G | 0.2930 | 0.3240 |
| seed agreement / benefit corr | 0.9929 / 0.9991 | 0.9955 / 0.9991 |
| predicted benefit mean / IQR / nontrivial | -5.133 / 15.050 / 0.9971 | -4.197 / 14.415 / 0.9971 |
| genomic split fraction obs / ACT | 0.3794 / 0.3622 | 0.3852 / 0.3314 |
| arm support obs (pt/ev) | 661.5 / 278.25 | 661.5 / 278.25 |
| arm support ACT (pt/ev) | 114.0 / 70.5 | 114.0 / 70.5 |
| raw propensity overlap / IPTW ESS | 0.8723 / 358.5 | 0.8694 / 366.3 |

Diagnostic score -6.448 (robust LCB -5.914, repeat range **0.5341**).
Gates 11/13 — `genomic_increment_positive`, `genomic_value_at_least_clinical`.

**Scorecard: the variance prediction came true and bought nothing.**

- *Repeat range below 0.6355*: **confirmed**, 0.5341 — the lowest of the whole
  search. The variance mechanism was real.
- *Both increments stay positive*: **falsified.** They collapsed from
  +2.574/+1.939 to +0.452/-0.082. That is a drop of roughly 2.1 and 2.0
  months, far outside the ~0.8-month noise floor I set in advance, so this is
  a real effect and not partition noise.
- *Seed diagnostics unmoved*: confirmed (0.9991 both repeats).

**What this establishes.** The ACT arm's genomic resolution is load-bearing,
despite that arm having only ~114 patients and ~70 events. Coarsening it
traded almost exactly two months of incremental alignment for one tenth of a
month of repeat range — a catastrophic exchange rate given the reward is
`LCB - range`. Note that the ACT genomic split *fraction* actually rose
(0.285 -> 0.362) while its contribution collapsed: with coarser trees there
are fewer splits overall and the one continuous column takes a larger share
of them. Split fraction measures how often the module is used, not how much
it resolves, and the two came apart here for the first time.

Per my prespecified rule, ACT geometry returns to locked clinical.

**The more useful reading is the inverse.** Reducing resolution cost two
months of increment. Nothing in this search has ever tested the other
direction: every candidate has used the locked clinical geometry or something
coarser. If benefit discrimination is resolution-limited rather than
noise-limited, a *finer* C+G forest should move the increment the other way.
That is the hypothesis for experiment 8, and it is the first one in several
runs that is derived from a measured effect rather than from a diagnosis of
what went wrong.

## Experiment 8 (prespecified before running)

**Candidate** `tlearner_fine_resolution_module16`

The eligible experiment 3 configuration — single sign-coherent 16-gene
module, 5 anchors of 8 genes, threshold 0.0 — with **finer** forest geometry
in both arms of the C+G panel:

- observation arm depth 9 -> **12**, leaf 8 -> **5**, split 16 -> **10**
- ACT arm depth 7 -> **8**, leaf 12 -> **10**, split 24 -> **20**

`max_features` stays 1.0 in both arms; the clinical comparator is untouched.

**Hypothesis.** Experiment 7 measured the exchange rate in the coarse
direction: two months of incremental alignment lost for a tenth of a month of
repeat range gained. That is the signature of a **resolution-limited**
problem, not a noise-limited one. Every candidate in this search has sat at
the locked clinical geometry or coarser, so the fine direction has never been
sampled. If the exchange rate is even roughly symmetric near the locked
point, a finer C+G forest should raise both increments above experiment 3's
+2.574 / +1.939 at a modest cost in repeat range, and the reward is
`LCB - range` with the LCB term about 2.8 times more sensitive.

**Why the asymmetry with the clinical comparator is legitimate and what it
costs interpretively.** Only the clinical T-learner is locked; the C+G panel's
geometry is explicitly editable, and red line 9 forbids weakening or
threshold-tuning the *clinical* learner, which this does not do. But it does
mean the C+G panel now has more capacity on the clinical features too, so any
gain is no longer attributable to the module alone. If this run wins, the
honest statement is that the winning policy uses both a Reactome module and a
finer forest, and the module's isolated contribution is the one measured by
experiment 3, where the geometries matched exactly.

**Prespecified predictions and what would falsify them.**

- At least one increment exceeds +2.574 and both stay positive, with the
  repeat range staying under about 1.0, giving a reward above -3.9105 by more
  than the ~0.8-month noise floor.
- The C-index drop stays inside 0.03. A depth-12, leaf-5 forest on ~662
  patients is the most overfittable configuration the arena permits, and the
  C-index gate is the natural place for that to show first.
- If the increments *fall*, then the locked clinical geometry is at or near
  the optimum rather than on the coarse side of it, experiment 7's exchange
  rate was one-sided, and the forest-geometry axis closes in both directions.
  With every other axis already closed by a prespecified rule, that would
  leave experiment 3 as the search's final answer and the remaining slots
  would go to confirming it rather than to further tuning.

**Result.** (to be appended after the run)

**Result — run_008_20260824T003334Z, 256 s wall. NOT ELIGIBLE, reward -1e6.**

| field | repeat 1 | repeat 2 |
|---|---|---|
| incremental alignment (months) | **+2.687** | **+1.401** |
| selection LCB | -2.518 | -4.401 |
| bootstrap mean / CI95 | +2.640 / [-0.890, +6.118] | +1.361 / [-2.468, +4.924] |
| alignment C / C+G | 5.098 / 7.785 | 5.361 / 6.762 |
| value C / C+G / best constant | 48.152 / 48.538 / 45.513 | 48.526 / **48.432** / 45.484 |
| Harrell C, C / C+G (drop) | 0.6561 / 0.6506 (0.0055) | 0.6641 / 0.6570 (0.0071) |
| ACT recommended fraction C+G | 0.3104 | 0.3656 |
| seed agreement / benefit corr | 0.9916 / 0.9989 | 0.9903 / 0.9987 |
| predicted benefit mean / IQR / nontrivial | -5.533 / 16.225 / 0.9942 | -3.951 / 16.365 / 0.9971 |
| genomic split fraction obs / ACT | 0.4036 / 0.2830 | 0.4094 / 0.2880 |
| arm support obs (pt/ev) | 661.5 / 278.25 | 661.5 / 278.25 |
| arm support ACT (pt/ev) | 114.0 / 70.5 | 114.0 / 70.5 |
| raw propensity overlap / IPTW ESS | 0.8723 / 358.5 | 0.8694 / 366.3 |

Diagnostic score -5.686 (robust LCB -4.401, repeat range 1.2851).
Jaccard 0.1417. Gates **12/13** — the only failure is
`genomic_value_at_least_clinical`, in repeat 2, by **0.094 months**.

**Scorecard.**

- *An increment above +2.574, both positive, range under 1.0, reward better
  than -3.9105*: **falsified.** Repeat 1 did reach +2.687, but by 0.11 months
  — inside the design-sensitivity scale — while repeat 2 fell to +1.401 and
  the range rose to 1.2851. Not eligible, so no reward.
- *C-index drop inside 0.03*: **confirmed**, and more comfortably than
  expected (0.0055 / 0.0071). A depth-12, leaf-5 forest on ~662 patients is
  the most overfittable configuration the arena permits and it did not
  degrade discrimination at all. Whatever limits this problem, it is not
  observation-arm overfitting.

**The geometry axis closes, from both sides.** Coarser (experiment 7) bought
0.10 months of repeat range for 2.0 months of increment. Finer (this run)
bought at most 0.11 months of increment for 0.65 months of repeat range.
Neither direction pays, and the locked clinical geometry sits at or very near
the optimum for the C+G panel as well. That is a mildly interesting result in
its own right: the geometry chosen for the locked comparator turns out to be
about right for a panel with one extra column.

**A correction to how I have been describing the "noise floor".** This
pipeline is fully deterministic — fixed fold seeds, fixed forest seeds, a
fixed bootstrap seed, and a deterministic selector — so re-running any
candidate reproduces its numbers exactly. The ~0.8-month spread I have been
attributing to noise is therefore not sampling variability and cannot be
averaged away. It is *sensitivity to scientifically arbitrary design
choices*: quota 5x8 versus 8x3 is not a hypothesis about biology, yet it
moves the increment by that much. For model selection it functions like
noise, which is exactly what the 0.25% selection LCB is priced to punish, but
it should be named accurately.

**What the last three runs jointly imply.** Experiments 3, 5 and 8 share a
selector and differ only in arbitrary design details, and their value-gate
margins were +0.44/+0.35, -0.48/-0.21, and +0.39/-0.09. The design family
clears 12 gates dependably and the thirteenth — value — sits within about
half a month of the line either way. Experiment 3's eligibility is real and
exactly reproducible, but it is not a comfortable margin, and further tuning
of arbitrary knobs would be selection on that margin rather than improvement
of it.

## Experiment 9 (prespecified before running)

**Candidate** `tlearner_precise_anchor_module16`

The eligible experiment 3 configuration in every respect — single
sign-coherent 16-gene module, 5 anchors of 8 genes, both arms at locked
clinical geometry, threshold 0.0 — with one change:

- `MIN_PATHWAY_MEMBERS` 25 -> **50**

**This is a different knob from the one experiment 5 closed.** Experiment 5's
rule retired *quota shape* — how many anchors contribute and how many genes
each contributes. This changes *which pathways are eligible to be anchors at
all*, leaving the 5x8 quota exactly as experiment 3 had it. I want that
distinction on the record rather than glossed, because it is the difference
between honouring a prespecified closure and quietly reopening it.

**Hypothesis.** A pathway's mean signed score has standard error falling as
1/sqrt(members), so a floor of 50 makes every eligible anchor's score about
1.4 times more precise and the ranking among anchors correspondingly more
reliable. Experiment 3 works because its anchors are consistent across folds;
experiment 5 failed because adding anchors reached down into the noise floor.
Raising the precision bar attacks the same mechanism from the opposite side —
better anchors rather than more of them — and should widen the value-gate
margin, which is the gate the design family keeps failing by half a month or
less.

**What I am no longer expecting.** After eight runs the increment for this
selector family sits between +1.1 and +2.7 regardless of what I change, and
the value margin sits within about half a month of the line either way. I do
not expect this run to transform either. It is worth a slot because it is the
last untested knob with a real statistical argument behind it, and because
the value margin is the only thing standing between this family and a
comfortable rather than a marginal win.

**Prespecified predictions and what would falsify them.**

- Both value margins are positive and larger than experiment 3's +0.44 /
  +0.35, and the run is eligible.
- Jaccard holds near or above 0.14, since restricting to large pathways
  should if anything concentrate selection further.
- If the value margins do not improve, then anchor precision is not what
  limits this design, the selector axis is exhausted along with every other,
  and experiment 3 stands as the search's answer. In that case the remaining
  slots should not be spent on further knob-turning, which would be selection
  on a half-month margin rather than improvement of it. I would stop and
  report, rather than consume the budget for its own sake.

**Result.** (to be appended after the run)

**Result — run_009_20260824T003957Z, 233 s wall. NOT ELIGIBLE, reward -1e6.**

| field | repeat 1 | repeat 2 |
|---|---|---|
| incremental alignment (months) | **+4.224** | **+0.902** |
| selection LCB | **-0.695** | -5.535 |
| bootstrap mean / CI95 | +4.231 / **[+0.694, +8.134]** | +0.875 / [-3.504, +4.667] |
| alignment C / C+G | 5.098 / **9.322** | 5.361 / 6.263 |
| value C / C+G / best constant | 48.152 / 48.934 / 45.380 | 48.526 / **48.011** / 45.430 |
| value margin vs clinical | **+0.781** | **-0.515** |
| Harrell C, C / C+G (drop) | 0.6561 / 0.6534 (0.0027) | 0.6641 / 0.6547 (0.0094) |
| ACT recommended fraction C+G | 0.3056 | 0.3656 |
| seed agreement / benefit corr | 0.9939 / 0.9992 | 0.9945 / 0.9991 |
| predicted benefit mean / IQR / nontrivial | -5.293 / 15.581 / 0.9932 | -4.276 / 16.086 / 0.9952 |
| genomic split fraction obs / ACT | 0.3814 / 0.2750 | 0.3882 / 0.3130 |
| arm support obs (pt/ev) | 661.5 / 278.25 | 661.5 / 278.25 |
| arm support ACT (pt/ev) | 114.0 / 70.5 | 114.0 / 70.5 |
| raw propensity overlap / IPTW ESS | 0.8723 / 358.5 | 0.8694 / 366.3 |

Diagnostic score -8.856 (robust LCB -5.535, repeat range **3.3216**).
Jaccard 0.1378. Gates 12/13 — `genomic_value_at_least_clinical`, repeat 2.

**Scorecard: falsified, but not in the way the rule anticipated.**

- *Both value margins positive and above +0.44/+0.35*: **falsified.** Repeat 1
  improved markedly (+0.781 against +0.44) while repeat 2 went negative
  (-0.515). The margins did not fail to move — they moved *apart*.
- *Jaccard near 0.14*: confirmed, 0.1378.

Repeat 1 produced the best single-repeat result of the entire search:
increment +4.224, alignment 9.322 against the clinical 5.098, and a bootstrap
CI95 of **[+0.694, +8.134]** that excludes zero — the only repeat in nine
experiments whose ordinary 95% interval is entirely positive. Repeat 2 gave
+0.902. The repeat range tripled to 3.3216.

**The mechanism is the opposite of the one I proposed, and it matters.** I
argued that more precise anchor scores would stabilize the panel. Instead,
raising the membership floor to 50 *shrank the eligible anchor set*, leaving
only the largest gene families to compete. Selection concentrated rather than
sharpened, and a concentrated panel tracks whatever its single family happens
to be worth in a given partition — spectacular in repeat 1, mediocre in
repeat 2. **Anchor diversity buys stability; anchor precision buys
sensitivity.** That is now the third independent confirmation of the same
mechanism, from a third direction: experiment 3 gained by spreading across
anchors, experiment 5 lost by spreading into anchors too weak to help, and
experiment 9 lost by removing the spread altogether.

**Why this does not trigger my stop rule, and what it points at instead.**
The rule was to stop if the value margins "do not improve". Repeat 1's
improved substantially; the failure was divergence, not inertia. And the
diagnosis raises a specific, previously unexamined question: *are experiment
3's five anchors actually five distinct gene sets?* Reactome is hierarchical,
and a family such as olfactory signalling appears as a chain of nested
pathways — the receptor set, its GPCR parent, that parent's signalling
parent, and so on — every one of which is large enough to pass a member
floor and all of which score alike because they contain largely the same
genes. If the top five anchors are nested relatives, then the quota has been
buying redundancy rather than diversity all along, and every result in this
family has been produced by an effective anchor count closer to one than to
five. That would explain why raising the floor to 50 concentrated things so
sharply: it retains exactly the large nested parents.

Experiment 10 tests it directly by de-duplicating overlapping anchors, which
is selector logic rather than another knob.

## Experiment 10 (prespecified before running)

**Candidate** `tlearner_distinct_anchor_module16`

The eligible experiment 3 configuration — single sign-coherent 16-gene
module, 5 anchors of 8 genes, `MIN_PATHWAY_MEMBERS` back to 25, both arms at
locked clinical geometry, threshold 0.0 — with one addition to the selector:

- **`MAX_ANCHOR_OVERLAP = 0.50`**: a pathway is skipped as an anchor when its
  development-present membership overlaps an already-accepted anchor by more
  than Jaccard 0.5. The quota then backfills from further down the ranking
  until five genuinely distinct anchors are held.

**Hypothesis.** Reactome is hierarchical. A receptor set, its GPCR parent and
that parent's signalling parent contain largely the same genes, so they score
almost identically and a naive "top 5 pathways" rule can return five nested
relatives instead of five distinct gene sets. If that is what has been
happening, every result in this family has been produced by an effective
anchor count near one, and the diversification mechanism that experiments 3,
5 and 9 all point at has never actually been switched on. De-duplicating
anchors should raise the *effective* diversity that experiment 3 only
appeared to have, narrowing the repeat range and lifting the binding repeat-2
LCB above -3.2751.

**Verified before launch.** On a synthetic nested chain constructed so the
family sweeps the ranking, the naive rule takes 5 of 5 anchors from the
nested family; de-duplication takes 1 and backfills four distinct sets. The
mechanism fires as intended, and all three fallback paths still return 16
valid unique genes.

**Prespecified predictions and what would falsify them.**

- The full-development panel is no longer dominated by a single gene family.
  If it comes back looking like experiment 3's — `OR*`/`KCN*` throughout —
  then the top anchors were already distinct pathways that merely happen to
  share a biology, the nesting diagnosis is wrong, and this run is a null
  that costs a slot and settles the question.
- Repeat range below 0.6355 and repeat-2 LCB above -3.2751, giving a reward
  better than -3.9105.
- If the panel does diversify but the reward does not improve, then effective
  anchor diversity is not what limits this design either. Combined with the
  closure of the capacity, representation, direction, quota, signal and
  geometry axes, that would exhaust every mechanism I can motivate from the
  evidence, and I will stop and report experiment 3 rather than spend the
  remaining slots on knobs whose effects are smaller than the design
  sensitivity documented under experiment 8.

**Result.** (to be appended after the run)

**Result — run_010_20260824T004651Z, 233 s wall. NOT ELIGIBLE, reward -1e6.**

| field | repeat 1 | repeat 2 |
|---|---|---|
| incremental alignment (months) | **+3.403** | **-0.327** |
| selection LCB | -1.882 | -6.484 |
| bootstrap mean / CI95 | +3.355 / [-0.162, +6.734] | -0.344 / [-4.584, +3.339] |
| alignment C / C+G | 5.098 / 8.502 | 5.361 / 5.034 |
| value C / C+G / best constant | 48.152 / 48.521 / 45.543 | 48.526 / 47.713 / 45.446 |
| Harrell C, C / C+G (drop) | 0.6561 / 0.6528 (0.0034) | 0.6641 / 0.6543 (0.0097) |
| ACT recommended fraction C+G | 0.2911 | 0.3395 |
| seed agreement / benefit corr | 0.9929 / 0.9991 | 0.9936 / 0.9989 |
| predicted benefit mean / IQR / nontrivial | -6.352 / 15.670 / 0.9952 | -4.036 / 14.852 / 0.9990 |
| genomic split fraction obs / ACT | 0.3892 / 0.2742 | 0.3837 / 0.2945 |
| arm support obs (pt/ev) | 661.5 / 278.25 | 661.5 / 278.25 |
| arm support ACT (pt/ev) | 114.0 / 70.5 | 114.0 / 70.5 |
| raw propensity overlap / IPTW ESS | 0.8723 / 358.5 | 0.8694 / 366.3 |

Diagnostic score -10.215 (robust LCB -6.484, repeat range 3.7305).
Jaccard 0.1412. Gates 11/13.

**Scorecard, and an error in my own test design.** The prediction was that the
panel would stop being dominated by one gene family. The full-development
panel came back **byte-identical to experiment 3's** — same sixteen genes in
the same order — and five of the eight folds were unchanged. By the letter of
my rule this is the null branch. But the reason is a flaw in how I measured
overlap, and it should be recorded as mine rather than as a fact about the
data: **Jaccard is the wrong statistic for detecting nesting.** A 60-member
child inside a 400-member parent has Jaccard 60/400 = 0.15, far below my 0.50
threshold, so the guard almost never fired on exactly the hierarchy it was
built to catch. Containment, `|A and B| / min(|A|, |B|)`, would have been the
right measure. My synthetic test missed this because the nested chain I built
had members of similar size, where Jaccard and containment nearly coincide.

**What the run nevertheless establishes, and it is the most useful result of
the second half of the search.** De-duplication fired in only three of eight
folds and changed between two and five genes in each:

- r1f3: 13 of 16 shared with experiment 3
- r2f1: 14 of 16 shared
- r2f2: 11 of 16 shared

Those small perturbations — under 10% of the total gene-slots across the
search — moved the increments from +2.574/+1.939 to +3.403/-0.327 and the
repeat range from 0.6355 to **3.7305**. A five-fold swing in the term the
reward subtracts, from changing a handful of genes in three folds.

This quantifies the design-sensitivity point from experiment 8 far more
sharply than experiment 8 did, and it settles the question of whether to keep
tuning. The objective surface is dominated by sensitivity to essentially
arbitrary micro-variation in gene selection, not by any mechanism that can be
steered. It also means a containment-based de-duplication is not worth a
slot: where the guard did fire it made things worse, and a stricter measure
would only fire more often.

**Why a positive reward is out of reach for this family.** Reward is
`LCB - range`, and the LCB is the 0.25% bootstrap quantile, about 2.84 SDs
below the mean. A positive LCB therefore needs `mean/SD > 2.84`. Experiment
3's best repeat runs at 2.543/1.74 = 1.46, and the ratio scales as
`sqrt(k) * E[tau|disagree] / sd(tau|disagree)` over the k patients where the
two policies disagree. Doubling it needs four times as many disagreeing
patients at the same per-patient signal, or twice the per-patient signal.
Neither is reachable by any knob in the editable surface. The eligible reward
is negative by construction of the multiplicity correction, not by a failure
of the search.

## Search status after 10 of 20 experiments

| run | candidate | elig | score | LCB | range | gates | jaccard | increments |
|---|---|---|---|---|---|---|---|---|
| 001 | raw4 dr_gene locked geometry | no | -16.027 | -13.116 | 2.911 | 9/13 | 0.0697 | -1.660 / -4.571 |
| 002 | pathway benefit module16 | no | -7.764 | -6.142 | 1.622 | 11/13 | 0.1462 | +1.470 / -0.153 |
| **003** | **bounded pool module16** | **YES** | **-3.911** | **-3.275** | **0.635** | **13/13** | **0.1417** | **+2.574 / +1.939** |
| 004 | bidirectional pool module16 | no | -13.214 | -10.342 | 2.871 | 10/13 | 0.0784 | -2.003 / +0.868 |
| 005 | wide anchor module16 | no | -5.446 | -4.568 | 0.878 | 12/13 | 0.1270 | +2.023 / +1.145 |
| 006 | dynamic range module16 | no | -15.601 | -12.702 | 2.899 | 9/13 | 0.0762 | -5.092 / -2.194 |
| 007 | regularized ACT module16 | no | -6.448 | -5.914 | 0.534 | 11/13 | 0.1417 | +0.452 / -0.082 |
| 008 | fine resolution module16 | no | -5.686 | -4.401 | 1.285 | 12/13 | 0.1417 | +2.686 / +1.401 |
| 009 | precise anchor module16 | no | -8.856 | -5.535 | 3.322 | 12/13 | 0.1378 | +4.224 / +0.902 |
| 010 | distinct anchor module16 | no | -10.215 | -6.484 | 3.731 | 11/13 | 0.1412 | +3.403 / -0.327 |

**Every axis is now closed by a prespecified rule**: capacity and column
count (001, 004), representation and sign coherence (002), quota shape (005),
which signal is tracked (006), ACT-arm and observation-arm geometry in both
directions (007, 008), anchor precision (009), anchor de-duplication (010),
and `benefit_threshold_months` by arithmetic, twice.

**Recommendation: stop and freeze `run_003_20260823T191920Z`.** It is the
only eligible run, it holds both `best_run.txt` and `diagnostic_leader.txt`,
and it is exactly reproducible. The ten unspent experiments should not be
consumed for their own sake: experiment 10 showed the objective moves by
several months on gene-level perturbations that carry no hypothesis, so
further attempts would be selection on a margin the 0.25% LCB exists to
punish, and each one risks nothing but also teaches nothing.

**Limitations carried forward to any test-set decision.** The panel is
dominated by olfactory-receptor and potassium-channel genes whose measured
variation sits near the array detection floor (experiment 2), the benefit
direction is asymmetric in a way consistent with a technical rather than
biological axis (experiment 4), and every gene set with real expression
dynamic range failed to reproduce across folds (experiment 6). The candidate
clears the arena's thirteen gates; it should not be described as a
biologically interpretable pathway result. Test data remains untouched.

## Experiment 11 (prespecified before running)

**Goal change.** The directive is now an eligible candidate with a *positive*
reward, not merely an eligible one. That requires the 0.25% bootstrap
quantile of the increment to be positive in **both** repeats, and then to
exceed the repeat range. Concretely: a mean increment near 5 months at
SD ~1.7 in both repeats, against a best-ever repeat-2 increment of +1.939
(experiment 3). Experiment 9's repeat 1 is the closest any run has come —
mean +4.224, LCB -0.695, with an ordinary CI95 of [+0.694, +8.134] that
already excludes zero — so the target is not obviously unreachable, but it
needs both repeats to behave like that one.

**Candidate** `tlearner_stability_selected_module16`

Experiment 3's configuration — single sign-coherent 16-gene module,
5 anchors of 8 genes, `MIN_PATHWAY_MEMBERS` 25, both arms at locked clinical
geometry, threshold 0.0 — with the single-shot ranking replaced by
**stability selection**:

- the cross-fitted DR benefit pseudo-outcome is computed once per fitting
  partition, as before
- the entire pathway-anchored, sign-coherent pick is then repeated on
  **40 random 80% subsamples** of that partition
- the 16 genes chosen most often across those 40 picks form the panel, with
  mean signed score breaking ties

Fit-only throughout, fixed seed, no gene hard-coded.

**Hypothesis.** Every one of the ten completed experiments points at the same
binding defect, and none of them has addressed it directly. Experiment 1
measured Jaccard 0.070 for a top-4 of 8,647 genes. Experiment 6 showed that
each fold lands on a different plausible biology when the pool is
biologically structured. Experiment 10 showed that changing 2-5 genes in 3 of
8 folds swings the repeat range by 3.1 months. The panel is being chosen by a
statistic whose sampling variability dominates its signal. Selection
frequency across subsamples is a far lower-variance statistic than any single
ranking — that is precisely what stability selection is for — so it should
raise Jaccard well above the 0.14 ceiling this family has never beaten, pull
the two repeats toward each other, and lift the binding repeat-2 LCB.

**A calibration I obtained before launching, which changes how I read
Jaccard.** Running the single-draw pathway-anchored pick on synthetic pure
noise at full scale (8,647 genes, 1,839 pathways) gives a pairwise panel
Jaccard of **0.0809**. Pathway anchoring manufactures apparent agreement even
when there is no signal at all, because the pool restriction alone forces
overlap. So 0.08 is this selector's *noise floor*, not zero, and the 0.10
gate sits barely above it. That reframes several earlier results: experiment
4 (0.0784), experiment 6 (0.0762) and experiment 10's perturbed folds were
not merely unstable, they were **indistinguishable from pure noise
selection**. It also means experiment 3's 0.1417 is a real but modest signal —
roughly 1.75x the noise floor — and that a large Jaccard gain here would be
strong evidence the method is working rather than a cosmetic improvement.

**Timing.** Measured at full scale on synthetic data: 5.0 s per selector
call, so about 45 s for the nine calls, against a 2,100 s budget with roughly
150 s of forests. No timeout risk.

**Prespecified predictions and what would falsify them.**

- Jaccard rises above 0.20, comfortably clear of the 0.08 noise floor. If it
  does not move above 0.14, subsample-to-subsample variation is not what
  drives fold-to-fold panel disagreement — the folds differ by patients, not
  by sampling noise within a fixed partition — and stability selection cannot
  fix it.
- Both increments positive with the repeat range below 0.6355, and the reward
  above -3.9105. This is the eligibility-preserving prediction.
- Reward positive. I state this as the goal rather than as an expectation:
  it needs roughly a 2.5x improvement in the binding repeat-2 increment, and
  nothing in ten experiments has moved that quantity by more than about a
  month. If Jaccard rises sharply and the increments do not follow, the
  conclusion is that panel instability was never the limiting factor for the
  *increment* — only for the range — and the positive-reward target is out of
  reach for this estimand at n=1034 with 152 treated patients.

**Result.** (to be appended after the run)

**Result — run_011_20260824T005914Z, 257 s wall. NOT ELIGIBLE, reward -1e6.**

| field | repeat 1 | repeat 2 |
|---|---|---|
| incremental alignment (months) | **+3.587** | **-0.104** |
| selection LCB | -1.140 | -6.077 |
| bootstrap mean / CI95 | +3.558 / **[+0.483, +6.664]** | -0.161 / [-4.268, +3.343] |
| alignment C / C+G | 5.098 / 8.685 | 5.361 / 5.257 |
| value C / C+G / best constant | 48.152 / 48.298 / 45.506 | 48.526 / 47.934 / 45.430 |
| Harrell C, C / C+G (drop) | 0.6561 / 0.6569 (**-0.0008**) | 0.6641 / 0.6496 (0.0145) |
| ACT recommended fraction C+G | 0.3008 | 0.3443 |
| seed agreement / benefit corr | 0.9939 / 0.9992 | 0.9923 / 0.9989 |
| predicted benefit mean / IQR / nontrivial | -5.861 / 15.940 / 0.9971 | -3.942 / 15.663 / 0.9942 |
| genomic split fraction obs / ACT | 0.3810 / 0.2761 | 0.3857 / 0.2785 |
| arm support obs (pt/ev) | 661.5 / 278.25 | 661.5 / 278.25 |
| arm support ACT (pt/ev) | 114.0 / 70.5 | 114.0 / 70.5 |
| raw propensity overlap / IPTW ESS | 0.8723 / 358.5 | 0.8694 / 366.3 |

Diagnostic score -9.768 (robust LCB -6.077, repeat range 3.6910).
Jaccard **0.1581 — the highest of the search.** Gates 11/13.

**Scorecard.**

- *Jaccard above 0.20*: **falsified**, but it did rise to 0.1581, the best
  eleven experiments have produced and about 1.95x the 0.0809 noise floor.
  Stability selection does what it claims; the ceiling on panel
  reproducibility in this cohort is simply lower than I hoped.
- *Both increments positive, range below 0.6355, reward above -3.9105*:
  **falsified.** Repeat 1 rose to +3.587 with a CI95 of [+0.483, +6.664] that
  excludes zero; repeat 2 fell to -0.104 and the range blew out to 3.691.
- *Reward positive*: **falsified.**

**The structural finding, which is now the main obstacle.** Within the
single-positive-module family — experiments 3, 5, 8, 9, 10 and 11, which
differ only in how the panel is picked — the two repeats respond in
*opposite* directions:

| | 003 | 005 | 008 | 009 | 010 | 011 |
|---|---|---|---|---|---|---|
| repeat 1 | +2.574 | +2.023 | +2.686 | +4.224 | +3.403 | +3.587 |
| repeat 2 | +1.939 | +1.145 | +1.401 | +0.902 | -0.327 | -0.104 |

Correlation across those six designs is **-0.535**. Every change that
sharpens the panel raises repeat 1 and lowers repeat 2. The two repeats use
the same 1,034 patients and differ only in the fold seed, so a genuine
biological effect modifier should help both partitions roughly equally. A
signal that improves one partition while degrading the other is
partition-specific — which is the third independent line of evidence, after
the detection-floor panel composition and the tail asymmetry of experiment 4,
that this module tracks a technical rather than biological axis.

The antagonism is also doubly punishing for the objective: pushing repeat 1
up pushes repeat 2 down, and the reward takes the *minimum* LCB and then
subtracts a range that the divergence itself inflates.

**Where that leaves the positive-reward target.** A positive reward needs the
binding repeat's LCB above zero, which at the observed bootstrap SD of about
1.7 means a mean increment near **+4.83**. The best repeat-2 increment in
eleven experiments is **+1.939**, and the antagonism means the designs that
reach +4 in repeat 1 are precisely those that drive repeat 2 negative.

One structural lever remains untried, and it is the only one that acts on the
LCB rather than on the mean. The increment's signal-to-noise scales roughly
as `sqrt(k) * E[tau | disagree] / sd(tau | disagree)` over the k patients
where the two policies disagree, because the mean grows linearly in k while
the bootstrap SD grows as its square root. Every candidate so far has *shrunk*
genomic influence, and therefore k. Experiment 12 goes the other way, on a
panel that is now demonstrably the most reproducible the search has produced.

## Experiment 12 (prespecified before running)

**Candidate** `tlearner_stability_two_module16`

Experiment 11's stability selector, unchanged — 40 subsampled repeats of the
pathway-anchored sign-coherent pick, 16 genes by selection frequency, both
arms at locked clinical geometry, threshold 0.0 — with one change:

- `module_count` 1 -> **2**, both modules from the benefit-**increasing**
  direction (ranks 1-8 and 9-16 by selection frequency; verified sign-pure
  positive on synthetic data before launch)

**This is not experiment 4.** That run's second module was the
benefit-decreasing tail, which the fold listings showed had no stable anchor
— mitochondrial, then ribosomal, then ER-glycosylation complexes in
successive folds. It was a noise module bolted onto a signal module. Here
both modules come from the same reproducible positive tail, and the panel
they are built from is the most stable the search has produced
(Jaccard 0.1581).

**Hypothesis, and why it targets the LCB rather than the mean.** Every
candidate since experiment 1 has *reduced* genomic influence, on the reading
that added columns cost more than they return. That reading is right about
the mean when the added column is noise, but it is the wrong objective. The
increment's signal-to-noise scales roughly as

    sqrt(k) * E[tau | disagree] / sd(tau | disagree)

over the k patients where the C+G and clinical policies disagree, because the
bootstrap mean grows linearly in k while the bootstrap SD grows as sqrt(k).
Doubling correctly-signed disagreement therefore buys about a factor of 1.41
in the ratio that the selection LCB actually depends on. Experiment 3 sits at
mean/SD = 2.543/1.74 = 1.46 in its better repeat and needs 2.84 for a
positive LCB. No lever tried so far moves that ratio structurally; this one
does, provided the extra influence is correctly signed — which is exactly
what experiment 4 failed to provide and experiment 11's panel now supplies.

**What I expect it to cost.** Genomic split fraction should rise from ~0.38
to roughly 0.5, as it did in experiment 4. The C-index drop has never
exceeded 0.0166 even at two columns, so the discrimination gate is not the
binding risk; the risk is that the second module adds influence without
adding independent information, since ranks 9-16 are correlated with ranks
1-8 by construction.

**Prespecified predictions and what would falsify them.**

- Both repeats' mean/SD ratios improve over experiment 11's, whatever happens
  to the means themselves. This is the direct test of the mechanism and I
  will compute it explicitly from the reported CI95 widths.
- Both increments positive, and the binding LCB above experiment 3's -3.275.
- If the ratios do *not* improve, then the extra column is adding influence
  that is correlated with what the first module already provides — k grows
  but `E[tau | disagree]` falls in proportion — and the disagreement-count
  lever is closed along with every other. At that point I will have tested
  every mechanism that acts on the mean, the variance, the panel, the
  representation and the geometry, and the honest conclusion is that a
  positive reward is unreachable for this estimand at n=1034 with 152 treated
  patients, with experiment 3 standing as the best eligible candidate.

**Result.** (to be appended after the run)

**Result — run_012_20260824T010720Z, 324 s wall. NOT ELIGIBLE, reward -1e6.**

| field | repeat 1 | repeat 2 |
|---|---|---|
| incremental alignment (months) | **+5.132** | **-0.268** |
| selection LCB | **-0.464** | -7.109 |
| bootstrap mean / CI95 | +5.109 / **[+1.488, +8.936]** | -0.314 / [-4.646, +3.548] |
| alignment C / C+G | 5.098 / **10.230** | 5.361 / 5.093 |
| value C / C+G / best constant | 48.152 / 48.933 / 45.468 | 48.526 / 47.685 / 45.502 |
| Harrell C, C / C+G (drop) | 0.6561 / 0.6588 (**-0.0027**) | 0.6641 / 0.6486 (0.0155) |
| ACT recommended fraction C+G | 0.2911 | 0.3366 |
| seed agreement / benefit corr | 0.9916 / 0.9990 | 0.9945 / 0.9988 |
| predicted benefit mean / IQR / nontrivial | -6.158 / 15.581 / 0.9923 | -4.313 / 15.371 / 0.9990 |
| genomic split fraction obs / ACT | 0.4926 / 0.3585 | 0.5090 / 0.3633 |
| arm support obs (pt/ev) | 661.5 / 278.25 | 661.5 / 278.25 |
| arm support ACT (pt/ev) | 114.0 / 70.5 | 114.0 / 70.5 |
| raw propensity overlap / IPTW ESS | 0.8723 / 358.5 | 0.8694 / 366.3 |

Diagnostic score -12.509 (robust LCB -7.109, repeat range **5.4004**).
Gates 11/13.

**The prespecified mechanism check, computed from the CI95 widths.**

| candidate | repeat | mean | SD | **mean/SD** | LCB |
|---|---|---|---|---|---|
| 003 bounded pool | 1 | +2.543 | 1.737 | 1.464 | -2.406 |
| 003 bounded pool | 2 | +1.904 | 1.735 | 1.097 | -3.275 |
| 011 stability | 1 | +3.558 | 1.577 | 2.257 | -1.140 |
| 011 stability | 2 | -0.161 | 1.942 | -0.083 | -6.077 |
| 012 stability + 2 modules | 1 | **+5.109** | 1.900 | **2.689** | **-0.464** |
| 012 stability + 2 modules | 2 | -0.313 | 2.090 | -0.150 | -7.109 |

**The mechanism is real and it worked.** Repeat 1's ratio went 1.464 -> 2.257
-> 2.689 across exactly the sequence of changes predicted to raise it, and
2.689 is within striking distance of the 2.84 a positive LCB requires. Repeat
1 produced the best numbers of the entire search by a wide margin: increment
+5.132, alignment 10.230 against the clinical 5.098, a C-index that
*improved* by 0.0027, and a bootstrap CI95 of [+1.488, +8.936]. More
correctly-signed genomic influence does raise the quantity the selection LCB
depends on, as predicted.

**And it changes nothing, because the antagonism is stronger.** The same
change drove repeat 2 from -0.083 to -0.150 and the repeat range to 5.4004.
Prediction 1 asked for *both* ratios to improve; one did and one did not, so
the prediction is falsified as stated.

**The decisive measurement.** The repeat range is now **5.40 months while the
bootstrap SD is 1.9-2.1**. The two repeats estimate the same population
quantity on the same 1,034 patients and differ only in the fold seed, so the
partition contributes roughly two and a half times more uncertainty than
patient resampling does. **The dominant source of error in this estimand is
cross-fitting variability, which the bootstrap does not see and the reward's
range term exists to charge for.** That reframes the whole optimisation: the
binding constraint was never the mean, the panel, or the variance the
bootstrap measures. It is that a policy improvement demonstrable in one
partition of these patients is not demonstrable in another.

**Consequence for the target, stated plainly.** Maximising repeat 1 is worth
nothing — the reward takes the *minimum* LCB and then subtracts a range that
strengthening inflates. Every remaining slot should go to raising the
*worst* repeat, which means trading panel strength for panel constancy across
folds. Experiment 13 does that: it stops selecting genes by rank, which is the
statistic that has proven unstable at every level of this search, and selects
both the anchor and its members by *frequency* statistics instead.

## Experiment 13 (prespecified before running)

**Candidate** `tlearner_modal_anchor_module16`

- `module_count` back to **1** (experiment 12 showed a second column inflates
  the repeat range, 3.69 -> 5.40, and the reward takes the *minimum* LCB)
- `POOL_PATHWAYS` **2**, `PER_PATHWAY_TOP` 8, `MIN_PATHWAY_MEMBERS` 25
- both arms at locked clinical geometry, threshold 0.0
- **rank is removed from both selection decisions:**
  1. anchors are the pathways appearing most *often* in the top 2 across the
     40 subsamples, not the ones ranking highest in any single fit
  2. within an anchor, members are ordered by how often their signed score is
     *positive* across those subsamples — sign consistency is bounded in
     [0, 1] and far more stable than score magnitude, while remaining a
     benefit-relevant criterion rather than an arbitrary one

**Hypothesis.** Experiment 12 identified the real binding constraint: the
repeat range is 5.40 months against a bootstrap SD of 1.9-2.1, so the fold
partition contributes about two and a half times more uncertainty than
patient resampling. The reward takes the minimum repeat LCB and subtracts
that range, so the only move with any value left is raising the *worst*
repeat by making the panel more nearly constant across folds. Rank is the
statistic that has been unstable at every level of this search — across folds
(experiment 1, Jaccard 0.070), across biologically structured pools
(experiment 6), and under 2-5 gene perturbations (experiment 10, range swing
of 3.1 months). Replacing it with bounded frequency counts at both decision
points should raise cross-fold panel constancy and pull the two repeats
together.

**A calibration that makes this a harder test than it looks.** Measured on
synthetic pure noise at full scale, cross-fold panel Jaccard over four
overlapping 3/4 partitions is **0.0330** for experiment 11's selector and
**0.0111** for this one. Pooling five anchors manufactures agreement through
the pool restriction alone; two modal anchors drawn from 1,839 candidates
almost never coincide by chance. So this design starts from a noise floor
three times lower and must clear the same 0.10 gate on real signal. If it
does clear it, the Jaccard is far more attributable to genuine pathway
reproducibility than any previous run's — and if it fails the gate, that is
itself the cleanest evidence yet about how little reproducible pathway-level
benefit signal this cohort contains.

**Prespecified predictions and what would falsify them.**

- Repeat range below experiment 3's 0.6355, and the *minimum* repeat
  increment above +1.939, which is the best any repeat-2 has achieved.
- Jaccard clears 0.10 despite the threefold lower noise floor.
- Reward above -3.9105, and eligible.
- If the range does not fall, then cross-fold panel constancy is not what
  drives the partition sensitivity either — the two partitions would be
  disagreeing about the *effect* of a stable panel rather than about which
  panel to use. That would mean the reproducible component of this signal is
  smaller than the cross-fitting noise at n=1034 with 152 treated patients,
  the positive-reward target is unreachable on this estimand, and I will
  report that conclusion with experiment 3 as the standing eligible
  candidate rather than spend the remaining slots re-testing it.

**Result.** (to be appended after the run)

**Result — run_013_20260824T011659Z, 260 s wall. NOT ELIGIBLE, reward -1e6.**

| field | repeat 1 | repeat 2 |
|---|---|---|
| incremental alignment (months) | **+0.749** | **+0.568** |
| selection LCB | -3.889 | -3.650 |
| bootstrap mean / CI95 | +0.746 / [-2.166, +3.479] | +0.573 / [-2.468, +3.205] |
| alignment C / C+G | 5.098 / 5.847 | 5.361 / 5.929 |
| value C / C+G / best constant | 48.152 / 47.496 / 45.392 | 48.526 / 48.096 / 45.477 |
| Harrell C, C / C+G (drop) | 0.6561 / 0.6503 (0.0058) | 0.6641 / 0.6547 (0.0094) |
| ACT recommended fraction C+G | 0.3008 | 0.3598 |
| seed agreement / benefit corr | 0.9968 / 0.9991 | 0.9952 / 0.9991 |
| predicted benefit mean / IQR / nontrivial | -4.856 / 15.498 / 0.9990 | -3.864 / 16.176 / 0.9952 |
| genomic split fraction obs / ACT | 0.3774 / 0.2296 | 0.3702 / 0.2878 |
| arm support obs (pt/ev) | 661.5 / 278.25 | 661.5 / 278.25 |
| arm support ACT (pt/ev) | 114.0 / 70.5 | 114.0 / 70.5 |
| raw propensity overlap / IPTW ESS | 0.8723 / 358.5 | 0.8694 / 366.3 |

Diagnostic score -4.070 (robust LCB -3.889, repeat range **0.1806**).
Jaccard **0.1089**. Gates 12/13 — `genomic_value_at_least_clinical`, and this
time in *both* repeats (-0.657 and -0.430).

**Scorecard.**

- *Repeat range below 0.6355*: **confirmed, decisively.** 0.1806, less than a
  third of the previous best (experiment 7's 0.5341) and thirty times smaller
  than experiment 12's. The two increments came in at +0.749 and +0.568 —
  the two repeats finally agree with each other.
- *Jaccard clears 0.10 from a threefold lower noise floor*: **confirmed**,
  0.1089 against a measured pure-noise floor of 0.0111. This is the most
  honestly-earned Jaccard of the search: roughly ten times its noise floor,
  where experiment 3's 0.1417 was 1.75x its own.
- *Minimum repeat increment above +1.939*: **falsified**, +0.568.
- *Reward above -3.9105 and eligible*: **falsified**, -4.070 and not eligible.

**The frontier is now mapped, and it is the answer.** Three designs bracket
it:

| | repeat range | increments | min mean | diagnostic |
|---|---|---|---|---|
| 012 max strength | 5.4004 | +5.132 / -0.268 | -0.313 | -12.509 |
| **003 balanced** | 0.6355 | +2.574 / +1.939 | +1.904 | **-3.911** |
| 013 max constancy | **0.1806** | +0.749 / +0.568 | +0.573 | -4.070 |

Panel constancy controls the repeat range almost perfectly, and it costs
increment magnitude almost exactly in proportion. Experiment 3 sits near the
interior optimum of that trade-off by accident rather than design.

**Why the trade-off is itself the diagnosis, and what it implies for the
positive-reward target.** If a large, genuinely reproducible benefit signal
existed in this cohort, strength and constancy would *not* trade off: every
fold would find the same strong genes, and an adaptive selector would be both
powerful and stable. Observing a clean inverse relationship between them
means the strength on display in experiment 12's repeat 1 — increment +5.132,
alignment 10.230, a CI95 excluding zero — is fitted to that partition's
particular patients. It is not a discovery that failed to replicate; it is
the thing that replication is designed to catch.

A positive reward requires a minimum-repeat mean near **+4.83** at the
observed SDs. The frontier says the best achievable minimum mean is about
**+1.9**, and that pushing toward +4.8 in one repeat drives the other to zero
or below while inflating the range that gets subtracted. The two requirements
are in direct opposition on this estimand.

One point on the frontier remains untested and is worth a slot: experiment 13
achieved its constancy by removing rank from *both* the anchor choice and the
within-anchor ordering. The anchor choice is what produces the large
fold-level jumps; the within-anchor ordering is a much smaller perturbation.
Restoring signal strength to the ordering while keeping the anchor modal
should land between experiments 3 and 13 — and experiment 13's bootstrap SD
was also the lowest seen (about 1.45 against experiment 3's 1.74), so a
design with experiment 3's mean and experiment 13's SD and range would score
near -1.8 rather than -3.9.

## Experiment 14 (prespecified before running)

**Candidate** `tlearner_modal_anchor_scored_module16`

Experiment 13's design — modal anchors chosen by appearance frequency across
40 subsamples, 2 anchors of 8 genes, one module, locked clinical geometry,
threshold 0.0 — with exactly one change:

- within an anchor, members are ordered by **mean signed score** rather than
  by sign-consistency count

**Hypothesis.** Experiment 13 removed rank from *both* selection decisions
and landed at the constancy end of the frontier: repeat range 0.1806 but a
minimum increment of only +0.568. The two decisions are not equally
responsible for instability. The *anchor* choice is what produces the large
fold-level jumps — an entire panel switching from olfactory receptors to
histones — while the *within-anchor* ordering only reshuffles which members
of the same pathway are averaged, a far smaller perturbation of the module.
Keeping the anchor modal while restoring score-based ordering should recover
increment magnitude at a small cost in range.

**The arithmetic that makes this worth a slot.** Experiment 13 also produced
the lowest bootstrap SD of the search, about **1.45** against experiment 3's
1.74, because a more constant panel yields a more stable per-patient
contribution. Reward is roughly `min_mean - 2.84*SD - range`. Experiment 3
scores -3.911 from (1.904, 1.74, 0.636); experiment 13 scores -4.070 from
(0.573, 1.45, 0.181). A design landing at experiment 3's minimum mean with
experiment 13's SD and range would score about **-1.8**, a full two months
better than the standing best, without needing any new signal — purely from
occupying a better point on the frontier already mapped.

**On the positive-reward target, stated honestly before the run.** This will
not reach a positive reward and I do not want the prediction record to
pretend otherwise. Positive requires a minimum-repeat mean near +4.83; the
frontier's best minimum mean is +1.904, and experiments 9 through 13
established that pushing one repeat toward +5 drives the other to zero while
inflating the subtracted range. What this run can do is improve the standing
eligible reward, and demonstrate that the improvement comes from the
variance and range terms rather than from finding more signal.

**Prespecified predictions and what would falsify them.**

- Both increments positive and the minimum above experiment 13's +0.568,
  with the repeat range staying below experiment 3's 0.6355.
- Eligible, with a reward better than -3.9105.
- If the minimum increment does not rise above +0.568, then the within-anchor
  ordering carries no usable signal at all and essentially everything this
  selector family achieves comes from the anchor choice alone. Combined with
  the frontier already mapped, that would close the last untested point on it
  and I will report the positive-reward target as unreachable for this
  estimand, with the best eligible candidate standing.

**Result.** (to be appended after the run)

**Result — run_014_20260824T012332Z, 260 s wall. NOT ELIGIBLE, reward -1e6.**

| field | repeat 1 | repeat 2 |
|---|---|---|
| incremental alignment (months) | **+1.203** | **-0.644** |
| selection LCB | -3.033 | -5.816 |
| bootstrap mean / CI95 | +1.202 / [-1.639, +3.891] | -0.680 / [-4.270, +2.314] |
| alignment C / C+G | 5.098 / 6.301 | 5.361 / 4.717 |
| value C / C+G / best constant | 48.152 / 47.647 / 45.439 | 48.526 / 47.595 / 45.466 |
| Harrell C, C / C+G (drop) | 0.6561 / 0.6508 (0.0053) | 0.6641 / 0.6547 (0.0094) |
| ACT recommended fraction C+G | 0.3056 | 0.3627 |
| seed agreement / benefit corr | 0.9958 / 0.9991 | 0.9958 / 0.9991 |
| predicted benefit mean / IQR / nontrivial | -4.876 / 15.578 / 0.9981 | -3.606 / 16.263 / 0.9971 |
| genomic split fraction obs / ACT | 0.3754 / 0.2314 | 0.3750 / 0.2799 |
| arm support obs (pt/ev) | 661.5 / 278.25 | 661.5 / 278.25 |
| arm support ACT (pt/ev) | 114.0 / 70.5 | 114.0 / 70.5 |
| raw propensity overlap / IPTW ESS | 0.8723 / 358.5 | 0.8694 / 366.3 |

Diagnostic score -7.663 (robust LCB -5.816, repeat range 1.8463). Gates 11/13.

**Scorecard: falsified, and the manner of failure is the finding.** Restoring
score-based ordering *within* a modal anchor raised repeat 1 from +0.749 to
+1.203 and drove repeat 2 from +0.568 to **-0.644**, with the range going
0.1806 -> 1.8463. The minimum increment fell rather than rose, so the
prespecified rule fires.

I had argued that within-anchor ordering would be "a far smaller
perturbation" than anchor choice. It is smaller in magnitude — the swing here
is about 1.2 months where anchor-level changes swung 5 — but it is **the same
antagonism, in the same direction, at a smaller scale**. Repeat 1 up, repeat
2 down, range up. The strength/constancy trade-off is not a property of one
particular decision in the selector; it holds at every level of the selection
hierarchy that has been tested: which genes (experiment 1), which pool
(experiments 2-3), how many anchors (experiments 5, 9), which anchors
(experiment 10), how many columns (experiments 4, 12), and now which members
within a fixed anchor. **There is no level at which strength is free.**

## Conclusion: no eligible candidate with a positive reward exists, and the target is unreachable on this estimand

**The identification result.** Across 14 experiments, exactly one run is
eligible — `run_003_20260823T191920Z`, 13/13 gates, reward **-3.9105**. No
run is both eligible and positively scored, and none is close: the best
diagnostic score among non-eligible runs is -4.070.

**Why this is a property of the estimand rather than of the search.** Write
`f` for the fraction of patients where the C+G and clinical policies
disagree, and `tau` for the per-patient AIPW benefit score. Only disagreeing
patients contribute to the increment, so

    increment mean  ~  2 * f * E[tau | disagree]
    bootstrap SD    ~  2 * sqrt(f) * sqrt(E[tau^2 | disagree]) / sqrt(n)
    =>  mean / SD   ~  sqrt(f * n) * (mu_tau / sigma_tau)

A positive selection LCB needs `mean/SD > 2.84`, the 0.25% normal quantile
that the alpha/20 multiplicity correction imposes. Calibrating from
experiment 3 (mean 2.543, SD 1.737, n 1034) gives `mu_tau / sigma_tau ~ 0.118`
— the AIPW scores carry a standard deviation near 70 months because the
propensity clip admits weights up to 20, against a mean signed benefit of
about 8.5 months on the disagreement set. The requirement becomes

    sqrt(f * 1034) * 0.118 > 2.84   =>   f > 0.56

Even at **f = 1**, with every patient's recommendation flipped relative to
the clinical rule, the requirement is `mu_tau > 6.4` months of correctly
signed benefit *averaged over the whole cohort* — against a cohort ATE of
about **-0.2 months**. Since `n`, the propensity clipping, the IPCW estimand
and the 0.25% quantile are all locked, no edit to `train.py` can reach it.
The negative reward is imposed by the multiplicity correction acting on a
weak effect in a small treated arm, not by a failure to find the right panel.

**And the empirical frontier says the same thing from the other side.**

| | range | increments | min mean | mean/SD (best repeat) | diagnostic |
|---|---|---|---|---|---|
| 012 max strength | 5.4004 | +5.132 / -0.268 | -0.313 | 2.689 | -12.509 |
| 009 | 3.3216 | +4.224 / +0.902 | +0.902 | 2.228 | -8.856 |
| **003 balanced** | 0.6355 | +2.574 / +1.939 | **+1.904** | 1.464 | **-3.911** |
| 014 | 1.8463 | +1.203 / -0.644 | -0.680 | 0.847 | -7.663 |
| 013 max constancy | **0.1806** | +0.749 / +0.568 | +0.573 | 0.516 | -4.070 |

Every design that raises one repeat lowers the other; the best attainable
minimum increment is **+1.904**, against the **+4.83** a positive LCB needs.
Strength and constancy trade off cleanly, which is precisely what one
observes when the apparent signal is fitted to a partition rather than
present in the population — a genuinely reproducible effect would let a
selector be powerful and stable at once.

**Standing answer to the question asked.** Eligible candidates: one
(`run_003_20260823T191920Z`). Eligible candidates with positive score: none,
and none is attainable within the arena's fixed budget, estimand and
multiplicity correction. Six experiments remain unspent; I am not consuming
them, because every mechanism that acts on the mean, the variance, the panel,
the representation, the geometry and the disagreement fraction has now been
tested and the remaining gap is a factor of 2.5 in a quantity the editable
surface does not control. Test data remains untouched.
