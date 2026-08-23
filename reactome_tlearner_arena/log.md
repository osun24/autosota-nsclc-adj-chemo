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
