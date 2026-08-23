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
