# HTE screen with `grf`

`grf_hte_test.py` asks one question: **is there heterogeneous treatment effect
(HTE) signal for adjuvant chemotherapy in this cohort?** It drives R's
[`grf`](https://grf-labs.github.io/grf/) through `rpy2` and runs the standard
screen — causal forest, overlap check, calibration test, RATE/AUTOC.

## Running it

```bash
/Library/Frameworks/Python.framework/Versions/3.11/bin/python3 test_for_hte/grf_hte_test.py --features genomic
```

Useful variations:

```bash
# clinical covariates only (p = 18)
... --features clinical

# genes + clinical covariates
... --features both

# the classic causal_forest recipe instead of the survival forest
... --forest causal

# cheap smoke run
... --num-trees 500 --tune none --rate-repeats 0
```

A full genomic run (n = 1034, p = 13055, 4000 trees, `tune.parameters="all"`,
3 split-sample RATE repeats) takes about 90 seconds on this machine.

Each run writes to `test_for_hte/results/<features>_<forest>_<splits>/`:

| file | contents |
| --- | --- |
| `report.txt` | the printed report |
| `summary.json` | every number, machine-readable |
| `oob_predictions.csv` | per-patient `W.hat` and out-of-bag `tau_hat` |
| `propensity_overlap.png` | `hist(cf$W.hat)`, overall and by arm |
| `toc_in_sample.png` | TOC curve, priorities = the forest's own OOB CATE |
| `toc_split_sample.png` | TOC curve, priorities from a forest fit on the other half |

`affyfRMATest.csv` is never opened — `--splits` only accepts `train` and
`validation` (see `red_lines.md`).

## What the script does, and why it deviates from the plain recipe

**Censoring.** The outcome is censored overall survival, so the default forest
is `causal_survival_forest(X, Y, W, D, target = "RMST", horizon = 60)` rather
than `causal_forest(X, Y, W)`; 60 months matches the repository's tau. Passing
raw `OS_MONTHS` as `Y` to `causal_forest` would treat every censoring time as a
death. `--forest causal` still runs the plain causal forest, but by default on
an IPCW restricted-mean pseudo-outcome (Kaplan–Meier censoring weights,
patients censored before the horizon dropped). `--causal-outcome raw-time`
gives the literal `causal_forest(X, Y, W)` recipe and is descriptive only.

**Calibration test.** `grf::test_calibration` supports `causal_forest` but
raises *"Calibration check not supported for this type of forest"* for
`causal_survival_forest`. For the survival forest the script computes the same
best linear projection on the forest's doubly-robust scores
(`grf::get_scores`), with the same `sandwich` HC3 variance and the same
one-sided p-value convention grf uses, so the two output tables read
identically. `differential.forest.prediction` is the omnibus HTE test.

**Overlap.** This is observational data with ~15% treated, so the report leads
with the `W.hat` distribution and flags propensity mass near 0 or 1 before the
calibration result. `--propensity-features clinical` fits `W.hat` from a
regression forest on the clinical covariates only, if a 13k-gene propensity
model is not something you want to lean on.

**RATE.** `rank_average_treatment_effect(cf, priorities = predict(cf)$predictions)`
ranks and evaluates on the same rows, which is optimistic. The script reports
that number (AUTOC and QINI) *and* an honest version: fit a forest on one
stratified half, use its predictions as priorities for a forest fit on the
other half, both directions, `--rate-repeats` times. Read the split-sample
column.

## Environment notes

- **R**: the script embeds R **4.5.1** with `grf` 2.6.1, not the default R 4.6.1
  on this machine. Embedding R 4.6.1 through `rpy2` segfaults inside
  `setup_Rmainloop`: the embedded interpreter sees an empty `environ`, so
  `base::Sys.getenv()` fails with *"invalid substring arguments"* and R aborts.
  `Rscript` 4.6.1 is unaffected — this only bites the embedded path.
  `grf_hte_test.py` probes candidate `R_HOME`s in a child process and picks the
  first that can actually be embedded with `grf` loadable, so it recovers on
  its own; `--r-home` overrides the choice. `grf` was installed into
  `~/Library/R/arm64/4.5/library` for this purpose.
- **Python**: needs `numpy`, `pandas`, and `rpy2` (3.6.7 here). Of the
  interpreters on this machine only
  `/Library/Frameworks/Python.framework/Versions/3.11/bin/python3` has all
  three; the conda base env has the survival stack but no `rpy2`. The script
  deliberately does not import `clinical_data.py`, so it needs no `lifelines`
  or `scikit-survival`.
