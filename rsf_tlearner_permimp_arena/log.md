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
