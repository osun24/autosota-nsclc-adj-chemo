"""Test whether the NSCLC adjuvant-chemotherapy data carry heterogeneous
treatment effect (HTE) signal, using R's ``grf`` through ``rpy2``.

The script follows the standard grf HTE-screening recipe:

1. Fit a causal forest (nuisances estimated internally, honest, out-of-bag).
2. Check overlap on the estimated propensity scores ``W.hat``.  The cohort is
   observational and heavily imbalanced (about 15% treated), so a calibration
   test run without an overlap check is not interpretable.
3. Run the calibration test on out-of-bag predictions with HC3 SEs.  The
   ``differential.forest.prediction`` coefficient is the omnibus HTE test.
4. Cross-check with RATE/AUTOC, both in-sample (as commonly written, and
   optimistically biased) and with honest train/evaluate splits.

Outcome handling
----------------
The outcome here is censored overall survival, not a plain numeric ``Y``, so
``causal_survival_forest`` (target ``RMST``, horizon 60 months to match the
repository's tau) is the default and the statistically correct choice.
``--forest causal`` runs the plain ``causal_forest`` from the classic recipe on
an IPCW restricted-mean pseudo-outcome (or, with
``--causal-outcome raw-time``, on the raw follow-up time, which ignores
censoring and is reported for reference only).

``grf::test_calibration`` supports ``causal_forest`` but not
``causal_survival_forest``; for the survival forest the identical best-linear-
projection test is computed here on the forest's doubly-robust scores, with the
same ``sandwich`` HC3 variance and the same one-sided p-value convention grf
uses.

Sealed test set
---------------
``affyfRMATest.csv`` is never read (see ``red_lines.md``); only the train and
validation splits are available to this script.

Usage
-----
    python test_for_hte/grf_hte_test.py --features genomic
    python test_for_hte/grf_hte_test.py --features clinical --tune all
    python test_for_hte/grf_hte_test.py --features both --top-var-genes 2000

Run it with an interpreter that has ``rpy2`` installed; see
``test_for_hte/README.md``.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent

GENOMIC_CSV = {
    "train": REPO_ROOT / "affyfRMATrain.csv",
    "validation": REPO_ROOT / "affyfRMAValidation.csv",
}
CLINICAL_CSV = {
    "train": REPO_ROOT / "clinicalTrain.csv",
    "validation": REPO_ROOT / "clinicalValidation.csv",
}

TREATMENT = "Adjuvant Chemo"
TIME = "OS_MONTHS"
EVENT = "OS_STATUS"

# Mirrors clinical_data.CLINICAL_VARS minus the treatment column.  Duplicated
# rather than imported so this script only needs numpy/pandas/rpy2.
CLINICAL_COVARIATES = [
    "Age",
    "IS_MALE",
    "Stage_IA",
    "Stage_IB",
    "Stage_II",
    "Stage_III",
    "Histology_Adenocarcinoma",
    "Histology_Adenosquamous Carcinoma",
    "Histology_Large Cell Carcinoma",
    "Histology_Squamous Cell Carcinoma",
    "Race_African American",
    "Race_Asian",
    "Race_Caucasian",
    "Race_Native Hawaiian or Other Pacific Islander",
    "Race_Unknown",
    "Smoked?_No",
    "Smoked?_Unknown",
    "Smoked?_Yes",
]
NON_GENE_COLUMNS = {"Unnamed: 0", TREATMENT, "Age", EVENT, "IS_MALE", TIME}


# --------------------------------------------------------------------------
# R environment bootstrap
# --------------------------------------------------------------------------
def _r_libpaths(r_home: str) -> list[str]:
    """Ask an R installation for its own library search path."""
    for exe in (Path(r_home) / "bin" / "exec" / "R", Path(r_home) / "bin" / "R"):
        if not exe.exists():
            continue
        env = dict(os.environ, R_HOME=r_home)
        env.pop("R_LIBS", None)
        env.pop("R_LIBS_USER", None)
        try:
            out = subprocess.run(
                [str(exe), "--vanilla", "--no-echo", "-e", 'cat(.libPaths(), sep="\\n")'],
                capture_output=True, text=True, env=env, timeout=120,
            )
        except (OSError, subprocess.TimeoutExpired):
            continue
        paths = [line.strip() for line in out.stdout.splitlines() if line.strip()]
        if paths:
            return paths
    return []


def _candidate_r_homes(explicit: str | None) -> list[str]:
    candidates: list[str] = []
    for value in (explicit, os.environ.get("HTE_R_HOME"), os.environ.get("R_HOME")):
        if value:
            candidates.append(value)
    framework = Path("/Library/Frameworks/R.framework/Versions")
    if framework.is_dir():
        # Newest first, but "Current" last: the newest R is not necessarily the
        # one that can be embedded (see README on the R 4.6.1 crash).
        versions = sorted(
            (p for p in framework.iterdir() if p.is_dir() and p.name != "Current"),
            key=lambda p: p.name,
            reverse=True,
        )
        candidates.extend(str(p / "Resources") for p in versions)
        candidates.append(str(framework / "Current" / "Resources"))
    rbin = shutil.which("R")
    if rbin:
        try:
            out = subprocess.run([rbin, "RHOME"], capture_output=True, text=True, timeout=120)
            if out.returncode == 0 and out.stdout.strip():
                candidates.append(out.stdout.strip())
        except (OSError, subprocess.TimeoutExpired):
            pass
    seen: set[str] = set()
    ordered = []
    for candidate in candidates:
        resolved = str(Path(candidate))
        if resolved not in seen and Path(resolved).is_dir():
            seen.add(resolved)
            ordered.append(resolved)
    return ordered


_PROBE = (
    "import rpy2.robjects as ro;"
    "ro.r('suppressMessages(library(grf))');"
    "print('RPROBE|' + ro.r('R.version.string')[0] +"
    "'|' + ro.r('as.character(packageVersion(\"grf\"))')[0])"
)


def _probe_r_home(r_home: str) -> tuple[bool, str, dict[str, str]]:
    """Return (works, message, env) for embedding this R with grf loadable.

    Embedding R is done in a child process because a mismatched R crashes the
    interpreter outright (R 4.6.1 on macOS segfaults inside setup_Rmainloop).
    """
    libpaths = _r_libpaths(r_home)
    if not libpaths:
        return False, "could not query .libPaths()", {}
    env_overrides = {"R_HOME": r_home, "R_LIBS": os.pathsep.join(libpaths)}
    env = dict(os.environ, **env_overrides)
    env.pop("R_LIBS_USER", None)
    try:
        out = subprocess.run(
            [sys.executable, "-c", _PROBE],
            capture_output=True, text=True, env=env, timeout=300,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, f"probe failed: {exc}", {}
    line = next((ln for ln in out.stdout.splitlines() if ln.startswith("RPROBE|")), None)
    if out.returncode != 0 or line is None:
        detail = (out.stderr or out.stdout).strip().splitlines()
        reason = detail[-1] if detail else f"exit {out.returncode}"
        if out.returncode < 0:
            reason = f"crashed with signal {-out.returncode}; {reason}"
        return False, reason, {}
    _, r_version, grf_version = line.split("|")
    env_overrides["_R_VERSION"] = r_version
    env_overrides["_GRF_VERSION"] = grf_version
    return True, f"{r_version}, grf {grf_version}", env_overrides


def bootstrap_r(explicit_r_home: str | None, verbose: bool = True) -> dict[str, str]:
    """Point rpy2 at an R installation that can actually be embedded."""
    tried: list[str] = []
    for r_home in _candidate_r_homes(explicit_r_home):
        works, message, env_overrides = _probe_r_home(r_home)
        if works:
            if verbose:
                print(f"[R] using R_HOME={r_home} ({message})")
            os.environ["R_HOME"] = env_overrides["R_HOME"]
            os.environ["R_LIBS"] = env_overrides["R_LIBS"]
            os.environ.pop("R_LIBS_USER", None)
            return env_overrides
        tried.append(f"  {r_home}: {message}")
        if explicit_r_home and r_home == str(Path(explicit_r_home)):
            break
    raise SystemExit(
        "No usable R installation found (needs an embeddable R with the grf "
        "package installed).\nTried:\n" + "\n".join(tried) +
        "\n\nInstall grf with:\n"
        "  R -e 'install.packages(\"grf\", repos=\"https://cloud.r-project.org\")'\n"
        "or point at a specific R with --r-home."
    )


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------
def _assert_not_test_path(path: Path) -> None:
    assert "test" not in path.name.lower(), f"Test data is sealed: {path}"


def _map_treatment(series: pd.Series) -> np.ndarray:
    mapped = series.map({"OBS": 0, "ACT": 1})
    mapped = mapped.where(mapped.notna(), series)
    return pd.to_numeric(mapped, errors="raise").astype(int).to_numpy()


@dataclass
class Cohort:
    X: np.ndarray
    feature_names: list[str]
    Y: np.ndarray          # follow-up time in months
    D: np.ndarray          # 1 = death observed, 0 = censored
    W: np.ndarray          # 1 = adjuvant chemotherapy, 0 = observation
    split: np.ndarray      # split label per row
    clinical_index: list[int] = field(default_factory=list)

    def subset(self, keep: np.ndarray) -> "Cohort":
        return Cohort(
            X=np.ascontiguousarray(self.X[keep]),
            feature_names=self.feature_names,
            Y=self.Y[keep], D=self.D[keep], W=self.W[keep], split=self.split[keep],
            clinical_index=self.clinical_index,
        )


def load_cohort(features: str, splits: list[str], top_var_genes: int) -> Cohort:
    frames_clin: list[pd.DataFrame] = []
    frames_gene: list[pd.DataFrame] = []
    split_labels: list[np.ndarray] = []

    for split in splits:
        clinical_path = CLINICAL_CSV[split]
        _assert_not_test_path(clinical_path)
        clinical = pd.read_csv(clinical_path)
        frames_clin.append(clinical)
        split_labels.append(np.full(len(clinical), split, dtype=object))

        if features in {"genomic", "both"}:
            genomic_path = GENOMIC_CSV[split]
            _assert_not_test_path(genomic_path)
            print(f"[data] reading {genomic_path.name} ...", flush=True)
            genomic = pd.read_csv(genomic_path)
            if len(genomic) != len(clinical):
                raise ValueError(
                    f"{genomic_path.name} has {len(genomic)} rows but "
                    f"{clinical_path.name} has {len(clinical)}"
                )
            # The two files are row-aligned; verify before relying on it.
            aligned = np.allclose(
                genomic[TIME].to_numpy(float), clinical[TIME].to_numpy(float)
            ) and np.array_equal(
                _map_treatment(genomic[TREATMENT]), _map_treatment(clinical[TREATMENT])
            )
            if not aligned:
                raise ValueError(
                    f"{genomic_path.name} and {clinical_path.name} are not row-aligned"
                )
            gene_cols = [c for c in genomic.columns if c not in NON_GENE_COLUMNS]
            frames_gene.append(genomic[gene_cols])

    clinical = pd.concat(frames_clin, ignore_index=True)
    split = np.concatenate(split_labels)

    Y = pd.to_numeric(clinical[TIME], errors="raise").to_numpy(float)
    D = pd.to_numeric(clinical[EVENT], errors="raise").to_numpy(int)
    W = _map_treatment(clinical[TREATMENT])

    blocks: list[np.ndarray] = []
    names: list[str] = []
    clinical_index: list[int] = []
    if features in {"clinical", "both"}:
        missing = [c for c in CLINICAL_COVARIATES if c not in clinical.columns]
        if missing:
            raise ValueError(f"missing clinical covariates: {missing}")
        block = clinical[CLINICAL_COVARIATES].apply(pd.to_numeric).to_numpy(float)
        clinical_index = list(range(block.shape[1]))
        blocks.append(block)
        names.extend(CLINICAL_COVARIATES)
    if features in {"genomic", "both"}:
        genes = pd.concat(frames_gene, ignore_index=True)
        gene_matrix = genes.to_numpy(dtype=float)
        gene_names = list(genes.columns)
        keep = gene_matrix.std(axis=0) > 0
        if not keep.all():
            print(f"[data] dropping {int((~keep).sum())} zero-variance genes")
            gene_matrix = gene_matrix[:, keep]
            gene_names = [n for n, k in zip(gene_names, keep) if k]
        if top_var_genes and top_var_genes < gene_matrix.shape[1]:
            # Unsupervised filter: uses neither the outcome nor the treatment.
            order = np.argsort(gene_matrix.var(axis=0))[::-1][:top_var_genes]
            order = np.sort(order)
            gene_matrix = gene_matrix[:, order]
            gene_names = [gene_names[i] for i in order]
            print(f"[data] keeping the {top_var_genes} highest-variance genes")
        blocks.append(gene_matrix)
        names.extend(gene_names)

    X = np.ascontiguousarray(np.column_stack(blocks), dtype=float)

    usable = np.isfinite(Y) & (Y > 0) & np.isfinite(X).all(axis=1)
    if not usable.all():
        print(f"[data] dropping {int((~usable).sum())} rows with non-positive or missing values")
        X, Y, D, W, split = X[usable], Y[usable], D[usable], W[usable], split[usable]

    return Cohort(X=X, feature_names=names, Y=Y, D=D, W=W, split=split,
                  clinical_index=clinical_index)


def ipcw_rmst_outcome(Y: np.ndarray, D: np.ndarray, horizon: float):
    """Restricted-mean pseudo-outcome with Kaplan-Meier censoring weights.

    Returns (Y_rmst, weights, kept).  Patients censored before the horizon get
    weight zero; everyone else is up-weighted by 1 / G(t-), the KM estimate of
    the censoring survival function.
    """
    order = np.argsort(Y, kind="mergesort")
    times = Y[order]
    censored = 1 - D[order]
    n = len(Y)
    at_risk = n - np.arange(n)
    surv = 1.0
    grid_t: list[float] = []
    grid_g: list[float] = []
    for i in range(n):
        if censored[i] == 1 and at_risk[i] > 0:
            surv *= 1.0 - 1.0 / at_risk[i]
        grid_t.append(times[i])
        grid_g.append(surv)
    grid_t_arr = np.asarray(grid_t)
    grid_g_arr = np.asarray(grid_g)

    def g_at(t: np.ndarray) -> np.ndarray:
        # G(t-): last censoring-survival value strictly before t.
        idx = np.searchsorted(grid_t_arr, t, side="left") - 1
        out = np.ones_like(t, dtype=float)
        valid = idx >= 0
        out[valid] = grid_g_arr[idx[valid]]
        return out

    eval_time = np.minimum(Y, horizon)
    kept = (D == 1) | (Y >= horizon)
    g = g_at(np.where(Y >= horizon, horizon, eval_time))
    g = np.clip(g, 1e-3, 1.0)
    weights = np.where(kept, 1.0 / g, 0.0)
    return eval_time, weights, kept


# --------------------------------------------------------------------------
# R plumbing
# --------------------------------------------------------------------------
class RSession:
    def __init__(self) -> None:
        import rpy2.rinterface_lib.callbacks as callbacks
        import rpy2.robjects as ro
        from rpy2.robjects import default_converter, numpy2ri
        from rpy2.robjects.conversion import localconverter

        callbacks.consolewrite_print = lambda s: sys.stdout.write(s)
        callbacks.consolewrite_warnerror = lambda s: sys.stderr.write(s)

        self.ro = ro
        self._converter = default_converter + numpy2ri.converter
        self._localconverter = localconverter
        ro.r('suppressMessages({library(grf); library(lmtest); library(sandwich)})')

    def assign(self, name: str, value: Any) -> None:
        # numpy2ri turns a 1-D array into an R array with a `dim` attribute,
        # which grf rejects ("Observations must be vectors"), so pass those
        # through as plain R numeric vectors.
        if isinstance(value, np.ndarray) and value.ndim == 1:
            self.ro.globalenv[name] = self.ro.FloatVector(np.asarray(value, dtype=float))
            return
        with self._localconverter(self._converter):
            self.ro.globalenv[name] = value

    def assign_strings(self, name: str, values: list[str]) -> None:
        self.ro.globalenv[name] = self.ro.StrVector(values)

    def r(self, code: str):
        return self.ro.r(code)

    def num(self, code: str) -> float:
        return float(np.asarray(self.ro.r(code))[0])

    def vec(self, code: str) -> np.ndarray:
        return np.asarray(self.ro.r(code), dtype=float)

    def strs(self, code: str) -> list[str]:
        return [str(v) for v in self.ro.r(code)]


CALIBRATION_HELPER = """
# Best-linear-projection calibration test on doubly-robust scores.  This is
# grf::test_calibration for forests it does not natively support (notably
# causal_survival_forest): same regressors, same sandwich variance, same
# one-sided p-value convention.
dr_calibration <- function(forest, vcov.type = "HC3") {
  preds <- predict(forest)$predictions
  ow <- if (length(forest$sample.weights) > 0) forest$sample.weights
        else rep(1, length(preds))
  mean.pred <- weighted.mean(preds, ow)
  DF <- data.frame(
    target = grf::get_scores(forest),
    mean.forest.prediction = mean.pred,
    differential.forest.prediction = preds - mean.pred
  )
  blp <- lm(target ~ mean.forest.prediction + differential.forest.prediction + 0,
            weights = ow, data = DF)
  out <- lmtest::coeftest(blp, vcov = sandwich::vcovCL, type = vcov.type,
                          cluster = seq_len(nrow(DF)))
  out[, 4] <- ifelse(out[, 3] < 0, 1 - out[, 4] / 2, out[, 4] / 2)
  out
}
"""


ROW_SUBSETTED_ARGS = {"sample.weights", "W.hat"}


def forest_call(kind: str, args: dict[str, str], subset: str | None = None) -> str:
    """Build the grf forest call, optionally restricted to a subset of rows."""
    if subset is None:
        matrix, vector = "X", lambda name: name
    else:
        matrix = f"X[{subset}, , drop = FALSE]"
        vector = lambda name: f"{name}[{subset}]"  # noqa: E731
    parts = [matrix, vector("Y"), vector("W")]
    if kind == "survival":
        parts.append(vector("D"))
    for key, value in args.items():
        if subset is not None and key in ROW_SUBSETTED_ARGS:
            value = f"{value}[{subset}]"
        parts.append(f"{key} = {value}")
    function = "causal_survival_forest" if kind == "survival" else "causal_forest"
    return f"{function}({', '.join(parts)})"


def stratified_halves(W: np.ndarray, D: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Split rows in half, keeping treatment/event composition balanced."""
    rng = np.random.default_rng(seed)
    first: list[np.ndarray] = []
    second: list[np.ndarray] = []
    for w in (0, 1):
        for d in (0, 1):
            idx = np.flatnonzero((W == w) & (D == d))
            rng.shuffle(idx)
            cut = len(idx) // 2
            first.append(idx[:cut])
            second.append(idx[cut:])
    return np.sort(np.concatenate(first)), np.sort(np.concatenate(second))


# --------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------
def run(args: argparse.Namespace) -> dict[str, Any]:
    r_info = bootstrap_r(args.r_home)
    cohort = load_cohort(args.features, args.splits, args.top_var_genes)
    n, p = cohort.X.shape
    print(f"[data] n = {n}, p = {p}, treated = {int(cohort.W.sum())}, "
          f"events = {int(cohort.D.sum())}")
    if p > n:
        print(f"[data] p >> n regime ({p} features vs {n} patients)")

    # Outcome construction happens before anything is handed to R: the IPCW
    # pseudo-outcome drops patients censored before the horizon, which changes
    # the analytic cohort.
    sample_weights = None
    if args.forest == "survival":
        outcome = cohort.Y
        outcome_label = f"RMST at {args.horizon:g} months (causal_survival_forest)"
    elif args.causal_outcome == "ipcw-rmst":
        outcome, weights, kept = ipcw_rmst_outcome(cohort.Y, cohort.D, args.horizon)
        dropped = int((~kept).sum())
        if dropped:
            print(f"[data] IPCW: dropping {dropped} patients censored before "
                  f"{args.horizon:g} months")
        cohort = cohort.subset(kept)
        outcome, sample_weights = outcome[kept], weights[kept]
        outcome_label = (
            f"IPCW restricted mean survival time at {args.horizon:g} months "
            f"({dropped} patients censored before the horizon were dropped)"
        )
    else:
        outcome = cohort.Y
        outcome_label = "raw follow-up time (censoring ignored - biased)"
        print("[warn] --causal-outcome raw-time ignores censoring; "
              "treat the estimates as descriptive only")

    n, p = cohort.X.shape

    outdir = Path(args.outdir) if args.outdir else HERE / "results" / (
        f"{args.features}_{args.forest}_{'-'.join(args.splits)}"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    session = RSession()
    session.r(CALIBRATION_HELPER)
    session.assign("X", cohort.X)
    session.assign_strings("feature_names", cohort.feature_names)
    session.r("colnames(X) <- feature_names")
    session.assign("Y", outcome.astype(float))
    session.assign("W", cohort.W.astype(float))
    session.assign("D", cohort.D.astype(float))

    results: dict[str, Any] = {
        "r_version": r_info.get("_R_VERSION"),
        "grf_version": r_info.get("_GRF_VERSION"),
        "features": args.features,
        "splits": args.splits,
        "forest": args.forest,
        "outcome": outcome_label,
        "n": n,
        "p": p,
        "n_treated": int(cohort.W.sum()),
        "n_events": int(cohort.D.sum()),
        "horizon_months": args.horizon,
        "seed": args.seed,
    }

    forest_args: dict[str, str] = {
        "num.trees": str(args.num_trees),
        "min.node.size": str(args.min_node_size),
        "honesty": "TRUE",
        "tune.parameters": f'"{args.tune}"',
        "seed": str(args.seed),
    }
    if args.threads:
        forest_args["num.threads"] = str(args.threads)
    if args.forest == "survival":
        forest_args["target"] = '"RMST"'
        forest_args["horizon"] = str(args.horizon)
    if sample_weights is not None:
        session.assign("ipcw", sample_weights)
        forest_args["sample.weights"] = "ipcw"

    # Optional clinical-only propensity model, for when a 13k-gene propensity
    # forest is not something you want to lean on.
    if args.propensity_features == "clinical":
        if not cohort.clinical_index:
            raise SystemExit("--propensity-features clinical needs --features clinical or both")
        cols = ",".join(str(i + 1) for i in cohort.clinical_index)
        session.r(f"Xclin <- X[, c({cols}), drop = FALSE]")
        session.r(f"w.forest <- regression_forest(Xclin, W, num.trees = 2000, seed = {args.seed})")
        session.r("W.hat.fixed <- predict(w.forest)$predictions")
        forest_args["W.hat"] = "W.hat.fixed"
        results["propensity_model"] = "regression_forest on clinical covariates only"
    else:
        results["propensity_model"] = "grf internal regression forest on all features"

    call = forest_call(args.forest, forest_args)
    print(f"[grf] fitting {call.split('(')[0]} "
          f"(num.trees={args.num_trees}, tune.parameters={args.tune}) ...", flush=True)
    started = time.time()
    session.r(f"set.seed({args.seed}); cf <- {call}")
    results["fit_seconds"] = round(time.time() - started, 1)
    print(f"[grf] fitted in {results['fit_seconds']:.1f}s")

    if args.tune != "none":
        # causal_forest reports what it tuned; causal_survival_forest does not.
        names = session.strs(
            "if (is.null(cf$tunable.params)) character(0) else names(cf$tunable.params)"
        )
        if names:
            values = session.vec("as.numeric(unlist(cf$tunable.params))")
            results["tuned_params"] = dict(zip(names, (float(v) for v in values)))
        else:
            results["tuned_params"] = None
            print("[grf] this forest does not report its tuned parameters; the "
                  "split-sample forests will re-tune inside each half")

    # ---- overlap -------------------------------------------------------
    w_hat = session.vec("cf$W.hat")
    treated = cohort.W == 1
    overlap = {
        "min": float(w_hat.min()),
        "max": float(w_hat.max()),
        "quantiles": {
            q: float(np.quantile(w_hat, float(q))) for q in ("0.01", "0.05", "0.5", "0.95", "0.99")
        },
        "share_below_0.05": float((w_hat < 0.05).mean()),
        "share_above_0.95": float((w_hat > 0.95).mean()),
        "share_outside_0.05_0.95": float(((w_hat < 0.05) | (w_hat > 0.95)).mean()),
        "mean_treated": float(w_hat[treated].mean()),
        "mean_control": float(w_hat[~treated].mean()),
        "treated_min": float(w_hat[treated].min()),
        "control_max": float(w_hat[~treated].max()),
    }
    results["overlap"] = overlap
    session.r(
        f'png("{outdir / "propensity_overlap.png"}", width = 1500, height = 1000, res = 150);'
        'op <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 1));'
        'br <- seq(min(cf$W.hat), max(cf$W.hat), length.out = 51);'
        'hist(cf$W.hat, breaks = br, col = "grey80", border = "white",'
        '     main = "Estimated propensity scores (all patients)", xlab = "W.hat");'
        'h1 <- hist(cf$W.hat[W == 1], breaks = br, plot = FALSE);'
        'h0 <- hist(cf$W.hat[W == 0], breaks = br, plot = FALSE);'
        'plot(h0, col = "#2980b980", border = "white",'
        '     ylim = c(0, max(h0$counts, h1$counts)),'
        '     main = "Treated (red) vs control (blue)", xlab = "W.hat");'
        'plot(h1, col = "#c0392b80", border = "white", add = TRUE);'
        'par(op); invisible(dev.off())'
    )

    # ---- ATE -----------------------------------------------------------
    try:
        ate = session.vec("as.numeric(average_treatment_effect(cf))")
        results["ate"] = {"estimate": float(ate[0]), "std_err": float(ate[1])}
    except Exception as exc:
        results["ate"] = {"error": str(exc)}

    # ---- calibration ---------------------------------------------------
    if args.forest == "causal":
        session.r(f'tc <- test_calibration(cf, vcov.type = "{args.vcov_type}")')
        results["calibration_method"] = f"grf::test_calibration (vcov.type={args.vcov_type})"
    else:
        session.r(f'tc <- dr_calibration(cf, vcov.type = "{args.vcov_type}")')
        results["calibration_method"] = (
            f"best linear projection on doubly-robust scores, sandwich {args.vcov_type} "
            "(grf::test_calibration does not support causal_survival_forest)"
        )
    tc = session.vec("as.numeric(as.matrix(tc))").reshape(2, 4, order="F")
    rownames = session.strs("rownames(tc)")
    calibration = {
        name: {
            "estimate": float(tc[i, 0]),
            "std_err": float(tc[i, 1]),
            "t_value": float(tc[i, 2]),
            "p_value_one_sided": float(tc[i, 3]),
        }
        for i, name in enumerate(rownames)
    }
    results["calibration"] = calibration

    # ---- CATE distribution --------------------------------------------
    tau_hat = session.vec("predict(cf)$predictions")
    results["tau_hat"] = {
        "mean": float(tau_hat.mean()),
        "sd": float(tau_hat.std(ddof=1)),
        "min": float(tau_hat.min()),
        "q25": float(np.quantile(tau_hat, 0.25)),
        "median": float(np.median(tau_hat)),
        "q75": float(np.quantile(tau_hat, 0.75)),
        "max": float(tau_hat.max()),
        "share_positive": float((tau_hat > 0).mean()),
    }
    pd.DataFrame({
        "split": cohort.split,
        "W": cohort.W,
        "OS_MONTHS": cohort.Y,
        "OS_STATUS": cohort.D,
        "W_hat": w_hat,
        "tau_hat": tau_hat,
    }).to_csv(outdir / "oob_predictions.csv", index=False)

    # ---- variable importance ------------------------------------------
    if args.top_importance:
        session.r("vi <- variable_importance(cf)")
        importance = session.vec("as.numeric(vi)")
        order = np.argsort(importance)[::-1][: args.top_importance]
        results["variable_importance_top"] = [
            {"feature": cohort.feature_names[i], "importance": float(importance[i])}
            for i in order
        ]

    # ---- RATE / AUTOC --------------------------------------------------
    rate: dict[str, Any] = {}
    for target in ("AUTOC", "QINI"):
        session.r(
            f'set.seed({args.seed}); rate <- rank_average_treatment_effect('
            f'cf, priorities = predict(cf)$predictions, target = "{target}", R = {args.rate_boot})'
        )
        est = session.num("rate$estimate")
        se = session.num("rate$std.err")
        rate[target] = _rate_row(est, se)
        if target == "AUTOC":
            # plot() draws the TOC curve, which is the same object for both
            # targets; only the weighting of its area differs.
            session.r(
                f'png("{outdir / "toc_in_sample.png"}", width = 1400, height = 1000, res = 150);'
                'plot(rate, main = "TOC, priorities = out-of-bag CATE (in-sample)");'
                'invisible(dev.off())'
            )
    results["rate_in_sample"] = rate
    results["rate_in_sample_note"] = (
        "Priorities are the forest's own out-of-bag predictions evaluated on the "
        "same rows, so these are optimistically biased; use the split-sample RATE."
    )

    # Honest version: rank on a forest fit to one half, evaluate on the other.
    split_rate: list[dict[str, Any]] = []
    if args.rate_repeats > 0:
        split_args = dict(forest_args)
        # Reuse the main forest's tuned hyperparameters when grf reports them,
        # otherwise let each half tune the same way the main forest did.
        tuned = results.get("tuned_params") or {}
        if tuned:
            split_args["tune.parameters"] = '"none"'
            for key, value in tuned.items():
                if key in {"mtry", "min.node.size"}:
                    split_args[key] = str(int(round(value)))
                elif key == "honesty.prune.leaves":
                    split_args[key] = "TRUE" if value >= 0.5 else "FALSE"
                elif key in {"sample.fraction", "honesty.fraction", "alpha", "imbalance.penalty"}:
                    split_args[key] = f"{value:.6g}"
        for repeat in range(args.rate_repeats):
            seed = args.seed + 1000 * (repeat + 1)
            first, second = stratified_halves(cohort.W, cohort.D, seed)
            session.assign("idx.a", (first + 1).astype(float))
            session.assign("idx.b", (second + 1).astype(float))
            print(f"[grf] split-sample RATE repeat {repeat + 1}/{args.rate_repeats} ...", flush=True)
            fit_args = dict(split_args, seed=str(seed))
            session.r(f"set.seed({seed}); f.a <- {forest_call(args.forest, fit_args, 'idx.a')}")
            session.r(f"set.seed({seed}); f.b <- {forest_call(args.forest, fit_args, 'idx.b')}")
            session.r("prio.b <- predict(f.a, X[idx.b, , drop = FALSE])$predictions")
            session.r("prio.a <- predict(f.b, X[idx.a, , drop = FALSE])$predictions")
            for label, forest, priorities in (("a_on_b", "f.b", "prio.b"), ("b_on_a", "f.a", "prio.a")):
                session.r(
                    f'set.seed({seed}); rate <- rank_average_treatment_effect('
                    f'{forest}, priorities = {priorities}, target = "AUTOC", R = {args.rate_boot})'
                )
                row = _rate_row(session.num("rate$estimate"), session.num("rate$std.err"))
                row.update({"repeat": repeat + 1, "direction": label, "seed": seed})
                split_rate.append(row)
                if repeat == 0 and label == "a_on_b":
                    session.r(
                        f'png("{outdir / "toc_split_sample.png"}", '
                        'width = 1400, height = 1000, res = 150);'
                        'plot(rate, main = "TOC, split-sample priorities"); invisible(dev.off())'
                    )
        estimates = np.array([row["estimate"] for row in split_rate])
        results["rate_split_sample"] = {
            "repeats": split_rate,
            "mean_estimate": float(estimates.mean()),
            "median_p_value": float(np.median([row["p_value_two_sided"] for row in split_rate])),
            "share_p_below_0.05": float(
                np.mean([row["p_value_two_sided"] < 0.05 for row in split_rate])
            ),
            "note": "Halves are re-drawn per repeat, so repeats are correlated; "
                    "read the spread, not a pooled p-value.",
        }

    results["verdict"] = verdict(results)
    (outdir / "summary.json").write_text(json.dumps(results, indent=2))
    report = format_report(results)
    (outdir / "report.txt").write_text(report)
    print("\n" + report)
    print(f"[out] {outdir}")
    return results


def _rate_row(estimate: float, std_err: float) -> dict[str, Any]:
    from math import erfc, sqrt
    z = estimate / std_err if std_err > 0 else 0.0
    return {
        "estimate": estimate,
        "std_err": std_err,
        "z": z,
        "p_value_two_sided": float(erfc(abs(z) / sqrt(2.0))),
    }


def verdict(results: dict[str, Any]) -> str:
    lines = []
    overlap = results["overlap"]
    treated_share = results["n_treated"] / results["n"]
    if overlap["share_outside_0.05_0.95"] > 0.10:
        side = "below 0.05" if overlap["share_below_0.05"] >= overlap["share_above_0.95"] \
            else "above 0.95"
        lines.append(
            f"OVERLAP: THIN - {overlap['share_outside_0.05_0.95']:.1%} of patients have "
            f"propensity {side} (treated share is {treated_share:.1%}, so W.hat is "
            f"low throughout; W.hat maxes out at {overlap['max']:.2f}, with no mass "
            "near 1). Effects for those patients rest on few comparable treated cases."
        )
    elif overlap["min"] < 0.02 or overlap["max"] > 0.98:
        lines.append("OVERLAP: MARGINAL - a few patients sit near 0 or 1.")
    else:
        lines.append("OVERLAP: OK - no propensity mass near 0 or 1.")

    diff = results["calibration"].get("differential.forest.prediction")
    if diff:
        p = diff["p_value_one_sided"]
        if p < 0.05:
            lines.append(
                f"CALIBRATION: HTE signal (differential.forest.prediction = "
                f"{diff['estimate']:.3f}, one-sided p = {p:.4f})."
            )
        else:
            lines.append(
                f"CALIBRATION: no HTE signal (differential.forest.prediction = "
                f"{diff['estimate']:.3f}, one-sided p = {p:.4f})."
            )
    split = results.get("rate_split_sample")
    if split:
        lines.append(
            f"RATE (split-sample AUTOC): mean estimate {split['mean_estimate']:.3f}, "
            f"median p = {split['median_p_value']:.3f}, "
            f"{split['share_p_below_0.05']:.0%} of splits significant at 0.05."
        )
    return "\n".join(lines)


def format_report(r: dict[str, Any]) -> str:
    out: list[str] = []
    add = out.append
    add("=" * 78)
    add("grf heterogeneous-treatment-effect screen")
    add("=" * 78)
    add(f"{r['r_version']}, grf {r['grf_version']}")
    add(f"cohort            : {', '.join(r['splits'])} (test set never read)")
    add(f"features          : {r['features']} (p = {r['p']}, n = {r['n']})")
    add(f"treated / control : {r['n_treated']} / {r['n'] - r['n_treated']}")
    add(f"deaths observed   : {r['n_events']}")
    add(f"outcome           : {r['outcome']}")
    add(f"propensity        : {r['propensity_model']}")
    add(f"fit time          : {r['fit_seconds']}s")
    if r.get("tuned_params"):
        tuned = ", ".join(f"{k}={v:g}" for k, v in r["tuned_params"].items())
        add(f"tuned params      : {tuned}")

    o = r["overlap"]
    add("")
    add("-- Overlap (estimated propensity W.hat) " + "-" * 38)
    add(f"  range              : [{o['min']:.4f}, {o['max']:.4f}]")
    add("  quantiles          : " + ", ".join(f"q{q}={v:.4f}" for q, v in o["quantiles"].items()))
    add(f"  outside [.05,.95]  : {o['share_outside_0.05_0.95']:.2%} "
        f"(below {o['share_below_0.05']:.2%}, above {o['share_above_0.95']:.2%})")
    add(f"  mean W.hat treated : {o['mean_treated']:.4f}   control: {o['mean_control']:.4f}")

    if isinstance(r.get("ate"), dict) and "estimate" in r["ate"]:
        add("")
        add("-- Average treatment effect " + "-" * 50)
        est, se = r["ate"]["estimate"], r["ate"]["std_err"]
        add(f"  ATE = {est:+.4f} (SE {se:.4f}, 95% CI [{est - 1.96 * se:+.4f}, "
            f"{est + 1.96 * se:+.4f}]), in months of survival restricted to "
            f"{r['horizon_months']:g} months")

    add("")
    add("-- Calibration test " + "-" * 58)
    add(f"  {r['calibration_method']}")
    add(f"  {'coefficient':<34}{'estimate':>11}{'SE':>10}{'t':>8}{'p(1-sided)':>13}")
    for name, row in r["calibration"].items():
        add(f"  {name:<34}{row['estimate']:>11.4f}{row['std_err']:>10.4f}"
            f"{row['t_value']:>8.2f}{row['p_value_one_sided']:>13.4f}")
    add("  mean.forest.prediction ~ 1 means the forest's average effect is calibrated;")
    add("  differential.forest.prediction > 0 with a small p-value is the HTE evidence.")
    if abs(r["tau_hat"]["mean"]) < 0.25 * r["tau_hat"]["sd"]:
        add("  NOTE: the mean CATE is near zero, so the mean.forest.prediction coefficient")
        add("  (a ratio with that mean in the denominator) is numerically unstable.  The")
        add("  two regressors are orthogonal, so the differential row is unaffected.")

    t = r["tau_hat"]
    add("")
    add("-- Out-of-bag CATE estimates " + "-" * 49)
    add(f"  mean {t['mean']:+.4f}  sd {t['sd']:.4f}  "
        f"[min {t['min']:+.4f}, q25 {t['q25']:+.4f}, med {t['median']:+.4f}, "
        f"q75 {t['q75']:+.4f}, max {t['max']:+.4f}]")
    add(f"  share with predicted benefit from chemotherapy: {t['share_positive']:.1%}")

    add("")
    add("-- RATE " + "-" * 70)
    for target, row in r["rate_in_sample"].items():
        add(f"  {target:<6} in-sample   : {row['estimate']:+.4f} (SE {row['std_err']:.4f}, "
            f"two-sided p = {row['p_value_two_sided']:.4f})")
    add("  in-sample RATE reuses the same rows for ranking and evaluation (optimistic)")
    split = r.get("rate_split_sample")
    if split:
        for row in split["repeats"]:
            add(f"  AUTOC  split {row['repeat']}/{row['direction']:<6}: "
                f"{row['estimate']:+.4f} (SE {row['std_err']:.4f}, "
                f"two-sided p = {row['p_value_two_sided']:.4f})")
        add(f"  mean split-sample AUTOC : {split['mean_estimate']:+.4f}, "
            f"median p = {split['median_p_value']:.4f}")

    if r.get("variable_importance_top"):
        add("")
        add("-- Top split-frequency importance " + "-" * 44)
        add("  " + ", ".join(f"{d['feature']} ({d['importance']:.4f})"
                             for d in r["variable_importance_top"]))

    add("")
    add("-- Verdict " + "-" * 67)
    for line in r["verdict"].splitlines():
        add(f"  {line}")
    add("=" * 78)
    return "\n".join(out)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--features", choices=["genomic", "clinical", "both"], default="genomic",
                        help="covariate matrix X (default: genomic)")
    parser.add_argument("--splits", default="train,validation",
                        help="comma-separated development splits to pool "
                             "(default: train,validation; the test set is sealed)")
    parser.add_argument("--forest", choices=["survival", "causal"], default="survival",
                        help="causal_survival_forest (default, handles censoring) or causal_forest")
    parser.add_argument("--causal-outcome", choices=["ipcw-rmst", "raw-time"], default="ipcw-rmst",
                        help="outcome for --forest causal (default: IPCW restricted mean)")
    parser.add_argument("--horizon", type=float, default=60.0,
                        help="RMST horizon in months (default: 60, the repository's tau)")
    parser.add_argument("--num-trees", type=int, default=4000,
                        help="more trees stabilize OOB predictions when p >> n (default: 4000)")
    parser.add_argument("--min-node-size", type=int, default=5)
    parser.add_argument("--tune", choices=["all", "none"], default="all",
                        help="grf tune.parameters (default: all)")
    parser.add_argument("--top-var-genes", type=int, default=0,
                        help="keep only the K highest-variance genes (0 = keep all, default)")
    parser.add_argument("--propensity-features", choices=["all", "clinical"], default="all",
                        help="features for the W.hat model (default: all, grf's internal forest)")
    parser.add_argument("--vcov-type", default="HC3",
                        help="sandwich variance type for the calibration test (default: HC3)")
    parser.add_argument("--rate-repeats", type=int, default=3,
                        help="split-sample RATE repeats, 0 to skip (default: 3)")
    parser.add_argument("--rate-boot", type=int, default=500,
                        help="bootstrap replications inside rank_average_treatment_effect")
    parser.add_argument("--top-importance", type=int, default=25,
                        help="how many features to list by split-frequency importance")
    parser.add_argument("--threads", type=int, default=0, help="grf num.threads (0 = all cores)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--r-home", default=None,
                        help="R installation to embed; auto-detected if omitted")
    args = parser.parse_args(argv)
    args.splits = [s.strip() for s in args.splits.split(",") if s.strip()]
    unknown = [s for s in args.splits if s not in CLINICAL_CSV]
    if unknown:
        parser.error(f"unknown split(s) {unknown}; the test set is sealed (see red_lines.md)")
    return args


if __name__ == "__main__":
    run(parse_args())
