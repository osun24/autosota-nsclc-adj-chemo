"""Locked evaluation for the Reactome-wide RSF permutation-importance arena.

The former train and validation cohorts are pooled adaptive development data.
All model scores are out-of-fold. The sealed test set is outside this module's
path contract.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import pickle
import tempfile
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rsf_pfi_mpl"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import optuna
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sksurv.ensemble import RandomSurvivalForest
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv

from reactome_rsf.config import MSIGDB_FILENAME
from reactome_rsf.msigdb import parse_gmt, retrieve_collection, sha256_file


ARENA_DIR = Path(__file__).resolve().parent
REPO_ROOT = ARENA_DIR.parent
TRAIN_CSV = REPO_ROOT / "affyfRMATrain.csv"
VALID_CSV = REPO_ROOT / "affyfRMAValidation.csv"
REACTOME_DATA_DIR = REPO_ROOT / "reactome_rsf" / "data"
BUDGET_PATH = ARENA_DIR / "budget.json"

OUTCOME_COLUMNS = ["OS_STATUS", "OS_MONTHS"]
TREATMENT = "Adjuvant Chemo"
CLINICAL_COLUMNS = [
    TREATMENT,
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
NUISANCE_COLUMNS = [
    "Age",
    "IS_MALE",
    "Stage_IB",
    "Stage_II",
    "Stage_III",
    "Histology_Adenosquamous Carcinoma",
    "Histology_Large Cell Carcinoma",
    "Histology_Squamous Cell Carcinoma",
    "Race_African American",
    "Race_Asian",
    "Race_Native Hawaiian or Other Pacific Islander",
    "Race_Unknown",
    "Smoked?_Unknown",
    "Smoked?_Yes",
]

SEARCH_KEYS = {
    "n_estimators", "max_depth", "min_samples_leaf",
    "split_leaf_multiplier", "max_features",
}
PROPENSITY_CLIP = (0.05, 0.95)
CENSORING_SURVIVAL_FLOOR = 0.05
WEIGHT_CLIP = (0.1, 5.0)
NUMERICAL_TOLERANCE = 1e-6


def load_budget() -> dict:
    return json.loads(BUDGET_PATH.read_text(encoding="utf-8"))


BUDGET = load_budget()
TAU = float(BUDGET["tau_months"])


@dataclass(frozen=True)
class CandidateSpec:
    name: str
    screening_space: dict
    reduced_space: dict
    plot_top_k: int

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "screening_space": self.screening_space,
            "reduced_space": self.reduced_space,
            "plot_top_k": self.plot_top_k,
        }


@dataclass
class FoldData:
    fold: int
    fit_index: np.ndarray
    assess_index: np.ndarray
    fit: pd.DataFrame
    assess: pd.DataFrame
    phi0: np.ndarray
    phi1: np.ndarray


class MedianTransformer:
    """Fit-only median imputation with stable feature ordering."""

    def __init__(self, feature_names: list[str]):
        self.feature_names = list(feature_names)
        self.medians_: np.ndarray | None = None

    def fit(self, frame: pd.DataFrame) -> "MedianTransformer":
        values = frame[self.feature_names].to_numpy(dtype=np.float32)
        with np.errstate(all="ignore"):
            medians = np.nanmedian(values, axis=0)
        if not np.isfinite(medians).all():
            bad = [
                name for name, value in zip(self.feature_names, medians)
                if not np.isfinite(value)
            ]
            raise ValueError(f"Features have no finite fitting values: {bad[:5]}")
        self.medians_ = medians.astype(np.float32)
        return self

    def transform(self, frame: pd.DataFrame) -> np.ndarray:
        if self.medians_ is None:
            raise RuntimeError("MedianTransformer must be fit before transform")
        values = frame[self.feature_names].to_numpy(dtype=np.float32)
        bad = ~np.isfinite(values)
        if bad.any():
            values[bad] = np.take(self.medians_, np.where(bad)[1])
        return np.ascontiguousarray(values, dtype=np.float32)


def _validate_space(space: object, label: str) -> dict:
    if not isinstance(space, dict) or set(space) != SEARCH_KEYS:
        raise ValueError(f"{label} must contain exactly {sorted(SEARCH_KEYS)}")
    normalized: dict[str, list] = {}
    for key in SEARCH_KEYS:
        raw = space[key]
        if not isinstance(raw, list) or not raw:
            raise ValueError(f"{label}.{key} must be a nonempty list")
        normalized[key] = list(raw)
    trees = [int(value) for value in normalized["n_estimators"]]
    depths = [int(value) for value in normalized["max_depth"]]
    leaves = [int(value) for value in normalized["min_samples_leaf"]]
    multipliers = [int(value) for value in normalized["split_leaf_multiplier"]]
    if min(trees) < 50 or max(trees) > int(BUDGET["max_trees_per_forest"]):
        raise ValueError(f"{label}.n_estimators is outside the fixed budget")
    if min(depths) < 2 or max(depths) > 16:
        raise ValueError(f"{label}.max_depth must be in [2, 16]")
    if min(leaves) < 3 or max(leaves) > 150:
        raise ValueError(f"{label}.min_samples_leaf must be in [3, 150]")
    if min(multipliers) < 2 or max(multipliers) > 8:
        raise ValueError(f"{label}.split_leaf_multiplier must be in [2, 8]")
    for value in normalized["max_features"]:
        if isinstance(value, str):
            if value not in {"sqrt", "log2"}:
                raise ValueError(f"Unsupported {label}.max_features value: {value}")
        elif not 0.001 <= float(value) <= 1.0:
            raise ValueError(f"{label}.max_features fractions must be in [0.001, 1]")
    return {
        "n_estimators": trees,
        "max_depth": depths,
        "min_samples_leaf": leaves,
        "split_leaf_multiplier": multipliers,
        "max_features": normalized["max_features"],
    }


def validate_candidate(raw: dict) -> CandidateSpec:
    expected = {"name", "screening_space", "reduced_space", "plot_top_k"}
    if not isinstance(raw, dict) or set(raw) != expected:
        raise ValueError(f"Candidate must contain exactly {sorted(expected)}")
    name = raw["name"]
    if not isinstance(name, str) or not name.strip():
        raise ValueError("Candidate name must be a nonempty string")
    plot_top_k = int(raw["plot_top_k"])
    if not 5 <= plot_top_k <= 100:
        raise ValueError("plot_top_k must be in [5, 100]")
    return CandidateSpec(
        name=name.strip(),
        screening_space=_validate_space(raw["screening_space"], "screening_space"),
        reduced_space=_validate_space(raw["reduced_space"], "reduced_space"),
        plot_top_k=plot_top_k,
    )


def _assert_allowed_search_path(path: str | Path) -> None:
    resolved = Path(path).resolve()
    if resolved not in {TRAIN_CSV.resolve(), VALID_CSV.resolve()}:
        raise ValueError(f"Search loader accepts only train and validation CSVs: {resolved}")
    if "test" in resolved.name.lower():
        raise ValueError(f"Test paths are forbidden during search: {resolved}")


def _prepare_frame(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    mapped = result[TREATMENT].map({"OBS": 0, "ACT": 1})
    result[TREATMENT] = mapped.where(mapped.notna(), result[TREATMENT])
    for column in result.columns:
        result[column] = pd.to_numeric(result[column], errors="coerce")
    required = OUTCOME_COLUMNS + CLINICAL_COLUMNS
    if result[required].isna().any().any():
        bad = result[required].columns[result[required].isna().any()].tolist()
        raise ValueError(f"Missing required analysis values: {bad}")
    result[TREATMENT] = result[TREATMENT].astype(np.int8)
    result["OS_STATUS"] = result["OS_STATUS"].astype(np.int8)
    if not set(result[TREATMENT].unique()).issubset({0, 1}):
        raise ValueError("Adjuvant Chemo must be binary")
    if not set(result["OS_STATUS"].unique()).issubset({0, 1}):
        raise ValueError("OS_STATUS must be binary")
    if (result["OS_MONTHS"] <= 0).any():
        raise ValueError("OS_MONTHS must be positive")
    return result


def load_reactome_sets(data_dir: str | Path = REACTOME_DATA_DIR) -> dict[str, tuple[str, ...]]:
    retrieve_collection(data_dir)
    return parse_gmt(Path(data_dir) / MSIGDB_FILENAME)


def load_development(
    train_csv: str | Path = TRAIN_CSV,
    valid_csv: str | Path = VALID_CSV,
    data_dir: str | Path = REACTOME_DATA_DIR,
) -> tuple[pd.DataFrame, list[str]]:
    _assert_allowed_search_path(train_csv)
    _assert_allowed_search_path(valid_csv)
    pathways = load_reactome_sets(data_dir)
    train_header = set(pd.read_csv(train_csv, nrows=0).columns)
    valid_header = set(pd.read_csv(valid_csv, nrows=0).columns)
    reactome = {gene for members in pathways.values() for gene in members}
    genes = sorted(reactome.intersection(train_header).intersection(valid_header))
    requested = OUTCOME_COLUMNS + CLINICAL_COLUMNS + genes
    train = _prepare_frame(pd.read_csv(train_csv, usecols=requested))
    valid = _prepare_frame(pd.read_csv(valid_csv, usecols=requested))
    development = pd.concat([train, valid], ignore_index=True)
    values = development[genes].to_numpy(dtype=np.float32)
    with np.errstate(all="ignore"):
        spans = np.nanmax(values, axis=0) - np.nanmin(values, axis=0)
    genes = [gene for gene, span in zip(genes, spans) if np.isfinite(span) and span > 0]
    if not genes:
        raise RuntimeError("No nonconstant Reactome genes are available")
    return development[OUTCOME_COLUMNS + CLINICAL_COLUMNS + genes], genes


def fixed_folds(frame: pd.DataFrame, n_splits: int, seed: int):
    strata = frame["OS_STATUS"].astype(str) + "_" + frame[TREATMENT].astype(str)
    splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    yield from splitter.split(np.zeros(len(frame)), strata)


def _outcome(frame: pd.DataFrame, *, censoring: bool = False) -> np.ndarray:
    event = 1 - frame["OS_STATUS"].to_numpy(int) if censoring else frame["OS_STATUS"].to_numpy(int)
    return Surv.from_arrays(event=event.astype(bool), time=frame["OS_MONTHS"].to_numpy(float))


def _nuisance_matrix(
    frame: pd.DataFrame,
    *,
    treatment: np.ndarray | None = None,
    interactions: bool = False,
) -> np.ndarray:
    x = frame[NUISANCE_COLUMNS].to_numpy(dtype=float)
    a = frame[TREATMENT].to_numpy(dtype=float) if treatment is None else np.asarray(treatment, dtype=float)
    blocks = [a[:, None], x]
    if interactions:
        blocks.append(x * a[:, None])
    return np.column_stack(blocks)


def _fit_propensity(fit: pd.DataFrame, assess: pd.DataFrame) -> np.ndarray:
    model = make_pipeline(
        StandardScaler(), LogisticRegression(max_iter=2000, solver="lbfgs", C=0.5)
    )
    model.fit(fit[NUISANCE_COLUMNS].to_numpy(float), fit[TREATMENT].to_numpy(int))
    raw = model.predict_proba(assess[NUISANCE_COLUMNS].to_numpy(float))[:, 1]
    return np.clip(raw, *PROPENSITY_CLIP)


def _fit_training_weights(fit: pd.DataFrame) -> np.ndarray:
    propensity = _fit_propensity(fit, fit)
    treatment = fit[TREATMENT].to_numpy(int)
    prevalence = float(treatment.mean())
    weights = np.where(
        treatment == 1,
        prevalence / propensity,
        (1.0 - prevalence) / (1.0 - propensity),
    )
    return np.clip(weights, *WEIGHT_CLIP).astype(float)


def _integrate_step_functions(functions: Iterable, tau: float = TAU) -> np.ndarray:
    output = []
    for function in functions:
        knots = np.asarray(function.x, dtype=float)
        internal = knots[(knots > 0) & (knots < tau)]
        bounds = np.concatenate(([0.0], internal, [float(tau)]))
        left = bounds[:-1]
        survival = np.ones_like(left)
        positive = left > 0
        if positive.any():
            survival[positive] = np.asarray(function(left[positive]), dtype=float)
        output.append(float(np.sum(np.diff(bounds) * survival)))
    return np.asarray(output, dtype=float)


def _clinical_outcome_predictions(
    fit: pd.DataFrame, assess: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray]:
    model = make_pipeline(StandardScaler(), CoxPHSurvivalAnalysis(alpha=1.0, ties="efron"))
    model.fit(_nuisance_matrix(fit, interactions=True), _outcome(fit))
    zeros = np.zeros(len(assess))
    ones = np.ones(len(assess))
    sf0 = model.predict_survival_function(
        _nuisance_matrix(assess, treatment=zeros, interactions=True)
    )
    sf1 = model.predict_survival_function(
        _nuisance_matrix(assess, treatment=ones, interactions=True)
    )
    return _integrate_step_functions(sf0), _integrate_step_functions(sf1)


def _ipcw_restricted_time(
    fit: pd.DataFrame, assess: pd.DataFrame, tau: float = TAU
) -> np.ndarray:
    model = make_pipeline(StandardScaler(), CoxPHSurvivalAnalysis(alpha=1.0, ties="efron"))
    model.fit(_nuisance_matrix(fit), _outcome(fit, censoring=True))
    functions = model.predict_survival_function(_nuisance_matrix(assess))
    observed = assess["OS_MONTHS"].to_numpy(float)
    result = np.empty(len(assess), dtype=float)
    for index, (function, observed_time) in enumerate(zip(functions, observed)):
        limit = float(min(observed_time, tau))
        knots = np.asarray(function.x, dtype=float)
        internal = knots[(knots > 0) & (knots < limit)]
        bounds = np.concatenate(([0.0], internal, [limit]))
        left = bounds[:-1]
        censor_survival = np.ones_like(left)
        positive = left > 0
        if positive.any():
            censor_survival[positive] = np.asarray(function(left[positive]), dtype=float)
        result[index] = float(np.sum(
            np.diff(bounds) / np.maximum(censor_survival, CENSORING_SURVIVAL_FLOOR)
        ))
    return result


def aipw_arm_scores(
    treatment: np.ndarray,
    propensity: np.ndarray,
    ipcw_time: np.ndarray,
    mu0: np.ndarray,
    mu1: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    treatment = np.asarray(treatment, dtype=int)
    propensity = np.asarray(propensity, dtype=float)
    phi1 = mu1 + treatment / propensity * (ipcw_time - mu1)
    phi0 = mu0 + (1 - treatment) / (1 - propensity) * (ipcw_time - mu0)
    return phi0, phi1


def make_fold_data(frame: pd.DataFrame, n_folds: int) -> list[FoldData]:
    folds: list[FoldData] = []
    seed = int(BUDGET["fold_seed"])
    for fold, (fit_index, assess_index) in enumerate(
        fixed_folds(frame, n_folds, seed), start=1
    ):
        fit = frame.iloc[fit_index].reset_index(drop=True)
        assess = frame.iloc[assess_index].reset_index(drop=True)
        propensity = _fit_propensity(fit, assess)
        ipcw_time = _ipcw_restricted_time(fit, assess)
        mu0, mu1 = _clinical_outcome_predictions(fit, assess)
        phi0, phi1 = aipw_arm_scores(
            assess[TREATMENT].to_numpy(int), propensity, ipcw_time, mu0, mu1
        )
        folds.append(FoldData(
            fold=fold,
            fit_index=np.asarray(fit_index),
            assess_index=np.asarray(assess_index),
            fit=fit,
            assess=assess,
            phi0=phi0,
            phi1=phi1,
        ))
    return folds


def _fit_rsf(
    fit: pd.DataFrame,
    feature_names: list[str],
    params: dict,
    seed: int,
) -> tuple[RandomSurvivalForest, MedianTransformer]:
    transformer = MedianTransformer(feature_names).fit(fit)
    model = RandomSurvivalForest(
        **params,
        bootstrap=True,
        oob_score=False,
        n_jobs=-1,
        random_state=int(seed),
        low_memory=False,
    )
    model.fit(
        transformer.transform(fit),
        _outcome(fit),
        sample_weight=_fit_training_weights(fit),
    )
    return model, transformer


def _rmst_from_survival_array(
    survival: np.ndarray, times: np.ndarray, tau: float = TAU
) -> np.ndarray:
    times = np.asarray(times, dtype=float)
    survival = np.asarray(survival, dtype=float)
    indices = np.flatnonzero((times > 0) & (times < tau))
    internal = times[indices]
    bounds = np.concatenate(([0.0], internal, [float(tau)]))
    values = np.column_stack([
        np.ones(survival.shape[0], dtype=float),
        survival[:, indices],
    ])
    return np.sum(values * np.diff(bounds)[None, :], axis=1)


def _predict_panel(
    model: RandomSurvivalForest,
    transformer: MedianTransformer,
    assess: pd.DataFrame,
) -> dict[str, np.ndarray]:
    observed = transformer.transform(assess)
    risk = np.asarray(model.predict(observed), dtype=float)
    rmst = []
    for treatment in (0, 1):
        counterfactual = assess.copy()
        counterfactual[TREATMENT] = treatment
        survival = model.predict_survival_function(
            transformer.transform(counterfactual), return_array=True
        )
        rmst.append(_rmst_from_survival_array(survival, model.unique_times_))
    return {"risk": risk, "rmst0": rmst[0], "rmst1": rmst[1]}


def recommendations(rmst0: np.ndarray, rmst1: np.ndarray) -> np.ndarray:
    benefit = np.asarray(rmst1, dtype=float) - np.asarray(rmst0, dtype=float)
    return (benefit > NUMERICAL_TOLERANCE).astype(np.int8)


def _policy_contribution(
    phi0: np.ndarray, phi1: np.ndarray, policy: np.ndarray
) -> np.ndarray:
    policy = np.asarray(policy, dtype=np.int8)
    chosen = np.where(policy == 1, phi1, phi0)
    anti = np.where(policy == 1, phi0, phi1)
    return chosen - anti


def _score_panel(fold: FoldData, panel: dict[str, np.ndarray]) -> dict:
    policy = recommendations(panel["rmst0"], panel["rmst1"])
    contribution = _policy_contribution(fold.phi0, fold.phi1, policy)
    cindex = float(concordance_index_censored(
        fold.assess["OS_STATUS"].to_numpy(bool),
        fold.assess["OS_MONTHS"].to_numpy(float),
        panel["risk"],
    )[0])
    return {
        "rmst_difference_months": float(np.mean(contribution)),
        "harrell_cindex": cindex,
        "act_recommended_fraction": float(np.mean(policy)),
        "contribution": contribution,
    }


def _evaluate_cv(
    frame: pd.DataFrame,
    folds: list[FoldData],
    features: list[str],
    params: dict,
) -> dict:
    risk = np.zeros(len(frame), dtype=float)
    contribution = np.zeros(len(frame), dtype=float)
    fold_results = []
    for fold in folds:
        model, transformer = _fit_rsf(
            fold.fit, features, params, int(BUDGET["forest_seed"]) + fold.fold
        )
        panel = _predict_panel(model, transformer, fold.assess)
        score = _score_panel(fold, panel)
        risk[fold.assess_index] = panel["risk"]
        contribution[fold.assess_index] = score.pop("contribution")
        fold_results.append({"fold": fold.fold, **score})
    cindex = float(concordance_index_censored(
        frame["OS_STATUS"].to_numpy(bool),
        frame["OS_MONTHS"].to_numpy(float),
        risk,
    )[0])
    return {
        "rmst_difference_months": float(np.mean(contribution)),
        "harrell_cindex": cindex,
        "act_recommended_fraction": float(np.mean([
            item["act_recommended_fraction"] for item in fold_results
        ])),
        "folds": fold_results,
    }


def _suggest_params(trial: optuna.Trial, space: dict, *, smoke: bool) -> dict:
    if smoke:
        return {
            "n_estimators": 20,
            "max_depth": 3,
            "min_samples_leaf": 10,
            "min_samples_split": 20,
            "max_features": "sqrt",
        }
    trees = int(trial.suggest_categorical("n_estimators", space["n_estimators"]))
    depth = int(trial.suggest_categorical("max_depth", space["max_depth"]))
    leaf = int(trial.suggest_categorical("min_samples_leaf", space["min_samples_leaf"]))
    multiplier = int(trial.suggest_categorical(
        "split_leaf_multiplier", space["split_leaf_multiplier"]
    ))
    max_features = trial.suggest_categorical("max_features", space["max_features"])
    return {
        "n_estimators": trees,
        "max_depth": depth,
        "min_samples_leaf": leaf,
        "min_samples_split": leaf * multiplier,
        "max_features": max_features,
    }


def _lexicographic_best(study: optuna.Study) -> optuna.trial.FrozenTrial:
    complete = [
        trial for trial in study.trials
        if trial.state == optuna.trial.TrialState.COMPLETE and trial.values is not None
    ]
    if not complete:
        raise RuntimeError("No completed optimization trials")
    best_primary = max(float(trial.values[0]) for trial in complete)
    tied = [
        trial for trial in complete
        if float(trial.values[0]) >= best_primary - NUMERICAL_TOLERANCE
    ]
    return max(tied, key=lambda trial: (float(trial.values[1]), -trial.number))


def _params_from_trial(trial: optuna.trial.FrozenTrial) -> dict:
    return dict(trial.user_attrs["rsf_params"])


def _study_rows(study: optuna.Study, stage: str) -> list[dict]:
    rows = []
    for trial in study.trials:
        if trial.state != optuna.trial.TrialState.COMPLETE or trial.values is None:
            continue
        rows.append({
            "stage": stage,
            "trial": int(trial.number),
            "top_n": trial.user_attrs.get("top_n"),
            "rmst_difference_months": float(trial.values[0]),
            "harrell_cindex": float(trial.values[1]),
            **{f"rsf_{key}": value for key, value in trial.user_attrs["rsf_params"].items()},
        })
    return rows


def _optimize_screening(
    frame: pd.DataFrame,
    folds: list[FoldData],
    genes: list[str],
    spec: CandidateSpec,
    *,
    smoke: bool,
) -> tuple[optuna.Study, optuna.trial.FrozenTrial]:
    features = CLINICAL_COLUMNS + genes
    sampler = optuna.samplers.NSGAIISampler(seed=int(BUDGET["forest_seed"]))
    study = optuna.create_study(directions=["maximize", "maximize"], sampler=sampler)

    def objective(trial: optuna.Trial) -> tuple[float, float]:
        params = _suggest_params(trial, spec.screening_space, smoke=smoke)
        result = _evaluate_cv(frame, folds, features, params)
        trial.set_user_attr("rsf_params", params)
        trial.set_user_attr("folds", result["folds"])
        print(
            f"[screening trial {trial.number:03d}] "
            f"RMST={result['rmst_difference_months']:.4f} "
            f"C={result['harrell_cindex']:.4f} params={params}",
            flush=True,
        )
        return result["rmst_difference_months"], result["harrell_cindex"]

    trials = 1 if smoke else int(BUDGET["screening_trials"])
    study.optimize(objective, n_trials=trials, gc_after_trial=True)
    return study, _lexicographic_best(study)


def _used_feature_indices(model: RandomSurvivalForest) -> set[int]:
    used: set[int] = set()
    for estimator in model.estimators_:
        features = np.asarray(estimator.tree_.feature, dtype=int)
        used.update(int(value) for value in features[features >= 0])
    return used


def _permutation_importance(
    folds: list[FoldData],
    feature_names: list[str],
    params: dict,
    repeats: int,
) -> tuple[pd.DataFrame, dict]:
    n_features = len(feature_names)
    count = len(folds) * repeats
    rmst_sum = np.zeros(n_features, dtype=float)
    rmst_sumsq = np.zeros(n_features, dtype=float)
    cindex_sum = np.zeros(n_features, dtype=float)
    cindex_sumsq = np.zeros(n_features, dtype=float)
    used_folds = np.zeros(n_features, dtype=int)
    baseline_contribution = []
    baseline_risk = []
    baseline_event = []
    baseline_time = []

    for fold in folds:
        model, transformer = _fit_rsf(
            fold.fit, feature_names, params,
            int(BUDGET["forest_seed"]) + fold.fold,
        )
        baseline_panel = _predict_panel(model, transformer, fold.assess)
        baseline = _score_panel(fold, baseline_panel)
        baseline_contribution.append(baseline.pop("contribution"))
        baseline_risk.append(baseline_panel["risk"])
        baseline_event.append(fold.assess["OS_STATUS"].to_numpy(bool))
        baseline_time.append(fold.assess["OS_MONTHS"].to_numpy(float))
        used = _used_feature_indices(model)
        used_folds[list(used)] += 1
        print(
            f"[PFI fold {fold.fold}] {len(used)}/{n_features} features used by a split",
            flush=True,
        )
        for feature_index in sorted(used):
            name = feature_names[feature_index]
            for repeat in range(repeats):
                rng = np.random.default_rng(
                    int(BUDGET["permutation_seed"])
                    + fold.fold * 10_000_000
                    + feature_index * 100
                    + repeat
                )
                permuted = fold.assess.copy()
                order = rng.permutation(len(permuted))
                permuted[name] = fold.assess[name].to_numpy()[order]
                score = _score_panel(fold, _predict_panel(model, transformer, permuted))
                score.pop("contribution")
                rmst_drop = baseline["rmst_difference_months"] - score["rmst_difference_months"]
                cindex_drop = baseline["harrell_cindex"] - score["harrell_cindex"]
                rmst_sum[feature_index] += rmst_drop
                rmst_sumsq[feature_index] += rmst_drop * rmst_drop
                cindex_sum[feature_index] += cindex_drop
                cindex_sumsq[feature_index] += cindex_drop * cindex_drop
            if (feature_index + 1) % 250 == 0:
                print(
                    f"[PFI fold {fold.fold}] processed through feature {feature_index + 1}",
                    flush=True,
                )

    rmst_mean = rmst_sum / count
    cindex_mean = cindex_sum / count
    denominator = max(1, count - 1)
    rmst_sd = np.sqrt(np.maximum(0.0, (rmst_sumsq - count * rmst_mean ** 2) / denominator))
    cindex_sd = np.sqrt(np.maximum(0.0, (cindex_sumsq - count * cindex_mean ** 2) / denominator))
    clinical = set(CLINICAL_COLUMNS)
    table = pd.DataFrame({
        "feature": feature_names,
        "feature_type": ["clinical" if name in clinical else "gene" for name in feature_names],
        "rmst_importance_mean_months": rmst_mean,
        "rmst_importance_sd_months": rmst_sd,
        "cindex_importance_mean": cindex_mean,
        "cindex_importance_sd": cindex_sd,
        "split_used_fold_fraction": used_folds / len(folds),
    })
    table = table.sort_values(
        ["rmst_importance_mean_months", "cindex_importance_mean", "feature"],
        ascending=[False, False, True],
        kind="stable",
    ).reset_index(drop=True)
    table.insert(0, "overall_rank", np.arange(1, len(table) + 1))
    gene_mask = table["feature_type"].eq("gene")
    table["gene_rank"] = pd.array([pd.NA] * len(table), dtype="Int64")
    table.loc[gene_mask, "gene_rank"] = np.arange(1, int(gene_mask.sum()) + 1)
    pooled_contribution = np.concatenate(baseline_contribution)
    pooled_risk = np.concatenate(baseline_risk)
    baseline_summary = {
        "rmst_difference_months": float(np.mean(pooled_contribution)),
        "harrell_cindex": float(concordance_index_censored(
            np.concatenate(baseline_event),
            np.concatenate(baseline_time),
            pooled_risk,
        )[0]),
    }
    return table, baseline_summary


def rank_genes_from_importance(table: pd.DataFrame) -> list[str]:
    required = {
        "feature", "feature_type", "rmst_importance_mean_months",
        "cindex_importance_mean",
    }
    if required - set(table):
        raise ValueError(f"Importance table is missing: {sorted(required - set(table))}")
    genes = table.loc[table["feature_type"].eq("gene")].sort_values(
        ["rmst_importance_mean_months", "cindex_importance_mean", "feature"],
        ascending=[False, False, True],
        kind="stable",
    )
    return genes["feature"].astype(str).tolist()


def _optimize_reduced(
    frame: pd.DataFrame,
    folds: list[FoldData],
    ranked_genes: list[str],
    spec: CandidateSpec,
    *,
    smoke: bool,
) -> tuple[optuna.Study, optuna.trial.FrozenTrial]:
    max_n = min(2 if smoke else int(BUDGET["max_top_genes"]), len(ranked_genes))
    n_trials = max_n + 1 if smoke else int(BUDGET["reduced_trials"])
    if n_trials < max_n + 1:
        raise RuntimeError("Reduced trial budget must evaluate every top-N value")
    sampler = optuna.samplers.NSGAIISampler(seed=int(BUDGET["forest_seed"]) + 1)
    study = optuna.create_study(directions=["maximize", "maximize"], sampler=sampler)
    reference = {
        "n_estimators": spec.reduced_space["n_estimators"][0],
        "max_depth": spec.reduced_space["max_depth"][0],
        "min_samples_leaf": spec.reduced_space["min_samples_leaf"][0],
        "split_leaf_multiplier": spec.reduced_space["split_leaf_multiplier"][0],
        "max_features": spec.reduced_space["max_features"][0],
    }
    for top_n in range(max_n + 1):
        # The first pass is a controlled N curve: every panel uses identical
        # forest geometry. Remaining trials jointly refine N and geometry.
        study.enqueue_trial({"top_n": top_n, **({} if smoke else reference)})

    def objective(trial: optuna.Trial) -> tuple[float, float]:
        top_n = int(trial.suggest_int("top_n", 0, max_n))
        params = _suggest_params(trial, spec.reduced_space, smoke=smoke)
        features = CLINICAL_COLUMNS + ranked_genes[:top_n]
        result = _evaluate_cv(frame, folds, features, params)
        trial.set_user_attr("top_n", top_n)
        trial.set_user_attr("rsf_params", params)
        trial.set_user_attr("genes", ranked_genes[:top_n])
        trial.set_user_attr("folds", result["folds"])
        print(
            f"[reduced trial {trial.number:03d}] N={top_n:02d} "
            f"RMST={result['rmst_difference_months']:.4f} "
            f"C={result['harrell_cindex']:.4f} params={params}",
            flush=True,
        )
        return result["rmst_difference_months"], result["harrell_cindex"]

    study.optimize(objective, n_trials=n_trials, gc_after_trial=True)
    observed = {
        int(trial.user_attrs["top_n"])
        for trial in study.trials
        if trial.state == optuna.trial.TrialState.COMPLETE
    }
    missing = sorted(set(range(max_n + 1)) - observed)
    if missing:
        raise RuntimeError(f"Reduced optimization did not evaluate top-N values: {missing}")
    return study, _lexicographic_best(study)


def _plot_importance(table: pd.DataFrame, path: Path, top_k: int) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(15, max(6, 0.28 * top_k + 2)))
    panels = [
        ("rmst_importance_mean_months", "rmst_importance_sd_months", "RMST importance (months)"),
        ("cindex_importance_mean", "cindex_importance_sd", "C-index importance"),
    ]
    colors = {"clinical": "#d95f02", "gene": "#1b9e77"}
    for axis, (mean_col, sd_col, title) in zip(axes, panels):
        tie_break = (
            "cindex_importance_mean"
            if mean_col == "rmst_importance_mean_months"
            else "rmst_importance_mean_months"
        )
        panel = table.sort_values(
            [mean_col, tie_break, "feature"],
            ascending=[False, False, True],
            kind="stable",
        ).head(top_k).iloc[::-1]
        y = np.arange(len(panel))
        axis.barh(
            y,
            panel[mean_col],
            xerr=panel[sd_col],
            color=[colors[value] for value in panel["feature_type"]],
            alpha=0.85,
            capsize=2,
        )
        axis.set_yticks(y, panel["feature"])
        axis.axvline(0.0, color="black", linewidth=0.8)
        axis.set_xlabel("Score loss after held-out permutation")
        axis.set_title(title)
        axis.grid(axis="x", alpha=0.2)
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=colors[label], label=label.title())
        for label in ("gene", "clinical")
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.005),
        ncol=2,
        frameon=False,
    )
    fig.suptitle("Reactome-wide RSF permutation feature importance", y=0.995)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _json_safe(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate_candidate(
    raw_candidate: dict,
    *,
    artifact_dir: str | Path,
    smoke: bool = False,
) -> dict:
    spec = validate_candidate(raw_candidate)
    artifact_dir = Path(artifact_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    frame, all_genes = load_development()
    genes = all_genes[:64] if smoke else all_genes
    n_folds = 2 if smoke else int(BUDGET["cv_folds"])
    repeats = 1 if smoke else int(BUDGET["permutation_repeats"])
    folds = make_fold_data(frame, n_folds)

    screening_study, screening_best = _optimize_screening(
        frame, folds, genes, spec, smoke=smoke
    )
    screening_params = _params_from_trial(screening_best)
    feature_names = CLINICAL_COLUMNS + genes
    importance, screening_oof = _permutation_importance(
        folds, feature_names, screening_params, repeats
    )
    ranked_genes = rank_genes_from_importance(importance)
    reduced_study, reduced_best = _optimize_reduced(
        frame, folds, ranked_genes, spec, smoke=smoke
    )
    final_top_n = int(reduced_best.user_attrs["top_n"])
    final_genes = list(reduced_best.user_attrs["genes"])
    final_params = _params_from_trial(reduced_best)
    final_features = CLINICAL_COLUMNS + final_genes

    final_model, final_transformer = _fit_rsf(
        frame, final_features, final_params, int(BUDGET["forest_seed"])
    )
    importance_path = artifact_dir / "permutation_importance.csv"
    plot_path = artifact_dir / "permutation_importance.png"
    top_n_path = artifact_dir / "top_n_results.csv"
    final_spec_path = artifact_dir / "final_spec.json"
    covariates_path = artifact_dir / "final_covariates.txt"
    medians_path = artifact_dir / "training_imputation_medians.csv"
    model_path = artifact_dir / "final_rsf.pkl"

    importance.to_csv(importance_path, index=False)
    _plot_importance(importance, plot_path, min(spec.plot_top_k, len(importance)))
    pd.DataFrame(_study_rows(reduced_study, "reduced")).to_csv(top_n_path, index=False)
    covariates_path.write_text("\n".join(final_features) + "\n", encoding="utf-8")
    pd.DataFrame({
        "feature": final_features,
        "median": final_transformer.medians_,
    }).to_csv(medians_path, index=False)
    with model_path.open("wb") as stream:
        pickle.dump(final_model, stream, protocol=pickle.HIGHEST_PROTOCOL)

    final_cv = {
        "rmst_difference_months": float(reduced_best.values[0]),
        "harrell_cindex": float(reduced_best.values[1]),
        "folds": reduced_best.user_attrs["folds"],
    }
    final_spec = {
        "top_n": final_top_n,
        "genes": final_genes,
        "clinical_covariates": CLINICAL_COLUMNS,
        "covariates": final_features,
        "rsf_params": final_params,
        "development_cv": final_cv,
    }
    final_spec_path.write_text(
        json.dumps(_json_safe(final_spec), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    artifacts = {
        path.name: {"sha256": _sha256(path), "bytes": path.stat().st_size}
        for path in (
            importance_path, plot_path, top_n_path, final_spec_path,
            covariates_path, medians_path, model_path,
        )
    }
    max_n = min(2 if smoke else int(BUDGET["max_top_genes"]), len(ranked_genes))
    result = {
        "schema_version": 1,
        "phase": "smoke" if smoke else "pooled_development_cv_search",
        "candidate": spec.as_dict(),
        "objective": {
            "primary": "60-month cross-fitted IPCW-AIPW policy-versus-anti-policy RMST difference",
            "secondary": "Harrell C-index",
            "selection": "lexicographic with RMST first",
        },
        "screening": {
            "n_reactome_genes": len(genes),
            "n_features": len(feature_names),
            "best_trial": int(screening_best.number),
            "best_params": screening_params,
            "optimization_rmst_difference_months": float(screening_best.values[0]),
            "optimization_harrell_cindex": float(screening_best.values[1]),
            "pfi_refit_oof": screening_oof,
            "trials": _study_rows(screening_study, "screening"),
        },
        "permutation_importance": {
            "repeats": repeats,
            "rank_rule": "RMST importance descending, C-index importance descending, symbol ascending",
            "top_32_genes": ranked_genes[:32],
            "features_with_positive_rmst_importance": int(
                (importance["rmst_importance_mean_months"] > 0).sum()
            ),
            "features_used_in_any_fold": int(
                (importance["split_used_fold_fraction"] > 0).sum()
            ),
        },
        "reduced_search": {
            "evaluated_top_n": list(range(max_n + 1)),
            "best_trial": int(reduced_best.number),
            "trials": _study_rows(reduced_study, "reduced"),
        },
        "final": final_spec,
        "reward": float(reduced_best.values[0]),
        "secondary_score": float(reduced_best.values[1]),
        "eligible": True,
        "data": {
            "development_n": int(len(frame)),
            "development_events": int(frame["OS_STATUS"].sum()),
            "development_act": int(frame[TREATMENT].sum()),
            "reactome_genes_before_smoke_reduction": len(all_genes),
            "reactome_genes_used": len(genes),
            "former_validation_used_as_development": True,
            "confirmatory_validation_used": False,
            "test_used": False,
            "train_sha256": sha256_file(TRAIN_CSV),
            "validation_sha256": sha256_file(VALID_CSV),
            "reactome_collection_sha256": sha256_file(REACTOME_DATA_DIR / MSIGDB_FILENAME),
        },
        "budget": {**BUDGET, "smoke_reductions_applied": bool(smoke)},
        "artifacts": artifacts,
        "limitations": [
            "Train and former validation are pooled adaptive development data; the reported CV scores are not external confirmation.",
            "The PFI ranking and reduced-model search reuse development folds, so final CV performance is selection-adaptive.",
            "Permutation importance is conditional on the fitted RSF and can divide credit unpredictably among correlated genes.",
            "Zero importance for a feature unused by every assessment-fold forest is exact for those fitted forests, not proof of no biological effect.",
            "IPCW-AIPW validity requires adequate measured-confounder adjustment, positivity, consistency, and conditionally independent censoring.",
            "Random row-level CV does not establish transportability across source studies or expression platforms.",
        ],
    }
    return _json_safe(result)
