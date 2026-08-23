"""Locked evaluator for the Reactome RSF T-learner arena.

The arena reuses the sealed data, nuisance-model, and policy-value primitives
from ``reactome_rsf_arena_v2``.  Policy predictions are different: two
survival forests are fit in every training fold, one in each treatment arm,
and ACT is excluded from both feature matrices.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored

from reactome_rsf.msigdb import sha256_file
from reactome_rsf_arena_v2 import prepare as base


ARENA_DIR = Path(__file__).resolve().parent
REPO_ROOT = ARENA_DIR.parent
BUDGET_PATH = ARENA_DIR / "budget.json"

TRAIN_CSV = base.TRAIN_CSV
VALID_CSV = base.VALID_CSV
REACTOME_DATA_DIR = base.REACTOME_DATA_DIR
OUTCOME_COLUMNS = base.OUTCOME_COLUMNS
TREATMENT = base.TREATMENT
CLINICAL_COLUMNS = base.CLINICAL_COLUMNS
PRETREATMENT_COLUMNS = base.PRETREATMENT_COLUMNS
NUISANCE_COLUMNS = base.NUISANCE_COLUMNS
PROPENSITY_CLIP = base.PROPENSITY_CLIP
NUMERICAL_BENEFIT_TOLERANCE_MONTHS = base.NUMERICAL_BENEFIT_TOLERANCE_MONTHS
FOLD_SEED = base.FOLD_SEED

load_development = base.load_development
fixed_folds = base.fixed_folds
cross_fitted_benefit_pseudo_outcome = base.cross_fitted_benefit_pseudo_outcome
_gene_effect_scores = base._gene_effect_scores
_fit_propensity = base._fit_propensity
_fit_training_weights = base._fit_training_weights
_ipcw_restricted_time = base._ipcw_restricted_time
_integrate_step_functions = base._integrate_step_functions
_outcome = base._outcome
aipw_arm_scores = base.aipw_arm_scores
recommendations = base.recommendations
policy_summary = base.policy_summary
_bootstrap_increment = base._bootstrap_increment
_pairwise_jaccard = base._pairwise_jaccard


def load_budget() -> dict:
    return json.loads(BUDGET_PATH.read_text(encoding="utf-8"))


BUDGET = load_budget()
TAU = float(BUDGET["tau_months"])


@dataclass(frozen=True)
class ArmForestSpec:
    n_estimators: int
    max_depth: int
    min_samples_leaf: int
    min_samples_split: int
    max_features: float

    def as_dict(self) -> dict:
        return {
            "n_estimators": self.n_estimators,
            "max_depth": self.max_depth,
            "min_samples_leaf": self.min_samples_leaf,
            "min_samples_split": self.min_samples_split,
            "max_features": self.max_features,
        }


@dataclass(frozen=True)
class CandidateSpec:
    name: str
    selector: str
    n_genes: int
    representation: str
    module_count: int
    benefit_threshold_months: float
    observation: ArmForestSpec
    act: ArmForestSpec

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "selector": self.selector,
            "n_genes": self.n_genes,
            "representation": self.representation,
            "module_count": self.module_count,
            "benefit_threshold_months": self.benefit_threshold_months,
            "tlearner": {
                "observation": self.observation.as_dict(),
                "act": self.act.as_dict(),
            },
        }


def _validate_arm(raw: dict, label: str) -> ArmForestSpec:
    expected = {
        "n_estimators", "max_depth", "min_samples_leaf",
        "min_samples_split", "max_features",
    }
    if not isinstance(raw, dict) or set(raw) != expected:
        raise ValueError(f"{label} has missing or unrecognized RSF keys")
    trees = int(raw["n_estimators"])
    depth = int(raw["max_depth"])
    leaf = int(raw["min_samples_leaf"])
    split = int(raw["min_samples_split"])
    mtry = float(raw["max_features"])
    if not 100 <= trees <= int(BUDGET["max_trees_per_forest"]):
        raise ValueError(f"{label} n_estimators is outside the fixed budget")
    if not 3 <= depth <= 12:
        raise ValueError(f"{label} max_depth must be in [3, 12]")
    if not 5 <= leaf <= 40:
        raise ValueError(f"{label} min_samples_leaf must be in [5, 40]")
    if split < 2 * leaf or split > 100:
        raise ValueError(f"{label} min_samples_split must be >= 2*leaf and <= 100")
    if not 0.20 <= mtry <= 1.0:
        raise ValueError(f"{label} max_features must be in [0.20, 1.0]")
    return ArmForestSpec(trees, depth, leaf, split, mtry)


def validate_candidate(raw: dict) -> CandidateSpec:
    expected = {
        "name", "selector", "n_genes", "representation", "module_count",
        "benefit_threshold_months", "tlearner",
    }
    if set(raw) != expected:
        raise ValueError("Candidate has missing or unrecognized top-level keys")
    if not isinstance(raw["name"], str) or not raw["name"].strip():
        raise ValueError("Candidate name must be a non-empty string")
    selector = str(raw["selector"])
    if selector not in base.ALLOWED_SELECTORS:
        raise ValueError(f"selector must be one of {sorted(base.ALLOWED_SELECTORS)}")
    n_genes = int(raw["n_genes"])
    if not 4 <= n_genes <= int(BUDGET["max_selected_genes"]):
        raise ValueError("n_genes is outside the fixed arena budget")
    representation = str(raw["representation"])
    if representation not in {"raw", "module"}:
        raise ValueError("representation must be raw or module")
    module_count = int(raw["module_count"])
    if representation == "raw" and module_count != n_genes:
        raise ValueError("raw representation requires module_count == n_genes")
    if representation == "module" and not 1 <= module_count <= min(4, n_genes):
        raise ValueError("module representation requires 1-4 nonempty modules")
    threshold = float(raw["benefit_threshold_months"])
    if not 0.0 <= threshold <= 3.0:
        raise ValueError("benefit_threshold_months must be between 0 and 3")
    learner = raw["tlearner"]
    if not isinstance(learner, dict) or set(learner) != {"observation", "act"}:
        raise ValueError("tlearner must contain exactly observation and act")
    added = n_genes if representation == "raw" else module_count
    if len(PRETREATMENT_COLUMNS) + added > int(BUDGET["max_total_features"]):
        raise ValueError("Candidate exceeds max_total_features")
    return CandidateSpec(
        name=raw["name"].strip(), selector=selector, n_genes=n_genes,
        representation=representation, module_count=module_count,
        benefit_threshold_months=threshold,
        observation=_validate_arm(learner["observation"], "observation"),
        act=_validate_arm(learner["act"], "act"),
    )


def select_genes_by_dr_benefit(
    fit: pd.DataFrame,
    pathways: dict[str, tuple[str, ...]],
    genes: list[str],
    spec: CandidateSpec,
    *,
    smoke: bool = False,
) -> list[str]:
    """Delegate to the v2 fit-only DR selector under the same sealed inputs."""
    return base.select_genes_by_dr_benefit(fit, pathways, genes, spec, smoke=smoke)


class FeatureTransformer:
    """Fit-only pretreatment imputation and optional gene-module compression."""

    def __init__(self, selected_genes: list[str], representation: str, module_count: int) -> None:
        self.selected_genes = list(selected_genes)
        self.representation = representation
        self.module_count = int(module_count)
        self.clinical_imputer = SimpleImputer(strategy="median")
        self.gene_imputer = SimpleImputer(strategy="median") if selected_genes else None
        self.gene_scaler = StandardScaler() if representation == "module" else None
        if not selected_genes:
            self.names = list(PRETREATMENT_COLUMNS)
            self.groups: list[np.ndarray] = []
        elif representation == "raw":
            self.names = list(PRETREATMENT_COLUMNS) + list(selected_genes)
            self.groups = []
        else:
            self.names = list(PRETREATMENT_COLUMNS) + [
                f"REACTOME_MODULE_{index + 1}" for index in range(module_count)
            ]
            self.groups = [
                np.asarray(group, dtype=int)
                for group in np.array_split(np.arange(len(selected_genes)), module_count)
            ]

    def fit(self, frame: pd.DataFrame) -> "FeatureTransformer":
        self.clinical_imputer.fit(frame[PRETREATMENT_COLUMNS].to_numpy(float))
        if self.selected_genes:
            genes = self.gene_imputer.fit_transform(frame[self.selected_genes].to_numpy(float))
            if self.gene_scaler is not None:
                self.gene_scaler.fit(genes)
        return self

    def transform(self, frame: pd.DataFrame) -> np.ndarray:
        clinical = self.clinical_imputer.transform(frame[PRETREATMENT_COLUMNS].to_numpy(float))
        if not self.selected_genes:
            return clinical
        genes = self.gene_imputer.transform(frame[self.selected_genes].to_numpy(float))
        if self.representation == "raw":
            return np.column_stack([clinical, genes])
        standardized = self.gene_scaler.transform(genes)
        modules = np.column_stack([
            standardized[:, group].mean(axis=1) for group in self.groups
        ])
        return np.column_stack([clinical, modules])


def _forest(parameters: ArmForestSpec, seed: int, smoke: bool) -> RandomSurvivalForest:
    return RandomSurvivalForest(
        n_estimators=min(int(BUDGET["smoke_trees"]), parameters.n_estimators) if smoke else parameters.n_estimators,
        max_depth=parameters.max_depth,
        min_samples_leaf=parameters.min_samples_leaf,
        min_samples_split=parameters.min_samples_split,
        max_features=parameters.max_features,
        bootstrap=True,
        oob_score=False,
        n_jobs=1,
        random_state=int(seed),
        low_memory=False,
    )


def _locked_clinical_parameters(arm: str) -> ArmForestSpec:
    raw = BUDGET["clinical_tlearner"][arm]
    return ArmForestSpec(
        int(raw["n_estimators"]), int(raw["max_depth"]),
        int(raw["min_samples_leaf"]), int(raw["min_samples_split"]),
        float(raw["max_features"]),
    )


def _fit_tlearner(
    fit: pd.DataFrame,
    selected: list[str],
    spec: CandidateSpec,
    seed: int,
    *,
    clinical_only: bool,
    smoke: bool,
) -> tuple[RandomSurvivalForest, RandomSurvivalForest, FeatureTransformer, dict]:
    transformer = FeatureTransformer(
        [] if clinical_only else selected,
        "raw" if clinical_only else spec.representation,
        0 if clinical_only else spec.module_count,
    ).fit(fit)
    matrix = transformer.transform(fit)
    treatment = fit[TREATMENT].to_numpy(int)
    weights = _fit_training_weights(fit)
    models: list[RandomSurvivalForest] = []
    support: dict[str, dict] = {}
    for arm, label in ((0, "observation"), (1, "act")):
        mask = treatment == arm
        events = int(fit.loc[mask, "OS_STATUS"].sum())
        if int(mask.sum()) < 2 or events < 1:
            raise RuntimeError(f"T-learner {label} arm has insufficient fold support")
        params = (
            _locked_clinical_parameters(label)
            if clinical_only else (spec.observation if arm == 0 else spec.act)
        )
        model = _forest(params, seed + arm * 100_003, smoke)
        model.fit(matrix[mask], _outcome(fit.loc[mask]), sample_weight=weights[mask])
        models.append(model)
        support[label] = {"patients": int(mask.sum()), "events": events}
    if TREATMENT in transformer.names:
        raise RuntimeError("ACT leaked into a T-learner feature matrix")
    return models[0], models[1], transformer, support


def _genomic_split_fraction(model: RandomSurvivalForest, first_genomic: int) -> float:
    genomic = 0
    total = 0
    for estimator in model.estimators_:
        features = estimator.tree_.feature
        used = features[features >= 0]
        total += int(used.size)
        genomic += int(np.sum(used >= first_genomic))
    return float(genomic / total) if total else 0.0


def _tlearner_predictions(
    model0: RandomSurvivalForest,
    model1: RandomSurvivalForest,
    transformer: FeatureTransformer,
    assess: pd.DataFrame,
    support: dict,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    matrix = transformer.transform(assess)
    risk0 = np.asarray(model0.predict(matrix), dtype=float)
    risk1 = np.asarray(model1.predict(matrix), dtype=float)
    treatment = assess[TREATMENT].to_numpy(int)
    risk = np.where(treatment == 1, risk1, risk0)
    rmst0 = _integrate_step_functions(model0.predict_survival_function(matrix), TAU)
    rmst1 = _integrate_step_functions(model1.predict_survival_function(matrix), TAU)
    first_genomic = len(PRETREATMENT_COLUMNS)
    usage = {
        "support": support,
        "genomic_split_fraction": {
            "observation": _genomic_split_fraction(model0, first_genomic),
            "act": _genomic_split_fraction(model1, first_genomic),
        },
    }
    return risk, rmst0, rmst1, usage


def _mean_pairwise_correlation(rows: list[np.ndarray]) -> float:
    values: list[float] = []
    for left in range(len(rows)):
        for right in range(left + 1, len(rows)):
            a, b = rows[left], rows[right]
            if np.std(a) <= 1e-12 or np.std(b) <= 1e-12:
                values.append(1.0 if np.allclose(a, b) else 0.0)
            else:
                values.append(float(np.corrcoef(a, b)[0, 1]))
    return float(np.mean(values)) if values else 1.0


def _mean_nested(rows: list[dict]) -> dict:
    output: dict = {}
    for key in rows[0]:
        values = [row[key] for row in rows]
        if isinstance(values[0], dict):
            output[key] = _mean_nested(values)
        elif isinstance(values[0], str):
            if len(set(values)) != 1:
                raise RuntimeError(f"Cannot aggregate inconsistent diagnostic {key}")
            output[key] = values[0]
        elif isinstance(values[0], (bool, np.bool_)):
            output[key] = bool(all(values))
        else:
            output[key] = float(np.mean(values))
    return output


def _combine_usage(rows: list[dict]) -> dict:
    return _mean_nested(rows)


def _empty_predictions(n: int) -> dict:
    return {
        name: {
            "risk": np.zeros(n), "rmst0": np.zeros(n), "rmst1": np.zeros(n),
            "seed_agreement": [], "seed_benefit_correlation": [], "usage": [],
        }
        for name in ("clinical", "clinical_plus_genomic")
    }


def _validate_selected(selected: list[str], genes: list[str], spec: CandidateSpec) -> None:
    if len(selected) != spec.n_genes or len(set(selected)) != len(selected):
        raise ValueError("Gene selector must return exactly n_genes unique symbols")
    if any(gene not in genes for gene in selected):
        raise ValueError("Gene selector returned a non-Reactome/non-development gene")


def _fit_predict_fold(
    fit: pd.DataFrame,
    assess: pd.DataFrame,
    selected: list[str],
    spec: CandidateSpec,
    seeds: tuple[int, ...],
    *,
    smoke: bool,
) -> tuple[dict, np.ndarray, np.ndarray, np.ndarray]:
    predictions = _empty_predictions(len(assess))
    for model_name, clinical_only in (("clinical", True), ("clinical_plus_genomic", False)):
        threshold = 0.0 if clinical_only else spec.benefit_threshold_months
        seed_risk, seed_rmst0, seed_rmst1, seed_rec, seed_benefit, seed_usage = [], [], [], [], [], []
        for seed in seeds:
            model0, model1, transformer, support = _fit_tlearner(
                fit, selected, spec, seed, clinical_only=clinical_only, smoke=smoke
            )
            risk, rmst0, rmst1, usage = _tlearner_predictions(
                model0, model1, transformer, assess, support
            )
            benefit = rmst1 - rmst0
            seed_risk.append(risk)
            seed_rmst0.append(rmst0)
            seed_rmst1.append(rmst1)
            seed_benefit.append(benefit)
            seed_rec.append(recommendations(benefit, threshold))
            seed_usage.append(usage)
        panel = predictions[model_name]
        panel["risk"] = np.mean(seed_risk, axis=0)
        panel["rmst0"] = np.mean(seed_rmst0, axis=0)
        panel["rmst1"] = np.mean(seed_rmst1, axis=0)
        ensemble_rec = recommendations(panel["rmst1"] - panel["rmst0"], threshold)
        panel["seed_agreement"].append(float(np.mean(np.asarray(seed_rec) == ensemble_rec[None, :])))
        panel["seed_benefit_correlation"].append(_mean_pairwise_correlation(seed_benefit))
        panel["usage"].append(_combine_usage(seed_usage))
    raw_propensity, propensity, _ = _fit_propensity(fit, assess)
    return predictions, raw_propensity, propensity, _ipcw_restricted_time(fit, assess)


def _summarize_cohort(
    frame: pd.DataFrame,
    predictions: dict,
    raw_propensity: np.ndarray,
    propensity: np.ndarray,
    ipcw_time: np.ndarray,
    threshold: float,
) -> tuple[dict, dict, dict]:
    treatment = frame[TREATMENT].to_numpy(int)
    model_results: dict[str, dict] = {}
    phi: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for name in ("clinical", "clinical_plus_genomic"):
        panel = predictions[name]
        phi0, phi1 = aipw_arm_scores(treatment, propensity, ipcw_time, panel["rmst0"], panel["rmst1"])
        benefit = panel["rmst1"] - panel["rmst0"]
        model_threshold = 0.0 if name == "clinical" else threshold
        rec = recommendations(benefit, model_threshold)
        phi[name] = (phi0, phi1, rec)
        usage = _combine_usage(panel["usage"])
        model_results[name] = {
            **policy_summary(phi0, phi1, rec),
            "harrell_cindex": float(concordance_index_censored(
                frame["OS_STATUS"].astype(bool).to_numpy(),
                frame["OS_MONTHS"].to_numpy(float), panel["risk"],
            )[0]),
            "mean_predicted_benefit_months": float(np.mean(benefit)),
            "predicted_benefit_iqr_months": float(np.percentile(benefit, 75) - np.percentile(benefit, 25)),
            "median_absolute_predicted_benefit_months": float(np.median(np.abs(benefit))),
            "nontrivial_benefit_fraction": float(np.mean(np.abs(benefit) > max(model_threshold, 0.10))),
            "seed_agreement": float(np.mean(panel["seed_agreement"])),
            "seed_benefit_correlation": float(np.mean(panel["seed_benefit_correlation"])),
            "act_mechanism": {
                "role": "arm assignment only; ACT absent from X",
                "treatment_absent_from_features": True,
            },
            "arm_support": usage["support"],
            "genomic_split_fraction": usage["genomic_split_fraction"],
        }
    treatment_weight = np.where(treatment == 1, 1 / propensity, 1 / (1 - propensity))
    overlap = {
        "raw_propensity_in_0_05_0_95_fraction": float(np.mean(
            (raw_propensity >= PROPENSITY_CLIP[0]) & (raw_propensity <= PROPENSITY_CLIP[1])
        )),
        "iptw_effective_sample_size": float(treatment_weight.sum() ** 2 / np.square(treatment_weight).sum()),
    }
    return model_results, phi, overlap


GeneSelector = Callable[[pd.DataFrame, dict[str, tuple[str, ...]], list[str], CandidateSpec], list[str]]


def _crossfit_repeat(
    frame: pd.DataFrame,
    pathways: dict[str, tuple[str, ...]],
    genes: list[str],
    spec: CandidateSpec,
    selector_fn: GeneSelector,
    seeds: tuple[int, ...],
    n_folds: int,
    repeat: int,
    draws: int,
    *,
    smoke: bool,
) -> tuple[dict, list[list[str]]]:
    predictions = _empty_predictions(len(frame))
    raw_propensity = np.zeros(len(frame))
    propensity = np.zeros(len(frame))
    ipcw_time = np.zeros(len(frame))
    selected_lists: list[list[str]] = []
    fold_seed = FOLD_SEED + 10_000 * repeat
    for fold, (fit_index, assess_index) in enumerate(fixed_folds(frame, n_folds, seed=fold_seed), start=1):
        fit = frame.iloc[fit_index].reset_index(drop=True)
        assess = frame.iloc[assess_index].reset_index(drop=True)
        selected = selector_fn(fit, pathways, genes, spec)
        _validate_selected(selected, genes, spec)
        selected_lists.append(list(selected))
        fold_predictions, raw_e, e, fold_ipcw = _fit_predict_fold(
            fit, assess, selected, spec, seeds, smoke=smoke
        )
        for model_name in ("clinical", "clinical_plus_genomic"):
            for metric in ("risk", "rmst0", "rmst1"):
                predictions[model_name][metric][assess_index] = fold_predictions[model_name][metric]
            for metric in ("seed_agreement", "seed_benefit_correlation", "usage"):
                predictions[model_name][metric].extend(fold_predictions[model_name][metric])
        raw_propensity[assess_index], propensity[assess_index], ipcw_time[assess_index] = raw_e, e, fold_ipcw
        print(f"[repeat {repeat} fold {fold}/{n_folds}] selected={','.join(selected[:5])}...", flush=True)
    models, phi, overlap = _summarize_cohort(
        frame, predictions, raw_propensity, propensity, ipcw_time,
        spec.benefit_threshold_months,
    )
    c0, c1, rec_c = phi["clinical"]
    g0, g1, rec_g = phi["clinical_plus_genomic"]
    increment = models["clinical_plus_genomic"]["alignment_months"] - models["clinical"]["alignment_months"]
    return {
        "repeat": int(repeat),
        "fold_seed": int(fold_seed),
        "incremental_alignment_months": float(increment),
        "paired_bootstrap": _bootstrap_increment(c0, c1, rec_c, g0, g1, rec_g, draws=draws),
        "models": models,
        "constant_policy_values": {
            "all_observation_months": float(np.mean(g0)),
            "all_act_months": float(np.mean(g1)),
        },
        "overlap": overlap,
    }, selected_lists


def evaluate_candidate(
    raw_candidate: dict,
    *,
    selector: GeneSelector | None = None,
    smoke: bool = False,
) -> dict:
    spec = validate_candidate(raw_candidate)
    frame, pathways, genes = load_development()
    n_folds = 2 if smoke else int(BUDGET["outer_folds"])
    n_repeats = 1 if smoke else int(BUDGET["outer_repeats"])
    seeds = tuple(BUDGET["forest_seeds"][:1] if smoke else BUDGET["forest_seeds"])
    draws = 100 if smoke else int(BUDGET["patient_bootstraps"])
    selector_fn = selector or (
        lambda fit, pathway_map, available, candidate: select_genes_by_dr_benefit(
            fit, pathway_map, available, candidate, smoke=smoke
        )
    )
    repeat_results: list[dict] = []
    selections_by_repeat: list[list[list[str]]] = []
    for repeat in range(1, n_repeats + 1):
        result, selections = _crossfit_repeat(
            frame, pathways, genes, spec, selector_fn, seeds, n_folds,
            repeat, draws, smoke=smoke,
        )
        repeat_results.append(result)
        selections_by_repeat.append(selections)
    selected_lists = [item for rows in selections_by_repeat for item in rows]
    selected_full = selector_fn(frame, pathways, genes, spec)
    _validate_selected(selected_full, genes, spec)
    selected_modules = (
        [[selected_full[int(index)] for index in group]
         for group in np.array_split(np.arange(len(selected_full)), spec.module_count)]
        if spec.representation == "module" else [[gene] for gene in selected_full]
    )
    development_models = _mean_nested([item["models"] for item in repeat_results])
    development_overlap = _mean_nested([item["overlap"] for item in repeat_results])
    development_constants = _mean_nested([item["constant_policy_values"] for item in repeat_results])
    increments = np.asarray([item["incremental_alignment_months"] for item in repeat_results])
    increment_range = float(np.ptp(increments)) if len(increments) > 1 else 0.0
    jaccard = _pairwise_jaccard(selected_lists)
    gates = {
        "all_repeat_genomic_increment_positive": bool(np.all(increments > 0)),
        "all_repeat_genomic_alignment_positive": all(item["models"]["clinical_plus_genomic"]["alignment_months"] > 0 for item in repeat_results),
        "all_repeat_genomic_value_at_least_clinical": all(item["models"]["clinical_plus_genomic"]["value_months"] >= item["models"]["clinical"]["value_months"] for item in repeat_results),
        "all_repeat_genomic_value_at_least_best_constant": all(item["models"]["clinical_plus_genomic"]["value_months"] >= max(item["constant_policy_values"].values()) for item in repeat_results),
        "all_repeat_cindex_drop_no_more_than_0_03": all(item["models"]["clinical_plus_genomic"]["harrell_cindex"] >= item["models"]["clinical"]["harrell_cindex"] - 0.03 for item in repeat_results),
        "all_repeat_genomic_seed_agreement_at_least_0_85": all(item["models"]["clinical_plus_genomic"]["seed_agreement"] >= 0.85 for item in repeat_results),
        "all_repeat_genomic_seed_benefit_correlation_at_least_0_50": all(item["models"]["clinical_plus_genomic"]["seed_benefit_correlation"] >= 0.50 for item in repeat_results),
        "all_repeat_treatment_absent_from_features": all(item["models"]["clinical_plus_genomic"]["act_mechanism"]["treatment_absent_from_features"] for item in repeat_results),
        "all_repeat_mean_genomic_split_fraction_at_least_0_01": all(np.mean(list(item["models"]["clinical_plus_genomic"]["genomic_split_fraction"].values())) >= 0.01 for item in repeat_results),
        "gene_selection_jaccard_at_least_0_10": jaccard >= 0.10,
        "all_repeat_nontrivial_benefit_fraction_at_least_0_10": all(item["models"]["clinical_plus_genomic"]["nontrivial_benefit_fraction"] >= 0.10 for item in repeat_results),
        "all_repeat_raw_propensity_overlap_at_least_0_80": all(item["overlap"]["raw_propensity_in_0_05_0_95_fraction"] >= 0.80 for item in repeat_results),
        "all_repeat_iptw_effective_sample_size_at_least_0_30n": all(item["overlap"]["iptw_effective_sample_size"] >= 0.30 * len(frame) for item in repeat_results),
    }
    eligible = all(gates.values())
    robust_lcb = min(item["paired_bootstrap"]["selection_lcb"] for item in repeat_results)
    score = float(robust_lcb - increment_range)
    counts: dict[str, int] = {}
    for selected in selected_lists:
        for gene in selected:
            counts[gene] = counts.get(gene, 0) + 1
    return {
        "schema_version": 1,
        "phase": "smoke" if smoke else "pooled_development_repeated_crossfit",
        "candidate": spec.as_dict(),
        "objective": "worst repeated-development OOF selection LCB minus repeat range for A60(C+G)-A60(C)",
        "estimator": "RSF T-learner policies graded by outer-cross-fitted IPCW-AIPW restricted-time value",
        "reward": score if eligible else -1_000_000.0,
        "eligible": bool(eligible),
        "gates": gates,
        "generalization": {
            "robust_selection_lcb": float(robust_lcb),
            "repeat_increment_range_months": increment_range,
            "score_before_eligibility": score,
        },
        "development_oof": {
            "incremental_alignment_months": float(np.mean(increments)),
            "models": development_models,
            "constant_policy_values": development_constants,
            "overlap": development_overlap,
            "repeat_increment_min_months": float(np.min(increments)),
            "repeat_increment_max_months": float(np.max(increments)),
        },
        "repeat_results": repeat_results,
        "gene_selection": {
            "mean_pairwise_jaccard": jaccard,
            "per_repeat_outer_fold": selections_by_repeat,
            "full_development_selection": selected_full,
            "full_development_modules": selected_modules,
            "frequency": [
                {"gene": gene, "count": count}
                for gene, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
            ],
        },
        "data": {
            "development_n": int(len(frame)),
            "development_events": int(frame["OS_STATUS"].sum()),
            "development_act": int(frame[TREATMENT].sum()),
            "outer_folds": n_folds,
            "outer_repeats": n_repeats,
            "forest_seeds": list(seeds),
            "former_validation_used_as_development": True,
            "confirmatory_validation_used": False,
            "test_used": False,
            "train_sha256": sha256_file(TRAIN_CSV),
            "validation_sha256": sha256_file(VALID_CSV),
            "reactome_collection_sha256": sha256_file(base.REACTOME_DATA_DIR / base.MSIGDB_FILENAME),
        },
        "budget": {**BUDGET, "smoke_reductions_applied": bool(smoke)},
        "limitations": [
            "The ACT-arm learners have substantially less training support than the observation-arm learners.",
            "Search uncertainty bootstraps fixed OOF predictions and does not refit the full pipeline in each draw.",
            "Repeated cross-fits reuse patients and are stability views, not independent validation cohorts.",
            "IPCW-AIPW validity requires adequate measured-confounder adjustment, positivity, and conditionally independent censoring.",
            "Random row-level cross-validation does not establish transportability across studies or expression platforms.",
        ],
    }

