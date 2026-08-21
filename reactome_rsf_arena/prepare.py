"""Frozen data, feature-selection, fitting, and policy-value utilities.

The autonomous search calls :func:`evaluate_candidate`, which may load the
training and validation splits but never a test split.  Training predictions
are out-of-fold; validation predictions come from a full-training refit.  The
primary objective rewards the worse conservative IPCW-AIPW 60-month policy
increment across training and validation and subtracts their generalization
gap.

This module is part of the arena's locked surface.  Agents edit ``train.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression, Ridge
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
PRETREATMENT_COLUMNS = [column for column in CLINICAL_COLUMNS if column != TREATMENT]

# Reference levels are omitted from nuisance Cox/logistic models.  The RSFs
# retain the repository's complete prespecified clinical feature set.
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

ALLOWED_SELECTORS = {"dr_gene", "dr_pathway", "dr_hybrid"}
FOLD_SEED = 20260821
BOOTSTRAP_SEED = 20260822
PROPENSITY_CLIP = (0.05, 0.95)
CENSORING_SURVIVAL_FLOOR = 0.05
WEIGHT_CLIP = (0.1, 5.0)
NUMERICAL_BENEFIT_TOLERANCE_MONTHS = 1e-6


def load_budget() -> dict:
    return json.loads(BUDGET_PATH.read_text(encoding="utf-8"))


BUDGET = load_budget()
TAU = float(BUDGET["tau_months"])


@dataclass(frozen=True)
class CandidateSpec:
    name: str
    selector: str
    n_genes: int
    benefit_threshold_months: float
    n_estimators: int
    max_depth: int
    min_samples_leaf: int
    min_samples_split: int
    max_features: float

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "selector": self.selector,
            "n_genes": self.n_genes,
            "benefit_threshold_months": self.benefit_threshold_months,
            "rsf": {
                "n_estimators": self.n_estimators,
                "max_depth": self.max_depth,
                "min_samples_leaf": self.min_samples_leaf,
                "min_samples_split": self.min_samples_split,
                "max_features": self.max_features,
            },
        }


def validate_candidate(raw: dict) -> CandidateSpec:
    if set(raw) != {
        "name", "selector", "n_genes", "benefit_threshold_months", "rsf",
    }:
        raise ValueError("Candidate has missing or unrecognized top-level keys")
    if not isinstance(raw["name"], str) or not raw["name"].strip():
        raise ValueError("Candidate name must be a non-empty string")
    selector = str(raw["selector"])
    if selector not in ALLOWED_SELECTORS:
        raise ValueError(f"selector must be one of {sorted(ALLOWED_SELECTORS)}")
    n_genes = int(raw["n_genes"])
    if not 4 <= n_genes <= int(BUDGET["max_selected_genes"]):
        raise ValueError("n_genes is outside the fixed arena budget")
    threshold = float(raw["benefit_threshold_months"])
    if not 0.0 <= threshold <= 3.0:
        raise ValueError("benefit_threshold_months must be between 0 and 3")
    rsf = raw["rsf"]
    if set(rsf) != {
        "n_estimators", "max_depth", "min_samples_leaf",
        "min_samples_split", "max_features",
    }:
        raise ValueError("Candidate has missing or unrecognized RSF keys")
    trees = int(rsf["n_estimators"])
    depth = int(rsf["max_depth"])
    leaf = int(rsf["min_samples_leaf"])
    split = int(rsf["min_samples_split"])
    mtry = float(rsf["max_features"])
    if not 100 <= trees <= int(BUDGET["max_trees_per_forest"]):
        raise ValueError("n_estimators is outside the fixed arena budget")
    if not 3 <= depth <= 12:
        raise ValueError("max_depth must be in [3, 12]")
    if not 8 <= leaf <= 40:
        raise ValueError("min_samples_leaf must be in [8, 40]")
    if split < 2 * leaf or split > 100:
        raise ValueError("min_samples_split must be >= 2*min_samples_leaf and <= 100")
    if not 0.20 <= mtry <= 1.0:
        raise ValueError("max_features must be in [0.20, 1.0]")
    total_features = len(CLINICAL_COLUMNS) + n_genes
    if total_features > int(BUDGET["max_total_features"]):
        raise ValueError("Candidate exceeds max_total_features")
    return CandidateSpec(
        name=raw["name"].strip(), selector=selector, n_genes=n_genes,
        benefit_threshold_months=threshold, n_estimators=trees,
        max_depth=depth, min_samples_leaf=leaf,
        min_samples_split=split, max_features=mtry,
    )


def _assert_allowed_search_path(path: str | Path) -> None:
    resolved = Path(path).resolve()
    allowed = {TRAIN_CSV.resolve(), VALID_CSV.resolve()}
    if resolved not in allowed:
        raise ValueError(f"Search loader accepts only train and validation CSVs: {resolved}")
    if "test" in resolved.name.lower():
        raise ValueError(f"Test paths are forbidden during search: {resolved}")


def _prepare_frame(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    mapped = result[TREATMENT].map({"OBS": 0, "ACT": 1})
    result[TREATMENT] = mapped.where(mapped.notna(), result[TREATMENT])
    for column in result.columns:
        result[column] = pd.to_numeric(result[column], errors="coerce")
    if result[OUTCOME_COLUMNS + CLINICAL_COLUMNS].isna().any().any():
        bad = result[OUTCOME_COLUMNS + CLINICAL_COLUMNS].columns[
            result[OUTCOME_COLUMNS + CLINICAL_COLUMNS].isna().any()
        ].tolist()
        raise ValueError(f"Missing required analysis values: {bad}")
    result[TREATMENT] = result[TREATMENT].astype(int)
    result["OS_STATUS"] = result["OS_STATUS"].astype(int)
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


def load_train_valid(
    train_csv: str | Path = TRAIN_CSV,
    valid_csv: str | Path = VALID_CSV,
    data_dir: str | Path = REACTOME_DATA_DIR,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, tuple[str, ...]], list[str]]:
    _assert_allowed_search_path(train_csv)
    _assert_allowed_search_path(valid_csv)
    pathways = load_reactome_sets(data_dir)
    train_header = set(pd.read_csv(train_csv, nrows=0).columns)
    valid_header = set(pd.read_csv(valid_csv, nrows=0).columns)
    all_reactome = {gene for members in pathways.values() for gene in members}
    genes = sorted(all_reactome.intersection(train_header).intersection(valid_header))
    requested = OUTCOME_COLUMNS + CLINICAL_COLUMNS + genes
    train = _prepare_frame(pd.read_csv(train_csv, usecols=requested))
    valid = _prepare_frame(pd.read_csv(valid_csv, usecols=requested))
    values = train[genes].to_numpy(dtype=np.float32)
    with np.errstate(all="ignore"):
        spans = np.nanmax(values, axis=0) - np.nanmin(values, axis=0)
    genes = [gene for gene, span in zip(genes, spans) if np.isfinite(span) and span > 0]
    pathways = {
        name: tuple(gene for gene in members if gene in set(genes))
        for name, members in pathways.items()
    }
    pathways = {name: members for name, members in pathways.items() if members}
    columns = OUTCOME_COLUMNS + CLINICAL_COLUMNS + genes
    return train[columns], valid[columns], pathways, genes


def fixed_folds(frame: pd.DataFrame, n_splits: int, seed: int = FOLD_SEED):
    strata = frame["OS_STATUS"].astype(str) + "_" + frame[TREATMENT].astype(str)
    splitter = StratifiedKFold(n_splits=int(n_splits), shuffle=True, random_state=int(seed))
    yield from splitter.split(np.zeros(len(frame)), strata)


def _outcome(frame: pd.DataFrame, *, censoring: bool = False) -> np.ndarray:
    event = 1 - frame["OS_STATUS"].to_numpy(int) if censoring else frame["OS_STATUS"].to_numpy(int)
    return Surv.from_arrays(event=event.astype(bool), time=frame["OS_MONTHS"].to_numpy(float))


def _nuisance_matrix(frame: pd.DataFrame, *, treatment: np.ndarray | None = None, interactions: bool = False) -> np.ndarray:
    x = frame[NUISANCE_COLUMNS].to_numpy(dtype=float)
    a = frame[TREATMENT].to_numpy(dtype=float) if treatment is None else np.asarray(treatment, dtype=float)
    blocks = [a[:, None], x]
    if interactions:
        blocks.append(x * a[:, None])
    return np.column_stack(blocks)


def _fit_propensity(fit: pd.DataFrame, assess: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    model = make_pipeline(
        StandardScaler(), LogisticRegression(max_iter=2000, solver="lbfgs", C=0.5)
    )
    model.fit(fit[NUISANCE_COLUMNS].to_numpy(float), fit[TREATMENT].to_numpy(int))
    raw = model.predict_proba(assess[NUISANCE_COLUMNS].to_numpy(float))[:, 1]
    clipped = np.clip(raw, *PROPENSITY_CLIP)
    treatment = assess[TREATMENT].to_numpy(int)
    prevalence = float(fit[TREATMENT].mean())
    stabilized = np.where(treatment == 1, prevalence / clipped, (1 - prevalence) / (1 - clipped))
    return raw, clipped, np.clip(stabilized, *WEIGHT_CLIP)


def _fit_training_weights(fit: pd.DataFrame) -> np.ndarray:
    _, clipped, _ = _fit_propensity(fit, fit)
    treatment = fit[TREATMENT].to_numpy(int)
    prevalence = float(treatment.mean())
    weights = np.where(treatment == 1, prevalence / clipped, (1 - prevalence) / (1 - clipped))
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


def _clinical_outcome_predictions(fit: pd.DataFrame, assess: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    model = make_pipeline(
        StandardScaler(), CoxPHSurvivalAnalysis(alpha=1.0, ties="efron")
    )
    model.fit(_nuisance_matrix(fit, interactions=True), _outcome(fit))
    zeros = np.zeros(len(assess))
    ones = np.ones(len(assess))
    sf0 = model.predict_survival_function(_nuisance_matrix(assess, treatment=zeros, interactions=True))
    sf1 = model.predict_survival_function(_nuisance_matrix(assess, treatment=ones, interactions=True))
    return _integrate_step_functions(sf0), _integrate_step_functions(sf1)


def _ipcw_restricted_time(fit: pd.DataFrame, assess: pd.DataFrame, tau: float = TAU) -> np.ndarray:
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
        censor_survival = np.maximum(censor_survival, CENSORING_SURVIVAL_FLOOR)
        result[index] = float(np.sum(np.diff(bounds) / censor_survival))
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
    ipcw_time = np.asarray(ipcw_time, dtype=float)
    mu0 = np.asarray(mu0, dtype=float)
    mu1 = np.asarray(mu1, dtype=float)
    phi1 = mu1 + treatment / propensity * (ipcw_time - mu1)
    phi0 = mu0 + (1 - treatment) / (1 - propensity) * (ipcw_time - mu0)
    return phi0, phi1


def cross_fitted_benefit_pseudo_outcome(frame: pd.DataFrame, n_splits: int) -> np.ndarray:
    gamma = np.zeros(len(frame), dtype=float)
    for fit_index, assess_index in fixed_folds(frame, n_splits, seed=FOLD_SEED + 101):
        fit = frame.iloc[fit_index].reset_index(drop=True)
        assess = frame.iloc[assess_index].reset_index(drop=True)
        _, propensity, _ = _fit_propensity(fit, assess)
        ipcw_time = _ipcw_restricted_time(fit, assess)
        mu0, mu1 = _clinical_outcome_predictions(fit, assess)
        phi0, phi1 = aipw_arm_scores(
            assess[TREATMENT].to_numpy(int), propensity, ipcw_time, mu0, mu1
        )
        gamma[assess_index] = phi1 - phi0
    return gamma


def _gene_effect_scores(frame: pd.DataFrame, genes: list[str], gamma: np.ndarray) -> dict[str, float]:
    clinical = frame[NUISANCE_COLUMNS].to_numpy(dtype=float)
    clinical = StandardScaler().fit_transform(clinical)
    design = np.column_stack([np.ones(len(frame)), clinical])
    gamma_residual = gamma - Ridge(alpha=1.0).fit(clinical, gamma).predict(clinical)
    gamma_norm = max(float(np.linalg.norm(gamma_residual)), 1e-12)
    scores: dict[str, float] = {}
    chunk_size = 512
    for start in range(0, len(genes), chunk_size):
        chunk = genes[start:start + chunk_size]
        values = frame[chunk].to_numpy(dtype=float)
        medians = np.nanmedian(values, axis=0)
        bad = ~np.isfinite(values)
        if bad.any():
            values[bad] = np.take(medians, np.where(bad)[1])
        coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
        residual = values - design @ coefficients
        denominators = np.linalg.norm(residual, axis=0) * gamma_norm
        numerators = np.abs(gamma_residual @ residual)
        correlations = np.divide(numerators, denominators, out=np.zeros_like(numerators), where=denominators > 0)
        scores.update({gene: float(score) for gene, score in zip(chunk, correlations)})
    return scores


def select_genes_by_dr_benefit(
    fit: pd.DataFrame,
    pathways: dict[str, tuple[str, ...]],
    genes: list[str],
    spec: CandidateSpec,
    *,
    smoke: bool = False,
) -> list[str]:
    available = genes[:160] if smoke else genes
    inner_folds = 2 if smoke else int(BUDGET["inner_folds"])
    gamma = cross_fitted_benefit_pseudo_outcome(fit, inner_folds)
    scores = _gene_effect_scores(fit, available, gamma)
    direct = sorted(available, key=lambda gene: (-scores[gene], gene))
    if spec.selector == "dr_gene":
        return direct[: spec.n_genes]

    available_set = set(available)
    pathway_rank = []
    for name, members in pathways.items():
        present = [gene for gene in members if gene in available_set]
        if not present:
            continue
        top = sorted((scores[gene] for gene in present), reverse=True)[: min(3, len(present))]
        pathway_rank.append((float(np.mean(top)), name, present))
    pathway_rank.sort(key=lambda item: (-item[0], item[1]))
    pathway_genes: list[str] = []
    for _, _, members in pathway_rank:
        for gene in sorted(members, key=lambda item: (-scores[item], item)):
            if gene not in pathway_genes:
                pathway_genes.append(gene)
                break
        if len(pathway_genes) >= spec.n_genes:
            break
    if spec.selector == "dr_pathway":
        selected = pathway_genes
    else:
        selected = []
        half = (spec.n_genes + 1) // 2
        for gene in direct[:half] + pathway_genes:
            if gene not in selected:
                selected.append(gene)
    for gene in direct:
        if len(selected) >= spec.n_genes:
            break
        if gene not in selected:
            selected.append(gene)
    return selected[: spec.n_genes]


def feature_names(selected_genes: list[str]) -> tuple[list[str], list[str]]:
    clinical = list(CLINICAL_COLUMNS)
    genomic = clinical + selected_genes
    return clinical, genomic


def build_matrix(frame: pd.DataFrame, names: list[str]) -> np.ndarray:
    return frame[names].to_numpy(dtype=float)


def _fit_rsf(
    fit: pd.DataFrame,
    names: list[str],
    spec: CandidateSpec,
    seed: int,
    *,
    smoke: bool,
) -> tuple[RandomSurvivalForest, SimpleImputer]:
    imputer = SimpleImputer(strategy="median")
    matrix = imputer.fit_transform(build_matrix(fit, names))
    model = RandomSurvivalForest(
        n_estimators=min(40, spec.n_estimators) if smoke else spec.n_estimators,
        max_depth=spec.max_depth,
        min_samples_leaf=spec.min_samples_leaf,
        min_samples_split=spec.min_samples_split,
        max_features=spec.max_features,
        bootstrap=True,
        oob_score=False,
        n_jobs=1,
        random_state=int(seed),
        low_memory=False,
    )
    model.fit(matrix, _outcome(fit), sample_weight=_fit_training_weights(fit))
    return model, imputer


def _rsf_predictions(
    model: RandomSurvivalForest,
    imputer: SimpleImputer,
    assess: pd.DataFrame,
    names: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    observed = imputer.transform(build_matrix(assess, names))
    untreated = assess.copy()
    untreated[TREATMENT] = 0
    treated = assess.copy()
    treated[TREATMENT] = 1
    x0 = imputer.transform(build_matrix(untreated, names))
    x1 = imputer.transform(build_matrix(treated, names))
    risk = np.asarray(model.predict(observed), dtype=float)
    rmst0 = _integrate_step_functions(model.predict_survival_function(x0))
    rmst1 = _integrate_step_functions(model.predict_survival_function(x1))
    act_usage = _act_usage_diagnostics(model, x0, x1, names)
    return risk, rmst0, rmst1, act_usage


def _act_usage_diagnostics(
    model: RandomSurvivalForest,
    x0: np.ndarray,
    x1: np.ndarray,
    names: list[str],
) -> dict:
    """Measure whether fitted trees operationally use the ACT feature.

    Path traversal is averaged over both members of each patient's ACT=0/ACT=1
    counterfactual pair. Terminal-node difference is patient specific and
    compares the two counterfactual copies in each tree.
    """
    act_feature = names.index(TREATMENT)
    n_patients = x0.shape[0]
    path_counts = np.zeros(n_patients, dtype=float)
    terminal_difference_counts = np.zeros(n_patients, dtype=float)
    trees_containing_act = 0
    # Direct tree traversal follows scikit-learn's float32 tree API even though
    # forest-level prediction accepts the imputer's float64 output.
    counterfactual_pair = np.asarray(np.vstack([x0, x1]), dtype=np.float32)

    for estimator in model.estimators_:
        act_nodes = np.flatnonzero(estimator.tree_.feature == act_feature)
        if act_nodes.size:
            trees_containing_act += 1
            path = estimator.decision_path(counterfactual_pair)
            traversed = np.asarray(path[:, act_nodes].sum(axis=1)).ravel() > 0
            path_counts += (
                traversed[:n_patients].astype(float)
                + traversed[n_patients:].astype(float)
            )
        leaves = estimator.apply(counterfactual_pair)
        terminal_difference_counts += leaves[:n_patients] != leaves[n_patients:]

    n_trees = len(model.estimators_)
    if n_trees == 0:
        raise RuntimeError("ACT diagnostics require at least one fitted tree")
    return {
        "tree_contains_act_split_fraction": float(trees_containing_act / n_trees),
        "patient_tree_path_act_fraction": path_counts / (2 * n_trees),
        "patient_terminal_difference_tree_fraction": terminal_difference_counts / n_trees,
    }


def _combine_seed_act_usage(seed_diagnostics: list[dict]) -> dict:
    if not seed_diagnostics:
        raise ValueError("At least one seed diagnostic is required")
    return {
        "tree_contains_act_split_fraction": float(np.mean([
            item["tree_contains_act_split_fraction"] for item in seed_diagnostics
        ])),
        "patient_tree_path_act_fraction": np.mean([
            item["patient_tree_path_act_fraction"] for item in seed_diagnostics
        ], axis=0),
        "patient_terminal_difference_tree_fraction": np.mean([
            item["patient_terminal_difference_tree_fraction"] for item in seed_diagnostics
        ], axis=0),
    }


def _summarize_act_usage(fold_diagnostics: list[dict]) -> dict:
    if not fold_diagnostics:
        raise ValueError("At least one fold diagnostic is required")
    path_fraction = np.concatenate([
        item["patient_tree_path_act_fraction"] for item in fold_diagnostics
    ])
    terminal_fraction = np.concatenate([
        item["patient_terminal_difference_tree_fraction"] for item in fold_diagnostics
    ])
    return {
        "tree_contains_act_split_fraction": float(np.mean([
            item["tree_contains_act_split_fraction"] for item in fold_diagnostics
        ])),
        "patient_tree_paths_traversing_act_fraction": float(np.mean(path_fraction)),
        "per_patient_different_terminal_tree_fraction": {
            "mean": float(np.mean(terminal_fraction)),
            "median": float(np.median(terminal_fraction)),
            "p10": float(np.percentile(terminal_fraction, 10)),
            "p90": float(np.percentile(terminal_fraction, 90)),
            "nonzero_patient_fraction": float(np.mean(terminal_fraction > 0)),
        },
    }


def recommendations(benefit: np.ndarray, threshold: float) -> np.ndarray:
    benefit = np.asarray(benefit, dtype=float)
    effective_threshold = max(float(threshold), NUMERICAL_BENEFIT_TOLERANCE_MONTHS)
    return (benefit > effective_threshold).astype(np.int8)


def _policy_value(phi0: np.ndarray, phi1: np.ndarray, policy: np.ndarray) -> float:
    policy = np.asarray(policy, dtype=np.int8)
    return float(np.mean(np.where(policy == 1, phi1, phi0)))


def policy_summary(phi0: np.ndarray, phi1: np.ndarray, policy: np.ndarray) -> dict:
    value = _policy_value(phi0, phi1, policy)
    anti_value = _policy_value(phi0, phi1, 1 - np.asarray(policy, dtype=np.int8))
    return {
        "value_months": value,
        "anti_policy_value_months": anti_value,
        "alignment_months": value - anti_value,
        "act_recommended_fraction": float(np.mean(policy)),
    }


def _pairwise_jaccard(gene_lists: list[list[str]]) -> float:
    values = []
    for left in range(len(gene_lists)):
        for right in range(left + 1, len(gene_lists)):
            a, b = set(gene_lists[left]), set(gene_lists[right])
            values.append(len(a & b) / len(a | b) if a | b else 1.0)
    return float(np.mean(values)) if values else 1.0


def _bootstrap_increment(
    phi_c0: np.ndarray,
    phi_c1: np.ndarray,
    rec_c: np.ndarray,
    phi_g0: np.ndarray,
    phi_g1: np.ndarray,
    rec_g: np.ndarray,
    *,
    draws: int,
) -> dict:
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    n = len(rec_c)
    estimates = np.empty(draws, dtype=float)
    for draw in range(draws):
        index = rng.integers(0, n, size=n)
        clinical = policy_summary(phi_c0[index], phi_c1[index], rec_c[index])["alignment_months"]
        genomic = policy_summary(phi_g0[index], phi_g1[index], rec_g[index])["alignment_months"]
        estimates[draw] = genomic - clinical
    ordinary_low, ordinary_high = np.percentile(estimates, [2.5, 97.5])
    selection_percent = 100.0 * float(BUDGET["multiplicity_alpha"]) / float(BUDGET["max_experiments"])
    return {
        "draws": int(draws),
        "mean": float(estimates.mean()),
        "sd": float(estimates.std(ddof=1)),
        "ci95_low": float(ordinary_low),
        "ci95_high": float(ordinary_high),
        "selection_lcb": float(np.percentile(estimates, selection_percent)),
        "selection_lcb_quantile_percent": selection_percent,
    }


GeneSelector = Callable[[pd.DataFrame, dict[str, tuple[str, ...]], list[str], CandidateSpec], list[str]]


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
        phi0, phi1 = aipw_arm_scores(
            treatment, propensity, ipcw_time, panel["rmst0"], panel["rmst1"]
        )
        benefit = panel["rmst1"] - panel["rmst0"]
        rec = recommendations(benefit, threshold)
        cindex = float(concordance_index_censored(
            frame["OS_STATUS"].astype(bool).to_numpy(),
            frame["OS_MONTHS"].to_numpy(float), panel["risk"],
        )[0])
        phi[name] = (phi0, phi1, rec)
        model_results[name] = {
            **policy_summary(phi0, phi1, rec),
            "harrell_cindex": cindex,
            "mean_predicted_benefit_months": float(benefit.mean()),
            "median_absolute_predicted_benefit_months": float(np.median(np.abs(benefit))),
            "nontrivial_benefit_fraction": float(
                np.mean(np.abs(benefit) > max(threshold, 0.10))
            ),
            "seed_agreement": float(np.mean(panel["seed_agreement"])),
            "act_usage": _summarize_act_usage(panel["act_usage"]),
        }
    treatment_weight = np.where(treatment == 1, 1 / propensity, 1 / (1 - propensity))
    overlap = {
        "raw_propensity_in_0_05_0_95_fraction": float(np.mean(
            (raw_propensity >= PROPENSITY_CLIP[0]) & (raw_propensity <= PROPENSITY_CLIP[1])
        )),
        "iptw_effective_sample_size": float(
            treatment_weight.sum() ** 2 / np.square(treatment_weight).sum()
        ),
    }
    return model_results, phi, overlap


def _fit_predict_holdout(
    fit: pd.DataFrame,
    assess: pd.DataFrame,
    selected: list[str],
    spec: CandidateSpec,
    seeds: tuple[int, ...],
    *,
    smoke: bool,
) -> tuple[dict, np.ndarray, np.ndarray, np.ndarray]:
    predictions = {
        name: {
            "risk": np.zeros(len(assess)), "rmst0": np.zeros(len(assess)),
            "rmst1": np.zeros(len(assess)), "seed_agreement": [], "act_usage": [],
        }
        for name in ("clinical", "clinical_plus_genomic")
    }
    clinical_names, genomic_names = feature_names(selected)
    for model_name, names in (("clinical", clinical_names), ("clinical_plus_genomic", genomic_names)):
        seed_risk, seed_rmst0, seed_rmst1, seed_rec, seed_act_usage = [], [], [], [], []
        for seed in seeds:
            model, imputer = _fit_rsf(fit, names, spec, seed, smoke=smoke)
            risk, rmst0, rmst1, act_usage = _rsf_predictions(
                model, imputer, assess, names
            )
            seed_risk.append(risk)
            seed_rmst0.append(rmst0)
            seed_rmst1.append(rmst1)
            seed_rec.append(recommendations(rmst1 - rmst0, spec.benefit_threshold_months))
            seed_act_usage.append(act_usage)
        predictions[model_name]["risk"] = np.mean(seed_risk, axis=0)
        predictions[model_name]["rmst0"] = np.mean(seed_rmst0, axis=0)
        predictions[model_name]["rmst1"] = np.mean(seed_rmst1, axis=0)
        ensemble_rec = recommendations(
            predictions[model_name]["rmst1"] - predictions[model_name]["rmst0"],
            spec.benefit_threshold_months,
        )
        predictions[model_name]["seed_agreement"].append(float(
            np.mean(np.asarray(seed_rec) == ensemble_rec[None, :])
        ))
        predictions[model_name]["act_usage"].append(
            _combine_seed_act_usage(seed_act_usage)
        )
    raw_propensity, propensity, _ = _fit_propensity(fit, assess)
    ipcw_time = _ipcw_restricted_time(fit, assess)
    return predictions, raw_propensity, propensity, ipcw_time


def evaluate_candidate(
    raw_candidate: dict,
    *,
    selector: GeneSelector | None = None,
    smoke: bool = False,
) -> dict:
    spec = validate_candidate(raw_candidate)
    frame, valid, pathways, genes = load_train_valid()
    n_folds = 2 if smoke else int(BUDGET["outer_folds"])
    seeds = tuple(BUDGET["forest_seeds"][:1] if smoke else BUDGET["forest_seeds"])
    draws = 100 if smoke else int(BUDGET["patient_bootstraps"])
    selector_fn = selector or (
        lambda fit, pathway_map, available, candidate: select_genes_by_dr_benefit(
            fit, pathway_map, available, candidate, smoke=smoke
        )
    )

    n = len(frame)
    predictions = {
        name: {
            "risk": np.zeros(n), "rmst0": np.zeros(n), "rmst1": np.zeros(n),
            "seed_agreement": [], "act_usage": [],
        }
        for name in ("clinical", "clinical_plus_genomic")
    }
    propensity = np.zeros(n)
    raw_propensity = np.zeros(n)
    ipcw_time = np.zeros(n)
    selected_lists: list[list[str]] = []

    for fold, (fit_index, assess_index) in enumerate(fixed_folds(frame, n_folds), start=1):
        fit = frame.iloc[fit_index].reset_index(drop=True)
        assess = frame.iloc[assess_index].reset_index(drop=True)
        selected = selector_fn(fit, pathways, genes, spec)
        if len(selected) != spec.n_genes or len(set(selected)) != len(selected):
            raise ValueError("Gene selector must return exactly n_genes unique symbols")
        if any(gene not in genes for gene in selected):
            raise ValueError("Gene selector returned a non-Reactome/non-training gene")
        selected_lists.append(list(selected))
        clinical_names, genomic_names = feature_names(selected)

        for model_name, names in (("clinical", clinical_names), ("clinical_plus_genomic", genomic_names)):
            seed_risk, seed_rmst0, seed_rmst1, seed_rec, seed_act_usage = [], [], [], [], []
            for seed in seeds:
                model, imputer = _fit_rsf(fit, names, spec, seed, smoke=smoke)
                risk, rmst0, rmst1, act_usage = _rsf_predictions(
                    model, imputer, assess, names
                )
                seed_risk.append(risk)
                seed_rmst0.append(rmst0)
                seed_rmst1.append(rmst1)
                seed_rec.append(recommendations(rmst1 - rmst0, spec.benefit_threshold_months))
                seed_act_usage.append(act_usage)
            ensemble_risk = np.mean(seed_risk, axis=0)
            ensemble_rmst0 = np.mean(seed_rmst0, axis=0)
            ensemble_rmst1 = np.mean(seed_rmst1, axis=0)
            ensemble_rec = recommendations(
                ensemble_rmst1 - ensemble_rmst0, spec.benefit_threshold_months
            )
            predictions[model_name]["risk"][assess_index] = ensemble_risk
            predictions[model_name]["rmst0"][assess_index] = ensemble_rmst0
            predictions[model_name]["rmst1"][assess_index] = ensemble_rmst1
            agreement = np.mean(np.asarray(seed_rec) == ensemble_rec[None, :])
            predictions[model_name]["seed_agreement"].append(float(agreement))
            predictions[model_name]["act_usage"].append(
                _combine_seed_act_usage(seed_act_usage)
            )

        raw_e, e, _ = _fit_propensity(fit, assess)
        raw_propensity[assess_index] = raw_e
        propensity[assess_index] = e
        ipcw_time[assess_index] = _ipcw_restricted_time(fit, assess)
        print(f"[fold {fold}/{n_folds}] selected={','.join(selected[:5])}...", flush=True)

    train_models, train_phi, train_overlap = _summarize_cohort(
        frame, predictions, raw_propensity, propensity, ipcw_time,
        spec.benefit_threshold_months,
    )
    tc0, tc1, train_rec_c = train_phi["clinical"]
    tg0, tg1, train_rec_g = train_phi["clinical_plus_genomic"]
    train_clinical = train_models["clinical"]
    train_genomic = train_models["clinical_plus_genomic"]
    train_increment = train_genomic["alignment_months"] - train_clinical["alignment_months"]
    train_bootstrap = _bootstrap_increment(
        tc0, tc1, train_rec_c, tg0, tg1, train_rec_g, draws=draws
    )

    # Validation is deliberately available to the search, but every transform,
    # gene choice, nuisance model, and RSF is fitted on training rows only.
    selected_full = selector_fn(frame, pathways, genes, spec)
    if len(selected_full) != spec.n_genes or len(set(selected_full)) != len(selected_full):
        raise ValueError("Full-training selector must return exactly n_genes unique symbols")
    if any(gene not in genes for gene in selected_full):
        raise ValueError("Full-training selector returned a non-Reactome/non-training gene")
    valid_predictions, valid_raw_e, valid_e, valid_ipcw = _fit_predict_holdout(
        frame, valid, selected_full, spec, seeds, smoke=smoke
    )
    valid_models, valid_phi, valid_overlap = _summarize_cohort(
        valid, valid_predictions, valid_raw_e, valid_e, valid_ipcw,
        spec.benefit_threshold_months,
    )
    vc0, vc1, valid_rec_c = valid_phi["clinical"]
    vg0, vg1, valid_rec_g = valid_phi["clinical_plus_genomic"]
    valid_clinical = valid_models["clinical"]
    valid_genomic = valid_models["clinical_plus_genomic"]
    valid_increment = valid_genomic["alignment_months"] - valid_clinical["alignment_months"]
    valid_bootstrap = _bootstrap_increment(
        vc0, vc1, valid_rec_c, vg0, vg1, valid_rec_g, draws=draws
    )

    jaccard = _pairwise_jaccard(selected_lists)
    generalization_gap = abs(train_increment - valid_increment)
    gates = {
        "train_genomic_alignment_positive": train_genomic["alignment_months"] > 0,
        "validation_genomic_alignment_positive": valid_genomic["alignment_months"] > 0,
        "train_genomic_value_at_least_clinical": train_genomic["value_months"] >= train_clinical["value_months"],
        "validation_genomic_value_at_least_clinical": valid_genomic["value_months"] >= valid_clinical["value_months"],
        "train_genomic_value_at_least_best_constant": train_genomic["value_months"] >= max(float(np.mean(tg0)), float(np.mean(tg1))),
        "validation_genomic_value_at_least_best_constant": valid_genomic["value_months"] >= max(float(np.mean(vg0)), float(np.mean(vg1))),
        "train_cindex_drop_no_more_than_0_03": train_genomic["harrell_cindex"] >= train_clinical["harrell_cindex"] - 0.03,
        "validation_cindex_drop_no_more_than_0_03": valid_genomic["harrell_cindex"] >= valid_clinical["harrell_cindex"] - 0.03,
        "train_genomic_seed_agreement_at_least_0_85": train_genomic["seed_agreement"] >= 0.85,
        "validation_genomic_seed_agreement_at_least_0_85": valid_genomic["seed_agreement"] >= 0.85,
        "gene_selection_jaccard_at_least_0_10": jaccard >= 0.10,
        "train_nontrivial_benefit_fraction_at_least_0_10": train_genomic["nontrivial_benefit_fraction"] >= 0.10,
        "validation_nontrivial_benefit_fraction_at_least_0_10": valid_genomic["nontrivial_benefit_fraction"] >= 0.10,
        "train_raw_propensity_overlap_at_least_0_80": train_overlap["raw_propensity_in_0_05_0_95_fraction"] >= 0.80,
        "validation_raw_propensity_overlap_at_least_0_80": valid_overlap["raw_propensity_in_0_05_0_95_fraction"] >= 0.80,
        "train_iptw_effective_sample_size_at_least_0_30n": train_overlap["iptw_effective_sample_size"] >= 0.30 * len(frame),
        "validation_iptw_effective_sample_size_at_least_0_30n": valid_overlap["iptw_effective_sample_size"] >= 0.30 * len(valid),
    }
    eligible = all(gates.values())
    robust_lcb = min(train_bootstrap["selection_lcb"], valid_bootstrap["selection_lcb"])
    generalization_score = robust_lcb - generalization_gap
    reward = generalization_score if eligible else -1_000_000.0
    counts: dict[str, int] = {}
    for selected in selected_lists:
        for gene in selected:
            counts[gene] = counts.get(gene, 0) + 1
    return {
        "schema_version": 1,
        "phase": "smoke" if smoke else "train_crossfit_plus_validation_search",
        "candidate": spec.as_dict(),
        "objective": (
            "min(train OOF LCB, validation LCB) minus the absolute train-validation "
            "gap for A60(C+G)-A60(C), where A60(d)=V60(d)-V60(1-d)"
        ),
        "estimator": "outer-cross-fitted IPCW-AIPW restricted-time policy value",
        "reward": float(reward),
        "eligible": bool(eligible),
        "gates": gates,
        "generalization": {
            "robust_selection_lcb": float(robust_lcb),
            "absolute_increment_gap_months": float(generalization_gap),
            "score_before_eligibility": float(generalization_score),
        },
        "train_oof": {
            "incremental_alignment_months": float(train_increment),
            "paired_bootstrap": train_bootstrap,
            "models": train_models,
            "constant_policy_values": {
                "all_observation_months": float(np.mean(tg0)),
                "all_act_months": float(np.mean(tg1)),
            },
            "overlap": train_overlap,
        },
        "validation": {
            "incremental_alignment_months": float(valid_increment),
            "paired_bootstrap": valid_bootstrap,
            "models": valid_models,
            "constant_policy_values": {
                "all_observation_months": float(np.mean(vg0)),
                "all_act_months": float(np.mean(vg1)),
            },
            "overlap": valid_overlap,
        },
        "gene_selection": {
            "mean_pairwise_jaccard": jaccard,
            "per_outer_fold": selected_lists,
            "full_training_selection": selected_full,
            "frequency": [
                {"gene": gene, "count": count}
                for gene, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
            ],
        },
        "data": {
            "train_n": int(n), "train_events": int(frame["OS_STATUS"].sum()),
            "train_act": int(frame[TREATMENT].sum()),
            "validation_n": int(len(valid)),
            "validation_events": int(valid["OS_STATUS"].sum()),
            "validation_act": int(valid[TREATMENT].sum()),
            "outer_folds": n_folds,
            "forest_seeds": list(seeds), "validation_used": True,
            "test_used": False,
            "train_sha256": sha256_file(TRAIN_CSV),
            "validation_sha256": sha256_file(VALID_CSV),
            "reactome_collection_sha256": sha256_file(REACTOME_DATA_DIR / MSIGDB_FILENAME),
        },
        "budget": {
            **BUDGET,
            "smoke_reductions_applied": bool(smoke),
        },
        "limitations": [
            "Search uncertainty uses paired bootstraps of fixed train-OOF and validation predictions; it does not refit the full pipeline in each bootstrap draw.",
            "Validation is adaptive search feedback, not untouched confirmation; the fixed experiment budget and multiplicity-adjusted LCB limit but do not eliminate validation overfitting.",
            "IPCW-AIPW validity requires adequate measured-confounder adjustment, treatment positivity, and conditionally independent censoring.",
            "Random row-level cross-validation does not establish transportability across source studies or expression platforms.",
        ],
    }
