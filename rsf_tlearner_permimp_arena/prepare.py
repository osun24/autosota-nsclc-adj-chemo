"""Locked repeated-cross-fit evaluator with fit-only, arm-specific policy PFI."""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import tempfile
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rsf_tlearner_permimp_mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored

from reactome_rsf.msigdb import sha256_file
from reactome_rsf_arena_v2 import prepare as base
from reactome_tlearner_arena import prepare as locked_t
from . import integrity


ARENA_DIR = Path(__file__).resolve().parent
REPO_ROOT = ARENA_DIR.parent
BUDGET = json.loads((ARENA_DIR / "budget.json").read_text(encoding="utf-8"))
TAU = float(BUDGET["tau_months"])
TREATMENT = base.TREATMENT
PRETREATMENT_COLUMNS = base.PRETREATMENT_COLUMNS
CLINICAL_COLUMNS = base.CLINICAL_COLUMNS
PROPENSITY_CLIP = base.PROPENSITY_CLIP


@dataclass(frozen=True)
class ForestSpec:
    n_estimators: int
    max_depth: int
    min_samples_leaf: int
    min_samples_split: int
    max_features: float | str

    def as_dict(self) -> dict:
        return self.__dict__.copy()


@dataclass(frozen=True)
class CandidateSpec:
    name: str
    parent_run_id: str | None
    max_panel_genes: dict[str, int]
    next_pool_size: dict[str, int]
    benefit_threshold_months: float
    screening: dict[str, ForestSpec]
    final: dict[str, ForestSpec]

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "parent_run_id": self.parent_run_id,
            "max_panel_genes": self.max_panel_genes,
            "next_pool_size": self.next_pool_size,
            "benefit_threshold_months": self.benefit_threshold_months,
            "screening_tlearner": {key: value.as_dict() for key, value in self.screening.items()},
            "tlearner": {key: value.as_dict() for key, value in self.final.items()},
        }


def _validate_forest(raw: dict, label: str, *, screening: bool) -> ForestSpec:
    expected = {"n_estimators", "max_depth", "min_samples_leaf", "min_samples_split", "max_features"}
    if not isinstance(raw, dict) or set(raw) != expected:
        raise ValueError(f"{label} has missing or unrecognized forest keys")
    trees, depth = int(raw["n_estimators"]), int(raw["max_depth"])
    leaf, split = int(raw["min_samples_leaf"]), int(raw["min_samples_split"])
    mtry = raw["max_features"]
    lower = 40 if screening else 100
    if not lower <= trees <= int(BUDGET["max_trees_per_forest"]):
        raise ValueError(f"{label} n_estimators must be in [{lower}, {BUDGET['max_trees_per_forest']}]")
    if not 2 <= depth <= 12 or not 3 <= leaf <= 60 or split < 2 * leaf or split > 150:
        raise ValueError(f"{label} has invalid shallow-forest geometry")
    if isinstance(mtry, str):
        if mtry not in {"sqrt", "log2"}:
            raise ValueError(f"{label} max_features string must be sqrt or log2")
    elif not 0.01 <= float(mtry) <= 1.0:
        raise ValueError(f"{label} max_features must be in [0.01, 1.0]")
    return ForestSpec(trees, depth, leaf, split, mtry if isinstance(mtry, str) else float(mtry))


def validate_candidate(raw: dict) -> CandidateSpec:
    expected = {"name", "parent_run_id", "max_panel_genes", "next_pool_size", "benefit_threshold_months", "screening_tlearner", "tlearner"}
    if not isinstance(raw, dict) or set(raw) != expected:
        raise ValueError("Candidate has missing or unrecognized top-level keys")
    if not isinstance(raw["name"], str) or not raw["name"].strip():
        raise ValueError("Candidate name must be nonempty")
    parent = raw["parent_run_id"]
    if parent is not None and (not isinstance(parent, str) or not parent.startswith("run_")):
        raise ValueError("parent_run_id must be null or a run_* identifier")
    panels = raw["max_panel_genes"]
    pools = raw["next_pool_size"]
    if not isinstance(panels, dict) or set(panels) != {"observation", "act"}:
        raise ValueError("max_panel_genes must contain observation and act")
    if not isinstance(pools, dict) or set(pools) != {"observation", "act"}:
        raise ValueError("next_pool_size must contain observation and act")
    panels = {arm: int(value) for arm, value in panels.items()}
    pools = {arm: int(value) for arm, value in pools.items()}
    if any(value < 0 or value > int(BUDGET["max_panel_genes_per_arm"]) for value in panels.values()) or not any(panels.values()):
        raise ValueError("Arm panel caps must be 0-32 with at least one nonzero")
    if any(value < 1 or value > int(BUDGET["max_pool_genes_per_arm"]) for value in pools.values()):
        raise ValueError("Arm next-pool sizes are outside the arena budget")
    threshold = float(raw["benefit_threshold_months"])
    if not 0 <= threshold <= 3:
        raise ValueError("benefit_threshold_months must be in [0, 3]")
    output: dict[str, dict[str, ForestSpec]] = {}
    for key, is_screen in (("screening_tlearner", True), ("tlearner", False)):
        panel = raw[key]
        if not isinstance(panel, dict) or set(panel) != {"observation", "act"}:
            raise ValueError(f"{key} must contain exactly observation and act")
        output[key] = {arm: _validate_forest(panel[arm], f"{key}.{arm}", screening=is_screen) for arm in panel}
    return CandidateSpec(raw["name"].strip(), parent, panels, pools, threshold, output["screening_tlearner"], output["tlearner"])


def paired_permutation_indices(n: int, *, context: int, gene_index: int, repeat: int) -> np.ndarray:
    seed = int(BUDGET["permutation_seed"]) + context * 10_000_000 + gene_index * 100 + repeat
    return np.random.default_rng(seed).permutation(n)


def component_drop(*, used: bool, baseline: float, permuted: float) -> float:
    """Return an exact zero for a component whose arm forest never used a gene."""
    return float(baseline - permuted) if used else 0.0


def rank_arm_importance(table: pd.DataFrame, arm: str) -> pd.DataFrame:
    """Complete deterministic arm ranking; negative scores remain visible."""
    needed = {"gene", f"{arm}_selection_score", f"{arm}_policy_pfi_mean", "joint_policy_pfi_mean", f"{arm}_cindex_pfi_mean", f"{arm}_positive_fraction"}
    if needed - set(table):
        raise ValueError(f"PFI table is missing {sorted(needed - set(table))}")
    in_pool = f"{arm}_in_pool"
    ranked = table.loc[table[in_pool].astype(bool)].copy() if in_pool in table else table.copy()
    ranked["_positive_score"] = ranked[f"{arm}_selection_score"].clip(lower=0)
    ranked = ranked.sort_values(
        ["_positive_score", f"{arm}_policy_pfi_mean", "joint_policy_pfi_mean", f"{arm}_cindex_pfi_mean", f"{arm}_positive_fraction", "gene"],
        ascending=[False, False, False, False, False, True], kind="stable",
    ).drop(columns="_positive_score").reset_index(drop=True)
    ranked.insert(0, f"{arm}_rank", np.arange(1, len(ranked) + 1))
    return ranked


def select_panels(table: pd.DataFrame, spec: CandidateSpec) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    panels, pools = {}, {}
    for arm in ("observation", "act"):
        ranked = rank_arm_importance(table, arm)
        positive = ranked.loc[ranked[f"{arm}_selection_score"] > 0, "gene"].tolist()
        panels[arm] = positive[: spec.max_panel_genes[arm]]
        pools[arm] = ranked["gene"].tolist()[: spec.next_pool_size[arm]]
    return panels, pools


def _forest(spec: ForestSpec, seed: int, smoke: bool) -> RandomSurvivalForest:
    return RandomSurvivalForest(
        n_estimators=min(int(BUDGET["smoke_trees"]), spec.n_estimators) if smoke else spec.n_estimators,
        max_depth=spec.max_depth, min_samples_leaf=spec.min_samples_leaf,
        min_samples_split=spec.min_samples_split, max_features=spec.max_features,
        bootstrap=True, oob_score=False, n_jobs=1, random_state=int(seed), low_memory=False,
    )


def _fit_arm(fit: pd.DataFrame, genes: list[str], arm: int, spec: ForestSpec, seed: int, smoke: bool):
    transformer = locked_t.FeatureTransformer(genes, "raw", len(genes)).fit(fit)
    mask = fit[TREATMENT].to_numpy(int) == arm
    if int(mask.sum()) < 2 or int(fit.loc[mask, "OS_STATUS"].sum()) < 1:
        raise RuntimeError(f"Insufficient support for arm {arm}")
    model = _forest(spec, seed + arm * 100_003, smoke)
    weights = base._fit_training_weights(fit)
    model.fit(transformer.transform(fit)[mask], base._outcome(fit.loc[mask]), sample_weight=weights[mask])
    if TREATMENT in transformer.names:
        raise RuntimeError("ACT leaked into a feature matrix")
    return model, transformer, {"patients": int(mask.sum()), "events": int(fit.loc[mask, "OS_STATUS"].sum())}


def _arm_predictions(model, transformer, frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    matrix = transformer.transform(frame)
    risk = np.asarray(model.predict(matrix), float)
    rmst = base._integrate_step_functions(model.predict_survival_function(matrix), TAU)
    return risk, rmst


def _used_genes(model, transformer) -> set[str]:
    indices: set[int] = set()
    for estimator in model.estimators_:
        values = np.asarray(estimator.tree_.feature, int)
        indices.update(int(value) for value in values[values >= len(PRETREATMENT_COLUMNS)])
    return {transformer.names[index] for index in indices if index < len(transformer.names)}


def _safe_cindex(frame: pd.DataFrame, risk: np.ndarray) -> float:
    try:
        return float(concordance_index_censored(frame["OS_STATUS"].to_numpy(bool), frame["OS_MONTHS"].to_numpy(float), risk)[0])
    except ValueError:
        return 0.5


def _policy_nuisance(fit: pd.DataFrame, assess: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    _, propensity, _ = base._fit_propensity(fit, assess)
    ipcw = base._ipcw_restricted_time(fit, assess)
    mu0, mu1 = base._clinical_outcome_predictions(fit, assess)
    return base.aipw_arm_scores(assess[TREATMENT].to_numpy(int), propensity, ipcw, mu0, mu1)[:2]


def _importance_table(fit: pd.DataFrame, pools: dict[str, list[str]], spec: CandidateSpec, *, context: int, smoke: bool) -> pd.DataFrame:
    union = sorted(set(pools["observation"]) | set(pools["act"]))
    index_of = {gene: index for index, gene in enumerate(union)}
    repeats = int(BUDGET["smoke_permutation_repeats"] if smoke else BUDGET["permutation_repeats"])
    inner_folds = int(BUDGET["smoke_inner_folds"] if smoke else BUDGET["inner_pfi_folds"])
    metrics = ["observation_policy", "act_policy", "joint_policy", "observation_cindex", "act_cindex"]
    values = {metric: {gene: [] for gene in union} for metric in metrics}
    used_count = {arm: {gene: 0 for gene in union} for arm in ("observation", "act")}
    folds = list(base.fixed_folds(fit, inner_folds, seed=base.FOLD_SEED + context))
    for inner, (train_index, assess_index) in enumerate(folds, 1):
        train = fit.iloc[train_index].reset_index(drop=True)
        assess = fit.iloc[assess_index].reset_index(drop=True)
        models, transforms, baseline_risk, baseline_rmst, used = {}, {}, {}, {}, {}
        for arm_name, arm_value in (("observation", 0), ("act", 1)):
            models[arm_name], transforms[arm_name], _ = _fit_arm(train, pools[arm_name], arm_value, spec.screening[arm_name], context + inner * 1_009, smoke)
            baseline_risk[arm_name], baseline_rmst[arm_name] = _arm_predictions(models[arm_name], transforms[arm_name], assess)
            used[arm_name] = _used_genes(models[arm_name], transforms[arm_name])
            for gene in used[arm_name]:
                used_count[arm_name][gene] += 1
        phi0, phi1 = _policy_nuisance(train, assess)
        baseline_policy = base.recommendations(baseline_rmst["act"] - baseline_rmst["observation"], spec.benefit_threshold_months)
        baseline_alignment = base.policy_summary(phi0, phi1, baseline_policy)["alignment_months"]
        masks = {"observation": assess[TREATMENT].to_numpy(int) == 0, "act": assess[TREATMENT].to_numpy(int) == 1}
        baseline_c = {arm: _safe_cindex(assess.loc[masks[arm]], baseline_risk[arm][masks[arm]]) for arm in masks}
        for gene in union:
            if gene not in used["observation"] and gene not in used["act"]:
                for _ in range(repeats):
                    for metric in metrics:
                        values[metric][gene].append(0.0)
                continue
            for permutation_repeat in range(repeats):
                changed_rmst = {arm: baseline_rmst[arm] for arm in ("observation", "act")}
                changed_risk = {arm: baseline_risk[arm] for arm in ("observation", "act")}
                order = paired_permutation_indices(len(assess), context=context * 100 + inner, gene_index=index_of[gene], repeat=permutation_repeat)
                permuted = assess.copy()
                permuted[gene] = assess[gene].to_numpy()[order]
                for arm in ("observation", "act"):
                    if gene in used[arm]:
                        changed_risk[arm], changed_rmst[arm] = _arm_predictions(models[arm], transforms[arm], permuted)
                obs_policy = base.recommendations(baseline_rmst["act"] - changed_rmst["observation"], spec.benefit_threshold_months)
                act_policy = base.recommendations(changed_rmst["act"] - baseline_rmst["observation"], spec.benefit_threshold_months)
                joint_policy = base.recommendations(changed_rmst["act"] - changed_rmst["observation"], spec.benefit_threshold_months)
                values["observation_policy"][gene].append(component_drop(used=gene in used["observation"], baseline=baseline_alignment, permuted=base.policy_summary(phi0, phi1, obs_policy)["alignment_months"]))
                values["act_policy"][gene].append(component_drop(used=gene in used["act"], baseline=baseline_alignment, permuted=base.policy_summary(phi0, phi1, act_policy)["alignment_months"]))
                values["joint_policy"][gene].append(component_drop(used=gene in used["observation"] or gene in used["act"], baseline=baseline_alignment, permuted=base.policy_summary(phi0, phi1, joint_policy)["alignment_months"]))
                for arm in ("observation", "act"):
                    drop = component_drop(used=gene in used[arm], baseline=baseline_c[arm], permuted=_safe_cindex(assess.loc[masks[arm]], changed_risk[arm][masks[arm]]))
                    values[f"{arm}_cindex"][gene].append(drop)
    rows = []
    for gene in union:
        row: dict[str, Any] = {
            "gene": gene,
            "observation_in_pool": gene in pools["observation"],
            "act_in_pool": gene in pools["act"],
        }
        for metric in metrics:
            array = np.asarray(values[metric][gene], float)
            mean = float(array.mean()) if len(array) else 0.0
            sd = float(array.std(ddof=1)) if len(array) > 1 else 0.0
            prefix = metric.replace("_policy", "_policy_pfi").replace("_cindex", "_cindex_pfi")
            row[f"{prefix}_mean"] = mean
            row[f"{prefix}_sd"] = sd
            row[f"{prefix}_se"] = sd / np.sqrt(max(1, len(array)))
            row[f"{prefix}_positive_fraction"] = float(np.mean(array > 0)) if len(array) else 0.0
        row["observation_positive_fraction"] = row["observation_policy_pfi_positive_fraction"]
        row["act_positive_fraction"] = row["act_policy_pfi_positive_fraction"]
        row["observation_selection_score"] = row["observation_policy_pfi_mean"] - row["observation_policy_pfi_se"]
        row["act_selection_score"] = row["act_policy_pfi_mean"] - row["act_policy_pfi_se"]
        row["observation_split_use_fraction"] = used_count["observation"][gene] / len(folds)
        row["act_split_use_fraction"] = used_count["act"][gene] / len(folds)
        rows.append(row)
    return pd.DataFrame(rows)


def _index_hash(frame: pd.DataFrame) -> str:
    return integrity.canonical_sha256(frame.index.astype(int).tolist())


def _pool_hash(pool: list[str]) -> str:
    return integrity.canonical_sha256(pool)


def load_parent_pools(parent_run_dir: Path, expected_parent: str, genes: list[str]) -> tuple[dict, str]:
    if parent_run_dir.name != expected_parent:
        raise RuntimeError("Parent directory does not match candidate parent_run_id")
    result_path, lineage_path = parent_run_dir / "result.json", parent_run_dir / "pool_lineage.json"
    if not result_path.exists() or not lineage_path.exists():
        raise RuntimeError("Parent is not a completed run with pool lineage")
    result = json.loads(result_path.read_text(encoding="utf-8"))
    lineage_sha = integrity.sha256_file(lineage_path)
    if result.get("artifacts", {}).get("pool_lineage_sha256") != lineage_sha:
        raise RuntimeError("Parent pool-lineage artifact hash mismatch")
    lineage = json.loads(lineage_path.read_text(encoding="utf-8"))
    if lineage.get("run_id") != expected_parent or lineage.get("eligible_gene_universe_sha256") != _pool_hash(genes):
        raise RuntimeError("Parent lineage identity or gene universe is invalid")
    allowed = set(genes)
    for context, payload in lineage.get("contexts", {}).items():
        for arm in ("observation", "act"):
            pool = payload.get("next_pools", {}).get(arm)
            if not isinstance(pool, list) or len(pool) != len(set(pool)) or any(gene not in allowed for gene in pool):
                raise RuntimeError(f"Parent {context}/{arm} pool is invalid")
            if payload.get("next_pool_sha256", {}).get(arm) != _pool_hash(pool):
                raise RuntimeError(f"Parent {context}/{arm} pool hash mismatch")
    return lineage, lineage_sha


def _clinical_params(arm: str) -> ForestSpec:
    raw = BUDGET["clinical_tlearner"][arm]
    return ForestSpec(int(raw["n_estimators"]), int(raw["max_depth"]), int(raw["min_samples_leaf"]), int(raw["min_samples_split"]), float(raw["max_features"]))


def _empty_predictions(n: int) -> dict:
    return {name: {metric: np.zeros(n) for metric in ("risk", "rmst0", "rmst1")} | {"seed_agreement": [], "seed_benefit_correlation": [], "usage": []} for name in ("clinical", "clinical_plus_genomic")}


def _correlation(rows: list[np.ndarray]) -> float:
    output = []
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            output.append(1.0 if np.allclose(rows[i], rows[j]) else (0.0 if np.std(rows[i]) < 1e-12 or np.std(rows[j]) < 1e-12 else float(np.corrcoef(rows[i], rows[j])[0, 1])))
    return float(np.mean(output)) if output else 1.0


def _fit_predict(fit: pd.DataFrame, assess: pd.DataFrame, panels: dict[str, list[str]], spec: CandidateSpec, seeds: tuple[int, ...], smoke: bool) -> dict:
    output = _empty_predictions(len(assess))
    for name, clinical in (("clinical", True), ("clinical_plus_genomic", False)):
        risks, zeros, ones, recs, benefits, usage = [], [], [], [], [], []
        threshold = 0.0 if clinical else spec.benefit_threshold_months
        for seed in seeds:
            genes0, genes1 = ([], []) if clinical else (panels["observation"], panels["act"])
            params0, params1 = (_clinical_params("observation"), _clinical_params("act")) if clinical else (spec.final["observation"], spec.final["act"])
            model0, transform0, support0 = _fit_arm(fit, genes0, 0, params0, seed, smoke)
            model1, transform1, support1 = _fit_arm(fit, genes1, 1, params1, seed, smoke)
            risk0, rmst0 = _arm_predictions(model0, transform0, assess)
            risk1, rmst1 = _arm_predictions(model1, transform1, assess)
            treatment = assess[TREATMENT].to_numpy(int)
            risk = np.where(treatment == 1, risk1, risk0)
            benefit = rmst1 - rmst0
            risks.append(risk); zeros.append(rmst0); ones.append(rmst1); benefits.append(benefit)
            recs.append(base.recommendations(benefit, threshold))
            usage.append({"support": {"observation": support0, "act": support1}, "genomic_split_fraction": {"observation": locked_t._genomic_split_fraction(model0, len(PRETREATMENT_COLUMNS)), "act": locked_t._genomic_split_fraction(model1, len(PRETREATMENT_COLUMNS))}})
        panel = output[name]
        panel["risk"] = np.mean(risks, axis=0); panel["rmst0"] = np.mean(zeros, axis=0); panel["rmst1"] = np.mean(ones, axis=0)
        ensemble = base.recommendations(panel["rmst1"] - panel["rmst0"], threshold)
        panel["seed_agreement"].append(float(np.mean(np.asarray(recs) == ensemble[None, :])))
        panel["seed_benefit_correlation"].append(_correlation(benefits))
        panel["usage"].extend(usage)
    return output


def _mean_nested(rows: list[dict]) -> dict:
    output = {}
    for key in rows[0]:
        vals = [row[key] for row in rows]
        if isinstance(vals[0], dict): output[key] = _mean_nested(vals)
        elif isinstance(vals[0], str): output[key] = vals[0]
        elif isinstance(vals[0], (bool, np.bool_)): output[key] = bool(all(vals))
        else: output[key] = float(np.mean(vals))
    return output


def _summarize(frame: pd.DataFrame, predictions: dict, raw_e: np.ndarray, e: np.ndarray, ipcw: np.ndarray, threshold: float):
    treatment = frame[TREATMENT].to_numpy(int)
    models, phis = {}, {}
    for name in ("clinical", "clinical_plus_genomic"):
        panel = predictions[name]
        phi0, phi1 = base.aipw_arm_scores(treatment, e, ipcw, panel["rmst0"], panel["rmst1"])
        benefit = panel["rmst1"] - panel["rmst0"]
        current_threshold = 0.0 if name == "clinical" else threshold
        rec = base.recommendations(benefit, current_threshold)
        phis[name] = (phi0, phi1, rec)
        use = _mean_nested(panel["usage"])
        models[name] = {**base.policy_summary(phi0, phi1, rec), "harrell_cindex": _safe_cindex(frame, panel["risk"]), "mean_predicted_benefit_months": float(np.mean(benefit)), "predicted_benefit_iqr_months": float(np.ptp(np.percentile(benefit, [25, 75]))), "median_absolute_predicted_benefit_months": float(np.median(np.abs(benefit))), "nontrivial_benefit_fraction": float(np.mean(np.abs(benefit) > max(current_threshold, .10))), "seed_agreement": float(np.mean(panel["seed_agreement"])), "seed_benefit_correlation": float(np.mean(panel["seed_benefit_correlation"])), "act_mechanism": {"role": "arm assignment only; ACT absent from X", "treatment_absent_from_features": True}, "arm_support": use["support"], "genomic_split_fraction": use["genomic_split_fraction"]}
    weights = np.where(treatment == 1, 1 / e, 1 / (1 - e))
    overlap = {"raw_propensity_in_0_05_0_95_fraction": float(np.mean((raw_e >= PROPENSITY_CLIP[0]) & (raw_e <= PROPENSITY_CLIP[1]))), "iptw_effective_sample_size": float(weights.sum() ** 2 / np.square(weights).sum())}
    return models, phis, overlap


def _nonempty_jaccard(rows: list[list[str]]) -> float:
    return base._pairwise_jaccard([row for row in rows if row])


def _plot_pfi(table: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric, title in zip(axes, ["observation_policy_pfi_mean", "act_policy_pfi_mean", "joint_policy_pfi_mean"], ["OBS component", "ACT component", "Joint policy"]):
        shown = table.sort_values([metric, "gene"], ascending=[False, True]).head(20).iloc[::-1]
        ax.barh(shown["gene"], shown[metric]); ax.set_title(title); ax.set_xlabel("held-out policy alignment drop (months)")
    fig.tight_layout(); fig.savefig(path, dpi=160, bbox_inches="tight"); plt.close(fig)


def _plot_single_pfi(table: pd.DataFrame, metric: str, title: str, path: Path) -> None:
    shown = table.sort_values([metric, "gene"], ascending=[False, True]).head(25).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.barh(shown["gene"], shown[metric])
    ax.set_xlabel("held-out policy alignment drop (months)")
    ax.set_title(title)
    fig.tight_layout(); fig.savefig(path, dpi=160, bbox_inches="tight"); plt.close(fig)


def evaluate_candidate(raw_candidate: dict, *, artifact_dir: str | Path, parent_run_dir: str | Path | None = None, smoke: bool = False) -> dict:
    spec = validate_candidate(raw_candidate)
    artifact_dir = Path(artifact_dir); artifact_dir.mkdir(parents=True, exist_ok=True)
    frame, _, all_genes = base.load_development()
    genes = list(all_genes[: int(BUDGET["smoke_genes"])]) if smoke else list(all_genes)
    parent_lineage, parent_sha = (None, None)
    if spec.parent_run_id:
        if smoke:
            raise ValueError("Smoke must exercise a root candidate")
        if parent_run_dir is None:
            raise ValueError("A parent run directory is required")
        parent_lineage, parent_sha = load_parent_pools(Path(parent_run_dir), spec.parent_run_id, genes)
    elif parent_run_dir is not None:
        raise ValueError("Root candidate cannot receive a parent directory")
    n_folds = int(BUDGET["smoke_outer_folds"] if smoke else BUDGET["outer_folds"])
    n_repeats = int(BUDGET["smoke_outer_repeats"] if smoke else BUDGET["outer_repeats"])
    seeds = tuple(BUDGET["forest_seeds"][:1] if smoke else BUDGET["forest_seeds"])
    draws = 100 if smoke else int(BUDGET["patient_bootstraps"])
    repeat_results, fold_panels, lineage_contexts, pfi_frames = [], {"observation": [], "act": []}, {}, []
    for repeat in range(1, n_repeats + 1):
        predictions = _empty_predictions(len(frame)); raw_e = np.zeros(len(frame)); e = np.zeros(len(frame)); ipcw = np.zeros(len(frame))
        fold_seed = base.FOLD_SEED + repeat * 10_000
        for fold, (fit_index, assess_index) in enumerate(base.fixed_folds(frame, n_folds, seed=fold_seed), 1):
            fit, assess = frame.iloc[fit_index].reset_index(drop=True), frame.iloc[assess_index].reset_index(drop=True)
            key = f"repeat_{repeat}_fold_{fold}"
            expected_fit_hash = integrity.canonical_sha256([int(value) for value in fit_index])
            if parent_lineage and parent_lineage["contexts"].get(key, {}).get("fit_index_sha256") != expected_fit_hash:
                raise RuntimeError(f"Parent cross-fold lineage mismatch for {key}")
            pools = {arm: list(parent_lineage["contexts"][key]["next_pools"][arm]) if parent_lineage else list(genes) for arm in ("observation", "act")}
            table = _importance_table(fit, pools, spec, context=repeat * 100 + fold, smoke=smoke)
            panels, next_pools = select_panels(table, spec)
            for arm in panels: fold_panels[arm].append(panels[arm])
            labeled = table.copy(); labeled.insert(0, "context", key); pfi_frames.append(labeled)
            lineage_contexts[key] = {"fit_index_sha256": expected_fit_hash, "input_pool_sha256": {arm: _pool_hash(pools[arm]) for arm in pools}, "selected_panels": panels, "selected_panel_sha256": {arm: _pool_hash(panels[arm]) for arm in panels}, "next_pools": next_pools, "next_pool_sha256": {arm: _pool_hash(next_pools[arm]) for arm in next_pools}}
            fold_predictions = _fit_predict(fit, assess, panels, spec, seeds, smoke)
            for name in predictions:
                for metric in ("risk", "rmst0", "rmst1"): predictions[name][metric][assess_index] = fold_predictions[name][metric]
                for metric in ("seed_agreement", "seed_benefit_correlation", "usage"): predictions[name][metric].extend(fold_predictions[name][metric])
            fold_raw, fold_e, _ = base._fit_propensity(fit, assess)
            raw_e[assess_index], e[assess_index], ipcw[assess_index] = fold_raw, fold_e, base._ipcw_restricted_time(fit, assess)
            print(f"[repeat {repeat} fold {fold}/{n_folds}] OBS={len(panels['observation'])} ACT={len(panels['act'])}", flush=True)
        models, phis, overlap = _summarize(frame, predictions, raw_e, e, ipcw, spec.benefit_threshold_months)
        c0, c1, cr = phis["clinical"]; g0, g1, gr = phis["clinical_plus_genomic"]
        increment = models["clinical_plus_genomic"]["alignment_months"] - models["clinical"]["alignment_months"]
        repeat_results.append({"repeat": repeat, "fold_seed": fold_seed, "incremental_alignment_months": float(increment), "paired_bootstrap": base._bootstrap_increment(c0, c1, cr, g0, g1, gr, draws=draws), "models": models, "constant_policy_values": {"all_observation_months": float(np.mean(g0)), "all_act_months": float(np.mean(g1))}, "overlap": overlap})
    full_key = "full_development"
    full_fit_hash = integrity.canonical_sha256(list(range(len(frame))))
    if parent_lineage and parent_lineage["contexts"].get(full_key, {}).get("fit_index_sha256") != full_fit_hash:
        raise RuntimeError("Parent full-development lineage mismatch")
    full_pools = {arm: list(parent_lineage["contexts"][full_key]["next_pools"][arm]) if parent_lineage else list(genes) for arm in ("observation", "act")}
    full_table = _importance_table(frame.reset_index(drop=True), full_pools, spec, context=999, smoke=smoke)
    full_panels, full_next = select_panels(full_table, spec)
    labeled = full_table.copy(); labeled.insert(0, "context", full_key); pfi_frames.append(labeled)
    lineage_contexts[full_key] = {"fit_index_sha256": full_fit_hash, "input_pool_sha256": {arm: _pool_hash(full_pools[arm]) for arm in full_pools}, "selected_panels": full_panels, "selected_panel_sha256": {arm: _pool_hash(full_panels[arm]) for arm in full_panels}, "next_pools": full_next, "next_pool_sha256": {arm: _pool_hash(full_next[arm]) for arm in full_next}}
    lineage = {"schema_version": 1, "run_id": None, "parent_run_id": spec.parent_run_id, "parent_pool_lineage_sha256": parent_sha, "eligible_gene_universe_sha256": _pool_hash(genes), "contexts": lineage_contexts}
    lineage_path = artifact_dir / "pool_lineage.json"; lineage_path.write_text(json.dumps(lineage, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    combined = pd.concat(pfi_frames, ignore_index=True)
    for arm, columns in (("observation", ["context", "gene", "observation_policy_pfi_mean", "observation_policy_pfi_sd", "observation_policy_pfi_se", "observation_selection_score", "observation_positive_fraction", "observation_cindex_pfi_mean", "observation_split_use_fraction"]), ("act", ["context", "gene", "act_policy_pfi_mean", "act_policy_pfi_sd", "act_policy_pfi_se", "act_selection_score", "act_positive_fraction", "act_cindex_pfi_mean", "act_split_use_fraction"]), ("joint", ["context", "gene", "joint_policy_pfi_mean", "joint_policy_pfi_sd", "joint_policy_pfi_se", "joint_policy_pfi_positive_fraction"])):
        combined[columns].to_csv(artifact_dir / f"{arm}_pfi.csv", index=False)
    _plot_pfi(full_table, artifact_dir / "permutation_importance.png")
    _plot_single_pfi(full_table, "observation_policy_pfi_mean", "OBS-component policy PFI", artifact_dir / "observation_pfi.png")
    _plot_single_pfi(full_table, "act_policy_pfi_mean", "ACT-component policy PFI", artifact_dir / "act_pfi.png")
    _plot_single_pfi(full_table, "joint_policy_pfi_mean", "Joint policy PFI", artifact_dir / "joint_pfi.png")
    (artifact_dir / "selected_panels.json").write_text(json.dumps({"per_fold": fold_panels, "full_development": full_panels}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    increments = np.asarray([row["incremental_alignment_months"] for row in repeat_results])
    development_models = _mean_nested([row["models"] for row in repeat_results]); development_overlap = _mean_nested([row["overlap"] for row in repeat_results]); constants = _mean_nested([row["constant_policy_values"] for row in repeat_results])
    obs_j, act_j = _nonempty_jaccard(fold_panels["observation"]), _nonempty_jaccard(fold_panels["act"])
    selected_scores = []
    for arm in ("observation", "act"):
        selected_scores.extend(full_table.set_index("gene").loc[full_panels[arm], f"{arm}_selection_score"].tolist() if full_panels[arm] else [])
    gates = {
        "all_repeat_genomic_increment_positive": bool(np.all(increments > 0)),
        "all_repeat_genomic_alignment_positive": all(row["models"]["clinical_plus_genomic"]["alignment_months"] > 0 for row in repeat_results),
        "all_repeat_genomic_value_at_least_clinical": all(row["models"]["clinical_plus_genomic"]["value_months"] >= row["models"]["clinical"]["value_months"] for row in repeat_results),
        "all_repeat_genomic_value_at_least_best_constant": all(row["models"]["clinical_plus_genomic"]["value_months"] >= max(row["constant_policy_values"].values()) for row in repeat_results),
        "all_repeat_cindex_drop_no_more_than_0_03": all(row["models"]["clinical_plus_genomic"]["harrell_cindex"] >= row["models"]["clinical"]["harrell_cindex"] - .03 for row in repeat_results),
        "all_repeat_genomic_seed_agreement_at_least_0_85": all(row["models"]["clinical_plus_genomic"]["seed_agreement"] >= .85 for row in repeat_results),
        "all_repeat_genomic_seed_benefit_correlation_at_least_0_50": all(row["models"]["clinical_plus_genomic"]["seed_benefit_correlation"] >= .50 for row in repeat_results),
        "all_repeat_treatment_absent_from_features": True,
        "all_repeat_mean_genomic_split_fraction_at_least_0_01": all(np.mean(list(row["models"]["clinical_plus_genomic"]["genomic_split_fraction"].values())) >= .01 for row in repeat_results),
        "observation_panel_jaccard_at_least_0_10_when_nonempty": obs_j >= .10,
        "act_panel_jaccard_at_least_0_10_when_nonempty": act_j >= .10,
        "all_repeat_nontrivial_benefit_fraction_at_least_0_10": all(row["models"]["clinical_plus_genomic"]["nontrivial_benefit_fraction"] >= .10 for row in repeat_results),
        "all_repeat_raw_propensity_overlap_at_least_0_80": all(row["overlap"]["raw_propensity_in_0_05_0_95_fraction"] >= .80 for row in repeat_results),
        "all_repeat_iptw_effective_sample_size_at_least_0_30n": all(row["overlap"]["iptw_effective_sample_size"] >= .30 * len(frame) for row in repeat_results),
        "fit_only_parent_pool_lineage_valid": True,
        "panel_caps_respected": all(len(panel) <= 32 for arm in fold_panels.values() for panel in arm) and all(len(panel) <= 32 for panel in full_panels.values()),
        "at_least_one_selected_genomic_feature": bool(selected_scores),
        "all_selected_genes_have_positive_pfi_selection_score": bool(selected_scores) and all(value > 0 for value in selected_scores),
    }
    eligible = all(gates.values())
    robust_lcb = min(row["paired_bootstrap"]["selection_lcb"] for row in repeat_results)
    increment_range = float(np.ptp(increments)) if len(increments) > 1 else 0.0
    score = float(robust_lcb - increment_range)
    result = {
        "schema_version": 1, "phase": "smoke" if smoke else "pooled_development_repeated_crossfit", "candidate": spec.as_dict(),
        "objective": "worst repeated-development OOF selection LCB minus repeat range for A60(C+G)-A60(C)",
        "estimator": "arm-specific RSF T-learners with held-out policy permutation importance and IPCW-AIPW value",
        "reward": score if eligible else -1_000_000.0, "eligible": bool(eligible), "gates": gates,
        "beats_run20_reward": bool(eligible and score > float(BUDGET["run20_reward_to_beat"])),
        "generalization": {"robust_selection_lcb": float(robust_lcb), "repeat_increment_range_months": increment_range, "score_before_eligibility": score},
        "development_oof": {"incremental_alignment_months": float(np.mean(increments)), "models": development_models, "constant_policy_values": constants, "overlap": development_overlap, "repeat_increment_min_months": float(np.min(increments)), "repeat_increment_max_months": float(np.max(increments))},
        "repeat_results": repeat_results,
        "gene_selection": {"observation_mean_pairwise_jaccard": obs_j, "act_mean_pairwise_jaccard": act_j, "per_repeat_outer_fold": fold_panels, "full_development_selection": full_panels, "full_development_next_pools": full_next},
        "permutation_importance": {"selection_score": "mean arm policy PFI minus one standard error", "unused_arm_importance": "exactly zero", "negative_importance_selection_contribution": 0, "shared_permutation_indices": True, "rank_rule": "positive selection score, arm PFI, joint PFI, arm C-index PFI, positive fraction, gene symbol"},
        "diagnostics": {"panel_stability": {"observation_jaccard": obs_j, "act_jaccard": act_j}, "arm_support": development_models["clinical_plus_genomic"]["arm_support"], "split_use": development_models["clinical_plus_genomic"]["genomic_split_fraction"], "seed_stability": {key: development_models["clinical_plus_genomic"][key] for key in ("seed_agreement", "seed_benefit_correlation")}, "policy": {"threshold_months": spec.benefit_threshold_months}},
        "data": {"development_n": int(len(frame)), "development_events": int(frame["OS_STATUS"].sum()), "development_act": int(frame[TREATMENT].sum()), "eligible_reactome_genes": len(genes), "outer_folds": n_folds, "outer_repeats": n_repeats, "forest_seeds": list(seeds), "former_validation_used_as_development": True, "confirmatory_validation_used": False, "test_used": False, "train_sha256": sha256_file(base.TRAIN_CSV), "validation_sha256": sha256_file(base.VALID_CSV), "reactome_collection_sha256": sha256_file(base.REACTOME_DATA_DIR / base.MSIGDB_FILENAME)},
        "environment": integrity.current_environment(), "budget": {**BUDGET, "smoke_reductions_applied": bool(smoke)},
        "artifacts": {"observation_pfi": "observation_pfi.csv", "act_pfi": "act_pfi.csv", "joint_pfi": "joint_pfi.csv", "observation_pfi_plot": "observation_pfi.png", "act_pfi_plot": "act_pfi.png", "joint_pfi_plot": "joint_pfi.png", "pfi_plot": "permutation_importance.png", "selected_panels": "selected_panels.json", "diagnostics": "diagnostics.json", "pool_lineage": "pool_lineage.json", "pool_lineage_sha256": integrity.sha256_file(lineage_path)},
        "limitations": ["PFI is conditional on fitted screening forests and may divide credit among correlated genes.", "Repeated cross-fits are stability views, not independent validation cohorts.", "The sealed test set is inaccessible to this arena."],
    }
    (artifact_dir / "candidate.json").write_text(json.dumps(spec.as_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (artifact_dir / "diagnostics.json").write_text(json.dumps(result["diagnostics"], indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result
