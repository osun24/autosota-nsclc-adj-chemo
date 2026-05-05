"""Agent-editable Random Survival Forest arena entry point."""

from __future__ import annotations

from pathlib import Path
import inspect
import json
import os
import pickle
import random
import time
import warnings

import numpy as np
import optuna
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.utils import restricted_mean_survival_time
from optuna.samplers import NSGAIISampler
from sksurv.ensemble import RandomSurvivalForest
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv

try:
    from . import prepare
except ImportError:
    import prepare

warnings.filterwarnings("ignore", message="Ties in event time detected; using efron's method to handle ties.")
optuna.logging.set_verbosity(optuna.logging.WARNING)

np.random.seed(42)
random.seed(42)

ARENA_DIR = Path(__file__).resolve().parent
RUNS_DIR = ARENA_DIR / "runs"
DEFAULT_N_TRIALS = int(os.environ.get("RSF_ARENA_N_TRIALS", "30"))
DEFAULT_BOOTSTRAPS = int(os.environ.get("RSF_ARENA_BOOTSTRAPS", "3"))
DEFAULT_SEEDS = int(os.environ.get("RSF_ARENA_SEEDS", "3"))
BOOTSTRAP_BASE_SEED = 31415
FEAT_EVENT_FRACTION = 0.50


def rank_genes_univariate(train_df: pd.DataFrame, gene_cols: list[str]) -> list[str]:
    y = Surv.from_arrays(
        event=train_df["OS_STATUS"].astype(bool).values,
        time=train_df["OS_MONTHS"].values.astype(float),
    )
    ranks = []
    for g in gene_cols:
        Xg = train_df[[g]].to_numpy(dtype=np.float64)
        try:
            model = CoxPHSurvivalAnalysis(alpha=1e-12)
            model.fit(Xg, y)
            pred = model.predict(Xg)
            ci = concordance_index_censored(y["event"], y["time"], pred)[0]
            ranks.append((g, float(ci)))
        except Exception:
            ranks.append((g, 0.5))
    ranks.sort(key=lambda z: z[1], reverse=True)
    return [g for g, _ in ranks]


def rank_genes_stability(train_df: pd.DataFrame, gene_cols: list[str], n_bootstraps: int = 25, seed: int = 42) -> list[str]:
    """Bootstrap stability selection: rank genes by how consistently they rank in the top half across half-samples."""
    rng = np.random.default_rng(seed)
    n = len(train_df)
    sub_n = max(10, n // 2)
    top_k = max(1, len(gene_cols) // 2)
    scores: dict[str, float] = {g: 0.0 for g in gene_cols}
    for _ in range(n_bootstraps):
        idx = rng.choice(n, size=sub_n, replace=False)
        sub_df = train_df.iloc[idx]
        sub_ranked = rank_genes_univariate(sub_df, gene_cols)
        for g in sub_ranked[:top_k]:
            scores[g] += 1.0 / n_bootstraps
    return sorted(gene_cols, key=lambda g: scores[g], reverse=True)


def build_features_with_interactions(
    df: pd.DataFrame,
    main_genes: list[str],
    inter_genes: list[str],
    act_col: str = "Adjuvant Chemo",
    dup_inter: int = 1,
    clin_cols: list[str] | None = None,
) -> tuple[np.ndarray, list[str]]:
    if clin_cols is None:
        clin_cols = prepare.CLINICAL_VARS
    base_cols = list(clin_cols) + list(main_genes)
    X_base = df[base_cols].to_numpy(dtype=np.float64)
    A = df[act_col].to_numpy(dtype=np.float64).reshape(-1, 1)
    names = list(base_cols)
    blocks = [X_base]
    if inter_genes:
        X_int = df[list(inter_genes)].to_numpy(dtype=np.float64) * A
        blocks.append(X_int)
        names += [f"{g}*ACT" for g in inter_genes]
        if int(dup_inter) > 1:
            for d in range(1, int(dup_inter)):
                blocks.append(X_int.copy())
                names += [f"{g}*ACT#dup{d}" for g in inter_genes]
    return (np.concatenate(blocks, axis=1) if len(blocks) > 1 else X_base).astype(np.float64), names


def make_rsf(**params) -> RandomSurvivalForest:
    sig = inspect.signature(RandomSurvivalForest.__init__)
    allowed = set(sig.parameters.keys()) - {"self"}
    return RandomSurvivalForest(**{k: v for k, v in params.items() if k in allowed})


def build_trial_mats_for_splits(
    train_fit_df: pd.DataFrame,
    valid_eval_df: pd.DataFrame,
    w_fit: np.ndarray | None,
    k_main: int,
    k_int: int,
    dup_inter: int,
    gene_rank: list[str],
    clin_cols: list[str],
):
    genes_main = gene_rank[:k_main]
    genes_inter = genes_main[:k_int]
    Xtr, feat_names = build_features_with_interactions(
        train_fit_df, genes_main, genes_inter, dup_inter=dup_inter, clin_cols=clin_cols
    )
    Xva, _ = build_features_with_interactions(valid_eval_df, genes_main, genes_inter, dup_inter=dup_inter, clin_cols=clin_cols)
    ytr = Surv.from_arrays(
        event=train_fit_df["OS_STATUS"].astype(bool).values,
        time=train_fit_df["OS_MONTHS"].astype(float).values,
    )
    return Xtr, ytr, w_fit, Xva, feat_names, genes_main, genes_inter


def suggest_hparams(trial: optuna.Trial, feat_budget: int, clin_cols: list[str], max_genes: int):
    max_nonclin = max(8, feat_budget - len(clin_cols))
    # Reserve budget headroom for interaction terms: main + inter <= max_nonclin
    base_main = [16, 32, 64, 96, 128, 192, max_genes]
    topk_main_choices = tuple(sorted({k for k in base_main if 1 <= k <= min(max_genes, max_nonclin)}))
    if not topk_main_choices:
        topk_main_choices = (min(max_genes, max_nonclin),)
    k_main = int(trial.suggest_categorical("top_k_genes", topk_main_choices))
    # Allow optimizer to choose interaction count within remaining budget
    max_k_int = min(k_main, max(0, max_nonclin - k_main), 32)
    k_int = int(trial.suggest_int("k_int", 0, max_k_int)) if max_k_int > 0 else 0
    dup_inter = 1
    mf_mode = trial.suggest_categorical("max_features_mode", ["sqrt", "frac"])
    max_features = trial.suggest_float("max_features_frac", 0.05, 0.40) if mf_mode == "frac" else mf_mode
    params = {
        "n_estimators": trial.suggest_int("n_estimators", 200, 1000, step=100),
        "max_depth": trial.suggest_int("max_depth", 2, 8),
        "min_samples_split": trial.suggest_int("min_samples_split", 2, 40),
        "min_samples_leaf": trial.suggest_int("min_samples_leaf", 10, 120),
        "max_features": max_features,
        "bootstrap": True,
        "oob_score": False,
        "n_jobs": -1,
        "random_state": 42,
        "verbose": 0,
    }
    return k_main, k_int, dup_inter, params


def _select_pareto_compromise(study: optuna.Study, w_ci: float = 0.40, w_rmst: float = 0.60) -> optuna.trial.FrozenTrial:
    """Select from Pareto front using normalized weighted score over *all* completed trials.

    Normalization pools all trials (not just Pareto) so that two-member fronts don't
    degenerate to a coin-flip tie. RMST is weighted higher than CI since it is the
    harder, clinically-meaningful HTE objective.
    """
    all_complete = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE and t.values is not None]
    if not all_complete:
        raise RuntimeError("No completed multi-objective trials found.")
    pareto = [t for t in study.best_trials if t.state == optuna.trial.TrialState.COMPLETE and t.values is not None]
    if not pareto:
        pareto = all_complete

    # Normalise using global (all-trial) range so Pareto ties are broken consistently.
    ci_all = np.array([float(t.values[0]) for t in all_complete])
    rmst_all = np.array([float(t.values[1]) for t in all_complete])
    ci_lo, ci_hi = float(ci_all.min()), float(ci_all.max())
    rmst_lo, rmst_hi = float(rmst_all.min()), float(rmst_all.max())

    def norm(x, lo, hi):
        return 0.0 if hi <= lo else (x - lo) / (hi - lo)

    scored = []
    for t in pareto:
        score = w_ci * norm(float(t.values[0]), ci_lo, ci_hi) + w_rmst * norm(float(t.values[1]), rmst_lo, rmst_hi)
        scored.append((score, t))
    return sorted(scored, key=lambda z: z[0], reverse=True)[0][1]


def _rsf_params_from_chosen(chosen: optuna.trial.FrozenTrial, random_state: int = 7) -> dict:
    mf_mode = chosen.params.get("max_features_mode", "sqrt")
    max_features = chosen.params.get("max_features_frac", 0.5) if mf_mode == "frac" else mf_mode
    return {
        "n_estimators": int(chosen.params.get("n_estimators", 300)),
        "max_depth": int(chosen.params.get("max_depth", 5)),
        "min_samples_split": int(chosen.params.get("min_samples_split", 10)),
        "min_samples_leaf": int(chosen.params.get("min_samples_leaf", 30)),
        "max_features": max_features,
        "bootstrap": True,
        "oob_score": False,
        "n_jobs": -1,
        "random_state": int(random_state),
        "verbose": 0,
    }


def _save_run_artifacts(result: dict, model, metadata: dict, feature_names: list[str], genes_main: list[str], genes_inter: list[str]) -> Path:
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    run_dir = RUNS_DIR / time.strftime("smoke_%Y%m%d_%H%M%S")
    run_dir.mkdir(parents=True, exist_ok=False)
    with open(run_dir / "rsf_model.pkl", "wb") as f:
        pickle.dump(model, f)
    (run_dir / "feature_names.txt").write_text("\n".join(feature_names) + "\n")
    (run_dir / "genes_main.txt").write_text("\n".join(genes_main) + "\n")
    (run_dir / "genes_inter.txt").write_text("\n".join(genes_inter) + "\n")
    (run_dir / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    (run_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return run_dir


def fit_model_from_metadata(metadata: dict, train_df: pd.DataFrame, valid_df: pd.DataFrame, refit_train_valid: bool = False):
    clin_cols = list(metadata["clin_cols"])
    clin_pretx = [c for c in clin_cols if c != "Adjuvant Chemo"]
    genes_main = list(metadata["genes_main"])
    genes_inter = list(metadata.get("genes_inter", []))
    dup_inter = int(metadata.get("dup_inter", 1))
    params = dict(metadata["rsf_params"])
    params["n_jobs"] = 1  # deterministic final fit: random_state fully controls tree ordering
    fit_df = pd.concat([train_df, valid_df], axis=0, ignore_index=True) if refit_train_valid else train_df
    eval_df = valid_df
    fit_df = fit_df.sort_values(["OS_MONTHS", "OS_STATUS"], ascending=[False, False]).reset_index(drop=True)
    w_fit, _, _ = prepare.compute_iptw(fit_df, covariate_cols=clin_pretx, act_col="Adjuvant Chemo")
    X_fit, feat_names = build_features_with_interactions(fit_df, genes_main, genes_inter, dup_inter=dup_inter, clin_cols=clin_cols)
    y_fit = Surv.from_arrays(event=fit_df["OS_STATUS"].astype(bool).values, time=fit_df["OS_MONTHS"].astype(float).values)
    model = make_rsf(**params)
    model.fit(X_fit, y_fit, sample_weight=w_fit)
    return model, feat_names, genes_main, genes_inter, clin_cols


def run(n_trials: int = DEFAULT_N_TRIALS, bootstrap_n: int = DEFAULT_BOOTSTRAPS, seed_eval_n: int = DEFAULT_SEEDS, save_artifacts: bool = True) -> dict:
    start = time.time()
    train_df, valid_df = prepare.load_train_valid()
    clin_cols, clin_pretx, gene_feats = prepare.clinical_and_gene_columns(train_df, valid_df)

    gene_rank = rank_genes_stability(train_df, gene_feats, n_bootstraps=25, seed=42)
    max_genes = len(gene_rank)
    feat_budget = max(24, int(FEAT_EVENT_FRACTION * int(train_df["OS_STATUS"].sum())))
    print(f"[Gene Ranking] Stability-ranked {max_genes} genes on full TRAIN (25 half-sample bootstraps)")
    print(f"[Budgets] feature budget <= {feat_budget}")
    print(f"Starting bootstrap optimization: {n_trials} trials x {bootstrap_n} bootstraps x {seed_eval_n} seeds")

    def objective(trial: optuna.Trial):
        k_main, k_int, dup_inter, rsf_params = suggest_hparams(trial, feat_budget, clin_cols, max_genes)
        boot_cis, boot_rmsts = [], []
        n_features = None
        for b in range(int(bootstrap_n)):
            boot_seed = BOOTSTRAP_BASE_SEED + trial.number * 1000 + b
            boot_train_df = prepare.bootstrap_resample_df(train_df, boot_seed, require_two_arms=True, require_event=True)
            w_boot, _, _ = prepare.compute_iptw(boot_train_df, covariate_cols=clin_pretx, act_col="Adjuvant Chemo")
            Xtr, ytr, wtr, Xva, feat_names, genes_main, genes_inter = build_trial_mats_for_splits(
                boot_train_df, valid_df, w_boot, k_main, k_int, dup_inter, gene_rank, clin_cols
            )
            seed_cis, seed_rmsts = [], []
            for s in range(int(seed_eval_n)):
                params = dict(rsf_params)
                params["random_state"] = int(boot_seed * 100 + s)
                model = make_rsf(**params)
                model.fit(Xtr, ytr, sample_weight=wtr)
                scored = prepare.evaluate_on_valid(model, valid_df, genes_main, genes_inter, dup_inter, feat_names, clin_cols)
                seed_cis.append(float(scored["val_ci"]))
                seed_rmsts.append(float(scored["val_rmst_diff"]))
                n_features = int(scored["n_features"])
            boot_cis.append(float(np.median(seed_cis)))
            boot_rmsts.append(float(np.median(seed_rmsts)))
        val_ci = float(np.median(boot_cis))
        rmst = float(np.median(boot_rmsts))
        val_ci_se = float(np.std(boot_cis, ddof=1)) if len(boot_cis) > 1 else 0.0
        rmst_iqr = float(np.percentile(boot_rmsts, 75) - np.percentile(boot_rmsts, 25)) if len(boot_rmsts) > 1 else 0.0
        trial.set_user_attr("n_features", int(n_features))
        trial.set_user_attr("k_main", int(k_main))
        trial.set_user_attr("k_int", int(k_int))
        trial.set_user_attr("dup_inter", int(dup_inter))
        trial.set_user_attr("val_ci_boot_se", val_ci_se)
        trial.set_user_attr("rmst_diff_boot_iqr", rmst_iqr)
        print(f"[Trial {trial.number:03d}] Boot Val CI={val_ci:.4f} (SE={val_ci_se:.4f}), RMST={rmst:.4f} (IQR={rmst_iqr:.4f}), N_feats={n_features}")
        return val_ci, rmst

    study = optuna.create_study(directions=["maximize", "maximize"], sampler=NSGAIISampler(seed=42))
    study.optimize(objective, n_trials=int(n_trials), gc_after_trial=True)
    chosen = _select_pareto_compromise(study)
    k_main = int(chosen.user_attrs["k_main"])
    k_int = int(chosen.user_attrs["k_int"])
    dup_inter = int(chosen.user_attrs.get("dup_inter", 1))
    genes_main = gene_rank[:k_main]
    genes_inter = genes_main[:k_int]
    rsf_params = _rsf_params_from_chosen(chosen, random_state=7)
    metadata = {
        "arena": "rsf",
        "clin_cols": clin_cols,
        "clin_pretx": clin_pretx,
        "genes_main": genes_main,
        "genes_inter": genes_inter,
        "dup_inter": dup_inter,
        "rsf_params": rsf_params,
        "chosen_trial_params": dict(chosen.params),
    }
    # Seed ensemble: average risk predictions across 10 seeds for stable CI and RMST.
    FINAL_SEEDS = [7, 13, 21, 37, 53, 71, 89, 97, 113, 127]
    ens_risks, ens_risks_treated, ens_risks_untreated = [], [], []
    seed_cis, seed_rmsts = [], []
    last_model, feat_names = None, None
    valid_tr = valid_df.copy(); valid_tr["Adjuvant Chemo"] = 1
    valid_co = valid_df.copy(); valid_co["Adjuvant Chemo"] = 0
    for fs in FINAL_SEEDS:
        meta_seed = dict(metadata)
        meta_seed["rsf_params"] = dict(rsf_params)
        meta_seed["rsf_params"]["random_state"] = int(fs)
        m, fn, _, _, _ = fit_model_from_metadata(meta_seed, train_df, valid_df, refit_train_valid=False)
        X_v = prepare.build_matrix_from_feature_names(valid_df, fn)
        X_tr = prepare.build_matrix_from_feature_names(valid_tr, fn)
        X_co = prepare.build_matrix_from_feature_names(valid_co, fn)
        r = prepare.predict_rsf_risk(m, X_v)
        r_tr = prepare.predict_rsf_risk(m, X_tr)
        r_co = prepare.predict_rsf_risk(m, X_co)
        ens_risks.append(r); ens_risks_treated.append(r_tr); ens_risks_untreated.append(r_co)
        seed_cis.append(prepare.cindex(r, valid_df["OS_MONTHS"].to_numpy(float), valid_df["OS_STATUS"].to_numpy(int)))
        vr = prepare.evaluate_on_valid(m, valid_df, genes_main, genes_inter, dup_inter, fn, clin_cols)
        seed_rmsts.append(float(vr["val_rmst_diff"]))
        last_model, feat_names = m, fn
    # Ensemble predictions
    ens_r = np.mean(ens_risks, axis=0)
    ens_r_tr = np.mean(ens_risks_treated, axis=0)
    ens_r_co = np.mean(ens_risks_untreated, axis=0)
    ens_ci = prepare.cindex(ens_r, valid_df["OS_MONTHS"].to_numpy(float), valid_df["OS_STATUS"].to_numpy(int))
    model_rec = np.where(ens_r_tr < ens_r_co, 1, 0)
    alignment = valid_df["Adjuvant Chemo"].to_numpy(int) == model_rec
    if int(alignment.sum()) > 0 and int((~alignment).sum()) > 0:
        km_a = KaplanMeierFitter().fit(valid_df.loc[alignment, "OS_MONTHS"], event_observed=valid_df.loc[alignment, "OS_STATUS"])
        km_n = KaplanMeierFitter().fit(valid_df.loc[~alignment, "OS_MONTHS"], event_observed=valid_df.loc[~alignment, "OS_STATUS"])
        ens_rmst = float(restricted_mean_survival_time(km_a, t=60) - restricted_mean_survival_time(km_n, t=60))
    else:
        ens_rmst = 0.0
    print(f"[Seed Panel] val_ci per seed:    {[round(v,4) for v in seed_cis]}")
    print(f"[Seed Panel] val_rmst per seed:  {[round(v,3) for v in seed_rmsts]}")
    print(f"[Ensemble ] val_ci={ens_ci:.4f}, val_rmst_diff={ens_rmst:.3f}")
    model = last_model
    result = {
        "val_ci": float(ens_ci),
        "val_rmst_diff": float(ens_rmst),
        "val_ci_se": float(np.std(seed_cis, ddof=1)),
        "rmst_iqr": float(np.percentile(seed_rmsts, 75) - np.percentile(seed_rmsts, 25)),
        "seed_panel_cis": seed_cis,
        "seed_panel_rmsts": seed_rmsts,
        "ensemble_ci": float(ens_ci),
        "ensemble_rmst_diff": float(ens_rmst),
        "n_features": len(feat_names),
        "k_main": k_main,
        "k_int": k_int,
        "dup_inter": dup_inter,
        "bootstrap_n": int(bootstrap_n),
        "seed_eval_n": int(seed_eval_n),
        "n_trials": int(n_trials),
        "chosen_trial": int(chosen.number),
        "notes": "Phase 1 RSF scaffold: existing RSF logic with reduced Optuna/bootstrap/seed budget, no test access.",
        "wall_clock_sec": round(time.time() - start, 2),
    }
    metadata["feature_names"] = feat_names
    metadata["result"] = result
    if save_artifacts:
        run_dir = _save_run_artifacts(result, model, metadata, feat_names, genes_main, genes_inter)
        result["run_dir"] = str(run_dir.relative_to(ARENA_DIR.parent))
        print(f"[Artifacts] Saved RSF artifacts to {run_dir}")
    print("\n[Result Row]")
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


if __name__ == "__main__":
    run()
