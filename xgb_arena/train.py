"""Agent-editable XGBoost-Cox arena entry point.

The autonomous loop edits this file. It may import frozen utilities from
prepare.py, but it must not load or reference the test set.
"""

from __future__ import annotations

from pathlib import Path
import json
import os
import random
import time
import warnings

import numpy as np
import optuna
import pandas as pd
import xgboost as xgb
from optuna.samplers import NSGAIISampler
from sklearn.model_selection import train_test_split
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv

try:
    from . import prepare
except ImportError:
    import prepare

warnings.filterwarnings(
    "ignore",
    message="Ties in event time detected; using efron's method to handle ties.",
)
optuna.logging.set_verbosity(optuna.logging.WARNING)

np.random.seed(42)
random.seed(42)

ARENA_DIR = Path(__file__).resolve().parent
RUNS_DIR = ARENA_DIR / "runs"

DEFAULT_N_TRIALS = int(os.environ.get("XGB_ARENA_N_TRIALS", "10"))
DEFAULT_BOOTSTRAPS = int(os.environ.get("XGB_ARENA_BOOTSTRAPS", "2"))
BOOTSTRAP_BASE_SEED = 42
FEAT_EVENT_FRACTION = 0.50


def rank_genes_univariate(train_df: pd.DataFrame, gene_cols: list[str]) -> list[str]:
    y = Surv.from_arrays(
        event=train_df["OS_STATUS"].astype(bool).values,
        time=train_df["OS_MONTHS"].values.astype(float),
    )
    ranks = []
    for g in gene_cols:
        Xg = train_df[[g]].to_numpy(dtype=np.float32)
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


def stability_selection_genes(
    train_df: pd.DataFrame,
    gene_cols: list[str],
    n_subsamples: int = 100,
    subsample_frac: float = 0.5,
    stability_threshold: float = 0.6,
    target_selections: int = 50,
    seed: int = 42,
) -> list[str]:
    """Rank genes by LASSO-Cox stability selection frequency over half-subsamples."""
    from sksurv.linear_model import CoxnetSurvivalAnalysis

    y = Surv.from_arrays(
        event=train_df["OS_STATUS"].astype(bool).values,
        time=train_df["OS_MONTHS"].values.astype(float),
    )
    X = train_df[gene_cols].to_numpy(dtype=np.float64)
    X_mean = X.mean(axis=0)
    X_std = X.std(axis=0)
    X_std = np.where(X_std > 0, X_std, 1.0)
    X_scaled = ((X - X_mean) / X_std).astype(np.float64)

    n = len(train_df)
    n_sub = int(n * subsample_frac)
    rng = np.random.default_rng(seed)

    try:
        cox_full = CoxnetSurvivalAnalysis(
            l1_ratio=1.0, fit_baseline_model=False,
            max_iter=500, n_alphas=30,
        )
        cox_full.fit(X_scaled, y)
        alphas = cox_full.alphas_
        coef_path = cox_full.coef_  # shape: (n_features, n_alphas)
        n_selected = (np.abs(coef_path) > 1e-8).sum(axis=0)
        candidates = np.where(n_selected >= target_selections)[0]
        alpha_idx = int(candidates[0]) if len(candidates) > 0 else int(len(alphas) - 1)
        target_alpha = float(alphas[alpha_idx])
        print(f"[StabSel] alpha={target_alpha:.4g} -> {int(n_selected[alpha_idx])} selected on full data")
    except Exception as e:
        print(f"[StabSel] Fallback to univariate: {e}")
        return rank_genes_univariate(train_df, gene_cols)

    selection_counts = np.zeros(len(gene_cols), dtype=int)
    failed = 0
    for i in range(n_subsamples):
        idx = rng.choice(n, size=n_sub, replace=False)
        try:
            cox = CoxnetSurvivalAnalysis(
                alphas=[target_alpha], l1_ratio=1.0,
                fit_baseline_model=False, max_iter=300,
            )
            cox.fit(X_scaled[idx], y[idx])
            selection_counts += (np.abs(cox.coef_.ravel()) > 1e-8).astype(int)
        except Exception:
            failed += 1

    selection_probs = selection_counts / max(n_subsamples - failed, 1)
    n_stable = int((selection_probs >= stability_threshold).sum())
    print(
        f"[StabSel] {n_stable}/{len(gene_cols)} genes stable (>={stability_threshold*100:.0f}%), "
        f"{failed} subsamples failed"
    )
    order = np.argsort(-selection_probs)
    return [gene_cols[int(i)] for i in order]


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
    X_base = df[base_cols].to_numpy(dtype=np.float32)
    A = df[act_col].to_numpy(dtype=np.float32).reshape(-1, 1)

    names = list(base_cols)
    blocks = [X_base]
    if len(inter_genes) > 0:
        X_int = df[list(inter_genes)].to_numpy(dtype=np.float32) * A
        blocks.append(X_int)
        names += [f"{g}*ACT" for g in inter_genes]
        if int(dup_inter) > 1:
            for d in range(1, int(dup_inter)):
                blocks.append(X_int.copy())
                names += [f"{g}*ACT#dup{d}" for g in inter_genes]
    X = np.concatenate(blocks, axis=1) if len(blocks) > 1 else X_base
    return X.astype(np.float32), names


def make_dmatrix(
    X: np.ndarray,
    time_vals: np.ndarray,
    events: np.ndarray,
    weight: np.ndarray | None = None,
    feature_names: list[str] | None = None,
) -> xgb.DMatrix:
    y_signed = prepare.pack_cox_labels(time_vals, events)
    return xgb.DMatrix(X.astype(np.float32), label=y_signed, weight=weight, feature_names=feature_names)


def xgb_cindex_eval(predt: np.ndarray, dtrain: xgb.DMatrix) -> tuple[str, float]:
    y_signed = dtrain.get_label()
    times = np.abs(y_signed)
    events = (y_signed > 0.0).astype(int)
    return "cindex", prepare.cindex(predt, times, events)


def train_xgb_cox(
    dtrain: xgb.DMatrix,
    dvalid: xgb.DMatrix,
    params: dict,
    num_boost_round: int,
    early_stopping_rounds: int,
) -> tuple[xgb.Booster, dict]:
    evals_result = {}
    booster = xgb.train(
        params=params,
        dtrain=dtrain,
        num_boost_round=int(num_boost_round),
        evals=[(dtrain, "train"), (dvalid, "valid")],
        custom_metric=xgb_cindex_eval,
        evals_result=evals_result,
        early_stopping_rounds=int(early_stopping_rounds),
        maximize=True,
        verbose_eval=False,
    )
    return booster, evals_result


def build_trial_mats_for_splits(
    train_fit_df: pd.DataFrame,
    valid_eval_df: pd.DataFrame,
    w_fit: np.ndarray | None,
    k_main: int,
    k_int: int,
    dup_inter: int,
    gene_rank: list[str],
    clin_cols: list[str],
) -> tuple[xgb.DMatrix, xgb.DMatrix, list[str], list[str], list[str]]:
    genes_main = gene_rank[:k_main]
    genes_inter = genes_main[:k_int]
    Xtr, feat_names = build_features_with_interactions(
        train_fit_df, genes_main, genes_inter, dup_inter=dup_inter, clin_cols=clin_cols
    )
    Xva, _ = build_features_with_interactions(
        valid_eval_df, genes_main, genes_inter, dup_inter=dup_inter, clin_cols=clin_cols
    )
    dtr = make_dmatrix(
        Xtr,
        train_fit_df["OS_MONTHS"].values,
        train_fit_df["OS_STATUS"].values,
        weight=w_fit,
        feature_names=feat_names,
    )
    dva = make_dmatrix(
        Xva,
        valid_eval_df["OS_MONTHS"].values,
        valid_eval_df["OS_STATUS"].values,
        weight=None,
        feature_names=feat_names,
    )
    return dtr, dva, feat_names, genes_main, genes_inter


def suggest_hparams(
    trial: optuna.Trial,
    feat_budget: int,
    clin_cols: list[str],
    max_genes: int,
) -> tuple[int, int, int, dict, int, int]:
    max_nonclin = max(8, feat_budget - len(clin_cols))
    max_main = max(1, max_nonclin - 1)

    base_main = [16, 32, 64, 96, 128, 192, 256, 384, 512, max_genes]
    topk_main_choices = tuple(sorted({k for k in base_main if 1 <= k <= min(max_genes, max_main)}))
    if len(topk_main_choices) == 0:
        topk_main_choices = (min(max_genes, max_main),)
    top_k_genes = int(trial.suggest_categorical("top_k_genes", topk_main_choices))

    k_main = int(min(top_k_genes, max_main))
    k_int_cap = int(max(0, min(k_main, max_nonclin - k_main)))
    base_inter = [0, 8, 16, 32, 64]
    topk_inter_choices = tuple(sorted({k for k in base_inter if 0 <= k <= k_int_cap}))
    if len(topk_inter_choices) == 0:
        topk_inter_choices = (0,)
    top_k_inter = int(trial.suggest_categorical(f"top_k_inter_cap_{k_int_cap}", topk_inter_choices))
    k_int = int(min(top_k_inter, k_int_cap))
    dup_inter = 1

    params = {
        "objective": "survival:cox",
        "booster": "gbtree",
        "tree_method": "hist",
        "disable_default_eval_metric": True,
        "eta": trial.suggest_float("eta", 0.01, 0.12, log=True),
        "max_depth": trial.suggest_int("max_depth", 3, 6),
        "min_child_weight": trial.suggest_float("min_child_weight", 2.0, 30.0, log=True),
        "gamma": trial.suggest_float("gamma", 0.0, 6.0),
        "subsample": trial.suggest_float("subsample", 0.65, 0.95),
        "colsample_bytree": trial.suggest_float("colsample_bytree", 0.5, 0.9),
        "colsample_bylevel": trial.suggest_float("colsample_bylevel", 0.6, 1.0),
        "colsample_bynode": trial.suggest_float("colsample_bynode", 0.6, 1.0),
        "reg_lambda": trial.suggest_float("reg_lambda", 3.0, 60.0, log=True),
        "reg_alpha": trial.suggest_float("reg_alpha", 0.0, 3.0),
        "max_bin": trial.suggest_int("max_bin", 128, 512),
        "seed": 42,
        "nthread": -1,
    }
    num_boost_round = trial.suggest_int("num_boost_round", 600, 3200, step=100)
    early_stopping_rounds = trial.suggest_int("early_stopping_rounds", 75, 175, step=25)
    return k_main, k_int, dup_inter, params, num_boost_round, early_stopping_rounds


def _select_pareto_compromise(study: optuna.Study) -> optuna.trial.FrozenTrial:
    candidates = [
        t for t in study.best_trials
        if t.state == optuna.trial.TrialState.COMPLETE and t.values is not None
    ]
    if not candidates:
        candidates = [
            t for t in study.trials
            if t.state == optuna.trial.TrialState.COMPLETE and t.values is not None
        ]
    if not candidates:
        raise RuntimeError("No completed multi-objective trials found.")

    ci_vals = np.array([float(t.values[0]) for t in candidates], dtype=float)
    rmst_vals = np.array([float(t.values[1]) for t in candidates], dtype=float)
    ci_min, ci_max = float(ci_vals.min()), float(ci_vals.max())
    rmst_min, rmst_max = float(rmst_vals.min()), float(rmst_vals.max())

    def norm(x: float, lo: float, hi: float) -> float:
        return 0.0 if hi <= lo else (x - lo) / (hi - lo)

    scored = []
    for t in candidates:
        score = norm(float(t.values[0]), ci_min, ci_max) + norm(float(t.values[1]), rmst_min, rmst_max)
        scored.append((score, t))
    scored.sort(key=lambda z: z[0], reverse=True)
    return scored[0][1]


def _xgb_params_from_trial_attrs(chosen: optuna.trial.FrozenTrial, seed: int = 7) -> tuple[dict, int, int]:
    params = {
        "objective": "survival:cox",
        "booster": "gbtree",
        "tree_method": "hist",
        "disable_default_eval_metric": True,
        "seed": int(seed),
        "nthread": -1,
    }
    params.update(
        {
            k: v
            for k, v in chosen.params.items()
            if k
            in {
                "eta",
                "max_depth",
                "min_child_weight",
                "gamma",
                "subsample",
                "colsample_bytree",
                "colsample_bylevel",
                "colsample_bynode",
                "reg_lambda",
                "reg_alpha",
                "max_bin",
            }
        }
    )
    return params, int(chosen.params.get("num_boost_round", 1600)), int(chosen.params.get("early_stopping_rounds", 125))


def fit_model_from_metadata(
    metadata: dict,
    train_df: pd.DataFrame,
    valid_df: pd.DataFrame,
    refit_train_valid: bool = False,
) -> tuple[xgb.Booster, int, list[str], list[str], list[str], list[str]]:
    clin_cols = list(metadata["clin_cols"])
    clin_pretx = [c for c in clin_cols if c != "Adjuvant Chemo"]
    genes_main = list(metadata["genes_main"])
    genes_inter = list(metadata["genes_inter"])
    dup_inter = int(metadata.get("dup_inter", 1))
    params = dict(metadata["xgb_params"])
    num_boost_round = int(metadata["num_boost_round"])
    early_stopping_rounds = int(metadata["early_stopping_rounds"])

    if refit_train_valid:
        trainval_df = pd.concat([train_df, valid_df], axis=0, ignore_index=True)
        trainval_df = trainval_df.sort_values(
            by=["OS_MONTHS", "OS_STATUS"], ascending=[False, False]
        ).reset_index(drop=True)
        X_all, feat_names = build_features_with_interactions(
            trainval_df, genes_main, genes_inter, dup_inter=dup_inter, clin_cols=clin_cols
        )
        y_all = prepare.pack_cox_labels(trainval_df["OS_MONTHS"].values, trainval_df["OS_STATUS"].values)
        w_all, _, _ = prepare.compute_iptw(trainval_df, covariate_cols=clin_pretx)
        idx = np.arange(len(y_all))
        tr_idx, va_idx = train_test_split(
            idx,
            test_size=0.25,
            random_state=42,
            stratify=(y_all > 0).astype(int),
        )
        dtr = xgb.DMatrix(X_all[tr_idx], label=y_all[tr_idx], weight=w_all[tr_idx], feature_names=feat_names)
        dva = xgb.DMatrix(X_all[va_idx], label=y_all[va_idx], weight=w_all[va_idx], feature_names=feat_names)
    else:
        w_tr, _, _ = prepare.compute_iptw(train_df, covariate_cols=clin_pretx)
        dtr, dva, feat_names, _, _ = build_trial_mats_for_splits(
            train_df,
            valid_df,
            w_fit=w_tr,
            k_main=len(genes_main),
            k_int=len(genes_inter),
            dup_inter=dup_inter,
            gene_rank=genes_main,
            clin_cols=clin_cols,
        )

    booster, _ = train_xgb_cox(dtr, dva, params, num_boost_round, early_stopping_rounds)
    best_ntree = booster.best_iteration + 1 if booster.best_iteration is not None else num_boost_round
    booster_best = prepare.slice_booster_to_best_iteration(booster, best_ntree)
    return booster_best, int(best_ntree), feat_names, genes_main, genes_inter, clin_cols


def _evaluate_tlearner_on_valid(
    booster0: xgb.Booster,
    booster1: xgb.Booster,
    valid_df: pd.DataFrame,
    feat_names: list[str],
    best_ntree0: int,
    best_ntree1: int,
    tau: float = 60,
) -> dict:
    """Evaluate T-learner pair: CI from model_0 risk, RMST from counterfactual recs."""
    from lifelines import KaplanMeierFitter
    from lifelines.utils import restricted_mean_survival_time

    X_val = valid_df[feat_names].to_numpy(dtype=np.float32)
    d = xgb.DMatrix(X_val, feature_names=feat_names)
    risk_0 = booster0.predict(d, iteration_range=(0, int(best_ntree0)), output_margin=True)
    risk_1 = booster1.predict(d, iteration_range=(0, int(best_ntree1)), output_margin=True)

    val_ci = prepare.cindex(
        risk_0,
        valid_df["OS_MONTHS"].to_numpy(dtype=float),
        valid_df["OS_STATUS"].to_numpy(dtype=int),
    )

    model_rec = (risk_1 < risk_0).astype(int)
    actual = valid_df["Adjuvant Chemo"].astype(int).to_numpy()
    alignment = actual == model_rec
    if int(alignment.sum()) == 0 or int((~alignment).sum()) == 0:
        val_rmst_diff = 0.0
    else:
        kmf_a = KaplanMeierFitter()
        kmf_b = KaplanMeierFitter()
        kmf_a.fit(valid_df.loc[alignment, "OS_MONTHS"], event_observed=valid_df.loc[alignment, "OS_STATUS"])
        kmf_b.fit(valid_df.loc[~alignment, "OS_MONTHS"], event_observed=valid_df.loc[~alignment, "OS_STATUS"])
        val_rmst_diff = float(
            restricted_mean_survival_time(kmf_a, t=tau) - restricted_mean_survival_time(kmf_b, t=tau)
        )

    n_rec_act = int(model_rec.sum())
    print(f"  [T-learner] Recommends ACT for {n_rec_act}/{len(model_rec)} val patients")
    return {"val_ci": float(val_ci), "val_rmst_diff": float(val_rmst_diff), "n_features": len(feat_names)}


def _evaluate_ensemble_on_valid(
    booster_s: xgb.Booster,
    feat_names_s: list[str],
    best_ntree_s: int,
    booster0: xgb.Booster,
    booster1: xgb.Booster,
    feat_names_t: list[str],
    best_ntree0: int,
    best_ntree1: int,
    valid_df: pd.DataFrame,
    tau: float = 60,
) -> dict:
    """Ensemble: CI from S-learner risk; RMST from averaged ITE_S and ITE_T."""
    from lifelines import KaplanMeierFitter
    from lifelines.utils import restricted_mean_survival_time

    # S-learner: CI on actual treatment assignment
    X_obs = prepare.build_matrix_from_feature_names(valid_df, feat_names_s)
    risk_s_obs = prepare.predict_xgb_risk(booster_s, X_obs, feat_names_s, best_ntree_s)
    val_ci = prepare.cindex(
        risk_s_obs,
        valid_df["OS_MONTHS"].to_numpy(dtype=float),
        valid_df["OS_STATUS"].to_numpy(dtype=int),
    )

    # S-learner counterfactual ITE
    vdf1 = valid_df.copy(); vdf1["Adjuvant Chemo"] = 1
    vdf0 = valid_df.copy(); vdf0["Adjuvant Chemo"] = 0
    X_s1 = prepare.build_matrix_from_feature_names(vdf1, feat_names_s)
    X_s0 = prepare.build_matrix_from_feature_names(vdf0, feat_names_s)
    ite_s = (
        prepare.predict_xgb_risk(booster_s, X_s1, feat_names_s, best_ntree_s)
        - prepare.predict_xgb_risk(booster_s, X_s0, feat_names_s, best_ntree_s)
    )

    # T-learner ITE
    X_t = valid_df[feat_names_t].to_numpy(dtype=np.float32)
    d_t = xgb.DMatrix(X_t, feature_names=feat_names_t)
    risk_t0 = booster0.predict(d_t, iteration_range=(0, int(best_ntree0)), output_margin=True)
    risk_t1 = booster1.predict(d_t, iteration_range=(0, int(best_ntree1)), output_margin=True)
    ite_t = risk_t1 - risk_t0

    # Ensemble recommendation: ACT if avg ITE < 0
    ite_ens = (ite_s + ite_t) / 2
    model_rec = (ite_ens < 0).astype(int)

    actual = valid_df["Adjuvant Chemo"].astype(int).to_numpy()
    alignment = actual == model_rec
    if int(alignment.sum()) == 0 or int((~alignment).sum()) == 0:
        val_rmst_diff = 0.0
    else:
        kmf_a = KaplanMeierFitter()
        kmf_b = KaplanMeierFitter()
        kmf_a.fit(valid_df.loc[alignment, "OS_MONTHS"], event_observed=valid_df.loc[alignment, "OS_STATUS"])
        kmf_b.fit(valid_df.loc[~alignment, "OS_MONTHS"], event_observed=valid_df.loc[~alignment, "OS_STATUS"])
        val_rmst_diff = float(
            restricted_mean_survival_time(kmf_a, t=tau) - restricted_mean_survival_time(kmf_b, t=tau)
        )

    n_rec_act = int(model_rec.sum())
    print(f"  [Ensemble] Rec ACT {n_rec_act}/{len(model_rec)}, ITE_S mean={ite_s.mean():.3f}, ITE_T mean={ite_t.mean():.3f}")
    n_feats = len(feat_names_s) + len(feat_names_t)
    return {"val_ci": float(val_ci), "val_rmst_diff": float(val_rmst_diff), "n_features": n_feats}


def _save_run_artifacts(
    result: dict,
    booster: xgb.Booster,
    metadata: dict,
    feature_names: list[str],
    genes_main: list[str],
    genes_inter: list[str],
) -> Path:
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    run_id = time.strftime("smoke_%Y%m%d_%H%M%S")
    run_dir = RUNS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    booster.save_model(run_dir / "xgb_model.json")
    (run_dir / "feature_names.txt").write_text("\n".join(feature_names) + "\n")
    (run_dir / "genes_main.txt").write_text("\n".join(genes_main) + "\n")
    (run_dir / "genes_inter.txt").write_text("\n".join(genes_inter) + "\n")
    (run_dir / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    (run_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return run_dir


def _save_tlearner_artifacts(
    result: dict,
    booster0: xgb.Booster,
    booster1: xgb.Booster,
    metadata: dict,
    feat_names: list[str],
    genes_main: list[str],
) -> Path:
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    run_id = time.strftime("smoke_%Y%m%d_%H%M%S")
    run_dir = RUNS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    booster0.save_model(run_dir / "xgb_model_arm0.json")
    booster1.save_model(run_dir / "xgb_model_arm1.json")
    (run_dir / "feature_names.txt").write_text("\n".join(feat_names) + "\n")
    (run_dir / "genes_main.txt").write_text("\n".join(genes_main) + "\n")
    (run_dir / "genes_inter.txt").write_text("\n")
    (run_dir / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    (run_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return run_dir


def _save_ensemble_artifacts(
    result: dict,
    booster_s: xgb.Booster,
    booster0: xgb.Booster,
    booster1: xgb.Booster,
    metadata: dict,
    feat_names_s: list[str],
    feat_names_t: list[str],
    genes_main: list[str],
    genes_inter: list[str],
) -> Path:
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    run_id = time.strftime("smoke_%Y%m%d_%H%M%S")
    run_dir = RUNS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    booster_s.save_model(run_dir / "xgb_model_s.json")
    booster0.save_model(run_dir / "xgb_model_arm0.json")
    booster1.save_model(run_dir / "xgb_model_arm1.json")
    (run_dir / "feat_names_s.txt").write_text("\n".join(feat_names_s) + "\n")
    (run_dir / "feat_names_t.txt").write_text("\n".join(feat_names_t) + "\n")
    (run_dir / "genes_main.txt").write_text("\n".join(genes_main) + "\n")
    (run_dir / "genes_inter.txt").write_text("\n".join(genes_inter) + "\n")
    (run_dir / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    (run_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return run_dir


def run(n_trials: int = DEFAULT_N_TRIALS, bootstrap_n: int = DEFAULT_BOOTSTRAPS, save_artifacts: bool = True) -> dict:
    start = time.time()
    train_df, valid_df = prepare.load_train_valid()
    clin_cols, clin_pretx, gene_feats = prepare.clinical_and_gene_columns(train_df, valid_df)
    gene_rank = stability_selection_genes(train_df, gene_feats)
    max_genes = len(gene_rank)

    n_events_arm1 = int(train_df.loc[train_df["Adjuvant Chemo"].astype(int) == 1, "OS_STATUS"].sum())
    n_events_arm0 = int(train_df.loc[train_df["Adjuvant Chemo"].astype(int) == 0, "OS_STATUS"].sum())
    feat_budget = max(24, int(FEAT_EVENT_FRACTION * (n_events_arm1 + n_events_arm0)))
    print(f"[Gene Ranking] Ranked {max_genes} genes on TRAIN")
    print(f"[Ensemble] arm0 events={n_events_arm0}, arm1 events={n_events_arm1}, feat_budget={feat_budget}")
    print(f"Starting S+T Ensemble bootstrap optimization: {n_trials} trials x {bootstrap_n} bootstraps/trial")

    def objective(trial: optuna.Trial) -> tuple[float, float]:
        # S-learner uses clin_cols (includes ACT) + genes + gene*ACT interactions
        k_main, k_int, dup_inter, params, num_boost_round, esr = suggest_hparams(
            trial, feat_budget=feat_budget, clin_cols=clin_cols, max_genes=max_genes
        )
        genes_main_s = gene_rank[:k_main]
        genes_inter_s = genes_main_s[:k_int]
        # T-learner uses clin_pretx (no ACT) + same genes, no interactions
        feat_names_t = list(clin_pretx) + list(genes_main_s)

        boot_val_cis: list[float] = []
        boot_rmst_diffs: list[float] = []
        boot_best_ntrees_s: list[int] = []
        boot_best_ntrees_t: list[int] = []

        for b in range(int(bootstrap_n)):
            boot_seed = BOOTSTRAP_BASE_SEED + trial.number * 1000 + b
            boot_df = prepare.bootstrap_resample_df(
                train_df, seed=boot_seed, require_two_arms=True, require_event=True,
            )
            boot_arm0 = boot_df[boot_df["Adjuvant Chemo"].astype(int) == 0].reset_index(drop=True)
            boot_arm1 = boot_df[boot_df["Adjuvant Chemo"].astype(int) == 1].reset_index(drop=True)

            w_full, _, _ = prepare.compute_iptw(boot_df, covariate_cols=clin_pretx, act_col="Adjuvant Chemo")
            act_vals = boot_df["Adjuvant Chemo"].astype(int).values
            w_arm0 = w_full[act_vals == 0]
            w_arm1 = w_full[act_vals == 1]

            # S-learner: full boot_df with ACT feature + interactions, IPTW-weighted
            X_s_tr, feat_names_s = build_features_with_interactions(
                boot_df, genes_main_s, genes_inter_s, dup_inter=dup_inter, clin_cols=clin_cols
            )
            X_s_va, _ = build_features_with_interactions(
                valid_df, genes_main_s, genes_inter_s, dup_inter=dup_inter, clin_cols=clin_cols
            )
            p_s = dict(params); p_s["seed"] = int(boot_seed)
            dtr_s = make_dmatrix(X_s_tr, boot_df["OS_MONTHS"].values, boot_df["OS_STATUS"].values, w_full, feat_names_s)
            dva_s = make_dmatrix(X_s_va, valid_df["OS_MONTHS"].values, valid_df["OS_STATUS"].values, None, feat_names_s)
            booster_s, _ = train_xgb_cox(dtr_s, dva_s, p_s, num_boost_round, esr)
            best_s = booster_s.best_iteration + 1 if booster_s.best_iteration is not None else num_boost_round

            # T-learner: arm0/arm1 separately, no ACT feature
            X_arm0 = boot_arm0[feat_names_t].to_numpy(dtype=np.float32)
            X_arm1 = boot_arm1[feat_names_t].to_numpy(dtype=np.float32)
            X_val_t = valid_df[feat_names_t].to_numpy(dtype=np.float32)
            p0 = dict(params); p0["seed"] = int(boot_seed) + 1
            p1 = dict(params); p1["seed"] = int(boot_seed) + 2
            dtr0 = make_dmatrix(X_arm0, boot_arm0["OS_MONTHS"].values, boot_arm0["OS_STATUS"].values, w_arm0, feat_names_t)
            dtr1 = make_dmatrix(X_arm1, boot_arm1["OS_MONTHS"].values, boot_arm1["OS_STATUS"].values, w_arm1, feat_names_t)
            dva_t = make_dmatrix(X_val_t, valid_df["OS_MONTHS"].values, valid_df["OS_STATUS"].values, None, feat_names_t)
            booster0, _ = train_xgb_cox(dtr0, dva_t, p0, num_boost_round, esr)
            booster1, _ = train_xgb_cox(dtr1, dva_t, p1, num_boost_round, esr)
            best0 = booster0.best_iteration + 1 if booster0.best_iteration is not None else num_boost_round
            best1 = booster1.best_iteration + 1 if booster1.best_iteration is not None else num_boost_round

            scored = _evaluate_ensemble_on_valid(
                booster_s, feat_names_s, best_s,
                booster0, booster1, feat_names_t, best0, best1,
                valid_df,
            )
            boot_val_cis.append(float(scored["val_ci"]))
            boot_rmst_diffs.append(float(scored["val_rmst_diff"]))
            boot_best_ntrees_s.append(int(best_s))
            boot_best_ntrees_t.append((int(best0) + int(best1)) // 2)

        median_val_ci = float(np.median(boot_val_cis))
        median_rmst_diff = float(np.median(boot_rmst_diffs))
        val_ci_se = float(np.std(boot_val_cis, ddof=1)) if len(boot_val_cis) > 1 else 0.0
        rmst_iqr = (
            float(np.percentile(boot_rmst_diffs, 75) - np.percentile(boot_rmst_diffs, 25))
            if len(boot_rmst_diffs) > 1 else 0.0
        )
        best_ntree_s_med = int(round(np.median(boot_best_ntrees_s))) if boot_best_ntrees_s else num_boost_round
        best_ntree_t_med = int(round(np.median(boot_best_ntrees_t))) if boot_best_ntrees_t else num_boost_round

        trial.set_user_attr("n_features", len(feat_names_s) + len(feat_names_t))
        trial.set_user_attr("k_main", int(k_main))
        trial.set_user_attr("k_int", int(k_int))
        trial.set_user_attr("dup_inter", int(dup_inter))
        trial.set_user_attr("best_ntree_s", int(best_ntree_s_med))
        trial.set_user_attr("best_ntree_t", int(best_ntree_t_med))
        trial.set_user_attr("val_ci_boot_se", float(val_ci_se))
        trial.set_user_attr("rmst_diff_boot_iqr", float(rmst_iqr))
        trial.set_user_attr("bootstrap_n", int(bootstrap_n))
        print(
            f"[Trial {trial.number:03d}] Boot Val CI median({bootstrap_n})={median_val_ci:.4f} "
            f"(SE={val_ci_se:.4f}), RMST diff median({bootstrap_n})={median_rmst_diff:.4f} "
            f"(IQR={rmst_iqr:.4f}), N_feats_s={len(feat_names_s)}, N_feats_t={len(feat_names_t)}, "
            f"K_main={k_main}, K_int={k_int}, N_tree_s={best_ntree_s_med}, N_tree_t={best_ntree_t_med}"
        )
        return median_val_ci, median_rmst_diff

    study = optuna.create_study(
        directions=["maximize", "maximize"],
        sampler=NSGAIISampler(seed=42),
    )
    study.optimize(objective, n_trials=int(n_trials), gc_after_trial=True)
    chosen = _select_pareto_compromise(study)
    params_fin, num_boost_round_fin, esr_fin = _xgb_params_from_trial_attrs(chosen, seed=7)

    k_main = int(chosen.user_attrs["k_main"])
    k_int = int(chosen.user_attrs.get("k_int", 0))
    dup_inter = int(chosen.user_attrs.get("dup_inter", 1))
    genes_main = gene_rank[:k_main]
    genes_inter = genes_main[:k_int]
    feat_names_t = list(clin_pretx) + list(genes_main)

    # Final S-learner on full train_df
    X_s_tr_fin, feat_names_s = build_features_with_interactions(
        train_df, genes_main, genes_inter, dup_inter=dup_inter, clin_cols=clin_cols
    )
    X_s_va_fin, _ = build_features_with_interactions(
        valid_df, genes_main, genes_inter, dup_inter=dup_inter, clin_cols=clin_cols
    )
    w_s_fin, _, _ = prepare.compute_iptw(train_df, covariate_cols=clin_pretx)
    dtr_s_fin = make_dmatrix(X_s_tr_fin, train_df["OS_MONTHS"].values, train_df["OS_STATUS"].values, w_s_fin, feat_names_s)
    dva_s_fin = make_dmatrix(X_s_va_fin, valid_df["OS_MONTHS"].values, valid_df["OS_STATUS"].values, None, feat_names_s)
    pf_s = dict(params_fin); pf_s["seed"] = 7
    booster_s_fin, _ = train_xgb_cox(dtr_s_fin, dva_s_fin, pf_s, num_boost_round_fin, esr_fin)
    best_s_fin = booster_s_fin.best_iteration + 1 if booster_s_fin.best_iteration is not None else num_boost_round_fin
    booster_s_best = prepare.slice_booster_to_best_iteration(booster_s_fin, best_s_fin)

    # Final T-learner arm0/arm1 on full train_df
    train_arm0 = train_df[train_df["Adjuvant Chemo"].astype(int) == 0].reset_index(drop=True)
    train_arm1 = train_df[train_df["Adjuvant Chemo"].astype(int) == 1].reset_index(drop=True)
    w_full_fin, _, _ = prepare.compute_iptw(train_df, covariate_cols=clin_pretx, act_col="Adjuvant Chemo")
    act_vals_tr = train_df["Adjuvant Chemo"].astype(int).values
    w_arm0_fin = w_full_fin[act_vals_tr == 0]
    w_arm1_fin = w_full_fin[act_vals_tr == 1]

    X_arm0_fin = train_arm0[feat_names_t].to_numpy(dtype=np.float32)
    X_arm1_fin = train_arm1[feat_names_t].to_numpy(dtype=np.float32)
    X_val_t_fin = valid_df[feat_names_t].to_numpy(dtype=np.float32)
    dtr0_fin = make_dmatrix(X_arm0_fin, train_arm0["OS_MONTHS"].values, train_arm0["OS_STATUS"].values, w_arm0_fin, feat_names_t)
    dtr1_fin = make_dmatrix(X_arm1_fin, train_arm1["OS_MONTHS"].values, train_arm1["OS_STATUS"].values, w_arm1_fin, feat_names_t)
    dva_t_fin = make_dmatrix(X_val_t_fin, valid_df["OS_MONTHS"].values, valid_df["OS_STATUS"].values, None, feat_names_t)
    pf0 = dict(params_fin); pf0["seed"] = 8
    pf1 = dict(params_fin); pf1["seed"] = 9
    booster0_fin, _ = train_xgb_cox(dtr0_fin, dva_t_fin, pf0, num_boost_round_fin, esr_fin)
    booster1_fin, _ = train_xgb_cox(dtr1_fin, dva_t_fin, pf1, num_boost_round_fin, esr_fin)
    best0_fin = booster0_fin.best_iteration + 1 if booster0_fin.best_iteration is not None else num_boost_round_fin
    best1_fin = booster1_fin.best_iteration + 1 if booster1_fin.best_iteration is not None else num_boost_round_fin
    booster0_best = prepare.slice_booster_to_best_iteration(booster0_fin, best0_fin)
    booster1_best = prepare.slice_booster_to_best_iteration(booster1_fin, best1_fin)

    valid_result = _evaluate_ensemble_on_valid(
        booster_s_best, feat_names_s, best_s_fin,
        booster0_best, booster1_best, feat_names_t, best0_fin, best1_fin,
        valid_df,
    )

    result = {
        "val_ci": float(valid_result["val_ci"]),
        "val_rmst_diff": float(valid_result["val_rmst_diff"]),
        "val_ci_se": float(chosen.user_attrs.get("val_ci_boot_se", 0.0)),
        "rmst_iqr": float(chosen.user_attrs.get("rmst_diff_boot_iqr", 0.0)),
        "n_features": int(valid_result["n_features"]),
        "k_main": int(k_main),
        "k_int": int(k_int),
        "dup_inter": int(dup_inter),
        "best_ntree_s": int(best_s_fin),
        "best_ntree_arm0": int(best0_fin),
        "best_ntree_arm1": int(best1_fin),
        "bootstrap_n": int(bootstrap_n),
        "n_trials": int(n_trials),
        "chosen_trial": int(chosen.number),
        "notes": "iter_007: s_t_ensemble — CI from S-learner; RMST from averaged S+T ITE.",
        "wall_clock_sec": round(time.time() - start, 2),
    }
    metadata = {
        "arena": "xgb_s_t_ensemble",
        "created_by": "xgb_arena/train.py",
        "sealed_test_policy": "train.py uses train/validation data only; use finalize.py manually.",
        "clin_cols": list(clin_cols),
        "clin_pretx": clin_pretx,
        "genes_main": list(genes_main),
        "genes_inter": list(genes_inter),
        "dup_inter": int(dup_inter),
        "feat_names_s": list(feat_names_s),
        "feat_names_t": list(feat_names_t),
        "xgb_params": params_fin,
        "num_boost_round": int(num_boost_round_fin),
        "early_stopping_rounds": int(esr_fin),
        "best_ntree_s": int(best_s_fin),
        "best_ntree_arm0": int(best0_fin),
        "best_ntree_arm1": int(best1_fin),
        "result": result,
        "chosen_trial_params": dict(chosen.params),
    }
    if save_artifacts:
        run_dir = _save_ensemble_artifacts(
            result, booster_s_best, booster0_best, booster1_best,
            metadata, feat_names_s, feat_names_t, genes_main, genes_inter,
        )
        result["run_dir"] = str(run_dir.relative_to(ARENA_DIR.parent))
        print(f"[Artifacts] Saved ensemble artifacts to {run_dir}")

    print("\n[Result Row]")
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


if __name__ == "__main__":
    run()
