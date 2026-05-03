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


def run(n_trials: int = DEFAULT_N_TRIALS, bootstrap_n: int = DEFAULT_BOOTSTRAPS, save_artifacts: bool = True) -> dict:
    start = time.time()
    train_df, valid_df = prepare.load_train_valid()
    clin_cols, clin_pretx, gene_feats = prepare.clinical_and_gene_columns(train_df, valid_df)
    gene_rank = rank_genes_univariate(train_df, gene_feats)
    max_genes = len(gene_rank)
    n_events_tr = int(train_df["OS_STATUS"].sum())
    feat_budget = max(24, int(FEAT_EVENT_FRACTION * n_events_tr))
    print(f"[Gene Ranking] Ranked {max_genes} genes on TRAIN")
    print(f"[Budgets] events(train)={n_events_tr} -> feature budget <= {feat_budget}")
    print(f"Starting bootstrap optimization: {n_trials} trials x {bootstrap_n} bootstraps/trial")

    def objective(trial: optuna.Trial) -> tuple[float, float]:
        k_main, k_int, dup_inter, params, num_boost_round, esr = suggest_hparams(
            trial, feat_budget=feat_budget, clin_cols=clin_cols, max_genes=max_genes
        )
        boot_val_cis = []
        boot_rmst_diffs = []
        boot_best_ntrees = []
        n_features = None

        for b in range(int(bootstrap_n)):
            boot_seed = BOOTSTRAP_BASE_SEED + trial.number * 1000 + b
            boot_train_df = prepare.bootstrap_resample_df(
                train_df,
                seed=boot_seed,
                require_two_arms=True,
                require_event=True,
            )
            w_tr_boot, _, _ = prepare.compute_iptw(
                boot_train_df,
                covariate_cols=clin_pretx,
                act_col="Adjuvant Chemo",
            )
            dtr, dva, feat_names, genes_main, genes_inter = build_trial_mats_for_splits(
                boot_train_df,
                valid_df,
                w_fit=w_tr_boot,
                k_main=k_main,
                k_int=k_int,
                dup_inter=dup_inter,
                gene_rank=gene_rank,
                clin_cols=clin_cols,
            )
            params_boot = dict(params)
            params_boot["seed"] = int(boot_seed)
            booster, _ = train_xgb_cox(dtr, dva, params_boot, num_boost_round, esr)
            best_ntree = booster.best_iteration + 1 if booster.best_iteration is not None else num_boost_round
            scored = prepare.evaluate_on_valid(
                booster,
                valid_df,
                genes_main=genes_main,
                genes_inter=genes_inter,
                dup_inter=dup_inter,
                feature_names=feat_names,
                best_ntree=best_ntree,
                clin_cols=clin_cols,
            )
            boot_val_cis.append(float(scored["val_ci"]))
            boot_rmst_diffs.append(float(scored["val_rmst_diff"]))
            boot_best_ntrees.append(int(best_ntree))
            n_features = int(scored["n_features"])

        median_val_ci = float(np.median(boot_val_cis))
        median_rmst_diff = float(np.median(boot_rmst_diffs))
        val_ci_se = float(np.std(boot_val_cis, ddof=1)) if len(boot_val_cis) > 1 else 0.0
        rmst_iqr = (
            float(np.percentile(boot_rmst_diffs, 75) - np.percentile(boot_rmst_diffs, 25))
            if len(boot_rmst_diffs) > 1
            else 0.0
        )
        best_ntree_median = int(round(np.median(boot_best_ntrees))) if boot_best_ntrees else int(num_boost_round)

        trial.set_user_attr("n_features", int(n_features))
        trial.set_user_attr("k_main", int(k_main))
        trial.set_user_attr("k_int", int(k_int))
        trial.set_user_attr("dup_inter", int(dup_inter))
        trial.set_user_attr("best_ntree", int(best_ntree_median))
        trial.set_user_attr("val_ci_boot_se", float(val_ci_se))
        trial.set_user_attr("rmst_diff_boot_iqr", float(rmst_iqr))
        trial.set_user_attr("bootstrap_n", int(bootstrap_n))
        print(
            f"[Trial {trial.number:03d}] Boot Val CI median({bootstrap_n})={median_val_ci:.4f} "
            f"(SE={val_ci_se:.4f}), RMST diff median({bootstrap_n})={median_rmst_diff:.4f} "
            f"(IQR={rmst_iqr:.4f}), N_feats={n_features}, K_main={k_main}, K_int={k_int}, "
            f"N_tree={best_ntree_median}"
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
    k_int = int(chosen.user_attrs["k_int"])
    dup_inter = int(chosen.user_attrs.get("dup_inter", 1))
    genes_main = gene_rank[:k_main]
    genes_inter = genes_main[:k_int]
    w_tr, _, _ = prepare.compute_iptw(train_df, covariate_cols=clin_pretx, act_col="Adjuvant Chemo")
    dtr, dva, feat_names, _, _ = build_trial_mats_for_splits(
        train_df,
        valid_df,
        w_fit=w_tr,
        k_main=k_main,
        k_int=k_int,
        dup_inter=dup_inter,
        gene_rank=gene_rank,
        clin_cols=clin_cols,
    )
    booster, _ = train_xgb_cox(dtr, dva, params_fin, num_boost_round_fin, esr_fin)
    best_ntree = booster.best_iteration + 1 if booster.best_iteration is not None else num_boost_round_fin
    booster_best = prepare.slice_booster_to_best_iteration(booster, best_ntree)
    valid_result = prepare.evaluate_on_valid(
        booster_best,
        valid_df,
        genes_main=genes_main,
        genes_inter=genes_inter,
        dup_inter=dup_inter,
        feature_names=feat_names,
        best_ntree=best_ntree,
        clin_cols=clin_cols,
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
        "best_ntree": int(best_ntree),
        "bootstrap_n": int(bootstrap_n),
        "n_trials": int(n_trials),
        "chosen_trial": int(chosen.number),
        "notes": "Phase 1 smoke-test scaffold: existing XGB logic, 10-trial/2-bootstrap default budget, no test access.",
        "wall_clock_sec": round(time.time() - start, 2),
    }
    metadata = {
        "arena": "xgb",
        "created_by": "xgb_arena/train.py",
        "sealed_test_policy": "train.py uses train/validation data only; use finalize.py manually.",
        "clin_cols": clin_cols,
        "clin_pretx": clin_pretx,
        "genes_main": genes_main,
        "genes_inter": genes_inter,
        "dup_inter": int(dup_inter),
        "feature_names": feat_names,
        "xgb_params": params_fin,
        "num_boost_round": int(num_boost_round_fin),
        "early_stopping_rounds": int(esr_fin),
        "best_ntree": int(best_ntree),
        "result": result,
        "chosen_trial_params": dict(chosen.params),
    }
    if save_artifacts:
        run_dir = _save_run_artifacts(result, booster_best, metadata, feat_names, genes_main, genes_inter)
        result["run_dir"] = str(run_dir.relative_to(ARENA_DIR.parent))
        print(f"[Artifacts] Saved smoke-test artifacts to {run_dir}")

    print("\n[Result Row]")
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


if __name__ == "__main__":
    run()
