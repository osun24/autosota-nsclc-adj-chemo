"""Agent-editable DeepSurv arena entry point."""

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
import torch
import torch.nn as nn
from optuna.samplers import NSGAIISampler
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
torch.manual_seed(42)

ARENA_DIR = Path(__file__).resolve().parent
RUNS_DIR = ARENA_DIR / "runs"
DEFAULT_N_TRIALS = int(os.environ.get("DEEPSURV_ARENA_N_TRIALS", "10"))
DEFAULT_BOOTSTRAPS = int(os.environ.get("DEEPSURV_ARENA_BOOTSTRAPS", "2"))
DEFAULT_MAX_EPOCHS = int(os.environ.get("DEEPSURV_ARENA_EPOCHS", "80"))
BOOTSTRAP_BASE_SEED = 27182
FEAT_EVENT_FRACTION = 0.35


class DeepSurvMLP(nn.Module):
    def __init__(self, in_features: int, hidden_layers: list[int], dropout: float = 0.0):
        super().__init__()
        layers: list[nn.Module] = []
        d = int(in_features)
        for h in hidden_layers:
            layers.append(nn.Linear(d, int(h)))
            layers.append(nn.ReLU())
            if dropout > 0:
                layers.append(nn.Dropout(float(dropout)))
            d = int(h)
        layers.append(nn.Linear(d, 1))
        self.model = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.model(x)


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
    if inter_genes:
        X_int = df[list(inter_genes)].to_numpy(dtype=np.float32) * A
        blocks.append(X_int)
        names += [f"{g}*ACT" for g in inter_genes]
        if int(dup_inter) > 1:
            for d in range(1, int(dup_inter)):
                blocks.append(X_int.copy())
                names += [f"{g}*ACT#dup{d}" for g in inter_genes]
    return (np.concatenate(blocks, axis=1) if len(blocks) > 1 else X_base).astype(np.float32), names


def cox_breslow_weighted(pred: torch.Tensor, event: torch.Tensor, time_vals: torch.Tensor, weight: torch.Tensor) -> torch.Tensor:
    eta = pred.reshape(-1)
    e = event.float().reshape(-1)
    t = time_vals.reshape(-1)
    w = weight.reshape(-1)
    order = torch.argsort(t, descending=True)
    t, e, eta, w = t[order], e[order], eta[order], w[order]
    exp_eta = torch.exp(torch.clamp(eta, -20, 20))
    cum_w_exp = torch.cumsum(w * exp_eta, dim=0)
    event_mask = e > 0.5
    if event_mask.sum() == 0:
        return eta.sum() * 0.0
    log_risk = torch.log(cum_w_exp + 1e-8)
    loss = -((w * event_mask.float()) * (eta - log_risk)).sum()
    return loss / (w[event_mask].sum() + 1e-8)


def train_deepsurv_model(
    X: np.ndarray,
    time_vals: np.ndarray,
    events: np.ndarray,
    weights: np.ndarray,
    params: dict,
    device: torch.device,
    max_epochs: int,
) -> DeepSurvMLP:
    torch.manual_seed(int(params.get("seed", 42)))
    model = DeepSurvMLP(X.shape[1], params["hidden_layers"], params["dropout"]).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=float(params["lr"]), weight_decay=float(params["weight_decay"]))
    x_t = torch.tensor(X, dtype=torch.float32, device=device)
    t_t = torch.tensor(time_vals.astype(np.float32), dtype=torch.float32, device=device)
    e_t = torch.tensor(events.astype(bool), dtype=torch.bool, device=device)
    w_t = torch.tensor(weights.astype(np.float32), dtype=torch.float32, device=device)
    best_state = None
    best_loss = np.inf
    patience = int(params.get("patience", 20))
    stale = 0
    for epoch in range(int(max_epochs)):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        pred = torch.clamp(model(x_t), -20, 20)
        loss = cox_breslow_weighted(pred, e_t, t_t, w_t)
        if params.get("l1", 0.0) > 0:
            first = next(m for m in model.modules() if isinstance(m, nn.Linear))
            loss = loss + float(params["l1"]) * first.weight.abs().sum()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        optimizer.step()
        loss_val = float(loss.detach().cpu().item())
        if loss_val < best_loss - 1e-5:
            best_loss = loss_val
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
        if stale >= patience:
            break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model


def suggest_hparams(trial: optuna.Trial, feat_budget: int, clin_cols: list[str], max_genes: int) -> tuple[int, int, int, dict]:
    max_nonclin = max(8, feat_budget - len(clin_cols))
    base_main = [16, 32, 64, 96, 128, max_genes]
    topk_main_choices = tuple(sorted({k for k in base_main if 1 <= k <= min(max_genes, max_nonclin)}))
    if not topk_main_choices:
        topk_main_choices = (min(max_genes, max_nonclin),)
    k_main = int(trial.suggest_categorical("top_k_genes", topk_main_choices))
    k_int_cap = int(max(0, min(k_main, max_nonclin - k_main)))
    inter_choices = tuple(sorted({k for k in [0, 8, 16, 32] if k <= k_int_cap})) or (0,)
    k_int = int(trial.suggest_categorical(f"top_k_inter_cap_{k_int_cap}", inter_choices))
    dup_inter = 1
    arch = trial.suggest_categorical("arch", ["32", "64", "64-32"])
    params = {
        "hidden_layers": [int(x) for x in arch.split("-")],
        "arch": arch,
        "dropout": trial.suggest_float("dropout", 0.0, 0.35),
        "lr": trial.suggest_float("lr", 1e-4, 3e-3, log=True),
        "weight_decay": trial.suggest_float("weight_decay", 1e-5, 1e-2, log=True),
        "l1": trial.suggest_float("l1", 0.0, 1e-5),
        "patience": 20,
        "seed": 42,
    }
    return k_main, k_int, dup_inter, params


def _select_pareto_compromise(study: optuna.Study) -> optuna.trial.FrozenTrial:
    candidates = [t for t in study.best_trials if t.state == optuna.trial.TrialState.COMPLETE and t.values is not None]
    if not candidates:
        candidates = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE and t.values is not None]
    if not candidates:
        raise RuntimeError("No completed multi-objective trials found.")
    ci_vals = np.array([float(t.values[0]) for t in candidates])
    rmst_vals = np.array([float(t.values[1]) for t in candidates])

    def norm(x, lo, hi):
        return 0.0 if hi <= lo else (x - lo) / (hi - lo)

    scored = []
    for t in candidates:
        scored.append((
            norm(float(t.values[0]), float(ci_vals.min()), float(ci_vals.max())) +
            norm(float(t.values[1]), float(rmst_vals.min()), float(rmst_vals.max())),
            t,
        ))
    return sorted(scored, key=lambda z: z[0], reverse=True)[0][1]


def _params_from_chosen(chosen: optuna.trial.FrozenTrial, seed: int = 7) -> dict:
    arch = chosen.params.get("arch", "64")
    return {
        "hidden_layers": [int(x) for x in arch.split("-")],
        "arch": arch,
        "dropout": float(chosen.params.get("dropout", 0.1)),
        "lr": float(chosen.params.get("lr", 1e-3)),
        "weight_decay": float(chosen.params.get("weight_decay", 1e-4)),
        "l1": float(chosen.params.get("l1", 0.0)),
        "patience": 20,
        "seed": int(seed),
    }


def fit_model_from_metadata(metadata: dict, train_df: pd.DataFrame, valid_df: pd.DataFrame, refit_train_valid: bool = False, device: torch.device | None = None):
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    clin_cols = list(metadata["clin_cols"])
    clin_pretx = [c for c in clin_cols if c != "Adjuvant Chemo"]
    genes_main = list(metadata["genes_main"])
    genes_inter = list(metadata.get("genes_inter", []))
    dup_inter = int(metadata.get("dup_inter", 1))
    fit_df = pd.concat([train_df, valid_df], axis=0, ignore_index=True) if refit_train_valid else train_df
    fit_df = fit_df.sort_values(["OS_MONTHS", "OS_STATUS"], ascending=[False, False]).reset_index(drop=True)
    X_raw, feature_names = build_features_with_interactions(fit_df, genes_main, genes_inter, dup_inter=dup_inter, clin_cols=clin_cols)
    transform = prepare.fit_tabular_transform(X_raw)
    X = prepare.apply_tabular_transform(X_raw, transform)
    weights, _, _ = prepare.compute_iptw(fit_df, covariate_cols=clin_pretx, act_col="Adjuvant Chemo")
    model = train_deepsurv_model(
        X,
        fit_df["OS_MONTHS"].to_numpy(np.float32),
        fit_df["OS_STATUS"].to_numpy(int),
        weights,
        dict(metadata["deepsurv_params"]),
        device,
        int(metadata.get("max_epochs", DEFAULT_MAX_EPOCHS)),
    )
    return model, transform, feature_names, genes_main, genes_inter, clin_cols


def _save_run_artifacts(result: dict, model: nn.Module, metadata: dict, transform: dict[str, np.ndarray], feature_names: list[str], genes_main: list[str], genes_inter: list[str]) -> Path:
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    run_dir = RUNS_DIR / time.strftime("smoke_%Y%m%d_%H%M%S")
    run_dir.mkdir(parents=True, exist_ok=False)
    torch.save(model.state_dict(), run_dir / "deepsurv_model.pt")
    np.savez(run_dir / "tabular_transform.npz", **transform)
    (run_dir / "feature_names.txt").write_text("\n".join(feature_names) + "\n")
    (run_dir / "genes_main.txt").write_text("\n".join(genes_main) + "\n")
    (run_dir / "genes_inter.txt").write_text("\n".join(genes_inter) + "\n")
    (run_dir / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    (run_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return run_dir


def run(n_trials: int = DEFAULT_N_TRIALS, bootstrap_n: int = DEFAULT_BOOTSTRAPS, max_epochs: int = DEFAULT_MAX_EPOCHS, save_artifacts: bool = True) -> dict:
    start = time.time()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")
    train_df, valid_df = prepare.load_train_valid()
    clin_cols, clin_pretx, gene_feats = prepare.clinical_and_gene_columns(train_df, valid_df)
    gene_rank = rank_genes_univariate(train_df, gene_feats)
    max_genes = len(gene_rank)
    feat_budget = max(24, int(FEAT_EVENT_FRACTION * int(train_df["OS_STATUS"].sum())))
    print(f"[Gene Ranking] Ranked {max_genes} genes on TRAIN")
    print(f"[Budgets] feature budget <= {feat_budget}")
    print(f"Starting bootstrap optimization: {n_trials} trials x {bootstrap_n} bootstraps x {max_epochs} epochs")

    def objective(trial: optuna.Trial):
        k_main, k_int, dup_inter, params = suggest_hparams(trial, feat_budget, clin_cols, max_genes)
        boot_cis, boot_rmsts = [], []
        n_features = None
        for b in range(int(bootstrap_n)):
            boot_seed = BOOTSTRAP_BASE_SEED + trial.number * 1000 + b
            boot_df = prepare.bootstrap_resample_df(train_df, boot_seed, require_two_arms=True, require_event=True)
            genes_main = gene_rank[:k_main]
            genes_inter = genes_main[:k_int]
            X_raw, feat_names = build_features_with_interactions(boot_df, genes_main, genes_inter, dup_inter=dup_inter, clin_cols=clin_cols)
            transform = prepare.fit_tabular_transform(X_raw)
            X = prepare.apply_tabular_transform(X_raw, transform)
            weights, _, _ = prepare.compute_iptw(boot_df, covariate_cols=clin_pretx, act_col="Adjuvant Chemo")
            params_b = dict(params)
            params_b["seed"] = int(boot_seed)
            model = train_deepsurv_model(
                X,
                boot_df["OS_MONTHS"].to_numpy(np.float32),
                boot_df["OS_STATUS"].to_numpy(int),
                weights,
                params_b,
                device,
                int(max_epochs),
            )
            scored = prepare.evaluate_on_valid(model, valid_df, feat_names, transform, device)
            boot_cis.append(float(scored["val_ci"]))
            boot_rmsts.append(float(scored["val_rmst_diff"]))
            n_features = int(scored["n_features"])
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
    params = _params_from_chosen(chosen, seed=7)
    metadata = {
        "arena": "deepsurv",
        "device_policy": "cuda if available, else cpu",
        "clin_cols": clin_cols,
        "clin_pretx": clin_pretx,
        "genes_main": genes_main,
        "genes_inter": genes_inter,
        "dup_inter": dup_inter,
        "deepsurv_params": params,
        "max_epochs": int(max_epochs),
        "chosen_trial_params": dict(chosen.params),
    }
    model, transform, feat_names, _, _, _ = fit_model_from_metadata(metadata, train_df, valid_df, refit_train_valid=False, device=device)
    valid_result = prepare.evaluate_on_valid(model, valid_df, feat_names, transform, device)
    result = {
        "val_ci": float(valid_result["val_ci"]),
        "val_rmst_diff": float(valid_result["val_rmst_diff"]),
        "val_ci_se": float(chosen.user_attrs.get("val_ci_boot_se", 0.0)),
        "rmst_iqr": float(chosen.user_attrs.get("rmst_diff_boot_iqr", 0.0)),
        "n_features": int(valid_result["n_features"]),
        "k_main": k_main,
        "k_int": k_int,
        "dup_inter": dup_inter,
        "bootstrap_n": int(bootstrap_n),
        "n_trials": int(n_trials),
        "max_epochs": int(max_epochs),
        "chosen_trial": int(chosen.number),
        "notes": "Phase 1 DeepSurv scaffold: notebook logic converted to script with reduced search budget, no test access.",
        "wall_clock_sec": round(time.time() - start, 2),
    }
    metadata["feature_names"] = feat_names
    metadata["result"] = result
    if save_artifacts:
        run_dir = _save_run_artifacts(result, model, metadata, transform, feat_names, genes_main, genes_inter)
        result["run_dir"] = str(run_dir.relative_to(ARENA_DIR.parent))
        print(f"[Artifacts] Saved DeepSurv artifacts to {run_dir}")
    print("\n[Result Row]")
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


if __name__ == "__main__":
    run()
