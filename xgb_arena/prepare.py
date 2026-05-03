"""Frozen data loading, IPTW, and validation scoring utilities for XGB arena.

This module intentionally has no test-set path constant and no function that
loads test data. The autonomous loop imports this file read-only.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable
import random
import warnings

import numpy as np
import pandas as pd
import xgboost as xgb
from lifelines import KaplanMeierFitter
from lifelines.utils import restricted_mean_survival_time
from sklearn.linear_model import LogisticRegression
from sksurv.metrics import concordance_index_censored
from threadpoolctl import threadpool_limits

warnings.filterwarnings(
    "ignore",
    message="Ties in event time detected; using efron's method to handle ties.",
)

np.random.seed(42)
random.seed(42)

ARENA_DIR = Path(__file__).resolve().parent
REPO_ROOT = ARENA_DIR.parent

TRAIN_CSV = REPO_ROOT / "affyfRMATrain.csv"
VALID_CSV = REPO_ROOT / "affyfRMAValidation.csv"
GENES_CSV = REPO_ROOT / "LOOCV_Genes2.csv"

CLINICAL_VARS = [
    "Adjuvant Chemo",
    "Age",
    "IS_MALE",
    "Stage_IA",
    "Stage_IB",
    "Stage_II",
    "Stage_III",
    "Histology_Adenocarcinoma",
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
CLIN_FEATS_PRETX = [c for c in CLINICAL_VARS if c != "Adjuvant Chemo"]

BOOTSTRAP_MAX_RESAMPLE_TRIES = 256


def _assert_not_test_path(path: str | Path) -> None:
    assert "test" not in str(path).lower(), f"Test path is sealed during loop: {path}"


def load_genes_list(genes_csv: str | Path = GENES_CSV) -> list[str]:
    _assert_not_test_path(genes_csv)
    g = pd.read_csv(genes_csv)
    if "Prop" not in g.columns or "Gene" not in g.columns:
        raise ValueError("LOOCV_Genes2.csv must have columns 'Gene' and 'Prop'.")
    g["Prop"] = pd.to_numeric(g["Prop"], errors="coerce").fillna(0)
    genes = g.loc[g["Prop"] == 1, "Gene"].astype(str).tolist()
    print(f"[Genes] Selected {len(genes)} genes with Prop == 1")
    return genes


def preprocess_split(df: pd.DataFrame, clinical_vars: Iterable[str], gene_names: Iterable[str]) -> pd.DataFrame:
    df = df.copy()
    if "Adjuvant Chemo" in df.columns:
        mapped_act = df["Adjuvant Chemo"].map({"OBS": 0, "ACT": 1})
        df["Adjuvant Chemo"] = mapped_act.where(mapped_act.notna(), df["Adjuvant Chemo"])
    for col in ["Adjuvant Chemo", "IS_MALE"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
    if "OS_STATUS" in df.columns:
        df["OS_STATUS"] = pd.to_numeric(df["OS_STATUS"], errors="coerce").fillna(0).astype(int)
    if "OS_MONTHS" in df.columns:
        df["OS_MONTHS"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce").astype(float)

    keep_cols = [c for c in clinical_vars if c in df.columns] + [g for g in gene_names if g in df.columns]
    cols = ["OS_STATUS", "OS_MONTHS"] + keep_cols
    return df[cols].copy()


def load_train_valid(
    train_csv: str | Path = TRAIN_CSV,
    valid_csv: str | Path = VALID_CSV,
    genes_csv: str | Path = GENES_CSV,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    _assert_not_test_path(train_csv)
    _assert_not_test_path(valid_csv)
    gene_list = load_genes_list(genes_csv)

    train_raw = pd.read_csv(train_csv)
    valid_raw = pd.read_csv(valid_csv)

    print("Train OS_STATUS value counts:")
    print(train_raw["OS_STATUS"].value_counts())
    print("Train Adjuvant Chemo value counts:")
    print(train_raw["Adjuvant Chemo"].value_counts())
    print("Valid OS_STATUS value counts:")
    print(valid_raw["OS_STATUS"].value_counts())
    print("Valid Adjuvant Chemo value counts:")
    print(valid_raw["Adjuvant Chemo"].value_counts())

    train_df = preprocess_split(train_raw, CLINICAL_VARS, gene_list)
    valid_df = preprocess_split(valid_raw, CLINICAL_VARS, gene_list)

    feat_candidates = [c for c in (CLINICAL_VARS + gene_list) if c in train_df.columns and c in valid_df.columns]
    keep_cols = ["OS_STATUS", "OS_MONTHS"] + feat_candidates
    train_df = train_df[keep_cols].sort_values(
        by=["OS_MONTHS", "OS_STATUS"], ascending=[False, False]
    ).reset_index(drop=True)
    valid_df = valid_df[keep_cols].sort_values(
        by=["OS_MONTHS", "OS_STATUS"], ascending=[False, False]
    ).reset_index(drop=True)
    return train_df, valid_df


def clinical_and_gene_columns(train_df: pd.DataFrame, valid_df: pd.DataFrame) -> tuple[list[str], list[str], list[str]]:
    feat_candidates = [c for c in train_df.columns if c in valid_df.columns and c not in {"OS_STATUS", "OS_MONTHS"}]
    clin_cols = [c for c in CLINICAL_VARS if c in feat_candidates]
    gene_cols = [c for c in feat_candidates if c not in set(CLINICAL_VARS)]
    clin_pretx = [c for c in clin_cols if c != "Adjuvant Chemo"]
    return clin_cols, clin_pretx, gene_cols


def pack_cox_labels(time: np.ndarray, event: np.ndarray) -> np.ndarray:
    time = np.asarray(time, dtype=np.float32)
    event = np.asarray(event, dtype=int)
    return np.where(event == 1, time, -time).astype(np.float32)


def compute_iptw(
    df: pd.DataFrame,
    covariate_cols: list[str],
    act_col: str = "Adjuvant Chemo",
    ps_clip: tuple[float, float] = (0.05, 0.95),
    w_clip: tuple[float, float] = (0.1, 10.0),
    ref_prev: float | None = None,
    model: LogisticRegression | None = None,
) -> tuple[np.ndarray, LogisticRegression, float]:
    A = df[act_col].astype(int).values
    X = df[covariate_cols].astype(float).values
    if model is None:
        model = LogisticRegression(max_iter=2000, solver="lbfgs", class_weight="balanced")
        model.fit(X, A)
    ps = np.clip(model.predict_proba(X)[:, 1], ps_clip[0], ps_clip[1])
    if ref_prev is None:
        ref_prev = float(A.mean())
    w = np.where(A == 1, ref_prev / ps, (1 - ref_prev) / (1 - ps))
    return np.clip(w, w_clip[0], w_clip[1]).astype(np.float32), model, float(ref_prev)


def cindex(pred: np.ndarray, time: np.ndarray, event: np.ndarray) -> float:
    return float(concordance_index_censored(event.astype(bool), time.astype(float), pred)[0])


def build_matrix_from_feature_names(
    df: pd.DataFrame,
    feature_names: list[str],
    act_col: str = "Adjuvant Chemo",
) -> np.ndarray:
    df = df.copy()
    df[act_col] = df[act_col].astype(int)
    A = df[act_col].to_numpy(dtype=np.float32)
    cols = []
    for feat in feature_names:
        if "*ACT" in feat:
            base = feat.split("*ACT", 1)[0]
            vals = df[base].to_numpy(dtype=np.float32) * A
        else:
            vals = df[feat].to_numpy(dtype=np.float32)
        cols.append(vals.reshape(-1, 1))
    return np.hstack(cols).astype(np.float32)


def predict_xgb_risk(booster: xgb.Booster, X: np.ndarray, feature_names: list[str], best_ntree: int) -> np.ndarray:
    dmat = xgb.DMatrix(X.astype(np.float32), feature_names=feature_names)
    return booster.predict(dmat, iteration_range=(0, int(best_ntree)), output_margin=True)


def slice_booster_to_best_iteration(booster: xgb.Booster, best_ntree: int | None) -> xgb.Booster:
    if best_ntree is None:
        return booster
    try:
        return booster[: int(best_ntree)]
    except Exception:
        return booster


def compute_alignment_rmst_diff_xgb(
    booster: xgb.Booster,
    df: pd.DataFrame,
    genes_main: list[str],
    genes_inter: list[str],
    dup_inter: int,
    feature_names: list[str],
    best_ntree: int,
    clin_cols: list[str],
    tau: float = 60,
    time_col: str = "OS_MONTHS",
    event_col: str = "OS_STATUS",
) -> float:
    del genes_main, genes_inter, dup_inter, clin_cols
    df = df.copy()
    df["Adjuvant Chemo"] = df["Adjuvant Chemo"].astype(int)
    df_treated = df.copy()
    df_treated["Adjuvant Chemo"] = 1
    df_untreated = df.copy()
    df_untreated["Adjuvant Chemo"] = 0

    X_treated = build_matrix_from_feature_names(df_treated, feature_names)
    X_untreated = build_matrix_from_feature_names(df_untreated, feature_names)
    with threadpool_limits(limits=None):
        risk_treated = predict_xgb_risk(booster, X_treated, feature_names, best_ntree)
        risk_untreated = predict_xgb_risk(booster, X_untreated, feature_names, best_ntree)

    model_rec = np.where(risk_treated < risk_untreated, 1, 0)
    actual = df["Adjuvant Chemo"].to_numpy(int)
    alignment = actual == model_rec
    mask_aligned = alignment
    mask_not_aligned = ~alignment
    if int(mask_aligned.sum()) == 0 or int(mask_not_aligned.sum()) == 0:
        return 0.0

    kmf_aligned = KaplanMeierFitter()
    kmf_not_aligned = KaplanMeierFitter()
    kmf_aligned.fit(df.loc[mask_aligned, time_col], event_observed=df.loc[mask_aligned, event_col])
    kmf_not_aligned.fit(df.loc[mask_not_aligned, time_col], event_observed=df.loc[mask_not_aligned, event_col])
    rmst_aligned = float(restricted_mean_survival_time(kmf_aligned, t=tau))
    rmst_not_aligned = float(restricted_mean_survival_time(kmf_not_aligned, t=tau))
    return float(rmst_aligned - rmst_not_aligned)


def evaluate_on_valid(
    booster: xgb.Booster,
    valid_df: pd.DataFrame,
    genes_main: list[str],
    genes_inter: list[str],
    dup_inter: int,
    feature_names: list[str],
    best_ntree: int,
    clin_cols: list[str],
) -> dict[str, float | int]:
    X_valid = build_matrix_from_feature_names(valid_df, feature_names)
    pred = predict_xgb_risk(booster, X_valid, feature_names, best_ntree)
    val_ci = cindex(
        pred,
        valid_df["OS_MONTHS"].to_numpy(dtype=float),
        valid_df["OS_STATUS"].to_numpy(dtype=int),
    )
    val_rmst_diff = compute_alignment_rmst_diff_xgb(
        booster,
        valid_df,
        genes_main=genes_main,
        genes_inter=genes_inter,
        dup_inter=dup_inter,
        feature_names=feature_names,
        best_ntree=best_ntree,
        clin_cols=clin_cols,
        tau=60,
    )
    return {
        "val_ci": float(val_ci),
        "val_rmst_diff": float(val_rmst_diff),
        "n_features": int(len(feature_names)),
    }


def bootstrap_resample_df(
    df: pd.DataFrame,
    seed: int,
    require_two_arms: bool = False,
    require_event: bool = True,
) -> pd.DataFrame:
    rng = np.random.default_rng(int(seed))
    n = len(df)
    for _ in range(BOOTSTRAP_MAX_RESAMPLE_TRIES):
        idx = rng.integers(0, n, size=n)
        boot = df.iloc[idx].copy()
        if require_event and int(boot["OS_STATUS"].sum()) == 0:
            continue
        if require_two_arms and boot["Adjuvant Chemo"].nunique() < 2:
            continue
        return boot.sort_values(by=["OS_MONTHS", "OS_STATUS"], ascending=[False, False]).reset_index(drop=True)
    raise RuntimeError("Could not draw a valid bootstrap sample with events and both treatment arms.")
