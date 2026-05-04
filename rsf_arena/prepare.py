"""Frozen data loading, IPTW, and validation scoring utilities for RSF arena."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable
import random
import warnings

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.utils import restricted_mean_survival_time
from sklearn.linear_model import LogisticRegression
from sksurv.metrics import concordance_index_censored
from threadpoolctl import threadpool_limits

warnings.filterwarnings("ignore", message="Ties in event time detected; using efron's method to handle ties.")

np.random.seed(42)
random.seed(42)

ARENA_DIR = Path(__file__).resolve().parent
REPO_ROOT = ARENA_DIR.parent
TRAIN_CSV = REPO_ROOT / "affyfRMATrain.csv"
VALID_CSV = REPO_ROOT / "affyfRMAValidation.csv"
GENES_CSV = REPO_ROOT / "LOOCV_Genes2.csv"

CLINICAL_VARS = [
    "Adjuvant Chemo", "Age", "IS_MALE",
    "Stage_IA", "Stage_IB", "Stage_II", "Stage_III",
    "Histology_Adenocarcinoma", "Histology_Large Cell Carcinoma", "Histology_Squamous Cell Carcinoma",
    "Race_African American", "Race_Asian", "Race_Caucasian", "Race_Native Hawaiian or Other Pacific Islander", "Race_Unknown",
    "Smoked?_No", "Smoked?_Unknown", "Smoked?_Yes",
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
        mapped = df["Adjuvant Chemo"].map({"OBS": 0, "ACT": 1})
        df["Adjuvant Chemo"] = mapped.where(mapped.notna(), df["Adjuvant Chemo"])
    for col in ["Adjuvant Chemo", "IS_MALE"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
    df["OS_STATUS"] = pd.to_numeric(df["OS_STATUS"], errors="coerce").fillna(0).astype(int)
    df["OS_MONTHS"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce").astype(float)
    keep_cols = [c for c in clinical_vars if c in df.columns] + [g for g in gene_names if g in df.columns]
    return df[["OS_STATUS", "OS_MONTHS"] + keep_cols].copy()


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
    train_df = train_df[keep_cols].sort_values(["OS_MONTHS", "OS_STATUS"], ascending=[False, False]).reset_index(drop=True)
    valid_df = valid_df[keep_cols].sort_values(["OS_MONTHS", "OS_STATUS"], ascending=[False, False]).reset_index(drop=True)
    return train_df, valid_df


def clinical_and_gene_columns(train_df: pd.DataFrame, valid_df: pd.DataFrame) -> tuple[list[str], list[str], list[str]]:
    feat_candidates = [c for c in train_df.columns if c in valid_df.columns and c not in {"OS_STATUS", "OS_MONTHS"}]
    clin_cols = [c for c in CLINICAL_VARS if c in feat_candidates]
    gene_cols = [c for c in feat_candidates if c not in set(CLINICAL_VARS)]
    return clin_cols, [c for c in clin_cols if c != "Adjuvant Chemo"], gene_cols


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
    return np.clip(w, w_clip[0], w_clip[1]).astype(np.float64), model, float(ref_prev)


def cindex(pred: np.ndarray, time: np.ndarray, event: np.ndarray) -> float:
    return float(concordance_index_censored(event.astype(bool), time.astype(float), pred)[0])


def build_matrix_from_feature_names(df: pd.DataFrame, feature_names: list[str], act_col: str = "Adjuvant Chemo") -> np.ndarray:
    df = df.copy()
    df[act_col] = df[act_col].astype(int)
    A = df[act_col].to_numpy(dtype=np.float64)
    cols = []
    for feat in feature_names:
        if "*ACT" in feat:
            base = feat.split("*ACT", 1)[0]
            vals = df[base].to_numpy(dtype=np.float64) * A
        else:
            vals = df[feat].to_numpy(dtype=np.float64)
        cols.append(vals.reshape(-1, 1))
    return np.hstack(cols).astype(np.float64)


def predict_rsf_risk(model, X: np.ndarray) -> np.ndarray:
    with threadpool_limits(limits=None):
        return np.asarray(model.predict(X.astype(np.float64)), dtype=float)


def compute_alignment_rmst_diff_rsf(
    model,
    df: pd.DataFrame,
    genes_main: list[str],
    genes_inter: list[str],
    dup_inter: int,
    feature_names: list[str],
    clin_cols: list[str],
    tau: float = 60,
) -> float:
    del genes_main, genes_inter, dup_inter, clin_cols
    df = df.copy()
    df["Adjuvant Chemo"] = df["Adjuvant Chemo"].astype(int)
    df_treated = df.copy()
    df_treated["Adjuvant Chemo"] = 1
    df_untreated = df.copy()
    df_untreated["Adjuvant Chemo"] = 0
    risk_treated = predict_rsf_risk(model, build_matrix_from_feature_names(df_treated, feature_names))
    risk_untreated = predict_rsf_risk(model, build_matrix_from_feature_names(df_untreated, feature_names))
    model_rec = np.where(risk_treated < risk_untreated, 1, 0)
    alignment = df["Adjuvant Chemo"].to_numpy(int) == model_rec
    if int(alignment.sum()) == 0 or int((~alignment).sum()) == 0:
        return 0.0
    km_a = KaplanMeierFitter().fit(df.loc[alignment, "OS_MONTHS"], event_observed=df.loc[alignment, "OS_STATUS"])
    km_n = KaplanMeierFitter().fit(df.loc[~alignment, "OS_MONTHS"], event_observed=df.loc[~alignment, "OS_STATUS"])
    return float(restricted_mean_survival_time(km_a, t=tau) - restricted_mean_survival_time(km_n, t=tau))


def evaluate_on_valid(
    model,
    valid_df: pd.DataFrame,
    genes_main: list[str],
    genes_inter: list[str],
    dup_inter: int,
    feature_names: list[str],
    clin_cols: list[str],
) -> dict[str, float | int]:
    X_valid = build_matrix_from_feature_names(valid_df, feature_names)
    pred = predict_rsf_risk(model, X_valid)
    val_ci = cindex(pred, valid_df["OS_MONTHS"].to_numpy(float), valid_df["OS_STATUS"].to_numpy(int))
    val_rmst_diff = compute_alignment_rmst_diff_rsf(
        model, valid_df, genes_main, genes_inter, dup_inter, feature_names, clin_cols, tau=60
    )
    return {"val_ci": float(val_ci), "val_rmst_diff": float(val_rmst_diff), "n_features": int(len(feature_names))}


def bootstrap_resample_df(df: pd.DataFrame, seed: int, require_two_arms: bool = False, require_event: bool = True) -> pd.DataFrame:
    rng = np.random.default_rng(int(seed))
    n = len(df)
    for _ in range(BOOTSTRAP_MAX_RESAMPLE_TRIES):
        idx = rng.integers(0, n, size=n)
        boot = df.iloc[idx].copy()
        if require_event and int(boot["OS_STATUS"].sum()) == 0:
            continue
        if require_two_arms and boot["Adjuvant Chemo"].nunique() < 2:
            continue
        return boot.sort_values(["OS_MONTHS", "OS_STATUS"], ascending=[False, False]).reset_index(drop=True)
    raise RuntimeError("Could not draw a valid bootstrap sample with events and both treatment arms.")
