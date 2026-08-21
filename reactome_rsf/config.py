"""Frozen inputs and model defaults for the Reactome RSF analysis."""

from __future__ import annotations

from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_DIR.parent

# This is intentionally a literal release, not a "latest" lookup. Updating the
# collection requires changing the version, URL, filename, and verified digest.
MSIGDB_VERSION = "2026.1.Hs"
MSIGDB_COLLECTION = "C2:CP:REACTOME"
MSIGDB_GENE_ID_TYPE = "HGNC gene symbols"
MSIGDB_FILENAME = f"c2.cp.reactome.v{MSIGDB_VERSION}.symbols.gmt"
MSIGDB_URL = (
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
    f"{MSIGDB_VERSION}/{MSIGDB_FILENAME}"
)
MSIGDB_SHA256 = "5d61f289a2400cddfbb3a3353829fd2284a360bbe50f2093b566c4b7bea93341"
MSIGDB_RELEASE_CATALOG_URL = (
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
    "msigdb_releases.json"
)
MSIGDB_RELEASE_NOTES_URL = (
    "https://docs.gsea-msigdb.org/MSigDB/Release_Notes/"
    "MSigDB_2026.1.Hs/"
)

DEFAULT_DATA_DIR = PACKAGE_DIR / "data"
DEFAULT_RUNS_DIR = PACKAGE_DIR / "runs"
DEFAULT_TRAIN_CSV = REPO_ROOT / "affyfRMATrain.csv"
DEFAULT_VALID_CSV = REPO_ROOT / "affyfRMAValidation.csv"

OUTCOME_COLUMNS = ["OS_STATUS", "OS_MONTHS"]
CLINICAL_COVARIATES = [
    "Adjuvant Chemo",
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

DEFAULT_RSF_PARAMS = {
    "n_estimators": 500,
    "min_samples_split": 20,
    "min_samples_leaf": 10,
    "max_features": "sqrt",
    "max_depth": None,
    "n_jobs": -1,
    "random_state": 42,
    "oob_score": False,
    "low_memory": True,
}
