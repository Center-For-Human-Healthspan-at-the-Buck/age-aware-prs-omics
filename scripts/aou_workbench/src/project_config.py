"""Repository-relative configuration for the archived All of Us workflow.

The original analysis ran in an authorized All of Us Controlled Tier workspace.
Set GOOGLE_PROJECT and WORKSPACE_CDR through environment variables before running.
"""
from __future__ import annotations
from pathlib import Path
import os

PROJECT_ROOT = Path(__file__).resolve().parents[3]
WORKSPACE_ROOT = Path(os.getenv("WORKSPACE_ROOT", PROJECT_ROOT / "aou_workspace")).expanduser()
CONFIG_DIR = PROJECT_ROOT / "scripts" / "aou_workbench" / "config"
TARGET_DIR = WORKSPACE_ROOT / "targets"
SHARED_DIR = WORKSPACE_ROOT / "shared"
PROTEOMICS_DIR = WORKSPACE_ROOT / "proteomics"
CHEMISTRY_DIR = WORKSPACE_ROOT / "chemistries"
RESULTS_DIR = WORKSPACE_ROOT / "results"
PLOTS_DIR = WORKSPACE_ROOT / "plots"
LOG_DIR = WORKSPACE_ROOT / "logs"
WORKING_DIR = WORKSPACE_ROOT / "working"
PREP = WORKSPACE_ROOT / "prepared"

for path in [WORKSPACE_ROOT, TARGET_DIR, SHARED_DIR, PROTEOMICS_DIR, CHEMISTRY_DIR,
             RESULTS_DIR, PLOTS_DIR, LOG_DIR, WORKING_DIR, PREP]:
    path.mkdir(parents=True, exist_ok=True)

GOOGLE_PROJECT = os.getenv("GOOGLE_PROJECT", "")
CDR = os.getenv("WORKSPACE_CDR", os.getenv("CDR", ""))
PRIMARY_ANALYSIS_NAME = os.getenv("AOU_PRIMARY_ANALYSIS_NAME", "primary_main_effect")
AGE_CUTOFF = float(os.getenv("AOU_AGE_CUTOFF", "45"))
N_GENETIC_PCS = int(os.getenv("AOU_N_GENETIC_PCS", "8"))
PC_COLUMNS = [f"PC{i}" for i in range(1, N_GENETIC_PCS + 1)]
MIN_CASES_CHEMISTRY = int(os.getenv("AOU_MIN_CASES_CHEMISTRY", "20"))
MIN_CONTROLS_CHEMISTRY = int(os.getenv("AOU_MIN_CONTROLS_CHEMISTRY", "100"))
MIN_CASES_PROTEOMICS = int(os.getenv("AOU_MIN_CASES_PROTEOMICS", "10"))
MIN_CONTROLS_PROTEOMICS = int(os.getenv("AOU_MIN_CONTROLS_PROTEOMICS", "50"))
FORCE_REBUILD = os.getenv("AOU_FORCE_REBUILD", "false").lower() in {"1", "true", "yes"}

DATA_ROOT = PROJECT_ROOT / "data" / "aou"
LINEAR_PROTEIN_INPUT = DATA_ROOT / "targets" / "ols_proteins.tsv"
LINEAR_CHEMISTRY_INPUT = DATA_ROOT / "targets" / "ols_chemistries.tsv"
VAE_PROTEIN_INPUT_CANDIDATES = [DATA_ROOT / "targets" / "vae_proteins.tsv"]
VAE_CHEMISTRY_INPUT_CANDIDATES = [DATA_ROOT / "targets" / "vae_chemistries.tsv"]

# Controlled-workspace files that are not distributed in this repository.
PROTEOMICS_FILE_CANDIDATES = [
    Path(p) for p in os.getenv("AOU_PROTEOMICS_FILES", "").split(os.pathsep) if p
]
PROTEOMICS_MANIFEST_CANDIDATES = [
    Path(p) for p in os.getenv("AOU_PROTEOMICS_MANIFESTS", "").split(os.pathsep) if p
]
ANCESTRY_FILE_CANDIDATES = [
    Path(p) for p in os.getenv("AOU_ANCESTRY_FILES", "").split(os.pathsep) if p
]
LEGACY_CHEMISTRY_STRATA_CANDIDATES = [
    Path(p) for p in os.getenv("AOU_CHEMISTRY_STRATA_FILES", "").split(os.pathsep) if p
]


def first_existing(candidates, required=False):
    for candidate in candidates:
        candidate = Path(candidate)
        if candidate.exists():
            return candidate
    if required:
        raise FileNotFoundError("None of the candidate files exists: " + ", ".join(map(str, candidates)))
    return None


def bigquery_configuration_summary():
    return {
        "google_project": GOOGLE_PROJECT,
        "cdr": CDR,
        "workspace_root": str(WORKSPACE_ROOT),
        "n_genetic_pcs": N_GENETIC_PCS,
    }
