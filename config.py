# config.py


# ============================================================
# PROJECT INFORMATION
# ============================================================

PROJECT_NAME = "Compound Triage Workbench"
VERSION = "0.1.0"


# ============================================================
# FILE PATHS
# ============================================================

INPUT_SDF = "data/input_compounds.sdf"

OUTPUT_SMILES = "data/output_compounds.smi"

INPUT_SMILES = "data/output_compounds.smi"

STANDARDIZED_SMILES = "data/standardized_compounds.smi"

# Input to filter.py
PROFILE_OUTPUT = "data/standardized_smiles_profiled.csv"

# Full filtered table
FILTERED_OUTPUT = "results/standardized_smiles_filtered.csv"

# Compounds that pass both structural and ADMET filters
FILTERED_SELECT_OUTPUT = "results/standardized_smiles_filtered_select.csv"

# SDF containing compounds selected for docking
FILTERED_COMPOUNDS_SDF = "results/filtered_compounds.sdf"

# Profile output
PROFILE_OUTPUT = "results/compound_profiles.csv"

# ADMET endpoints
ADMET_ENDPOINTS = [
    "Solubility_AqSolDB",
    "hERG",
    "Caco2_Wang",
    "BBB_Martins",
    "Clearance_Hepatocyte_AZ",
    "Clearance_Microsome_AZ",
    "PAMPA_NCATS",
    "PPBR_AZ"
]

# ============================================================
# STRUCTURE STANDARDIZATION
# ============================================================

REMOVE_SALTS = True
REMOVE_DUPLICATES = True

# Set to True only if standardization.py implements neutralization
NEUTRALIZE = False


# ============================================================
# STRUCTURAL LIABILITY FILTERS
# ============================================================

# Structural alerts
REJECT_PAINS = True
REJECT_BRENK = True
REJECT_NIH = True

# Synthetic accessibility
SA_SCORE_MAX = 6.0

# Molecular weight
MW_MAX = 500.0

# Topological polar surface area
TPSA_MAX = 140.0

# Hydrogen bond donors
HBD_MAX = 5

# Hydrogen bond acceptors
HBA_MAX = 10

# Aromatic rings
AROMATIC_RING_NUMBER_MAX = 7

# Rotatable bonds
ROTB_MAX = 10

# Lipophilicity
LOGP_MIN = 0.0
LOGP_MAX = 5.0


# ============================================================
# ADMET LIABILITY FILTERS
# ============================================================

# hERG
HERG_MAX = 0.6

# Caco2_Wang
CACO2_WANG_MAX = 6.0

# Hepatocyte clearance
CLEARANCE_HEPATOCYTE_AZ_MAX = 15.0

# Microsome clearance
CLEARANCE_MICROSOME_AZ_MAX = 50.0

# PAMPA_NCATS
PAMPA_NCATS_MIN = 0.5

# Plasma protein binding
PPBR_AZ_MAX = 90.0

# Aqueous solubility
SOLUBILITY_AQSOLDB_MIN = -4.0


# ============================================================
# WORKFLOW SETTINGS
# ============================================================

SAVE_PROFILE = True
SAVE_FILTERED = True

# Maximum number of compounds to send to docking.
# None = send all compounds that pass the filters.
MAX_DOCKING_COMPOUNDS = None