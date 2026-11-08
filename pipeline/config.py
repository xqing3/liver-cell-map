"""Shared paths, sample sheets and QC thresholds for the preprocessing notebooks.

Every notebook does `from config import *`. Change DATA_ROOT here if the data moves.
"""
import os

# Raw downloads and .h5ad outputs live in <project>/data/CSE6242_data (moved off the external drive 2026-11-07).
_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_ROOT = os.environ.get("LCM_DATA_ROOT", os.path.join(_PROJECT, "data", "CSE6242_data"))
RAW = os.path.join(DATA_ROOT, "raw")
PROCESSED = os.path.join(DATA_ROOT, "processed")
EXPORT = os.path.join(DATA_ROOT, "export")

# DEV_SUBSET: None = full run. An int = use only the first N samples of each dataset,
# for a fast logic test (same idea as DEV_SUBSET in the Leone Lab DataProcessing.ipynb).
DEV_SUBSET = None

GEO_SAMPLE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/samples/{prefix}nnn/{gsm}/suppl/{fname}"
MGI_ORTHOLOGS_URL = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"

# --- Xiong et al. 2019, mouse healthy (chow) vs NASH, GEO GSE129516 ---
XIONG = [  # (gsm, sample, condition)
    ("GSM3714747", "chow1", "healthy"), ("GSM3714748", "chow2", "healthy"), ("GSM3714749", "chow3", "healthy"),
    ("GSM3714750", "nash1", "disease"), ("GSM3714751", "nash2", "disease"), ("GSM3714752", "nash3", "disease"),
]
XIONG_FILE = "{gsm}_{sample}_filtered_gene_bc_matrices.csv.gz"

# --- Ramachandran et al. 2019, human healthy vs cirrhotic liver, GEO GSE136103 ---
# Liver samples only (blood1-4 and the two mouse_* samples in the same series are skipped).
# Each donor was sorted into CD45+ (immune) and CD45- (non-immune) fractions; A/B are repeat runs.
RAMACHANDRAN = [  # (gsm, sample)
    ("GSM4041150", "healthy1_cd45+"), ("GSM4041151", "healthy1_cd45-A"), ("GSM4041152", "healthy1_cd45-B"),
    ("GSM4041153", "healthy2_cd45+"), ("GSM4041154", "healthy2_cd45-"),
    ("GSM4041155", "healthy3_cd45+"), ("GSM4041156", "healthy3_cd45-A"), ("GSM4041157", "healthy3_cd45-B"),
    ("GSM4041158", "healthy4_cd45+"), ("GSM4041159", "healthy4_cd45-"),
    ("GSM4041160", "healthy5_cd45+"),
    ("GSM4041161", "cirrhotic1_cd45+"), ("GSM4041162", "cirrhotic1_cd45-A"), ("GSM4041163", "cirrhotic1_cd45-B"),
    ("GSM4041164", "cirrhotic2_cd45+"), ("GSM4041165", "cirrhotic2_cd45-"),
    ("GSM4041166", "cirrhotic3_cd45+"), ("GSM4041167", "cirrhotic3_cd45-"),
    ("GSM4041168", "cirrhotic4_cd45+"),
    ("GSM4041169", "cirrhotic5_cd45+"),
]
RAMACHANDRAN_FILES = ["{gsm}_{sample}_barcodes.tsv.gz", "{gsm}_{sample}_genes.tsv.gz", "{gsm}_{sample}_matrix.mtx.gz"]

# --- QC gene families per species (Leone Lab convention: flag mito, ribo, hemoglobin) ---
QC_PATTERNS = {
    "mouse": {"mt": r"^mt-", "ribo": r"^Rp[sl]", "hb": r"^Hb[ab]-"},
    "human": {"mt": r"^MT-", "ribo": r"^RP[SL]", "hb": r"^HB[AB]\d?$"},
}

# --- QC cutoffs: starting values; tune them against each notebook's percentile diagnostics ---
QC = {
    "xiong":        {"min_genes": 500, "max_genes": 6000, "max_pct_mt": 10},   # paper: >500 genes; mito cut is ours
    "ramachandran": {"min_genes": 300, "max_genes": 6000, "max_pct_mt": 30},   # paper: >=300 genes, mito <=30%
}
MIN_CELLS = 3          # min cells per gene
N_NEIGHBORS = 15       # kNN graph size (project spec)
N_PCS = 50
N_HVG = 3000
