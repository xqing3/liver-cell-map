# Preprocessing pipeline

Turns the two public liver datasets into one mouse + human cell map and the team data contract
(`cells.parquet`, `graph.npz`, `signaling.csv`) that the web tool reads. It is modeled on the Leone Lab
scanpy pipeline (`DataProcessing.ipynb` → `Annotation.ipynb`), adapted for two species.

| Notebook | What it does | Output (under `DATA_ROOT`) |
|---|---|---|
| `00_download.ipynb` | Downloads the GEO files and the MGI mouse–human gene table, checking file sizes | `raw/` |
| `01_xiong_qc.ipynb` | Mouse (GSE129516): load → QC → Scrublet → CellTypist `Healthy_Mouse_Liver` | `processed/xiong_qc.h5ad` |
| `02_ramachandran_qc.ipynb` | Human (GSE136103, 20 liver samples): same steps, CellTypist `Healthy_Human_Liver` | `processed/ramachandran_qc.h5ad` |
| `03_merge_embed.ipynb` | One-to-one gene matching, merge, HVG/PCA, Harmony (by donor), kNN (k=15), Leiden, UMAP, PAGA, shared cell-type names | `processed/merged.h5ad` |
| `04_check.ipynb` | Biology checks: species mixing, donor islands, TREM2+CD9+ macrophages in disease, cells kept per step | `figures/04_check/` |
| `05_export.ipynb` | Writes the data contract, plus a 5,000-cell subsample for the prototype | `export/full/`, `export/n5000/`, `../prototype/data_real.json` |

`config.py` holds the paths, the sample sheets and the QC cutoffs. `lcm.py` holds helpers shared by 01 and 02.
Figures and per-step cell counts go to `figures/` (in this folder, small enough to share).

## Setup

Uses the existing `scrnaseq` conda env (Python 3.11, scanpy 1.11). Two extra packages were added:
```bash
conda run -n scrnaseq python -m pip install pyarrow scikit-misc
```
Data lives in `data/CSE6242_data` at the project root (about 0.4 GB raw), found relative to `config.py`.
A teammate with a different location sets `LCM_DATA_ROOT=/their/path` before starting Jupyter.

## Run

Open the notebooks with the **Python (scrnaseq)** kernel and run them in order, or headless:
```bash
for nb in 00_download 01_xiong_qc 02_ramachandran_qc 03_merge_embed 04_check 05_export; do
  conda run -n scrnaseq jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=scrnaseq --ExecutePreprocessor.timeout=7200 $nb.ipynb
done
```
For a quick logic test, set `DEV_SUBSET = 2` in `config.py` (two samples per dataset).

## Choices worth knowing

- **QC cutoffs** follow each paper where it says: Ramachandran ≥ 300 genes and < 30% mitochondrial reads; Xiong > 500 genes.
  The upper gene cap (6,000) and Xiong's 10% mito cap are ours. Check the percentile table in 01/02 before changing them.
- **Gene matching** keeps only one-to-one mouse–human pairs from MGI, so each gene means the same thing in both species.
- **Variable genes** exclude mitochondrial, ribosomal and hemoglobin genes (Leone Lab convention). These genes stay in `.X`.
- **Harmony corrects by donor**, not by sequencing run. Ramachandran sorted each donor's cells into immune (CD45+) and
  non-immune (CD45−) runs, so runs differ in cell make-up on purpose and should not be forced together.
- **Cell-type labels** come from CellTypist models trained on *healthy* liver, so disease-only states get the nearest
  healthy label. `03` maps both label sets to one shared list (`figures/03_merge/celltype_mapping.csv`).
- **Not done here:** MELD disease scores (`meld_score` is empty) and CellPhoneDB. `signaling.csv` is a placeholder
  (mean ligand × mean receptor expression per cell type).
