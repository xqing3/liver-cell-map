# Guilliams human data: provisional handoff

Prepared for team handoff by Lei. This is a basic-QC and provisional-annotation checkpoint, not a final validated cell atlas.

## Scope and verification

- GEO source: GSE192742_RAW.tar.
- Selected human single-cell subset: 20 libraries from 12 donors (9 healthy, 3 disease under the team-provided donor mapping).
- 154,158 cells retained across 20 readable H5AD files; 32,643 unique gene symbols per file.
- Batch reported zero failures and count-matrix read-back verification. Handoff inspection also checked dimensions and required annotation/metadata fields in all 20 files.
- Other organisms, human nuclei/spatial libraries, and samples with unknown condition are outside this selected subset, not evidence of missing source data.

## Files

- `GSM*_provisional_annotation.h5ad`: one file per library.
- `annotation_summary.csv`: library metadata, cell totals, label changes, and local file paths.
- `celltype_counts_by_sample.csv`: majority-voting label counts.
- `annotation_failures.csv`: no failed libraries in this run.
- Notebook sources: `liver-cell-map/pipeline/00_inspect_guilliams.ipynb`, `02b_guilliams_human_qc.ipynb`, and `02c_guilliams_human_annotation.ipynb` under the same pipeline directory. Save current editor changes before sharing notebooks.

## Matrix and metadata

`.X` contains raw RNA counts, summed across duplicate gene symbols; it is not normalized or log-transformed. `var_names` are gene symbols and `var.source_gene_ids` preserves contributing Ensembl IDs. Cells retain sample-prefixed IDs, barcode, sample, donor, condition, and QC fields. TXT gene IDs were inferred from validated ordered correspondence to the CS41 H5 reference; additional non-RNA rows were excluded.

`celltypist_label` is the individual model prediction; `celltype_native` is the majority-voting prediction. `celltypist_score` is the maximum individual probability score; `celltypist_conf` is the score for the voted label (not necessarily the maximum). Neither is a calibrated annotation accuracy. Preserve both labels because voting absorbs some small predicted populations. `label_changed_by_voting` records changes.

`needs_marker_review` is a review flag, not a confirmed error: all cells in the 19 batch-annotated libraries are flagged pending review; CS41 flags the specifically questioned original or voted Neutrophils/Fibroblasts labels. A false flag in CS41 does not mean independent validation.

Only CS41 has an exploratory UMAP in this handoff. Do not interpret it as a joint embedding across libraries. Donor IDs such as H06 are study-local and should be prefixed by dataset when combining studies.

## Processing

1. RNA-only reading and per-library feature validation.
2. Provisional human cell QC: detected gene rows >=300 and <6000, mitochondrial RNA percentage <30, finite mitochondrial percentage. QC was computed before duplicate-symbol aggregation.
3. Scrublet per library: simulation ratio 1.0, expected doublet rate 0.05, random seed 0, exact neighbors to avoid a sparse-array compatibility error.
4. CellTypist 1.7.1, Healthy_Human_Liver.pkl v1: prediction on a separate log1p-normalized copy (10,000 counts per cell), with all genes retained. For voting, 2,000 highly variable genes, 30 PCs, and 15-neighbor graphs were used. Scanpy 1.11.5 and AnnData 0.11.4.
5. No cells were removed by doublet prediction or annotation score. No additional gene-frequency filter was applied to these exports.

## Unresolved items for the receiving team

- Scrublet automatic thresholds remain unapproved. Scores and automatic flags are included; these files are not confirmed singlet-only data.
- CS41 Neutrophils predictions have limited support from the inspected marker panel. No manual relabeling was performed. Rare classes and voting changes also require appropriate review.
- The CellTypist model originates from the same Guilliams study. Model agreement is not independent validation, and reference overlap should be considered in any evaluation.
- Healthy/disease status follows the team-provided donor mapping and is distinct from cell type; retain its source and confirm the grouping for downstream comparisons.
- Cross-dataset integration, final annotation adjudication, disease comparisons, and frontend outputs have not been performed.

## Supporting results

Relative to the Project folder: `guilliams_metadata/human_sample_review.csv`, `guilliams_qc/human/qc_summary.csv`, `guilliams_qc/human/basic_qc/`, and `guilliams_qc/human/doublet_review/` contain sample selection, basic QC, and doublet review outputs. Paths in summaries and source-file metadata are relative to the shared project folder, not to the CSV location. Resolve them with PROJECT_DIR / Path(value).

## Agreed stopping point

Basic QC, per-library doublet scoring, provisional annotation of all 20 selected libraries, and this handoff documentation are complete. Further biological adjudication and integration are separate follow-up work.

## Portable folder layout

Keep this layout when sharing (the enclosing folder can have any name):

```text
project/
  liver-cell-map/pipeline/     # the three Guilliams notebooks
  guilliams_metadata/
  guilliams_qc/
  guilliams_annotation/
  GSE192742_RAW.tar           # needed only to rerun inspection/basic QC
```

Run notebooks from the project folder, repository root, or pipeline folder. Setup searches parent folders for guilliams_metadata and uses relative paths. For another layout, set GUILLIAMS_PROJECT_DIR to the shared project folder. Run the setup cell after restarting the kernel so old absolute path variables are replaced. Installation logs were cleared for portability; scientific outputs were retained with displayed paths shortened. The Python environment itself must be installed on the receiving machine, not copied from .venv.
