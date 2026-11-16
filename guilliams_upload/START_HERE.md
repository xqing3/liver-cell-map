# Start here: Guilliams handoff

20 annotated H5AD files contain 154,158 cells from 20 libraries and 12 donors. Labels are provisional; doublet scores are provided but doublets were not removed. Read README_handoff.md before downstream use.

The three notebooks are available separately for inspection. Guilliams_notebooks_and_support.zip contains the same notebooks in their relative directory layout, sample metadata, per-cell QC tables, QC summaries, doublet review tables and figures, and the exploratory marker review.

To reconstruct the project layout, extract the ZIP into a new project directory and place the 20 H5AD files and four handoff documentation/CSV files in guilliams_annotation/human/handoff/ under that directory. The local source folder was renamed Guilliams, but the notebook output location and recorded relative paths use handoff.

The ZIP excludes raw source archives and duplicate intermediate H5AD matrices. To rerun inspection/basic QC, obtain GSE192742_RAW.tar and place it in the project directory. Recreate intermediate basic_qc and doublet_review matrices with the notebooks before rerunning annotation. Install the required Python dependencies on the receiving computer; the virtual environment is not included.

No cross-dataset integration or final biological adjudication is included. This is the agreed preparation and provisional annotation checkpoint.
