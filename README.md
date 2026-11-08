# Liver Cell Map

An interactive graph of liver cells in health and disease, from mouse and human. Each cell is a node linked
to the 15 cells with the most similar gene activity. A user selects cells in the browser, finds similar
cells by a random walk with restart, and compares disease scores and cell-to-cell signals between healthy
and diseased livers.

CSE6242 Data and Visual Analytics, Fall 2026, Team 186: Crystal Xin Qing, Ying Lei, Yalin Zhu, Jiancheng Zhu.

## Status (progress report, Nov 2026)

- Two of three datasets are merged: Xiong et al. 2019 (mouse NASH, GEO GSE129516) and Ramachandran et al.
  2019 (human cirrhosis, GEO GSE136103), 93,077 cells and 1,149,880 graph edges. The Liver Cell Atlas
  (Guilliams et al. 2022) is next.
- The browser prototype runs on all 93,077 real cells, a 5,000-cell sample, and synthetic data.
- Not done yet: MELD/Milo disease scores, and CellChat/CellPhoneDB signalling (the current signalling table
  uses 10 hand-picked ligand–receptor pairs).

## Layout

| Folder | Contents |
| --- | --- |
| `pipeline/` | Notebooks `00_download` … `05_export`, shared settings (`config.py`), helpers (`lcm.py`), QC and check figures. See `pipeline/README.md`. |
| `src/` | `make_fake_data.py` (synthetic data with a planted answer) and `export_for_web.py` (data contract → browser JSON). |
| `prototype/` | `index.html` (D3 7.9 + canvas) and its data files. |

## Run the prototype

```bash
cd prototype
python3 -m http.server 8000
```

Open `http://localhost:8000/` for all real cells, `#sample` for the 5,000-cell sample, or `#synthetic` for
fake data.

## Rebuild the data

1. Create the `scrnaseq` conda env (Python 3.11, scanpy 1.11, harmonypy, scrublet, celltypist, pyarrow,
   scikit-misc).
2. Run the notebooks in `pipeline/` in order (headless commands in `pipeline/README.md`). Raw data
   (about 0.4 GB) downloads to `data/CSE6242_data/`; set `LCM_DATA_ROOT` to put it elsewhere.
3. `05_export` writes the 5,000-cell `prototype/data_real.json` itself. Convert the full export for the
   browser with:

```bash
python3 src/export_for_web.py data/CSE6242_data/export/full prototype/data_full.json --graph-bin
```

Data files are not in this repository; all sources are public direct downloads.
