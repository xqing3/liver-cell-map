"""Helpers shared by the preprocessing notebooks (QC, diagnostics, annotation, safe save).

Patterns follow the Leone Lab pipeline (PreProcessing/DataProcessing.ipynb, Annotation.ipynb).
"""
import os

import numpy as np
import pandas as pd
import psutil
import scanpy as sc

from config import QC_PATTERNS

QC_COLS = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"]


def ram(label=""):
    used = psutil.Process().memory_info().rss / 1e9
    print(f"[RAM] {label} process uses {used:.1f} GB")
    return used


def add_qc_metrics(adata, species):
    """Flag mito / ribosomal / hemoglobin genes and compute per-cell QC metrics."""
    for key, pat in QC_PATTERNS[species].items():
        adata.var[key] = adata.var_names.str.match(pat)
    print({k: int(adata.var[k].sum()) for k in QC_PATTERNS[species]}, "genes flagged")
    sc.pp.calculate_qc_metrics(adata, qc_vars=list(QC_PATTERNS[species]), percent_top=None, log1p=False, inplace=True)


def qc_diagnostics(adata, cut):
    """Percentiles of each QC metric and how many cells each cutoff would remove (read-only)."""
    rows = []
    for col in QC_COLS:
        q = np.percentile(adata.obs[col], [1, 5, 25, 50, 75, 95, 99])
        rows.append([col] + [round(v, 1) for v in q])
    print(pd.DataFrame(rows, columns=["metric", "p1", "p5", "p25", "p50", "p75", "p95", "p99"]).to_string(index=False))
    o = adata.obs
    lo = o.n_genes_by_counts < cut["min_genes"]
    hi = o.n_genes_by_counts >= cut["max_genes"]
    mt = o.pct_counts_mt >= cut["max_pct_mt"]
    print(f"\n< {cut['min_genes']} genes: {int(lo.sum()):,} | >= {cut['max_genes']} genes: {int(hi.sum()):,} | "
          f">= {cut['max_pct_mt']}% mito: {int(mt.sum()):,} | removed by any: {int((lo | hi | mt).sum()):,} of {adata.n_obs:,}")


def apply_qc(adata, cut):
    keep = ((adata.obs.n_genes_by_counts >= cut["min_genes"]) & (adata.obs.n_genes_by_counts < cut["max_genes"])
            & (adata.obs.pct_counts_mt < cut["max_pct_mt"]))
    return adata[keep].copy()


def per_sample_table(stages):
    """stages: dict of name -> AnnData; returns cells per sample at each stage."""
    return pd.DataFrame({k: v.obs["sample"].value_counts() for k, v in stages.items()}).fillna(0).astype(int)


def remove_doublets(adata):
    """Scrublet per sample (a doublet forms within one 10x run), as in Leone Lab."""
    sc.pp.scrublet(adata, batch_key="sample", sim_doublet_ratio=1.0)
    print(adata.obs.groupby("sample", observed=True)["predicted_doublet"].mean().round(3).rename("doublet rate").to_string())
    return adata[~adata.obs["predicted_doublet"]].copy()


def celltypist_labels(adata, model_name):
    """Annotate on a log-normalized copy (CellTypist expects log1p to 1e4); returns obs columns."""
    import celltypist
    from celltypist import models

    models.download_models(model=model_name, force_update=False)
    a = adata.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    med = float(np.median(np.asarray(a.X.expm1().sum(axis=1)).ravel()))
    assert 9_000 < med < 11_000, f"X not log1p-normalized to 1e4 (median expm1 sum {med})"
    sc.pp.highly_variable_genes(a, n_top_genes=2000)
    sc.tl.pca(a, n_comps=30, mask_var="highly_variable")
    sc.pp.neighbors(a, n_neighbors=15)
    pred = celltypist.annotate(a, model=model_name, majority_voting=True)
    out = pred.predicted_labels.rename(columns={"predicted_labels": "celltypist_label",
                                                "majority_voting": "celltype_native", "over_clustering": "celltypist_cluster"})
    out["celltypist_conf"] = pred.probability_matrix.max(axis=1)
    return out[["celltypist_label", "celltype_native", "celltypist_conf"]]


def safe_write(adata, path):
    """Write, then reopen the file to prove it is readable (Leone Lab 'robust save')."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    adata.write_h5ad(path, compression="gzip")
    chk = sc.read_h5ad(path, backed="r")
    assert chk.shape == adata.shape, "written file has the wrong shape"
    chk.file.close()
    print(f"wrote {path} | {adata.n_obs:,} cells x {adata.n_vars:,} genes | {os.path.getsize(path) / 1e9:.2f} GB")
