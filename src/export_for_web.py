"""Convert the data-contract files into one JSON file for the browser prototype.

Usage:
  python3 src/export_for_web.py data/fake/n5000 prototype/data.json
  python3 src/export_for_web.py data/CSE6242_data/export/full prototype/data_full.json --graph-bin

--graph-bin keeps large graphs out of the JSON: the graph goes to <dst stem>_graph.txt as the upper
triangle in CSR form (uint32 indptr of length N+1, then uint32 column indices, little-endian),
base64-encoded because artifact hosting serves text but not raw binary. The browser decodes it and
rebuilds the symmetric adjacency. Donors are stored as indices into a donor list.
"""
import base64
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def main(src, dst, graph_bin=False):
    src = Path(src)
    cells = pd.read_parquet(src / "cells.parquet").sort_values("id")
    graph = sp.load_npz(src / "graph.npz").tocsr()
    signaling = pd.read_csv(src / "signaling.csv")

    clusters = sorted(cells.cluster.unique())
    cl_idx = cells.cluster.map({c: i for i, c in enumerate(clusters)}).to_numpy()

    # Cluster summary for the zoomed-out view, per species and merged.
    summary = []
    for species in ("all", "mouse", "human"):
        sub = cells if species == "all" else cells[cells.species == species]
        for i, c in enumerate(clusters):
            g = sub[sub.cluster == c]
            if len(g):
                summary.append({"cluster": i, "species": species, "x": round(float(g.x.mean()), 3),
                                "y": round(float(g.y.mean()), 3), "n": int(len(g)),
                                "score": None if g.meld_score.isna().all() else round(float(g.meld_score.mean()), 3)})

    # PAGA-style links: number of graph edges between each pair of clusters.
    coo = sp.triu(graph, k=1).tocoo()
    a, b = cl_idx[coo.row], cl_idx[coo.col]
    m = a != b
    pairs = pd.Series(1, index=pd.MultiIndex.from_arrays([np.minimum(a[m], b[m]), np.maximum(a[m], b[m])])).groupby(level=[0, 1]).sum()
    links = [{"source": int(s), "target": int(t), "edges": int(n)} for (s, t), n in pairs.items()]

    meta_path = src / "meta.json"
    donors = sorted(cells.donor.unique())
    out = {
        "meta": json.loads(meta_path.read_text()) if meta_path.exists() else {},
        "clusters": clusters,
        "species": ["mouse", "human"],
        "conditions": ["healthy", "disease"],
        "cells": {
            # float64 first, so float32 noise (-1.8860000371932983) does not bloat the JSON
            "x": cells.x.astype("float64").round(3).tolist(),
            "y": cells.y.astype("float64").round(3).tolist(),
            "cluster": cl_idx.tolist(),
            "species": (cells.species == "human").astype(int).tolist(),
            "condition": (cells.condition == "disease").astype(int).tolist(),
            "donor": cells.donor.map({d: i for i, d in enumerate(donors)}).tolist() if graph_bin else cells.donor.tolist(),
            # NaN (MELD not run yet) becomes null, since JSON has no NaN
            "score": [None if np.isnan(v) else round(float(v), 3) for v in cells.meld_score],
        },
        "graph": {"indptr": graph.indptr.tolist(), "indices": graph.indices.tolist()},
        "clusterSummary": summary,
        "clusterLinks": links,
        "signaling": signaling.to_dict(orient="split")["data"],
        "signalingColumns": list(signaling.columns),
    }
    if graph_bin:
        out["donors"] = donors
        up = sp.triu(graph, k=1).tocsr()
        up.sort_indices()
        bin_path = Path(dst).with_name(Path(dst).stem + "_graph.txt")
        bin_path.parent.mkdir(parents=True, exist_ok=True)
        raw = up.indptr.astype("<u4").tobytes() + up.indices.astype("<u4").tobytes()
        bin_path.write_bytes(base64.b64encode(raw))
        out["graph"] = {"format": "upper-csr-u32-base64", "file": bin_path.name, "n": int(graph.shape[0]), "edges": int(up.nnz)}
        print(f"wrote {bin_path}: {bin_path.stat().st_size / 1e6:.2f} MB, {up.nnz} undirected edges")
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    Path(dst).write_text(json.dumps(out, separators=(",", ":")))
    print(f"wrote {dst}: {Path(dst).stat().st_size / 1e6:.2f} MB, {len(cells)} cells, {len(links)} cluster links")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    main(args[0], args[1], graph_bin="--graph-bin" in sys.argv)
