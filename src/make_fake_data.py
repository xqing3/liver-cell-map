"""Generate a synthetic Liver Cell Map dataset that follows the team data contract.

Outputs (in --out):
  cells.parquet   id, x, y, cluster, species, condition, donor, meld_score
  graph.npz       scipy sparse CSR, n x n, symmetric 0/1 kNN connectivity (k=15)
  signaling.csv   sender, receiver, ligand, receptor, score, condition, species

Planted answer, so search and plots can be checked:
  - "Disease macrophage" cells come almost only from disease donors, in both species.
  - meld_score is higher in disease, highest in disease macrophages, and carries a
    small per-donor offset (to mimic the donor-vs-disease risk).
  - Disease macrophage -> Stellate signaling (SPP1-CD44, PDGFB-PDGFRB, TGFB1-TGFBR1)
    is boosted in disease.

Usage:
  python3 src/make_fake_data.py --n-cells 5000 --out data/fake/n5000
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors

CLUSTERS = [
    "Hepatocyte", "Cholangiocyte", "Endothelial", "Stellate",
    "Kupffer cell", "Disease macrophage", "T cell", "B cell",
]
DISEASE_MAC = CLUSTERS.index("Disease macrophage")

# Ring order puts related cell types next to each other. Cluster proportions per condition; disease macrophages appear almost only in disease.
PROPS = {
    "healthy": np.array([0.40, 0.05, 0.15, 0.08, 0.15, 0.01, 0.10, 0.06]),
    "disease": np.array([0.30, 0.05, 0.13, 0.12, 0.10, 0.15, 0.10, 0.05]),
}
DONORS_PER_GROUP = {"mouse": 3, "human": 5}  # matches the real study design

LR_PAIRS = [
    ("SPP1", "CD44"), ("TGFB1", "TGFBR1"), ("PDGFB", "PDGFRB"), ("CCL2", "CCR2"),
    ("TNF", "TNFRSF1A"), ("IL1B", "IL1R1"), ("CXCL12", "CXCR4"), ("VEGFA", "KDR"),
    ("APOE", "TREM2"), ("CSF1", "CSF1R"),
]
BOOSTED = {("SPP1", "CD44"), ("PDGFB", "PDGFRB"), ("TGFB1", "TGFBR1")}


def make_cells(n, rng, latent_dim=10):
    species = rng.choice(["mouse", "human"], size=n)
    condition = rng.choice(["healthy", "disease"], size=n)

    donor = np.empty(n, dtype=object)
    for s, k in DONORS_PER_GROUP.items():
        for c in ("healthy", "disease"):
            m = (species == s) & (condition == c)
            donor[m] = [f"{s}_{c}_{i + 1}" for i in rng.integers(0, k, m.sum())]

    cluster = np.empty(n, dtype=int)
    for c, p in PROPS.items():
        m = condition == c
        cluster[m] = rng.choice(len(CLUSTERS), size=m.sum(), p=p / p.sum())

    # 2D layout centers on a ring, plus extra latent dims with their own centers.
    # kNN is built in the full latent space, so the graph is close to, not identical to, the layout.
    angles = np.linspace(0, 2 * np.pi, len(CLUSTERS), endpoint=False)
    centers_2d = np.c_[np.cos(angles), np.sin(angles)] * 6
    centers_rest = rng.normal(0, 0.8, size=(len(CLUSTERS), latent_dim - 2))
    # Disease macrophages sit between Kupffer cells and their own spot, so the boundary is fuzzy.
    kupffer = CLUSTERS.index("Kupffer cell")
    centers_2d[DISEASE_MAC] = 0.55 * centers_2d[kupffer] + 0.45 * centers_2d[DISEASE_MAC]
    centers_rest[DISEASE_MAC] = 0.6 * centers_rest[kupffer] + 0.4 * centers_rest[DISEASE_MAC]
    noise = rng.normal(0, 1.3, size=(n, latent_dim))
    latent = np.c_[centers_2d[cluster], centers_rest[cluster]] + noise
    latent[:, 0] += np.where(species == "human", 0.6, -0.6)  # mild species shift

    # Drawn layout: same ring, spread wider with less noise, so cell types read as islands like a UMAP.
    layout = 1.7 * centers_2d[cluster] + 0.55 * noise[:, :2]
    layout[:, 0] += np.where(species == "human", 0.4, -0.4)

    donor_offset = {d: rng.normal(0, 0.05) for d in np.unique(donor)}
    score = (
        0.2
        + 0.2 * (condition == "disease")
        + 0.4 * ((cluster == DISEASE_MAC) & (condition == "disease"))
        + np.array([donor_offset[d] for d in donor])
        + rng.normal(0, 0.08, n)
    )

    cells = pd.DataFrame({
        "id": np.arange(n),
        "x": layout[:, 0].astype(np.float32),
        "y": layout[:, 1].astype(np.float32),
        "cluster": [CLUSTERS[i] for i in cluster],
        "species": species,
        "condition": condition,
        "donor": donor.astype(str),
        "meld_score": np.clip(score, 0, 1).astype(np.float32),
    })
    return cells, latent


def make_graph(latent, k=15):
    nn = NearestNeighbors(n_neighbors=k + 1, n_jobs=-1).fit(latent)
    _, idx = nn.kneighbors(latent)
    n = latent.shape[0]
    rows = np.repeat(np.arange(n), k)
    cols = idx[:, 1:].ravel()  # drop self
    a = sp.csr_matrix((np.ones_like(rows, dtype=np.float32), (rows, cols)), shape=(n, n))
    return a.maximum(a.T).tocsr()


def make_signaling(rng):
    rows = []
    for species in ("mouse", "human"):
        for condition in ("healthy", "disease"):
            for s in CLUSTERS:
                for r in CLUSTERS:
                    for lig, rec in LR_PAIRS:
                        score = rng.lognormal(-1.5, 0.6)
                        if (condition == "disease" and s == "Disease macrophage"
                                and r == "Stellate" and (lig, rec) in BOOSTED):
                            score += 1.5
                        rows.append((s, r, lig, rec, round(score, 4), condition, species))
    return pd.DataFrame(rows, columns=["sender", "receiver", "ligand", "receptor",
                                       "score", "condition", "species"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-cells", type=int, default=5000)
    ap.add_argument("--k", type=int, default=15)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    out = args.out or Path("data/fake") / f"n{args.n_cells}"
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    cells, latent = make_cells(args.n_cells, rng)
    graph = make_graph(latent, args.k)
    signaling = make_signaling(rng)

    cells.to_parquet(out / "cells.parquet", index=False)
    sp.save_npz(out / "graph.npz", graph)
    signaling.to_csv(out / "signaling.csv", index=False)
    (out / "meta.json").write_text(json.dumps(vars(args) | {"out": str(out), "synthetic": True}, indent=2))

    print(f"wrote {out}: {len(cells)} cells, {graph.nnz // 2} undirected edges, {len(signaling)} signaling rows")


if __name__ == "__main__":
    main()
