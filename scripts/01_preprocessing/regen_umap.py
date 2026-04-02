#!/usr/bin/env python
"""Regenerate the 4-panel integration UMAP from the saved atlas."""
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

FIG = Path("out/01_preprocessing/figures")
FIG.mkdir(parents=True, exist_ok=True)

print("Loading atlas (backed)...")
adata = ad.read_h5ad("data/processed/integrated_io_atlas.h5ad", backed="r")
umap = np.array(adata.obsm["X_umap"])
obs = adata.obs.copy()
print(f"  {len(obs):,} cells, UMAP shape: {umap.shape}")

fig, axes = plt.subplots(2, 2, figsize=(14, 13))

# Panel A: Dataset
ax = axes[0, 0]
ds_list = sorted(obs["dataset"].dropna().unique())
cmap_ds = plt.cm.Set2(np.linspace(0, 1, len(ds_list)))
for i, ds in enumerate(ds_list):
    idx = np.where(obs["dataset"].values == ds)[0]
    if len(idx) > 0:
        ax.scatter(umap[idx, 0], umap[idx, 1],
                   s=0.5, alpha=0.2, c=[cmap_ds[i]], label=ds, rasterized=True)
        print(f"  A: {ds} -> {len(idx):,} cells")
ax.legend(fontsize=7, markerscale=10)
ax.set_title("A. Dataset", fontweight="bold")
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

# Panel B: Study
ax = axes[0, 1]
studies = sorted(obs["study"].dropna().unique())
cmap_st = plt.cm.tab10(np.linspace(0, 1, len(studies)))
for i, st in enumerate(studies):
    idx = np.where(obs["study"].values == st)[0]
    if len(idx) > 0:
        ax.scatter(umap[idx, 0], umap[idx, 1],
                   s=0.5, alpha=0.2, c=[cmap_st[i]], label=st, rasterized=True)
        print(f"  B: {st} -> {len(idx):,} cells")
ax.legend(fontsize=8, markerscale=10)
ax.set_title("B. Study", fontweight="bold")
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

# Panel C: Tissue
ax = axes[1, 0]
tissues = sorted(obs["tissue"].dropna().unique())
cmap_ti = plt.cm.Paired(np.linspace(0, 1, max(len(tissues), 2)))
for i, ti in enumerate(tissues):
    idx = np.where(obs["tissue"].values == ti)[0]
    if len(idx) > 0:
        ax.scatter(umap[idx, 0], umap[idx, 1],
                   s=0.5, alpha=0.2, c=[cmap_ti[i]], label=ti, rasterized=True)
        print(f"  C: {ti} -> {len(idx):,} cells")
ax.legend(fontsize=8, markerscale=10)
ax.set_title("C. Tissue", fontweight="bold")
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

# Panel D: IO Response
ax = axes[1, 1]
resp_colors = {"Responder": "#2ecc71", "Non-Responder": "#e74c3c"}
for resp, color in resp_colors.items():
    idx = np.where(obs["response"].values == resp)[0]
    if len(idx) > 0:
        ax.scatter(umap[idx, 0], umap[idx, 1],
                   s=0.5, alpha=0.2, c=color, label=f"{resp} ({len(idx):,})",
                   rasterized=True)
        print(f"  D: {resp} -> {len(idx):,} cells")
idx_na = np.where(obs["response"].isna().values)[0]
if len(idx_na) > 0:
    ax.scatter(umap[idx_na, 0], umap[idx_na, 1],
               s=0.3, alpha=0.05, c="gray", label=f"No label ({len(idx_na):,})",
               rasterized=True)
    print(f"  D: No label -> {len(idx_na):,} cells")
ax.legend(fontsize=8, markerscale=10)
ax.set_title("D. IO Response (GSE206325)", fontweight="bold")
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

plt.suptitle(f"Integrated HCC IO Atlas \u2014 {len(obs):,} cells\n"
             f"Harmony batch correction, {len(ds_list)} datasets",
             fontweight="bold", fontsize=13)
plt.tight_layout()
for ext in ["png", "pdf"]:
    plt.savefig(FIG / f"01_integration_umap.{ext}", dpi=140, bbox_inches="tight")
plt.close()
print(f"\nSaved to {FIG / '01_integration_umap.png'}")
