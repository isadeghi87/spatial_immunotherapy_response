#!/usr/bin/env python3
"""02b_visium_scvi_patch_tumor_and_replot.py

Post-hoc fix: load the trained scVI-integrated h5ad, merge the correct
CNV-derived tumor_status (from out/02b_spatial_biomarkers/cnv/tables/...),
re-save and regenerate the composite UMAP figure.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.colors as mcolors

PROJECT = Path(__file__).resolve().parents[2]
H5      = PROJECT / "data" / "processed" / "spatial_integrated_scvi.h5ad"
CNV_TAB = PROJECT / "out" / "02b_spatial_biomarkers" / "cnv" / "tables" / "spatial_cnv_universal_spotmeta.csv"
FIG_DIR = PROJECT / "out" / "02b_spatial_biomarkers" / "figures"
TAB_DIR = PROJECT / "out" / "02b_spatial_biomarkers" / "tables"

print("Loading integrated h5ad …")
adata = sc.read_h5ad(H5)
print(adata.shape)

cnv = pd.read_csv(CNV_TAB)
cnv["key"] = cnv["sample_id"] + "__" + cnv["spot_id"]
tumor_map = dict(zip(cnv["key"], cnv["tumor_status"]))
score_map = dict(zip(cnv["key"], cnv["cnv_score"]))

adata.obs["tumor_status"] = [tumor_map.get(b, "unknown") for b in adata.obs_names]
adata.obs["cnv_score"]    = [float(score_map.get(b, np.nan)) for b in adata.obs_names]
print("tumor_status:", adata.obs["tumor_status"].value_counts().to_dict())

adata.write_h5ad(H5, compression="gzip")
print(f"  → {H5.name} re-saved")

# ─────────────────────────────────────────────────────────────────────────────
# Composite figure
# ─────────────────────────────────────────────────────────────────────────────
def palette_n(n, cmap="tab20"):
    cm = plt.get_cmap(cmap, n)
    return [mcolors.to_hex(cm(i)) for i in range(n)]

X = adata.obsm["X_umap"]

fig = plt.figure(figsize=(18, 12), dpi=180)
gs = GridSpec(2, 3, figure=fig, wspace=0.18, hspace=0.22)

def panel(ax, color_by, title, palette=None, legend_loc="right",
          continuous=False, cmap="viridis", legend_size=6):
    vals = adata.obs[color_by].values
    if continuous:
        sc_ = ax.scatter(X[:, 0], X[:, 1], c=vals, cmap=cmap,
                         s=1.5, alpha=0.7, rasterized=True, linewidths=0)
        cax = fig.colorbar(sc_, ax=ax, fraction=0.04, pad=0.02)
        cax.ax.tick_params(labelsize=6)
    else:
        cats = pd.Series(vals).astype(str)
        unique = sorted(cats.unique())
        if palette is None:
            palette = dict(zip(unique, palette_n(len(unique))))
        for cat in unique:
            m = (cats == cat).values
            ax.scatter(X[m, 0], X[m, 1], c=palette.get(cat, "#cccccc"),
                       s=1.5, alpha=0.7, rasterized=True,
                       linewidths=0, label=cat)
        if legend_loc:
            ax.legend(
                bbox_to_anchor=(1.02, 1), loc="upper left",
                fontsize=legend_size, frameon=False,
                markerscale=3, handletextpad=0.4, labelspacing=0.3,
                ncol=1 if len(unique) <= 14 else 2,
            )
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1", fontsize=8); ax.set_ylabel("UMAP2", fontsize=8)
    ax.set_title(title, fontsize=11)

ax = fig.add_subplot(gs[0, 0])
panel(ax, "sample_id", f"Sample (n={adata.obs['sample_id'].nunique()})", legend_size=5)

ax = fig.add_subplot(gs[0, 1])
panel(ax, "dataset", "Dataset", palette={"GSE238264": "#1ABC9C", "GSE281759": "#E67E22"})

ax = fig.add_subplot(gs[0, 2])
panel(ax, "response", "Response",
      palette={"Responder": "#2E86C1", "Non_Responder": "#C0392B", "NA": "#bbbbbb"})

ct_pal = {
    "CD8 T":"#1F77B4","CD4 T":"#17BECF","Treg":"#9467BD","NK":"#2CA02C",
    "B":"#8C564B","Macrophage":"#D62728","DC":"#FF7F0E","Endothelial":"#7F7F7F",
    "Fibroblast":"#BCBD22","Hepatocyte":"#E5C07B","Malignant":"#7B241C",
    "Other":"#CCCCCC",
}
ax = fig.add_subplot(gs[1, 0])
panel(ax, "cell_type", "Cell type (harmonized)", palette=ct_pal, legend_size=7)

ax = fig.add_subplot(gs[1, 1])
panel(ax, "tumor_status", "Tumor status (CNV-derived)",
      palette={"tumor": "#7B241C", "non_tumor": "#85929E", "unknown": "#dddddd"})

ax = fig.add_subplot(gs[1, 2])
panel(ax, "immunotype", "Spatial immunotype")

n_genes = adata.shape[1]
fig.suptitle(
    f"scVI-integrated Visium atlas — {adata.shape[0]:,} spots × "
    f"{adata.obs['sample_id'].nunique()} sections × {n_genes:,} HVG-pool genes",
    fontsize=13, y=0.98,
)

png = FIG_DIR / "spatial_scvi_umap_composite.png"
pdf = FIG_DIR / "spatial_scvi_umap_composite.pdf"
fig.savefig(png, bbox_inches="tight", dpi=200)
fig.savefig(pdf, bbox_inches="tight")
plt.close(fig)
print(f"  → {png.name}")
print(f"  → {pdf.name}")
print("Done.")
