#!/usr/bin/env python3
"""02b_visium_scvi_integration.py

scVI integration of all Visium sections (GSE238264 + GSE281759).

Pipeline:
  1. Load all per-sample h5ads (raw counts), make var_names unique.
  2. Outer-merge on gene-symbol intersection; build sample/dataset/response obs.
  3. Merge cell-type annotations (harmonized) and CNV-derived tumor_status.
  4. Select 2,000 HVG with seurat_v3 flavor across batches.
  5. Train scVI with batch_key='sample_id'.
  6. kNN + UMAP on latent.
  7. 6-panel UMAP composite (sample, dataset, response, cell type, tumor, immunotype).
  8. Save integrated h5ad and figure.
"""
from __future__ import annotations

from pathlib import Path
import gc
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import scvi
import torch

PROJECT = Path(__file__).resolve().parents[2]
SD      = PROJECT / "data" / "spatial_data"
ANN_DIR = PROJECT / "out" / "02b_spatial_biomarkers" / "tables"
OUT_DIR = PROJECT / "out" / "02b_spatial_biomarkers"
H5_DIR  = PROJECT / "data" / "processed"
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"
H5_DIR.mkdir(parents=True, exist_ok=True)

GSE238 = SD / "GSE238264_HCC_IO_Visium"
GSE281 = SD / "GSE281759_HCC_MultiRegion"

GSE238_RESP = {
    "HCC1R": "Responder", "HCC2R": "Responder", "HCC3R": "Responder",
    "HCC4R": "Responder", "HCC5NR": "Non_Responder", "HCC6NR": "Non_Responder",
    "HCC7NR": "Non_Responder",
}

# ─────────────────────────────────────────────────────────────────────────────
# Cell-type harmonization
# ─────────────────────────────────────────────────────────────────────────────
HARMONIZE = {
    "B cell": "B", "B_cell": "B",
    "CD4+ T cell": "CD4 T", "CD4_T": "CD4 T",
    "CD8+ T cell": "CD8 T", "CD8_T": "CD8 T",
    "NK cell": "NK", "NK": "NK",
    "Treg": "Treg",
    "Macrophage": "Macrophage",
    "DC": "DC",
    "Endothelial": "Endothelial",
    "Fibroblast": "Fibroblast",
    "Hepatocyte": "Hepatocyte",
    "Tumor": "Malignant",
    "Plasma": "B",
}

def harmonize(x):
    return HARMONIZE.get(str(x), "Other")

# ─────────────────────────────────────────────────────────────────────────────
# Load per-sample h5ads
# ─────────────────────────────────────────────────────────────────────────────
print("Loading per-sample h5ads …")
adatas = []

# GSE238264
for sample, resp in GSE238_RESP.items():
    p = GSE238 / f"{sample}.h5ad"
    if not p.exists(): continue
    a = sc.read_h5ad(p)
    a.var_names_make_unique()
    a.obs["sample_id"] = sample
    a.obs["dataset"]   = "GSE238264"
    a.obs["response"]  = resp
    a.obs_names = [f"{sample}__{b}" for b in a.obs_names]
    adatas.append(a)
    print(f"  GSE238264/{sample}: {a.shape}")

# GSE281759
for p in sorted(GSE281.glob("*.h5ad")):
    sample = p.stem
    a = sc.read_h5ad(p)
    a.var_names_make_unique()
    a.obs["sample_id"] = sample
    a.obs["dataset"]   = "GSE281759"
    a.obs["response"]  = "NA"
    a.obs_names = [f"{sample}__{b}" for b in a.obs_names]
    adatas.append(a)
    print(f"  GSE281759/{sample}: {a.shape}")

# ─────────────────────────────────────────────────────────────────────────────
# Inner-join on gene symbols
# ─────────────────────────────────────────────────────────────────────────────
gene_sets = [set(a.var_names) for a in adatas]
common = sorted(set.intersection(*gene_sets))
print(f"\nCommon genes across {len(adatas)} samples: {len(common):,}")

adatas = [a[:, common].copy() for a in adatas]
adata = ad.concat(adatas, axis=0, join="outer", label=None, merge="unique")
del adatas
gc.collect()
print(f"Combined: {adata.shape}")

adata.layers["counts"] = adata.X.copy()

# ─────────────────────────────────────────────────────────────────────────────
# Merge cell-type annotations
# ─────────────────────────────────────────────────────────────────────────────
ann1 = pd.read_csv(ANN_DIR / "gse238264_celltype_annotations.csv", index_col=0)
ann1.index = [f"{s}__{b.split('__')[-1]}" if '__' not in b else b
              for s, b in zip(ann1["sample"], ann1.index)]
# original idx in ann1 is "{sample}_{barcode-suffix}" – matches old obs_names.
# Build a fresh map from sample+barcode using the celltype column:
ann1_map = {}
for idx, row in ann1.iterrows():
    sid = row["sample"]
    # original idx like "HCC1R_AAACAAGTATCTCCCA-1-0"
    barcode = idx[len(sid) + 1:] if idx.startswith(sid + "_") else idx
    ann1_map[f"{sid}__{barcode}"] = harmonize(row["celltype"])

ann2 = pd.read_csv(ANN_DIR / "gse281759_celltype_annotations.csv")
ann2_map = {f"{r['sample_id']}__{r['barcodes']}": harmonize(r["dominant_cell_type"])
            for _, r in ann2.iterrows()}
imm2_map = {f"{r['sample_id']}__{r['barcodes']}": str(r["spatial_immunotype"])
            for _, r in ann2.iterrows()}

adata.obs["cell_type"]   = [ann1_map.get(b, ann2_map.get(b, "Other")) for b in adata.obs_names]
adata.obs["immunotype"]  = [imm2_map.get(b, "NA") for b in adata.obs_names]

# CNV-derived tumor status
cnv_meta_path = TAB_DIR / "spatial_cnv_universal_spotmeta.csv"
if cnv_meta_path.exists():
    cnv_m = pd.read_csv(cnv_meta_path)
    cnv_m["key"] = cnv_m["sample_id"] + "__" + cnv_m["spot_id"]
    tumor_map = dict(zip(cnv_m["key"], cnv_m["tumor_status"]))
    score_map = dict(zip(cnv_m["key"], cnv_m["cnv_score"]))
else:
    tumor_map, score_map = {}, {}
adata.obs["tumor_status"] = [tumor_map.get(b, "unknown") for b in adata.obs_names]
adata.obs["cnv_score"]    = [float(score_map.get(b, np.nan)) for b in adata.obs_names]

print("Annotation summary:")
print("  cell_type:", adata.obs["cell_type"].value_counts().head(10).to_dict())
print("  tumor_status:", adata.obs["tumor_status"].value_counts().to_dict())

# ─────────────────────────────────────────────────────────────────────────────
# QC filter and HVG
# ─────────────────────────────────────────────────────────────────────────────
sc.pp.filter_cells(adata, min_counts=200)
sc.pp.filter_genes(adata, min_cells=10)
print(f"After QC: {adata.shape}")

sc.pp.highly_variable_genes(
    adata, layer="counts", flavor="seurat_v3",
    n_top_genes=2000, batch_key="sample_id", subset=False,
)
print(f"HVG: {int(adata.var['highly_variable'].sum())}")

# ─────────────────────────────────────────────────────────────────────────────
# scVI
# ─────────────────────────────────────────────────────────────────────────────
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
adata_hvg.X = adata_hvg.layers["counts"].copy()
scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts", batch_key="sample_id")

print("\nTraining scVI on CPU …")
torch.set_num_threads(8)
model = scvi.model.SCVI(adata_hvg, n_layers=2, n_latent=30, gene_likelihood="zinb")
model.train(
    max_epochs=80,
    accelerator="cpu",
    early_stopping=True,
    check_val_every_n_epoch=5,
    plan_kwargs={"lr": 1e-3},
    enable_progress_bar=True,
)

adata.obsm["X_scVI"] = model.get_latent_representation()
print(f"X_scVI shape: {adata.obsm['X_scVI'].shape}")

# ─────────────────────────────────────────────────────────────────────────────
# kNN + UMAP
# ─────────────────────────────────────────────────────────────────────────────
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=20)
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, resolution=0.6, key_added="leiden_scvi")

# ─────────────────────────────────────────────────────────────────────────────
# Save integrated h5ad
# ─────────────────────────────────────────────────────────────────────────────
out_h5 = H5_DIR / "spatial_integrated_scvi.h5ad"
# strip layers to keep file small
adata_save = adata.copy()
del adata_save.layers["counts"]
adata_save.write_h5ad(out_h5, compression="gzip")
print(f"  → {out_h5.name} ({adata.shape})")

# ─────────────────────────────────────────────────────────────────────────────
# 6-panel composite UMAP
# ─────────────────────────────────────────────────────────────────────────────
import matplotlib.colors as mcolors

def palette_n(n, cmap="tab20"):
    cm = plt.get_cmap(cmap, n)
    return [mcolors.to_hex(cm(i)) for i in range(n)]

fig = plt.figure(figsize=(18, 12), dpi=180)
gs = GridSpec(2, 3, figure=fig, wspace=0.18, hspace=0.22)

X = adata.obsm["X_umap"]

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
            leg = ax.legend(
                bbox_to_anchor=(1.02, 1), loc="upper left",
                fontsize=legend_size, frameon=False,
                markerscale=3, handletextpad=0.4, labelspacing=0.3,
                ncol=1 if len(unique) <= 14 else 2,
            )
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1", fontsize=8); ax.set_ylabel("UMAP2", fontsize=8)
    ax.set_title(title, fontsize=11)

# 1. Sample
ax = fig.add_subplot(gs[0, 0])
panel(ax, "sample_id", f"Sample (n={adata.obs['sample_id'].nunique()})", legend_size=5)

# 2. Dataset
ax = fig.add_subplot(gs[0, 1])
panel(ax, "dataset", "Dataset", palette={"GSE238264": "#1ABC9C", "GSE281759": "#E67E22"})

# 3. Response
ax = fig.add_subplot(gs[0, 2])
panel(ax, "response", "Response",
      palette={"Responder": "#2E86C1", "Non_Responder": "#C0392B", "NA": "#bbbbbb"})

# 4. Cell type
ax = fig.add_subplot(gs[1, 0])
ct_pal = {
    "CD8 T":"#1F77B4","CD4 T":"#17BECF","Treg":"#9467BD","NK":"#2CA02C",
    "B":"#8C564B","Macrophage":"#D62728","DC":"#FF7F0E","Endothelial":"#7F7F7F",
    "Fibroblast":"#BCBD22","Hepatocyte":"#E5C07B","Malignant":"#7B241C",
    "Other":"#CCCCCC",
}
panel(ax, "cell_type", "Cell type (harmonized)", palette=ct_pal, legend_size=7)

# 5. Tumor status
ax = fig.add_subplot(gs[1, 1])
panel(ax, "tumor_status", "Tumor status (CNV-derived)",
      palette={"tumor": "#7B241C", "non_tumor": "#85929E", "unknown": "#dddddd"})

# 6. Immunotype
ax = fig.add_subplot(gs[1, 2])
panel(ax, "immunotype", "Spatial immunotype")

fig.suptitle(
    f"scVI-integrated Visium atlas — {adata.shape[0]:,} spots × "
    f"{adata.obs['sample_id'].nunique()} sections × {len(common):,} genes",
    fontsize=13, y=0.98,
)

png = FIG_DIR / "spatial_scvi_umap_composite.png"
pdf = FIG_DIR / "spatial_scvi_umap_composite.pdf"
fig.savefig(png, bbox_inches="tight", dpi=200)
fig.savefig(pdf, bbox_inches="tight")
plt.close(fig)
print(f"  → {png.name}")
print(f"  → {pdf.name}")

# ─────────────────────────────────────────────────────────────────────────────
# Lightweight integration metric: mean kNN purity by sample (lower = more mixed)
# ─────────────────────────────────────────────────────────────────────────────
from sklearn.neighbors import NearestNeighbors

n_eval = min(15000, adata.shape[0])
idx = np.random.default_rng(0).choice(adata.shape[0], n_eval, replace=False)
Z = adata.obsm["X_scVI"][idx]
samples = adata.obs["sample_id"].values[idx]

nn = NearestNeighbors(n_neighbors=20).fit(Z)
_, knn_idx = nn.kneighbors(Z)
purity = np.mean([np.mean(samples[knn_idx[i]] == samples[i]) for i in range(n_eval)])
print(f"Mean within-sample kNN purity (lower=better mix): {purity:.3f}")

with open(TAB_DIR / "spatial_scvi_integration_metrics.txt", "w") as fh:
    fh.write(f"n_spots\t{adata.shape[0]}\n")
    fh.write(f"n_samples\t{adata.obs['sample_id'].nunique()}\n")
    fh.write(f"n_genes\t{adata.shape[1]}\n")
    fh.write(f"latent_dim\t{adata.obsm['X_scVI'].shape[1]}\n")
    fh.write(f"mean_within_sample_knn_purity\t{purity:.4f}\n")
print("Done.")
