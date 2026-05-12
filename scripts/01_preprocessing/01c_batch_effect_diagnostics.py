#!/usr/bin/env python3
"""
01c_batch_effect_diagnostics.py
===============================

Quantitative assessment of batch correction on the integrated atlas,
using the primary integrated embedding saved in the atlas and optionally
training scVI on the subsample if no integrated latent is available.

Metrics reported
----------------
- Per-dataset silhouette (PC before vs integrated embedding) — batch mixing
- Per-cell-type silhouette (biological conservation)
- Average-silhouette-width (ASW) batch and bio, and their balance
- kBET-style batch acceptance rate (proxy implementation over kNN graph):
  for each cell, test whether the k nearest neighbours form a chi-square
  null composition matching the global batch distribution.
- LISI (local inverse Simpson's index) over local neighbourhoods
  (iLISI for integration, cLISI for cell-type purity).
- Optional scVI latent embedding on a 60k-cell subsample if the atlas
    does not already contain an scVI latent.

Outputs
-------
out/01_preprocessing/tables/
    batch_effect_metrics.csv             -- headline metrics table
    scvi_vs_harmony_metrics.csv          -- if scVI run
    per_dataset_silhouette.csv
out/01_preprocessing/figures/
    01c_batch_metrics_summary.png
    01c_umap_integration_comparison.png
    01c_lisi_histograms.png
"""
from __future__ import annotations
import warnings
import sys
from pathlib import Path

warnings.filterwarnings("ignore")

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy.stats import chisquare
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

from config.publication_style import apply_style, PAL
apply_style()

OUT = PROJECT / "out" / "01_preprocessing"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

ATLAS = PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad"

SAMPLE_N = 60_000
SEED = 0
N_NN = 30

print("=" * 72)
print("[01c] Batch-effect diagnostics (scVI-first)")
print("=" * 72)

# ----------------------------------------------------------------------
# 1. Load atlas and subsample
# ----------------------------------------------------------------------
print("\n[1] Loading atlas …")
adata = sc.read_h5ad(ATLAS)
print(f"   shape: {adata.shape}")
rng = np.random.default_rng(SEED)
if adata.n_obs > SAMPLE_N:
    idx = rng.choice(adata.n_obs, SAMPLE_N, replace=False)
    adata = adata[idx].copy()
    print(f"   subsampled to: {adata.n_obs}")

batch_col = next((c for c in ("dataset", "batch", "sample", "donor")
                  if c in adata.obs.columns), None)
ct_col = next((c for c in ("cell_type", "celltype", "cell_type_original")
               if c in adata.obs.columns), None)
if batch_col is None:
    raise RuntimeError("No batch column found in adata.obs")
if ct_col is None:
    raise RuntimeError("No cell-type column found in adata.obs")

print(f"   batch column: {batch_col!r} ({adata.obs[batch_col].nunique()} batches)")
print(f"   cell-type column: {ct_col!r} "
      f"({adata.obs[ct_col].nunique()} types)")

# ----------------------------------------------------------------------
# 2. Build PC (pre-integration) and fetch primary integrated embedding
# ----------------------------------------------------------------------
print("\n[2] Computing PCA pre-integration …")
if "X_pca" not in adata.obsm:
    sc.pp.pca(adata, n_comps=50)
pca = adata.obsm["X_pca"]
print(f"   PCA shape: {pca.shape}")

hkey = None
for cand in ("X_scVI", "X_integrated", "X_pca_harmony", "X_harmony", "X_umap_harmony"):
    if cand in adata.obsm:
        hkey = cand
        break
if hkey is None:
    print("   [info] Integrated embedding not found in .obsm; training scVI on the subsample.")
    try:
        import scvi  # type: ignore
        scvi_ad = adata.copy()
        sc.pp.highly_variable_genes(scvi_ad, n_top_genes=3000,
                                    flavor="seurat_v3", batch_key=batch_col)
        scvi_ad = scvi_ad[:, scvi_ad.var["highly_variable"]].copy()
        scvi.model.SCVI.setup_anndata(
            scvi_ad,
            layer="counts" if "counts" in scvi_ad.layers else None,
            batch_key=batch_col,
        )
        model = scvi.model.SCVI(scvi_ad, n_layers=2, n_latent=30)
        model.train(max_epochs=30, batch_size=512, early_stopping=True,
                    plan_kwargs=dict(lr=1e-3))
        adata.obsm["X_scVI"] = model.get_latent_representation(scvi_ad)
        hkey = "X_scVI"
    except Exception as e:
        print(f"   [warn] scVI fallback failed: {e}; using PCA as integrated proxy")
        adata.obsm["X_integrated"] = pca.copy()
        hkey = "X_integrated"
integrated = adata.obsm[hkey]
if hkey == "X_scVI":
    embedding_label = "scVI"
elif hkey in {"X_pca_harmony", "X_harmony", "X_umap_harmony"}:
    embedding_label = "Harmony"
else:
    embedding_label = str(adata.uns.get("integration_method", "Integrated"))
print(f"   Integrated embedding ({embedding_label}): {integrated.shape}")


# ----------------------------------------------------------------------
# 3. Silhouette (ASW) — batch vs bio on PCA vs Harmony
# ----------------------------------------------------------------------
def _safe_sil(X, labels, sample=min(20_000, SAMPLE_N), seed=0):
    labels = np.asarray(labels)
    keep = pd.Series(labels).value_counts()
    keep = keep[keep >= 10].index
    m = np.isin(labels, keep)
    X = X[m]
    labels = labels[m]
    if len(labels) > sample:
        rs = np.random.default_rng(seed)
        sel = rs.choice(len(labels), sample, replace=False)
        X = X[sel]
        labels = labels[sel]
    if len(np.unique(labels)) < 2:
        return np.nan
    return float(silhouette_score(X, labels, metric="euclidean"))


print("\n[3] ASW silhouette: batch (lower=better mixing) / bio (higher=better) …")
sil_rows = []
for emb_name, emb in [("PCA", pca), (embedding_label, integrated)]:
    sb = _safe_sil(emb, adata.obs[batch_col].values)
    sc_ = _safe_sil(emb, adata.obs[ct_col].astype(str).values)
    asw_batch = (sb + 1) / 2  # scale to [0,1]
    asw_bio = (sc_ + 1) / 2
    mix = 1.0 - asw_batch
    bal = (mix + asw_bio) / 2
    sil_rows.append(dict(embedding=emb_name,
                         asw_batch=asw_batch,
                         asw_bio=asw_bio,
                         batch_mixing=mix,
                         asw_balance=bal))
    print(f"   {emb_name:10s}  asw_batch={asw_batch:.3f}  "
          f"asw_bio={asw_bio:.3f}  mixing={mix:.3f}  balance={bal:.3f}")
sil_df = pd.DataFrame(sil_rows)


# ----------------------------------------------------------------------
# 4. Per-dataset silhouette
# ----------------------------------------------------------------------
print("\n[4] Per-dataset silhouette …")
pd_rows = []
for emb_name, emb in [("PCA", pca), (embedding_label, integrated)]:
    for dset in adata.obs[batch_col].unique():
        mask = adata.obs[batch_col] == dset
        if mask.sum() < 50:
            continue
        # silhouette of this dataset vs "others"
        lbl = np.where(mask, "this", "other")
        if len(np.unique(lbl)) < 2:
            continue
        s = _safe_sil(emb, lbl, sample=10_000)
        pd_rows.append(dict(embedding=emb_name, dataset=dset,
                            n_cells=int(mask.sum()),
                            sil_vs_rest=(s + 1) / 2 if s is not None else np.nan))
per_dataset = pd.DataFrame(pd_rows)
per_dataset.to_csv(TAB / "per_dataset_silhouette.csv", index=False)


# ----------------------------------------------------------------------
# 5. kBET-proxy and LISI on kNN graph
# ----------------------------------------------------------------------
def _knn_indices(X, k=N_NN):
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="auto").fit(X)
    return nn.kneighbors(X, return_distance=False)[:, 1:]


def _kbet(emb, batch, k=N_NN, sample=5000, seed=0):
    rng = np.random.default_rng(seed)
    n = emb.shape[0]
    sel = rng.choice(n, min(sample, n), replace=False)
    knn = _knn_indices(emb[sel], k=k)
    cats, cat_idx = np.unique(batch, return_inverse=True)
    global_prop = np.bincount(cat_idx, minlength=len(cats)) / len(cat_idx)
    reject = 0
    for i, neigh in enumerate(knn):
        obs_counts = np.bincount(cat_idx[sel[neigh]], minlength=len(cats))
        exp_counts = global_prop * k
        mask = exp_counts > 1
        if mask.sum() < 2:
            continue
        try:
            _, p = chisquare(obs_counts[mask], exp_counts[mask])
        except Exception:
            continue
        if p < 0.05:
            reject += 1
    return 1.0 - reject / len(knn)


def _lisi(emb, labels, k=N_NN, sample=5000, seed=0):
    rng = np.random.default_rng(seed)
    n = emb.shape[0]
    sel = rng.choice(n, min(sample, n), replace=False)
    knn = _knn_indices(emb[sel], k=k)
    lab = np.asarray(labels)[sel]
    lab_u, lab_i = np.unique(lab, return_inverse=True)
    neigh_lab = lab_i[knn]
    lisi_vals = np.empty(len(sel))
    for i in range(len(sel)):
        counts = np.bincount(neigh_lab[i], minlength=len(lab_u))
        p = counts / counts.sum()
        p = p[p > 0]
        lisi_vals[i] = 1.0 / np.sum(p ** 2)  # effective diversity
    return lisi_vals


print("\n[5] kBET-proxy + LISI …")
kbet_rows = []
lisi_rows = {"PCA": {}, embedding_label: {}}
for emb_name, emb in [("PCA", pca), (embedding_label, integrated)]:
    kb = _kbet(emb, adata.obs[batch_col].values, k=N_NN)
    ilisi = _lisi(emb, adata.obs[batch_col].values, k=N_NN)
    clisi = _lisi(emb, adata.obs[ct_col].astype(str).values, k=N_NN)
    kbet_rows.append(dict(embedding=emb_name,
                          kBET_accept=kb,
                          iLISI_mean=float(np.mean(ilisi)),
                          iLISI_median=float(np.median(ilisi)),
                          cLISI_mean=float(np.mean(clisi)),
                          cLISI_median=float(np.median(clisi))))
    lisi_rows[emb_name] = {"iLISI": ilisi, "cLISI": clisi}
    print(f"   {emb_name:10s} kBET_accept={kb:.3f} "
          f"iLISI={np.median(ilisi):.2f} cLISI={np.median(clisi):.2f}")
kbet_df = pd.DataFrame(kbet_rows)


# ----------------------------------------------------------------------
# 6. Optional scVI rerun for orthogonal check
# ----------------------------------------------------------------------
scvi_metrics = None
if embedding_label != "scVI":
    try:
        import scvi  # type: ignore
        print("\n[6] Running scVI on the subsample (may be slow, no GPU) …")
        scvi_ad = adata.copy()
        sc.pp.highly_variable_genes(scvi_ad, n_top_genes=3000,
                                    flavor="seurat_v3", batch_key=batch_col)
        scvi_ad = scvi_ad[:, scvi_ad.var["highly_variable"]].copy()
        scvi.model.SCVI.setup_anndata(
            scvi_ad,
            layer="counts" if "counts" in scvi_ad.layers else None,
            batch_key=batch_col,
        )
        model = scvi.model.SCVI(scvi_ad, n_layers=2, n_latent=30)
        model.train(max_epochs=30, batch_size=512, early_stopping=True,
                    plan_kwargs=dict(lr=1e-3))
        scvi_ad.obsm["X_scVI"] = model.get_latent_representation()
        sv = scvi_ad.obsm["X_scVI"]
        sb = _safe_sil(sv, scvi_ad.obs[batch_col].values)
        sc_ = _safe_sil(sv, scvi_ad.obs[ct_col].astype(str).values)
        kb = _kbet(sv, scvi_ad.obs[batch_col].values, k=N_NN)
        ilisi = _lisi(sv, scvi_ad.obs[batch_col].values, k=N_NN)
        clisi = _lisi(sv, scvi_ad.obs[ct_col].astype(str).values, k=N_NN)
        scvi_metrics = pd.DataFrame([dict(
            embedding="scVI",
            asw_batch=(sb + 1) / 2,
            asw_bio=(sc_ + 1) / 2,
            batch_mixing=1.0 - (sb + 1) / 2,
            asw_balance=((1.0 - (sb + 1) / 2) + (sc_ + 1) / 2) / 2,
            kBET_accept=kb,
            iLISI_median=float(np.median(ilisi)),
            cLISI_median=float(np.median(clisi)),
        )])
        scvi_metrics.to_csv(TAB / "scvi_vs_harmony_metrics.csv", index=False)
    except Exception as e:
        print(f"\n[6] scVI benchmark skipped: {e}")
else:
    print("\n[6] Primary atlas embedding is already scVI; no secondary scVI rerun needed.")

# ----------------------------------------------------------------------
# 7. Combine metrics and save
# ----------------------------------------------------------------------
metrics = sil_df.merge(kbet_df, on="embedding", how="outer")
if scvi_metrics is not None:
    metrics = pd.concat([metrics, scvi_metrics], ignore_index=True)
metrics.to_csv(TAB / "batch_effect_metrics.csv", index=False)
print(f"\n   metrics saved to {TAB/'batch_effect_metrics.csv'}")
print(metrics.to_string(index=False, float_format=lambda x: f"{x:.3f}"))


# ----------------------------------------------------------------------
# 8. Figures
# ----------------------------------------------------------------------
def _save(fig, stem):
    fig.savefig(FIG / f"{stem}.png", dpi=300, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print(f"   saved {stem}.png")


print("\n[7] Figures …")

# metrics summary bar
mm = metrics.melt(id_vars=["embedding"],
                  value_vars=["batch_mixing", "asw_bio", "kBET_accept",
                              "iLISI_median", "cLISI_median"],
                  var_name="metric", value_name="value")
fig, ax = plt.subplots(figsize=(9, 4.2), constrained_layout=True)
pal = {"PCA": "#B0B0B0", "Harmony": "#117733", "scVI": "#332288"}
pal[embedding_label] = "#117733"
sns.barplot(data=mm, x="metric", y="value", hue="embedding",
            palette=pal, ax=ax, edgecolor="black", linewidth=0.5)
ax.set_title("Batch-effect diagnostics summary",
             fontsize=10.5, fontweight="bold")
ax.set_xlabel("")
ax.tick_params(axis="x", rotation=20)
ax.legend(title="Embedding", frameon=False, fontsize=8, loc="upper right")
_save(fig, "01c_batch_metrics_summary")

# LISI histograms
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
for ax, key, ttl in [(axes[0], "iLISI",
                      "iLISI (higher = better batch mixing)"),
                     (axes[1], "cLISI",
                      "cLISI (lower = better cell-type purity)")]:
    for emb_name, color in [("PCA", "#B0B0B0"), (embedding_label, "#117733")]:
        vals = lisi_rows[emb_name][key]
        ax.hist(vals, bins=40, alpha=0.55, color=color,
                density=True, label=emb_name, edgecolor="black",
                linewidth=0.4)
    ax.set_title(ttl, fontsize=10, fontweight="bold")
    ax.set_xlabel("LISI")
    ax.set_ylabel("density")
    ax.legend(frameon=False, fontsize=8)
_save(fig, "01c_lisi_histograms")

# UMAP comparison
try:
    print("   computing PCA-only UMAP for comparison …")
    ad_pca = adata.copy()
    sc.pp.neighbors(ad_pca, use_rep="X_pca", n_neighbors=15)
    sc.tl.umap(ad_pca)
    ad_pca.obsm["X_umap_pca"] = ad_pca.obsm["X_umap"]
    sc.pp.neighbors(adata, use_rep=hkey, n_neighbors=15)
    sc.tl.umap(adata)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5),
                              constrained_layout=True)
    for ax, emb_xy, emb_name in [
        (axes[0], ad_pca.obsm["X_umap_pca"], "PCA (pre-integration)"),
        (axes[1], adata.obsm["X_umap"], f"{embedding_label} (post-integration)")]:
        batches = adata.obs[batch_col].astype(str).values
        for b in sorted(set(batches)):
            m = batches == b
            ax.scatter(emb_xy[m, 0], emb_xy[m, 1],
                       s=1.2, alpha=0.3, label=b, edgecolor="none")
        ax.set_title(emb_name, fontsize=11, fontweight="bold")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_aspect("equal")
    axes[1].legend(markerscale=6, fontsize=7, frameon=False,
                   loc="upper left", bbox_to_anchor=(1.02, 1.0))
    _save(fig, "01c_umap_integration_comparison")
except Exception as e:
    print(f"   UMAP figure skipped: {e}")

print("\n[01c] DONE")
