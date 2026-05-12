#!/usr/bin/env python3
"""
01d_scvi_benchmark.py
=====================

Orthogonal integration benchmark: scVI vs Harmony vs PCA on a 60k-cell
subsample of the integrated atlas. Runs in the isolated venv at
`envs/scvi_env/` (scvi-tools 1.1.6, jax 0.4.30, torch CUDA).

Outputs
-------
out/01_preprocessing/tables/
    scvi_vs_harmony_metrics.csv
out/01_preprocessing/figures/
    01d_scvi_umap.png
    01d_scvi_metric_bars.png

Notes
-----
Metrics are computed on 10k-cell subsamples for speed and parity with
01c. Random seed fixed to 0 for reproducibility.
"""
from __future__ import annotations
import os
# Hide GPU from torch because installed CUDA driver (12020) is older than
# bundled torch CUDA runtime; scVI is run on CPU for this benchmark.
os.environ["CUDA_VISIBLE_DEVICES"] = ""
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
import scanpy as sc
from scipy.stats import chisquare
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

OUT = PROJECT / "out" / "01_preprocessing"
FIG = OUT / "figures"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

ATLAS = PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad"
SAMPLE_N = 60_000
METRIC_N = 10_000
N_NN = 30
SEED = 0

print("=" * 72)
print("[01d] scVI orthogonal integration benchmark")
print("=" * 72)

# ------------------------------------------------------------------
# 1. Load atlas via h5py (anndata version mismatch on /uns/log1p);
#    we only need counts, obs, var
# ------------------------------------------------------------------
print("\n[1] Loading atlas via h5py (subsample, counts-only) …")
import h5py
import scipy.sparse as sp
import anndata as ad

def _read_df(g):
    cols = {}
    for k in g.keys():
        if k.startswith("_"): continue
        item = g[k]
        if isinstance(item, h5py.Group) and "codes" in item:
            codes = item["codes"][()]
            cats = item["categories"][()]
            cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats])
            vals = pd.Categorical.from_codes(codes, categories=cats)
            cols[k] = vals
        else:
            arr = item[()]
            if arr.dtype.kind == "O" or arr.dtype.kind == "S":
                arr = np.array([x.decode() if isinstance(x, bytes) else x for x in arr])
            cols[k] = arr
    idx_name = g.attrs.get("_index", "_index")
    if isinstance(idx_name, bytes): idx_name = idx_name.decode()
    idx = cols.pop(idx_name, None)
    df = pd.DataFrame(cols)
    if idx is not None:
        df.index = np.array([x.decode() if isinstance(x, bytes) else x for x in idx])
    return df

rng = np.random.default_rng(SEED)
with h5py.File(ATLAS, "r") as f:
    n_obs = int(f["obs"].attrs.get("column-order", f["obs"].keys()) and
                f[f"obs/{list(f['obs'].keys())[0]}"].shape[0]
                if not hasattr(f['obs'], 'attrs') else 0)
    # robust n_obs
    first_key = [k for k in f['obs'].keys() if not k.startswith("_")][0]
    first_item = f[f'obs/{first_key}']
    if isinstance(first_item, h5py.Group) and 'codes' in first_item:
        n_obs = first_item['codes'].shape[0]
    else:
        n_obs = first_item.shape[0]
    print(f"   full n_obs: {n_obs}")
    sample_n = min(SAMPLE_N, n_obs)
    idx = np.sort(rng.choice(n_obs, sample_n, replace=False))
    # Counts layer as CSR
    c = f["layers/counts"]
    data = c["data"][()]
    indices = c["indices"][()]
    indptr = c["indptr"][()]
    n_var = f["var"].attrs.get("_index") and f[f"var/{[k for k in f['var'].keys() if not k.startswith('_')][0]}"].shape[0]
    first_var_key = [k for k in f['var'].keys() if not k.startswith("_")][0]
    fvi = f[f'var/{first_var_key}']
    n_var = fvi['codes'].shape[0] if (isinstance(fvi, h5py.Group) and 'codes' in fvi) else fvi.shape[0]
    X_full = sp.csr_matrix((data, indices, indptr), shape=(n_obs, n_var))
    X = X_full[idx].tocsr()
    obs = _read_df(f["obs"]).iloc[idx].reset_index(drop=True)
    var = _read_df(f["var"])

print(f"   subsampled counts shape: {X.shape}")
print(f"   X max={X.data.max():.2f}, min={X.data.min():.2f}, "
      f"integer-like={np.allclose(X.data, X.data.round())}")
adata = ad.AnnData(X=X, obs=obs, var=var)
adata.layers["counts"] = X.copy()

batch_col = next(c for c in ("dataset", "batch", "sample") if c in adata.obs.columns)
ct_col = next(c for c in ("cell_type", "celltype", "cell_type_original", "Cell.type")
              if c in adata.obs.columns)
print(f"   batch col = {batch_col!r} ({adata.obs[batch_col].nunique()} batches)")
print(f"   cell-type col = {ct_col!r} ({adata.obs[ct_col].nunique()} types)")

# ------------------------------------------------------------------
# 2. HVG selection + setup
# ------------------------------------------------------------------
print("\n[2] HVG selection (n=3000, Seurat v3 on counts) …")
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3",
                            batch_key=batch_col, subset=True)
print(f"   HVG-subset shape: {adata.shape}")

# ------------------------------------------------------------------
# 3. Train scVI
# ------------------------------------------------------------------
print("\n[3] Training scVI …")
import scvi
import torch
scvi.settings.seed = SEED
scvi.model.SCVI.setup_anndata(adata, batch_key=batch_col)
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30,
                        gene_likelihood="nb")
use_gpu = False  # CUDA driver on host is older than bundled torch; force CPU
print(f"   CUDA available but forcing CPU (driver mismatch): torch.cuda={torch.cuda.is_available()}")
model.train(max_epochs=25, batch_size=1024, early_stopping=True,
            accelerator="cpu",
            plan_kwargs=dict(lr=1e-3))
adata.obsm["X_scVI"] = model.get_latent_representation()
print(f"   scVI latent shape: {adata.obsm['X_scVI'].shape}")

# ------------------------------------------------------------------
# 4. Metric helpers (parity with 01c)
# ------------------------------------------------------------------
def _safe_sil(X, labels, sample=METRIC_N, seed=0):
    labels = np.asarray(labels)
    keep = pd.Series(labels).value_counts()
    keep = keep[keep >= 10].index
    m = np.isin(labels, keep)
    X = X[m]; labels = labels[m]
    if len(labels) > sample:
        rs = np.random.default_rng(seed)
        sel = rs.choice(len(labels), sample, replace=False)
        X = X[sel]; labels = labels[sel]
    if len(np.unique(labels)) < 2:
        return np.nan
    return float(silhouette_score(X, labels, metric="euclidean"))

def _knn(X, k=N_NN):
    nn = NearestNeighbors(n_neighbors=k+1, algorithm="auto").fit(X)
    return nn.kneighbors(X, return_distance=False)[:, 1:]

def _kbet(emb, batch, k=N_NN, sample=5000, seed=0):
    rng = np.random.default_rng(seed)
    n = emb.shape[0]
    sel = rng.choice(n, min(sample, n), replace=False)
    knn = _knn(emb[sel], k=k)
    cats, ci = np.unique(batch, return_inverse=True)
    gp = np.bincount(ci, minlength=len(cats)) / len(ci)
    reject = 0
    for i, nb in enumerate(knn):
        obs = np.bincount(ci[sel[nb]], minlength=len(cats))
        exp = gp * k
        mk = exp > 1
        if mk.sum() < 2: continue
        try:
            _, p = chisquare(obs[mk], exp[mk])
        except Exception:
            continue
        if p < 0.05: reject += 1
    return 1.0 - reject / len(knn)

def _lisi(emb, labels, k=N_NN, sample=5000, seed=0):
    rng = np.random.default_rng(seed)
    n = emb.shape[0]
    sel = rng.choice(n, min(sample, n), replace=False)
    knn = _knn(emb[sel], k=k)
    lab = np.asarray(labels)[sel]
    _, li = np.unique(lab, return_inverse=True)
    nl = li[knn]
    out = np.empty(len(sel))
    for i in range(len(sel)):
        c = np.bincount(nl[i], minlength=li.max()+1)
        p = c / c.sum(); p = p[p > 0]
        out[i] = 1.0 / np.sum(p**2)
    return out

# ------------------------------------------------------------------
# 5. Compute metrics on scVI latent
# ------------------------------------------------------------------
print("\n[4] Metrics on scVI latent …")
sv = adata.obsm["X_scVI"]
batch = adata.obs[batch_col].astype(str).values
ctype = adata.obs[ct_col].astype(str).values

sb = _safe_sil(sv, batch)
sc_ = _safe_sil(sv, ctype)
asw_batch = (sb + 1) / 2
asw_bio = (sc_ + 1) / 2
mixing = 1.0 - asw_batch
balance = (mixing + asw_bio) / 2
kb = _kbet(sv, batch)
ilisi = _lisi(sv, batch)
clisi = _lisi(sv, ctype)
print(f"   asw_batch={asw_batch:.3f} asw_bio={asw_bio:.3f} "
      f"mixing={mixing:.3f} balance={balance:.3f}")
print(f"   kBET={kb:.3f} iLISI_med={np.median(ilisi):.2f} "
      f"cLISI_med={np.median(clisi):.2f}")

row = dict(
    embedding="scVI",
    asw_batch=asw_batch,
    asw_bio=asw_bio,
    batch_mixing=mixing,
    asw_balance=balance,
    kBET_accept=kb,
    iLISI_mean=float(np.mean(ilisi)),
    iLISI_median=float(np.median(ilisi)),
    cLISI_mean=float(np.mean(clisi)),
    cLISI_median=float(np.median(clisi)),
)
df_scvi = pd.DataFrame([row])
df_scvi.to_csv(TAB / "scvi_vs_harmony_metrics.csv", index=False)
print(f"   wrote {TAB/'scvi_vs_harmony_metrics.csv'}")

# ------------------------------------------------------------------
# 6. Append to existing batch_effect_metrics.csv (non-destructive)
# ------------------------------------------------------------------
main = TAB / "batch_effect_metrics.csv"
if main.exists():
    m = pd.read_csv(main)
    if "scVI" not in m["embedding"].values:
        mnew = pd.concat([m, df_scvi], ignore_index=True)
        # Align columns
        mnew = mnew[[c for c in m.columns if c in mnew.columns] +
                    [c for c in mnew.columns if c not in m.columns]]
        mnew.to_csv(main, index=False)
        print(f"   appended scVI row to {main.name}")
    else:
        print(f"   scVI row already present in {main.name}; skipped append")

# ------------------------------------------------------------------
# 7. UMAP + metric bar figure
# ------------------------------------------------------------------
print("\n[5] UMAP on scVI …")
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
sc.tl.umap(adata)
fig, ax = plt.subplots(figsize=(6.2, 5.2), constrained_layout=True)
u = adata.obsm["X_umap"]
for b in sorted(adata.obs[batch_col].astype(str).unique()):
    m = adata.obs[batch_col].astype(str).values == b
    ax.scatter(u[m, 0], u[m, 1], s=1.2, alpha=0.35, label=b, edgecolor="none")
ax.set_title("scVI-integrated UMAP (60k subsample)",
             fontsize=10.5, fontweight="bold")
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.set_aspect("equal")
ax.legend(markerscale=6, fontsize=7, frameon=False,
          bbox_to_anchor=(1.02, 1.0), loc="upper left")
fig.savefig(FIG / "01d_scvi_umap.png", dpi=300, bbox_inches="tight",
            facecolor="white")
plt.close(fig)
print(f"   saved {FIG/'01d_scvi_umap.png'}")

# Combined bar with Harmony + PCA if available
try:
    main_df = pd.read_csv(TAB / "batch_effect_metrics.csv")
    keep_cols = ["embedding", "batch_mixing", "asw_bio",
                 "iLISI_median", "cLISI_median"]
    bar = main_df[keep_cols].melt(id_vars="embedding",
                                  var_name="metric", value_name="value")
    fig, ax = plt.subplots(figsize=(8.6, 4.2), constrained_layout=True)
    import seaborn as sns
    pal = {"PCA": "#B0B0B0", "Harmony": "#117733", "scVI": "#332288"}
    sns.barplot(data=bar, x="metric", y="value", hue="embedding",
                palette=pal, ax=ax, edgecolor="black", linewidth=0.5)
    ax.set_title("PCA vs Harmony vs scVI — batch-effect metrics",
                 fontsize=10.5, fontweight="bold")
    ax.tick_params(axis="x", rotation=20); ax.set_xlabel("")
    ax.legend(frameon=False, fontsize=8, title="Embedding")
    fig.savefig(FIG / "01d_scvi_metric_bars.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"   saved {FIG/'01d_scvi_metric_bars.png'}")
except Exception as e:
    print(f"   bar figure skipped: {e}")

print("\n[01d] DONE")
