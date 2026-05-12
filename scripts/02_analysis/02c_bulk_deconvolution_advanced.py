#!/usr/bin/env python3
"""02c_bulk_deconvolution_advanced.py

Advanced bulk-deconvolution analysis pack.

Inputs:
  - data/processed/integrated_io_atlas.h5ad           (sc atlas; 1.06M cells, 79 patients)
  - out/02c_bulk_deconvolution/tables/cibersortx_hcc_reference.tsv   (CIBERSORTx ref: 15,697 × 43)
  - data/clinical_trials/IMbrave150/imbrave150_fig3b_patient_data.csv (177 patients with CD8 % + survival)

Pipeline:
  1. Build per-patient pseudobulk from the atlas (counts summed → log1p+CPM).
  2. NNLS deconvolution of each patient pseudobulk against the 43-cell-type reference.
  3. Compare predicted vs ground-truth (atlas-derived) fractions per coarse lineage.
  4. Plot panel set (UMAP, stacked bars, boxplots, heatmap, radar, KM survival).

Outputs:
  out/02c_bulk_deconvolution/tables/pseudobulk_per_patient.csv
  out/02c_bulk_deconvolution/tables/deconv_fractions_per_patient.csv
  out/02c_bulk_deconvolution/tables/deconv_RvsNR_stats.csv
  out/02c_bulk_deconvolution/tables/imbrave150_km_summary.csv
  out/02c_bulk_deconvolution/figures/bulk_deconv_advanced_composite.{png,pdf}
"""
from __future__ import annotations

from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from scipy.optimize import nnls
from scipy.stats import mannwhitneyu
from sklearn.preprocessing import StandardScaler
import umap as umap_

PROJECT = Path(__file__).resolve().parents[2]
ATLAS   = PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad"
REF_TSV = PROJECT / "out" / "02c_bulk_deconvolution" / "tables" / "cibersortx_hcc_reference.tsv"
IMB     = PROJECT / "data" / "clinical_trials" / "IMbrave150" / "imbrave150_fig3b_patient_data.csv"
OUT_TAB = PROJECT / "out" / "02c_bulk_deconvolution" / "tables"
OUT_FIG = PROJECT / "out" / "02c_bulk_deconvolution" / "figures"
OUT_TAB.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

# Coarse lineage map for the 43 reference cell types
LINEAGE = {
    "B_Memory": "B", "B_Naive": "B", "B_Plasma": "B",
    "CD4_Effector-Memory": "CD4 T", "CD4_Memory": "CD4 T", "CD4_Tfh-like": "CD4 T",
    "CD4_Th1 GZMK": "CD4 T", "CD4_Th1/Th17": "CD4 T",
    "CD8_Cytotoxic": "CD8 T", "CD8_Effector": "CD8 T", "CD8_Memory": "CD8 T",
    "CD8_Dysfunctional-effector": "CD8 T", "CD8_Dysfunctional-Proliferating": "CD8 T",
    "CD8_Dysfunctional-progenitor": "CD8 T", "CD8_Dysfunctional-terminal": "CD8 T",
    "DC_cDC1": "DC", "DC_cDC2": "DC", "DC_mregDC": "DC", "DC_pDC": "DC",
    "Macrophage_FCN1/CXCL9": "Macrophage", "Macrophage_FOLR2": "Macrophage",
    "Macrophage_I": "Macrophage", "Macrophage_II": "Macrophage",
    "Mast": "Other", "MAIT": "T (other)", "ILC": "Other",
    "Monocyte_Classical": "Macrophage", "Monocyte_Non-classical": "Macrophage",
    "NK_CD16+": "NK", "NK_CD56+": "NK",
    "Treg": "Treg", "gdT": "T (other)", "iNKT": "T (other)",
    "CD8_Exhausted": "CD8 T", "CD8_GZMK": "CD8 T", "Neutrophil": "Other",
    "Endothelial": "Stroma", "Fibroblast": "Stroma", "Hepatocyte": "Hepatocyte",
    "Cholangiocyte": "Stroma", "Smooth_Muscle": "Stroma",
    "Malignant": "Malignant", "Unknown_Doublets": "Other",
}

# ─────────────────────────────────────────────────────────────────────────────
# Build per-patient pseudobulk
# ─────────────────────────────────────────────────────────────────────────────
print("Streaming atlas to build pseudobulk …")
adata = ad.read_h5ad(ATLAS, backed="r")
patients = adata.obs["patient_ID"].astype(str).values
responses = adata.obs["response"].astype(str).values
ct        = adata.obs["cell_type"].astype(str).values
genes     = adata.var_names.values

uniq_pat, pat_idx = np.unique(patients, return_inverse=True)
n_p = len(uniq_pat)
n_g = adata.shape[1]
print(f"  patients={n_p}  genes={n_g}")

# Stream chunks → sum per patient
PB_counts = np.zeros((n_p, n_g), dtype=np.float64)
PB_cells  = np.zeros(n_p, dtype=np.int64)
cache = OUT_TAB / "_pseudobulk_counts.npz"
if cache.exists():
    print("  loading cached pseudobulk counts …")
    z = np.load(cache, allow_pickle=True)
    if list(z["patients"]) == list(uniq_pat) and list(z["genes"]) == list(genes):
        PB_counts = z["counts"]
        PB_cells  = z["ncells"]
    else:
        print("  cache schema mismatch → recomputing")
        cache.unlink()
if not cache.exists():
    chunk = 50_000
    for s in range(0, adata.shape[0], chunk):
        e = min(s + chunk, adata.shape[0])
        X = adata.X[s:e]
        if sp.issparse(X):
            X = X.toarray()
        pi = pat_idx[s:e]
        np.add.at(PB_counts, pi, X)
        np.add.at(PB_cells, pi, 1)
        if (s // chunk) % 5 == 0:
            print(f"    [{e:,}/{adata.shape[0]:,}]", flush=True)
    np.savez_compressed(cache, counts=PB_counts, ncells=PB_cells,
                        patients=np.array(uniq_pat), genes=np.array(genes))
    print(f"  cached → {cache.name}")

PB = pd.DataFrame(PB_counts, index=uniq_pat, columns=genes)

# Patient-level meta (response: take majority, treatment: majority)
meta = pd.DataFrame(index=uniq_pat)
def majority(arr, idxs):
    vals, c = np.unique(arr[idxs], return_counts=True)
    return vals[np.argmax(c)]
resp_pat, ncells_pat = [], []
for p in uniq_pat:
    m = patients == p
    resp_pat.append(majority(responses, np.where(m)[0]))
    ncells_pat.append(int(m.sum()))
meta["response"] = resp_pat
meta["n_cells"]  = ncells_pat
PB.insert(0, "n_cells", meta["n_cells"].values)
PB.insert(0, "response", meta["response"].values)
PB.to_csv(OUT_TAB / "pseudobulk_per_patient.csv")
print(f"  → pseudobulk_per_patient.csv ({PB.shape})")

PB_mat = PB.iloc[:, 2:].values            # patients × genes (counts)
# CPM + log1p
lib = PB_mat.sum(axis=1, keepdims=True)
PB_cpm = np.log1p(PB_mat / np.maximum(lib, 1) * 1e6)

# Ground-truth coarse fractions (from atlas cell_type column)
gt_frac = pd.crosstab(patients, ct, normalize="index")
gt_frac = gt_frac.loc[uniq_pat]

# ─────────────────────────────────────────────────────────────────────────────
# NNLS deconvolution
# ─────────────────────────────────────────────────────────────────────────────
print("\nNNLS deconvolution …")
ref = pd.read_csv(REF_TSV, sep="\t", index_col=0)
common = [g for g in ref.index if g in PB.columns[2:]]
print(f"  ref overlap: {len(common):,} genes")
R = ref.loc[common].values                # genes × ctypes
B = np.log1p(PB_mat[:, [list(genes).index(g) for g in common]] /
             np.maximum(lib, 1) * 1e6)    # patients × genes (log1p CPM)

# z-scale rows of R and B (gene-wise) to balance contributions
sc_R = StandardScaler().fit(R.T).transform(R.T).T  # standardize across cell types
sc_B = (B - B.mean(axis=0)) / (B.std(axis=0) + 1e-6)

frac = np.zeros((n_p, R.shape[1]))
for i in range(n_p):
    sol, _ = nnls(sc_R, sc_B[i])
    if sol.sum() > 0:
        frac[i] = sol / sol.sum()
    else:
        # fallback: positive only NNLS on raw
        sol, _ = nnls(R, B[i])
        s = sol.sum()
        frac[i] = sol / s if s > 0 else 0
frac_df = pd.DataFrame(frac, index=uniq_pat, columns=ref.columns)
frac_df.insert(0, "response", meta["response"].values)
frac_df.to_csv(OUT_TAB / "deconv_fractions_per_patient.csv")
print(f"  → deconv_fractions_per_patient.csv ({frac_df.shape})")

# Aggregate to coarse lineage for benchmark
lin_pred = pd.DataFrame(index=uniq_pat)
lin_pred_cols = {}
for lin in set(LINEAGE.values()):
    cols = [c for c in ref.columns if LINEAGE.get(c, "Other") == lin]
    if not cols: continue
    lin_pred_cols[lin] = cols
    lin_pred[lin] = frac_df[cols].sum(axis=1)

# ─────────────────────────────────────────────────────────────────────────────
# R vs NR stats (use Responder / Non-Responder only)
# ─────────────────────────────────────────────────────────────────────────────
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); rk = p[o]
    adj = rk * n / (np.arange(n) + 1); adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(adj, 0, 1); return out

mask_r = meta["response"].values == "Responder"
mask_nr = meta["response"].values == "Non-Responder"
print(f"  R={mask_r.sum()}  NR={mask_nr.sum()}")

stat_rows = []
for col in ref.columns:
    r_v = frac_df.loc[mask_r, col].values
    nr_v = frac_df.loc[mask_nr, col].values
    if len(r_v) < 3 or len(nr_v) < 3 or (r_v.sum() + nr_v.sum()) == 0:
        continue
    try:
        U, p = mannwhitneyu(r_v, nr_v, alternative="two-sided")
    except ValueError:
        continue
    stat_rows.append(dict(
        cell_type=col, lineage=LINEAGE.get(col, "Other"),
        median_R=float(np.median(r_v)), median_NR=float(np.median(nr_v)),
        delta=float(np.median(r_v) - np.median(nr_v)),
        n_R=int(mask_r.sum()), n_NR=int(mask_nr.sum()), p=float(p),
    ))
stat_df = pd.DataFrame(stat_rows)
if len(stat_df):
    stat_df["padj"] = bh(stat_df["p"])
    stat_df = stat_df.sort_values("padj")
stat_df.to_csv(OUT_TAB / "deconv_RvsNR_stats.csv", index=False)
print(f"  → deconv_RvsNR_stats.csv ({len(stat_df)} cell types)")
if len(stat_df):
    print(stat_df.head(10).to_string(index=False))

# ─────────────────────────────────────────────────────────────────────────────
# IMbrave150 survival KM (CD8 H vs L) per arm — use provided columns
# ─────────────────────────────────────────────────────────────────────────────
imb = pd.read_csv(IMB)
km_rows = []
def km_curve(t, e):
    """Return arrays (sorted_t, surv, ci_lo, ci_hi) — Kaplan-Meier with Greenwood SE."""
    order = np.argsort(t)
    t = t[order]; e = e[order]
    n_at = len(t)
    surv = 1.0
    var  = 0.0
    out_t, out_s, out_lo, out_hi = [0.0], [1.0], [1.0], [1.0]
    i = 0
    while i < len(t):
        j = i
        while j < len(t) and t[j] == t[i]: j += 1
        d = e[i:j].sum()
        if n_at > 0 and d > 0:
            surv  *= (1 - d / n_at)
            var   += d / (n_at * (n_at - d)) if n_at != d else 0
        out_t.append(t[i]); out_s.append(surv)
        if surv > 0:
            ci = 1.96 * np.sqrt(var) * surv
        else:
            ci = 0
        out_lo.append(max(0.0, surv - ci)); out_hi.append(min(1.0, surv + ci))
        n_at -= (j - i)
        i = j
    return np.array(out_t), np.array(out_s), np.array(out_lo), np.array(out_hi)

# Logrank
def logrank(t1, e1, t2, e2):
    ts = np.unique(np.concatenate([t1[e1 == 1], t2[e2 == 1]]))
    O1, E1, V = 0.0, 0.0, 0.0
    for t in ts:
        n1 = (t1 >= t).sum(); n2 = (t2 >= t).sum()
        d1 = ((t1 == t) & (e1 == 1)).sum()
        d2 = ((t2 == t) & (e2 == 1)).sum()
        n = n1 + n2; d = d1 + d2
        if n <= 1: continue
        e = d * n1 / n
        v = d * (n - d) * n1 * n2 / (n * n * (n - 1))
        O1 += d1; E1 += e; V += v
    if V == 0: return 1.0, 0.0
    z = (O1 - E1) / np.sqrt(V)
    from math import erfc
    pval = erfc(abs(z) / np.sqrt(2))
    return pval, z

km_arm_data = {}
for arm in imb["arm"].unique():
    sub = imb[imb["arm"] == arm]
    h = sub[sub["cd8_hl"] == "High"]
    l = sub[sub["cd8_hl"] == "Low"]
    if len(h) < 5 or len(l) < 5: continue
    p, z = logrank(h["os_m"].values, h["os_event"].values, l["os_m"].values, l["os_event"].values)
    km_rows.append(dict(arm=arm, endpoint="OS", n_high=len(h), n_low=len(l),
                        events_high=int(h["os_event"].sum()), events_low=int(l["os_event"].sum()),
                        logrank_p=float(p), z=float(z)))
    km_arm_data[arm] = dict(
        h_t=h["os_m"].values, h_e=h["os_event"].values,
        l_t=l["os_m"].values, l_e=l["os_event"].values,
    )
km_df = pd.DataFrame(km_rows)
km_df.to_csv(OUT_TAB / "imbrave150_km_summary.csv", index=False)
print(f"  → imbrave150_km_summary.csv ({len(km_df)} arms)")
print(km_df.to_string(index=False))

# ─────────────────────────────────────────────────────────────────────────────
# UMAP on patient × cell-type fractions
# ─────────────────────────────────────────────────────────────────────────────
print("\nUMAP on deconvolved fractions …")
clrs = StandardScaler().fit_transform(frac)
reducer = umap_.UMAP(n_neighbors=10, min_dist=0.4, random_state=0)
emb = reducer.fit_transform(clrs)
emb_df = pd.DataFrame(emb, index=uniq_pat, columns=["UMAP1", "UMAP2"])
emb_df["response"] = meta["response"].values
emb_df.to_csv(OUT_TAB / "deconv_patient_umap.csv")

# ─────────────────────────────────────────────────────────────────────────────
# Composite figure (4×3)
# ─────────────────────────────────────────────────────────────────────────────
print("Building composite figure …")
fig = plt.figure(figsize=(20, 16), dpi=180)
gs = GridSpec(4, 3, figure=fig, wspace=0.32, hspace=0.55)

resp_pal = {"Responder": "#2E86C1", "Non-Responder": "#C0392B", "Unknown": "#bbbbbb"}

# A) Patient UMAP
ax = fig.add_subplot(gs[0, 0])
for r, c in resp_pal.items():
    m = emb_df["response"] == r
    if m.sum():
        ax.scatter(emb_df.loc[m, "UMAP1"], emb_df.loc[m, "UMAP2"],
                   c=c, s=44, edgecolor="black", linewidth=0.5, label=f"{r} (n={m.sum()})")
ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
ax.set_title("A. Patient pseudobulk UMAP\n(deconv-fraction space)", fontsize=10)
ax.legend(fontsize=7, frameon=False, loc="best")

# B) Stacked-bar fractions (top-15 cell types) ordered by R / NR
ax = fig.add_subplot(gs[0, 1:])
order_pat = list(np.where(meta["response"].values == "Responder")[0]) + \
            list(np.where(meta["response"].values == "Non-Responder")[0]) + \
            list(np.where(~np.isin(meta["response"].values, ["Responder", "Non-Responder"]))[0])
top15 = ref.columns[np.argsort(-frac.sum(axis=0))[:15]]
others = [c for c in ref.columns if c not in top15]
plot_mat = np.column_stack([frac_df[c].values for c in top15] + [frac_df[others].sum(axis=1).values])
plot_cols = list(top15) + ["Other"]
plot_mat = plot_mat[order_pat]
bottom = np.zeros(len(order_pat))
ct_pal = plt.cm.tab20(np.linspace(0, 1, len(plot_cols)))
for i, c in enumerate(plot_cols):
    ax.bar(np.arange(len(order_pat)), plot_mat[:, i],
           bottom=bottom, color=ct_pal[i], label=c, width=1.0,
           edgecolor="white", linewidth=0.1)
    bottom += plot_mat[:, i]
n_R = int((meta["response"].values == "Responder").sum())
n_NR = int((meta["response"].values == "Non-Responder").sum())
ax.axvline(n_R - 0.5, color="black", linewidth=0.7)
ax.axvline(n_R + n_NR - 0.5, color="black", linewidth=0.7, linestyle="--")
ax.set_xlim(-0.5, len(order_pat) - 0.5)
ax.set_ylim(0, 1)
ax.set_xticks([n_R / 2, n_R + n_NR / 2, n_R + n_NR + (n_p - n_R - n_NR) / 2])
ax.set_xticklabels(["Responder", "Non-Responder", "Unknown"])
ax.set_ylabel("Fraction")
ax.set_title(f"B. Per-patient cell-type composition (NNLS, top-15 + Other) — {n_p} patients", fontsize=10)
ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=6, frameon=False, ncol=1)

# C) Boxplot of top-8 differentially abundant cell types R vs NR
ax = fig.add_subplot(gs[1, :2])
top_diff = stat_df.sort_values("p").head(8)["cell_type"].values if len(stat_df) else top15[:8]
positions = np.arange(len(top_diff))
for i, ctn in enumerate(top_diff):
    r_v  = frac_df.loc[mask_r, ctn].values
    nr_v = frac_df.loc[mask_nr, ctn].values
    bp = ax.boxplot([r_v, nr_v], positions=[i - 0.18, i + 0.18], widths=0.32,
                    patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], [resp_pal["Responder"], resp_pal["Non-Responder"]]):
        patch.set_facecolor(c); patch.set_alpha(0.7); patch.set_edgecolor("black")
    ax.scatter(np.full_like(r_v, i - 0.18, dtype=float) + np.random.uniform(-0.06, 0.06, len(r_v)),
               r_v, c="black", s=8, alpha=0.5, zorder=3)
    ax.scatter(np.full_like(nr_v, i + 0.18, dtype=float) + np.random.uniform(-0.06, 0.06, len(nr_v)),
               nr_v, c="black", s=8, alpha=0.5, zorder=3)
ax.set_xticks(positions)
ax.set_xticklabels(top_diff, rotation=35, ha="right", fontsize=8)
ax.set_ylabel("Predicted fraction")
ax.set_title("C. Top-8 differentially abundant cell types (R vs NR, MWU)", fontsize=10)
handles = [Line2D([0],[0],marker="s",color="w",markerfacecolor=resp_pal["Responder"],
                  markeredgecolor="black",label="R",markersize=10),
           Line2D([0],[0],marker="s",color="w",markerfacecolor=resp_pal["Non-Responder"],
                  markeredgecolor="black",label="NR",markersize=10)]
ax.legend(handles=handles, fontsize=8, frameon=False, loc="best")

# D) Volcano (-log10 p vs delta)
ax = fig.add_subplot(gs[1, 2])
if len(stat_df):
    x = stat_df["delta"].values
    y = -np.log10(stat_df["p"].clip(lower=1e-10))
    sig = stat_df["padj"] < 0.1
    ax.scatter(x[~sig.values], y[~sig.values], s=18, c="#7f7f7f", alpha=0.6)
    ax.scatter(x[sig.values], y[sig.values], s=42, c="#C0392B", edgecolor="black", linewidth=0.4, label="BH<0.1")
    for _, row in stat_df.head(6).iterrows():
        ax.annotate(row["cell_type"], (row["delta"], -np.log10(row["p"])),
                    fontsize=6, xytext=(3, 2), textcoords="offset points")
    ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=0.6)
    ax.axvline(0, color="black", linewidth=0.4)
    ax.set_xlabel("Δ median fraction (R − NR)")
    ax.set_ylabel("−log10 p")
    ax.set_title("D. R vs NR volcano", fontsize=10)
    ax.legend(fontsize=7, frameon=False)

# E) Heatmap: cell-type × patient (clustered cols by response)
ax = fig.add_subplot(gs[2, :2])
order_idx = np.array(order_pat)
H = frac.T[:, order_idx]   # ctypes × patients
H = (H - H.mean(axis=1, keepdims=True)) / (H.std(axis=1, keepdims=True) + 1e-6)
im = ax.imshow(H, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2,
               interpolation="nearest", rasterized=True)
ax.set_yticks(np.arange(len(ref.columns)))
ax.set_yticklabels(ref.columns, fontsize=5)
ax.set_xticks([n_R / 2, n_R + n_NR / 2, n_R + n_NR + (n_p - n_R - n_NR) / 2])
ax.set_xticklabels(["R", "NR", "Unknown"], fontsize=8)
ax.axvline(n_R - 0.5, color="black", linewidth=0.7)
ax.axvline(n_R + n_NR - 0.5, color="black", linewidth=0.7, linestyle="--")
ax.set_title("E. Per-patient z-scored cell-type fractions (43 types × patients)", fontsize=10)
fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01, label="z(fraction)")

# F) Radar: mean coarse-lineage fractions R vs NR
ax = fig.add_subplot(gs[2, 2], projection="polar")
lin_keys = sorted(set(LINEAGE.values()))
mean_R  = [lin_pred.loc[mask_r, k].mean() if k in lin_pred.columns else 0.0 for k in lin_keys]
mean_NR = [lin_pred.loc[mask_nr, k].mean() if k in lin_pred.columns else 0.0 for k in lin_keys]
angles = np.linspace(0, 2 * np.pi, len(lin_keys), endpoint=False)
mean_R  = mean_R  + [mean_R[0]]
mean_NR = mean_NR + [mean_NR[0]]
angles  = np.concatenate([angles, [angles[0]]])
ax.plot(angles, mean_R, color=resp_pal["Responder"], linewidth=1.5, label="R")
ax.fill(angles, mean_R, color=resp_pal["Responder"], alpha=0.18)
ax.plot(angles, mean_NR, color=resp_pal["Non-Responder"], linewidth=1.5, label="NR")
ax.fill(angles, mean_NR, color=resp_pal["Non-Responder"], alpha=0.18)
ax.set_xticks(angles[:-1])
ax.set_xticklabels(lin_keys, fontsize=7)
ax.set_yticklabels([])
ax.set_title("F. Mean lineage fractions", fontsize=10, pad=15)
ax.legend(loc="lower right", bbox_to_anchor=(1.25, -0.05), fontsize=7, frameon=False)

# G) Cell-type co-abundance correlation heatmap
ax = fig.add_subplot(gs[3, 0])
corr = np.corrcoef(frac.T)
im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest", rasterized=True)
ax.set_xticks([]); ax.set_yticks([])
ax.set_title("G. Cell-type co-abundance correlations\n(43×43, across 79 patients)", fontsize=10)
fig.colorbar(im, ax=ax, fraction=0.045, label="Pearson r")

# H) Pred vs ground-truth scatter (lineage fractions)
ax = fig.add_subplot(gs[3, 1])
common_lin = [lin for lin in lin_pred.columns if lin in gt_frac.columns or lin in {"B","CD4 T","CD8 T","NK","Treg","Macrophage","DC","Hepatocyte","Malignant"}]
gt_coarse = pd.DataFrame(index=uniq_pat)
mapping = {
    "B": ["B_cell", "Plasma"], "CD4 T": ["CD4_T"], "CD8 T": ["CD8_T"],
    "NK": ["NK"], "Treg": ["Treg"], "Macrophage": ["Macrophage"],
    "DC": ["DC", "cDC1", "cDC2", "mregDC"],
    "Hepatocyte": ["Hepatocyte"], "Malignant": ["Malignant"],
}
for lin, atlas_keys in mapping.items():
    cols = [c for c in atlas_keys if c in gt_frac.columns]
    if cols:
        gt_coarse[lin] = gt_frac[cols].sum(axis=1)
    else:
        gt_coarse[lin] = 0
both_lins = [l for l in mapping if l in lin_pred.columns and l in gt_coarse.columns]
xs, ys, labels = [], [], []
for lin in both_lins:
    a = gt_coarse[lin].reindex(uniq_pat).fillna(0).values
    b = lin_pred[lin].reindex(uniq_pat).fillna(0).values
    xs.extend(a); ys.extend(b); labels.extend([lin] * len(a))
xs = np.array(xs); ys = np.array(ys); labels = np.array(labels)
lin_color = {l: c for l, c in zip(both_lins, plt.cm.tab10(np.linspace(0, 1, len(both_lins))))}
for lin in both_lins:
    m = labels == lin
    ax.scatter(xs[m], ys[m], c=[lin_color[lin]], s=22, alpha=0.7,
               edgecolor="black", linewidth=0.3, label=lin)
mx = max(xs.max(), ys.max()) * 1.05
ax.plot([0, mx], [0, mx], "k--", linewidth=0.6)
ax.set_xlabel("Ground truth fraction (atlas)")
ax.set_ylabel("Deconvolved fraction (NNLS)")
ax.set_xlim(0, mx); ax.set_ylim(0, mx)
r2 = np.corrcoef(xs, ys)[0, 1] ** 2
ax.set_title(f"H. Deconv vs ground truth (r²={r2:.2f})", fontsize=10)
ax.legend(fontsize=6, frameon=False, loc="best")

# I) IMbrave150 KM curves
ax = fig.add_subplot(gs[3, 2])
arm_styles = {"Atezo+Bev": "-", "Sorafenib": "--"}
for arm, d in km_arm_data.items():
    ls = arm_styles.get(arm, "-")
    t, s, lo, hi = km_curve(d["h_t"], d["h_e"])
    ax.step(t, s, where="post", color="#27AE60", linestyle=ls, linewidth=1.5,
            label=f"{arm} CD8-Hi (n={len(d['h_t'])})")
    t, s, lo, hi = km_curve(d["l_t"], d["l_e"])
    ax.step(t, s, where="post", color="#C0392B", linestyle=ls, linewidth=1.5,
            label=f"{arm} CD8-Lo (n={len(d['l_t'])})")
ax.set_xlabel("OS (months)"); ax.set_ylabel("Survival probability")
ax.set_ylim(0, 1.02)
title_p = "  ".join([f"{r['arm']} p={r['logrank_p']:.3f}" for _, r in km_df.iterrows()])
ax.set_title(f"I. IMbrave150 OS by CD8 stratum\n{title_p}", fontsize=9)
ax.legend(fontsize=6, frameon=False, loc="lower left")

fig.suptitle(
    f"Bulk-deconvolution advanced analyses — {n_p} pseudobulk patients × "
    f"{R.shape[1]} reference cell types  +  IMbrave150 (n={len(imb)})",
    fontsize=14, y=0.995,
)

png = OUT_FIG / "bulk_deconv_advanced_composite.png"
pdf = OUT_FIG / "bulk_deconv_advanced_composite.pdf"
fig.savefig(png, dpi=200, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
plt.close(fig)
print(f"  → {png.name}")
print(f"  → {pdf.name}")
print("Done.")
