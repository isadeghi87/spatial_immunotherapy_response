#!/usr/bin/env python3
"""02g_infercnv_spatial_universal_heatmap.py

Stack per-sample inferCNVpy npz outputs into a single universal heatmap:
  rows  = spots across all 23 Visium samples (≈43k)
  cols  = ~2,115 genome-window CNV bins (shared chr_pos axis)
  row tracks = response, dataset, sample, cell_type, tumor_status, cnv_score

Outputs:
  out/02b_spatial_biomarkers/cnv/figures/spatial_cnv_universal_heatmap.{png,pdf}
  out/02b_spatial_biomarkers/cnv/tables/spatial_cnv_universal_spotmeta.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec

PROJECT  = Path(__file__).resolve().parents[2]
NPZ_DIR  = PROJECT / "out" / "02b_spatial_biomarkers" / "cnv" / "npz"
FIG_DIR  = PROJECT / "out" / "02b_spatial_biomarkers" / "cnv" / "figures"
TAB_DIR  = PROJECT / "out" / "02b_spatial_biomarkers" / "cnv" / "tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)

CLIP = 0.3            # CNV value clip
TUMOR_PCT = 75.0      # within-sample percentile of cnv_score → tumor

# ─────────────────────────────────────────────────────────────────────────────
# Load and stack
# ─────────────────────────────────────────────────────────────────────────────
files = sorted(NPZ_DIR.glob("*_cnv.npz"))
print(f"Loading {len(files)} npz files …")

X_blocks, meta_rows = [], []
chr_keys, chr_vals, n_windows = None, None, None

for f in files:
    d = np.load(f, allow_pickle=True)
    X = d["X_cnv"]
    if chr_keys is None:
        chr_keys = list(d["chr_pos_keys"])
        chr_vals = list(d["chr_pos_values"])
        n_windows = X.shape[1]
    else:
        assert list(d["chr_pos_keys"]) == chr_keys, f"chr mismatch: {f.name}"
        assert X.shape[1] == n_windows, f"window mismatch: {f.name}"
    X_blocks.append(X)
    meta_rows.append(pd.DataFrame({
        "spot_id":    d["barcodes"],
        "sample_id":  str(d["sample_id"]),
        "dataset":    str(d["dataset"]),
        "response":   str(d["response"]),
        "cell_type":  d["cell_types"],
        "group":      d["groups"],
        "cnv_score":  d["cnv_score"].astype(float),
        "cnv_leiden": d["cnv_leiden"],
    }))

X_all  = np.vstack(X_blocks).astype(np.float32)
meta   = pd.concat(meta_rows, ignore_index=True)
print(f"Stacked: spots={X_all.shape[0]:,}  windows={X_all.shape[1]}")

# ─────────────────────────────────────────────────────────────────────────────
# Per-sample tumor classification (within-sample percentile of cnv_score on
# observation spots; reference spots stay non-tumor by definition).
# ─────────────────────────────────────────────────────────────────────────────
meta["tumor_status"] = "non_tumor"
for sid, sub in meta.groupby("sample_id"):
    obs = sub[sub["group"] == "observation"]
    if len(obs) == 0:
        continue
    thr = np.percentile(obs["cnv_score"].values, TUMOR_PCT)
    tumor_mask = (meta["sample_id"] == sid) & (meta["group"] == "observation") & (meta["cnv_score"] >= thr)
    meta.loc[tumor_mask, "tumor_status"] = "tumor"

# Save spot metadata
meta_out = TAB_DIR / "spatial_cnv_universal_spotmeta.csv"
meta.to_csv(meta_out, index=False)
print(f"  → {meta_out.name} ({len(meta):,} rows)")

# ─────────────────────────────────────────────────────────────────────────────
# Order rows: response (R, NR, NA) → tumor_status (tumor first) → sample → cnv_score
# ─────────────────────────────────────────────────────────────────────────────
resp_order  = {"Responder": 0, "Non_Responder": 1}
meta["_resp_rank"]  = meta["response"].map(resp_order).fillna(2).astype(int)
meta["_tum_rank"]   = (meta["tumor_status"] != "tumor").astype(int)  # tumor first
meta["_orig_index"] = np.arange(len(meta))
order = meta.sort_values(
    ["_resp_rank", "_tum_rank", "sample_id", "cnv_score"],
    ascending=[True, True, True, False],
).index.values

X_ord = X_all[order]
m_ord = meta.iloc[order].reset_index(drop=True)

# ─────────────────────────────────────────────────────────────────────────────
# Color encodings
# ─────────────────────────────────────────────────────────────────────────────
def encode(values, palette):
    cats = list(palette.keys())
    cat_to_int = {c: i for i, c in enumerate(cats)}
    arr = np.array([cat_to_int.get(v, len(cats)) for v in values])
    cmap = ListedColormap([palette[c] for c in cats] + ["#dddddd"])
    norm = BoundaryNorm(np.arange(len(cats) + 2) - 0.5, cmap.N)
    return arr, cmap, norm, cats

resp_pal = {"Responder": "#2E86C1", "Non_Responder": "#C0392B"}
tum_pal  = {"tumor": "#7B241C", "non_tumor": "#85929E"}

# Datasets
datasets = sorted(m_ord["dataset"].unique())
ds_pal   = {ds: c for ds, c in zip(datasets, plt.cm.Set2.colors)}

# Samples (one color per sample, alternating dataset hues)
samples = list(m_ord["sample_id"].drop_duplicates())
sample_pal = {s: c for s, c in zip(samples, plt.cm.tab20(np.linspace(0, 1, len(samples))))}

# Cell types — collapse to coarse grain to keep legend usable
def coarse(ct):
    s = str(ct)
    sl = s.lower().replace("+", "").replace(" ", "_")
    if "cd8" in sl: return "CD8 T"
    if "cd4" in sl: return "CD4 T"
    if "treg" in sl: return "Treg"
    if "nk"   in sl: return "NK"
    if "b_cell" in sl or sl == "b": return "B"
    if "macro" in sl or "kupffer" in sl: return "Macrophage"
    if "dc" in sl or "dendr" in sl: return "DC"
    if "endo" in sl: return "Endothelial"
    if "fibro" in sl or "stromal" in sl: return "Fibroblast"
    if "hepato" in sl: return "Hepatocyte"
    if "tumor"  in sl or "malign" in sl or "cancer" in sl: return "Malignant"
    if "chol" in sl: return "Cholangiocyte"
    return "Other"

m_ord["cell_coarse"] = m_ord["cell_type"].apply(coarse)
ct_pal = {
    "CD8 T":"#1F77B4","CD4 T":"#17BECF","Treg":"#9467BD","NK":"#2CA02C",
    "B":"#8C564B","Macrophage":"#D62728","DC":"#FF7F0E","Endothelial":"#7F7F7F",
    "Fibroblast":"#BCBD22","Hepatocyte":"#E5C07B","Malignant":"#7B241C",
    "Cholangiocyte":"#A569BD","Other":"#CCCCCC",
}

# ─────────────────────────────────────────────────────────────────────────────
# Figure
# ─────────────────────────────────────────────────────────────────────────────
n_spots = X_ord.shape[0]
fig_h = 12.0
fig_w = 16.0

fig = plt.figure(figsize=(fig_w, fig_h), dpi=200)
gs = GridSpec(
    2, 8,
    width_ratios=[0.18, 0.18, 0.18, 0.18, 0.18, 0.55, 16, 0.5],
    height_ratios=[0.5, 16],
    wspace=0.05, hspace=0.02,
)
# Track positions:  resp | tumor | dataset | sample | celltype | cnv_score | heatmap | colorbar
ax_resp  = fig.add_subplot(gs[1, 0])
ax_tum   = fig.add_subplot(gs[1, 1])
ax_ds    = fig.add_subplot(gs[1, 2])
ax_samp  = fig.add_subplot(gs[1, 3])
ax_ct    = fig.add_subplot(gs[1, 4])
ax_score = fig.add_subplot(gs[1, 5])
ax_main  = fig.add_subplot(gs[1, 6])
ax_cbar  = fig.add_subplot(gs[1, 7])
ax_chr   = fig.add_subplot(gs[0, 6])

def draw_track(ax, vals, palette, title):
    arr, cmap, norm, cats = encode(vals, palette)
    ax.imshow(arr.reshape(-1, 1), aspect="auto", cmap=cmap, norm=norm,
              interpolation="nearest", rasterized=True)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(title, rotation=90, fontsize=7, ha="center", va="top", labelpad=4)
    for s in ax.spines.values(): s.set_visible(False)

draw_track(ax_resp, m_ord["response"].values, resp_pal,        "Response")
draw_track(ax_tum,  m_ord["tumor_status"].values, tum_pal,     "Tumor")
draw_track(ax_ds,   m_ord["dataset"].values, ds_pal,           "Dataset")
draw_track(ax_samp, m_ord["sample_id"].values, sample_pal,     "Sample")
draw_track(ax_ct,   m_ord["cell_coarse"].values, ct_pal,       "Cell type")

# CNV score (continuous)
ax_score.imshow(m_ord["cnv_score"].values.reshape(-1, 1), aspect="auto",
                cmap="viridis", interpolation="nearest", rasterized=True)
ax_score.set_xticks([]); ax_score.set_yticks([])
ax_score.set_xlabel("CNV score", rotation=90, fontsize=7, ha="center", va="top", labelpad=4)
for s in ax_score.spines.values(): s.set_visible(False)

# Main heatmap — CNV matrix
im = ax_main.imshow(
    np.clip(X_ord, -CLIP, CLIP),
    aspect="auto", cmap="RdBu_r", vmin=-CLIP, vmax=CLIP,
    interpolation="nearest", rasterized=True,
)
ax_main.set_xticks([]); ax_main.set_yticks([])
ax_main.set_xlabel("Genome (chr1 → chrY)", fontsize=10)
ax_main.set_ylabel(f"Spots (n = {n_spots:,})", fontsize=10)

# Chromosome boundaries
for v in chr_vals:
    if v == 0: continue
    ax_main.axvline(v, color="black", linewidth=0.4, alpha=0.6)

# Chromosome labels track
ax_chr.set_xlim(0, n_windows)
ax_chr.set_ylim(0, 1)
boundaries = list(chr_vals) + [n_windows]
for i, chrom in enumerate(chr_keys):
    a, b = boundaries[i], boundaries[i + 1]
    color = "#EEEEEE" if i % 2 == 0 else "#CCCCCC"
    ax_chr.axvspan(a, b, facecolor=color, edgecolor="none")
    if (b - a) >= 25:
        ax_chr.text((a + b) / 2, 0.5, chrom.replace("chr", ""),
                    ha="center", va="center", fontsize=7)
ax_chr.set_xticks([]); ax_chr.set_yticks([])
for s in ax_chr.spines.values(): s.set_visible(False)

# Sample dividers on main heatmap
sample_runs = (m_ord["sample_id"].values[1:] != m_ord["sample_id"].values[:-1]).nonzero()[0]
for r in sample_runs:
    ax_main.axhline(r + 0.5, color="white", linewidth=0.2, alpha=0.4)

# Response dividers
resp_runs = (m_ord["response"].values[1:] != m_ord["response"].values[:-1]).nonzero()[0]
for r in resp_runs:
    ax_main.axhline(r + 0.5, color="black", linewidth=0.6, alpha=0.7)

# Colorbar
fig.colorbar(im, cax=ax_cbar, label="Inferred CNV (clipped ±0.3)")

# ─────────────────────────────────────────────────────────────────────────────
# External legend
# ─────────────────────────────────────────────────────────────────────────────
legend_handles = []
def add(title, palette):
    legend_handles.append(Patch(facecolor="none", edgecolor="none", label=f"$\\bf{{{title}}}$"))
    for k, v in palette.items():
        legend_handles.append(Patch(facecolor=v, edgecolor="black", linewidth=0.3, label=str(k)))

add("Response",  resp_pal)
add("Tumor",     tum_pal)
add("Dataset",   ds_pal)
add("Cell type", ct_pal)

fig.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(0.985, 0.5),
    frameon=False, fontsize=7, ncol=1,
    handlelength=1.0, labelspacing=0.4,
)

fig.suptitle(
    f"Universal spatial CNV heatmap — {n_spots:,} spots × {n_windows} windows × "
    f"{m_ord['sample_id'].nunique()} samples",
    fontsize=12, y=0.995,
)

png = FIG_DIR / "spatial_cnv_universal_heatmap.png"
pdf = FIG_DIR / "spatial_cnv_universal_heatmap.pdf"
fig.savefig(png, bbox_inches="tight", dpi=200)
fig.savefig(pdf, bbox_inches="tight")
plt.close(fig)
print(f"  → {png.name}")
print(f"  → {pdf.name}")
print("Done.")
