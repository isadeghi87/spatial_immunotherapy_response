#!/usr/bin/env python3
"""
05z_regenerate_spatial_plots.py — Regenerate per-spot spatial maps for all
GSE238264 Visium sections (4R top row, 3NR bottom row).

Reproduces the style of SPATA2 plotSurface:
  - rocket_r colormap (cream → salmon → magenta → dark-purple → near-black)
  - hexagonal-grid spots, white background, no axes
  - shared global vmin/vmax per feature (p2–p98 across all 7 sections)
  - single vertical colorbar on the right

Features plotted:
  Signatures  : Angiogenesis, IFN_gamma, Immune_checkpoint,
                T_cell_exhaustion, Tfh_like, TGF_beta
  L3 compo    : CD8_T, Plasma, cDC2, Macrophage  (cell2location q05)
  NMF programs: NMF4, NMF5, NMF10                (top-50 gene mean)

Outputs → out/05_multilayer/figures/05z_expression_spatial/
"""

import warnings
warnings.filterwarnings("ignore")

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.sparse import issparse
import seaborn as sns

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))

DATA    = PROJECT / "data" / "spatial_data" / "GSE238264_HCC_IO_Visium"
OUT     = PROJECT / "out" / "05_multilayer" / "figures" / "05z_expression_spatial"
TAB     = PROJECT / "out" / "05_multilayer" / "tables"
C2L_TAB = PROJECT / "out" / "02b_spatial_biomarkers" / "tables"
OUT.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 0

# ── Layout ─────────────────────────────────────────────────────────────
SAMPLES_R  = ["HCC1R", "HCC2R", "HCC3R", "HCC4R"]
SAMPLES_NR = ["HCC5NR", "HCC6NR", "HCC7NR"]
RESPONSE_LABEL = {s: "R" for s in SAMPLES_R}
RESPONSE_LABEL.update({s: "N" for s in SAMPLES_NR})
ALL_SAMPLES = SAMPLES_R + SAMPLES_NR

# Per-feature colormaps matching SPATA2 plotSurface viridis palette options:
#   A = inferno  (dark → bright orange/yellow)  — dense stromal / myeloid
#   B = magma    (dark → cream)                  — T-cell / immune infiltrate
#   C = plasma   (dark → bright yellow/magenta)  — DC / vascular  
#   D = viridis  (dark → bright green/yellow)    — general immune activation
#
# High values map to the bright end (matches R viridis default, NOT reversed).
FEATURE_CMAP: dict[str, str] = {
    # IO Signatures
    "sig_Angiogenesis":       "plasma",     # vascular / VEGF
    "sig_IFN_gamma":          "viridis",    # immune activation (D)
    "sig_Immune_checkpoint":  "magma",      # checkpoint (B)
    "sig_T_cell_exhaustion":  "inferno",    # exhaustion (A)
    "sig_Tfh_like":           "viridis",    # Tfh / TLS (D)
    "sig_TGF_beta":           "plasma",     # stromal / TGF (C)
    # L3 Cell-type composition
    "L3_CD8_T":               "magma",      # cytotoxic T (B)
    "L3_Plasma":              "viridis",    # B / Plasma (D)
    "L3_cDC2":                "plasma",     # myeloid DC (C)
    "L3_Macrophage":          "inferno",    # macrophage / myeloid (A)
    # NMF programs
    "NMF4":                   "magma",      # T-cell program (B)
    "NMF5":                   "viridis",    # cytotoxic program (D)
    "NMF10":                  "inferno",    # myeloid/stromal NR-associated (A)
}
DEFAULT_CMAP = "magma"  # fallback

# ── IO signature gene lists ────────────────────────────────────────────
IO_SIGS = {
    "sig_Angiogenesis":       ["VEGFA", "PECAM1", "KDR", "FLT1", "ANGPT2", "FGF2"],
    "sig_IFN_gamma":          ["IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "IDO1"],
    "sig_Immune_checkpoint":  ["CD274", "PDCD1LG2", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "VSIR"],
    "sig_T_cell_exhaustion":  ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "ENTPD1"],
    "sig_Tfh_like":           ["CXCL13", "CXCR5", "CR2", "CD79A", "MS4A1", "CD19", "IGHG1"],
    "sig_TGF_beta":           ["TGFB1", "TGFB2", "TGFB3", "ACTA2", "COL1A1", "SERPINE1"],
}

# ── Friendly display names ─────────────────────────────────────────────
DISPLAY = {
    "sig_Angiogenesis":      "Angiogenesis",
    "sig_IFN_gamma":         "IFN-γ",
    "sig_Immune_checkpoint": "Immune Checkpoint",
    "sig_T_cell_exhaustion": "T Cell Exhaustion",
    "sig_Tfh_like":          "Tfh-like / TLS",
    "sig_TGF_beta":          "TGF-β",
    "L3_CD8_T":              "L3: CD8+ T cells",
    "L3_Plasma":             "L3: Plasma/B cells",
    "L3_cDC2":               "L3: cDC2",
    "L3_Macrophage":         "L3: Macrophage",
    "NMF4":                  "NMF4 – T cell / Immune",
    "NMF5":                  "NMF5 – T cell / Cytotoxic",
    "NMF10":                 "NMF10 – Myeloid/Stromal (NR-associated)",
}

# ── Load cell2location per-spot abundances ─────────────────────────────
print("Loading cell2location abundances...")
c2l = pd.read_csv(C2L_TAB / "gse238264_cell2location_abundances.csv", index_col=0)
# c2l spot IDs have a section-specific trailing digit (e.g. "-0" for HCC1R,
# "-6" for HCC7NR). h5ad obs_names after renaming are "HCC1R_AAAC...-1".
# Strip the last "-\d+" suffix (anything after the final barcode "-1").
c2l.index = c2l.index.str.replace(r"-\d+$", "", regex=True)
# spot IDs like "HCC1R_AAACAAGTATCTCCCA-1"
c2l["section"] = c2l.index.str.split("_").str[0]

# ── Load NMF gene loadings → build top-50 gene lists ──────────────────
print("Loading NMF gene loadings...")
gl = pd.read_csv(TAB / "nmf_programs_gene_loadings.csv", index_col=0)
NMF_GENES = {}
for prog in ["NMF4", "NMF5", "NMF10"]:
    top50 = gl[prog].nlargest(50).index.tolist()
    NMF_GENES[prog] = top50

# ── Load and score all Visium sections ────────────────────────────────
print("Loading and scoring Visium sections...")
sections: dict[str, sc.AnnData] = {}
for name in ALL_SAMPLES:
    path = DATA / f"{name}.h5ad"
    if not path.exists():
        print(f"  [SKIP] {name}: not found")
        continue
    ad = sc.read_h5ad(str(path))
    ad.var_names_make_unique()
    ad.obs_names = [f"{name}_{b}" for b in ad.obs_names]
    ad.obsm["spatial"] = ad.obs[["pxl_col", "pxl_row"]].values.astype(float)
    ad.obs["section"] = name

    sc.pp.filter_genes(ad, min_cells=3)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    avail = set(ad.var_names)
    # Score IO signatures
    for sig, genes in IO_SIGS.items():
        present = [g for g in genes if g in avail]
        if len(present) >= 2:
            sc.tl.score_genes(ad, present, score_name=sig, use_raw=False)
        else:
            ad.obs[sig] = np.nan

    # Score NMF programs: mean of top-50 genes
    for prog, genes in NMF_GENES.items():
        g_ok = [g for g in genes if g in avail]
        if g_ok:
            X = ad[:, g_ok].X
            if issparse(X):
                X = X.toarray()
            ad.obs[prog] = np.clip(X.mean(axis=1), 0, None)
        else:
            ad.obs[prog] = np.nan

    # Attach cell2location L3 abundances
    sec_c2l = c2l[c2l["section"] == name].copy()
    for ct, col in [("CD8_T", "q05cell_abundance_w_sf_CD8_T"),
                    ("Plasma", "q05cell_abundance_w_sf_Plasma"),
                    ("cDC2",  "q05cell_abundance_w_sf_cDC2"),
                    ("Macrophage", "q05cell_abundance_w_sf_Macrophage")]:
        if col in sec_c2l.columns:
            mapped = sec_c2l[col].reindex(ad.obs_names)
            ad.obs[f"L3_{ct}"] = mapped.values
        else:
            ad.obs[f"L3_{ct}"] = np.nan

    sections[name] = ad
    print(f"  {name}: {ad.n_obs} spots")

# ── Spot-size estimation ───────────────────────────────────────────────
# Choose spot size so adjacent hexagonal grid points just-touch.
# Visium standard pitch is ~100 µm → pixel pitch varies per section.
# Estimate from median nearest-neighbour distance and convert to figure units.
def estimate_spot_size(ad: sc.AnnData, fig_w_in: float = 4.0, dpi: int = 300) -> float:
    """Return matplotlib scatter s so spots tile the Visium grid."""
    coords = ad.obsm["spatial"]
    if len(coords) < 2:
        return 14
    # sample distances
    sample_idx = np.random.choice(len(coords), min(200, len(coords)), replace=False)
    pts = coords[sample_idx]
    dists = []
    for i in range(len(pts)):
        d = np.sqrt(((pts - pts[i]) ** 2).sum(axis=1))
        d[i] = np.inf
        dists.append(d.min())
    pitch_px = float(np.median(dists))
    ax_range = coords[:, 0].ptp()
    ax_w_data = ax_range if ax_range > 0 else 1
    pts_per_data_unit = (fig_w_in * dpi) / ax_w_data
    spot_diameter_pts = pitch_px * pts_per_data_unit / dpi * 72  # convert to points
    return max(5.0, spot_diameter_pts ** 2)


# ── Core plot function ─────────────────────────────────────────────────
PANEL_W = 3.8   # inches per panel column
PANEL_H = 3.6   # inches per panel row
CBAR_W  = 0.55  # colorbar width strip

def make_spatial_grid(feature: str, vmin=None, vmax=None) -> None:
    """Render a 2-row × 4-col grid for feature across all 7 sections."""
    cmap = FEATURE_CMAP.get(feature, DEFAULT_CMAP)
    label = DISPLAY.get(feature, feature)

    # Collect all values for global scale
    all_vals = []
    for name in ALL_SAMPLES:
        ad = sections.get(name)
        if ad is None or feature not in ad.obs.columns:
            continue
        v = ad.obs[feature].dropna().values
        all_vals.extend(v.tolist())

    if not all_vals:
        print(f"  [SKIP] {feature}: no data")
        return

    all_arr = np.array(all_vals, dtype=float)
    if vmin is None:
        vmin = np.nanpercentile(all_arr, 2)
    if vmax is None:
        vmax = np.nanpercentile(all_arr, 98)
    if vmax <= vmin:
        vmax = all_arr.max(); vmin = all_arr.min()

    fig_w = PANEL_W * 4 + CBAR_W + 0.3
    fig_h = PANEL_H * 2
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

    # Leave right margin for colorbar
    left_w = PANEL_W * 4
    # [left, bottom, width, height] in figure fraction
    ax_gs = fig.add_gridspec(
        2, 4,
        left=0.01,
        right=left_w / fig_w - 0.01,
        top=0.90,
        bottom=0.03,
        wspace=0.04,
        hspace=0.10,
    )

    rows = [SAMPLES_R, SAMPLES_NR]
    for row_i, row_samples in enumerate(rows):
        n_cols = len(row_samples)
        for col_j in range(4):
            ax = fig.add_subplot(ax_gs[row_i, col_j])
            ax.set_facecolor("white")
            for sp in ax.spines.values():
                sp.set_visible(False)
            ax.set_xticks([]); ax.set_yticks([])

            if col_j >= n_cols:
                ax.axis("off")
                continue

            name = row_samples[col_j]
            ad = sections.get(name)
            if ad is None:
                ax.axis("off")
                continue

            coords = ad.obsm["spatial"]
            vals   = ad.obs[feature].values.astype(float)

            # Estimate spot size once per section
            ax_x_range = float(coords[:, 0].max() - coords[:, 0].min())
            ax_y_range = float(coords[:, 1].max() - coords[:, 1].min())

            # Spot size: map data-unit pitch → figure points, then to s (pts²).
            # PANEL_W inches × 72 pt/in = axis width in pts (approximate).
            ax_w_pts = PANEL_W * 72
            pts_per_data = ax_w_pts / (ax_x_range if ax_x_range > 0 else 1)

            # Estimate spot pitch from nearest-neighbour distances
            sample = np.random.choice(len(coords), min(150, len(coords)), replace=False)
            pts_sample = coords[sample]
            dists = np.array([
                sorted(np.sqrt(((pts_sample - p) ** 2).sum(axis=1)))[1]
                for p in pts_sample
            ])
            pitch_data = float(np.median(dists))
            # 0.72 gives a slight gap so adjacent spots do not fully occlude
            spot_diam_pts = pitch_data * pts_per_data * 0.72
            s = max(5.0, spot_diam_pts ** 2)

            sc_ = ax.scatter(
                coords[:, 0], coords[:, 1],
                c=vals, cmap=cmap,
                vmin=vmin, vmax=vmax,
                s=s,
                linewidths=0,
                rasterized=True,
            )
            ax.set_aspect("equal")
            ax.invert_yaxis()

            resp_tag = RESPONSE_LABEL[name]
            ax.set_title(f"{name}\n({resp_tag})", fontsize=9, fontweight="bold",
                         pad=3, color="#222222")

    # Shared vertical colorbar
    cbar_left = left_w / fig_w + 0.005
    cbar_ax = fig.add_axes([cbar_left, 0.13, CBAR_W / fig_w, 0.70])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cbar_ax)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_linewidth(0.5)
    cb.set_label(label, fontsize=9, rotation=270, labelpad=14, va="bottom")

    fig.suptitle(label, fontsize=13, fontweight="bold", color="#111111", y=0.97)

    for ext in ("png", "pdf"):
        out_path = OUT / f"05z_spatial_{feature}.{ext}"
        fig.savefig(out_path, dpi=280, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved 05z_spatial_{feature}.{{png,pdf}}")


# ── Run all features ───────────────────────────────────────────────────
ALL_FEATURES = list(IO_SIGS.keys()) + ["L3_CD8_T", "L3_Plasma", "L3_cDC2", "L3_Macrophage"] + ["NMF4", "NMF5", "NMF10"]

print("\nGenerating spatial plots...")
np.random.seed(42)
for feat in ALL_FEATURES:
    print(f"  → {feat}")
    make_spatial_grid(feat)

print(f"\nDone. Output: {OUT}")
