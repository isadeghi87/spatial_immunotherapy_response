#!/usr/bin/env python3
"""04h – High-contrast biomarker topology maps from real atlas data.

Computes gene-signature scores on the integrated IO atlas and visualises
them on UMAP coordinates:
  A. Hypoxia gradient
  B. CD8+ T-cell exhaustion hubs
  C. Macrophage M1 / M2 polarization
  D. Predicted IO response vulnerability
"""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
ATLAS   = PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad"
OUT     = PROJECT / "out" / "04_complementary" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# ── Curated gene signatures ─────────────────────────────────────────
SIGNATURES = {
    "hypoxia": ["HIF1A", "VEGFA", "CA9", "SLC2A1", "LDHA", "PGK1",
                "ENO1", "ALDOA", "PKM", "BNIP3", "NDRG1", "ADM"],
    "cd8_exhaustion": ["PDCD1", "HAVCR2", "LAG3", "TOX", "TIGIT",
                       "ENTPD1", "CTLA4", "LAYN", "CXCL13"],
    "m1_macrophage": ["TNF", "NOS2", "IL1B", "CXCL10", "CXCL9",
                      "IL12A", "CCL5", "IRF1", "STAT1"],
    "m2_macrophage": ["CD163", "MRC1", "MSR1", "TGFB1", "IL10",
                      "CCL18", "FOLR2", "MAF", "STAB1"],
    "io_response": ["CD8A", "CD8B", "GZMB", "PRF1", "IFNG", "CXCL13",
                    "CXCL9", "CXCL10", "IDO1", "STAT1", "GBP1", "PDCD1"],
}


def _score(adata, name, genes):
    present = [g for g in genes if g in adata.var_names]
    if len(present) < 3:
        print(f"  Warning: {name}: only {len(present)} genes found – skipping")
        return False
    sc.tl.score_genes(adata, gene_list=present, score_name=name, use_raw=False)
    print(f"  {name}: scored with {len(present)}/{len(genes)} genes")
    return True


def plot_spatial_maps(adata_path=None):
    if adata_path is None:
        adata_path = str(ATLAS)

    print("[04h] Loading atlas …")
    adata = sc.read_h5ad(adata_path)
    print(f"  Shape: {adata.shape}")

    # UMAP coordinates
    if "X_umap" not in adata.obsm:
        rep = "X_pca_harmony" if "X_pca_harmony" in adata.obsm else "X_pca"
        sc.pp.neighbors(adata, use_rep=rep, n_neighbors=15)
        sc.tl.umap(adata)

    # Subsample for visualisation (max 80k cells)
    if adata.n_obs > 80_000:
        idx = np.random.RandomState(42).choice(adata.n_obs, 80_000, replace=False)
        adata = adata[idx].copy()
        print(f"  Subsampled to {adata.n_obs} cells")

    umap = adata.obsm["X_umap"]

    # ── Score signatures ─────────────────────────────────────────────
    print("  Computing gene-signature scores …")
    for name, genes in SIGNATURES.items():
        _score(adata, name, genes)

    # M1 – M2 difference
    if "m1_macrophage" in adata.obs and "m2_macrophage" in adata.obs:
        adata.obs["m1_m2_ratio"] = adata.obs["m1_macrophage"] - adata.obs["m2_macrophage"]

    # ── 4-panel figure ───────────────────────────────────────────────
    rng = np.random.RandomState(0)
    order = rng.permutation(adata.n_obs)          # avoid overplot bias

    panels = [
        ("hypoxia",        "A. Hypoxia Gradient across TME",
         "Hypoxia Signature Score",  "magma"),
        ("cd8_exhaustion", "B. CD8+ T-Cell Exhaustion Hubs",
         "CD8 Exhaustion Score",     "inferno"),
        ("m1_m2_ratio",    "C. Macrophage Polarization (M1 vs M2)",
         "M1 - M2 Score",           "RdBu_r"),
        ("io_response",    "D. Predicted IO Response Vulnerability",
         "IO Response Score",       "viridis"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(18, 16))
    fig.suptitle("Biomarker Topology across the Tumour Micro-Environment",
                 fontsize=20, fontweight="bold", y=0.98)

    for ax, (key, title, cbar_label, cmap) in zip(axes.ravel(), panels):
        if key not in adata.obs:
            ax.set_title(title, fontweight="bold", fontsize=14)
            ax.text(0.5, 0.5, "score unavailable", ha="center", va="center",
                    transform=ax.transAxes, fontsize=12, color="grey")
            continue

        vals = adata.obs[key].values
        kw = dict(s=4, alpha=0.7, edgecolor="none", rasterized=True)

        if cmap == "RdBu_r":
            vmax = np.nanpercentile(np.abs(vals), 95)
            norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
            sc_pl = ax.scatter(umap[order, 0], umap[order, 1],
                               c=vals[order], cmap=cmap, norm=norm, **kw)
        else:
            sc_pl = ax.scatter(umap[order, 0], umap[order, 1],
                               c=vals[order], cmap=cmap, **kw)

        cb = plt.colorbar(sc_pl, ax=ax, shrink=0.7, pad=0.02)
        cb.set_label(cbar_label, fontweight="bold")
        cb.outline.set_visible(False)
        ax.set_title(title, fontweight="bold", fontsize=14)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    combined = OUT / "04h_spatial_high_contrast_maps.png"
    plt.savefig(combined, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved combined -> {combined}")

    # ── Individual panels ────────────────────────────────────────────
    names_map = {
        "hypoxia":        "spatial_1_hypoxia_gradient.png",
        "cd8_exhaustion": "spatial_2_cd8_exhaustion_hubs.png",
        "m1_m2_ratio":    "spatial_3_macrophage_polarization.png",
        "io_response":    "spatial_4_io_response_prediction.png",
    }
    for key, fname in names_map.items():
        if key not in adata.obs:
            continue
        vals = adata.obs[key].values
        fig_i, ax_i = plt.subplots(figsize=(8, 6), dpi=140)
        kw = dict(s=12, alpha=0.85, edgecolor="none", rasterized=True)
        cmap_i = {"hypoxia": "magma", "cd8_exhaustion": "inferno",
                  "m1_m2_ratio": "RdBu_r", "io_response": "viridis"}[key]
        if cmap_i == "RdBu_r":
            vmax = np.nanpercentile(np.abs(vals), 95)
            norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
            sc_i = ax_i.scatter(umap[order, 0], umap[order, 1],
                                c=vals[order], cmap=cmap_i, norm=norm, **kw)
        else:
            sc_i = ax_i.scatter(umap[order, 0], umap[order, 1],
                                c=vals[order], cmap=cmap_i, **kw)

        cb_i = plt.colorbar(sc_i, ax=ax_i, shrink=0.8)
        panel_info = [p for p in panels if p[0] == key][0]
        cb_i.set_label(panel_info[2], fontweight="bold")
        cb_i.outline.set_visible(False)
        ax_i.set_title(panel_info[1], fontsize=14, fontweight="bold")
        ax_i.set_xticks([]); ax_i.set_yticks([])
        for sp in ax_i.spines.values():
            sp.set_visible(False)
        plt.tight_layout()
        plt.savefig(OUT / fname, bbox_inches="tight")
        plt.close()
        print(f"  Saved -> {OUT / fname}")

    # ── Tissue-stratified spatial surface plots ──────────────────────
    # Resolve tissue column
    tissue_col = None
    for col in ["tissue", "Tissue", "site"]:
        if col in adata.obs.columns:
            vals = adata.obs[col].astype(str)
            if (vals != "nan").sum() > 1000:
                tissue_col = col
                break

    if tissue_col is not None:
        tissue_vals = adata.obs[tissue_col].astype(str)
        keep_tissues = ["Tumor", "Normal"]
        avail = [t for t in keep_tissues if (tissue_vals == t).sum() > 200]
        print(f"  Tissue-stratified panels for: {avail}")

        sig_keys = [k for k in ["hypoxia", "cd8_exhaustion", "m1_m2_ratio", "io_response"]
                    if k in adata.obs]
        if avail and sig_keys:
            n_sigs = len(sig_keys)
            n_tissues = len(avail) + 1  # +1 for "All"
            fig_t, axes_t = plt.subplots(n_sigs, n_tissues,
                                          figsize=(6 * n_tissues, 5 * n_sigs))
            if n_sigs == 1:
                axes_t = axes_t[np.newaxis, :]
            if n_tissues == 1:
                axes_t = axes_t[:, np.newaxis]

            sig_info = {
                "hypoxia":        ("Hypoxia",         "magma"),
                "cd8_exhaustion": ("CD8 Exhaustion",   "inferno"),
                "m1_m2_ratio":    ("M1 - M2",         "RdBu_r"),
                "io_response":    ("IO Response",      "viridis"),
            }
            col_labels = ["All Cells"] + avail

            for row, key in enumerate(sig_keys):
                label, cmap_s = sig_info[key]
                vals_s = adata.obs[key].values
                vmin_s, vmax_s = np.nanpercentile(vals_s, [2, 98])

                for col_idx, tissue_label in enumerate(col_labels):
                    ax_t = axes_t[row, col_idx]

                    if tissue_label == "All Cells":
                        mask = np.ones(adata.n_obs, dtype=bool)
                    else:
                        mask = tissue_vals.values == tissue_label

                    idx_m = np.where(mask)[0]
                    rng_t = np.random.RandomState(row + col_idx)
                    idx_m = rng_t.permutation(idx_m)

                    kw_t = dict(s=6, alpha=0.75, edgecolor="none", rasterized=True)

                    # Background: grey for non-selected cells
                    if tissue_label != "All Cells":
                        bg = np.where(~mask)[0]
                        ax_t.scatter(umap[bg, 0], umap[bg, 1],
                                     c="#e0e0e0", s=1, alpha=0.15, rasterized=True)

                    if cmap_s == "RdBu_r":
                        vabs = max(abs(vmin_s), abs(vmax_s))
                        norm_t = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
                        sc_t = ax_t.scatter(umap[idx_m, 0], umap[idx_m, 1],
                                            c=vals_s[idx_m], cmap=cmap_s, norm=norm_t, **kw_t)
                    else:
                        sc_t = ax_t.scatter(umap[idx_m, 0], umap[idx_m, 1],
                                            c=vals_s[idx_m], cmap=cmap_s,
                                            vmin=vmin_s, vmax=vmax_s, **kw_t)

                    if col_idx == n_tissues - 1:
                        cb_t = plt.colorbar(sc_t, ax=ax_t, shrink=0.7, pad=0.02)
                        cb_t.outline.set_visible(False)

                    if row == 0:
                        ax_t.set_title(tissue_label, fontweight="bold", fontsize=13)
                    if col_idx == 0:
                        ax_t.set_ylabel(label, fontweight="bold", fontsize=12)

                    n_cells = mask.sum() if tissue_label != "All Cells" else adata.n_obs
                    med = np.nanmedian(vals_s[mask] if tissue_label != "All Cells" else vals_s)
                    ax_t.text(0.02, 0.02, f"n={n_cells:,}\nmed={med:.3f}",
                              transform=ax_t.transAxes, fontsize=8, va="bottom",
                              bbox=dict(fc="white", ec="none", alpha=0.7))

                    ax_t.set_xticks([]); ax_t.set_yticks([])
                    for sp in ax_t.spines.values():
                        sp.set_visible(False)

            fig_t.suptitle("Tissue-Stratified Biomarker Maps (Tumor vs Normal)",
                           fontsize=18, fontweight="bold", y=1.01)
            plt.tight_layout()
            out_tissue = OUT / "04h_tissue_stratified_spatial.png"
            plt.savefig(out_tissue, dpi=250, bbox_inches="tight")
            plt.close()
            print(f"  Saved tissue-stratified -> {out_tissue}")
    else:
        print("  No tissue column found – skipping tissue-stratified plots.")

    print("[04h] Done.")


if __name__ == "__main__":
    plot_spatial_maps()
