#!/usr/bin/env python3
"""
02k_cappuyns_spatial_integration.py — Map Cappuyns signatures onto Visium spatial data.

Scores the 7 significant Cappuyns signatures on GSE238264 HCC IO Visium spots
and generates spatial heatmaps overlaid on tissue sections.

Requires:
  - data/references/cappuyns2025_mmc2_signatures.xlsx
  - data/spatial_data/GSE238264_HCC_IO_Visium/*.h5ad
  - out/02i_cappuyns_signatures/tables/cappuyns_signatures_parsed.json
"""
import warnings
warnings.filterwarnings("ignore")

import json
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import sys; sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
from config.publication_style import apply_style, PAL, CMAP, label_clean
apply_style()
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr

PROJECT = Path(__file__).resolve().parents[2]
SPATIAL = PROJECT / "data" / "spatial_data" / "GSE238264_HCC_IO_Visium"
OUT     = PROJECT / "out" / "02k_cappuyns_spatial"
FIG     = OUT / "figures"
TAB     = OUT / "tables"
for d in [FIG, TAB]:
    d.mkdir(parents=True, exist_ok=True)

# 7 significant signatures from Cappuyns analysis (p_adj < 0.05)
SIGNIFICANT_SIGS = [
    "Macro CXCL10", "CD8 Tex", "CD4 CXCL13",
    "Mono CD14", "Kupffer cells", "Macro SPP1", "Fibroblasts",
]


def load_signatures():
    """Load parsed Cappuyns signatures from JSON."""
    json_path = PROJECT / "out" / "02i_cappuyns_signatures" / "tables" / "cappuyns_signatures_parsed.json"
    if not json_path.exists():
        print(f"  ERROR: {json_path} not found — run 02i first")
        return None
    with open(json_path) as f:
        sigs = json.load(f)
    return sigs


def score_spatial_section(adata, signatures, sig_names):
    """Score selected Cappuyns signatures on a spatial section."""
    available = set(adata.var_names)
    scored_cols = []

    for sig_name in sig_names:
        if sig_name not in signatures:
            continue
        up_genes = [g for g in signatures[sig_name].get("up", []) if g in available]
        down_genes = [g for g in signatures[sig_name].get("down", []) if g in available]

        col = f"cap_{sig_name.replace(' ', '_')}"

        if len(up_genes) >= 3:
            sc.tl.score_genes(adata, up_genes, score_name=f"{col}_up")
            if len(down_genes) >= 3:
                sc.tl.score_genes(adata, down_genes, score_name=f"{col}_down")
                adata.obs[col] = adata.obs[f"{col}_up"] - adata.obs[f"{col}_down"]
            else:
                adata.obs[col] = adata.obs[f"{col}_up"]
            scored_cols.append(col)

    return scored_cols


def main():
    print("=" * 70)
    print("Cappuyns Signatures — Spatial Integration (GSE238264 Visium)")
    print("=" * 70)

    signatures = load_signatures()
    if signatures is None:
        return

    h5_files = sorted(SPATIAL.glob("*.h5ad"))
    if len(h5_files) == 0:
        print("  ERROR: No Visium h5ad files found")
        return

    print(f"  Found {len(h5_files)} Visium sections")

    # ── Score all sections ─────────────────────────────────────────
    all_section_scores = []
    section_data = {}

    for h5_path in h5_files:
        sample_name = h5_path.stem
        response = "Responder" if "R" in sample_name and "NR" not in sample_name else "Non-Responder"
        print(f"\n  [{sample_name}] ({response}) Loading...")

        adata = sc.read_h5ad(h5_path)
        adata.var_names_make_unique()
        print(f"    {adata.shape[0]} spots × {adata.shape[1]} genes")

        scored_cols = score_spatial_section(adata, signatures, SIGNIFICANT_SIGS)
        print(f"    Scored {len(scored_cols)} signatures")

        # Store for plotting
        section_data[sample_name] = {
            "adata": adata,
            "scored_cols": scored_cols,
            "response": response,
        }

        # Aggregate per-section means
        for col in scored_cols:
            all_section_scores.append({
                "section": sample_name,
                "response": response,
                "signature": col.replace("cap_", "").replace("_", " "),
                "mean_score": adata.obs[col].mean(),
                "median_score": adata.obs[col].median(),
                "std_score": adata.obs[col].std(),
            })

    scores_df = pd.DataFrame(all_section_scores)
    scores_df.to_csv(TAB / "spatial_cappuyns_section_scores.csv", index=False)

    # ═══════════════════════════════════════════════════════════════
    # FIGURE 1: Spatial heatmaps of top 3 signatures per section
    # ═══════════════════════════════════════════════════════════════
    top_sigs = ["cap_CD8_Tex", "cap_Macro_CXCL10", "cap_CD4_CXCL13"]
    # Verify available
    top_sigs = [s for s in top_sigs if any(s in sd["scored_cols"] for sd in section_data.values())]

    if top_sigs:
        n_sections = len(section_data)
        n_sigs = len(top_sigs)
        fig, axes = plt.subplots(n_sigs, n_sections,
                                 figsize=(3.2 * n_sections, 3.5 * n_sigs),
                                 squeeze=False)

        for j, (sample_name, sd) in enumerate(sorted(section_data.items())):
            ad = sd["adata"]
            resp_label = "R" if sd["response"] == "Responder" else "NR"

            for i, sig_col in enumerate(top_sigs):
                ax = axes[i][j]
                if sig_col not in ad.obs.columns:
                    ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                            transform=ax.transAxes)
                    ax.set_title(f"{sample_name}", fontsize=8)
                    ax.axis("off")
                    continue

                # Get spatial coordinates
                if "spatial" in ad.obsm:
                    coords = ad.obsm["spatial"]
                elif "pxl_col" in ad.obs.columns:
                    coords = ad.obs[["pxl_col", "pxl_row"]].values.astype(float)
                else:
                    ax.text(0.5, 0.5, "No spatial", ha="center", va="center",
                            transform=ax.transAxes)
                    ax.axis("off")
                    continue

                vals = ad.obs[sig_col].values
                vmin, vmax = np.nanpercentile(vals, [2, 98])

                sc_plot = ax.scatter(coords[:, 0], coords[:, 1],
                                     c=vals, cmap="RdBu_r", vmin=vmin, vmax=vmax,
                                     s=4, alpha=0.8, edgecolor="none", rasterized=True)
                ax.set_aspect("equal")
                ax.invert_yaxis()
                ax.set_xticks([]); ax.set_yticks([])

                if i == 0:
                    ax.set_title(f"{sample_name} ({resp_label})", fontsize=9, fontweight="bold")
                if j == 0:
                    sig_label = sig_col.replace("cap_", "").replace("_", " ")
                    ax.set_ylabel(sig_label, fontsize=9, fontweight="bold")

        plt.suptitle("Cappuyns Signature Spatial Maps — GSE238264 HCC IO Visium",
                     fontweight="bold", fontsize=12, y=1.01)
        plt.tight_layout()
        plt.savefig(FIG / "cappuyns_spatial_heatmaps.png", dpi=200, bbox_inches="tight")
        plt.close()
        print(f"\n  Saved: {FIG / 'cappuyns_spatial_heatmaps.png'}")

    # ═══════════════════════════════════════════════════════════════
    # FIGURE 2: Section-level signature comparison R vs NR
    # ═══════════════════════════════════════════════════════════════
    if len(scores_df) > 0:
        unique_sigs = scores_df["signature"].unique()
        n = len(unique_sigs)
        ncols = min(4, n)
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows), squeeze=False)
        axes_flat = axes.flatten()

        palette = {"Responder": "#0077BB", "Non-Responder": "#CC3311"}
        for i, sig in enumerate(unique_sigs):
            ax = axes_flat[i]
            sub = scores_df[scores_df["signature"] == sig]

            sns.stripplot(data=sub, x="response", y="mean_score", palette=palette,
                          size=10, alpha=0.8, ax=ax, jitter=0.1,
                          edgecolor="gray", linewidth=0.5)

            # Connect paired data isn't possible here, so just show means
            for resp in ["Responder", "Non-Responder"]:
                vals = sub.loc[sub["response"] == resp, "mean_score"]
                mean_val = vals.mean()
                x_pos = 0 if resp == "Responder" else 1
                ax.hlines(mean_val, x_pos - 0.2, x_pos + 0.2,
                          color="black", linewidth=2, zorder=5)

            r_vals = sub.loc[sub["response"] == "Responder", "mean_score"]
            nr_vals = sub.loc[sub["response"] == "Non-Responder", "mean_score"]
            if len(r_vals) >= 2 and len(nr_vals) >= 2:
                _, p = mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
                ax.set_title(f"{sig}\np = {p:.3g}", fontsize=9, fontweight="bold")
            else:
                ax.set_title(sig, fontsize=9, fontweight="bold")

            ax.set_xlabel("")
            ax.set_ylabel("Mean Spot Score" if i % ncols == 0 else "")

        for j in range(i + 1, len(axes_flat)):
            axes_flat[j].set_visible(False)

        plt.suptitle("Cappuyns Signatures in Visium Sections: R vs NR",
                     fontweight="bold", fontsize=12, y=1.01)
        plt.tight_layout()
        plt.savefig(FIG / "cappuyns_spatial_section_comparison.png", dpi=200, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {FIG / 'cappuyns_spatial_section_comparison.png'}")

    print(f"\n{'=' * 70}")
    print("Done.")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
