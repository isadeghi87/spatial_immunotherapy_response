#!/usr/bin/env python3
"""04d – Cell-cell communication analysis using real atlas data.

Computes ligand-receptor interaction scores between cell types using
a curated IO-relevant LR pair list and the real expression data.
"""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu

PROJECT = Path(__file__).resolve().parents[2]
ATLAS   = PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad"
OUT     = PROJECT / "out" / "04_complementary" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def _resolve_celltype_col(adata):
    for col in ("cell_type", "celltype_sub", "Cell.type", "celltype_global", "celltype"):
        if col in adata.obs.columns:
            return col
    raise KeyError("No cell-type column found")


# Curated IO-relevant ligand-receptor pairs
LR_PAIRS = [
    ("CD274", "PDCD1"),       # PD-L1 → PD-1
    ("CD80", "CD28"),         # CD80 → CD28
    ("CD80", "CTLA4"),        # CD80 → CTLA-4
    ("CD86", "CTLA4"),        # CD86 → CTLA-4
    ("CXCL10", "CXCR3"),     # CXCL10 → CXCR3
    ("CCL19", "CCR7"),        # CCL19 → CCR7
    ("VEGFA", "KDR"),         # VEGF → VEGFR2
    ("TGFB1", "TGFBR1"),     # TGF-β → TGFBR1
    ("TNF", "TNFRSF1A"),     # TNF-α → TNFR1
    ("FASLG", "FAS"),         # FasL → Fas
]


def plot_advanced_ccc(adata_path=None):
    if adata_path is None:
        adata_path = str(ATLAS)

    print("[04d] Loading atlas …")
    adata = sc.read_h5ad(adata_path)
    ct_col = _resolve_celltype_col(adata)
    print(f"  Shape: {adata.shape}, cell-type col: {ct_col}")

    # Filter to LR genes that exist in data
    all_genes = set(adata.var_names)
    valid_pairs = [(l, r) for l, r in LR_PAIRS if l in all_genes and r in all_genes]
    print(f"  Valid LR pairs: {len(valid_pairs)} / {len(LR_PAIRS)}")

    if len(valid_pairs) == 0:
        raise ValueError("No valid LR pairs found in the dataset")

    # Pick top cell types by abundance
    top_cts = adata.obs[ct_col].value_counts().nlargest(8).index.tolist()

    # Compute mean expression of ligand in sender and receptor in receiver
    # for each cell-type pair
    print("  Computing LR scores …")
    ct_means = {}
    for ct in top_cts:
        mask = adata.obs[ct_col] == ct
        subset = adata[mask]
        genes_needed = list(set([g for pair in valid_pairs for g in pair]))
        genes_present = [g for g in genes_needed if g in adata.var_names]
        if hasattr(subset[:, genes_present].X, "toarray"):
            expr = pd.DataFrame(subset[:, genes_present].X.toarray(),
                                columns=genes_present, index=subset.obs_names)
        else:
            expr = pd.DataFrame(np.array(subset[:, genes_present].X),
                                columns=genes_present, index=subset.obs_names)
        ct_means[ct] = expr.mean()

    # Build interaction score matrix: sender_ct × (LR pair)
    pair_labels = [f"{l}→{r}" for l, r in valid_pairs]
    # For each sender-receiver combo, compute product of mean ligand in sender × mean receptor in receiver
    # Use top 5 cell-type pairs for readability
    sender_receiver_pairs = []
    for s in top_cts[:5]:
        for r_ct in top_cts[:5]:
            if s != r_ct:
                sender_receiver_pairs.append((s, r_ct))
    sender_receiver_pairs = sender_receiver_pairs[:10]  # limit for plot

    scores = pd.DataFrame(0.0, index=[f"{s} → {r}" for s, r in sender_receiver_pairs],
                           columns=pair_labels)
    pvals_mat = pd.DataFrame(1.0, index=scores.index, columns=scores.columns)

    for (s_ct, r_ct) in sender_receiver_pairs:
        row_label = f"{s_ct} → {r_ct}"
        for (lig, rec) in valid_pairs:
            col_label = f"{lig}→{rec}"
            lig_expr = ct_means[s_ct].get(lig, 0)
            rec_expr = ct_means[r_ct].get(rec, 0)
            scores.loc[row_label, col_label] = lig_expr * rec_expr

    # Normalise scores for visualisation
    scores_norm = scores / (scores.max().max() + 1e-10)

    # ── Dotplot ──────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(12, 8))

    x_idx = {c: i for i, c in enumerate(scores_norm.columns)}
    y_idx = {r: i for i, r in enumerate(scores_norm.index)}

    for col in scores_norm.columns:
        for row in scores_norm.index:
            val = scores_norm.loc[row, col]
            raw = scores.loc[row, col]
            ax.scatter(x_idx[col], y_idx[row], s=val * 500 + 10, c=[raw],
                       cmap="RdBu_r", vmin=0, vmax=scores.values.max(), alpha=0.85,
                       edgecolor="black", linewidth=0.5)

    ax.set_xticks(range(len(pair_labels)))
    ax.set_xticklabels(pair_labels, rotation=45, ha="right", fontweight="bold", fontsize=9)
    ax.set_yticks(range(len(scores_norm.index)))
    ax.set_yticklabels(scores_norm.index, fontsize=9)
    ax.set_title("Ligand-Receptor Communication Scores (Real Data)", fontweight="bold", pad=15)

    sm = plt.cm.ScalarMappable(cmap="RdBu_r",
                                norm=plt.Normalize(vmin=0, vmax=scores.values.max()))
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label="Interaction Score", shrink=0.6)

    plt.tight_layout()
    plt.savefig(OUT / "04d_ccc_lr_dotplot.png", dpi=300, bbox_inches="tight")
    plt.savefig(OUT / "04d_ccc_lr_dotplot.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {OUT / '04d_ccc_lr_dotplot.png'}")


if __name__ == "__main__":
    plot_advanced_ccc()
