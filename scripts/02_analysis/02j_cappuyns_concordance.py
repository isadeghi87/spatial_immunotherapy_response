#!/usr/bin/env python3
"""
02j_cappuyns_concordance.py — Concordance between internal immune subtypes
and Cappuyns molecular subsets.

Reads:
  - out/04_complementary/tables/patient_immune_subtypes.csv   (internal Subtype_1/2/3)
  - out/02i_cappuyns_signatures/tables/patient_molecular_subsets.csv (Cappuyns subsets)

Produces:
  - Concordance heatmap (confusion matrix style)
  - Alluvial / Sankey-style summary of patient flow
  - Cohen's kappa and concordance statistics
"""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import sys; sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
from config.publication_style import apply_style, PAL, CMAP, label_clean
apply_style()
from pathlib import Path
from scipy.stats import chi2_contingency

PROJECT = Path(__file__).resolve().parents[2]
OUT = PROJECT / "out" / "02j_concordance"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in [FIG, TAB]:
    d.mkdir(parents=True, exist_ok=True)


def cohens_kappa(ct):
    """Compute Cohen's kappa from a confusion matrix (numpy array)."""
    n = ct.sum()
    p_o = np.diag(ct).sum() / n
    row_sums = ct.sum(axis=1)
    col_sums = ct.sum(axis=0)
    p_e = (row_sums * col_sums).sum() / (n * n)
    if p_e == 1.0:
        return 1.0
    return (p_o - p_e) / (1.0 - p_e)


def main():
    print("=" * 70)
    print("Cappuyns vs Internal Subtype Concordance Analysis")
    print("=" * 70)

    # ── Load internal subtypes ─────────────────────────────────────
    internal_path = PROJECT / "out" / "04_complementary" / "tables" / "patient_immune_subtypes.csv"
    if not internal_path.exists():
        print(f"  ERROR: {internal_path} not found")
        return
    internal = pd.read_csv(internal_path)
    print(f"  Internal subtypes: {len(internal)} patients")

    # ── Load Cappuyns subsets ──────────────────────────────────────
    cappuyns_path = PROJECT / "out" / "02i_cappuyns_signatures" / "tables" / "patient_molecular_subsets.csv"
    if not cappuyns_path.exists():
        print(f"  ERROR: {cappuyns_path} not found")
        return
    cappuyns = pd.read_csv(cappuyns_path, index_col=0)
    print(f"  Cappuyns subsets: {len(cappuyns)} patients")

    # ── Merge ──────────────────────────────────────────────────────
    # internal uses patient_ID column, cappuyns uses index
    internal["patient_ID"] = internal["patient_ID"].astype(float)
    cappuyns_sub = cappuyns[["response", "molecular_subset"]].copy()
    cappuyns_sub.index = cappuyns_sub.index.astype(float)
    cappuyns_sub = cappuyns_sub.rename(columns={"molecular_subset": "cappuyns_subset"})

    merged = internal[["patient_ID", "response", "immune_subtype"]].merge(
        cappuyns_sub[["cappuyns_subset"]], left_on="patient_ID", right_index=True, how="inner"
    )
    print(f"  Merged: {len(merged)} patients")

    if len(merged) < 5:
        print("  [SKIP] Too few common patients")
        return

    merged.to_csv(TAB / "concordance_merged.csv", index=False)

    # ── Relabel for clarity ────────────────────────────────────────
    subtype_labels = {
        "Subtype_1": "Immune-active (S1)",
        "Subtype_2": "Mixed (S2)",
        "Subtype_3": "Immune-cold (S3)",
    }
    merged["internal_label"] = merged["immune_subtype"].map(subtype_labels).fillna(merged["immune_subtype"])

    # ── Crosstab ───────────────────────────────────────────────────
    ct = pd.crosstab(merged["internal_label"], merged["cappuyns_subset"])
    ct.to_csv(TAB / "concordance_crosstab.csv")
    print(f"\n  Cross-tabulation:")
    print(ct.to_string())

    # ── Chi-squared test ───────────────────────────────────────────
    if ct.shape[0] >= 2 and ct.shape[1] >= 2 and ct.values.min() >= 0:
        chi2, p_chi, dof, expected = chi2_contingency(ct.values)
        print(f"\n  Chi-squared = {chi2:.2f}, p = {p_chi:.4f}, dof = {dof}")
    else:
        chi2, p_chi = np.nan, np.nan

    # ── Cohen's kappa (binarised: S1 ↔ Immune-comp, S3 ↔ Resistant) ──
    merged["internal_binary"] = merged["immune_subtype"].map(
        {"Subtype_1": "Immune-competent", "Subtype_3": "Resistant"})
    merged["cappuyns_binary"] = merged["cappuyns_subset"]
    binary = merged.dropna(subset=["internal_binary"])
    if len(binary) >= 5:
        ct_bin = pd.crosstab(binary["internal_binary"], binary["cappuyns_binary"])
        # Reorder to match for kappa
        labels = sorted(set(ct_bin.index) | set(ct_bin.columns))
        ct_bin = ct_bin.reindex(index=labels, columns=labels, fill_value=0)
        kappa = cohens_kappa(ct_bin.values)
        print(f"  Cohen's kappa (binary: S1↔IC, S3↔Resistant): {kappa:.3f}")
    else:
        kappa = np.nan

    # ═══════════════════════════════════════════════════════════════
    # FIGURE 1: Concordance heatmap
    # ═══════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [1.2, 1]})

    # Panel A: Count heatmap
    ax = axes[0]
    sns.heatmap(ct, annot=True, fmt="d", cmap="YlOrRd", linewidths=1, linecolor="white",
                ax=ax, cbar_kws={"label": "n patients", "shrink": 0.7},
                annot_kws={"fontsize": 14, "fontweight": "bold"})
    ax.set_xlabel("Cappuyns Molecular Subset", fontsize=11)
    ax.set_ylabel("Internal Immune Subtype", fontsize=11)
    ax.set_title("A. Cross-Classification Concordance", fontweight="bold", fontsize=12)
    ax.tick_params(axis="x", rotation=25, labelsize=9)
    ax.tick_params(axis="y", rotation=0, labelsize=9)

    # Panel B: Alluvial-style summary
    ax = axes[1]
    ax.axis("off")

    # Draw patient flow lines
    internal_cats = merged["internal_label"].unique()
    cappuyns_cats = merged["cappuyns_subset"].unique()
    internal_cats = sorted(internal_cats)
    cappuyns_cats = sorted(cappuyns_cats)

    left_y = {cat: i for i, cat in enumerate(internal_cats)}
    right_y = {cat: i for i, cat in enumerate(cappuyns_cats)}
    n_left = len(internal_cats)
    n_right = len(cappuyns_cats)

    colors_resp = {"Responder": "#0077BB", "Non-Responder": "#CC3311"}

    for _, row in merged.iterrows():
        y1 = left_y[row["internal_label"]]
        y2 = right_y[row["cappuyns_subset"]]
        c = colors_resp.get(row["response"], "#aaa")
        ax.plot([0, 1], [y1, y2], color=c, alpha=0.5, linewidth=2)

    for cat, y in left_y.items():
        ax.text(-0.05, y, cat, ha="right", va="center", fontsize=9, fontweight="bold")
    for cat, y in right_y.items():
        ax.text(1.05, y, cat, ha="left", va="center", fontsize=9, fontweight="bold")

    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.5, max(n_left, n_right) - 0.5)

    # Legend for response color
    from matplotlib.lines import Line2D
    legend_els = [Line2D([0], [0], color="#0077BB", lw=2, label="Responder"),
                  Line2D([0], [0], color="#CC3311", lw=2, label="Non-Responder")]
    ax.legend(handles=legend_els, loc="lower center", fontsize=9, ncol=2)

    stat_str = f"χ² = {chi2:.2f}, p = {p_chi:.4f}" if not np.isnan(chi2) else ""
    if not np.isnan(kappa):
        stat_str += f"  |  Cohen's κ = {kappa:.3f}"
    ax.set_title(f"B. Patient Flow  ({stat_str})", fontweight="bold", fontsize=11, pad=10)

    plt.tight_layout()
    plt.savefig(FIG / "cappuyns_internal_concordance.png", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {FIG / 'cappuyns_internal_concordance.png'}")

    # ═══════════════════════════════════════════════════════════════
    # Summary stats table
    # ═══════════════════════════════════════════════════════════════
    stats = {
        "n_patients": len(merged),
        "chi2": chi2, "p_chi2": p_chi,
        "cohens_kappa_binary": kappa,
        "internal_subtypes": merged["immune_subtype"].value_counts().to_dict(),
        "cappuyns_subsets": merged["cappuyns_subset"].value_counts().to_dict(),
    }
    pd.DataFrame([stats]).to_csv(TAB / "concordance_stats.csv", index=False)
    print(f"  Saved: {TAB / 'concordance_stats.csv'}")

    print(f"\n{'=' * 70}")
    print("Done.")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
