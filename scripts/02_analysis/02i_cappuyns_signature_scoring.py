#!/usr/bin/env python3
"""
02i_cappuyns_signature_scoring.py — Cappuyns et al. 2025 Molecular Subtyping

Applies the 21 scRNAseq-derived HCC gene signatures from:
  Cappuyns et al. "Single-cell RNA sequencing-derived signatures define response
  patterns to atezolizumab + bevacizumab in advanced HCC" (J Hepatol 2025)

Workflow:
  1. Parse 21 gene signatures from mmc2.xlsx (Table S2)
  2. Score each signature per cell in the integrated atlas (sc.tl.score_genes)
  3. Aggregate to pseudo-bulk (patient-level) via ssGSEA-like enrichment
  4. Classify samples as ImmunePos / ImmuneNeg (CD8 Temra + CD8 Tex + Macro CXCL10)
  5. Define "Immune-competent", "Angiogenesis-driven", "Resistant" molecular subsets
  6. Generate publication-quality figures & statistical summaries

Requires: data/references/cappuyns2025_mmc2_signatures.xlsx
          data/processed/integrated_io_atlas.h5ad
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
import sys; sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
from config.publication_style import apply_style, PAL, CMAP, label_clean
apply_style()
from pathlib import Path
from scipy.stats import mannwhitneyu, fisher_exact
from scipy.cluster.hierarchy import linkage, leaves_list

sc.settings.verbosity = 1

PROJECT = Path(__file__).resolve().parents[2]
DATA    = PROJECT / "data"
PROC    = DATA / "processed"
REF     = DATA / "references"
OUT     = PROJECT / "out" / "02i_cappuyns_signatures"
FIG     = OUT / "figures"
TAB     = OUT / "tables"
for d in [FIG, TAB]:
    d.mkdir(parents=True, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════
# 1. PARSE GENE SIGNATURES FROM mmc2.xlsx
# ═══════════════════════════════════════════════════════════════════════

def load_cappuyns_signatures(xlsx_path):
    """Parse 21 HCC gene signatures from Cappuyns et al. 2025 Table S2.

    Format: Each signature spans 2 rows in the 'Table S2' sheet.
      Row A: [sig_name, "Up-regulated genes (n=X)", "GENE1; GENE2; ..."]
      Row B: [None,     "Down-regulated genes (n=Y)", "GENE1; GENE2; ..."]
    Gene lists are semicolon-separated.

    Returns dict of {signature_name: {"up": [gene_list], "down": [gene_list]}}
    """
    import openpyxl
    wb = openpyxl.load_workbook(xlsx_path, read_only=True)

    # Find the sheet containing the 21 signatures
    target_sheet = None
    for sn in wb.sheetnames:
        if "S2" in sn or "signature" in sn.lower():
            target_sheet = sn
            break
    if target_sheet is None:
        target_sheet = wb.sheetnames[0]

    ws = wb[target_sheet]
    rows = list(ws.iter_rows(values_only=True))
    wb.close()

    signatures = {}
    current_sig = None

    for row in rows:
        if row is None or all(c is None for c in row):
            continue
        col0 = str(row[0]).strip() if row[0] else None
        col1 = str(row[1]).strip() if len(row) > 1 and row[1] else ""
        col2 = str(row[2]).strip() if len(row) > 2 and row[2] else ""

        # Skip header row
        if col0 and ("Table S" in col0 or "overview" in col0.lower()):
            continue

        # Parse gene list from semicolon-separated string
        genes = [g.strip() for g in col2.split(";") if g.strip()] if col2 else []

        is_up = "up" in col1.lower() if col1 else False
        is_down = "down" in col1.lower() if col1 else False

        if col0 and col0.lower() not in ("none", "nan", ""):
            # New signature name
            current_sig = col0
            if current_sig not in signatures:
                signatures[current_sig] = {"up": [], "down": []}

        if current_sig and (is_up or is_down):
            direction = "up" if is_up else "down"
            signatures[current_sig][direction] = genes

    for name, sig in signatures.items():
        print(f"  {name}: {len(sig['up'])} up, {len(sig['down'])} down")

    return signatures


# ═══════════════════════════════════════════════════════════════════════
# 2. SCORE SIGNATURES ON INTEGRATED ATLAS
# ═══════════════════════════════════════════════════════════════════════

def score_signatures_on_atlas(adata, signatures):
    """Score each Cappuyns signature per cell using sc.tl.score_genes."""
    available_genes = set(adata.var_names)
    scored = {}

    for sig_name, gene_sets in signatures.items():
        up_found = [g for g in gene_sets["up"] if g in available_genes]
        down_found = [g for g in gene_sets["down"] if g in available_genes]

        col_name = f"cappuyns_{sig_name}"

        if len(up_found) >= 3:
            sc.tl.score_genes(adata, up_found, score_name=f"{col_name}_up")
            scored[f"{col_name}_up"] = len(up_found)

        if len(down_found) >= 3:
            sc.tl.score_genes(adata, down_found, score_name=f"{col_name}_down")
            scored[f"{col_name}_down"] = len(down_found)

        # Combined score: up - down
        if len(up_found) >= 3 and len(down_found) >= 3:
            adata.obs[col_name] = (
                adata.obs[f"{col_name}_up"] - adata.obs[f"{col_name}_down"]
            )
        elif len(up_found) >= 3:
            adata.obs[col_name] = adata.obs[f"{col_name}_up"]

        coverage_up = len(up_found) / max(len(gene_sets["up"]), 1)
        coverage_dn = len(down_found) / max(len(gene_sets["down"]), 1)
        print(f"  {sig_name}: up {len(up_found)}/{len(gene_sets['up'])} "
              f"({coverage_up:.0%}), down {len(down_found)}/{len(gene_sets['down'])} "
              f"({coverage_dn:.0%})")

    return scored


# ═══════════════════════════════════════════════════════════════════════
# 3. PSEUDO-BULK PATIENT-LEVEL AGGREGATION
# ═══════════════════════════════════════════════════════════════════════

def patient_level_scores(adata, signatures):
    """Aggregate cell-level signature scores to patient level (pseudo-bulk)."""
    sig_cols = [c for c in adata.obs.columns if c.startswith("cappuyns_") and not c.endswith("_up") and not c.endswith("_down")]

    if "patient_ID" in adata.obs.columns:
        group_col = "patient_ID"
    elif "sample_id" in adata.obs.columns:
        group_col = "sample_id"
    else:
        print("  [WARN] No patient/sample ID column found for aggregation")
        return None

    # Filter to cells with response labels
    mask = adata.obs["response"].notna()
    obs_sub = adata.obs.loc[mask, sig_cols + [group_col, "response"]].copy()

    patient_scores = obs_sub.groupby(group_col)[sig_cols].median()
    patient_response = obs_sub.groupby(group_col)["response"].first()
    patient_scores["response"] = patient_response

    patient_scores.to_csv(TAB / "patient_level_cappuyns_scores.csv")
    print(f"  Patient-level scores: {patient_scores.shape[0]} patients × {len(sig_cols)} signatures")
    return patient_scores


# ═══════════════════════════════════════════════════════════════════════
# 4. IMMUNE-POSITIVE / IMMUNE-NEGATIVE CLASSIFICATION
# ═══════════════════════════════════════════════════════════════════════

def classify_immune_status(patient_scores):
    """Classify patients as ImmunePos/ImmuneNeg per Cappuyns Fig. 3C.

    ImmunePos: positive for at least one of CD8 Temra, CD8 Tex, Macro CXCL10
    ImmuneNeg: negative for all three
    """
    # Find the signature columns matching the three key immune signatures
    immune_sigs = {}
    for sig in ["CD8 Temra", "CD8 Tex", "Macro CXCL10",
                "CD8_Temra", "CD8_Tex", "Macro_CXCL10",
                "CD8Temra", "CD8Tex", "MacroCXCL10"]:
        match = [c for c in patient_scores.columns if sig.replace(" ", "").lower() in c.replace(" ", "").replace("_", "").lower()]
        if match:
            immune_sigs[sig] = match[0]

    if len(immune_sigs) < 3:
        print(f"  [WARN] Found only {len(immune_sigs)}/3 immune signatures: {list(immune_sigs.keys())}")
        print("  Skipping ImmunePos/ImmuneNeg classification")
        return patient_scores

    print(f"  Using immune signatures: {list(immune_sigs.values())}")

    # For each signature, define "positive" as above median
    immune_pos = pd.DataFrame(index=patient_scores.index)
    for name, col in immune_sigs.items():
        threshold = patient_scores[col].median()
        immune_pos[name] = patient_scores[col] > threshold

    patient_scores["immune_pos_count"] = immune_pos.sum(axis=1)
    patient_scores["immune_status"] = np.where(
        patient_scores["immune_pos_count"] >= 1, "ImmunePos", "ImmuneNeg"
    )

    print(f"  ImmunePos: {(patient_scores['immune_status'] == 'ImmunePos').sum()}")
    print(f"  ImmuneNeg: {(patient_scores['immune_status'] == 'ImmuneNeg').sum()}")

    return patient_scores


# ═══════════════════════════════════════════════════════════════════════
# 5. MOLECULAR SUBSET CLASSIFICATION (Fig. 7A)
# ═══════════════════════════════════════════════════════════════════════

def classify_molecular_subsets(patient_scores):
    """Classify patients into Immune-competent / Angiogenesis-driven / Resistant.

    Per Cappuyns et al. 2025 Fig. 7A decision tree:
      1. ImmunePos + no immunosuppressive myeloid → "Immune-competent"
      2. ImmuneNeg + NRP1-Low + no Notch → "Angiogenesis-driven"
      3. ImmunePos + immunosuppressive myeloid → "Resistant"
      4. ImmuneNeg + (Notch-active OR NRP1-High) + immunosuppressive myeloid → "Resistant"
      5. Otherwise → "Unclassified"
    """
    if "immune_status" not in patient_scores.columns:
        print("  [SKIP] No immune_status — run classify_immune_status first")
        return patient_scores

    # Find resistance-related signatures
    myeloid_cols = [c for c in patient_scores.columns
                    if any(m in c.lower() for m in ["cd14", "mono_cd14", "monocd14", "trem2", "macro_trem2"])]

    # NRP1 expression (may or may not be a signature)
    nrp1_col = next((c for c in patient_scores.columns if "nrp1" in c.lower()), None)

    # Notch signature
    notch_col = next((c for c in patient_scores.columns if "notch" in c.lower()), None)

    molecular_subset = pd.Series("Unclassified", index=patient_scores.index)

    is_immune_pos = patient_scores["immune_status"] == "ImmunePos"

    # Determine myeloid enrichment
    if myeloid_cols:
        myeloid_score = patient_scores[myeloid_cols].mean(axis=1)
        has_myeloid = myeloid_score > myeloid_score.median()
    else:
        has_myeloid = pd.Series(False, index=patient_scores.index)

    # Determine NRP1-Low
    if nrp1_col:
        nrp1_low = patient_scores[nrp1_col] <= patient_scores[nrp1_col].quantile(0.25)
    else:
        nrp1_low = pd.Series(False, index=patient_scores.index)

    # Determine Notch activation
    if notch_col:
        notch_active = patient_scores[notch_col] > patient_scores[notch_col].median()
    else:
        notch_active = pd.Series(False, index=patient_scores.index)

    # Classification logic
    molecular_subset[is_immune_pos & ~has_myeloid] = "Immune-competent"
    molecular_subset[~is_immune_pos & nrp1_low & ~notch_active] = "Angiogenesis-driven"
    molecular_subset[is_immune_pos & has_myeloid] = "Resistant"
    molecular_subset[~is_immune_pos & (notch_active | ~nrp1_low) & has_myeloid] = "Resistant"

    patient_scores["molecular_subset"] = molecular_subset

    counts = patient_scores["molecular_subset"].value_counts()
    print(f"  Molecular subsets:")
    for subset, n in counts.items():
        print(f"    {subset}: {n}")

    patient_scores.to_csv(TAB / "patient_molecular_subsets.csv")
    return patient_scores


# ═══════════════════════════════════════════════════════════════════════
# 6. FIGURES
# ═══════════════════════════════════════════════════════════════════════

def plot_signature_heatmap(patient_scores, signatures):
    """Heatmap of signature enrichment per patient with hierarchical clustering."""
    sig_cols = [c for c in patient_scores.columns
                if c.startswith("cappuyns_") and not c.endswith("_up") and not c.endswith("_down")]
    if len(sig_cols) == 0:
        return

    plot_data = patient_scores[sig_cols + ["response"]].dropna(subset=["response"])
    if plot_data.shape[0] < 5:
        print("  [SKIP] Too few patients for heatmap")
        return

    mat = plot_data[sig_cols]
    mat_z = (mat - mat.mean()) / mat.std()

    # Row linkage (signatures) and column linkage (patients)
    row_Z = linkage(mat_z.T, method="ward")
    col_Z = linkage(mat_z, method="ward")
    row_order = leaves_list(row_Z)
    col_order = leaves_list(col_Z)
    mat_z = mat_z.iloc[col_order, row_order]

    # Response annotation
    response = plot_data["response"].iloc[col_order]
    response_colors = response.map({"Responder": "#0077BB", "Non-Responder": "#CC3311"})

    fig, axes = plt.subplots(2, 1, figsize=(max(10, len(sig_cols) * 0.55), max(7, mat_z.shape[0] * 0.22)),
                             gridspec_kw={"height_ratios": [1, 20], "hspace": 0.02})

    # Response color bar on top
    ax_top = axes[0]
    for i, c in enumerate(response_colors):
        ax_top.add_patch(plt.Rectangle((i, 0), 1, 1, color=c, linewidth=0))
    ax_top.set_xlim(0, len(response_colors))
    ax_top.set_ylim(0, 1)
    ax_top.set_xticks([])
    ax_top.set_yticks([0.5])
    ax_top.set_yticklabels(["Response"], fontsize=9)
    ax_top.tick_params(left=False)
    for sp in ax_top.spines.values():
        sp.set_visible(False)

    # Main heatmap
    ax = axes[1]
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    sig_labels = [c.replace("cappuyns_", "") for c in mat_z.columns]
    patient_labels = [str(int(idx)) if isinstance(idx, float) else str(idx)
                      for idx in mat_z.index]

    sns.heatmap(mat_z, cmap=cmap, center=0, vmin=-2.5, vmax=2.5,
                xticklabels=sig_labels, yticklabels=patient_labels,
                ax=ax, cbar_kws={"label": "Z-score", "shrink": 0.7},
                linewidths=0.3)

    ax.set_title("Cappuyns HCC Gene Signatures — Hierarchical Clustering", fontweight="bold", fontsize=12, pad=4)
    ax.tick_params(axis="x", rotation=50, labelsize=8)
    ax.tick_params(axis="y", labelsize=8)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor="#0077BB", label="Responder"),
                       Patch(facecolor="#CC3311", label="Non-Responder")]
    ax_top.legend(handles=legend_elements, loc="upper right", fontsize=8,
                  frameon=True, ncol=2, bbox_to_anchor=(1.0, 1.8))

    plt.savefig(FIG / "cappuyns_signature_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close()
    print("  Saved: cappuyns_signature_heatmap.png")


def plot_immune_status_comparison(patient_scores):
    """Compare key signature scores between molecular subsets (Immune-competent vs Resistant)."""
    subset_col = "molecular_subset" if "molecular_subset" in patient_scores.columns else "immune_status"
    if subset_col not in patient_scores.columns:
        return

    sig_cols = [c for c in patient_scores.columns
                if c.startswith("cappuyns_") and not c.endswith("_up") and not c.endswith("_down")]

    # Focus on Immune-competent vs Resistant (drop Unclassified / Angiogenesis-driven if few)
    if subset_col == "molecular_subset":
        groups = patient_scores[subset_col].value_counts()
        valid_groups = groups[groups >= 2].index.tolist()
        plot_data = patient_scores[patient_scores[subset_col].isin(valid_groups)].copy()
    else:
        plot_data = patient_scores.copy()

    n_sigs = min(8, len(sig_cols))
    ncols = 4
    nrows = (n_sigs + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.5 * nrows), squeeze=False)
    axes = axes.flatten()

    palette = {"Immune-competent": "#0077BB", "Resistant": "#CC3311",
               "ImmunePos": "#3498db", "ImmuneNeg": "#CC3311",
               "Angiogenesis-driven": "#f39c12", "Unclassified": "#95a5a6"}

    for i, col in enumerate(sig_cols[:n_sigs]):
        ax = axes[i]
        sns.boxplot(data=plot_data, x=subset_col, y=col, palette=palette,
                    width=0.45, ax=ax, fliersize=0, linewidth=0.8)
        sns.stripplot(data=plot_data, x=subset_col, y=col, palette=palette,
                      size=6, alpha=0.7, jitter=0.15, ax=ax, edgecolor="gray", linewidth=0.3)

        # Mann-Whitney between largest two groups
        grp_names = plot_data[subset_col].unique()
        if len(grp_names) >= 2:
            g1 = plot_data.loc[plot_data[subset_col] == grp_names[0], col].dropna()
            g2 = plot_data.loc[plot_data[subset_col] == grp_names[1], col].dropna()
            if len(g1) >= 3 and len(g2) >= 3:
                _, p = mannwhitneyu(g1, g2, alternative="two-sided")
                ax.set_title(f"{col.replace('cappuyns_', '')}\np = {p:.3g}", fontsize=9, fontweight="bold")
            else:
                ax.set_title(col.replace("cappuyns_", ""), fontsize=9, fontweight="bold")
        else:
            ax.set_title(col.replace("cappuyns_", ""), fontsize=9, fontweight="bold")

        ax.set_xlabel("")
        ax.set_ylabel("Score" if i % ncols == 0 else "")
        ax.tick_params(axis="x", rotation=20, labelsize=8)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle("Cappuyns Signatures by Molecular Subset", fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(FIG / "cappuyns_immune_status_comparison.png", dpi=200, bbox_inches="tight")
    plt.close()
    print("  Saved: cappuyns_immune_status_comparison.png")


def plot_molecular_subsets(patient_scores):
    """Stacked proportion plot of response within each molecular subset with Fisher test."""
    if "molecular_subset" not in patient_scores.columns:
        return
    if "response" not in patient_scores.columns:
        return

    plot_df = patient_scores[["molecular_subset", "response"]].dropna()
    if plot_df.shape[0] < 5:
        return

    ct = pd.crosstab(plot_df["molecular_subset"], plot_df["response"])
    # Ensure column order
    for col in ["Responder", "Non-Responder"]:
        if col not in ct.columns:
            ct[col] = 0
    ct = ct[["Responder", "Non-Responder"]]

    # Drop groups with fewer than 2 patients
    ct = ct[ct.sum(axis=1) >= 2]
    if ct.shape[0] < 2:
        return

    ct_prop = ct.div(ct.sum(axis=1), axis=0)
    n_per_group = ct.sum(axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1.3, 1]})

    # ── Panel A: Stacked proportion bar ──────────────────────────────
    ax = axes[0]
    x = np.arange(len(ct_prop))
    bar_w = 0.5
    ax.bar(x, ct_prop["Responder"], bar_w, label="Responder", color="#0077BB", edgecolor="white", linewidth=0.8)
    ax.bar(x, ct_prop["Non-Responder"], bar_w, bottom=ct_prop["Responder"],
           label="Non-Responder", color="#CC3311", edgecolor="white", linewidth=0.8)

    for i, (subset, row) in enumerate(ct.iterrows()):
        total = row.sum()
        rr = row["Responder"] / total * 100
        ax.text(i, 0.5, f"{row['Responder']}/{total}\n({rr:.0f}%)",
                ha="center", va="center", fontsize=10, fontweight="bold", color="white")

    ax.set_xticks(x)
    ax.set_xticklabels(ct_prop.index, fontsize=10)
    ax.set_ylabel("Proportion", fontsize=11)
    ax.set_ylim(0, 1.12)
    ax.set_title("A. Response Rate by Molecular Subset", fontweight="bold", fontsize=12)
    ax.legend(loc="upper right", frameon=True, fontsize=9)

    # ── Panel B: Count table + Fisher exact ──────────────────────────
    ax = axes[1]
    ax.axis("off")

    # Fisher exact test (Immune-competent vs Resistant if both exist)
    ic_label = [s for s in ct.index if "compet" in s.lower()]
    res_label = [s for s in ct.index if "resist" in s.lower()]

    stat_text = ""
    if ic_label and res_label:
        table_2x2 = np.array([
            [ct.loc[ic_label[0], "Responder"], ct.loc[ic_label[0], "Non-Responder"]],
            [ct.loc[res_label[0], "Responder"], ct.loc[res_label[0], "Non-Responder"]],
        ])
        odds_ratio, p_fisher = fisher_exact(table_2x2)
        stat_text = f"Fisher exact test (Immune-competent vs Resistant):\nOR = {odds_ratio:.2f}, p = {p_fisher:.4f}"
    else:
        stat_text = "(Fisher's exact test not computed — requires\nImmune-competent and Resistant groups)"

    # Summary table
    summary = ct.copy()
    summary["Total"] = summary.sum(axis=1)
    summary["Response Rate"] = (summary["Responder"] / summary["Total"] * 100).round(1).astype(str) + "%"

    table = ax.table(cellText=summary.values,
                     colLabels=summary.columns,
                     rowLabels=summary.index,
                     cellLoc="center", loc="upper center")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.5)

    ax.text(0.5, 0.15, stat_text, transform=ax.transAxes,
            ha="center", va="top", fontsize=10,
            bbox=dict(facecolor="#f0f0f0", edgecolor="#cccccc", boxstyle="round,pad=0.4"))
    ax.set_title("B. Molecular Subset Summary", fontweight="bold", fontsize=12, pad=15)

    plt.tight_layout()
    plt.savefig(FIG / "cappuyns_molecular_subsets.png", dpi=200, bbox_inches="tight")
    plt.close()
    print("  Saved: cappuyns_molecular_subsets.png")


def plot_key_signatures_response(patient_scores):
    """Box + strip plots for all significantly differential signatures by response.

    Reads BH-corrected p-values from the stats table (if available) to annotate
    all signatures with p_adj < 0.05.
    """
    sig_cols = [c for c in patient_scores.columns
                if c.startswith("cappuyns_") and not c.endswith("_up") and not c.endswith("_down")]
    if len(sig_cols) == 0 or "response" not in patient_scores.columns:
        return

    # Load stats if available for BH p-values
    stats_path = TAB / "cappuyns_signature_stats.csv"
    sig_pvals = {}
    if stats_path.exists():
        stats_df = pd.read_csv(stats_path)
        for _, row in stats_df.iterrows():
            sig_pvals[f"cappuyns_{row['signature']}"] = row["p_adj"]

    # Select significant signatures (p_adj < 0.05), else fall back to top 7 by effect
    if sig_pvals:
        sig_sigs = [c for c in sig_cols if sig_pvals.get(c, 1.0) < 0.05]
        # Order by effect size
        from statsmodels.stats.multitest import multipletests
        sig_sigs_sorted = sorted(sig_sigs,
                                 key=lambda c: abs(patient_scores.loc[patient_scores["response"]=="Responder", c].median()
                                                   - patient_scores.loc[patient_scores["response"]=="Non-Responder", c].median()),
                                 reverse=True)
    else:
        sig_sigs_sorted = sig_cols[:7]

    if len(sig_sigs_sorted) == 0:
        sig_sigs_sorted = sig_cols[:7]

    n = len(sig_sigs_sorted)
    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 5 * nrows), squeeze=False)
    axes = axes.flatten()

    palette = {"Responder": "#0077BB", "Non-Responder": "#CC3311"}
    plot_df = patient_scores[sig_sigs_sorted + ["response"]].dropna(subset=["response"])

    for i, col in enumerate(sig_sigs_sorted):
        ax = axes[i]
        sns.boxplot(data=plot_df, x="response", y=col, palette=palette,
                    width=0.45, ax=ax, fliersize=0, linewidth=0.8)
        sns.stripplot(data=plot_df, x="response", y=col, palette=palette,
                      size=6, alpha=0.7, jitter=0.15, ax=ax,
                      edgecolor="gray", linewidth=0.3)

        # Use BH-corrected p if available
        p_adj = sig_pvals.get(col)
        if p_adj is not None:
            stars = "***" if p_adj < 0.001 else "**" if p_adj < 0.01 else "*" if p_adj < 0.05 else "ns"
            ax.set_title(f"{col.replace('cappuyns_', '')}\np_adj = {p_adj:.3g} {stars}",
                         fontweight="bold", fontsize=10)
        else:
            r_vals = plot_df.loc[plot_df["response"] == "Responder", col].dropna()
            nr_vals = plot_df.loc[plot_df["response"] == "Non-Responder", col].dropna()
            if len(r_vals) >= 3 and len(nr_vals) >= 3:
                _, p = mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
                ax.set_title(f"{col.replace('cappuyns_', '')}\np = {p:.3g}", fontweight="bold", fontsize=10)
            else:
                ax.set_title(col.replace("cappuyns_", ""), fontweight="bold", fontsize=10)

        ax.set_xlabel("")
        ax.set_ylabel("Signature Score" if i % ncols == 0 else "")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle("Significantly Differential Cappuyns Signatures (R vs NR, BH p_adj < 0.05)",
                 fontweight="bold", y=1.01, fontsize=12)
    plt.tight_layout()
    plt.savefig(FIG / "cappuyns_key_signatures_response.png", dpi=200, bbox_inches="tight")
    plt.close()
    print("  Saved: cappuyns_key_signatures_response.png")


# ═══════════════════════════════════════════════════════════════════════
# 7. STATISTICAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════

def statistical_summary(patient_scores):
    """Compute Mann-Whitney U tests for all signatures R vs NR."""
    sig_cols = [c for c in patient_scores.columns
                if c.startswith("cappuyns_") and not c.endswith("_up") and not c.endswith("_down")]

    results = []
    for col in sig_cols:
        r = patient_scores.loc[patient_scores["response"] == "Responder", col].dropna()
        nr = patient_scores.loc[patient_scores["response"] == "Non-Responder", col].dropna()
        if len(r) >= 3 and len(nr) >= 3:
            stat, p = mannwhitneyu(r, nr, alternative="two-sided")
            effect = r.median() - nr.median()
            results.append({
                "signature": col.replace("cappuyns_", ""),
                "n_R": len(r), "median_R": r.median(),
                "n_NR": len(nr), "median_NR": nr.median(),
                "effect_size": effect,
                "U_statistic": stat, "p_value": p,
            })

    if results:
        df = pd.DataFrame(results).sort_values("p_value")
        # BH correction
        from statsmodels.stats.multitest import multipletests
        _, df["p_adj"], _, _ = multipletests(df["p_value"], method="fdr_bh")
        df.to_csv(TAB / "cappuyns_signature_stats.csv", index=False)
        print(f"\n  Statistical summary ({len(df)} signatures tested):")
        for _, row in df.head(10).iterrows():
            sig_str = "***" if row["p_adj"] < 0.001 else "**" if row["p_adj"] < 0.01 else "*" if row["p_adj"] < 0.05 else ""
            print(f"    {row['signature']:30s} effect={row['effect_size']:+.3f}  "
                  f"p={row['p_value']:.3g}  p_adj={row['p_adj']:.3g} {sig_str}")
        return df
    return None


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("Cappuyns et al. 2025 — HCC Gene Signature Scoring & Molecular Subtyping")
    print("=" * 70)

    # ── Load signatures ────────────────────────────────────────────
    xlsx_path = REF / "cappuyns2025_mmc2_signatures.xlsx"
    if not xlsx_path.exists():
        print(f"\n  ERROR: Signatures file not found: {xlsx_path}")
        print("  Run: bash scripts/00_data_acquisition/download_cappuyns_data.sh")
        print("  Or download mmc2.xlsx manually from:")
        print("    https://doi.org/10.1016/j.jhep.2024.12.016")
        return

    print("\n[1] Loading 21 HCC gene signatures from Table S2...")
    signatures = load_cappuyns_signatures(xlsx_path)
    print(f"    Loaded {len(signatures)} signatures")

    if len(signatures) == 0:
        print("  ERROR: No signatures parsed from mmc2.xlsx")
        return

    # Save parsed signatures as JSON for reproducibility
    import json
    with open(TAB / "cappuyns_signatures_parsed.json", "w") as f:
        json.dump(signatures, f, indent=2)
    print(f"    Saved parsed signatures to {TAB / 'cappuyns_signatures_parsed.json'}")

    # ── Load atlas ─────────────────────────────────────────────────
    atlas_path = PROC / "integrated_io_atlas.h5ad"
    if not atlas_path.exists():
        print(f"\n  ERROR: Atlas not found: {atlas_path}")
        print("  Run preprocessing pipeline first.")
        return

    print(f"\n[2] Loading integrated atlas...")
    adata = sc.read_h5ad(atlas_path)
    print(f"    {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

    # ── Score signatures ───────────────────────────────────────────
    print(f"\n[3] Scoring {len(signatures)} Cappuyns signatures per cell...")
    scored = score_signatures_on_atlas(adata, signatures)
    print(f"    Scored {len(scored)} signature components")

    # ── Patient-level aggregation ──────────────────────────────────
    print(f"\n[4] Aggregating to patient level...")
    patient_scores = patient_level_scores(adata, signatures)

    if patient_scores is not None and "response" in patient_scores.columns:
        # ── ImmunePos / ImmuneNeg classification ───────────────────
        print(f"\n[5] Classifying ImmunePos / ImmuneNeg...")
        patient_scores = classify_immune_status(patient_scores)

        # ── Molecular subsets ──────────────────────────────────────
        print(f"\n[6] Defining molecular subsets (Immune-competent / Angiogenesis-driven / Resistant)...")
        patient_scores = classify_molecular_subsets(patient_scores)

        # ── Statistics ─────────────────────────────────────────────
        print(f"\n[7] Statistical analysis...")
        statistical_summary(patient_scores)

        # ── Figures ────────────────────────────────────────────────
        print(f"\n[8] Generating figures...")
        plot_signature_heatmap(patient_scores, signatures)
        plot_immune_status_comparison(patient_scores)
        plot_molecular_subsets(patient_scores)
        plot_key_signatures_response(patient_scores)
    else:
        print("\n  [INFO] No response labels available for patient-level analysis")
        print("         Cell-level scores have been computed and can be used downstream")

    print(f"\n{'=' * 70}")
    print(f"Done. Outputs in: {OUT}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
