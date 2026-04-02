#!/usr/bin/env python3
"""
04_complementary_analyses.py — Complementary Downstream Analyses

Additional analyses to strengthen the findings:
  1. Gene Set Enrichment Analysis (GSEA) on R vs NR DEGs per cell type
  2. Cell-cell communication inference (ligand-receptor analysis)
  3. CD8 T cell state characterization (exhaustion trajectory)
  4. Clonotype diversity analysis (Shannon, Simpson, Chao1)
  5. Immune cell ratio biomarkers (CD8/Treg, M1/M2, etc.)
  6. Patient stratification by immune subtypes
  7. Cross-dataset reproducibility of key signatures
  8. Multi-panel summary figure
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
from scipy.stats import mannwhitneyu, spearmanr, fisher_exact, entropy
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler
from collections import Counter
import re

sc.settings.verbosity = 1

PROJECT = Path(__file__).resolve().parents[2]
PROC    = PROJECT / "data" / "processed"
DATA    = PROJECT / "data"
OUT     = PROJECT / "out" / "04_complementary"
FIG     = OUT / "figures"
TAB     = OUT / "tables"
for d in [FIG, TAB]:
    d.mkdir(parents=True, exist_ok=True)

IO_SIGNATURES = {
    "T_cell_exhaustion":    ["PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","TOX","ENTPD1"],
    "T_cell_activation":    ["IFNG","GZMB","PRF1","TNF","CD69","IL2"],
    "T_cell_cytotoxicity":  ["GZMA","GZMB","GZMH","GZMK","PRF1","GNLY","NKG7"],
    "Progenitor_CD8":       ["TCF7","SLAMF6","LEF1","CCR7","IL7R","SELL"],
    "Terminal_exhaustion":  ["TOX","ENTPD1","HAVCR2","LAYN","TIGIT","CXCL13"],
    "Treg_signature":       ["FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"],
    "M1_macrophage":        ["NOS2","IL1B","IL6","TNF","CD80","CD86"],
    "M2_macrophage":        ["CD163","MRC1","MSR1","CD200R1","IL10"],
    "IFN_gamma":            ["IFNG","STAT1","IRF1","CXCL9","CXCL10","CXCL11","IDO1"],
    "Antigen_presentation": ["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","PSMB9"],
    "Immune_checkpoint":    ["CD274","PDCD1LG2","PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","VSIR"],
    "mregDC":               ["CCR7","LAMP3","CCL19","CCL22","CD80","CD83","FSCN1","MARCKSL1"],
    "Tfh_like":             ["CXCL13","CH25H","IL21","PDCD1","ICOS","BCL6","TOX2"],
}

# Curated pathway gene sets for GSEA-style analysis
HALLMARK_PATHWAYS = {
    "Interferon_alpha":     ["ISG15","MX1","MX2","IFIT1","IFIT2","IFIT3","OAS1","OAS2","IFI44","IFI44L","IFI6","RSAD2"],
    "Interferon_gamma":     ["IFNG","STAT1","IRF1","CXCL9","CXCL10","CXCL11","IDO1","GBP1","GBP2","GBP4","GBP5"],
    "TNF_via_NFkB":         ["NFKBIA","TNFAIP3","CXCL1","CXCL2","CCL2","CCL5","IL6","IL1B","ICAM1","VCAM1"],
    "IL6_JAK_STAT3":        ["IL6","IL6R","STAT3","SOCS3","BCL2L1","PIM1","MYC"],
    "Complement":           ["C1QA","C1QB","C1QC","C3","C4A","C4B","CFB","CFD","CFH","CFI","SERPING1"],
    "Coagulation":          ["FGA","FGB","FGG","F2","F5","F7","F10","SERPINC1","PROC","PROS1"],
    "Epithelial_mesenchymal": ["VIM","CDH2","FN1","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2","MMP2","MMP9"],
    "Glycolysis":           ["HK1","HK2","PKM","LDHA","GAPDH","ENO1","ALDOA","TPI1","PGK1","PFKP"],
    "Oxidative_phosphorylation": ["NDUFA1","NDUFB1","UQCRC1","COX5A","ATP5F1A","SDHB","CYC1"],
    "Fatty_acid_metabolism": ["ACACA","FASN","SCD","ACADM","ACADL","CPT1A","CPT2","HADHA"],
    "Apoptosis":            ["BAX","BAK1","BCL2","BCL2L1","MCL1","CASP3","CASP8","CASP9","BID","PARP1"],
    "Cell_cycle":           ["CDK1","CDK2","CDK4","CCNA2","CCNB1","CCND1","CCNE1","RB1","E2F1","MKI67"],
    "DNA_repair":           ["BRCA1","BRCA2","ATM","ATR","RAD51","MLH1","MSH2","XRCC1","PARP1"],
    "Hypoxia":              ["HIF1A","VEGFA","LDHA","PGK1","ENO1","ALDOA","SLC2A1","CA9","BNIP3"],
    "Angiogenesis":         ["VEGFA","VEGFB","VEGFC","FLT1","KDR","ANGPT1","ANGPT2","TEK","PDGFB"],
    "TGF_beta":             ["TGFB1","TGFB2","TGFB3","SMAD2","SMAD3","SMAD4","ACVR1","BMPR2"],
    "WNT_beta_catenin":     ["CTNNB1","APC","AXIN2","LGR5","WNT3A","WNT5A","FZD1","LEF1","TCF7L2"],
    "Notch_signaling":      ["NOTCH1","NOTCH2","JAG1","JAG2","DLL1","DLL4","HES1","HEY1"],
    "PI3K_AKT_mTOR":        ["PIK3CA","PIK3CB","AKT1","MTOR","PTEN","TSC1","TSC2","RPTOR","RICTOR"],
    "MYC_targets":          ["MYC","MYCN","ODC1","NME1","NPM1","NCL","FBL","RPS6","LDHA"],
}

# Ligand-receptor pairs (curated selection relevant to IO)
LR_PAIRS = [
    ("PDCD1", "CD274"),     # PD-1 / PD-L1
    ("PDCD1", "PDCD1LG2"),  # PD-1 / PD-L2
    ("CTLA4", "CD80"),      # CTLA-4 / B7-1
    ("CTLA4", "CD86"),      # CTLA-4 / B7-2
    ("HAVCR2", "LGALS9"),   # TIM-3 / Galectin-9
    ("LAG3", "HLA-DRA"),    # LAG-3 / MHC-II
    ("TIGIT", "PVR"),       # TIGIT / CD155
    ("TNFRSF9", "TNFSF9"),  # 4-1BB / 4-1BBL
    ("CD28", "CD80"),       # CD28 / B7-1
    ("CD28", "CD86"),       # CD28 / B7-2
    ("ICOS", "ICOSLG"),     # ICOS / ICOSL
    ("CXCR3", "CXCL9"),    # CXCR3 / CXCL9
    ("CXCR3", "CXCL10"),   # CXCR3 / CXCL10
    ("CCR5", "CCL5"),       # CCR5 / CCL5
    ("CXCR5", "CXCL13"),   # CXCR5 / CXCL13
    ("CCR7", "CCL19"),      # CCR7 / CCL19
    ("CCR7", "CCL21"),      # CCR7 / CCL21
    ("FASLG", "FAS"),       # FasL / Fas
    ("IFNG", "IFNGR1"),    # IFN-gamma / IFNGR1
    ("CD40LG", "CD40"),    # CD40L / CD40
    ("TGFB1", "TGFBR2"),  # TGFb / TGFbR
    ("IL10", "IL10RA"),    # IL-10 / IL-10R
    ("VEGFA", "KDR"),      # VEGF-A / VEGFR2
    ("TNF", "TNFRSF1A"),   # TNF / TNFR1
]


def recover_metadata(adata):
    """Recover GSE206325 metadata if missing."""
    if "response" in adata.obs.columns and adata.obs["response"].notna().sum() > 0:
        return adata
    print("[META] Recovering GSE206325 clinical metadata...")
    cell_meta = pd.read_csv(PROC / "GSE206325_metadata.csv", index_col=0)
    sample_ann = pd.read_csv(PROC / "GSE206325_sample_annotations.csv", index_col=0)
    sample_ann.index = sample_ann.index.astype(str)
    clust_ann = pd.read_csv(PROC / "GSE206325_cluster_annotations.csv", index_col=0)
    clust_ann["cluster_ID"] = clust_ann["cluster_ID"].astype(str)

    shared = adata.obs_names.intersection(cell_meta.index)
    print(f"  Matched {len(shared):,} GSE206325 barcodes")

    adata.obs["sample_id"] = np.nan
    adata.obs.loc[shared, "sample_id"] = cell_meta.loc[shared, "cell_to_sample_ID"].values.astype(str)
    adata.obs["cluster_id"] = np.nan
    adata.obs.loc[shared, "cluster_id"] = cell_meta.loc[shared, "cell_to_cluster"].values.astype(str)

    for col in ["patient_ID", "tissue", "treatment_Resp", "Immune_Infiltration"]:
        mapping = sample_ann[col].to_dict()
        adata.obs[col] = adata.obs["sample_id"].map(mapping)

    resp_map = {"anti-PD1_R": "Responder", "anti-PD1_NR": "Non-Responder"}
    adata.obs["response"] = adata.obs["treatment_Resp"].map(resp_map)

    ct_map = {}
    for _, row in clust_ann.iterrows():
        cid = str(row["cluster_ID"])
        subgroup = row["subgroup"] if pd.notna(row.get("subgroup")) else ""
        ct = row["type"] if pd.notna(row.get("type")) else ""
        group = row["group"] if pd.notna(row.get("group")) else ""
        if subgroup:
            ct_map[cid] = f"{ct}_{subgroup}" if ct else subgroup
        elif ct:
            ct_map[cid] = ct
        else:
            ct_map[cid] = group
    adata.obs["cell_type_original"] = adata.obs["cluster_id"].map(ct_map)
    return adata


# ═══════════════════════════════════════════════════════════════════════
# 1.  GSEA-STYLE PATHWAY SCORING
# ═══════════════════════════════════════════════════════════════════════

def pathway_enrichment_analysis(adata):
    """Score hallmark pathways and compare R vs NR."""
    print("\n[1] Hallmark Pathway Scoring (R vs NR)...")
    genes_available = set(adata.raw.var_names) if adata.raw else set(adata.var_names)

    scored_pathways = []
    for pw_name, genes in HALLMARK_PATHWAYS.items():
        present = [g for g in genes if g in genes_available]
        if len(present) >= 3:
            sc.tl.score_genes(adata, present, score_name=f"pw_{pw_name}", use_raw=True)
            scored_pathways.append(pw_name)

    print(f"  Scored {len(scored_pathways)}/{len(HALLMARK_PATHWAYS)} pathways")

    io_mask = adata.obs["response"].isin(["Responder", "Non-Responder"])
    io_cells = adata[io_mask]

    results = []
    for pw in scored_pathways:
        col = f"pw_{pw}"
        r_vals = io_cells.obs.loc[io_cells.obs["response"] == "Responder", col]
        nr_vals = io_cells.obs.loc[io_cells.obs["response"] == "Non-Responder", col]
        stat, pval = mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
        results.append({
            "pathway": pw, "mean_R": r_vals.mean(), "mean_NR": nr_vals.mean(),
            "effect": r_vals.mean() - nr_vals.mean(), "pvalue": pval,
        })

    pw_df = pd.DataFrame(results).sort_values("pvalue")
    pw_df["padj"] = np.minimum(pw_df["pvalue"] * len(pw_df), 1.0)
    pw_df.to_csv(TAB / "hallmark_pathway_R_vs_NR.csv", index=False)

    # Fig: Pathway effect sizes
    fig, ax = plt.subplots(figsize=(10, 8))
    pw_plot = pw_df.sort_values("effect")
    colors = ["#2ecc71" if e > 0 else "#e74c3c" for e in pw_plot["effect"]]
    stars = []
    for _, row in pw_plot.iterrows():
        s = "***" if row["padj"] < 0.001 else "**" if row["padj"] < 0.01 else "*" if row["padj"] < 0.05 else ""
        stars.append(s)
    ax.barh(range(len(pw_plot)), pw_plot["effect"], color=colors, alpha=0.85, height=0.5)
    ax.set_yticks(range(len(pw_plot)))
    ax.set_yticklabels(pw_plot["pathway"].str.replace("_", " "), fontsize=8)
    for i, (star, (_, row)) in enumerate(zip(stars, pw_plot.iterrows())):
        xoff = 0.002 * np.sign(row["effect"]) if row["effect"] != 0 else 0.002
        ax.text(row["effect"] + xoff, i, star, va="center", fontsize=10)
    ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Effect Size (R - NR)")
    ax.set_title("Hallmark Pathway Activity: Responder vs Non-Responder\n(* p<0.05, ** p<0.01, *** p<0.001)",
                 fontweight="bold")
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"04_hallmark_pathways.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    return pw_df, scored_pathways


# ═══════════════════════════════════════════════════════════════════════
# 2.  LIGAND-RECEPTOR ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

def ligand_receptor_analysis(adata):
    """Analyze ligand-receptor expression in R vs NR."""
    print("\n[2] Ligand-Receptor Communication Analysis...")
    genes_available = set(adata.raw.var_names) if adata.raw else set(adata.var_names)

    io_mask = adata.obs["response"].isin(["Responder", "Non-Responder"])
    io_cells = adata[io_mask]

    adata_raw = io_cells.raw.to_adata() if io_cells.raw else io_cells

    results = []
    for lig, rec in LR_PAIRS:
        if lig not in genes_available or rec not in genes_available:
            continue

        lig_idx = list(adata_raw.var_names).index(lig)
        rec_idx = list(adata_raw.var_names).index(rec)

        for resp in ["Responder", "Non-Responder"]:
            mask = io_cells.obs["response"] == resp
            X = adata_raw[mask].X
            if hasattr(X, "toarray"):
                X_dense = X.toarray()
            else:
                X_dense = X

            lig_expr = X_dense[:, lig_idx]
            rec_expr = X_dense[:, rec_idx]

            results.append({
                "ligand": lig, "receptor": rec,
                "response": resp,
                "lig_mean": lig_expr.mean(),
                "rec_mean": rec_expr.mean(),
                "lig_pct": (lig_expr > 0).mean() * 100,
                "rec_pct": (rec_expr > 0).mean() * 100,
                "lr_product": lig_expr.mean() * rec_expr.mean(),
            })

    lr_df = pd.DataFrame(results)
    lr_df.to_csv(TAB / "ligand_receptor_expression.csv", index=False)

    # Compute R vs NR log2FC of LR product
    lr_wide = lr_df.pivot_table(index=["ligand", "receptor"], columns="response",
                                 values="lr_product", aggfunc="first")
    if "Responder" in lr_wide.columns and "Non-Responder" in lr_wide.columns:
        lr_wide["log2FC"] = np.log2((lr_wide["Responder"] + 1e-8) / (lr_wide["Non-Responder"] + 1e-8))
        lr_wide = lr_wide.sort_values("log2FC")
        lr_wide.to_csv(TAB / "lr_interaction_scores.csv")

        # Fig: LR interaction comparison
        fig, ax = plt.subplots(figsize=(10, 8))
        labels = [f"{row.Index[0]}→{row.Index[1]}" for row in lr_wide.itertuples()]
        colors = ["#2ecc71" if fc > 0 else "#e74c3c" for fc in lr_wide["log2FC"]]
        ax.barh(range(len(lr_wide)), lr_wide["log2FC"], color=colors, alpha=0.85, height=0.5)
        ax.set_yticks(range(len(lr_wide)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
        ax.set_xlabel("log2 FC (R / NR) of LR Interaction Score")
        ax.set_title("Ligand-Receptor Interactions: Responder vs Non-Responder", fontweight="bold")
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"04_ligand_receptor.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

    # Fig: LR dotplot
    fig, axes = plt.subplots(1, 2, figsize=(14, 8))
    for idx, resp in enumerate(["Responder", "Non-Responder"]):
        ax = axes[idx]
        sub = lr_df[lr_df["response"] == resp].copy()
        sub["pair"] = sub["ligand"] + "→" + sub["receptor"]
        sub = sub.sort_values("lr_product", ascending=True)
        scatter = ax.scatter(sub["lig_pct"], range(len(sub)),
                            s=sub["lr_product"]*5000 + 5, c=sub["rec_pct"],
                            cmap="YlOrRd", alpha=0.8, edgecolors="gray", linewidths=0.5)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(sub["pair"], fontsize=7)
        ax.set_xlabel("Ligand % Expressing")
        ax.set_title(f"{resp}", fontweight="bold")
        plt.colorbar(scatter, ax=ax, label="Receptor % Expressing", shrink=0.7)
    plt.suptitle("Ligand-Receptor Communication: R vs NR", fontweight="bold", fontsize=13)
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"04_lr_dotplot.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    return lr_df


# ═══════════════════════════════════════════════════════════════════════
# 3.  CD8 T CELL STATE CHARACTERIZATION
# ═══════════════════════════════════════════════════════════════════════

def cd8_state_analysis(adata):
    """Detailed CD8 T cell state analysis: progenitor, effector, exhausted."""
    print("\n[3] CD8 T Cell State Characterization...")

    cd8_markers = ["CD8A", "CD8B", "CD3E"]
    genes_available = set(adata.var_names)

    # Identify CD8 T cells
    if "cell_type" in adata.obs.columns:
        cd8_types = [ct for ct in adata.obs["cell_type"].unique()
                     if any(k in str(ct).lower() for k in ["cd8", "cytotoxic", "ctl"])]
        if cd8_types:
            cd8_mask = adata.obs["cell_type"].isin(cd8_types)
        else:
            cd8_mask = adata.obs["cell_type"].str.contains("CD8|Cytotoxic", case=False, na=False)
    else:
        cd8_mask = pd.Series(False, index=adata.obs_names)

    io_mask = adata.obs["response"].isin(["Responder", "Non-Responder"])
    cd8_io = adata[cd8_mask & io_mask].copy()

    if cd8_io.shape[0] < 100:
        print("  [SKIP] Too few CD8 T cells with IO labels")
        return None

    print(f"  CD8 T cells with IO labels: {cd8_io.shape[0]:,}")

    # Define exhaustion trajectory markers
    states = {
        "Naive_like":      ["TCF7", "LEF1", "CCR7", "SELL", "IL7R"],
        "Progenitor_exh":  ["TCF7", "SLAMF6", "CXCR5", "PDCD1", "TOX"],
        "Effector":        ["GZMB", "PRF1", "GNLY", "NKG7", "IFNG", "GZMH"],
        "Transitional":    ["GZMK", "EOMES", "PDCD1", "CXCR3"],
        "Terminal_exh":    ["TOX", "ENTPD1", "HAVCR2", "LAYN", "TIGIT", "CXCL13"],
        "Tissue_resident": ["ITGAE", "ITGA1", "CD69", "CXCR6", "HOBIT"],
    }

    genes_raw = set(cd8_io.raw.var_names) if cd8_io.raw else set(cd8_io.var_names)
    scored_states = []
    for state, genes in states.items():
        present = [g for g in genes if g in genes_raw]
        if len(present) >= 2:
            sc.tl.score_genes(cd8_io, present, score_name=state, use_raw=True)
            scored_states.append(state)

    # Classify CD8 states
    if scored_states:
        score_cols = scored_states
        score_df = cd8_io.obs[score_cols]
        cd8_io.obs["cd8_state"] = score_df.idxmax(axis=1)

        state_counts = cd8_io.obs.groupby(["response", "cd8_state"]).size().reset_index(name="count")
        state_fracs = state_counts.copy()
        totals = state_counts.groupby("response")["count"].transform("sum")
        state_fracs["fraction"] = state_fracs["count"] / totals
        state_fracs.to_csv(TAB / "cd8_state_fractions.csv", index=False)

        # Fig: CD8 state composition R vs NR
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Panel A: Stacked bar
        ax = axes[0]
        pivot = state_fracs.pivot(index="response", columns="cd8_state", values="fraction").fillna(0)
        pivot.plot(kind="bar", stacked=True, ax=ax, colormap="Set2", width=0.5)
        ax.set_ylabel("Fraction")
        ax.set_title("A. CD8 State Composition", fontweight="bold")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
        ax.legend(fontsize=7, bbox_to_anchor=(1.05, 1), loc="upper left")

        # Panel B: Score violin per state
        ax = axes[1]
        melt = cd8_io.obs[["response"] + scored_states].melt(
            id_vars="response", var_name="State", value_name="Score")
        order = scored_states
        sns.boxplot(data=melt, x="State", y="Score", hue="response",
                    ax=ax, palette={"Responder": "#2ecc71", "Non-Responder": "#e74c3c"},
                    fliersize=0.5, linewidth=0.8)
        ax.set_title("B. CD8 State Scores: R vs NR", fontweight="bold")
        ax.tick_params(axis="x", rotation=30, labelsize=8)
        ax.legend(fontsize=8)

        plt.suptitle("CD8 T Cell State Analysis", fontweight="bold", fontsize=13)
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"04_cd8_states.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

        # Statistical tests
        test_results = []
        for state in scored_states:
            r_vals = cd8_io.obs.loc[cd8_io.obs["response"] == "Responder", state]
            nr_vals = cd8_io.obs.loc[cd8_io.obs["response"] == "Non-Responder", state]
            stat, pval = mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
            test_results.append({
                "state": state, "mean_R": r_vals.mean(), "mean_NR": nr_vals.mean(),
                "effect": r_vals.mean() - nr_vals.mean(), "pvalue": pval,
            })
        test_df = pd.DataFrame(test_results).sort_values("pvalue")
        test_df["padj"] = np.minimum(test_df["pvalue"] * len(test_df), 1.0)
        test_df.to_csv(TAB / "cd8_state_tests.csv", index=False)
        print(f"  CD8 states scored: {scored_states}")

    return cd8_io


# ═══════════════════════════════════════════════════════════════════════
# 4.  CLONOTYPE DIVERSITY
# ═══════════════════════════════════════════════════════════════════════

def clonotype_diversity(adata):
    """Analyze TCR clonotype diversity metrics."""
    print("\n[4] Clonotype Diversity Analysis...")
    tcr_path = PROC / "GSE206325_tcr.csv"
    if not tcr_path.exists():
        print("  [SKIP] TCR data not found")
        return None

    tcr = pd.read_csv(tcr_path)
    print(f"  TCR entries: {len(tcr):,}, columns: {list(tcr.columns)[:10]}")

    # Identify clonotype and barcode columns
    clone_col = None
    bc_col = None
    for c in tcr.columns:
        cl = c.lower()
        if "clone" in cl or "cdr3" in cl:
            clone_col = c
        if "barcode" in cl or "cell" in cl:
            bc_col = c
    if clone_col is None or bc_col is None:
        print("  [SKIP] Cannot identify clonotype/barcode columns")
        return None

    # Merge with atlas metadata
    if "response" in adata.obs.columns:
        shared = set(tcr[bc_col]) & set(adata.obs_names)
        if len(shared) < 10:
            print(f"  [SKIP] Only {len(shared)} shared barcodes")
            return None

        tcr_meta = tcr.set_index(bc_col)
        tcr_meta = tcr_meta.loc[tcr_meta.index.isin(adata.obs_names)]
        tcr_meta["response"] = adata.obs.loc[tcr_meta.index, "response"].values
        tcr_meta["patient_ID"] = adata.obs.loc[tcr_meta.index, "patient_ID"].values

        # Diversity metrics per patient
        diversity_results = []
        for pid in tcr_meta["patient_ID"].dropna().unique():
            p_tcr = tcr_meta[tcr_meta["patient_ID"] == pid]
            response = p_tcr["response"].iloc[0] if p_tcr["response"].notna().any() else None
            if response is None:
                continue

            clones = p_tcr[clone_col].dropna()
            if len(clones) < 10:
                continue

            clone_counts = clones.value_counts()
            n_clones = len(clone_counts)
            n_cells = clone_counts.sum()

            # Shannon diversity
            p_i = clone_counts / n_cells
            shannon = entropy(p_i)

            # Simpson diversity (1 - D)
            simpson = 1 - (p_i ** 2).sum()

            # Clonality (1 - normalized Shannon)
            max_entropy = np.log(n_clones) if n_clones > 1 else 1
            clonality = 1 - (shannon / max_entropy) if max_entropy > 0 else 0

            # Top clone fraction
            top1_frac = clone_counts.iloc[0] / n_cells

            # Expanded clones (>1 cell)
            n_expanded = (clone_counts > 1).sum()
            frac_expanded = n_expanded / n_clones

            diversity_results.append({
                "patient_ID": pid, "response": response,
                "n_cells_tcr": n_cells, "n_unique_clones": n_clones,
                "shannon_diversity": shannon, "simpson_diversity": simpson,
                "clonality": clonality, "top1_clone_frac": top1_frac,
                "n_expanded_clones": n_expanded, "frac_expanded": frac_expanded,
            })

        div_df = pd.DataFrame(diversity_results)
        if len(div_df) < 4:
            print("  [SKIP] Too few patients for diversity analysis")
            return None

        div_df.to_csv(TAB / "tcr_diversity_metrics.csv", index=False)

        # Fig: Diversity metrics R vs NR
        metrics = ["shannon_diversity", "clonality", "top1_clone_frac", "frac_expanded"]
        metric_labels = ["Shannon Diversity", "Clonality", "Top Clone Fraction", "Fraction Expanded"]
        fig, axes = plt.subplots(1, 4, figsize=(18, 5))
        for ax, metric, label in zip(axes, metrics, metric_labels):
            r_vals = div_df.loc[div_df["response"] == "Responder", metric]
            nr_vals = div_df.loc[div_df["response"] == "Non-Responder", metric]
            bp = ax.boxplot([r_vals, nr_vals], labels=["R", "NR"], patch_artist=True,
                           widths=0.5, showfliers=True)
            bp["boxes"][0].set_facecolor("#2ecc71")
            bp["boxes"][1].set_facecolor("#e74c3c")
            if len(r_vals) >= 3 and len(nr_vals) >= 3:
                stat, pval = mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
                ax.set_title(f"{label}\np={pval:.2e}", fontsize=9)
            else:
                ax.set_title(label, fontsize=9)
            # Overlay individual points
            for j, (vals, x) in enumerate([(r_vals, 1), (nr_vals, 2)]):
                ax.scatter(np.random.normal(x, 0.05, len(vals)), vals,
                          alpha=0.6, s=20, c=["#2ecc71", "#e74c3c"][j], zorder=5)
        plt.suptitle("TCR Clonotype Diversity: R vs NR", fontweight="bold", fontsize=13)
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"04_tcr_diversity.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

        return div_df
    return None


# ═══════════════════════════════════════════════════════════════════════
# 5.  IMMUNE CELL RATIO BIOMARKERS
# ═══════════════════════════════════════════════════════════════════════

def immune_ratio_biomarkers(adata):
    """Compute and test immune cell ratio biomarkers."""
    print("\n[5] Immune Cell Ratio Biomarkers...")

    io_mask = adata.obs["response"].isin(["Responder", "Non-Responder"])
    io_cells = adata[io_mask]

    if "patient_ID" not in io_cells.obs.columns:
        print("  [SKIP] No patient_ID column")
        return None

    # Define ratios to test
    ratio_defs = [
        ("CD8_Treg_ratio", ["CD8", "Cytotoxic", "CTL"], ["Treg", "Foxp3"]),
        ("M1_M2_ratio", ["M1", "Macro_M1"], ["M2", "Macro_M2", "CD163"]),
        ("Effector_Exhausted_ratio", ["Effector", "Cytotoxic"], ["Exhausted", "Terminal"]),
        ("T_Tumor_ratio", ["CD8", "CD4", "T_cell", "NK"], ["Tumor", "Hepato", "Malig"]),
        ("DC_Macro_ratio", ["DC", "cDC", "mregDC"], ["Macro", "Monocyte"]),
    ]

    all_cts = io_cells.obs["cell_type"].unique()

    results = []
    for ratio_name, numer_keys, denom_keys in ratio_defs:
        numer_cts = [ct for ct in all_cts if any(k.lower() in ct.lower() for k in numer_keys)]
        denom_cts = [ct for ct in all_cts if any(k.lower() in ct.lower() for k in denom_keys)]
        if not numer_cts or not denom_cts:
            continue

        for pid in io_cells.obs["patient_ID"].dropna().unique():
            p_cells = io_cells[io_cells.obs["patient_ID"] == pid]
            n_numer = p_cells.obs["cell_type"].isin(numer_cts).sum()
            n_denom = p_cells.obs["cell_type"].isin(denom_cts).sum()
            ratio_val = np.log2((n_numer + 1) / (n_denom + 1))
            results.append({
                "patient_ID": pid,
                "response": p_cells.obs["response"].iloc[0],
                "ratio": ratio_name,
                "value": ratio_val,
                "n_numerator": n_numer,
                "n_denominator": n_denom,
            })

    if not results:
        return None

    ratio_df = pd.DataFrame(results)
    ratio_df.to_csv(TAB / "immune_ratios.csv", index=False)

    # Fig: Ratio biomarkers boxplot
    unique_ratios = ratio_df["ratio"].unique()
    n_plots = len(unique_ratios)
    if n_plots == 0:
        return None

    fig, axes = plt.subplots(1, min(n_plots, 5), figsize=(4*min(n_plots, 5), 5))
    if n_plots == 1:
        axes = [axes]
    for i, ratio_name in enumerate(unique_ratios[:5]):
        ax = axes[i]
        sub = ratio_df[ratio_df["ratio"] == ratio_name]
        r_vals = sub.loc[sub["response"] == "Responder", "value"]
        nr_vals = sub.loc[sub["response"] == "Non-Responder", "value"]

        bp = ax.boxplot([r_vals, nr_vals], labels=["R", "NR"], patch_artist=True,
                       widths=0.5, showfliers=True)
        bp["boxes"][0].set_facecolor("#2ecc71")
        bp["boxes"][1].set_facecolor("#e74c3c")

        # Overlay points
        for j, (vals, x) in enumerate([(r_vals, 1), (nr_vals, 2)]):
            ax.scatter(np.random.normal(x, 0.05, len(vals)), vals,
                      alpha=0.6, s=20, c=["#2ecc71", "#e74c3c"][j], zorder=5)

        if len(r_vals) >= 3 and len(nr_vals) >= 3:
            stat, pval = mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
            ax.set_title(f"{ratio_name.replace('_', ' ')}\np={pval:.2e}", fontsize=9)
        else:
            ax.set_title(ratio_name.replace("_", " "), fontsize=9)
        ax.set_ylabel("log2 Ratio")

    plt.suptitle("Immune Cell Ratio Biomarkers: R vs NR", fontweight="bold", fontsize=13)
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"04_immune_ratios.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    # AUC for each ratio
    auc_results = []
    for ratio_name in unique_ratios:
        sub = ratio_df[ratio_df["ratio"] == ratio_name].dropna()
        if sub["response"].nunique() < 2:
            continue
        y = (sub["response"] == "Responder").astype(int).values
        vals = sub["value"].values
        if len(np.unique(y)) < 2:
            continue
        auc = roc_auc_score(y, vals)
        auc_results.append({"ratio": ratio_name, "AUC": auc})

    if auc_results:
        auc_df = pd.DataFrame(auc_results).sort_values("AUC", ascending=False)
        auc_df.to_csv(TAB / "immune_ratio_aucs.csv", index=False)
        print(f"  Ratio AUCs:")
        for _, row in auc_df.iterrows():
            print(f"    {row['ratio']}: AUC={row['AUC']:.3f}")

    return ratio_df


# ═══════════════════════════════════════════════════════════════════════
# 6.  PATIENT IMMUNE SUBTYPING
# ═══════════════════════════════════════════════════════════════════════

def patient_immune_subtyping(adata):
    """Classify patients into immune subtypes based on TME composition."""
    print("\n[6] Patient Immune Subtyping...")

    io_mask = adata.obs["response"].isin(["Responder", "Non-Responder"])
    io_cells = adata[io_mask]

    if "patient_ID" not in io_cells.obs.columns:
        return None

    # Build patient feature matrix
    sig_cols = [c for c in io_cells.obs.columns if c in IO_SIGNATURES]
    features = {}
    for pid in io_cells.obs["patient_ID"].dropna().unique():
        p_cells = io_cells[io_cells.obs["patient_ID"] == pid]
        if len(p_cells) < 50:
            continue
        row = {"patient_ID": pid, "response": p_cells.obs["response"].iloc[0]}
        for sig in sig_cols:
            row[sig] = p_cells.obs[sig].mean()
        # Add cell type fractions
        ct_fracs = p_cells.obs["cell_type"].value_counts(normalize=True)
        for ct, f in ct_fracs.items():
            row[f"frac_{ct}"] = f
        features[pid] = row

    patient_df = pd.DataFrame(features.values()).fillna(0)
    if len(patient_df) < 6:
        print("  [SKIP] Too few patients")
        return None

    # Hierarchical clustering of patients by signature profiles
    feat_cols = sig_cols + [c for c in patient_df.columns if c.startswith("frac_")]
    X = patient_df[feat_cols].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
    Z = linkage(X_scaled, method="ward")

    # Cut into 3 clusters
    patient_df["immune_subtype"] = fcluster(Z, t=3, criterion="maxclust")
    patient_df["immune_subtype"] = "Subtype_" + patient_df["immune_subtype"].astype(str)
    patient_df.to_csv(TAB / "patient_immune_subtypes.csv", index=False)

    # Fig: Patient heatmap with response and subtype annotations
    fig_width = max(10, len(patient_df) * 0.3)
    plot_data = patient_df.set_index("patient_ID")[sig_cols].T

    # Color bars for response and subtype
    col_colors_resp = patient_df.set_index("patient_ID")["response"].map(
        {"Responder": "#2ecc71", "Non-Responder": "#e74c3c"})
    subtype_colors = {"Subtype_1": "#3498db", "Subtype_2": "#f39c12", "Subtype_3": "#9b59b6"}
    col_colors_sub = patient_df.set_index("patient_ID")["immune_subtype"].map(subtype_colors)

    col_colors_df = pd.DataFrame({
        "Response": col_colors_resp,
        "Subtype": col_colors_sub,
    })

    g = sns.clustermap(plot_data, cmap="RdBu_r", center=0, figsize=(fig_width, 8),
                       method="ward", z_score=1, linewidths=0.3,
                       col_colors=col_colors_df,
                       xticklabels=True, yticklabels=True,
                       cbar_kws={"label": "Z-score"})
    g.fig.suptitle("Patient Immune Subtyping by IO Signatures", fontweight="bold", y=1.02)
    for ext in ["png", "pdf"]:
        g.savefig(FIG / f"04_patient_subtypes.{ext}", dpi=200, bbox_inches="tight")
    plt.close("all")

    # Response rates per subtype
    subtype_resp = patient_df.groupby("immune_subtype")["response"].value_counts().unstack(fill_value=0)
    subtype_resp.to_csv(TAB / "subtype_response_rates.csv")
    print("  Subtype response rates:")
    print(subtype_resp.to_string())

    return patient_df


# ═══════════════════════════════════════════════════════════════════════
# 7.  CROSS-DATASET SIGNATURE REPRODUCIBILITY
# ═══════════════════════════════════════════════════════════════════════

def cross_dataset_reproducibility(adata):
    """Test if key signatures are reproduced across datasets."""
    print("\n[7] Cross-Dataset Signature Reproducibility...")

    sig_cols = [c for c in adata.obs.columns if c in IO_SIGNATURES]
    if len(sig_cols) < 3:
        print("  [SKIP] too few signatures scored")
        return None

    datasets = adata.obs["dataset"].unique()
    if len(datasets) < 2:
        print("  [SKIP] Single dataset")
        return None

    # Compare mean signatures across datasets
    results = []
    for ds in sorted(datasets):
        ds_cells = adata[adata.obs["dataset"] == ds]
        row = {"dataset": ds, "n_cells": ds_cells.shape[0]}
        for sig in sig_cols:
            row[f"mean_{sig}"] = ds_cells.obs[sig].mean()
        results.append(row)

    cross_df = pd.DataFrame(results)
    cross_df.to_csv(TAB / "cross_dataset_signatures.csv", index=False)

    # Fig: Heatmap of signatures across datasets
    mean_cols = [f"mean_{s}" for s in sig_cols if f"mean_{s}" in cross_df.columns]
    plot_df = cross_df.set_index("dataset")[mean_cols]
    plot_df.columns = [c.replace("mean_", "") for c in plot_df.columns]

    fig, ax = plt.subplots(figsize=(12, max(4, len(datasets)*0.8)))
    # Z-score per signature for visualization
    plot_z = plot_df.sub(plot_df.mean(axis=0), axis=1).div(plot_df.std(axis=0) + 1e-10, axis=1)
    sns.heatmap(plot_z, cmap="RdBu_r", center=0, ax=ax, annot=True, fmt=".2f",
                linewidths=0.5, cbar_kws={"label": "Z-score"})
    ax.set_title("IO Signature Profiles Across Datasets", fontweight="bold")
    ax.set_xlabel("Signature")
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"04_cross_dataset_signatures.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    return cross_df


# ═══════════════════════════════════════════════════════════════════════
# 8.  MULTI-PANEL SUMMARY FIGURE
# ═══════════════════════════════════════════════════════════════════════

def summary_figure(adata, pw_df, lr_df):
    """Generate a multi-panel summary figure."""
    print("\n[8] Generating multi-panel summary figure...")

    io_mask = adata.obs["response"].isin(["Responder", "Non-Responder"])
    io_cells = adata[io_mask]

    fig, axes = plt.subplots(2, 3, figsize=(20, 13))

    # Panel A: UMAP by response
    ax = axes[0, 0]
    resp_colors = {"Responder": "#2ecc71", "Non-Responder": "#e74c3c"}
    for resp, color in resp_colors.items():
        mask = io_cells.obs["response"] == resp
        if mask.sum() > 0:
            ax.scatter(io_cells.obsm["X_umap"][mask, 0], io_cells.obsm["X_umap"][mask, 1],
                       s=0.3, alpha=0.1, c=color, label=resp, rasterized=True)
    ax.legend(fontsize=8, markerscale=10)
    ax.set_title("A. UMAP: IO Response", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    # Panel B: Top IO signatures effect (forest plot style)
    ax = axes[0, 1]
    sig_cols = [c for c in io_cells.obs.columns if c in IO_SIGNATURES]
    effects = []
    for sig in sig_cols:
        r_mean = io_cells.obs.loc[io_cells.obs["response"] == "Responder", sig].mean()
        nr_mean = io_cells.obs.loc[io_cells.obs["response"] == "Non-Responder", sig].mean()
        effects.append({"signature": sig, "effect": r_mean - nr_mean})
    eff_df = pd.DataFrame(effects).sort_values("effect")
    colors = ["#2ecc71" if e > 0 else "#e74c3c" for e in eff_df["effect"]]
    ax.barh(range(len(eff_df)), eff_df["effect"], color=colors, height=0.5)
    ax.set_yticks(range(len(eff_df)))
    ax.set_yticklabels(eff_df["signature"].str.replace("_", " "), fontsize=7)
    ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Effect Size (R - NR)")
    ax.set_title("B. IO Signature Differences", fontweight="bold")

    # Panel C: Top pathways
    ax = axes[0, 2]
    if pw_df is not None and len(pw_df) > 0:
        top_pw = pw_df.nlargest(8, "effect").append(pw_df.nsmallest(4, "effect")).sort_values("effect")
        colors_pw = ["#2ecc71" if e > 0 else "#e74c3c" for e in top_pw["effect"]]
        ax.barh(range(len(top_pw)), top_pw["effect"], color=colors_pw, height=0.5)
        ax.set_yticks(range(len(top_pw)))
        ax.set_yticklabels(top_pw["pathway"].str.replace("_", " "), fontsize=7)
        ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Effect Size (R - NR)")
    ax.set_title("C. Top Pathway Differences", fontweight="bold")

    # Panel D: Immune infiltration pattern
    ax = axes[1, 0]
    if "Immune_Infiltration" in io_cells.obs.columns:
        inf_resp = pd.crosstab(io_cells.obs["Immune_Infiltration"],
                               io_cells.obs["response"], normalize="index")
        if "Responder" in inf_resp.columns:
            inf_resp = inf_resp.sort_values("Responder")
            inf_resp.plot(kind="barh", stacked=True, ax=ax,
                         color=["#e74c3c", "#2ecc71"], width=0.5)
            ax.set_xlabel("Fraction")
            ax.legend(fontsize=8)
    ax.set_title("D. Response by Immune Infiltration", fontweight="bold")

    # Panel E: Checkpoint expression heatmap
    ax = axes[1, 1]
    checkpoints = ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "CD274"]
    genes_avail = set(io_cells.raw.var_names) if io_cells.raw else set(io_cells.var_names)
    cp_present = [g for g in checkpoints if g in genes_avail]
    if cp_present:
        adata_raw = io_cells.raw.to_adata() if io_cells.raw else io_cells
        cp_data = {}
        for resp in ["Responder", "Non-Responder"]:
            mask = io_cells.obs["response"] == resp
            sub = adata_raw[mask]
            X = sub[:, cp_present].X
            if hasattr(X, "toarray"):
                X = X.toarray()
            cp_data[resp] = X.mean(axis=0)
        cp_df = pd.DataFrame(cp_data, index=cp_present)
        sns.heatmap(cp_df, cmap="YlOrRd", ax=ax, annot=True, fmt=".3f", linewidths=0.5)
    ax.set_title("E. Checkpoint Expression", fontweight="bold")

    # Panel F: Response prediction (placeholder text)
    ax = axes[1, 2]
    # Load LASSO results if available
    lasso_path = PROJECT / "out" / "03_validation" / "tables" / "lasso_coefficients.csv"
    if lasso_path.exists():
        coefs = pd.read_csv(lasso_path)
        top = coefs[coefs["abs_coef"] > 0.01].head(10)
        colors_c = ["#2ecc71" if c > 0 else "#e74c3c" for c in top["coefficient"]]
        ax.barh(range(len(top)), top["coefficient"], color=colors_c, height=0.5)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top["feature"].str.replace("frac_", "% "), fontsize=7)
        ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
    ax.set_xlabel("LASSO Coefficient")
    ax.set_title("F. Top Response Predictors", fontweight="bold")

    plt.suptitle("HCC Immunotherapy Response — Complementary Analyses Summary",
                 fontweight="bold", fontsize=14)
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"04_summary_multipanel.{ext}", dpi=200, bbox_inches="tight")
    plt.close()


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("HCC IO Response — Complementary Downstream Analyses")
    print("=" * 70)

    # ── Load integrated atlas ──────────────────────────────────────
    h5_path = PROC / "integrated_io_atlas.h5ad"
    if not h5_path.exists():
        print("[ERROR] Integrated atlas not found.")
        return
    adata = sc.read_h5ad(h5_path)
    print(f"  {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

    # ── Recover metadata if needed ─────────────────────────────────
    adata = recover_metadata(adata)

    # ── Ensure signatures are scored ───────────────────────────────
    sig_cols = [c for c in adata.obs.columns if c in IO_SIGNATURES]
    if len(sig_cols) < 5:
        print("  Scoring IO signatures...")
        genes_available = set(adata.raw.var_names) if adata.raw else set(adata.var_names)
        for name, genes in IO_SIGNATURES.items():
            present = [g for g in genes if g in genes_available]
            if len(present) >= 2:
                sc.tl.score_genes(adata, present, score_name=name, use_raw=True)

    # ── Ensure cell types annotated ────────────────────────────────
    if "cell_type" not in adata.obs.columns:
        print("  [WARN] No cell_type column — using cell_type_original if available")
        if "cell_type_original" in adata.obs.columns:
            adata.obs["cell_type"] = adata.obs["cell_type_original"]
        else:
            print("  [ERROR] No cell type annotations. Run 02a first.")
            return

    # ═══════════════════════════════════════════════════════════════
    # RUN ALL ANALYSES
    # ═══════════════════════════════════════════════════════════════
    pw_df, scored_pw = pathway_enrichment_analysis(adata)
    lr_df = ligand_receptor_analysis(adata)
    cd8_state_analysis(adata)
    div_df = clonotype_diversity(adata)
    immune_ratio_biomarkers(adata)
    patient_immune_subtyping(adata)
    cross_dataset_reproducibility(adata)
    summary_figure(adata, pw_df, lr_df)

    # ═══════════════════════════════════════════════════════════════
    # DONE
    # ═══════════════════════════════════════════════════════════════
    n_figs = len(list(FIG.glob("*.png")))
    n_tabs = len(list(TAB.glob("*.csv")))
    print(f"\n{'='*70}")
    print(f"[DONE] Complementary Analyses")
    print(f"  Figures: {n_figs} in {FIG}")
    for f in sorted(FIG.glob("*.png")):
        print(f"    {f.name}")
    print(f"  Tables: {n_tabs} in {TAB}")
    for f in sorted(TAB.glob("*.csv")):
        print(f"    {f.name}")
    print(f"={'='*70}")


if __name__ == "__main__":
    main()
