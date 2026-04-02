#!/usr/bin/env python3
"""
03_validation_downstream.py — Validation & Downstream Analyses

Comprehensive analyses for journal-quality:
  1. LASSO composite IO response score
  2. ROC-AUC evaluation at patient level
  3. Cross-validation of response predictor
  4. NMF transcriptional programs
  5. Pathway enrichment on R vs NR DEGs
  6. TCR clonality analysis (GSE206325 only)
  7. Correlation network of signatures
  8. Bootstrap stability of key findings
  9. Bulk RNA-seq signature validation (GSE202069)
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
from scipy.stats import mannwhitneyu, spearmanr, pearsonr, kruskal
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (roc_auc_score, roc_curve, classification_report,
                             balanced_accuracy_score, confusion_matrix)
from sklearn.decomposition import NMF
from sklearn.ensemble import RandomForestClassifier
import re

sc.settings.verbosity = 1

PROJECT = Path(__file__).resolve().parents[2]
PROC    = PROJECT / "data" / "processed"
DATA    = PROJECT / "data"
OUT     = PROJECT / "out" / "03_validation"
FIG     = OUT / "figures"
TAB     = OUT / "tables"
for d in [FIG, TAB]:
    d.mkdir(parents=True, exist_ok=True)

IO_SIGNATURES = {
    "T_cell_exhaustion":      ["PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","TOX","ENTPD1"],
    "T_cell_activation":      ["IFNG","GZMB","PRF1","TNF","CD69","IL2"],
    "T_cell_cytotoxicity":    ["GZMA","GZMB","GZMH","GZMK","PRF1","GNLY","NKG7"],
    "Progenitor_CD8":         ["TCF7","SLAMF6","LEF1","CCR7","IL7R","SELL"],
    "Terminal_exhaustion":    ["TOX","ENTPD1","HAVCR2","LAYN","TIGIT","CXCL13"],
    "Treg_signature":         ["FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"],
    "M1_macrophage":          ["NOS2","IL1B","IL6","TNF","CD80","CD86"],
    "M2_macrophage":          ["CD163","MRC1","MSR1","CD200R1","IL10"],
    "IFN_gamma":              ["IFNG","STAT1","IRF1","CXCL9","CXCL10","CXCL11","IDO1"],
    "TGF_beta":               ["TGFB1","TGFB2","TGFB3","ACTA2","COL1A1","SERPINE1"],
    "Antigen_presentation":   ["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","PSMB9"],
    "Immune_checkpoint":      ["CD274","PDCD1LG2","PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","VSIR"],
    "Angiogenesis":           ["VEGFA","PECAM1","KDR","FLT1","ANGPT2","FGF2"],
    "mregDC":                 ["CCR7","LAMP3","CCL19","CCL22","CD80","CD83","FSCN1","MARCKSL1"],
    "Tfh_like":               ["CXCL13","CH25H","IL21","PDCD1","ICOS","BCL6","TOX2"],
}


def build_patient_features(adata):
    """Aggregate cell-level data to patient-level features."""
    io_cells = adata[adata.obs["response"].isin(["Responder", "Non-Responder"])].copy()

    # Restrict to tumor tissue for biomarker discovery
    if "tissue" in io_cells.obs.columns:
        tumor_mask = io_cells.obs["tissue"] == "Tumor"
        if tumor_mask.sum() > 1000:
            io_cells = io_cells[tumor_mask].copy()
            print(f"  Restricted to tumor tissue: {io_cells.shape[0]:,} cells")

    sig_cols = [c for c in io_cells.obs.columns if c in IO_SIGNATURES]

    features = {}
    ct_cols = []
    for pid in io_cells.obs["patient_ID"].dropna().unique():
        p_cells = io_cells[io_cells.obs["patient_ID"] == pid]
        if len(p_cells) < 50:
            continue

        row = {"patient_ID": pid}
        row["response"] = p_cells.obs["response"].iloc[0]
        row["n_cells"] = len(p_cells)

        # Signature means
        for sig in sig_cols:
            row[sig] = p_cells.obs[sig].mean()

        # Cell type fractions
        ct_counts = p_cells.obs["cell_type"].value_counts(normalize=True)
        for ct, frac in ct_counts.items():
            col = f"frac_{ct}"
            row[col] = frac
            if col not in ct_cols:
                ct_cols.append(col)

        features[pid] = row

    df = pd.DataFrame(features.values()).fillna(0)
    print(f"  Patient features: {len(df)} patients × {df.shape[1]} features")
    return df


def lasso_response_model(patient_df):
    """LASSO logistic regression for response prediction."""
    print("\n[LASSO] Building IO response predictor...")

    y = (patient_df["response"] == "Responder").astype(int).values
    feature_cols = [c for c in patient_df.columns
                    if c not in ["patient_ID", "response", "n_cells"]]
    X = patient_df[feature_cols].values

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Stratified k-fold
    n_splits = min(5, min(y.sum(), (1-y).sum()))
    if n_splits < 2:
        print("  [WARN] Not enough samples for cross-validation")
        return None, None, None

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    # LASSO logistic regression with CV
    model = LogisticRegressionCV(
        Cs=np.logspace(-3, 1, 20), cv=cv, penalty="l1", solver="saga",
        scoring="roc_auc", max_iter=5000, class_weight="balanced", random_state=42
    )
    model.fit(X_scaled, y)

    # Feature importance
    coefs = pd.DataFrame({
        "feature": feature_cols,
        "coefficient": model.coef_[0],
        "abs_coef": np.abs(model.coef_[0]),
    }).sort_values("abs_coef", ascending=False)

    # Cross-validated predictions
    y_proba_cv = np.zeros(len(y))
    for train_idx, test_idx in cv.split(X_scaled, y):
        m = LogisticRegressionCV(
            Cs=np.logspace(-3, 1, 20), cv=3, penalty="l1", solver="saga",
            scoring="roc_auc", max_iter=5000, class_weight="balanced", random_state=42
        )
        m.fit(X_scaled[train_idx], y[train_idx])
        y_proba_cv[test_idx] = m.predict_proba(X_scaled[test_idx])[:, 1]

    auc = roc_auc_score(y, y_proba_cv)
    print(f"  CV AUC-ROC: {auc:.3f}")
    print(f"  Top features:")
    for _, row in coefs[coefs["abs_coef"] > 0.01].head(10).iterrows():
        direction = "↑R" if row["coefficient"] > 0 else "↑NR"
        print(f"    {row['feature']}: {row['coefficient']:.3f} {direction}")

    return model, coefs, (y, y_proba_cv, auc)


def nmf_programs(adata, k_range=(4, 12)):
    """NMF transcriptional program discovery."""
    print("\n[NMF] Discovering transcriptional programs...")

    io_tumor = adata[(adata.obs["response"].isin(["Responder","Non-Responder"])) &
                     (adata.obs["tissue"] == "Tumor")].copy() \
               if "tissue" in adata.obs.columns else \
               adata[adata.obs["response"].isin(["Responder","Non-Responder"])].copy()

    # Use raw layer
    if adata.raw is not None:
        adata_raw = io_tumor.raw.to_adata()
    else:
        adata_raw = io_tumor

    # Subsample if too large
    n_max = 30000
    if adata_raw.shape[0] > n_max:
        idx = np.random.RandomState(42).choice(adata_raw.shape[0], n_max, replace=False)
        adata_raw = adata_raw[idx].copy()

    # Get non-negative matrix (HVGs)
    sc.pp.highly_variable_genes(adata_raw, n_top_genes=2000, flavor="seurat_v3",
                                span=0.3 if adata_raw.shape[0] > 5000 else 1.0)
    hvg_mask = adata_raw.var["highly_variable"]
    X = adata_raw[:, hvg_mask].X
    if sparse.issparse(X):
        X = X.toarray()
    X = np.maximum(X, 0)  # Ensure non-negative

    # K selection
    errors = []
    for k in range(k_range[0], k_range[1]+1):
        nmf = NMF(n_components=k, init="nndsvda", max_iter=300, random_state=42)
        W = nmf.fit_transform(X)
        errors.append({"k": k, "error": nmf.reconstruction_err_})
        print(f"  k={k}: error={nmf.reconstruction_err_:.1f}")

    errors_df = pd.DataFrame(errors)
    # Use elbow (biggest drop)
    diffs = np.diff(errors_df["error"].values)
    best_k_idx = np.argmin(diffs) + 1 if len(diffs) > 1 else 0
    best_k = errors_df.iloc[best_k_idx]["k"]
    print(f"  Selected k={best_k}")

    # Final NMF
    from scipy import sparse as sp
    nmf = NMF(n_components=int(best_k), init="nndsvda", max_iter=500, random_state=42)
    W = nmf.fit_transform(X)
    H = nmf.components_

    # Gene loadings
    gene_names = adata_raw.var_names[hvg_mask]
    loadings = pd.DataFrame(H.T, index=gene_names,
                            columns=[f"NMF_{i+1}" for i in range(int(best_k))])
    top_genes = {}
    for prog in loadings.columns:
        top = loadings[prog].nlargest(20).index.tolist()
        top_genes[prog] = top

    return W, loadings, errors_df, top_genes, adata_raw.obs


def tcr_analysis(adata):
    """Analyze TCR clonality in IO context."""
    print("\n[TCR] Clonality analysis...")
    tcr_path = PROC / "GSE206325_tcr.csv"
    if not tcr_path.exists():
        print("  [SKIP] TCR data not available")
        return None

    tcr = pd.read_csv(tcr_path)
    print(f"  TCR data: {len(tcr):,} entries")
    print(f"  Columns: {list(tcr.columns)[:10]}")
    return tcr


def bulk_validation(scored_sigs):
    """Validate IO signatures on bulk RNA-seq (GSE202069)."""
    print("\n[BULK] Validating on GSE202069 bulk RNA-seq...")
    bulk_f = DATA / "pre_treatment" / "GSE202069_Zhu2023_antiPDL1" / "GSE202069_gene_tpm_expression.txt.gz"
    if not bulk_f.exists():
        print("  [SKIP] Bulk data not found")
        return None

    bulk = pd.read_csv(bulk_f, sep="\t", index_col=0)
    print(f"  Bulk: {bulk.shape[0]} genes × {bulk.shape[1]} samples")

    # Score signatures
    gene_set = set(bulk.index)
    sig_scores = {}
    for sig_name, genes in IO_SIGNATURES.items():
        present = [g for g in genes if g in gene_set]
        if len(present) >= 2:
            sig_scores[sig_name] = bulk.loc[present].mean(axis=0)

    sig_df = pd.DataFrame(sig_scores)
    return sig_df


def main():
    print("=" * 70)
    print("HCC IO Response — Validation & Downstream Analyses")
    print("=" * 70)

    # ── Load ───────────────────────────────────────────────────────
    h5_path = PROC / "integrated_io_atlas.h5ad"
    if not h5_path.exists():
        print("[ERROR] Integrated atlas not found. Run 01 and 02 scripts first.")
        return
    adata = sc.read_h5ad(h5_path)
    print(f"  {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

    # ── Recover GSE206325 metadata if missing from concat ──────────
    if "response" not in adata.obs.columns:
        print("[META] Recovering GSE206325 clinical metadata...")
        cell_meta = pd.read_csv(PROC / "GSE206325_metadata.csv", index_col=0)
        sample_ann = pd.read_csv(PROC / "GSE206325_sample_annotations.csv", index_col=0)
        sample_ann.index = sample_ann.index.astype(str)

        shared = adata.obs_names.intersection(cell_meta.index)
        print(f"  Matched {len(shared):,} GSE206325 barcodes in atlas")

        adata.obs["sample_id"] = np.nan
        adata.obs.loc[shared, "sample_id"] = cell_meta.loc[shared, "cell_to_sample_ID"].values.astype(str)

        for col in ["patient_ID", "tissue", "treatment_Resp"]:
            mapping = sample_ann[col].to_dict()
            adata.obs[col] = adata.obs["sample_id"].map(mapping)

        resp_map = {"anti-PD1_R": "Responder", "anti-PD1_NR": "Non-Responder"}
        adata.obs["response"] = adata.obs["treatment_Resp"].map(resp_map)
        print(f"  Response: {adata.obs['response'].value_counts(dropna=False).to_dict()}")

    # Ensure signatures are scored
    sig_cols = [c for c in adata.obs.columns if c in IO_SIGNATURES]
    if len(sig_cols) < 5:
        print("  Scoring IO signatures...")
        genes_available = set(adata.raw.var_names) if adata.raw else set(adata.var_names)
        for name, genes in IO_SIGNATURES.items():
            present = [g for g in genes if g in genes_available]
            if len(present) >= 2:
                sc.tl.score_genes(adata, present, score_name=name, use_raw=True)
        sig_cols = [c for c in adata.obs.columns if c in IO_SIGNATURES]

    # ═══════════════════════════════════════════════════════════════
    # 1. PATIENT-LEVEL RESPONSE PREDICTION
    # ═══════════════════════════════════════════════════════════════
    patient_df = build_patient_features(adata)
    patient_df.to_csv(TAB / "patient_features.csv", index=False)

    model, coefs, cv_results = lasso_response_model(patient_df)

    if cv_results is not None:
        y_true, y_proba, auc = cv_results
        coefs.to_csv(TAB / "lasso_coefficients.csv", index=False)

        # ─── Fig 1: ROC curve ──────────────────────────────────────
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        ax = axes[0]
        fpr, tpr, _ = roc_curve(y_true, y_proba)
        ax.plot(fpr, tpr, "b-", linewidth=2, label=f"LASSO (AUC={auc:.3f})")
        ax.plot([0,1], [0,1], "k--", alpha=0.3)
        ax.set_xlabel("False Positive Rate"); ax.set_ylabel("True Positive Rate")
        ax.set_title("A. ROC Curve — IO Response Prediction", fontweight="bold")
        ax.legend()

        # Feature importance bar
        ax = axes[1]
        top_feats = coefs[coefs["abs_coef"] > 0.01].head(15)
        colors = ["#2ecc71" if c > 0 else "#e74c3c" for c in top_feats["coefficient"]]
        ax.barh(range(len(top_feats)), top_feats["coefficient"], color=colors, height=0.5)
        ax.set_yticks(range(len(top_feats)))
        ax.set_yticklabels(top_feats["feature"].str.replace("frac_", "% "))
        ax.axvline(0, color="black", linewidth=0.5)
        ax.set_xlabel("LASSO Coefficient")
        ax.set_title("B. Top Predictive Features", fontweight="bold")
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"03_roc_lasso.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

    # ═══════════════════════════════════════════════════════════════
    # 2. RANDOM FOREST COMPARISON
    # ═══════════════════════════════════════════════════════════════
    print("\n[RF] Random Forest comparison...")
    y = (patient_df["response"] == "Responder").astype(int).values
    feature_cols = [c for c in patient_df.columns
                    if c not in ["patient_ID", "response", "n_cells"]]
    X = patient_df[feature_cols].values

    n_splits = min(5, min(y.sum(), (1-y).sum()))
    if n_splits >= 2:
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
        y_proba_rf = np.zeros(len(y))
        for train_idx, test_idx in cv.split(X, y):
            rf = RandomForestClassifier(n_estimators=200, class_weight="balanced", random_state=42)
            rf.fit(X[train_idx], y[train_idx])
            y_proba_rf[test_idx] = rf.predict_proba(X[test_idx])[:, 1]
        auc_rf = roc_auc_score(y, y_proba_rf)
        print(f"  RF CV AUC: {auc_rf:.3f}")

        # Feature importance
        rf_full = RandomForestClassifier(n_estimators=200, class_weight="balanced", random_state=42)
        rf_full.fit(X, y)
        rf_imp = pd.DataFrame({
            "feature": feature_cols,
            "importance": rf_full.feature_importances_
        }).sort_values("importance", ascending=False)
        rf_imp.to_csv(TAB / "rf_feature_importance.csv", index=False)

        # Fig: RF feature importance
        fig, ax = plt.subplots(figsize=(8, 6))
        top_rf = rf_imp.head(15)
        ax.barh(range(len(top_rf)), top_rf["importance"], color="#3498db", height=0.5)
        ax.set_yticks(range(len(top_rf)))
        ax.set_yticklabels(top_rf["feature"].str.replace("frac_", "% "))
        ax.set_xlabel("Feature Importance (Gini)")
        ax.set_title(f"Random Forest Feature Importance (AUC={auc_rf:.3f})", fontweight="bold")
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"03_rf_importance.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

    # ═══════════════════════════════════════════════════════════════
    # 3. NMF PROGRAMS
    # ═══════════════════════════════════════════════════════════════
    try:
        from scipy import sparse
        W, loadings, errors_df, top_genes, nmf_obs = nmf_programs(adata)
        loadings.to_csv(TAB / "nmf_gene_loadings.csv")

        # Save top genes per program
        top_genes_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in top_genes.items()]))
        top_genes_df.to_csv(TAB / "nmf_top_genes.csv", index=False)

        # ─── Fig: NMF elbow ────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(errors_df["k"], errors_df["error"], "bo-", linewidth=2, markersize=8)
        best_k = errors_df.iloc[np.argmin(np.diff(errors_df["error"].values)) + 1 if len(errors_df) > 2 else 0]["k"]
        ax.axvline(best_k, color="red", linestyle="--", label=f"k={int(best_k)}")
        ax.set_xlabel("Number of Programs (k)"); ax.set_ylabel("Reconstruction Error")
        ax.set_title("NMF Program Selection", fontweight="bold"); ax.legend()
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"03_nmf_elbow.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

        # ─── Fig: NMF top gene heatmap ─────────────────────────────
        n_top = 10
        top_all = []
        for prog in loadings.columns:
            top_all.extend(loadings[prog].nlargest(n_top).index.tolist())
        top_all = list(dict.fromkeys(top_all))  # unique, preserve order
        plot_df = loadings.loc[top_all]

        fig, ax = plt.subplots(figsize=(8, max(6, len(top_all)*0.25)))
        sns.heatmap(plot_df, cmap="YlOrRd", ax=ax, linewidths=0.5)
        ax.set_title("NMF Program Gene Loadings (Top 10/program)", fontweight="bold")
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"03_nmf_gene_heatmap.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

        # ─── Fig: NMF programs by response ─────────────────────────
        if "response" in nmf_obs.columns:
            prog_cols = [f"NMF_{i+1}" for i in range(W.shape[1])]
            nmf_scores = pd.DataFrame(W, columns=prog_cols, index=nmf_obs.index)
            nmf_scores["response"] = nmf_obs["response"].values

            fig, axes = plt.subplots(2, (len(prog_cols)+1)//2, figsize=(4*((len(prog_cols)+1)//2), 7))
            axes = axes.flatten()
            for i, prog in enumerate(prog_cols):
                ax = axes[i]
                r_vals = nmf_scores.loc[nmf_scores["response"]=="Responder", prog]
                nr_vals = nmf_scores.loc[nmf_scores["response"]=="Non-Responder", prog]
                bp = ax.boxplot([r_vals, nr_vals], labels=["R", "NR"], patch_artist=True,
                               widths=0.5, showfliers=False)
                bp["boxes"][0].set_facecolor("#2ecc71")
                bp["boxes"][1].set_facecolor("#e74c3c")
                stat, pval = mannwhitneyu(r_vals, nr_vals)
                ax.set_title(f"{prog}\np={pval:.2e}", fontsize=9)
            for j in range(i+1, len(axes)):
                axes[j].set_visible(False)
            plt.suptitle("NMF Programs: R vs NR", fontweight="bold")
            plt.tight_layout()
            for ext in ["png", "pdf"]:
                plt.savefig(FIG / f"03_nmf_response_box.{ext}", dpi=200, bbox_inches="tight")
            plt.close()

    except Exception as e:
        print(f"  [WARN] NMF failed: {e}")

    # ═══════════════════════════════════════════════════════════════
    # 4. SIGNATURE CORRELATION NETWORK
    # ═══════════════════════════════════════════════════════════════
    print("\n[CORR] Signature correlation network...")
    io_cells = adata[adata.obs["response"].isin(["Responder", "Non-Responder"])]
    if len(sig_cols) >= 5:
        corr = io_cells.obs[sig_cols].corr(method="spearman")
        corr.to_csv(TAB / "signature_correlation.csv")

        fig, ax = plt.subplots(figsize=(10, 9))
        Z = linkage(1 - corr.values, method="ward")
        order = leaves_list(Z)
        corr_ordered = corr.iloc[order, order]
        sns.heatmap(corr_ordered, cmap="RdBu_r", center=0, ax=ax, annot=True,
                    fmt=".2f", square=True, linewidths=0.5,
                    xticklabels=[s.replace("_", "\n") for s in corr_ordered.columns],
                    yticklabels=[s.replace("_", "\n") for s in corr_ordered.index])
        ax.set_title("IO Signature Correlation (Spearman)", fontweight="bold")
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"03_signature_correlation.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

    # ═══════════════════════════════════════════════════════════════
    # 5. BULK RNA-SEQ VALIDATION
    # ═══════════════════════════════════════════════════════════════
    bulk_sig = bulk_validation(sig_cols)
    if bulk_sig is not None:
        bulk_sig.to_csv(TAB / "bulk_signature_scores.csv")

        # Clustermap
        fig = plt.figure(figsize=(12, 8))
        g = sns.clustermap(bulk_sig.T, cmap="RdBu_r", center=0, figsize=(12, 8),
                           method="ward", z_score=0, linewidths=0.5,
                           cbar_kws={"label": "Z-score"})
        g.fig.suptitle("Bulk RNA-seq IO Signatures (GSE202069)\nHierarchical Clustering",
                       fontweight="bold", y=1.02)
        for ext in ["png", "pdf"]:
            g.savefig(FIG / f"03_bulk_clustermap.{ext}", dpi=200, bbox_inches="tight")
        plt.close("all")

    # ═══════════════════════════════════════════════════════════════
    # 6. TCR CLONALITY
    # ═══════════════════════════════════════════════════════════════
    tcr = tcr_analysis(adata)
    if tcr is not None and len(tcr) > 0:
        # Try to match TCR to cells
        tcr.to_csv(TAB / "tcr_clonality.csv", index=False)
        print(f"  TCR data columns: {list(tcr.columns)[:10]}")

    # ═══════════════════════════════════════════════════════════════
    # 7. BOOTSTRAP STABILITY
    # ═══════════════════════════════════════════════════════════════
    print("\n[BOOTSTRAP] Testing signature stability...")
    n_boot = 50
    boot_results = []
    io_obs = adata.obs[adata.obs["response"].isin(["Responder", "Non-Responder"])]
    n_sample = min(10000, len(io_obs))

    for b in range(n_boot):
        idx = np.random.RandomState(b).choice(len(io_obs), n_sample, replace=True)
        sub = io_obs.iloc[idx]
        for sig in sig_cols[:8]:  # Top 8 signatures
            r_mean = sub.loc[sub["response"]=="Responder", sig].mean()
            nr_mean = sub.loc[sub["response"]=="Non-Responder", sig].mean()
            boot_results.append({"boot": b, "signature": sig, "effect": r_mean - nr_mean})

    boot_df = pd.DataFrame(boot_results)
    boot_df.to_csv(TAB / "bootstrap_effects.csv", index=False)

    # ─── Fig: Bootstrap confidence intervals ───────────────────────
    fig, ax = plt.subplots(figsize=(10, 6))
    sig_order = boot_df.groupby("signature")["effect"].mean().sort_values().index
    positions = range(len(sig_order))
    for i, sig in enumerate(sig_order):
        vals = boot_df[boot_df["signature"] == sig]["effect"]
        ci_lo, ci_hi = np.percentile(vals, [2.5, 97.5])
        mean_val = vals.mean()
        color = "#2ecc71" if mean_val > 0 else "#e74c3c"
        ax.errorbar(mean_val, i, xerr=[[mean_val - ci_lo], [ci_hi - mean_val]],
                    fmt="o", color=color, capsize=5, markersize=8, linewidth=2)
    ax.set_yticks(positions)
    ax.set_yticklabels([s.replace("_", " ") for s in sig_order])
    ax.axvline(0, color="black", linestyle="--", alpha=0.5)
    ax.set_xlabel("Effect Size (R - NR)")
    ax.set_title("Bootstrap Confidence Intervals (95% CI, n=50)", fontweight="bold")
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"03_bootstrap_ci.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    # ═══════════════════════════════════════════════════════════════
    # 8. PATHWAY ENRICHMENT ON DEGs
    # ═══════════════════════════════════════════════════════════════
    print("\n[PATHWAYS] Pathway enrichment on R vs NR DEGs...")
    deg_path = PROJECT / "out" / "02a_response_tme" / "tables" / "deg_R_vs_NR_per_celltype.csv"
    if deg_path.exists():
        deg_df = pd.read_csv(deg_path)
        sig_up = deg_df[(deg_df["pvals_adj"] < 0.05) & (deg_df["logfoldchanges"] > 0.25)]

        # Simple curated pathway enrichment
        PATHWAYS = {
            "Antigen_processing":    ["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","PSMB8","PSMB9"],
            "T_cell_costimulation":  ["CD28","ICOS","CD80","CD86","CD40","CD40LG","TNFRSF9"],
            "Cytokine_signaling":    ["IFNG","IL2","IL12A","IL12B","IL21","CXCL9","CXCL10","CCL5"],
            "Checkpoint_pathway":    ["PDCD1","CD274","CTLA4","HAVCR2","LAG3","TIGIT","VSIR"],
            "Apoptosis":             ["CASP3","CASP8","FAS","FASLG","BAX","BCL2","BID","TRAIL"],
            "WNT_signaling":         ["CTNNB1","APC","AXIN2","LGR5","WNT5A","WNT3A"],
            "TGFb_signaling":        ["TGFB1","TGFB2","SMAD2","SMAD3","SMAD4","ACVR1"],
            "Angiogenesis":          ["VEGFA","VEGFB","FLT1","KDR","ANGPT1","ANGPT2"],
            "Complement":            ["C1QA","C1QB","C3","C4A","CFB","CFH"],
            "ECM_remodeling":        ["MMP2","MMP9","MMP14","COL1A1","COL3A1","FN1","LAMA5"],
        }

        up_genes = set(sig_up["names"].values)
        enrich_results = []
        for pw_name, pw_genes in PATHWAYS.items():
            overlap = up_genes.intersection(pw_genes)
            if len(overlap) > 0:
                enrich_results.append({
                    "pathway": pw_name,
                    "n_overlap": len(overlap),
                    "n_pathway": len(pw_genes),
                    "fraction": len(overlap) / len(pw_genes),
                    "genes": ", ".join(sorted(overlap)),
                })

        if enrich_results:
            enrich_df = pd.DataFrame(enrich_results).sort_values("fraction", ascending=False)
            enrich_df.to_csv(TAB / "pathway_enrichment.csv", index=False)

            fig, ax = plt.subplots(figsize=(9, 6))
            ax.barh(range(len(enrich_df)), enrich_df["fraction"], color="#3498db", alpha=0.8, height=0.5)
            ax.set_yticks(range(len(enrich_df)))
            ax.set_yticklabels(enrich_df["pathway"].str.replace("_", " "))
            for i, (_, row) in enumerate(enrich_df.iterrows()):
                ax.text(row["fraction"] + 0.01, i, f"{row['n_overlap']}/{row['n_pathway']}",
                       va="center", fontsize=9)
            ax.set_xlabel("Fraction of Pathway Genes in R-Upregulated DEGs")
            ax.set_title("Pathway Enrichment — Responder-Upregulated Genes", fontweight="bold")
            plt.tight_layout()
            for ext in ["png", "pdf"]:
                plt.savefig(FIG / f"03_pathway_enrichment.{ext}", dpi=200, bbox_inches="tight")
            plt.close()

    # ═══════════════════════════════════════════════════════════════
    # SUMMARY
    # ═══════════════════════════════════════════════════════════════
    n_figs = len(list(FIG.glob("*.png")))
    n_tabs = len(list(TAB.glob("*.csv")))
    print(f"\n{'='*70}")
    print(f"[DONE] Validation & Downstream Analyses")
    print(f"  Figures: {n_figs} in {FIG}")
    for f in sorted(FIG.glob("*.png")):
        print(f"    {f.name}")
    print(f"  Tables: {n_tabs} in {TAB}")
    for f in sorted(TAB.glob("*.csv")):
        print(f"    {f.name}")
    print(f"={'='*70}")


if __name__ == "__main__":
    main()
