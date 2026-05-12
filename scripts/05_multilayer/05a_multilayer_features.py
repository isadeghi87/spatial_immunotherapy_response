#!/usr/bin/env python3
"""
05a_multilayer_features.py — Build a patient-level multi-layer feature matrix
by combining five orthogonal axes of variation:

  L1 — Tumour-intrinsic genomic state (CNV burden, CNV clone diversity,
       HLA-loss proxy, chr8q amplification proxy)
  L2 — TF-regulon activity (master-regulator signatures: IRF1/STAT1 axis vs
       MYC/CTNNB1 axis; computed as score_ratio per patient)
  L3 — Cell-composition / TME abundance (% per cell-type)
  L4 — IO signatures (16 curated mechanistic signatures)
  L5 — Spatial niche features (TLS score, spatial entropy — from Visium)

Inputs
------
  data/processed/integrated_io_atlas_cnv.h5ad   (48k cells w/ CNV + signatures)
  out/02a_response_tme/tables/celltype_abundance_R_vs_NR.csv
  out/02b_spatial_biomarkers/tables/gse238264_per_sample_stats.csv
  out/02b_spatial_biomarkers/tables/gse238264_spatial_entropy.csv

Outputs
-------
  out/05_multilayer/tables/patient_multilayer_features.csv
  out/05_multilayer/tables/cnv_response_stats.csv
  out/05_multilayer/figures/05a_cnv_RvsNR.png            (4-panel)
  out/05_multilayer/figures/05a_multilayer_correlations.png
"""
import warnings
warnings.filterwarnings("ignore")

import os
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
sys.path.insert(0, str(PROJECT / "scripts" / "05_multilayer"))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import stats

from config.publication_style import apply_style, PAL, CMAP
from geneformer_utils import DEFAULT_MODEL_VARIANT, sanitize_model_name
apply_style()

OUT = PROJECT / "out" / "05_multilayer"
FIG = OUT / "figures"; TAB = OUT / "tables"
for d in (FIG, TAB): d.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("[05a] Multi-layer feature engineering")
print("=" * 70)

# ── Load atlases ────────────────────────────────────────────────────────
print("\n[L1] Loading CNV-annotated atlas …")
adata_cnv = sc.read_h5ad(PROJECT / "data" / "processed" / "integrated_io_atlas_cnv.h5ad")
print(f"   CNV atlas shape: {adata_cnv.shape}")
mask_cnv = adata_cnv.obs["response"].isin(["Responder", "Non-Responder"])
adata_cnv = adata_cnv[mask_cnv].copy()
print(f"   CNV labelled subset: {adata_cnv.shape}")

print("\n[L2-L4] Loading full integrated atlas …")
adata = sc.read_h5ad(PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad")
print(f"   Full atlas shape: {adata.shape}")
mask = adata.obs["response"].isin(["Responder", "Non-Responder"])
adata = adata[mask].copy()
print(f"   Full labelled subset: {adata.shape}")

# ── Patient-level feature table ─────────────────────────────────────────
# Identify a usable patient_id column
pid = "patient_ID" if "patient_ID" in adata.obs.columns else "sample_id"
adata.obs[pid] = adata.obs[pid].astype(str)
patients = adata.obs.groupby(pid).agg(
    response=("response", lambda x: x.value_counts().idxmax()),
    dataset=("dataset", "first"),
).reset_index()
patients.columns = ["patient_id", "response", "dataset"]
print(f"   Patients with labels: {len(patients)}  "
      f"(R={(patients.response=='Responder').sum()}, "
      f"NR={(patients.response=='Non-Responder').sum()})")
print(f"   By dataset: {patients.groupby('dataset')['response'].value_counts().to_dict()}")

# ── L1: CNV features per patient ────────────────────────────────────────
print("\n[L1] Computing CNV features per patient …")
pid_cnv = pid if pid in adata_cnv.obs.columns else ("patient_ID" if "patient_ID" in adata_cnv.obs.columns else "sample_id")
adata_cnv.obs[pid_cnv] = adata_cnv.obs[pid_cnv].astype(str)
is_malig = adata_cnv.obs["malignant_status"].astype(str).str.lower().str.contains("malig")
adata_cnv.obs["is_malignant"] = is_malig.astype(int)

cnv_per_pat = adata_cnv.obs.groupby(pid_cnv).agg(
    cnv_score_mean=("cnv_score", "mean"),
    cnv_score_sd=("cnv_score", "std"),
    malig_frac=("is_malignant", "mean"),
    n_cells=("is_malignant", "size"),
)
# CNV clone diversity (Shannon over cnv_leiden among malignant cells)
def _shannon(x):
    p = x.value_counts(normalize=True)
    return float(-(p * np.log(p + 1e-12)).sum())
mal_cells = adata_cnv.obs[adata_cnv.obs["is_malignant"] == 1]
if len(mal_cells):
    clone_div = mal_cells.groupby(pid_cnv)["cnv_leiden"].apply(_shannon)
    cnv_per_pat["cnv_clone_diversity"] = clone_div
else:
    cnv_per_pat["cnv_clone_diversity"] = np.nan
cnv_per_pat = cnv_per_pat.reset_index().rename(columns={pid_cnv: "patient_id"})

# ── L4: IO signatures per patient (mean over all cells) ─────────────────
print("[L4] Aggregating IO signatures per patient …")
SIGS = ["T_cell_exhaustion", "T_cell_activation", "T_cell_cytotoxicity",
        "Progenitor_CD8", "Terminal_exhaustion", "Treg_signature",
        "M1_macrophage", "M2_macrophage", "IFN_gamma", "TGF_beta",
        "Antigen_presentation", "Immune_checkpoint", "Angiogenesis",
        "WNT_beta_catenin", "mregDC", "Tfh_like"]
sig_cols = [c for c in SIGS if c in adata.obs.columns]
sig_per_pat = adata.obs.groupby(pid)[sig_cols].mean().reset_index()
sig_per_pat = sig_per_pat.rename(columns={pid: "patient_id"})
sig_per_pat.columns = ["patient_id"] + [f"sig_{c}" for c in sig_cols]

# ── L2: TF-regulon proxies (gene-score differentials) ──────────────────
# Since we don't have a SCENIC AUC matrix, build two master-regulator
# proxy axes using canonical IFN-γ (IRF1/STAT1) and Wnt/MYC axes.
print("[L2] Building TF-regulon proxy scores …")
irf1_axis_genes = ["IRF1", "STAT1", "CIITA", "NLRC5", "GBP1", "PSMB8",
                   "PSMB9", "TAP1", "TAP2"]
myc_axis_genes  = ["MYC", "CTNNB1", "AXIN2", "LGR5", "CCND1", "NOTCH1"]

def score_axis(genes, name):
    g_in = [g for g in genes if g in adata.var_names]
    if len(g_in) < 2:
        adata.obs[name] = np.nan
        return len(g_in)
    sc.tl.score_genes(adata, gene_list=g_in, score_name=name, random_state=0)
    return len(g_in)

n_irf1 = score_axis(irf1_axis_genes, "regulon_IRF1_axis")
n_myc  = score_axis(myc_axis_genes,  "regulon_MYC_axis")
print(f"   IRF1-axis genes: {n_irf1}/{len(irf1_axis_genes)}; "
      f"MYC-axis genes: {n_myc}/{len(myc_axis_genes)}")

reg_per_pat = adata.obs.groupby(pid)[
    ["regulon_IRF1_axis", "regulon_MYC_axis"]].mean().reset_index()
reg_per_pat = reg_per_pat.rename(columns={pid: "patient_id"})
reg_per_pat["regulon_IRF1_over_MYC"] = (reg_per_pat["regulon_IRF1_axis"]
                                         - reg_per_pat["regulon_MYC_axis"])

# ── L3: Cell-composition per patient ────────────────────────────────────
# IMPORTANT: Harmonise to a coarse label taxonomy shared across all datasets.
# GSE206325 has fine-grained labels (>40 subtypes); Cappuyns2023 has only 15
# coarse labels. Without harmonisation, fine-grained features such as
# frac_B_Plasma and frac_CD4_Tfh-like are structurally zero for Cappuyns,
# making cross-cohort L3 validation impossible.
print("[L3] Computing per-patient cell-type fractions (coarse-harmonised) …")
ct_col = "cell_type" if "cell_type" in adata.obs.columns else "celltype"
adata.obs[ct_col] = adata.obs[ct_col].astype(str)

# ------ coarse label harmonisation map -----------------------------------
_COARSE_MAP = {
    # B cells
    "B_Memory": "B_cell", "B_Naive": "B_cell",
    "B_Plasma": "Plasma",
    # CD4 T cells
    "CD4_Effector-Memory": "CD4_T", "CD4_Memory": "CD4_T",
    "CD4_Tfh-like": "CD4_T", "CD4_Th1 GZMK": "CD4_T", "CD4_Th1/Th17": "CD4_T",
    # CD8 T cells
    "CD8_Cytotoxic": "CD8_T",
    "CD8_Dysfunctional-Proliferating": "CD8_T",
    "CD8_Dysfunctional-effector": "CD8_T",
    "CD8_Dysfunctional-progenitor": "CD8_T",
    "CD8_Dysfunctional-terminal": "CD8_T",
    "CD8_Effector": "CD8_T",
    "CD8_Memory": "CD8_T",
    # Other T
    "Treg_Treg": "Treg", "MAIT": "CD8_T", "gdT": "CD8_T",
    "iNKT": "NK", "ILC": "NK",
    # DCs
    "DC_cDC1": "cDC1", "DC_cDC2": "cDC2",
    "DC_mregDC": "mregDC", "DC_pDC": "mregDC",
    # Myeloid
    "Macrophage_FCN1/CXCL9": "Macrophage",
    "Macrophage_FOLR2": "Macrophage",
    "Macrophage_I": "Macrophage", "Macrophage_II": "Macrophage",
    "Macrophage_Kupffer": "Macrophage",
    "Macrophage_SPP1": "Macrophage", "Macrophage_TREM2": "Macrophage",
    "Monocyte_CD14": "Monocyte", "Monocyte_CD14 SPP1": "Monocyte",
    "Monocyte_CD16": "Monocyte", "Monocyte_Intermediate": "Monocyte",
    "Mast": "Mast",
    # Stromal
    "Stromal_Endothelial cells": "Endothelial",
    "Stromal_HSC": "Fibroblast",
    "Stromal_Hepatocytes": "Hepatocyte",
    # Ambiguous / technical
    "Unknown_Doublets": "Unknown", "Unknown_Empty": "Unknown",
    "Naive_Naive": "Unknown", "Proliferating_Proliferating": "Unknown",
}
adata.obs["cell_type_coarse"] = (
    adata.obs[ct_col].map(_COARSE_MAP).fillna(adata.obs[ct_col])
)
n_unique_before = adata.obs[ct_col].nunique()
n_unique_after  = adata.obs["cell_type_coarse"].nunique()
print(f"   Labels before harmonisation: {n_unique_before}  →  after: {n_unique_after}")

ct_pat = (adata.obs.groupby([pid, "cell_type_coarse"]).size()
          .unstack(fill_value=0))
ct_pat = ct_pat.div(ct_pat.sum(axis=1), axis=0)  # fractions
ct_pat.columns = [f"frac_{c}" for c in ct_pat.columns]
ct_pat = ct_pat.reset_index().rename(columns={pid: "patient_id"})

# Verify Cappuyns Plasma fraction is non-zero
cap_idx = ct_pat["patient_id"].isin(
    adata.obs.loc[adata.obs["dataset"] == "Cappuyns2023", pid]
          .astype(str).unique()
) if "dataset" in adata.obs.columns else pd.Series([False] * len(ct_pat))
if cap_idx.any() and "frac_Plasma" in ct_pat.columns:
    cap_plasma = ct_pat.loc[cap_idx, "frac_Plasma"]
    print(f"   Cappuyns frac_Plasma — mean={cap_plasma.mean():.4f}, "
          f"non-zero patients: {(cap_plasma > 0).sum()}/{cap_idx.sum()}")

# ── L5: Spatial niche features — from Visium summaries ─────────────────
print("[L5] Loading spatial niche features (Visium) …")
try:
    sp_stats = pd.read_csv(
        PROJECT / "out" / "02b_spatial_biomarkers" / "tables"
        / "gse238264_per_sample_stats.csv")
    # per-sample table; map back via response and use as cohort-level feature
    sp_r = sp_stats.groupby("response_group", observed=True).mean(numeric_only=True) \
        if "response_group" in sp_stats.columns else None
    print(f"   Visium per-sample: {sp_stats.shape[0]} rows, "
          f"{sp_stats.shape[1]} cols")
except Exception as e:
    print(f"   [WARN] Visium per-sample stats not usable: {e}")
    sp_stats = None

# ── Merge all layers ────────────────────────────────────────────────────
print("\n[MERGE] Combining all layers …")
feats = (patients
         .merge(cnv_per_pat, on="patient_id", how="left")
         .merge(reg_per_pat, on="patient_id", how="left")
         .merge(ct_pat, on="patient_id", how="left")
         .merge(sig_per_pat, on="patient_id", how="left"))

geneformer_env = os.environ.get("GENEFORMER_FEATURE_FILE")
if geneformer_env:
    geneformer_path = Path(geneformer_env)
    if not geneformer_path.is_absolute():
        geneformer_path = PROJECT / geneformer_path
else:
    geneformer_model_tag = sanitize_model_name(DEFAULT_MODEL_VARIANT)
    geneformer_path = TAB / f"geneformer_patient_features_{geneformer_model_tag}.csv"

if geneformer_path.exists():
    print(f"   Merging Geneformer L6: {geneformer_path}")
    geneformer = pd.read_csv(geneformer_path)
    geneformer_cols = [
        column for column in geneformer.columns
        if column == "patient_id" or column.startswith("gf_cls_") or column in {"n_cells_geneformer", "geneformer_model"}
    ]
    feats = feats.merge(geneformer[geneformer_cols], on="patient_id", how="left")
    print(f"   Geneformer dims : {len([c for c in feats.columns if c.startswith('gf_cls_')])}")
else:
    print(f"   Geneformer L6   : not merged ({geneformer_path.name} not found)")

print(f"   Feature matrix: {feats.shape[0]} patients × {feats.shape[1]} cols")

feats.to_csv(TAB / "patient_multilayer_features.csv", index=False)
print(f"   Saved: {TAB / 'patient_multilayer_features.csv'}")

disc_mask = feats["dataset"].astype(str) != "Cappuyns2023"
feats_disc = feats.loc[disc_mask].copy()
print(f"   Discovery modelling cohort: {len(feats_disc)} patients "
      f"(R={(feats_disc.response=='Responder').sum()}, "
      f"NR={(feats_disc.response=='Non-Responder').sum()})")

# ── Statistical test: CNV features R vs NR ─────────────────────────────
print("\n[STAT] CNV features R vs NR …")
from scipy.stats import mannwhitneyu
cnv_feats = ["cnv_score_mean", "cnv_score_sd", "malig_frac", "cnv_clone_diversity"]
stats_rows = []
for f in cnv_feats:
    r  = feats_disc.loc[feats_disc.response == "Responder",     f].dropna()
    nr = feats_disc.loc[feats_disc.response == "Non-Responder", f].dropna()
    if len(r) < 3 or len(nr) < 3:
        continue
    U, p = mannwhitneyu(r, nr, alternative="two-sided")
    eff = (r.mean() - nr.mean()) / (np.concatenate([r, nr]).std() + 1e-12)
    stats_rows.append(dict(feature=f,
                            n_R=len(r), n_NR=len(nr),
                            mean_R=float(r.mean()), mean_NR=float(nr.mean()),
                            effect_size=float(eff), U=float(U), p_value=float(p)))
cnv_stats = pd.DataFrame(stats_rows)
cnv_stats.to_csv(TAB / "cnv_response_stats.csv", index=False)
print(cnv_stats.to_string(index=False))

# ── Figure: CNV features R vs NR ────────────────────────────────────────
print("\n[FIG] CNV R-vs-NR comparison …")
fig, axes = plt.subplots(1, 4, figsize=(18, 4.5), constrained_layout=True)
for ax, f in zip(axes, cnv_feats):
    if f not in feats_disc.columns:
        ax.set_visible(False); continue
    sub = feats_disc[["response", f]].dropna()
    order = ["Non-Responder", "Responder"]
    sns.boxplot(data=sub, x="response", y=f, order=order, ax=ax,
                palette=[PAL["Non_Responder"], PAL["Responder"]],
                width=0.55, fliersize=0, linewidth=1)
    sns.stripplot(data=sub, x="response", y=f, order=order, ax=ax,
                   color="black", size=3.5, alpha=0.75, jitter=0.18)
    row = cnv_stats[cnv_stats.feature == f]
    p = row["p_value"].values[0] if len(row) else np.nan
    pretty = {"cnv_score_mean": "CNV score (mean)",
              "cnv_score_sd": "CNV score (SD)",
              "malig_frac": "Malignant fraction",
              "cnv_clone_diversity": "CNV clone diversity\n(Shannon)"}.get(f, f)
    ax.set_title(f"{pretty}\nMW p = {p:.2g}" if np.isfinite(p) else pretty,
                 fontsize=10)
    ax.set_xlabel("")
    ax.set_ylabel(pretty.split("\n")[0])
    ax.tick_params(axis="x", rotation=20)
fig.suptitle("Tumour-intrinsic genomic state — discovery cohort (GSE206325)",
             fontsize=12, fontweight="bold")
for fmt in ("png", "pdf"):
    fig.savefig(FIG / f"05a_cnv_RvsNR.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"   Saved: {FIG / '05a_cnv_RvsNR.png'}")

# ── Figure: multi-layer correlation heatmap ────────────────────────────
print("[FIG] Cross-layer correlation heatmap …")
num_cols = feats.select_dtypes(include=[np.number]).columns.tolist()
# keep only columns with <50% missing
keep = [c for c in num_cols if feats[c].notna().mean() > 0.5]
corr = feats[keep].corr(method="spearman")
fig, ax = plt.subplots(figsize=(12, 10), constrained_layout=True)
sns.heatmap(corr, cmap=CMAP["correlation"], center=0, vmin=-1, vmax=1,
            square=False, cbar_kws=dict(label="Spearman ρ", shrink=0.6),
            xticklabels=True, yticklabels=True, ax=ax, linewidths=0)
ax.tick_params(axis="x", labelsize=6, rotation=90)
ax.tick_params(axis="y", labelsize=6)
ax.set_title("Multi-layer feature correlation structure\n"
             "(CNV · regulon · TME composition · IO signatures)",
             fontsize=11, fontweight="bold")
for fmt in ("png", "pdf"):
    fig.savefig(FIG / f"05a_multilayer_correlations.{fmt}",
                dpi=300, bbox_inches="tight")
plt.close(fig)

# ── Sample-level prevalence: what fraction of patients carry each result? ─
print("\n[PREVALENCE] Per-feature carrier prevalence + carrier-vs-other comparison …")
from scripts._shared.sample_prevalence import prevalence_panel

# Pick the strongest discovery features per layer for the prevalence audit.
PREVALENCE_FEATURES = [c for c in [
    "cnv_clone_diversity", "cnv_score_sd", "malig_frac",   # L1
    "regulon_score_ratio",                                  # L2
    "frac_CD8+ T cell", "frac_Treg", "frac_Macrophage",     # L3 (if present)
    "IFN_gamma", "T_cell_exhaustion", "TLS_signature",      # L4 (if present)
    "tls_score", "spatial_entropy",                          # L5 (if present)
] if c in feats_disc.columns]

PREVALENCE_PROPERTIES = [c for c in [
    "malig_frac", "cnv_score_mean", "cnv_clone_diversity",
    "frac_CD8+ T cell", "frac_Treg", "frac_Macrophage",
    "IFN_gamma", "TLS_signature", "tls_score", "spatial_entropy",
    "regulon_score_ratio",
] if c in feats_disc.columns]

prev_summary = prevalence_panel(
    feats_disc,
    features=PREVALENCE_FEATURES,
    compare_cols=PREVALENCE_PROPERTIES,
    group_col="response",
    threshold="q75",
    direction="high",
    out_dir=TAB,
    fig_dir=FIG,
    tag="05a_multilayer",
    title="Per-feature carrier prevalence (Q3 cutoff) — discovery cohort (GSE206325)",
)
if not prev_summary.empty:
    print(f"   Prevalence rows: {len(prev_summary)} → {TAB / '05a_multilayer_prevalence.csv'}")

print("\n[DONE] Multi-layer features ready.")
