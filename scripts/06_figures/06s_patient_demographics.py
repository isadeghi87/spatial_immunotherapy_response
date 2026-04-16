#!/usr/bin/env python3
"""Build patient / sample demographics for supplementary results.

Outputs
-------
out/06_manuscript_figures/SupplementaryFig_S9_patient_demographics.png / .pdf
out/06_manuscript_figures/tables/S12_patient_demographics.csv
out/06_manuscript_figures/tables/S13_sample_overview.csv
out/06_manuscript_figures/tables/S14_cohort_summary.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.publication_style import apply_style, PAL  # type: ignore

apply_style()

ROOT = Path(__file__).resolve().parents[2]
OUT_FIG = ROOT / "out/06_manuscript_figures"
OUT_TAB = OUT_FIG / "tables"
OUT_TAB.mkdir(parents=True, exist_ok=True)

# ── 1. Per-patient GSE206325 demographics (primary IO cohort) ─────────
sample_csv = ROOT / "data/processed/GSE206325_sample_annotations.csv"
samp = pd.read_csv(sample_csv)

# Clean response label
samp["response"] = samp["treatment_Resp"].map(
    lambda x: "Responder" if str(x).endswith("_R") else ("Non-Responder" if str(x).endswith("_NR") else "Unknown")
)

# Per-patient aggregation
pats = (
    samp.groupby("patient_ID")
    .agg(
        n_samples=("sample_ID", "nunique"),
        tissues=("tissue", lambda s: ",".join(sorted(set(map(str, s))))),
        treatment=("treatment", lambda s: sorted(set(s))[0] if len(s) else "NA"),
        disease=("disease", lambda s: sorted(set(s))[0] if len(s) else "NA"),
        response=("response", lambda s: sorted(set(s))[0] if len(s) else "Unknown"),
        immune_infiltration=("Immune_Infiltration", lambda s: sorted(set(map(str, s)))[0] if len(s) else "NA"),
    )
    .reset_index()
)
pats.to_csv(OUT_TAB / "S12_patient_demographics.csv", index=False)

# Sample-level overview
samp_tab = samp[[
    "sample_ID", "patient_ID", "tissue", "treatment", "disease",
    "library_chemistry", "response", "Immune_Infiltration",
]].copy()
samp_tab.to_csv(OUT_TAB / "S13_sample_overview.csv", index=False)

# ── 2. Cross-cohort summary (5 scRNA studies + 2 Visium + 2 bulk) ──────
cohort_rows = [
    ("GSE206325", "scRNA-seq (10x)",   "Magen 2023 Nat Med",         30, 682_000, "anti-PD-1 R/NR"),
    ("GSE125449", "scRNA-seq (10x)",   "Ma 2019 Cancer Cell",         8,  10_000, "treatment-naive"),
    ("GSE140228", "scRNA-seq (10x)",   "Zhang 2019 Cell",             8,  67_000, "CD45+ immune"),
    ("GSE149614", "scRNA-seq (10x)",   "Lu 2022 Nat Commun",         10,  70_000, "multi-region"),
    ("E-MTAB-10553","scRNA-seq (10x)", "Pfister 2021 Nature",        12, 150_000, "NASH + anti-PD-1"),
    ("GSE173905", "Visium",            "Wu 2021 Sci Adv",             8,       0, "NC / ER / RS"),
    ("GSE238264", "Visium",            "Liu 2023 Cell Rep Med",       7,       0, "ICB-treated (4R / 3NR)"),
    ("IMbrave150","Bulk RNA-seq",      "Finn 2020 NEJM",            358,       0, "Atezo+Bev survival"),
    ("CheckMate-459","Bulk RNA-seq",   "Yau 2022 Lancet Oncol",     743,       0, "Nivo vs sorafenib"),
]
cohort = pd.DataFrame(cohort_rows, columns=[
    "Accession", "Modality", "Source", "Patients", "Cells", "Notes"])
cohort.to_csv(OUT_TAB / "S14_cohort_summary.csv", index=False)

# ── 3. Supplementary Figure S9 panels ──────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(14, 8))

# (a) Response distribution per patient
resp_counts = pats["response"].value_counts().reindex(
    ["Responder", "Non-Responder", "Unknown"]).fillna(0).astype(int)
colors = [PAL["Responder"], PAL["Non-Responder"], PAL["Unknown"]]
axes[0, 0].bar(resp_counts.index, resp_counts.values, color=colors, width=0.55, edgecolor="black")
for i, v in enumerate(resp_counts.values):
    axes[0, 0].text(i, v + 0.3, str(v), ha="center", fontsize=10)
axes[0, 0].set_title("a. Patient response labels (GSE206325)")
axes[0, 0].set_ylabel("n patients")
axes[0, 0].set_ylim(0, max(resp_counts.max(), 1) + 3)

# (b) Samples per tissue
tissue_counts = samp["tissue"].value_counts()
axes[0, 1].bar(tissue_counts.index, tissue_counts.values, color="#4477AA", width=0.55, edgecolor="black")
for i, v in enumerate(tissue_counts.values):
    axes[0, 1].text(i, v + 0.3, str(v), ha="center", fontsize=10)
axes[0, 1].set_title("b. Sample distribution by tissue")
axes[0, 1].set_ylabel("n samples")

# (c) Samples per patient histogram
sp = pats["n_samples"].values
axes[0, 2].hist(sp, bins=np.arange(sp.min(), sp.max() + 2) - 0.5,
                color="#4477AA", edgecolor="black")
axes[0, 2].set_xlabel("samples per patient")
axes[0, 2].set_ylabel("n patients")
axes[0, 2].set_title("c. Samples collected per patient")
_xt = np.arange(int(sp.min()), int(sp.max()) + 1)
axes[0, 2].set_xticks(_xt[::2])

# (d) Cross-cohort cells
bulk_mask = cohort["Cells"] > 0
c_sc = cohort[bulk_mask]
axes[1, 0].barh(c_sc["Accession"], c_sc["Cells"] / 1000, color="#4477AA", edgecolor="black")
axes[1, 0].set_xlabel("cells × 10³")
axes[1, 0].set_title("d. Cells per scRNA-seq cohort")
axes[1, 0].invert_yaxis()

# (e) Cross-cohort patients
axes[1, 1].barh(cohort["Accession"], cohort["Patients"], color="#228833", edgecolor="black")
axes[1, 1].set_xlabel("n patients")
axes[1, 1].set_title("e. Patients per cohort (all modalities)")
axes[1, 1].invert_yaxis()

# (f) Immune infiltration categories in labelled subset
inf = pats.loc[pats["response"].isin(["Responder", "Non-Responder"]), ["response", "immune_infiltration"]].copy()
inf["immune_infiltration"] = inf["immune_infiltration"].replace({"nan": "Unannotated"}).fillna("Unannotated")
ct = pd.crosstab(inf["immune_infiltration"], inf["response"]).reindex(
    columns=["Responder", "Non-Responder"]).fillna(0)
ct.plot(kind="bar", stacked=False,
        color=[PAL["Responder"], PAL["Non-Responder"]],
        ax=axes[1, 2], edgecolor="black", width=0.7)
axes[1, 2].set_title("f. Immune infiltration × response")
axes[1, 2].set_ylabel("n patients")
axes[1, 2].set_xlabel("")
axes[1, 2].tick_params(axis="x", rotation=20)
axes[1, 2].legend(title="", frameon=False, fontsize=9)

for ax in axes.ravel():
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=9)

fig.suptitle("Supplementary Figure S9 | Patient and cohort demographics",
             fontsize=13, y=1.01)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(OUT_FIG / f"SupplementaryFig_S9_patient_demographics.{ext}",
                dpi=300, bbox_inches="tight")
plt.close(fig)

print("Wrote:",
      OUT_TAB / "S12_patient_demographics.csv",
      OUT_TAB / "S13_sample_overview.csv",
      OUT_TAB / "S14_cohort_summary.csv",
      OUT_FIG / "SupplementaryFig_S9_patient_demographics.png",
      sep="\n  ")
