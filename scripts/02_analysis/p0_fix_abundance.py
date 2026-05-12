#!/usr/bin/env python3
"""P0.3 — Recompute celltype_abundance_R_vs_NR.csv from atlas obs.

The previous CSV had empty mean/p columns. Reads the integrated atlas's obs
table only (no expression matrix needed), aggregates per-patient cell-type
fractions, runs Mann–Whitney R vs NR with BH correction, writes table.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import anndata
from scipy.stats import mannwhitneyu

PROJECT = Path(__file__).resolve().parents[2]
ATLAS = PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad"
OUT = PROJECT / "out" / "02a_response_tme" / "tables" / "celltype_abundance_R_vs_NR.csv"

a = anndata.read_h5ad(ATLAS, backed="r")
obs = a.obs[["patient_ID", "response", "cell_type"]].copy()
print(f"[ATLAS] cells={len(obs):,}")
print(obs["response"].value_counts(dropna=False))
print(f"unique patients (all): {obs['patient_ID'].nunique()}")

io = obs[obs["response"].isin(["Responder", "Non-Responder"])].copy()
print(f"IO-labelled cells: {len(io):,}; patients: {io['patient_ID'].nunique()}")
print(io.groupby("response")["patient_ID"].nunique())

# Patient × cell_type counts → per-patient fraction
pcm = (io.groupby(["patient_ID", "response", "cell_type"], observed=True)
          .size().rename("count").reset_index())
pcm["fraction"] = pcm.groupby("patient_ID")["count"].transform(lambda s: s / s.sum())

# Test per cell type
results = []
for ct, sub in pcm.groupby("cell_type", observed=True):
    r  = sub.loc[sub["response"] == "Responder",     "fraction"].to_numpy()
    nr = sub.loc[sub["response"] == "Non-Responder", "fraction"].to_numpy()
    if len(r) < 3 or len(nr) < 3:
        continue
    stat, pval = mannwhitneyu(r, nr, alternative="two-sided")
    results.append({
        "cell_type":     ct,
        "mean_R":        float(np.mean(r)),
        "mean_NR":       float(np.mean(nr)),
        "median_R":      float(np.median(r)),
        "median_NR":     float(np.median(nr)),
        "log2FC":        float(np.log2((np.mean(r)+1e-6)/(np.mean(nr)+1e-6))),
        "U_stat":        float(stat),
        "pvalue":        float(pval),
        "n_R_patients":  int(len(r)),
        "n_NR_patients": int(len(nr)),
    })

df = pd.DataFrame(results).sort_values("pvalue").reset_index(drop=True)
# BH (Benjamini–Hochberg)
m = len(df)
df["padj"] = (df["pvalue"] * m / (df.index + 1)).clip(upper=1.0)
df["padj"] = df["padj"][::-1].cummin()[::-1]

OUT.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT, index=False)
print(f"[OUT] {OUT} ({len(df)} rows)")
print(df.to_string(index=False))
