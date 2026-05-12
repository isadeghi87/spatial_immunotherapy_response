#!/usr/bin/env python3
"""
03d_spatial_spot_qc.py  — T3.1

Per-spot QC across Visium cohorts (GSE238264 + GSE281759).
Computes:
- n_genes_by_counts, total_counts, pct_counts_mt per spot
- spot pass-rate under canonical thresholds
  (n_genes >= 200, total_counts >= 500, pct_mt <= 20)
- per-section summary (medians, IQR, pass-rate, n_spots)
- per-cohort aggregate

Outputs
-------
out/03_validation/tables/spatial_spot_qc.csv            (per-section summary)
out/03_validation/tables/spatial_spot_qc_cohort.csv     (per-cohort aggregate)
out/03_validation/figures/03d_spatial_spot_qc.{png,pdf}
"""
from __future__ import annotations
import sys, warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from config.publication_style import apply_style, PAL
apply_style()

TAB = PROJECT / "out" / "03_validation" / "tables"; TAB.mkdir(parents=True, exist_ok=True)
FIG = PROJECT / "out" / "03_validation" / "figures"; FIG.mkdir(parents=True, exist_ok=True)

COHORTS = {
    "GSE238264": PROJECT / "data" / "spatial_data" / "GSE238264_HCC_IO_Visium",
    "GSE281759": PROJECT / "data" / "spatial_data" / "GSE281759_HCC_MultiRegion",
}
# Thresholds (standard 10x Visium QC)
MIN_N_GENES, MIN_COUNTS, MAX_PCT_MT = 200, 500, 20.0

def qc_one(h5: Path, cohort: str):
    ad = sc.read_h5ad(h5)
    # ensure mt annotation — use np.asarray for safety across pandas versions
    var_names = np.asarray(ad.var_names).astype(str)
    mt_mask = np.array([v.upper().startswith(("MT-", "MT."))
                        for v in var_names])
    ad.var["mt"] = mt_mask
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True,
                               percent_top=None, log1p=False)
    obs = ad.obs.copy()
    ng = obs["n_genes_by_counts"].astype(float).values
    tc = obs["total_counts"].astype(float).values
    pm = obs["pct_counts_mt"].astype(float).values if "pct_counts_mt" in obs else np.zeros_like(ng)
    passed = (ng >= MIN_N_GENES) & (tc >= MIN_COUNTS) & (pm <= MAX_PCT_MT)
    section = h5.stem
    # responder label encoded in name for GSE238264 (HCCxR / HCCxNR)
    resp = "Responder" if section.endswith("R") and not section.endswith("NR") else (
           "Non-Responder" if section.endswith("NR") else "n/a")
    return dict(
        cohort=cohort, section=section, response=resp,
        n_spots=int(len(obs)),
        n_genes_median=float(np.median(ng)),
        n_genes_iqr_lo=float(np.percentile(ng, 25)),
        n_genes_iqr_hi=float(np.percentile(ng, 75)),
        total_counts_median=float(np.median(tc)),
        total_counts_iqr_lo=float(np.percentile(tc, 25)),
        total_counts_iqr_hi=float(np.percentile(tc, 75)),
        pct_mt_median=float(np.median(pm)),
        pct_mt_iqr_lo=float(np.percentile(pm, 25)),
        pct_mt_iqr_hi=float(np.percentile(pm, 75)),
        spots_pass=int(passed.sum()),
        spots_pass_frac=float(passed.mean()),
        min_n_genes_thr=MIN_N_GENES,
        min_counts_thr=MIN_COUNTS,
        max_pct_mt_thr=MAX_PCT_MT,
    ), ng, tc, pm

rows = []
raw_by_section = {}
for cohort, d in COHORTS.items():
    if not d.exists(): continue
    files = sorted(d.glob("*.h5ad"))
    print(f"[{cohort}] {len(files)} sections in {d}")
    for f in files:
        try:
            r, ng, tc, pm = qc_one(f, cohort)
            rows.append(r)
            raw_by_section[r["section"]] = (ng, tc, pm)
            print(f"   {r['section']:30s}  "
                  f"spots={r['n_spots']:5d}  "
                  f"n_genes={r['n_genes_median']:.0f}  "
                  f"counts={r['total_counts_median']:.0f}  "
                  f"mt%={r['pct_mt_median']:.1f}  "
                  f"pass={r['spots_pass_frac']*100:.1f}%")
        except Exception as e:
            print(f"   [WARN] {f.name}: {e}")

df = pd.DataFrame(rows)
df.to_csv(TAB / "spatial_spot_qc.csv", index=False)
print(f"\nwrote {TAB/'spatial_spot_qc.csv'} ({len(df)} rows)")

# Per-cohort aggregate (weighted by n_spots)
agg_rows = []
for coh, sub in df.groupby("cohort"):
    w = sub["n_spots"].values
    agg_rows.append(dict(
        cohort=coh,
        n_sections=len(sub),
        total_spots=int(sub["n_spots"].sum()),
        n_genes_median=float(np.average(sub["n_genes_median"], weights=w)),
        total_counts_median=float(np.average(sub["total_counts_median"], weights=w)),
        pct_mt_median=float(np.average(sub["pct_mt_median"], weights=w)),
        spots_pass_frac=float((sub["spots_pass"].sum()) / sub["n_spots"].sum()),
    ))
agg = pd.DataFrame(agg_rows)
agg.to_csv(TAB / "spatial_spot_qc_cohort.csv", index=False)
print(f"wrote {TAB/'spatial_spot_qc_cohort.csv'}")
print(agg.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

# ── Figure ──────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.8), constrained_layout=True)

# A — boxplot n_genes by section (GSE238264 only for readability)
axA = axes[0, 0]
g238 = df[df["cohort"] == "GSE238264"].sort_values("section")
sec_order = g238["section"].tolist()
data = [raw_by_section[s][0] for s in sec_order]
colors = [PAL.get("Responder", "#0077BB") if s.endswith("R") and not s.endswith("NR")
          else PAL.get("Non_Responder", "#CC3311") for s in sec_order]
bp = axA.boxplot(data, patch_artist=True, widths=0.7, showfliers=False,
                 medianprops=dict(color="black", linewidth=1.2))
for patch, c in zip(bp["boxes"], colors):
    patch.set_facecolor(c); patch.set_alpha(0.55); patch.set_edgecolor("black")
axA.set_xticks(np.arange(1, len(sec_order)+1))
axA.set_xticklabels(sec_order, rotation=35, ha="right", fontsize=8.5)
axA.axhline(MIN_N_GENES, ls="--", color="black", lw=0.7, label=f"threshold {MIN_N_GENES}")
axA.set_ylabel("n_genes_by_counts (per spot)", fontsize=9)
axA.set_title("A  GSE238264 per-section gene coverage",
              loc="left", fontsize=10.5, fontweight="bold")
axA.legend(frameon=False, fontsize=8)

# B — boxplot total_counts
axB = axes[0, 1]
data2 = [raw_by_section[s][1] for s in sec_order]
bp2 = axB.boxplot(data2, patch_artist=True, widths=0.7, showfliers=False,
                  medianprops=dict(color="black", linewidth=1.2))
for patch, c in zip(bp2["boxes"], colors):
    patch.set_facecolor(c); patch.set_alpha(0.55); patch.set_edgecolor("black")
axB.set_xticks(np.arange(1, len(sec_order)+1))
axB.set_xticklabels(sec_order, rotation=35, ha="right", fontsize=8.5)
axB.axhline(MIN_COUNTS, ls="--", color="black", lw=0.7, label=f"threshold {MIN_COUNTS}")
axB.set_ylabel("total_counts (per spot)", fontsize=9)
axB.set_yscale("log")
axB.set_title("B  GSE238264 per-section UMI depth (log)",
              loc="left", fontsize=10.5, fontweight="bold")
axB.legend(frameon=False, fontsize=8)

# C — pass-rate per section, both cohorts stacked
axC = axes[1, 0]
df_sorted = df.sort_values(["cohort", "section"]).reset_index(drop=True)
y = np.arange(len(df_sorted))
cols_c = [PAL.get("Responder", "#0077BB") if r == "Responder"
          else PAL.get("Non_Responder", "#CC3311") if r == "Non-Responder"
          else "#888888" for r in df_sorted["response"]]
axC.barh(y, df_sorted["spots_pass_frac"] * 100, color=cols_c,
         edgecolor="black", linewidth=0.4)
for i, r in df_sorted.iterrows():
    axC.text(r["spots_pass_frac"] * 100 + 0.5, i,
             f"{r['spots_pass']}/{r['n_spots']}",
             va="center", fontsize=7)
axC.set_yticks(y)
axC.set_yticklabels([f"{r['cohort']}·{r['section']}" for _, r in df_sorted.iterrows()],
                    fontsize=6.5)
axC.set_xlabel("% spots passing QC (n_genes≥200, counts≥500, mt≤20%)", fontsize=9)
axC.set_title("C  Spot-level pass rate (all sections, both cohorts)",
              loc="left", fontsize=10.5, fontweight="bold")
axC.set_xlim(0, 105)

# D — pct_mt distribution per section (violin)
axD = axes[1, 1]
g238_pm = [raw_by_section[s][2] for s in sec_order]
parts = axD.violinplot(g238_pm, positions=np.arange(len(sec_order)),
                       widths=0.75, showmedians=True, showextrema=False)
for body, c in zip(parts["bodies"], colors):
    body.set_facecolor(c); body.set_alpha(0.55); body.set_edgecolor("black")
axD.axhline(MAX_PCT_MT, ls="--", color="black", lw=0.7, label=f"threshold {MAX_PCT_MT}%")
axD.set_xticks(np.arange(len(sec_order)))
axD.set_xticklabels(sec_order, rotation=35, ha="right", fontsize=8.5)
axD.set_ylabel("pct_counts_mt (per spot)", fontsize=9)
axD.set_title("D  GSE238264 mitochondrial fraction",
              loc="left", fontsize=10.5, fontweight="bold")
axD.legend(frameon=False, fontsize=8)

fig.suptitle("Figure: Spatial spot-level QC across Visium cohorts",
             fontsize=12, fontweight="bold")
fig.savefig(FIG / "03d_spatial_spot_qc.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(FIG / "03d_spatial_spot_qc.pdf", bbox_inches="tight", facecolor="white")
print(f"\nwrote {FIG/'03d_spatial_spot_qc.png'}")
