#!/usr/bin/env python3
"""02g_infercnv_spatial_arm_burden.py

Per-chromosome and per-arm CNV burden in tumor spots:
  - Aggregate window-level CNV scores within each chromosome arm (p, q split
    on the median window position within chr; approximate but consistent).
  - Per-patient mean amplitude (signed) and burden (|CNV|) on tumor spots.
  - Mann–Whitney + BH (Responder vs Non_Responder).
  - Secondary: GSE281759 primary-vs-recurrence (sample-level T1P vs T1R).

Outputs:
  out/02b_spatial_biomarkers/cnv/tables/spatial_arm_burden_per_patient.csv
  out/02b_spatial_biomarkers/cnv/tables/spatial_arm_burden_RvsNR.csv
  out/02b_spatial_biomarkers/cnv/tables/spatial_arm_burden_PvsR_GSE281759.csv
  out/02b_spatial_biomarkers/cnv/figures/spatial_arm_burden_RvsNR.{png,pdf}
"""
from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

PROJECT  = Path(__file__).resolve().parents[2]
NPZ_DIR  = PROJECT / "out" / "02b_spatial_biomarkers" / "cnv" / "npz"
TAB_DIR  = PROJECT / "out" / "02b_spatial_biomarkers" / "cnv" / "tables"
FIG_DIR  = PROJECT / "out" / "02b_spatial_biomarkers" / "cnv" / "figures"
TAB_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

TUMOR_PCT = 75.0  # within-sample tumor cutoff on cnv_score (obs spots)

CHRS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# ─────────────────────────────────────────────────────────────────────────────
# Build per-(patient, arm) burden table
# ─────────────────────────────────────────────────────────────────────────────
files = sorted(NPZ_DIR.glob("*_cnv.npz"))
print(f"Loading {len(files)} npz files …")

records = []
for f in files:
    d = np.load(f, allow_pickle=True)
    X = d["X_cnv"]
    n_w = X.shape[1]
    chr_keys = list(d["chr_pos_keys"])
    chr_vals = list(d["chr_pos_values"]) + [n_w]
    obs_mask = (d["groups"] == "observation")
    if obs_mask.sum() == 0:
        continue
    cnv_score = d["cnv_score"][obs_mask].astype(float)
    if len(cnv_score) == 0:
        continue
    thr = np.percentile(cnv_score, TUMOR_PCT)
    full_obs = obs_mask.copy()
    tumor_mask = full_obs & (d["cnv_score"].astype(float) >= thr)
    if tumor_mask.sum() == 0:
        continue
    Xt = X[tumor_mask]
    sample_id = str(d["sample_id"])
    dataset   = str(d["dataset"])
    response  = str(d["response"])
    for i, chrom in enumerate(chr_keys):
        a, b = chr_vals[i], chr_vals[i + 1]
        if b - a < 4:
            continue
        block = Xt[:, a:b]
        if block.shape[1] == 0:
            continue
        mid = (a + b) // 2
        # arms (approx by window count)
        for arm, sl in [
            (chrom + "p", slice(a, mid)),
            (chrom + "q", slice(mid, b)),
            (chrom,       slice(a, b)),  # whole chr
        ]:
            sub = X[:, sl][tumor_mask]
            if sub.shape[1] == 0:
                continue
            mean_amp  = float(sub.mean())
            mean_burd = float(np.abs(sub).mean())
            records.append(dict(
                sample_id=sample_id, dataset=dataset, response=response,
                arm=arm, mean_amp=mean_amp, mean_burden=mean_burd,
                n_tumor_spots=int(tumor_mask.sum()),
            ))

per_pat = pd.DataFrame.from_records(records)
print(f"Per-patient × arm rows: {len(per_pat):,}")
out1 = TAB_DIR / "spatial_arm_burden_per_patient.csv"
per_pat.to_csv(out1, index=False)
print(f"  → {out1.name}")

# ─────────────────────────────────────────────────────────────────────────────
# R vs NR test (use GSE238264 — annotated response)
# ─────────────────────────────────────────────────────────────────────────────
def bh(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out

rvr = per_pat[per_pat["response"].isin(["Responder", "Non_Responder"])].copy()
rows = []
for arm, sub in rvr.groupby("arm"):
    r  = sub[sub["response"] == "Responder"]
    nr = sub[sub["response"] == "Non_Responder"]
    if len(r) < 2 or len(nr) < 2:
        continue
    for metric in ("mean_amp", "mean_burden"):
        try:
            U, p = mannwhitneyu(r[metric].values, nr[metric].values, alternative="two-sided")
        except ValueError:
            continue
        rows.append(dict(
            arm=arm, metric=metric,
            n_R=len(r), n_NR=len(nr),
            median_R=float(r[metric].median()),
            median_NR=float(nr[metric].median()),
            delta=float(r[metric].median() - nr[metric].median()),
            U=float(U), p=float(p),
        ))
res = pd.DataFrame(rows)
if len(res):
    for metric in res["metric"].unique():
        m = res["metric"] == metric
        res.loc[m, "padj"] = bh(res.loc[m, "p"].values)
    res = res.sort_values(["metric", "padj"]).reset_index(drop=True)
out2 = TAB_DIR / "spatial_arm_burden_RvsNR.csv"
res.to_csv(out2, index=False)
print(f"  → {out2.name} ({len(res)} tests)")
if len(res):
    sig = res[(res["metric"] == "mean_burden") & (res["padj"] < 0.1)]
    print("  BH<0.1 (mean_burden):", len(sig), "arms")
    print(res.head(10).to_string(index=False))

# ─────────────────────────────────────────────────────────────────────────────
# GSE281759 primary vs recurrence
# Sample naming: *_T1P (primary) vs *_T1R (recurrence) etc; "P" / "R" suffixes.
# ─────────────────────────────────────────────────────────────────────────────
def primary_recur(sid):
    m = re.search(r"T\d+([PR])$", sid)
    if not m:
        return None
    return "Primary" if m.group(1) == "P" else "Recurrence"

p_only = per_pat[per_pat["dataset"] == "GSE281759"].copy()
p_only["category"] = p_only["sample_id"].apply(primary_recur)
p_only = p_only[p_only["category"].notna()]
rows = []
for arm, sub in p_only.groupby("arm"):
    pr = sub[sub["category"] == "Primary"]
    re_ = sub[sub["category"] == "Recurrence"]
    if len(pr) < 2 or len(re_) < 2:
        continue
    for metric in ("mean_amp", "mean_burden"):
        try:
            U, p = mannwhitneyu(pr[metric].values, re_[metric].values, alternative="two-sided")
        except ValueError:
            continue
        rows.append(dict(
            arm=arm, metric=metric,
            n_P=len(pr), n_R=len(re_),
            median_P=float(pr[metric].median()),
            median_R=float(re_[metric].median()),
            delta=float(pr[metric].median() - re_[metric].median()),
            U=float(U), p=float(p),
        ))
res2 = pd.DataFrame(rows)
if len(res2):
    for metric in res2["metric"].unique():
        m = res2["metric"] == metric
        res2.loc[m, "padj"] = bh(res2.loc[m, "p"].values)
    res2 = res2.sort_values(["metric", "padj"]).reset_index(drop=True)
out3 = TAB_DIR / "spatial_arm_burden_PvsR_GSE281759.csv"
res2.to_csv(out3, index=False)
print(f"  → {out3.name} ({len(res2)} tests)")

# ─────────────────────────────────────────────────────────────────────────────
# Figure: forest of -log10(p) per arm for R vs NR (mean_burden)
# ─────────────────────────────────────────────────────────────────────────────
arms_only = pd.DataFrame()
if len(res) and "metric" in res.columns:
    arms_only = res[(res["metric"] == "mean_burden") & res["arm"].str.match(r"^chr.+[pq]$")].copy()
if len(arms_only):
    arms_only["chr_n"] = arms_only["arm"].str.extract(r"chr(\d+|X|Y)")
    chr_order = [str(i) for i in range(1, 23)] + ["X", "Y"]
    arms_only["chr_rank"] = arms_only["chr_n"].apply(lambda x: chr_order.index(x) if x in chr_order else 99)
    arms_only["arm_kind"] = arms_only["arm"].str[-1]
    arms_only = arms_only.sort_values(["chr_rank", "arm_kind"]).reset_index(drop=True)

    fig, axes = plt.subplots(1, 2, figsize=(13, 8), gridspec_kw={"width_ratios": [1.4, 1]})
    y = np.arange(len(arms_only))
    colors = ["#2E86C1" if d > 0 else "#C0392B" for d in arms_only["delta"]]

    ax = axes[0]
    ax.barh(y, arms_only["delta"], color=colors, edgecolor="black", linewidth=0.4)
    ax.set_yticks(y); ax.set_yticklabels(arms_only["arm"], fontsize=7)
    ax.axvline(0, color="black", linewidth=0.6)
    ax.set_xlabel("Δ median CNV burden (R − NR)")
    ax.set_title("Per-arm CNV burden difference")
    ax.invert_yaxis()

    ax = axes[1]
    neglog = -np.log10(arms_only["p"].clip(lower=1e-10))
    ax.scatter(neglog, y, c=colors, edgecolor="black", linewidth=0.3, s=22)
    ax.axvline(-np.log10(0.05), color="grey", linestyle="--", linewidth=0.7, label="p=0.05")
    sig_padj = arms_only["padj"] < 0.1
    ax.scatter(neglog[sig_padj.values], y[sig_padj.values],
               s=80, facecolor="none", edgecolor="black", linewidth=1.0, label="BH<0.1")
    ax.set_yticks([])
    ax.set_xlabel("−log10 p (MWU)")
    ax.set_title("Significance")
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    ax.invert_yaxis()

    fig.suptitle(
        "Spatial CNV burden: Responder vs Non-Responder (per chromosome arm, tumor spots)",
        fontsize=11, y=0.995,
    )
    fig.tight_layout()
    png = FIG_DIR / "spatial_arm_burden_RvsNR.png"
    pdf = FIG_DIR / "spatial_arm_burden_RvsNR.pdf"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {png.name}")
    print(f"  → {pdf.name}")
print("Done.")
