#!/usr/bin/env python3
"""04e – TCR clonotype dynamics from real GSE206325 TCR data.

Loads paired TCR repertoire data, joins with atlas metadata for tissue
type and response, then plots:
  A. Clonal expansion bar chart by tissue (top expanded clonotypes)
  B. Diversity metrics (Shannon / clonality) by tissue and response
"""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
ATLAS   = PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad"
TCR_CSV = PROJECT / "data" / "processed" / "GSE206325_tcr.csv"
OUT = PROJECT / "out" / "04_complementary"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def _load_tcr_with_metadata():
    """Load TCR repertoire and map subject_id → response from atlas."""
    tcr = pd.read_csv(TCR_CSV)

    # Load atlas metadata for patient → response mapping
    print("  Loading atlas metadata …")
    adata = sc.read_h5ad(str(ATLAS), backed="r")
    obs = adata.obs

    # Map TCR subject_id (e.g. 840001001) to atlas patient_ID (e.g. 1001.0)
    # by extracting the last 4 digits
    gse = obs[obs["dataset"] == "GSE206325"] if "dataset" in obs.columns else obs

    # Build patient → response
    resp_map = {}
    if "patient_ID" in gse.columns and "response" in gse.columns:
        for pid, grp in gse.groupby("patient_ID"):
            pid_str = str(int(float(pid))) if not pd.isna(pid) else None
            if pid_str:
                resp_map[pid_str] = grp["response"].astype(str).value_counts().index[0]

    # Map subject_id → patient → response
    tcr["patient_short"] = tcr["subject_id"].astype(str).str[-4:]
    tcr["response"] = tcr["patient_short"].map(resp_map).fillna("Unknown")

    # Use sample_id to infer per-sample clonal structure
    # (each sample = a tissue site, but we don't have explicit tissue labels in TCR file)
    print(f"  TCR records: {len(tcr)}, response: {tcr['response'].value_counts().to_dict()}")
    print(f"  Patients: {tcr['subject_id'].nunique()}, Samples: {tcr['sample_id'].nunique()}")
    return tcr


def plot_tcr_clonality():
    print("[04e] Generating TCR Clonotype Dynamics from real data …")
    tcr = _load_tcr_with_metadata()

    # Focus on TRB chain (more informative for clonality)
    trb = tcr[tcr["TR_chain"] == "TRB"].copy()
    if len(trb) < 100:
        trb = tcr.copy()

    # ── Per-sample expansion categories ──────────────────────────────
    expansion_rows = []
    for sid, grp in trb.groupby("sample_id"):
        top = grp.groupby("CDR3.aa")["Clones"].sum().sort_values(ascending=False)
        n = len(top)
        if n < 10:
            continue
        q99 = top.quantile(0.99)
        q95 = top.quantile(0.95)
        q80 = top.quantile(0.80)
        hyper  = (top > q99).sum() / n
        large  = ((top > q95) & (top <= q99)).sum() / n
        medium = ((top > q80) & (top <= q95)).sum() / n
        small  = 1.0 - hyper - large - medium
        resp = grp["response"].iloc[0]
        expansion_rows.append({
            "sample_id": sid, "response": resp,
            "Hyperexpanded": hyper, "Large": large,
            "Medium": medium, "Small/Singleton": small,
        })
    exp_df = pd.DataFrame(expansion_rows)

    # ── Per-sample diversity ─────────────────────────────────────────
    diversity_rows = []
    for sid, grp in trb.groupby("sample_id"):
        props = grp["Clones"] / grp["Clones"].sum()
        props = props[props > 0]
        shannon = -np.sum(props * np.log2(props))
        n_unique = len(props)
        clonality = 1 - (shannon / np.log2(n_unique)) if n_unique > 1 else 0
        resp = grp["response"].iloc[0]
        diversity_rows.append({
            "sample_id": sid, "response": resp,
            "shannon": shannon, "clonality": clonality,
            "n_clonotypes": n_unique,
        })
    div_df = pd.DataFrame(diversity_rows)

    # ── Figure: 3-panel ──────────────────────────────────────────────
    plt.style.use("seaborn-v0_8-whitegrid")
    fig = plt.figure(figsize=(20, 7))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.5, 1, 1], wspace=0.30)

    colors_exp = ["#c0392b", "#e67e22", "#f1c40f", "#95a5a6"]
    cats = ["Hyperexpanded", "Large", "Medium", "Small/Singleton"]
    resp_palette = {"Responder": "#2ecc71", "Non-Responder": "#e74c3c", "Unknown": "#95a5a6"}

    # Panel A: stacked expansion by sample, coloured by response status
    axA = fig.add_subplot(gs[0, 0])
    exp_sorted = exp_df.sort_values(["response", "Hyperexpanded"], ascending=[True, False])
    x_positions = np.arange(len(exp_sorted))
    bottom = np.zeros(len(exp_sorted))
    for cat, color in zip(cats, colors_exp):
        vals = exp_sorted[cat].values
        axA.bar(x_positions, vals, bottom=bottom, color=color,
                label=cat, edgecolor="white", lw=0.3, width=0.9)
        bottom += vals
    # Colour x-axis labels or add bar at bottom for response
    for i, (_, row) in enumerate(exp_sorted.iterrows()):
        rc = resp_palette.get(row["response"], "#95a5a6")
        axA.plot(i, -0.02, "s", color=rc, markersize=5, clip_on=False)
    axA.set_ylabel("Fraction of Unique Clonotypes", fontweight="bold")
    axA.set_title("A. Clonal Expansion per Sample", fontweight="bold", fontsize=14, pad=10)
    axA.set_xticks(x_positions)
    axA.set_xticklabels(exp_sorted["sample_id"].values, rotation=60, ha="right", fontsize=7)
    axA.legend(loc="upper right", fontsize=8, frameon=True, facecolor="white")
    axA.set_ylim(-0.04, 1.05)
    for sp in ("top", "right"):
        axA.spines[sp].set_visible(False)

    # Panel B: Shannon diversity by response
    axB = fig.add_subplot(gs[0, 1])
    resp_order = [r for r in ["Responder", "Non-Responder"]
                  if r in div_df["response"].unique()]
    for j, resp in enumerate(resp_order):
        sub = div_df[div_df["response"] == resp]
        c = resp_palette.get(resp, "#95a5a6")
        jitter = np.random.RandomState(42).uniform(-0.15, 0.15, len(sub))
        axB.scatter(j + jitter, sub["shannon"], color=c, s=60,
                    edgecolor="black", lw=0.5, alpha=0.8, zorder=3, label=resp)
        axB.errorbar(j, sub["shannon"].mean(), yerr=sub["shannon"].sem(),
                     fmt="D", color=c, markersize=10, markeredgecolor="black",
                     capsize=5, capthick=1.5, lw=1.5, zorder=4)
    axB.set_xticks(range(len(resp_order)))
    axB.set_xticklabels(resp_order, fontweight="bold")
    axB.set_ylabel("Shannon Diversity (bits)", fontweight="bold")
    axB.set_title("B. TCR Diversity by IO Response", fontweight="bold", fontsize=14, pad=10)
    for sp in ("top", "right"):
        axB.spines[sp].set_visible(False)

    # Panel C: Clonality by response
    axC = fig.add_subplot(gs[0, 2])
    for j, resp in enumerate(resp_order):
        sub = div_df[div_df["response"] == resp]
        c = resp_palette.get(resp, "#95a5a6")
        jitter = np.random.RandomState(43).uniform(-0.15, 0.15, len(sub))
        axC.scatter(j + jitter, sub["clonality"], color=c, s=60,
                    edgecolor="black", lw=0.5, alpha=0.8, zorder=3)
        axC.errorbar(j, sub["clonality"].mean(), yerr=sub["clonality"].sem(),
                     fmt="D", color=c, markersize=10, markeredgecolor="black",
                     capsize=5, capthick=1.5, lw=1.5, zorder=4)
    # Mann-Whitney between groups
    if len(resp_order) == 2:
        from scipy.stats import mannwhitneyu
        g1 = div_df[div_df["response"] == resp_order[0]]["clonality"]
        g2 = div_df[div_df["response"] == resp_order[1]]["clonality"]
        if len(g1) > 1 and len(g2) > 1:
            _, pval = mannwhitneyu(g1, g2, alternative="two-sided")
            axC.text(0.5, 0.95, f"p = {pval:.3g}", transform=axC.transAxes,
                     ha="center", fontsize=11, fontweight="bold",
                     bbox=dict(fc="white", ec="black", pad=3, boxstyle="round"))
    axC.set_xticks(range(len(resp_order)))
    axC.set_xticklabels(resp_order, fontweight="bold")
    axC.set_ylabel("Clonality (1 − norm. entropy)", fontweight="bold")
    axC.set_title("C. Clonality by IO Response", fontweight="bold", fontsize=14, pad=10)
    for sp in ("top", "right"):
        axC.spines[sp].set_visible(False)

    plt.tight_layout()
    plt.savefig(FIG / "04e_tcr_clonotype_alluvial.pdf", dpi=300, bbox_inches="tight")
    plt.savefig(FIG / "04e_tcr_clonotype_alluvial.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {FIG / '04e_tcr_clonotype_alluvial.png'}")


if __name__ == "__main__":
    plot_tcr_clonality()
