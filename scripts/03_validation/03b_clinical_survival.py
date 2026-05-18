#!/usr/bin/env python3
"""
03b_clinical_survival.py — Layer-stratified survival analysis (IMC validation cohort)

Design
------
Survival data is available only for the IMC validation cohort
(Salie et al. Gut 2025; n=42; n=35 survival-eligible, Child-Pugh < 7).

The scRNA-seq discovery cohort (GSE206325, n=24) and external validation
cohort (Cappuyns 2023, n=32) lack OS/PFS data and are therefore NOT
APPLICABLE for survival stratification.

Per-layer analysis:
  L3 (Composition):  IMC CD8 density (continuous → median split)
  L4 (IO functional): CD8_tex_like_frac (exhaustion proxy, median split)
                       CD8_gzmb_frac     (cytotoxicity proxy, median split)
  L1 / L2 / L5:     Not applicable — CNV, regulon, and NMF features are not
                     available in IMC imaging data.

Outputs
-------
  out/03_validation/figures/03b_clinical_survival_validation.{png,pdf}
  out/03_validation/tables/03b_layer_survival_km_summary.csv
"""

import warnings
warnings.filterwarnings("ignore")

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy.stats import norm

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from config.publication_style import apply_style, PAL, LAYER_PAL, save_fig
apply_style()

OUT  = PROJECT / "out" / "03_validation"
FIG  = OUT / "figures"
TAB  = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

IMC_PATH = (PROJECT / "out" / "03_validation" / "imc_validation"
            / "tables" / "imc_patient_features.csv")

# ── Layer definitions for survival stratification ─────────────────────
# Each entry: (feature_col, layer_label, panel_label, color_key)
LAYERS = [
    ("CD8_density_mm2",    "L3 · Composition",
     "CD8 density\n(cells/mm²)",      "L3_composition"),
    ("CD8_tex_like_frac",  "L4 · IO functional",
     "CD8 exhaustion-like\n(fraction PD-1⁺TOX⁺)",  "L4_IO_signature"),
    ("CD8_gzmb_frac",      "L4 · IO functional",
     "CD8 cytotoxicity\n(fraction GzmB⁺)",           "L4_IO_signature"),
]

NOT_APPLICABLE = {
    "L1 (CNV)":   "Genomic CNV features not available in IMC data.",
    "L2 (Regulon)": "Transcription-factor regulon axes not available in IMC data.",
    "L5 (NMF)":   "NMF gene programs not derivable from protein-level IMC data.",
    "Discovery (GSE206325)": "No OS/PFS data available for scRNA-seq discovery cohort.",
    "Validation (Cappuyns 2023)": "No OS/PFS data available for Cappuyns atezo+bev cohort.",
}

HI_COLOR = "#0077BB"
LO_COLOR = "#CC3311"


# ── Helper ────────────────────────────────────────────────────────────
def median_split(series: pd.Series) -> pd.Series:
    """Return 'High' / 'Low' based on the median; ties assigned to 'High'."""
    med = series.median()
    return series.apply(lambda v: "High" if v >= med else "Low")


def km_panel(ax: plt.Axes, df: pd.DataFrame,
             duration_col: str, event_col: str,
             group_col: str,
             endpoint_label: str,
             layer_color: str) -> dict:
    """
    Draw a single KM panel on ax.

    Returns a dict with summary statistics for the table.
    """
    hi = df[df[group_col] == "High"]
    lo = df[df[group_col] == "Low"]

    kmf_hi = KaplanMeierFitter()
    kmf_lo = KaplanMeierFitter()

    n_hi = len(hi); n_lo = len(lo)

    kmf_hi.fit(hi[duration_col], hi[event_col], label=f"High (n={n_hi})")
    kmf_lo.fit(lo[duration_col], lo[event_col], label=f"Low (n={n_lo})")

    kmf_hi.plot_survival_function(
        ax=ax, ci_show=True,
        color=HI_COLOR,
        ci_alpha=0.12,
        at_risk_counts=False,
    )
    kmf_lo.plot_survival_function(
        ax=ax, ci_show=True,
        color=LO_COLOR,
        ci_alpha=0.12,
        at_risk_counts=False,
    )

    # Censoring marks
    for kmf, col in [(kmf_hi, HI_COLOR), (kmf_lo, LO_COLOR)]:
        events = kmf.event_table
        cens_times = events.index[events["censored"] > 0]
        for t in cens_times:
            sf_val = kmf.predict(t)
            ax.plot(t, sf_val, "|", color=col, ms=5, mew=1.2, zorder=5)

    # Log-rank test
    results = logrank_test(
        hi[duration_col], lo[duration_col],
        event_observed_A=hi[event_col],
        event_observed_B=lo[event_col],
    )
    pval = results.p_value

    if pval < 0.001:
        pstr = "p<0.001"
    else:
        pstr = f"p={pval:.3f}"

    ax.text(
        0.97, 0.97, f"Log-rank {pstr}",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=7.5, fontstyle="italic",
        color="#c0392b" if pval < 0.05 else "0.45",
        fontweight="bold" if pval < 0.05 else "normal",
    )

    ax.set_ylim(-0.04, 1.09)
    ax.set_ylabel(endpoint_label, fontsize=8.5)
    ax.set_xlabel("Time (months)", fontsize=8)
    sns.despine(ax=ax)

    # Legend — place lower-left to avoid overlap with p-value annotation
    handles = [
        mpatches.Patch(facecolor=HI_COLOR, label=f"High (n={n_hi})"),
        mpatches.Patch(facecolor=LO_COLOR, label=f"Low (n={n_lo})"),
    ]
    ax.legend(handles=handles, fontsize=7, loc="lower left",
              frameon=False, handlelength=1.1, handleheight=0.9)

    med_hi = kmf_hi.median_survival_time_
    med_lo = kmf_lo.median_survival_time_

    return {
        "n_high": n_hi, "n_low": n_lo,
        "median_high": med_hi, "median_low": med_lo,
        "logrank_p": pval,
    }


# ── Main ──────────────────────────────────────────────────────────────
def run():
    print("[03b] Layer-stratified survival analysis — IMC validation cohort")

    # ── Load data ─────────────────────────────────────────────────────
    df_raw = pd.read_csv(IMC_PATH)

    # Censoring convention (from imc_05): censored_0yes=0 → event occurred
    df_raw["event_os"]  = (df_raw["censored_0yes"] == 0).astype(int)
    df_raw["event_pfs"] = (df_raw["censored_0yes"] == 0).astype(int)

    surv = df_raw[df_raw["survival_eligible"] == 1].copy()
    n_elig = len(surv)
    print(f"  Survival-eligible: {n_elig}  (R={surv.response_binary.sum()}, "
          f"NR={n_elig - surv.response_binary.sum()})")

    # ── Median-split all layer features ───────────────────────────────
    for feat, *_ in LAYERS:
        surv[f"{feat}_group"] = median_split(surv[feat])

    # ── Build figure: 2 endpoints × 3 layer features ─────────────────
    n_layers = len(LAYERS)
    n_endpoints = 2
    fig, axes = plt.subplots(
        n_endpoints, n_layers,
        figsize=(n_layers * 3.4, n_endpoints * 3.6),
        squeeze=False,
    )

    endpoints = [
        ("os_months",  "event_os",  "Overall survival probability"),
        ("pfs_months", "event_pfs", "Progression-free survival probability"),
    ]
    endpoint_labels = ["OS", "PFS"]

    km_records = []

    for row_i, (dur_col, evt_col, ep_label) in enumerate(endpoints):
        for col_i, (feat, layer_lbl, panel_lbl, color_key) in enumerate(LAYERS):
            ax = axes[row_i][col_i]
            grp_col = f"{feat}_group"
            stats = km_panel(
                ax, surv, dur_col, evt_col, grp_col,
                endpoint_label=ep_label,
                layer_color=LAYER_PAL.get(color_key, "#555555"),
            )
            # Column title (top row only)
            if row_i == 0:
                ax.set_title(
                    f"{panel_lbl}\n({layer_lbl})",
                    fontsize=8.5, fontweight="bold", pad=5, linespacing=1.5,
                )
            # Row label — only the endpoint abbreviation on left, full label from km_panel
            ax.set_ylabel(ep_label, fontsize=8)

            # Panel letter
            letter = chr(ord("A") + row_i * n_layers + col_i)
            ax.text(-0.18, 1.02, letter, transform=ax.transAxes,
                    fontsize=12, fontweight="bold", va="bottom")

            km_records.append({
                "endpoint":     endpoint_labels[row_i],
                "layer":        layer_lbl,
                "feature":      feat,
                "split":        "median",
                "n_high":       stats["n_high"],
                "n_low":        stats["n_low"],
                "median_high":  stats["median_high"],
                "median_low":   stats["median_low"],
                "logrank_p":    stats["logrank_p"],
            })

    fig.suptitle(
        f"Layer-stratified Kaplan–Meier survival  ·  IMC validation cohort  "
        f"(n={n_elig} survival-eligible)\n"
        "Discovery (GSE206325) and Cappuyns validation cohorts lack OS/PFS data "
        "(not applicable).  "
        "L1 (CNV), L2 (Regulon), L5 (NMF) not available in IMC imaging data.",
        fontsize=8, y=1.03, linespacing=1.6,
    )
    fig.tight_layout(h_pad=3.0, w_pad=2.0)

    # ── Save figure ───────────────────────────────────────────────────
    out_base = FIG / "03b_clinical_survival_validation"
    save_fig(fig, out_base)
    plt.close(fig)
    print(f"  Saved: {out_base}.png / .pdf")

    # ── Save table ────────────────────────────────────────────────────
    km_df = pd.DataFrame(km_records)
    # Flag significance
    km_df["sig_nominal"] = km_df["logrank_p"] < 0.05
    km_df["sig_label"] = km_df["logrank_p"].apply(
        lambda p: "***" if p < 0.001 else ("**" if p < 0.01 else
                  ("*"   if p < 0.05 else "ns"))
    )
    tab_path = TAB / "03b_layer_survival_km_summary.csv"
    km_df.to_csv(tab_path, index=False)
    print(f"  Saved table: {tab_path}")

    # ── Print summary ─────────────────────────────────────────────────
    print("\n  KM summary:")
    for _, row in km_df.iterrows():
        print(f"    {row.endpoint:3s}  {row.layer:22s}  {row.feature:24s}  "
              f"Hi_med={row.median_high:.1f}  Lo_med={row.median_low:.1f}  "
              f"p={row.logrank_p:.4f}  {row.sig_label}")

    # ── Print not-applicable note ──────────────────────────────────────
    print("\n  NOT APPLICABLE (no survival data or feature unavailable):")
    for key, reason in NOT_APPLICABLE.items():
        print(f"    {key}: {reason}")


if __name__ == "__main__":
    run()
