#!/usr/bin/env python3
"""
05z_patient_strips.py — Per-patient strip/violin plots for L3 composition,
L4 IO signatures, and L5 NMF programs.

Statistical design
------------------
Each cohort is analysed independently.

Discovery (GSE206325, n_R=7, n_NR=17):
  - Mann-Whitney U, BH-corrected across features within each layer.
  - Effect size: rank-biserial r = 2U/(n_R × n_NR) - 1.
  - Reported in annotation: stars + r above each strip.

Cappuyns2023 (n_R=19, n_NR=17):
  - Same statistical procedure, run independently on this cohort.

Inputs
------
  out/05_multilayer/tables/patient_multilayer_features.csv
  out/05_multilayer/tables/nmf_programs_patient_activity.csv

Outputs → out/05_multilayer/figures/05z_expression_spatial/
  05z_patient_strips_L3_Composition_discovery.{png,pdf}
  05z_patient_strips_L3_Composition_cappuyns.{png,pdf}
  05z_patient_strips_L4_IO_Signatures_discovery.{png,pdf}
  05z_patient_strips_L4_IO_Signatures_cappuyns.{png,pdf}
  05z_patient_strips_L5_NMF_Programs_discovery.{png,pdf}
  05z_patient_strips_L5_NMF_Programs_cappuyns.{png,pdf}
"""
from __future__ import annotations
import warnings
warnings.filterwarnings("ignore")

import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from config.publication_style import apply_style, PAL

apply_style()

TAB = PROJECT / "out" / "05_multilayer" / "tables"
OUT = PROJECT / "out" / "05_multilayer" / "figures" / "05z_expression_spatial"
OUT.mkdir(parents=True, exist_ok=True)

# ── Palettes ───────────────────────────────────────────────────────────
R_COL  = PAL["Responder"]       # blue
NR_COL = PAL["Non_Responder"]   # red
RESP_ORDER = ["Responder", "Non-Responder"]

PT_MARKER = "o"   # patient jitter marker (same for both cohorts)

# ── Feature definitions ─────────────────────────────────────────────────
L3_FEATURES = [
    ("frac_CD8_T",    "CD8+ T"),
    ("frac_Macrophage","Macrophage"),
    ("frac_NK",       "NK"),
    ("frac_Treg",     "Treg"),
    ("frac_Plasma",   "Plasma / B"),
    ("frac_B_cell",   "B cell"),
    ("frac_cDC2",     "cDC2"),
    ("frac_cDC1",     "cDC1"),
    ("frac_mregDC",   "mregDC"),
]

L4_FEATURES = [
    ("sig_T_cell_exhaustion",   "T-cell Exhaust."),
    ("sig_T_cell_activation",   "T-cell Activat."),
    ("sig_T_cell_cytotoxicity", "T-cell Cytotox."),
    ("sig_IFN_gamma",           "IFN-γ"),
    ("sig_Immune_checkpoint",   "Immune Checkpoint"),
    ("sig_TGF_beta",            "TGF-β"),
    ("sig_Angiogenesis",        "Angiogenesis"),
    ("sig_Tfh_like",            "Tfh-like"),
    ("sig_Antigen_presentation","Antigen Present."),
]


# ── Helpers ────────────────────────────────────────────────────────────
def rank_biserial(u_stat: float, n1: int, n2: int) -> float:
    """Rank-biserial correlation from MWU U-statistic."""
    return float(2 * u_stat / (n1 * n2) - 1)


def sig_stars(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def strip_panel(
    ax: plt.Axes,
    r_vals: np.ndarray,
    nr_vals: np.ndarray,
    label: str,
    p_bh: float,
    r_rb: float,
) -> None:
    """Draw one feature panel: half-violin + jitter + BH-corrected annotation."""
    is_sig = (not np.isnan(p_bh)) and p_bh < 0.05

    # Highlight significant panels with a pale amber background
    if is_sig:
        ax.set_facecolor("#fff8e1")

    for vals, xpos, color in [
        (r_vals,  0, R_COL),
        (nr_vals, 1, NR_COL),
    ]:
        if len(vals) == 0:
            continue
        parts = ax.violinplot(
            vals, positions=[xpos], widths=0.7,
            showmeans=False, showmedians=False, showextrema=False,
        )
        for pc in parts["bodies"]:
            pc.set_facecolor(color)
            pc.set_alpha(0.35)
            pc.set_edgecolor("none")

        # median line
        med = np.median(vals)
        ax.plot([xpos - 0.15, xpos + 0.15], [med, med],
                color=color, lw=1.5, zorder=4)

        # jitter
        rng = np.random.default_rng(seed=42 + xpos)
        jitter = rng.uniform(-0.12, 0.12, len(vals))
        ax.scatter(
            np.full(len(vals), xpos) + jitter, vals,
            marker=PT_MARKER, s=28, color=color, alpha=0.9,
            linewidths=0.4, edgecolors="white", zorder=5,
        )

    # Significance annotation
    all_vals = np.concatenate([v for v in (r_vals, nr_vals) if len(v) > 0])
    if len(all_vals) == 0:
        ax.axis("off")
        return
    y_max   = np.max(all_vals)
    y_min   = np.min(all_vals)
    y_range = max(y_max - y_min, 1e-9)

    bar_y2  = y_max + y_range * 0.08   # bracket tick ends
    bar_y   = y_max + y_range * 0.14   # bracket horizontal bar
    star_y  = bar_y + y_range * 0.02   # stars/ns text bottom
    r_y     = bar_y + y_range * 0.14   # r-value text bottom
    p_y     = bar_y + y_range * 0.26   # BH p-value text (sig only)

    # Extend y-axis to accommodate stacked annotation labels
    y_top = y_max + y_range * (0.68 if is_sig else 0.52)
    ax.set_ylim(y_min - y_range * 0.05, y_top)

    ax.plot([0, 0, 1, 1], [bar_y2, bar_y, bar_y, bar_y2],
            lw=0.9, color="0.30" if is_sig else "0.45")

    stars     = sig_stars(p_bh)
    sig_color = "#c0392b" if is_sig else "0.55"

    # Stars (or "ns") — larger and bold when significant
    ax.text(0.5, star_y, stars,
            ha="center", va="bottom",
            fontsize=10 if is_sig else 8,
            fontweight="bold" if is_sig else "normal",
            fontstyle="normal" if is_sig else "italic",
            color=sig_color,
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])

    # Rank-biserial correlation (separate, smaller)
    ax.text(0.5, r_y, f"r = {r_rb:+.2f}",
            ha="center", va="bottom",
            fontsize=7,
            color=sig_color if is_sig else "0.50",
            fontstyle="italic",
            path_effects=[pe.withStroke(linewidth=1.5, foreground="white")])

    # BH-corrected p-value label — only for significant panels
    if is_sig:
        p_str = f"p_BH < 0.001" if p_bh < 0.001 else f"p_BH = {p_bh:.3f}"
        ax.text(0.5, p_y, p_str,
                ha="center", va="bottom",
                fontsize=6.5,
                fontweight="semibold",
                color=sig_color,
                path_effects=[pe.withStroke(linewidth=1.5, foreground="white")])

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["R", "NR"], fontsize=9, fontweight="semibold")
    ax.set_title(label, fontsize=9.5, fontweight="bold", pad=4)
    ax.tick_params(axis="y", labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)

    # Colored border frame for significant panels
    if is_sig:
        for sname in ["left", "bottom"]:
            ax.spines[sname].set_edgecolor("#e67e22")
            ax.spines[sname].set_linewidth(1.8)


def make_cohort_figure(
    layer_name: str,
    features: list[tuple[str, str]],
    df: pd.DataFrame,
    cohort_label: str,
    ncols: int = 5,
) -> plt.Figure:
    """Build a multi-panel figure for one layer, one cohort."""
    n = len(features)
    nrows = int(np.ceil(n / ncols))
    fig_w = ncols * 2.2 + 0.6
    fig_h = nrows * 3.2 + 1.2

    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h))
    axes = np.array(axes).flatten()

    # Run MWU + BH within this cohort for this layer
    p_raw      = []
    stats_store = []
    for col, _ in features:
        if col not in df.columns:
            p_raw.append(np.nan)
            stats_store.append((np.nan, np.nan))
            continue
        r_v  = df.loc[df["response"] == "Responder",     col].dropna().values
        nr_v = df.loc[df["response"] == "Non-Responder", col].dropna().values
        if len(r_v) < 2 or len(nr_v) < 2:
            p_raw.append(np.nan)
            stats_store.append((np.nan, np.nan))
        else:
            u, p = mannwhitneyu(r_v, nr_v, alternative="two-sided")
            rb   = rank_biserial(u, len(r_v), len(nr_v))
            p_raw.append(p)
            stats_store.append((p, rb))

    # BH correction
    valid_mask = [not np.isnan(p) for p in p_raw]
    p_bh_full  = [np.nan] * n
    if any(valid_mask):
        p_valid = [p_raw[i] for i in range(n) if valid_mask[i]]
        _, p_corr, _, _ = multipletests(p_valid, method="fdr_bh")
        j = 0
        for i in range(n):
            if valid_mask[i]:
                p_bh_full[i] = p_corr[j]
                j += 1

    for idx, (col, label) in enumerate(features):
        ax = axes[idx]
        if col not in df.columns:
            ax.axis("off")
            continue

        r_v  = df.loc[df["response"] == "Responder",     col].dropna().values
        nr_v = df.loc[df["response"] == "Non-Responder", col].dropna().values

        p_bh = p_bh_full[idx] if not np.isnan(p_bh_full[idx]) else 1.0
        rb   = stats_store[idx][1] if not np.isnan(stats_store[idx][1]) else 0.0

        strip_panel(ax, r_v, nr_v, label, p_bh, rb)

    # Hide unused axes
    for idx in range(n, len(axes)):
        axes[idx].axis("off")

    # Figure-level legend
    n_r  = (df["response"] == "Responder").sum()
    n_nr = (df["response"] == "Non-Responder").sum()
    legend_elements = [
        Line2D([0], [0], marker=PT_MARKER, color="w",
               markerfacecolor=R_COL,  markersize=7, label=f"Responder (n={n_r})"),
        Line2D([0], [0], marker=PT_MARKER, color="w",
               markerfacecolor=NR_COL, markersize=7, label=f"Non-Responder (n={n_nr})"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=2,
        fontsize=7.5,
        frameon=True,
        bbox_to_anchor=(0.5, 0.01),
    )

    fig.suptitle(
        f"{layer_name}  ·  {cohort_label}"
        f"\nR n={n_r}   NR n={n_nr}   |   Mann–Whitney U + BH correction",
        fontsize=10, fontweight="bold", y=1.02,
        linespacing=1.6,
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.97])
    return fig


def save_fig(fig: plt.Figure, stem: str) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=300, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"  saved  {stem}")


# ── Load data ──────────────────────────────────────────────────────────
feats = pd.read_csv(TAB / "patient_multilayer_features.csv")
nmf   = pd.read_csv(TAB / "nmf_programs_patient_activity.csv")

# Merge NMF into feats table
nmf_cols = [c for c in nmf.columns if c.startswith("NMF")]
feats = feats.merge(nmf[["patient_id"] + nmf_cols], on="patient_id", how="left")

disc = feats.loc[feats["dataset"].astype(str) != "Cappuyns2023"].copy()
capp = feats.loc[feats["dataset"].astype(str) == "Cappuyns2023"].copy()

np.random.seed(42)

L5_FEATURES = [(f"NMF{i}", f"NMF {i}") for i in range(1, 11)]

COHORTS = [
    (disc, "GSE206325 (discovery)",  "discovery"),
    (capp, "Cappuyns2023 (validation)", "cappuyns"),
]

# ── L3 Composition ─────────────────────────────────────────────────────
for df, cohort_label, cohort_key in COHORTS:
    print(f"L3 Composition — {cohort_label} ...")
    fig = make_cohort_figure("L3  Cell composition", L3_FEATURES, df, cohort_label, ncols=5)
    save_fig(fig, f"05z_patient_strips_L3_Composition_{cohort_key}")

# ── L4 IO Signatures ───────────────────────────────────────────────────
for df, cohort_label, cohort_key in COHORTS:
    print(f"L4 IO Signatures — {cohort_label} ...")
    fig = make_cohort_figure("L4  IO signatures", L4_FEATURES, df, cohort_label, ncols=5)
    save_fig(fig, f"05z_patient_strips_L4_IO_Signatures_{cohort_key}")

# ── L5 NMF Programs ────────────────────────────────────────────────────
for df, cohort_label, cohort_key in COHORTS:
    print(f"L5 NMF Programs — {cohort_label} ...")
    fig = make_cohort_figure("L5  NMF programs", L5_FEATURES, df, cohort_label, ncols=5)
    save_fig(fig, f"05z_patient_strips_L5_NMF_Programs_{cohort_key}")

print(f"\nDone. Output: {OUT}")
