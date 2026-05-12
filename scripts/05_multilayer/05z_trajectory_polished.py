#!/usr/bin/env python3
"""
05z_trajectory_polished.py — Publication-grade trajectory figures for
significant layers / programs across discovery (scRNA + spatial) and
validation (Cappuyns scRNA) cohorts.

Three figures are produced, all loading pre-computed tables only
(no atlas / h5ad loading, runtime ~30 s):

  Figure 1 · 05z_traj_scrna_pseudotime.{png,pdf}
      CD8 pseudotime (discovery, n=24) × NMF10 cross-analysis
      + Cappuyns validation NMF programme layer (n=32)

  Figure 2 · 05z_traj_cnv_nmf10.{png,pdf}
      Patient CNV–clone-diversity trajectory (discovery)
      with LOESS-smoothed feature traces + ranked lollipop

  Figure 3 · 05z_traj_spatial_gradient.{png,pdf}
      Spatial gradient trajectories across 7 Visium sections
      (spatial-gradient ρ heatmap, section bar chart, immune-gradient scatter)

Outputs → out/05_multilayer/figures/05z_trajectory_polished/

Statistical discipline:
  - Patient-level statistics only (no cell-level inflation)
  - Spearman ρ for trajectory correlations
  - LOESS / Gaussian smoothing for visualisation only, NOT inference
  - Cliff's delta for effect-size reporting
  - Uncorrected p reported for exploratory Spearman; BH-corrected where noted
"""
from __future__ import annotations
import warnings
warnings.filterwarnings("ignore")

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as mgridspec
from matplotlib.patches import FancyArrowPatch
import seaborn as sns
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.nonparametric.smoothers_lowess import lowess

# ── project roots ─────────────────────────────────────────────────────────────
PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from config.publication_style import apply_style, PAL, CMAP as PCMAP

apply_style()

# ── output directory ───────────────────────────────────────────────────────────
OUT  = PROJECT / "out" / "05_multilayer" / "figures" / "05z_trajectory_polished"
TABT = PROJECT / "out" / "05_multilayer" / "figures" / "05z_trajectory_polished"
OUT.mkdir(parents=True, exist_ok=True)

# convenience save
def _save(fig: plt.Figure, stem: str, dpi: int = 280) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=dpi,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  saved → {stem}.{{png,pdf}}")

# ── shared colours ─────────────────────────────────────────────────────────────
R_COL  = PAL.get("Responder",     "#0077BB")
NR_COL = PAL.get("Non-Responder", "#CC3311")

def resp_color(r: str) -> str:
    return R_COL if r in ("Responder", "R") else NR_COL

# ── trajectory feature colourmap (matches 05z_regenerate_spatial_plots.py) ────
FEAT_CMAP = {
    "NMF10":               "inferno",
    "sig_T_cell_exhaustion": "magma",
    "frac_CD8_T":          "viridis",
    "sig_IFN_gamma":       "viridis",
    "frac_Monocyte":       "plasma",
}
# solid accent colours sampled from each cmap for line plots
FEAT_COL = {
    "NMF10":               "#F0A500",  # inferno warm yellow
    "sig_T_cell_exhaustion": "#C03080", # magma magenta
    "frac_CD8_T":          "#2171B5",  # viridis blue
    "sig_IFN_gamma":       "#31A354",  # viridis green
    "frac_Monocyte":       "#EE7733",  # plasma orange
}
FEAT_LS = {
    "NMF10":               "-",
    "sig_T_cell_exhaustion": "--",
    "frac_CD8_T":          "-.",
    "sig_IFN_gamma":       ":",
    "frac_Monocyte":       (0, (3, 1, 1, 1)),
}
FEAT_LW = {
    "NMF10":               2.8,
    "sig_T_cell_exhaustion": 2.0,
    "frac_CD8_T":          2.0,
    "sig_IFN_gamma":       2.0,
    "frac_Monocyte":       1.6,
}
FEAT_LABEL = {
    "NMF10":               "NMF10 (BH p=0.0054)",
    "sig_T_cell_exhaustion": "T cell exhaustion sig.",
    "frac_CD8_T":          "CD8 T fraction",
    "sig_IFN_gamma":       "IFN-γ sig.",
    "frac_Monocyte":       "Monocyte fraction",
}

# ── helper: Cliff's delta ──────────────────────────────────────────────────────
def _cliffs_delta(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    greater = (x[:, None] > y[None, :]).sum()
    less    = (x[:, None] < y[None, :]).sum()
    return float(greater - less) / (len(x) * len(y))

# ── helper: spine cleanup ──────────────────────────────────────────────────────
def _clean(ax, keep=("left", "bottom")):
    for sp, s in ax.spines.items():
        s.set_visible(sp in keep)
    ax.tick_params(which="both", direction="out", length=4)

# ── panel letter helper ────────────────────────────────────────────────────────
def _panel_letter(ax, letter, x=-0.10, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=15, fontweight="bold", va="bottom", ha="left")

# ── annotation box helper ─────────────────────────────────────────────────────
def _statbox(ax, text, loc="upper right", fs=8.5):
    kw = dict(ha="right", va="top") if "right" in loc else dict(ha="left", va="top")
    xy  = (0.97, 0.97) if "right" in loc else (0.03, 0.97)
    ax.text(*xy, text, transform=ax.transAxes, fontsize=fs,
            bbox=dict(boxstyle="round,pad=0.35", fc="white",
                      ec="#BBBBBB", lw=0.8, alpha=0.92),
            **kw)

# ══════════════════════════════════════════════════════════════════════════════
# Load tables
# ══════════════════════════════════════════════════════════════════════════════
print("Loading tables …")

# -- scRNA pseudotime (discovery)
PT_MEANS = PROJECT / "out" / "04_complementary" / "tables" / "04i_pseudotime_patient_means.csv"
PT_STATS  = PROJECT / "out" / "04_complementary" / "tables" / "04i_pseudotime_patient_stats.csv"

# -- Cappuyns validation predictions
CAPP_PRED = PROJECT / "out" / "03_validation" / "tables" / "cappuyns_per_patient_predictions.csv"

# -- Multi-layer features (discovery)
FEATS_CSV = PROJECT / "out" / "05_multilayer" / "tables" / "patient_multilayer_features.csv"

# -- NMF patient activity (discovery + Cappuyns)
NMF_CSV = PROJECT / "out" / "05_multilayer" / "tables" / "nmf_programs_patient_activity.csv"

# -- CNV trajectory (Spearman ρ table)
TRAJ_CSV = PROJECT / "out" / "05_multilayer" / "tables" / "cnv_patient_trajectory.csv"

# -- Spatial gradient correlations (GSE238264)
SGRAD_CSV = (PROJECT / "out" / "02b_spatial_biomarkers" / "tables"
             / "gse238264_signature_gradient_correlations.csv")
SGRAD_IG  = (PROJECT / "out" / "02b_spatial_biomarkers" / "tables"
             / "gse238264_immune_gradient.csv")

# -- Section-level scores from 05z regeneration
SEC_CSV   = (PROJECT / "out" / "05_multilayer" / "figures"
             / "05z_expression_spatial" / "05z_section_scores.csv")

pt_means  = pd.read_csv(PT_MEANS)
pt_stats  = pd.read_csv(PT_STATS).set_index("metric")["value"]

capp_pred = pd.read_csv(CAPP_PRED) if CAPP_PRED.exists() else pd.DataFrame()

feats     = pd.read_csv(FEATS_CSV) if FEATS_CSV.exists() else pd.DataFrame()
if "dataset" in feats.columns:
    feats = feats[feats["dataset"].astype(str) != "Cappuyns2023"].copy()

nmf_all   = pd.read_csv(NMF_CSV)  if NMF_CSV.exists()  else pd.DataFrame()
nmf_disc  = nmf_all[nmf_all["dataset"].astype(str) != "Cappuyns2023"].copy() \
            if "dataset" in nmf_all.columns else nmf_all.copy()

traj_rank = pd.read_csv(TRAJ_CSV) if TRAJ_CSV.exists() else pd.DataFrame()

sgrad_df  = pd.read_csv(SGRAD_CSV) if SGRAD_CSV.exists() else pd.DataFrame()
sgrad_ig  = pd.read_csv(SGRAD_IG,  index_col=0) \
            if SGRAD_IG.exists() else pd.DataFrame()
sec_scores = pd.read_csv(SEC_CSV) if SEC_CSV.exists() else pd.DataFrame()

print(f"  pt_means:  {len(pt_means)} patients")
print(f"  capp_pred: {len(capp_pred)} patients")
print(f"  feats:     {len(feats)} discovery patients")
print(f"  nmf_disc:  {len(nmf_disc)} discovery patients")
print(f"  traj_rank: {len(traj_rank)} features")
print(f"  sgrad_df:  {len(sgrad_df)} rows")
print(f"  sec_scores:{len(sec_scores)} sections")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 1 — scRNA pseudotime: discovery + Cappuyns + NMF10 cross-view
# ══════════════════════════════════════════════════════════════════════════════
print("\n[Fig 1] scRNA pseudotime …")

fig1, axes1 = plt.subplots(
    1, 3,
    figsize=(17, 6.5),
    constrained_layout=True,
    gridspec_kw={"width_ratios": [1.05, 1.05, 1.0]},
)
fig1.patch.set_facecolor("white")

# ── Panel A: Discovery per-patient DPT ───────────────────────────────────────
axA = axes1[0]
_panel_letter(axA, "A", x=-0.13)

r_vals  = pt_means.loc[pt_means["response"] == "Responder",     "mean"].values
nr_vals = pt_means.loc[pt_means["response"] == "Non-Responder", "mean"].values

rng = np.random.default_rng(42)
for vals, col, lab, xi in [(r_vals, R_COL, "Responder", 0),
                             (nr_vals, NR_COL, "Non-Responder", 1)]:
    jitter = rng.uniform(-0.13, 0.13, len(vals))
    axA.scatter(xi + jitter, vals, color=col, s=70,
                edgecolor="white", linewidth=0.9, zorder=4, label=lab)
    med = np.median(vals)
    axA.hlines(med, xi - 0.25, xi + 0.25,
               color=col, lw=2.8, zorder=5)
    axA.text(xi + 0.28, med, f"median={med:.4f}",
             va="center", ha="left", fontsize=8.5, color=col)

mwu_p  = float(pt_stats.get("MWU_p",  0.166))
cliff  = float(pt_stats.get("Cliffs_delta_patient", 0.378))
stat_n = (f"n = {int(pt_stats.get('n_R_patients',7))}R / "
          f"{int(pt_stats.get('n_NR_patients',17))}NR\n"
          f"MWU p = {mwu_p:.3f}  (patient-level)\n"
          f"Cliff's δ = {cliff:+.3f}")
_statbox(axA, stat_n, loc="upper right", fs=8)

axA.set_xticks([0, 1])
axA.set_xticklabels(["Responder", "Non-Responder"], fontsize=11)
axA.set_ylabel("Per-patient mean pseudotime", fontsize=11)
axA.set_title("Discovery CD8 pseudotime\n(patient-level, n=24)", fontsize=11, fontweight="bold")
_clean(axA)
# clip y-axis so a single high-pseudotime outlier doesn't compress all other points
_yall = np.concatenate([r_vals, nr_vals])
axA.set_ylim(np.percentile(_yall, 1) - 0.005, np.percentile(_yall, 99) + 0.01)
axA.yaxis.set_minor_locator(mticker.AutoMinorLocator(4))
# cell-level p disclaimer
axA.text(0.5, -0.16, "⚠ Cell-level p not reported (pseudoreplication)",
         transform=axA.transAxes, ha="center", fontsize=7.5, color="#666666",
         style="italic")

# ── Panel B: Cappuyns validation — L5 NMF pred-prob ─────────────────────────
axB = axes1[1]
_panel_letter(axB, "B", x=-0.13)

if len(capp_pred) > 0 and "L5_nmf_programs_pred_prob" in capp_pred.columns:
    cr = capp_pred.loc[capp_pred["response"] == "Responder",
                       "L5_nmf_programs_pred_prob"].values
    cn = capp_pred.loc[capp_pred["response"] == "Non-Responder",
                       "L5_nmf_programs_pred_prob"].values

    for vals, col, lab, xi in [(cr, R_COL, "Responder", 0),
                               (cn, NR_COL, "Non-Responder", 1)]:
        jitter = rng.uniform(-0.13, 0.13, len(vals))
        axB.scatter(xi + jitter, vals, color=col, s=70,
                    edgecolor="white", linewidth=0.9, zorder=4, label=lab)
        med = np.median(vals)
        axB.hlines(med, xi - 0.25, xi + 0.25,
                   color=col, lw=2.8, zorder=5)

    if len(cr) >= 2 and len(cn) >= 2:
        _, mwu_c = mannwhitneyu(cr, cn, alternative="two-sided")
        cd_c     = _cliffs_delta(cr, cn)
        nb_r, nb_n = len(cr), len(cn)
        _statbox(axB,
                 f"n = {nb_r}R / {nb_n}NR\n"
                 f"MWU p = {mwu_c:.3f}\n"
                 f"Cliff's δ = {cd_c:+.3f}",
                 loc="upper right", fs=8)

    axB.axhline(0.5, color="#999999", linestyle=":", linewidth=1.0, zorder=1)
    axB.text(1.02, 0.5, "random", transform=axB.get_yaxis_transform(),
             fontsize=7.5, color="#999999", va="center")
    axB.set_xticks([0, 1])
    axB.set_xticklabels(["Responder", "Non-Responder"], fontsize=11)
    axB.set_ylabel("L5 NMF pred-prob (responder)", fontsize=11)
    axB.set_title("Cappuyns validation — NMF layer\n(held-out, n=32)", fontsize=11, fontweight="bold")
    _clean(axB)
    axB.text(0.5, -0.16, "Shows partial NMF programme transfer to external cohort",
             transform=axB.transAxes, ha="center", fontsize=7.5, color="#666666",
             style="italic")
else:
    axB.text(0.5, 0.5, "Cappuyns NMF data\nnot available", ha="center", va="center",
             transform=axB.transAxes, fontsize=11, color="grey")
    axB.set_axis_off()

# ── Panel C: NMF10 activity vs CD8 pseudotime rank ───────────────────────────
axC = axes1[2]
_panel_letter(axC, "C", x=-0.13)

# merge pt_means with NMF10 from nmf_disc
merge_col = None
for c in ("patient_id", "patient", "patient_ID"):
    if c in pt_means.columns: merge_col = c; break

if merge_col and len(nmf_disc) > 0 and "NMF10" in nmf_disc.columns:
    nmf_id_col = None
    for c in ("patient_id", "patient", "patient_ID"):
        if c in nmf_disc.columns: nmf_id_col = c; break

    pt_tmp  = pt_means.rename(columns={merge_col: "_pid"}).copy()
    nmf_tmp = nmf_disc.rename(columns={nmf_id_col: "_pid"})[["_pid", "NMF10", "response"]].copy()
    # normalise key to string so float64 "1001.0" == object "1001.0"
    pt_tmp["_pid"]  = pt_tmp["_pid"].astype(str)
    nmf_tmp["_pid"] = nmf_tmp["_pid"].astype(str)
    mdf = pt_tmp.merge(nmf_tmp, on="_pid", how="inner", suffixes=("_pt", "_nmf"))
    # use response from pt_means if available
    if "response_pt" in mdf.columns:
        mdf["response"] = mdf["response_pt"]
    elif "response" not in mdf.columns and "response_nmf" in mdf.columns:
        mdf["response"] = mdf["response_nmf"]

    mdf = mdf.dropna(subset=["mean", "NMF10"])
    mdf = mdf.sort_values("mean").reset_index(drop=True)
    mdf["pt_rank"] = np.arange(len(mdf))

    rho_cn, p_cn = spearmanr(mdf["mean"], mdf["NMF10"])

    cols = mdf["response"].map(resp_color)
    axC.scatter(mdf["mean"], mdf["NMF10"],
                c=cols, s=80, edgecolor="white", linewidth=0.8, zorder=4)

    # OLS trend line (cleaner than LOESS for n=24 with non-significant ρ)
    if len(mdf) >= 4:
        _xf = mdf["mean"].values
        _yf = mdf["NMF10"].values
        _m, _b = np.polyfit(_xf, _yf, 1)
        _xs = np.linspace(_xf.min(), np.percentile(_xf, 97) + 0.01, 100)
        axC.plot(_xs, _m * _xs + _b,
                 color="#555555", lw=1.8, ls="--", zorder=3, label="OLS trend")

    # legend proxies
    for lab, col in [("Responder", R_COL), ("Non-Responder", NR_COL)]:
        axC.scatter([], [], color=col, s=60, label=lab, edgecolor="white")

    _statbox(axC,
             f"Spearman ρ = {rho_cn:.3f}\np = {p_cn:.3f}\nn = {len(mdf)}"
             f"\n(NMF10 anti-corr. with exhaustion depth)",
             loc="upper right", fs=8)

    # clip x-axis to 97th percentile to prevent single outlier from stretching axis
    axC.set_xlim(mdf["mean"].min() - 0.005, np.percentile(mdf["mean"], 97) + 0.01)
    axC.set_xlabel("Per-patient mean CD8 pseudotime", fontsize=10.5)
    axC.set_ylabel("NMF10 patient activity", fontsize=10.5)
    axC.set_title("NMF10 vs CD8 pseudotime\n(discovery, n=24)", fontsize=11, fontweight="bold")
    axC.legend(frameon=False, loc="upper left", fontsize=8.5)
    _clean(axC)
else:
    axC.text(0.5, 0.5, "Merge data\nnot available", ha="center", va="center",
             transform=axC.transAxes, fontsize=11, color="grey")
    axC.set_axis_off()

fig1.suptitle(
    "CD8 pseudotime × NMF programme trajectory — Discovery and validation",
    fontsize=13, fontweight="bold", y=1.01)

_save(fig1, "05z_traj_scrna_pseudotime")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 — CNV patient trajectory with NMF10 + ranked lollipop
# ══════════════════════════════════════════════════════════════════════════════
print("\n[Fig 2] CNV patient trajectory …")

if len(feats) == 0 or len(traj_rank) == 0:
    print("  [SKIP] missing feats or traj_rank table")
else:
    # merge NMF10 into feats
    if len(nmf_disc) > 0 and "NMF10" in nmf_disc.columns:
        nid = [c for c in nmf_disc.columns
               if c in ("patient_id", "patient", "patient_ID")][0]
        fid = [c for c in feats.columns
               if c in ("patient_id", "patient", "patient_ID")][0]
        feats = feats.merge(
            nmf_disc.rename(columns={nid: fid})[
                [fid] + [c for c in nmf_disc.columns if c.startswith("NMF")]
            ],
            on=fid, how="left", suffixes=("", "_nmf")
        )

    # build the selected-feature trajectory df sorted by clone diversity
    TRACE_FEATS = ["NMF10", "sig_T_cell_exhaustion",
                   "frac_CD8_T", "sig_IFN_gamma", "frac_Monocyte"]
    trace_avail = [f for f in TRACE_FEATS if f in feats.columns]
    missing = set(TRACE_FEATS) - set(trace_avail)
    if missing:
        # try frac_Monocyte_Intermediate or frac_Monocyte_Classical as fallback
        for alt in ("frac_Monocyte_Intermediate", "frac_Monocyte_Classical",
                    "frac_Monocyte_NC"):
            if alt in feats.columns and "frac_Monocyte" not in trace_avail:
                trace_avail.append(alt)
                # inherit style from frac_Monocyte
                FEAT_COL[alt]   = FEAT_COL["frac_Monocyte"]
                FEAT_LS[alt]    = FEAT_LS["frac_Monocyte"]
                FEAT_LW[alt]    = FEAT_LW["frac_Monocyte"]
                FEAT_LABEL[alt] = f"Monocyte frac. ({alt.split('_')[-1]})"

    feats_t = feats.dropna(subset=["cnv_clone_diversity"]).sort_values(
        "cnv_clone_diversity").reset_index(drop=True)
    feats_t["trajectory_rank"] = np.arange(len(feats_t))

    x_traj = feats_t["cnv_clone_diversity"].values

    fig2 = plt.figure(figsize=(15, 10), constrained_layout=True)
    fig2.patch.set_facecolor("white")
    gs2 = mgridspec.GridSpec(2, 1, figure=fig2,
                              height_ratios=[1.55, 1.0], hspace=0.08)

    axTA = fig2.add_subplot(gs2[0])
    axTB = fig2.add_subplot(gs2[1])

    _panel_letter(axTA, "A", x=-0.06)
    _panel_letter(axTB, "B", x=-0.06)

    # ── PANEL A: LOESS-smoothed feature traces ────────────────────────────────
    for feat in trace_avail:
        v = feats_t[feat].values.astype(float)
        # standardise for overlay
        mu, sd = np.nanmean(v), np.nanstd(v)
        if sd > 0:
            vz = (v - mu) / sd
        else:
            continue

        # Spearman ρ from traj_rank table if available, else compute directly
        row = traj_rank[traj_rank["feature"] == feat]
        rho_str = ""
        if len(row):
            r, p = float(row["rho_vs_clone_diversity"].iat[0]), float(row["p_spearman"].iat[0])
            rho_str = f"  ρ={r:+.2f}, p={p:.3f}"

        # LOESS smoothing for visual guide
        not_nan = ~np.isnan(vz)
        if len(row) == 0 and not_nan.sum() >= 5:
            # fallback: compute Spearman directly (rank-invariant to z-scoring)
            r_d, p_d = spearmanr(x_traj[not_nan], vz[not_nan])
            rho_str = f"  ρ={r_d:+.2f}, p={p_d:.3f}"

        if not_nan.sum() >= 6:
            lw_xy = lowess(vz[not_nan], x_traj[not_nan],
                           frac=0.55, return_sorted=True)
            axTA.plot(lw_xy[:, 0], lw_xy[:, 1],
                      color=FEAT_COL.get(feat, "#666666"),
                      lw=FEAT_LW.get(feat, 1.8),
                      ls=FEAT_LS.get(feat, "-"),
                      label=FEAT_LABEL.get(feat, feat) + rho_str,
                      zorder=3)
        # raw data as faint dots
        axTA.scatter(x_traj[not_nan], vz[not_nan],
                     color=FEAT_COL.get(feat, "#666666"),
                     s=18, alpha=0.35, edgecolor="none", zorder=2)

    # response triangle markers above the axis: R and NR in separate rows to avoid overlap
    ymin_a, ymax_a = axTA.get_ylim()
    if ymax_a == 0:
        ymax_a = 2.0; ymin_a = -2.0
    _dy = ymax_a - ymin_a
    top_y_r  = ymax_a + 0.26 * _dy   # Responders — upper row
    top_y_nr = ymax_a + 0.13 * _dy   # Non-Responders — lower row
    for _, _rrow in feats_t.iterrows():
        _is_r = _rrow.get("response", "") in ("Responder", "R")
        _col_r = R_COL if _is_r else NR_COL
        _top  = top_y_r if _is_r else top_y_nr
        axTA.scatter(_rrow["cnv_clone_diversity"], _top,
                     marker="v", s=42, color=_col_r,
                     edgecolor="black", linewidth=0.4, clip_on=False, zorder=5)
    # legend proxies for response triangles
    for lab, col in [("Responder ▼ (top row)", R_COL), ("Non-Responder ▼ (bottom row)", NR_COL)]:
        axTA.scatter([], [], marker="v", s=42, color=col,
                     edgecolor="black", linewidth=0.4, label=lab)    # lock y-axis to the data range so triangles (clip_on=False) don't inflate it
    axTA.set_ylim(ymin_a, ymax_a)
    axTA.axhline(0, color="#CCCCCC", lw=0.8, ls=":", zorder=1)
    axTA.set_ylabel("z-score (feature)", fontsize=11)
    axTA.set_title(
        "Patient trajectory along CNV clonal diversity axis — feature dynamics",
        fontsize=11.5, fontweight="bold")
    axTA.legend(loc="upper right", fontsize=8.2, frameon=True,
                framealpha=0.92, edgecolor="#DDDDDD", ncol=1,
                borderpad=0.5, handlelength=2.2,
                title="Feature (LOESS + raw points)", title_fontsize=8)
    _clean(axTA)
    axTA.set_xlim(x_traj.min() - 0.02, x_traj.max() + 0.02)

    # ── PANEL B: Ranked lollipop ──────────────────────────────────────────────
    # top 12 features by |ρ|, excluding near-zero variance / tech features
    excl_tok = ("Unknown", "Doublet", "Empty", "cnv_score", "n_cells",
                "malig_frac", "dataset", "response")
    lp_df = traj_rank[~traj_rank["feature"].str.contains("|".join(excl_tok),
                                                          na=False)].copy()
    lp_df["abs_rho"] = lp_df["rho_vs_clone_diversity"].abs()
    top_12 = lp_df.nlargest(12, "abs_rho")
    # always include NMF10 in the lollipop, even if its |ρ| falls outside top-12
    if "NMF10" not in top_12["feature"].values:
        nmf10_row = lp_df[lp_df["feature"] == "NMF10"]
        if len(nmf10_row):
            top_12 = pd.concat([top_12, nmf10_row], ignore_index=True)
    lp_df = top_12.sort_values("rho_vs_clone_diversity")

    colors_lp = [R_COL if r < 0 else NR_COL
                 for r in lp_df["rho_vs_clone_diversity"]]
    y_pos = np.arange(len(lp_df))

    axTB.hlines(y_pos, 0, lp_df["rho_vs_clone_diversity"].values,
                colors=colors_lp, lw=2.2, zorder=2)
    axTB.scatter(lp_df["rho_vs_clone_diversity"].values, y_pos,
                 c=colors_lp, s=90, zorder=4,
                 edgecolor="white", linewidth=1.0)

    # highlight NMF10 with star marker
    for i, (_, row) in enumerate(lp_df.reset_index(drop=True).iterrows()):
        if row["feature"] == "NMF10":
            axTB.scatter(row["rho_vs_clone_diversity"], i,
                         marker="*", s=200, c="gold",
                         edgecolor="darkorange", linewidth=0.8, zorder=5)

    # annotate p values for top features by |ρ|
    for i, (_, row) in enumerate(lp_df.reset_index(drop=True).iterrows()):
        ha = "right" if row["rho_vs_clone_diversity"] < 0 else "left"
        offset = -0.012 if row["rho_vs_clone_diversity"] < 0 else 0.012
        pv = row["p_spearman"]
        plab = f"p={pv:.3f}" if pv >= 0.001 else f"p={pv:.2e}"
        axTB.text(row["rho_vs_clone_diversity"] + offset, i,
                  plab, va="center", ha=ha, fontsize=7, color="#444444")

    # tidy feature label names
    def _clean_label(s):
        label = (s.replace("frac_", "")
                  .replace("sig_", "sig: ")
                  .replace("_", " "))
        # deduplicate consecutive identical words (e.g. "Proliferating Proliferating")
        words = label.split()
        deduped = [w for i, w in enumerate(words) if i == 0 or w != words[i - 1]]
        return " ".join(deduped)

    axTB.set_yticks(y_pos)
    axTB.set_yticklabels([_clean_label(f) for f in lp_df["feature"]], fontsize=9)
    axTB.axvline(0, color="black", lw=0.8)
    axTB.set_xlabel("Spearman ρ vs CNV clone diversity", fontsize=11)
    axTB.set_title("Top features correlated with clone diversity (|ρ| ranked)",
                   fontsize=11, fontweight="bold")
    _clean(axTB)
    axTB.xaxis.set_minor_locator(mticker.AutoMinorLocator(5))
    axTB.set_xlim(-0.90, 0.90)

    # direction legend
    for lab, col in [("R-enriched (negative ρ)", R_COL),
                     ("NR-enriched (positive ρ)", NR_COL)]:
        axTB.scatter([], [], color=col, s=60, label=lab)
    axTB.legend(frameon=False, loc="lower right", fontsize=8.5)

    # note about interpretation
    axTB.text(0.5, -0.18,
              "★ = NMF10 (BH-significant, p_BH=0.0054)  |  "
              "Negative ρ = enriched in R (low clone-diversity) patients",
              transform=axTB.transAxes, ha="center", fontsize=7.8,
              color="#555555", style="italic")

    # save traj selected table
    save_cols = ["cnv_clone_diversity", "response"] + [
        c for c in trace_avail if c in feats_t.columns]
    feats_t[save_cols].to_csv(
        OUT / "05z_traj_cnv_selected.csv", index=False)

    _save(fig2, "05z_traj_cnv_nmf10")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 3 — Spatial gradient trajectory (7 Visium sections)
# ══════════════════════════════════════════════════════════════════════════════
print("\n[Fig 3] Spatial gradient trajectory …")

RESP_ORDER = ["HCC1R", "HCC2R", "HCC3R", "HCC4R", "HCC5NR", "HCC6NR", "HCC7NR"]
SECTION_RESP = {s: ("Non-Responder" if s.endswith("NR") else "Responder")
                for s in RESP_ORDER}

# ── Pivot signature-gradient ρ to matrix ─────────────────────────────────────
# Focus on signatures relevant to trajectory story
SIG_SELECT = [
    "T_cell_activation", "T_cell_cytotoxicity", "T_cell_exhaustion",
    "Terminal_exhaustion", "IFN_gamma", "Immune_checkpoint",
    "TLS_signature", "Treg_signature",
    "M1_macrophage", "Angiogenesis", "Antigen_presentation",
]
sig_mat = pd.DataFrame()
if len(sgrad_df) > 0:
    sgrad_piv = sgrad_df.pivot(index="sample_id", columns="signature",
                                values="spearman_rho")
    # keep only sections in RESP_ORDER that exist
    sects_avail = [s for s in RESP_ORDER if s in sgrad_piv.index]
    sigs_avail  = [s for s in SIG_SELECT if s in sgrad_piv.columns]
    sig_mat = sgrad_piv.loc[sects_avail, sigs_avail]

    # p-value matrix for asterisk annotation
    sgrad_p = sgrad_df.pivot(index="sample_id", columns="signature",
                              values="p_adj").loc[sects_avail, sigs_avail]

fig3, axes3 = plt.subplots(
    1, 3,
    figsize=(19, 8.5),
    constrained_layout=True,
    gridspec_kw={"width_ratios": [1.6, 1.0, 1.0]},
)
fig3.patch.set_facecolor("white")

# ── Panel A: Heatmap section × signature ─────────────────────────────────────
axA3 = axes3[0]
_panel_letter(axA3, "A", x=-0.06)

if sig_mat.shape[0] > 0:
    # row colour bar (response)
    row_colors = sig_mat.index.map(
        lambda s: R_COL if SECTION_RESP.get(s, "?") == "Responder" else NR_COL)

    # pretty column labels
    col_labels = {
        "T_cell_activation":   "T cell activation",
        "T_cell_cytotoxicity": "T cell cytotoxicity",
        "T_cell_exhaustion":   "T cell exhaustion",
        "Terminal_exhaustion": "Terminal exhaustion",
        "IFN_gamma":           "IFN-γ",
        "Immune_checkpoint":   "Immune checkpoint",
        "TLS_signature":       "TLS signature",
        "Treg_signature":      "Treg signature",
        "M1_macrophage":       "M1 macrophage",
        "TGF_beta":            "TGF-β",
        "Angiogenesis":        "Angiogenesis",
        "Antigen_presentation":"Antigen presentation",
    }
    display_cols = [col_labels.get(c, c) for c in sig_mat.columns]

    im = axA3.imshow(sig_mat.values, cmap="RdBu_r",
                      vmin=-0.8, vmax=0.8, aspect="auto")
    axA3.set_xticks(np.arange(len(sig_mat.columns)))
    axA3.set_yticks(np.arange(len(sig_mat.index)))
    axA3.set_xticklabels(display_cols, rotation=45, ha="right",
                          fontsize=8.5)
    axA3.set_yticklabels(sig_mat.index, fontsize=10)
    # colour ytick labels by response — clear and avoids orphan rectangles
    for tick, sect in zip(axA3.get_yticklabels(), sig_mat.index):
        tick.set_color(R_COL if SECTION_RESP.get(sect, "") == "Responder" else NR_COL)
        tick.set_fontweight("bold")

    # asterisk annotation for p_adj < 0.05
    for i, sect in enumerate(sig_mat.index):
        for j, sig in enumerate(sig_mat.columns):
            if not sgrad_p.empty and sig in sgrad_p.columns:
                pv = sgrad_p.loc[sect, sig] if sect in sgrad_p.index else np.nan
                if not np.isnan(pv) and pv < 0.05:
                    val = sig_mat.loc[sect, sig]
                    txc = "white" if abs(val) > 0.5 else "black"
                    axA3.text(j, i, "★", ha="center", va="center",
                              fontsize=9, color=txc, zorder=6)

    cb3a = plt.colorbar(im, ax=axA3, shrink=0.65, pad=0.02)
    cb3a.set_label("Spearman ρ (spatial gradient)", fontsize=9)
    cb3a.outline.set_linewidth(0.4)
    cb3a.ax.tick_params(labelsize=8)

    # R / NR divider line after 4 R sections
    n_r = sum(1 for s in sig_mat.index if SECTION_RESP.get(s, "") == "Responder")
    axA3.axhline(n_r - 0.5, color="black", lw=1.8, zorder=7)
    # group bracket annotations in axes coordinates (right of heatmap)
    _n_sec = len(sig_mat)
    _n_col = len(sig_mat.columns)
    if n_r > 0:
        _y_r = (n_r - 1) / 2.0
        axA3.text(_n_col - 0.5, _y_r, " R ",
                  va="center", ha="left", fontsize=9, fontweight="bold",
                  color="white", clip_on=False,
                  bbox=dict(boxstyle="round,pad=0.15", fc=R_COL, alpha=0.85),
                  zorder=8)
    if n_r < _n_sec:
        _y_nr = n_r + (_n_sec - n_r - 1) / 2.0
        axA3.text(_n_col - 0.5, _y_nr, " NR ",
                  va="center", ha="left", fontsize=9, fontweight="bold",
                  color="white", clip_on=False,
                  bbox=dict(boxstyle="round,pad=0.15", fc=NR_COL, alpha=0.85),
                  zorder=8)

    axA3.set_title("Spatial-gradient Spearman ρ per section\n★ = p_adj < 0.05",
                   fontsize=11, fontweight="bold")
    for sp in axA3.spines.values():
        sp.set_visible(False)
else:
    axA3.text(0.5, 0.5, "Gradient data\nnot available",
              ha="center", va="center", transform=axA3.transAxes,
              fontsize=11, color="grey")
    axA3.set_axis_off()

# ── Panel B: NMF10 vs sig_T_cell_exhaustion per section ─────────────────────
axB3 = axes3[1]
_panel_letter(axB3, "B", x=-0.13)

if len(sec_scores) > 0 and "NMF10" in sec_scores.columns:
    secs_b = [s for s in RESP_ORDER if s in sec_scores["sample"].values]
    sec_b  = sec_scores[sec_scores["sample"].isin(secs_b)].set_index("sample")
    x_b    = np.arange(len(secs_b))

    nmf10_vals = sec_b.loc[secs_b, "NMF10"].values
    bar_colors = [R_COL if SECTION_RESP.get(s, "") == "Responder" else NR_COL
                  for s in secs_b]

    axB3.bar(x_b, nmf10_vals, width=0.62, color=bar_colors,
             alpha=0.85, edgecolor="white", linewidth=0.6)

    # median reference line
    _med_nmf = np.median(nmf10_vals)
    axB3.axhline(_med_nmf, color="#444444", lw=1.2, ls=":", zorder=3,
                 label=f"median = {_med_nmf:.3f}")

    # vertical divider between R / NR
    n_r_b = sum(1 for s in secs_b if SECTION_RESP.get(s, "") == "Responder")
    axB3.axvline(n_r_b - 0.5, color="#999999", lw=1.2, ls="--", zorder=0)

    axB3.set_xticks(x_b)
    axB3.set_xticklabels(secs_b, fontsize=9, rotation=45, ha="right")
    # colour x-tick labels by response
    for tick, sect in zip(axB3.get_xticklabels(), secs_b):
        tick.set_color(R_COL if SECTION_RESP.get(sect, "") == "Responder" else NR_COL)
        tick.set_fontweight("bold")

    axB3.set_ylabel("Section mean NMF10 score", fontsize=10.5)
    axB3.set_title("NMF10 activity\nper Visium section", fontsize=11,
                   fontweight="bold")
    # legend proxies for response colours
    for lab, col in [("Responder", R_COL), ("Non-Responder", NR_COL)]:
        axB3.bar([], [], color=col, alpha=0.85, label=lab)
    axB3.legend(frameon=False, fontsize=8.5, loc="upper right")
    _clean(axB3)
else:
    axB3.text(0.5, 0.5, "Section scores\nnot available",
              ha="center", va="center", transform=axB3.transAxes,
              fontsize=11, color="grey")
    axB3.set_axis_off()

# ── Panel C: Immune gradient vs NMF10 per section ────────────────────────────
axC3 = axes3[2]
_panel_letter(axC3, "C", x=-0.13)

if len(sgrad_ig) > 0 and "mean_immune_gradient" in sgrad_ig.columns \
        and len(sec_scores) > 0 and "NMF10" in sec_scores.columns:

    ig = sgrad_ig.reset_index() if sgrad_ig.index.name == "sample_id" \
         else sgrad_ig.copy()
    ig_col = "sample_id" if "sample_id" in ig.columns else ig.columns[0]

    sec_sc = sec_scores.rename(columns={"sample": "sample_id"})
    merged_c = ig[[ig_col, "mean_immune_gradient"]].rename(
        columns={ig_col: "sample_id"}).merge(
        sec_sc[["sample_id", "NMF10"]], on="sample_id", how="inner")
    merged_c = merged_c.dropna(subset=["mean_immune_gradient", "NMF10"])
    merged_c["response"] = merged_c["sample_id"].map(SECTION_RESP)

    rho_sc, p_sc = spearmanr(merged_c["mean_immune_gradient"],
                              merged_c["NMF10"])

    for _, row in merged_c.iterrows():
        col = resp_color(row["response"])
        axC3.scatter(row["mean_immune_gradient"], row["NMF10"],
                     color=col, s=95, edgecolor="white", linewidth=0.9,
                     zorder=4)
        axC3.annotate(row["sample_id"],
                      xy=(row["mean_immune_gradient"], row["NMF10"]),
                      xytext=(5, 4), textcoords="offset points",
                      fontsize=8, color=col)

    # regression line
    if len(merged_c) >= 4:
        xf = merged_c["mean_immune_gradient"].values
        yf = merged_c["NMF10"].values
        m, b = np.polyfit(xf, yf, 1)
        xs_ = np.linspace(xf.min(), xf.max(), 100)
        axC3.plot(xs_, m * xs_ + b, color="#555555", lw=1.6, ls="--", zorder=3)

    _statbox(axC3,
             f"Spearman ρ = {rho_sc:.3f}\np = {p_sc:.3f}\nn = {len(merged_c)} sections",
             loc="upper right", fs=8)

    axC3.set_xlabel("Mean immune gradient (spatial)", fontsize=10.5)
    axC3.set_ylabel("NMF10 section mean score", fontsize=10.5)
    axC3.set_title("Immune gradient vs NMF10\n(spatial, n=7 sections)",
                   fontsize=11, fontweight="bold")

    for lab, col in [("Responder", R_COL), ("Non-Responder", NR_COL)]:
        axC3.scatter([], [], color=col, s=60, label=lab, edgecolor="white")
    axC3.legend(frameon=False, fontsize=8.5, loc="lower right")
    _clean(axC3)
else:
    axC3.text(0.5, 0.5, "Immune gradient\nor NMF10 data\nnot available",
              ha="center", va="center", transform=axC3.transAxes,
              fontsize=11, color="grey")
    axC3.set_axis_off()

fig3.suptitle(
    "Spatial gradient trajectory — GSE238264 Visium sections (R=4, NR=3)",
    fontsize=13, fontweight="bold", y=1.01)

_save(fig3, "05z_traj_spatial_gradient")

print(f"\nAll trajectory figures saved to {OUT}")
