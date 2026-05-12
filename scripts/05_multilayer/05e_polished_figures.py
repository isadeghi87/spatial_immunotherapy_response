#!/usr/bin/env python3
"""
05e_polished_figures.py — Regenerate the key multi-layer figures with
stronger, publication-ready aesthetics. Replaces raw 05a/05b/05d panels.

Philosophy
----------
 - Show distribution shape, not just boxes (half-violin + strip + box).
 - Use a semantically meaningful layer palette (Tol muted / Okabe-Ito-ish).
 - Order heatmap rows and columns by the analytical message, not alphabet.
 - Separate Responders and Non-Responders with an annotation bar.
 - Remove chart-junk (inner gridlines, unused axes).
 - Annotate effect sizes and p-values on every comparison figure.

Inputs
------
  out/05_multilayer/tables/patient_multilayer_features.csv
  out/05_multilayer/tables/cnv_response_stats.csv
  out/05_multilayer/tables/nmf_programs_patient_activity.csv
  out/05_multilayer/tables/nmf_programs_response_stats.csv
  out/05_multilayer/tables/nmf_programs_spatial_activity.csv
  out/05_multilayer/tables/nmf_programs_gene_loadings.csv
  out/05_multilayer/tables/nmf_programs_bulk_activity.csv
  out/05_multilayer/tables/balanced_layer_aucs.csv
  out/05_multilayer/tables/balanced_selected_features.csv
  out/05_multilayer/tables/balanced_variance_partition.csv

Outputs  (all PNG + PDF, 300 dpi, tight bbox)
-------
  05e_cnv_RvsNR.png                 ← replaces 05a_cnv_RvsNR
  05e_nmf_program_vs_response.png   ← replaces 05b_nmf_program_vs_response
  05e_nmf_spatial_heatmap.png       ← replaces 05b_nmf_spatial_projection
  05e_nmf_bulk_vs_sc_concordance.png
  05e_nmf_top_genes_dot.png
  05e_roc_by_layer.png              ← replaces 05d_balanced_roc
  05e_variance_partition.png        ← replaces 05d_balanced_variance_partition
  05e_feature_effects.png           ← replaces 05d_balanced_feature_effects
  05e_auc_bars.png
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
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from scipy.stats import mannwhitneyu, ttest_ind, spearmanr

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score, roc_curve

from config.publication_style import apply_style, PAL, label_clean
apply_style()

OUT = PROJECT / "out" / "05_multilayer"
TAB = OUT / "tables"; FIG = OUT / "figures"

# ── Semantic palettes ──────────────────────────────────────────────────
RESP_PAL = {"Responder": PAL["Responder"], "Non-Responder": PAL["Non-Responder"]}

# Tol muted-ish layer palette — distinct & print-friendly
LAYER_COL = {
    "L1_tumour_intrinsic": "#882255",
    "L2_regulon":          "#DDCC77",
    "L3_composition":      "#117733",
    "L4_io_signatures":    "#332288",
    "L5_nmf_programs":     "#AA4499",
    "Balanced ensemble":   "#000000",
}
LAYER_LABEL = {
    "L1_tumour_intrinsic": "L1  Tumour-intrinsic",
    "L2_regulon":          "L2  Regulon (IRF1/MYC)",
    "L3_composition":      "L3  Cell composition",
    "L4_io_signatures":    "L4  IO signatures",
    "L5_nmf_programs":     "L5  NMF programs",
    "Balanced ensemble":   "Balanced ensemble",
}

# Perceptually smooth sequential
VIR = plt.get_cmap("magma")
# Diverging (white center) — Tol sunset inspired
DIV = LinearSegmentedColormap.from_list(
    "tol_div", ["#364B9A", "#6EA6CD", "#FAFBCD", "#F67E4B", "#A50026"])

FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

def _save(fig, stem):
    for fmt in ("png", "pdf"):
        fig.savefig(FIG / f"{stem}.{fmt}",
                     dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"   saved  {stem}")

# ── Load everything ────────────────────────────────────────────────────
feats = pd.read_csv(TAB / "patient_multilayer_features.csv")
if "dataset" not in feats.columns:
    feats["dataset"] = np.where(feats["patient_id"].astype(str).str.startswith("pat"),
                                 "Cappuyns2023", "GSE206325")
feats_disc = feats.loc[feats["dataset"].astype(str) != "Cappuyns2023"].copy()
cnv_stats = pd.read_csv(TAB / "cnv_response_stats.csv")
nmf_patient = pd.read_csv(TAB / "nmf_programs_patient_activity.csv")
nmf_resp = pd.read_csv(TAB / "nmf_programs_response_stats.csv")
nmf_spatial = pd.read_csv(TAB / "nmf_programs_spatial_activity.csv") \
    if (TAB / "nmf_programs_spatial_activity.csv").exists() else None
nmf_bulk = pd.read_csv(TAB / "nmf_programs_bulk_activity.csv", index_col=0) \
    if (TAB / "nmf_programs_bulk_activity.csv").exists() else None
nmf_loadings = pd.read_csv(TAB / "nmf_programs_gene_loadings.csv", index_col=0)
aucs = pd.read_csv(TAB / "balanced_layer_aucs.csv")
sel_feats = pd.read_csv(TAB / "balanced_selected_features.csv")
vp = pd.read_csv(TAB / "balanced_variance_partition.csv")

print("[05e] Polishing multi-layer figures")

# ══════════════════════════════════════════════════════════════════════
# 1. CNV R-vs-NR  —  violin + strip + box
# ══════════════════════════════════════════════════════════════════════
cnv_panels = [
    ("cnv_score_mean",     "CNV score  (mean)"),
    ("cnv_score_sd",       "CNV score  (SD, clonal heterogeneity)"),
    ("malig_frac",         "Malignant fraction"),
    ("cnv_clone_diversity","CNV clone diversity  (Shannon)"),
]
pvals = dict(zip(cnv_stats["metric"], cnv_stats["p_value"])) \
    if "metric" in cnv_stats.columns else {}

fig, axes = plt.subplots(1, 4, figsize=(14, 4.2), constrained_layout=True)
order = ["Responder", "Non-Responder"]
for ax, (col, ttl) in zip(axes, cnv_panels):
    d = feats_disc[["response", col]].dropna().copy()
    d["response"] = d["response"].astype(str)
    # half violin
    sns.violinplot(data=d, x="response", y=col, order=order, hue="response",
                   palette=RESP_PAL, inner=None, cut=0, ax=ax,
                   linewidth=0.8, legend=False, saturation=0.9)
    # translucent
    for coll in ax.collections:
        coll.set_alpha(0.35)
    # narrow box inside
    sns.boxplot(data=d, x="response", y=col, order=order, ax=ax,
                width=0.18, showfliers=False, boxprops=dict(facecolor="white",
                edgecolor="black", linewidth=1.0),
                medianprops=dict(color="black", linewidth=1.6),
                whiskerprops=dict(color="black", linewidth=0.9),
                capprops=dict(color="black", linewidth=0.9))
    # strip
    sns.stripplot(data=d, x="response", y=col, order=order, hue="response",
                  palette={k: v for k, v in RESP_PAL.items()},
                  jitter=0.12, size=4.2, edgecolor="black", linewidth=0.4,
                  alpha=0.95, ax=ax, legend=False)
    # stats
    r = d.loc[d.response == "Responder", col].dropna().values
    nr = d.loc[d.response == "Non-Responder", col].dropna().values
    if len(r) >= 2 and len(nr) >= 2:
        U, pv = mannwhitneyu(r, nr, alternative="two-sided")
        eff = (r.mean() - nr.mean()) / (np.concatenate([r, nr]).std() + 1e-12)
    else:
        pv, eff = np.nan, np.nan
    # annotation
    star = "***" if pv < 0.001 else "**" if pv < 0.01 \
        else "*" if pv < 0.05 else "ns"
    y_max = np.nanmax(d[col].values)
    y_min = np.nanmin(d[col].values)
    span = y_max - y_min
    h = y_max + span * 0.08
    ax.plot([0, 0, 1, 1], [h, h + span * 0.03, h + span * 0.03, h],
            color="black", linewidth=0.9)
    ax.text(0.5, h + span * 0.05,
            f"{star}   p={pv:.2g}   d={eff:+.2f}",
            ha="center", va="bottom", fontsize=8.5)
    ax.set_ylim(y_min - span * 0.08, h + span * 0.15)
    ax.set_title(ttl, fontsize=10.5, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    # bold significant panel
    if pv < 0.05:
        for side in ("top", "right", "bottom", "left"):
            ax.spines[side].set_linewidth(1.8)
            ax.spines[side].set_color("#B22222")
    ax.tick_params(axis="x", labelsize=9)

# annotations for axes
axes[0].set_ylabel("per-patient value", fontsize=9)
fig.suptitle("Tumour-intrinsic genomic state  —  discovery cohort only",
             fontsize=12, fontweight="bold")
_save(fig, "05e_cnv_RvsNR")


# ══════════════════════════════════════════════════════════════════════
# 2. NMF programs  vs  response — effect-size forest with top genes
# ══════════════════════════════════════════════════════════════════════
ps = nmf_resp.copy()
ps = ps.sort_values("effect_size")
fig, ax = plt.subplots(figsize=(10, max(3.5, len(ps) * 0.38)),
                        constrained_layout=True)
cols = [PAL["Responder"] if e > 0 else PAL["Non-Responder"]
        for e in ps["effect_size"]]
ax.barh(ps["program"], ps["effect_size"], color=cols, edgecolor="black",
        linewidth=0.6, height=0.7)
ax.axvline(0, color="black", linewidth=0.9)
# significance halos
for y, (_, r) in enumerate(ps.iterrows()):
    if r["p_value"] < 0.05:
        ax.barh(r["program"], r["effect_size"],
                color="none", edgecolor="#FFA500", linewidth=2.2, height=0.7)
# annotate: top genes + p
xmax = ps["effect_size"].abs().max() * 1.05
# anchor the gene-annotation column to a fixed x (right of data range)
x_anchor = xmax * 1.25
for y, (_, r) in enumerate(ps.iterrows()):
    ef = r["effect_size"]
    stars = "***" if r["p_value"] < 0.001 else "**" if r["p_value"] < 0.01 \
        else "*" if r["p_value"] < 0.05 else ""
    # p-value label just past the bar tip (always same side as bar)
    ax.text(ef + 0.03 * np.sign(ef or 1), y,
             f"p={r['p_value']:.2g} {stars}",
             va="center", ha="left" if ef > 0 else "right", fontsize=8)
    # top genes in a fixed right column, italic
    txt = ", ".join(t.strip() for t in r["top_genes"].split(",")[:5])
    ax.text(x_anchor, y, txt, va="center", ha="left",
             fontsize=7.5, color="#444444", style="italic")

ax.set_xlabel("Standardised mean difference  (Responder − Non-Responder)",
              fontsize=9.5)
ax.set_title("NMF programs vs response  (k = 10; marked * when p < 0.05)",
              fontsize=11, fontweight="bold")
# legend anchored outside data to avoid overlap
legend_elem = [
    Patch(facecolor=PAL["Responder"], label="Higher in Responders"),
    Patch(facecolor=PAL["Non-Responder"], label="Higher in Non-Responders"),
    Patch(facecolor="none", edgecolor="#FFA500",
           linewidth=2, label="p < 0.05  (Mann-Whitney)"),
]
ax.legend(handles=legend_elem, loc="upper left",
           bbox_to_anchor=(0.02, -0.08),
           ncol=3, frameon=False, fontsize=8)
# fix x-limits to show bars + gene annotations on one line each
ax.set_xlim(-xmax * 1.1, x_anchor + xmax * 1.4)
ax.axvline(x_anchor - xmax * 0.03, color="#DDDDDD", linewidth=0.6)
_save(fig, "05e_nmf_program_vs_response")


# ══════════════════════════════════════════════════════════════════════
# 3. NMF Visium spatial projection  — ordered, response-annotated
# ══════════════════════════════════════════════════════════════════════
if nmf_spatial is not None and len(nmf_spatial):
    RM = {"HCC1R": "R", "HCC2R": "R", "HCC3R": "R", "HCC4R": "R",
           "HCC5NR": "NR", "HCC6NR": "NR", "HCC7NR": "NR"}
    sp = nmf_spatial.copy()
    sp["resp"] = sp["section"].map(RM)
    mat = sp.pivot_table(index="program", columns="section",
                           values="mean", observed=True)
    # column order: R first then NR
    col_order = [c for c, r in RM.items() if r == "R" and c in mat.columns] + \
                [c for c, r in RM.items() if r == "NR" and c in mat.columns]
    mat = mat[col_order]
    # row order: by R - NR effect size
    eff = {}
    for p in mat.index:
        r = sp[(sp.program == p) & (sp.resp == "R")]["mean"].values
        nr = sp[(sp.program == p) & (sp.resp == "NR")]["mean"].values
        eff[p] = (r.mean() if len(r) else 0) - (nr.mean() if len(nr) else 0)
    mat = mat.loc[sorted(mat.index, key=lambda p: eff[p], reverse=True)]

    # z-score rows to make patterns visible
    matz = (mat.sub(mat.mean(axis=1), axis=0)
               .div(mat.std(axis=1) + 1e-9, axis=0))

    fig = plt.figure(figsize=(8.5, 6.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.07, 1.0],
                           width_ratios=[1.0, 0.04], hspace=0.04,
                           wspace=0.02)
    ax_ann = fig.add_subplot(gs[0, 0])
    ax_hm  = fig.add_subplot(gs[1, 0], sharex=ax_ann)
    ax_cb  = fig.add_subplot(gs[1, 1])
    # annotation bar
    ann_colors = [PAL["Responder"] if RM[c] == "R" else PAL["Non-Responder"]
                   for c in col_order]
    ax_ann.bar(range(len(col_order)), [1] * len(col_order),
                color=ann_colors, edgecolor="white", linewidth=1.2,
                width=1.0)
    ax_ann.set_xlim(-0.5, len(col_order) - 0.5)
    ax_ann.set_ylim(0, 1)
    ax_ann.set_axis_off()
    # heatmap
    im = ax_hm.imshow(matz.values, aspect="auto", cmap=DIV,
                       vmin=-2, vmax=2)
    ax_hm.set_yticks(range(matz.shape[0])); ax_hm.set_yticklabels(matz.index)
    ax_hm.set_xticks(range(matz.shape[1])); ax_hm.set_xticklabels(col_order,
                                                                     rotation=30, ha="right")
    # colour labels by response
    for i, lbl in enumerate(ax_hm.get_xticklabels()):
        lbl.set_color(PAL["Responder"] if RM[col_order[i]] == "R"
                       else PAL["Non-Responder"])
        lbl.set_fontweight("bold")
    # annotate top-genes on right
    for i, p in enumerate(matz.index):
        tg = nmf_resp.loc[nmf_resp["program"] == p, "top_genes"]
        if len(tg):
            s = ", ".join(tg.iloc[0].split(",")[:4]).strip()
            ax_hm.text(matz.shape[1] + 0.3, i, s, va="center", ha="left",
                        fontsize=7.5, color="#333333", style="italic")
    # highlight sig programs
    for i, p in enumerate(matz.index):
        pv = nmf_resp.loc[nmf_resp.program == p, "p_value"].iloc[0]
        if pv < 0.05:
            ax_hm.add_patch(plt.Rectangle((-0.5, i - 0.5),
                                           matz.shape[1], 1,
                                           fill=False, edgecolor="#FFA500",
                                           linewidth=2.2))
    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_label("program activity  (row z-score)", fontsize=8.5)
    # Legend for responders vs non-responders
    legend_elem = [Patch(color=PAL["Responder"], label="Responder"),
                    Patch(color=PAL["Non-Responder"], label="Non-Responder"),
                    Patch(facecolor="none", edgecolor="#FFA500", linewidth=2,
                          label="p < 0.05 (sc-atlas)")]
    ax_hm.legend(handles=legend_elem, loc="upper left",
                  bbox_to_anchor=(0.0, -0.12), ncol=3,
                  frameon=False, fontsize=8)
    fig.suptitle("NMF program activity in Visium sections "
                 "(rows sorted by Δ(R – NR); top genes at right)",
                 fontsize=11, fontweight="bold")
    _save(fig, "05e_nmf_spatial_heatmap")


# ══════════════════════════════════════════════════════════════════════
# 4. NMF bulk vs single-cell concordance
# ══════════════════════════════════════════════════════════════════════
if nmf_bulk is not None:
    # rank programs by |effect size| in sc-atlas
    rank_sc = nmf_resp.set_index("program")["effect_size"]
    # bulk per-sample mean by program
    bulk_mean = nmf_bulk.mean(axis=0)
    common = rank_sc.index.intersection(bulk_mean.index)
    x = rank_sc.loc[common].values
    y = bulk_mean.loc[common].values
    rho, pv = spearmanr(x, y)
    # concordance on programs that actually have bulk signal
    non_zero = [p for p in common if bulk_mean[p] > 0.01]
    if len(non_zero) >= 4:
        rho_nz, pv_nz = spearmanr(rank_sc.loc[non_zero].values,
                                     bulk_mean.loc[non_zero].values)
    else:
        rho_nz, pv_nz = np.nan, np.nan

    fig, ax = plt.subplots(figsize=(6.0, 5.4), constrained_layout=True)
    xy = [(p, rank_sc[p], bulk_mean[p]) for p in common]
    for p, x, y in xy:
        col = PAL["Responder"] if x > 0 else PAL["Non-Responder"]
        ax.scatter(x, y, s=140, color=col, edgecolor="black",
                    linewidth=0.9, alpha=0.92, zorder=3)
    # label offsets — alternate up/down for near-zero-y clusters
    ordered = sorted(xy, key=lambda t: (round(t[2], 2), t[1]))
    toggle = 1
    for p, x, y in ordered:
        if abs(y) < 0.01:     # programs that vanish in bulk
            dy = 12 if toggle > 0 else -14; toggle *= -1
        else:
            dy = 8
        dx = 8 if x >= 0 else -8
        ha = "left" if x >= 0 else "right"
        ax.annotate(p, (x, y), xytext=(dx, dy),
                     textcoords="offset points", fontsize=9,
                     ha=ha, color="#333333", fontweight="bold")
    ax.axvline(0, color="grey", linewidth=0.7, linestyle=":")
    ax.set_xlabel("sc-atlas effect size  (Responder − Non-Responder)")
    ax.set_ylabel("Bulk GSE202069  mean program activity")
    ttl = f"Cross-modality NMF concordance  (all ρ = {rho:+.2f}, p = {pv:.2g}"
    if not np.isnan(rho_nz):
        ttl += f";  nonzero ρ = {rho_nz:+.2f}, p = {pv_nz:.2g}"
    ttl += ")"
    ax.set_title(ttl, fontsize=10, fontweight="bold")
    _save(fig, "05e_nmf_bulk_vs_sc_concordance")


# ══════════════════════════════════════════════════════════════════════
# 5. NMF top-gene dot plot — biological interpretation
# ══════════════════════════════════════════════════════════════════════
top_genes_per_prog = {}
for p in nmf_loadings.columns:
    top_genes_per_prog[p] = nmf_loadings[p].nlargest(8).index.tolist()

all_genes = []
for p in nmf_loadings.columns:
    all_genes.extend(top_genes_per_prog[p])
# keep order but dedupe
seen = set(); ordered_genes = []
# order programs by sc effect size (descending; R-high at top)
prog_order = nmf_resp.sort_values("effect_size", ascending=False)["program"].tolist()
for p in prog_order:
    for g in top_genes_per_prog[p]:
        if g not in seen:
            seen.add(g); ordered_genes.append(g)

dot = nmf_loadings.loc[ordered_genes, prog_order]
# max-scale each gene
dot_z = dot.div(dot.max(axis=1) + 1e-9, axis=0)
fig, ax = plt.subplots(figsize=(max(5, 0.45 * len(prog_order)),
                                  max(5, 0.22 * len(ordered_genes))),
                          constrained_layout=True)
im = ax.imshow(dot_z.values, aspect="auto", cmap=VIR, vmin=0, vmax=1)
ax.set_xticks(range(len(prog_order))); ax.set_xticklabels(prog_order,
                                                            rotation=30, ha="right")
ax.set_yticks(range(len(ordered_genes))); ax.set_yticklabels(ordered_genes,
                                                                fontsize=7)
# colour program labels by R-vs-NR sign
for i, p in enumerate(prog_order):
    ef = nmf_resp.set_index("program").loc[p, "effect_size"]
    col = PAL["Responder"] if ef > 0 else PAL["Non-Responder"]
    ax.get_xticklabels()[i].set_color(col)
    ax.get_xticklabels()[i].set_fontweight("bold")
cb = fig.colorbar(im, ax=ax, shrink=0.4, pad=0.02)
cb.set_label("row-max scaled loading", fontsize=8.5)
ax.set_title("NMF programs — top 8 genes per program",
             fontsize=11, fontweight="bold")
_save(fig, "05e_nmf_top_genes_dot")


# ══════════════════════════════════════════════════════════════════════
# 6. ROC overlay — redo cleanly from saved features
# ══════════════════════════════════════════════════════════════════════
df = feats.merge(nmf_patient[["patient_id"] + [f"NMF{i+1}" for i in range(10)]]
                   .rename(columns={f"NMF{i+1}": f"nmf_NMF{i+1}"
                                      for i in range(10)}),
                   on="patient_id", how="left")
y = (df["response"] == "Responder").astype(int).values

layers = {
    "L1_tumour_intrinsic": [c for c in df.columns if c.startswith("cnv_")],
    "L2_regulon":          [c for c in df.columns if c.startswith("regulon_")],
    "L3_composition":      [c for c in df.columns if c.startswith("frac_")],
    "L4_io_signatures":    [c for c in df.columns if c.startswith("sig_")],
    "L5_nmf_programs":     [c for c in df.columns if c.startswith("nmf_")],
}
# Use the top-3 features per layer we already selected
top_cols = sel_feats.groupby("layer")["feature"].apply(list).to_dict()

def scale_impute(X):
    X = pd.DataFrame(X)
    X = X.fillna(X.median()).fillna(0)
    return StandardScaler().fit_transform(X.values)

def loocv_proba(X, y):
    loo = LeaveOneOut(); p = np.zeros(len(y))
    for tr, te in loo.split(X):
        m = LogisticRegression(penalty="l2", C=1.0, max_iter=3000,
                                class_weight="balanced")
        m.fit(X[tr], y[tr])
        p[te] = m.predict_proba(X[te])[0, 1]
    return p

roc_data = {}
for L, feats_L in top_cols.items():
    X = scale_impute(df[feats_L].values)
    p = loocv_proba(X, y)
    auc = roc_auc_score(y, p)
    fpr, tpr, _ = roc_curve(y, p)
    roc_data[L] = (fpr, tpr, auc)

# ensemble
all_top = sum(top_cols.values(), [])
X_all = scale_impute(df[all_top].values)
p_all = loocv_proba(X_all, y)
auc_all = roc_auc_score(y, p_all)
fpr, tpr, _ = roc_curve(y, p_all)
roc_data["Balanced ensemble"] = (fpr, tpr, auc_all)

fig, ax = plt.subplots(figsize=(5.6, 5.4), constrained_layout=True)
# order by AUC desc for legend
order = sorted(roc_data.keys(),
                key=lambda k: (-roc_data[k][2], k != "Balanced ensemble"))
for name in order:
    fpr, tpr, auc = roc_data[name]
    c = LAYER_COL.get(name, "#666")
    lw = 3.0 if name == "Balanced ensemble" else 1.8
    ls = "-" if name == "Balanced ensemble" else "--"
    alpha = 1.0 if name == "Balanced ensemble" else 0.85
    ax.plot(fpr, tpr, color=c, linewidth=lw, linestyle=ls, alpha=alpha,
             label=f"{LAYER_LABEL.get(name, name)}   AUC = {auc:.2f}")
ax.plot([0, 1], [0, 1], color="grey", linewidth=0.8, linestyle=":")
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title("LOOCV ROC  —  balanced top-3 features per layer",
             fontsize=10.5, fontweight="bold")
ax.legend(loc="lower right", frameon=False, fontsize=8)
ax.set_aspect("equal")
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
_save(fig, "05e_roc_by_layer")


# ══════════════════════════════════════════════════════════════════════
# 7. AUC horizontal bar — clearer than ROC for headline
# ══════════════════════════════════════════════════════════════════════
aa = pd.DataFrame({"layer": list(roc_data.keys()),
                   "auc":   [v[2] for v in roc_data.values()]})
aa = aa.sort_values("auc")
fig, ax = plt.subplots(figsize=(6.5, 3.5), constrained_layout=True)
cols = [LAYER_COL[L] for L in aa["layer"]]
ax.barh([LAYER_LABEL.get(L, L) for L in aa["layer"]], aa["auc"],
         color=cols, edgecolor="black", linewidth=0.5, height=0.6)
ax.axvline(0.5, color="grey", linewidth=0.7, linestyle=":")
for i, v in enumerate(aa["auc"]):
    ax.text(v + 0.008, i, f"{v:.2f}", va="center", fontsize=9,
             fontweight="bold")
ax.set_xlim(0, 1.0)
ax.set_xlabel("LOOCV AUC")
ax.set_title("Layer-wise predictive value  (balanced features)",
              fontsize=10.5, fontweight="bold")
_save(fig, "05e_auc_bars")


# ══════════════════════════════════════════════════════════════════════
# 8. Variance partition — clearer add-layer cumulative story
# ══════════════════════════════════════════════════════════════════════
# cumulative AUC: rank layers by single-layer AUC desc; add one at a time.
single_auc = {L: roc_data[L][2] for L in layers}
rank = sorted(single_auc, key=single_auc.get, reverse=True)
cum = []
picked = []
for L in rank:
    picked.append(L)
    feats_cum = sum((top_cols[P] for P in picked), [])
    Xc = scale_impute(df[feats_cum].values)
    pc = loocv_proba(Xc, y)
    ac = roc_auc_score(y, pc)
    cum.append((L, ac))

cum_df = pd.DataFrame(cum, columns=["added_layer", "cum_auc"])
cum_df.to_csv(TAB / "balanced_cumulative_auc.csv", index=False)

fig, ax = plt.subplots(figsize=(7.5, 3.8), constrained_layout=True)
cum_labels = ["+ " + LAYER_LABEL[L].split("  ", 1)[-1]
               for L in cum_df["added_layer"]]
ax.plot(range(len(cum_df)), cum_df["cum_auc"], marker="o",
         markersize=9, linewidth=2.2, color="black", zorder=3)
for i, (L, a) in enumerate(zip(cum_df["added_layer"], cum_df["cum_auc"])):
    ax.scatter(i, a, s=160, color=LAYER_COL[L], edgecolor="black",
                linewidth=1.0, zorder=4)
    ax.text(i, a + 0.015, f"{a:.2f}", ha="center", va="bottom",
             fontsize=9, fontweight="bold")
ax.set_xticks(range(len(cum_df)))
ax.set_xticklabels(cum_labels, rotation=18, ha="right", fontsize=8.5)
ax.set_ylim(0.45, 1.02)
ax.axhline(0.5, color="grey", linewidth=0.6, linestyle=":")
ax.set_ylabel("LOOCV AUC")
ax.set_title("Cumulative gain as layers are added in AUC-rank order",
              fontsize=10.5, fontweight="bold")
_save(fig, "05e_variance_partition")


# ══════════════════════════════════════════════════════════════════════
# 9. Feature effects forest — clean labels, layer colour stripe
# ══════════════════════════════════════════════════════════════════════
sel = sel_feats.copy().sort_values("t_stat")
# prettier labels
def _pretty(row):
    name = row["feature"]
    name = name.replace("frac_", "% ").replace("sig_", "sig: ") \
                 .replace("regulon_", "reg: ").replace("nmf_", "") \
                 .replace("cnv_", "CNV: ")
    name = name.replace("_", " ")
    return name
sel["lbl"] = sel.apply(_pretty, axis=1)
# map face (R/NR direction) and edge (layer) colours
face_cols = [PAL["Responder"] if t > 0 else PAL["Non-Responder"] for t in sel["t_stat"]]
edge_cols = [LAYER_COL[L] for L in sel["layer"]]
fig, ax = plt.subplots(figsize=(10.5, max(4.5, 0.42 * len(sel))),
                        constrained_layout=True)
bars = ax.barh(sel["lbl"], sel["t_stat"],
                color=face_cols, edgecolor=edge_cols, linewidth=2.2,
                height=0.7, zorder=3)
ax.axvline(0, color="black", linewidth=0.8)
xmax = sel["t_stat"].abs().max()
# p-value labels anchored outside the bar tip so they never overlap
for i, (_, r) in enumerate(sel.iterrows()):
    pv = r["p_value"]
    stars = "***" if pv < 0.001 else "**" if pv < 0.01 else "*" if pv < 0.05 else ""
    xlab = r["t_stat"] + (0.15 * np.sign(r["t_stat"]) if r["t_stat"] != 0 else 0.15)
    ax.text(xlab, i, f"p={pv:.2g} {stars}",
             va="center", ha="left" if r["t_stat"] >= 0 else "right",
             fontsize=8.5)
ax.set_xlim(-xmax * 1.45, xmax * 1.55)
ax.set_xlabel("Welch's t-statistic   (positive → higher in Responders)")
ax.set_title("Top-3 features per layer  —  stratified by response  "
              "(bar edge colour = layer)",
              fontsize=11, fontweight="bold")
# legend of layer colour stripes
legend = [Patch(color=LAYER_COL[L], label=LAYER_LABEL[L]) for L in
           ["L1_tumour_intrinsic", "L2_regulon", "L3_composition",
             "L4_io_signatures", "L5_nmf_programs"]]
ax.legend(handles=legend, loc="lower right", frameon=False, fontsize=7.5,
           title="Layer", title_fontsize=8)
ax.tick_params(axis="y", labelsize=8.5)
_save(fig, "05e_feature_effects")

print("[05e] DONE — all polished figures in", FIG)
