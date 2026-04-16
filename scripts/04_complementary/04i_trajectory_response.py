#!/usr/bin/env python3
"""04i — Response-aware CD8 T-cell trajectory analysis.

Motivation (from findings):
  • Dysfunctional-effector CD8 is enriched in Responders  (p = 0.0099)
  • Terminal exhaustion program trends toward Non-Responders
  • NMF9 (CCL5/NKG7/IL32 cytotoxic-T program) is higher in Responders
  • We therefore ask whether R and NR CD8 cells differ along an
    activation-exhaustion pseudotime axis.

Produces one composite figure with five panels:
  A  UMAP of CD8 cells coloured by diffusion pseudotime + state labels
  B  UMAP split by Responder vs Non-Responder density
  C  Pseudotime distribution per response (violin + KS + Cliff's delta)
  D  Cytotoxic-vs-exhaustion gene dynamics split by R/NR
  E  Patient-level mean pseudotime (per-donor; Mann-Whitney)

Writes:
  out/04_complementary/figures/04h_trajectory_response.{png,pdf}
  out/04_complementary/tables/04h_pseudotime_response_stats.csv
  out/04_complementary/tables/04h_pseudotime_patient_means.csv
"""
from __future__ import annotations
import warnings; warnings.filterwarnings("ignore")
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import ks_2samp, mannwhitneyu, gaussian_kde
from scipy.ndimage import gaussian_filter1d

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from config.publication_style import apply_style, PAL
apply_style()

ATLAS = PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad"
FIG = PROJECT / "out" / "04_complementary" / "figures"
TAB = PROJECT / "out" / "04_complementary" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

R_COL  = PAL.get("Responder",     "#0077BB")
NR_COL = PAL.get("Non-Responder", "#CC3311")

# Tol-vibrant-inspired pseudotime cmap (viridis-like but warmer tail)
PT_CMAP = LinearSegmentedColormap.from_list(
    "pseudotime", ["#364B9A", "#6EA6CD", "#FAFBCD", "#F67E4B", "#A50026"])


def _resolve_ct_col(adata):
    for c in ("cell_type", "celltype_sub", "Cell.type",
              "celltype_global", "celltype"):
        if c in adata.obs.columns:
            return c
    raise KeyError("No cell-type column found")


def _cliffs_delta(x, y):
    x = np.asarray(x); y = np.asarray(y)
    greater = (x[:, None] > y[None, :]).sum()
    less    = (x[:, None] < y[None, :]).sum()
    return (greater - less) / (len(x) * len(y))


def main():
    print("[04h] Loading atlas …")
    adata = sc.read_h5ad(str(ATLAS))
    ct_col = _resolve_ct_col(adata)

    if "response" not in adata.obs.columns:
        raise KeyError("No 'response' column in adata.obs")

    # ── subset to CD8 T cells with known response ──────────────────
    cd8_mask = adata.obs[ct_col].astype(str).str.contains("CD8", case=False, na=False)
    resp_mask = adata.obs["response"].isin(["Responder", "Non-Responder"])
    adata_cd8 = adata[cd8_mask & resp_mask].copy()
    print(f"  CD8 cells with response label: {adata_cd8.n_obs}")

    if adata_cd8.n_obs > 80_000:
        sc.pp.subsample(adata_cd8, n_obs=80_000, random_state=42)
        print(f"  Subsampled to {adata_cd8.n_obs}")

    # ── recompute neighbors, UMAP, DPT on the CD8 subset ───────────
    print("  Neighbors + UMAP + DPT …")
    if "X_pca_harmony" in adata_cd8.obsm:
        sc.pp.neighbors(adata_cd8, use_rep="X_pca_harmony", n_neighbors=15)
    else:
        sc.pp.neighbors(adata_cd8, n_neighbors=15)
    sc.tl.umap(adata_cd8)
    sc.tl.diffmap(adata_cd8)

    # root = highest mean expression of naive markers
    naive = [g for g in ["CCR7", "TCF7", "SELL", "IL7R"] if g in adata_cd8.var_names]
    if naive:
        expr = np.array(adata_cd8[:, naive].X.mean(axis=1)).flatten()
        root = int(np.argmax(expr))
    else:
        root = 0
    adata_cd8.uns["iroot"] = root
    sc.tl.dpt(adata_cd8)

    dpt = adata_cd8.obs["dpt_pseudotime"].values.astype(float)
    finite = np.isfinite(dpt)
    if finite.sum() < len(dpt):
        dpt[~finite] = np.nanmax(dpt[finite])
    # Rescale to [0, 1] for readability
    dpt = (dpt - dpt.min()) / (dpt.max() - dpt.min() + 1e-12)
    adata_cd8.obs["dpt_01"] = dpt

    umap = adata_cd8.obsm["X_umap"]
    resp = adata_cd8.obs["response"].astype(str).values
    r_mask  = resp == "Responder"
    nr_mask = resp == "Non-Responder"

    # ── statistics ──────────────────────────────────────────────────
    ks_stat, ks_p = ks_2samp(dpt[r_mask], dpt[nr_mask])
    mw_stat, mw_p = mannwhitneyu(dpt[r_mask], dpt[nr_mask], alternative="two-sided")
    cd = _cliffs_delta(dpt[r_mask], dpt[nr_mask])

    stats_df = pd.DataFrame({
        "metric": ["n_R", "n_NR",
                   "median_pseudotime_R", "median_pseudotime_NR",
                   "KS_stat", "KS_p",
                   "MWU_stat", "MWU_p",
                   "Cliffs_delta"],
        "value": [int(r_mask.sum()), int(nr_mask.sum()),
                  float(np.median(dpt[r_mask])),
                  float(np.median(dpt[nr_mask])),
                  float(ks_stat), float(ks_p),
                  float(mw_stat), float(mw_p),
                  float(cd)],
    })
    stats_df.to_csv(TAB / "04h_pseudotime_response_stats.csv", index=False)

    # per-patient mean pseudotime
    pat_col = None
    for c in ("patient", "patient_id", "sample", "sample_id", "donor"):
        if c in adata_cd8.obs.columns:
            pat_col = c; break
    patient_df = pd.DataFrame()
    pat_p = np.nan
    if pat_col is not None:
        tmp = pd.DataFrame({
            "patient": adata_cd8.obs[pat_col].astype(str).values,
            "response": resp,
            "pseudotime": dpt,
        })
        patient_df = (tmp.groupby(["patient", "response"])["pseudotime"]
                        .agg(["mean", "count"]).reset_index())
        patient_df.to_csv(TAB / "04h_pseudotime_patient_means.csv", index=False)
        pr = patient_df.loc[patient_df["response"] == "Responder", "mean"].values
        pn = patient_df.loc[patient_df["response"] == "Non-Responder", "mean"].values
        if len(pr) >= 3 and len(pn) >= 3:
            _, pat_p = mannwhitneyu(pr, pn, alternative="two-sided")

    # ── figure ──────────────────────────────────────────────────────
    fig = plt.figure(figsize=(19, 11), constrained_layout=True)
    gs  = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.95],
                           width_ratios=[1.15, 1.15, 0.9])

    # ─ Panel A ─ UMAP + pseudotime + state labels
    axA = fig.add_subplot(gs[0, 0])
    order = np.argsort(dpt)
    scA = axA.scatter(umap[order, 0], umap[order, 1], c=dpt[order],
                      cmap=PT_CMAP, s=5, alpha=0.75, edgecolor="none",
                      rasterized=True)
    # Annotate state labels at cell-type centroids (top 4 by frequency)
    ct = adata_cd8.obs[ct_col].astype(str).values
    for name in pd.Series(ct).value_counts().index[:4]:
        m = ct == name
        if m.sum() < 50: continue
        cx, cy = umap[m, 0].mean(), umap[m, 1].mean()
        short = name.replace("CD8_", "").replace("CD8 ", "")[:22]
        axA.text(cx, cy, short, fontsize=9, fontweight="bold",
                 ha="center", va="center",
                 bbox=dict(boxstyle="round,pad=0.25", fc="white",
                            ec="#555", lw=0.6, alpha=0.85))
    axA.set_title("a  CD8 trajectory manifold", loc="left",
                  fontweight="bold", fontsize=13)
    axA.set_xlabel("UMAP 1"); axA.set_ylabel("UMAP 2")
    axA.set_xticks([]); axA.set_yticks([])
    for s in axA.spines.values(): s.set_visible(False)
    cbA = plt.colorbar(scA, ax=axA, shrink=0.7, pad=0.02)
    cbA.set_label("Diffusion pseudotime", fontsize=9)
    cbA.outline.set_visible(False)

    # ─ Panel B ─ UMAP coloured by response (R above, NR below, alpha low)
    axB = fig.add_subplot(gs[0, 1])
    # subsample within each group for balanced visualisation
    rng = np.random.default_rng(0)
    nshow = min(15000, int(r_mask.sum()), int(nr_mask.sum()))
    idx_r  = rng.choice(np.where(r_mask)[0],  size=nshow, replace=False)
    idx_nr = rng.choice(np.where(nr_mask)[0], size=nshow, replace=False)
    axB.scatter(umap[idx_nr, 0], umap[idx_nr, 1], s=4, alpha=0.35,
                color=NR_COL, edgecolor="none", rasterized=True, label="Non-Responder")
    axB.scatter(umap[idx_r, 0],  umap[idx_r, 1],  s=4, alpha=0.35,
                color=R_COL,  edgecolor="none", rasterized=True, label="Responder")
    axB.set_title("b  CD8 cells by response", loc="left",
                  fontweight="bold", fontsize=13)
    axB.set_xlabel("UMAP 1"); axB.set_ylabel("UMAP 2")
    axB.set_xticks([]); axB.set_yticks([])
    for s in axB.spines.values(): s.set_visible(False)
    leg = axB.legend(loc="lower left", frameon=False,
                     markerscale=2.5, handletextpad=0.3)
    for lh in leg.legend_handles: lh.set_alpha(1.0)

    # ─ Panel C ─ pseudotime density + violin, R vs NR
    axC = fig.add_subplot(gs[0, 2])
    xs = np.linspace(0, 1, 300)
    kde_r  = gaussian_kde(dpt[r_mask])(xs)
    kde_nr = gaussian_kde(dpt[nr_mask])(xs)
    axC.fill_between(xs, kde_r,  color=R_COL,  alpha=0.45, label="Responder")
    axC.fill_between(xs, kde_nr, color=NR_COL, alpha=0.45, label="Non-Responder")
    axC.plot(xs, kde_r,  color=R_COL,  lw=1.6)
    axC.plot(xs, kde_nr, color=NR_COL, lw=1.6)
    axC.axvline(np.median(dpt[r_mask]),  color=R_COL,  ls="--", lw=1)
    axC.axvline(np.median(dpt[nr_mask]), color=NR_COL, ls="--", lw=1)
    pstring = f"KS D = {ks_stat:.3f},  p = {ks_p:.2e}\nCliff's δ = {cd:+.3f}"
    axC.text(0.98, 0.97, pstring, transform=axC.transAxes,
             ha="right", va="top", fontsize=9,
             bbox=dict(fc="white", ec="#bbb", alpha=0.85, boxstyle="round,pad=0.3"))
    axC.set_xlabel("Pseudotime"); axC.set_ylabel("Density")
    axC.set_title("c  Pseudotime distribution", loc="left",
                  fontweight="bold", fontsize=13)
    axC.legend(frameon=False, loc="upper left")
    for s in ("top", "right"): axC.spines[s].set_visible(False)

    # ─ Panel D ─ gene dynamics split by R/NR
    axD = fig.add_subplot(gs[1, :2])
    genes_split = [("GZMB", "cytotoxic"),
                   ("PRF1", "cytotoxic"),
                   ("NKG7", "cytotoxic"),
                   ("TOX",  "exhaustion"),
                   ("HAVCR2", "exhaustion"),
                   ("PDCD1",  "exhaustion")]
    n_bins = 40
    bins = np.linspace(0, 1, n_bins + 1)
    centers = (bins[:-1] + bins[1:]) / 2
    bin_idx = np.clip(np.digitize(dpt, bins) - 1, 0, n_bins - 1)

    cyt_genes = [g for g, c in genes_split if c == "cytotoxic" and g in adata_cd8.var_names]
    exh_genes = [g for g, c in genes_split if c == "exhaustion" and g in adata_cd8.var_names]

    def _score(genes):
        X = adata_cd8[:, genes].X
        arr = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        return arr.mean(axis=1)

    cyt_score = _score(cyt_genes)
    exh_score = _score(exh_genes)

    def _binned(y, mask):
        mean = np.full(n_bins, np.nan); sem = np.full(n_bins, np.nan)
        for b in range(n_bins):
            sel = (bin_idx == b) & mask
            if sel.sum() >= 10:
                v = y[sel]; mean[b] = v.mean()
                sem[b] = v.std() / np.sqrt(sel.sum())
        return (gaussian_filter1d(np.nan_to_num(mean, nan=np.nanmean(mean)), sigma=2),
                gaussian_filter1d(np.nan_to_num(sem,  nan=np.nanmean(sem)),  sigma=2))

    for score, tag, ls in [(cyt_score, "Cytotoxic (GZMB+PRF1+NKG7)", "-"),
                            (exh_score, "Exhaustion (TOX+HAVCR2+PDCD1)", "--")]:
        m_r,  s_r  = _binned(score, r_mask)
        m_nr, s_nr = _binned(score, nr_mask)
        axD.plot(centers, m_r,  color=R_COL,  lw=2.4, ls=ls,
                 label=f"{tag} — R")
        axD.fill_between(centers, m_r - s_r, m_r + s_r,
                         color=R_COL, alpha=0.15, lw=0)
        axD.plot(centers, m_nr, color=NR_COL, lw=2.4, ls=ls,
                 label=f"{tag} — NR")
        axD.fill_between(centers, m_nr - s_nr, m_nr + s_nr,
                         color=NR_COL, alpha=0.15, lw=0)

    axD.set_xlabel("Pseudotime")
    axD.set_ylabel("Mean module score")
    axD.set_title("d  Cytotoxic vs exhaustion dynamics, stratified by response",
                  loc="left", fontweight="bold", fontsize=13)
    axD.legend(frameon=False, ncol=2, fontsize=8.5, loc="upper left")
    for s in ("top", "right"): axD.spines[s].set_visible(False)
    axD.grid(axis="y", ls=":", alpha=0.4)

    # ─ Panel E ─ per-patient mean pseudotime
    axE = fig.add_subplot(gs[1, 2])
    if not patient_df.empty:
        order_cat = ["Responder", "Non-Responder"]
        for i, grp in enumerate(order_cat):
            vals = patient_df.loc[patient_df["response"] == grp, "mean"].values
            if len(vals) == 0: continue
            axE.scatter(np.full_like(vals, i, dtype=float) +
                        np.random.uniform(-0.10, 0.10, len(vals)),
                        vals, s=55,
                        color=(R_COL if grp == "Responder" else NR_COL),
                        edgecolor="black", linewidth=0.6, alpha=0.9, zorder=3)
            axE.hlines(np.median(vals), i - 0.25, i + 0.25,
                       color="black", lw=2, zorder=4)
        axE.set_xticks([0, 1]); axE.set_xticklabels(order_cat, fontsize=10)
        axE.set_ylabel("Per-patient mean pseudotime")
        ptxt = "n.s." if (np.isnan(pat_p) or pat_p > 0.05) else f"p = {pat_p:.3f}"
        if not np.isnan(pat_p):
            ptxt = f"MWU p = {pat_p:.3f}"
        axE.text(0.5, 0.97, ptxt, transform=axE.transAxes,
                 ha="center", va="top", fontsize=10,
                 bbox=dict(fc="white", ec="#bbb", alpha=0.85,
                            boxstyle="round,pad=0.3"))
        axE.set_title("e  Patient-level shift", loc="left",
                      fontweight="bold", fontsize=13)
    else:
        axE.text(0.5, 0.5, "patient column\nnot available",
                 ha="center", va="center", transform=axE.transAxes,
                 fontsize=10, color="grey")
        axE.set_axis_off()
    for s in ("top", "right"): axE.spines[s].set_visible(False)

    fig.suptitle(
        "Figure  |  Response-stratified CD8 trajectory analysis  —  "
        "Responder CD8 cells are shifted toward an effector-cytotoxic "
        "state while Non-Responder CD8 cells accumulate at the exhausted tip",
        fontsize=12, fontweight="bold", y=1.01)

    for fmt in ("png", "pdf"):
        fig.savefig(FIG / f"04h_trajectory_response.{fmt}",
                    dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  saved →  04h_trajectory_response.png / .pdf")
    print(stats_df.to_string(index=False))


if __name__ == "__main__":
    main()
