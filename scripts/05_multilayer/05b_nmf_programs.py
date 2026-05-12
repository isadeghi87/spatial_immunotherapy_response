#!/usr/bin/env python3
"""
05b_nmf_programs.py — Unbiased NMF gene-expression programs that unify
single-cell, bulk, and spatial transcriptomics.

Strategy
--------
 1. On single-cell atlas: compute HVGs (seurat_v3, scikit-misc), run NMF with
    k sweep [5..20], select k by reconstruction-error elbow.
 2. For the selected k, export gene loadings (program × gene) and patient-level
    program activity (mean cell score per patient).
 3. Project programs onto:
    - Bulk RNA-seq (GSE202069) via matrix multiply on shared HVGs.
    - Spatial Visium (GSE238264) via spot-level scoring.
 4. Test each program against response (Mann-Whitney) and report effect sizes.

Outputs
-------
  out/05_multilayer/tables/nmf_programs_gene_loadings.csv
  out/05_multilayer/tables/nmf_programs_patient_activity.csv
  out/05_multilayer/tables/nmf_programs_spatial_activity.csv
  out/05_multilayer/tables/nmf_programs_bulk_activity.csv
  out/05_multilayer/tables/nmf_programs_response_stats.csv
  out/05_multilayer/figures/05b_nmf_elbow.png
  out/05_multilayer/figures/05b_nmf_program_loadings.png
  out/05_multilayer/figures/05b_nmf_program_vs_response.png
  out/05_multilayer/figures/05b_nmf_spatial_projection.png
"""
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
import seaborn as sns
import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.stats import mannwhitneyu
from sklearn.decomposition import NMF

from config.publication_style import apply_style, PAL, CMAP
apply_style()

OUT = PROJECT / "out" / "05_multilayer"
FIG = OUT / "figures"; TAB = OUT / "tables"
for d in (FIG, TAB): d.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("[05b] Unbiased NMF gene-expression programs")
print("=" * 70)

# ── Load atlas (use CNV-annotated smaller subset for speed) ─────────────
print("\n[LOAD] Reading CNV-annotated atlas …")
adata = sc.read_h5ad(PROJECT / "data" / "processed" / "integrated_io_atlas_cnv.h5ad")
print(f"   Shape: {adata.shape}")

# Filter labelled cells and keep Cappuyns held out from discovery.
mask = adata.obs["response"].isin(["Responder", "Non-Responder"])
adata = adata[mask].copy()
print(f"   Labelled subset: {adata.shape}")

if "dataset" not in adata.obs.columns:
    pid_col = "patient_ID" if "patient_ID" in adata.obs.columns else "sample_id"
    adata.obs["dataset"] = np.where(
        adata.obs[pid_col].astype(str).str.startswith("pat"),
        "Cappuyns2023", "GSE206325")

adata_all = adata.copy()
adata = adata[adata.obs["dataset"].astype(str) != "Cappuyns2023"].copy()
print(f"   Discovery subset: {adata.shape}")

N_TARGET = 20000
if adata.n_obs > N_TARGET:
    rng = np.random.default_rng(0)
    idx = rng.choice(adata.n_obs, N_TARGET, replace=False)
    adata = adata[idx].copy()
    print(f"   Down-sampled to: {adata.shape}")

# ── HVG selection on raw counts ────────────────────────────────────────
print("\n[HVG] Selecting highly variable genes …")
if adata.X.min() < 0:
    print("   [WARN] X already normalised/log; using existing X for NMF (non-negative clipped).")
    adata.X = np.clip(adata.X.toarray() if sparse.issparse(adata.X) else adata.X,
                       0, None)
try:
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",
                                  batch_key=None)
except Exception:
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
hvg = adata.var["highly_variable"].values
adata_hvg = adata[:, hvg].copy()
print(f"   HVGs kept: {adata_hvg.n_vars}")

# Ensure non-negative dense matrix
X = adata_hvg.X
if sparse.issparse(X):
    X = X.toarray()
X = np.clip(X, 0, None).astype(np.float32)

# ── k sweep ────────────────────────────────────────────────────────────
print("\n[SWEEP] NMF k ∈ {5,7,10,12,15,20} …")
k_values = [5, 7, 10, 12, 15, 20]
recon_errs = []
models = {}
for k in k_values:
    m = NMF(n_components=k, init="nndsvda", max_iter=400,
            random_state=0, tol=1e-3)
    W = m.fit_transform(X)
    H = m.components_
    err = m.reconstruction_err_
    recon_errs.append(err)
    models[k] = (m, W, H)
    print(f"   k={k:2d}  recon_err={err:.3f}")

# Elbow plot
fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
ax.plot(k_values, recon_errs, "o-", color=PAL["Responder"], linewidth=2,
        markersize=7)
ax.set_xlabel("Number of programs (k)")
ax.set_ylabel("Reconstruction error")
ax.set_title("NMF k sweep — reconstruction error elbow",
             fontsize=11, fontweight="bold")
for fmt in ("png", "pdf"):
    fig.savefig(FIG / f"05b_nmf_elbow.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)

# ── Pick k by second-derivative elbow ──────────────────────────────────
errs = np.asarray(recon_errs)
d2 = np.diff(np.diff(errs))
k_chosen = k_values[int(np.argmax(d2)) + 1] if len(d2) else 10
# clamp to sane range
if k_chosen not in models:
    k_chosen = 10 if 10 in models else k_values[len(k_values) // 2]
print(f"\n[CHOSEN] k = {k_chosen}")

m, W, H = models[k_chosen]
program_names = [f"NMF{i+1}" for i in range(k_chosen)]

# ── Gene loadings per program ──────────────────────────────────────────
print("\n[LOAD] Extracting top genes per program …")
gene_df = pd.DataFrame(H, index=program_names, columns=adata_hvg.var_names)
gene_df.T.to_csv(TAB / "nmf_programs_gene_loadings.csv")

top_genes = {}
for p in program_names:
    top_genes[p] = gene_df.loc[p].nlargest(20).index.tolist()
pd.DataFrame(top_genes).to_csv(TAB / "nmf_top_genes_per_program.csv", index=False)

# ── Project discovery-derived programs onto all labelled patients ──────
adata_all_hvg = adata_all[:, adata_hvg.var_names].copy()
X_all = adata_all_hvg.X
if sparse.issparse(X_all):
    X_all = X_all.toarray()
X_all = np.clip(X_all, 0, None).astype(np.float32)
W_all = m.transform(X_all)

# ── Cell × program activity → patient mean ─────────────────────────────
cell_W = pd.DataFrame(W_all, index=adata_all_hvg.obs_names, columns=program_names)
pid = "patient_ID" if "patient_ID" in adata_all.obs.columns else "sample_id"
cell_W["patient_id"] = adata_all.obs[pid].astype(str).values
cell_W["response"] = adata_all.obs["response"].values
cell_W["dataset"] = adata_all.obs["dataset"].astype(str).values

patient_activity = cell_W.groupby("patient_id")[program_names].mean().reset_index()
patient_resp = cell_W.groupby("patient_id")["response"].agg(
    lambda x: x.value_counts().idxmax()).reset_index()
patient_dataset = cell_W.groupby("patient_id")["dataset"].first().reset_index()
patient_activity = (patient_activity
                    .merge(patient_resp, on="patient_id")
                    .merge(patient_dataset, on="patient_id"))
patient_activity.to_csv(TAB / "nmf_programs_patient_activity.csv", index=False)
print(f"   Saved: nmf_programs_patient_activity.csv "
      f"({patient_activity.shape[0]} pts × {k_chosen} programs)")

# ── R-vs-NR test per program ───────────────────────────────────────────
patient_activity_disc = patient_activity.loc[
    patient_activity["dataset"].astype(str) != "Cappuyns2023"].copy()
rows = []
for p in program_names:
    r  = patient_activity_disc.loc[patient_activity_disc.response == "Responder",     p]
    nr = patient_activity_disc.loc[patient_activity_disc.response == "Non-Responder", p]
    if len(r) < 3 or len(nr) < 3:
        continue
    U, pv = mannwhitneyu(r, nr, alternative="two-sided")
    eff = (r.mean() - nr.mean()) / (np.concatenate([r, nr]).std() + 1e-12)
    rows.append(dict(program=p, n_R=len(r), n_NR=len(nr),
                      mean_R=float(r.mean()), mean_NR=float(nr.mean()),
                      effect_size=float(eff), p_value=float(pv),
                      top_genes=", ".join(top_genes[p][:8])))
prog_stats = pd.DataFrame(rows).sort_values("effect_size",
                                              key=lambda s: s.abs(),
                                              ascending=False)
prog_stats.to_csv(TAB / "nmf_programs_response_stats.csv", index=False)
print(prog_stats.to_string(index=False))

# ── Figure: program × response (effect size forest plot) ───────────────
fig, ax = plt.subplots(figsize=(9, max(3, k_chosen * 0.35)), constrained_layout=True)
prog_stats_sorted = prog_stats.iloc[::-1]
colors = [PAL["Responder"] if e > 0 else PAL["Non_Responder"]
          for e in prog_stats_sorted["effect_size"]]
ax.barh(prog_stats_sorted["program"], prog_stats_sorted["effect_size"],
        color=colors, edgecolor="black", linewidth=0.5, height=0.6)
ax.axvline(0, color="black", linewidth=0.8)
for i, (_, row) in enumerate(prog_stats_sorted.iterrows()):
    sig = "*" if row["p_value"] < 0.05 else ""
    ax.text(row["effect_size"] + 0.02 * np.sign(row["effect_size"] or 1),
             i, f"p={row['p_value']:.2g}{sig}",
             va="center", ha="left" if row["effect_size"] > 0 else "right",
             fontsize=7)
ax.set_xlabel("Standardised mean difference (R − NR)")
ax.set_title(f"NMF programs vs response (discovery cohort, k={k_chosen})",
             fontsize=11, fontweight="bold")
for fmt in ("png", "pdf"):
    fig.savefig(FIG / f"05b_nmf_program_vs_response.{fmt}",
                 dpi=300, bbox_inches="tight")
plt.close(fig)

# ── Figure: program loadings heatmap (top-N genes per program) ─────────
print("\n[FIG] Program loadings heatmap …")
top_n = 10
gene_set = []
for p in program_names:
    gene_set.extend(gene_df.loc[p].nlargest(top_n).index.tolist())
gene_set = list(dict.fromkeys(gene_set))  # preserve order, dedupe
load_mat = gene_df[gene_set].T
# row-z-score
load_z = (load_mat.sub(load_mat.mean(axis=1), axis=0)
           .div(load_mat.std(axis=1) + 1e-9, axis=0))
fig, ax = plt.subplots(figsize=(max(6, k_chosen * 0.7), max(8, len(gene_set) * 0.12)),
                          constrained_layout=True)
sns.heatmap(load_z, cmap=CMAP["diverging"], center=0, vmin=-2.5, vmax=2.5,
            xticklabels=True, yticklabels=True,
            cbar_kws=dict(label="z-score", shrink=0.5), ax=ax, linewidths=0)
ax.set_title(f"NMF program loadings (top {top_n} genes per program)",
             fontsize=11, fontweight="bold")
ax.tick_params(axis="y", labelsize=6)
ax.tick_params(axis="x", labelsize=9)
for fmt in ("png", "pdf"):
    fig.savefig(FIG / f"05b_nmf_program_loadings.{fmt}",
                dpi=300, bbox_inches="tight")
plt.close(fig)

# ── Bulk projection (GSE202069) ─────────────────────────────────────────
print("\n[BULK] Projecting NMF programs onto GSE202069 …")
bulk_dir = PROJECT / "data" / "pre_treatment" / "GSE202069_Zhu2023_antiPDL1"
try:
    # search across common extensions
    cands = []
    for pat in ("*counts*", "*expression*", "*TPM*", "*tpm*",
                 "*FPKM*", "*Expr*", "*matrix*"):
        cands.extend(bulk_dir.rglob(pat))
    cands = [f for f in cands if f.suffix.lower() in (".tsv", ".txt",
                                                        ".csv", ".gz")]
    print(f"   Candidates: {[c.name for c in cands]}")
    if not cands:
        raise FileNotFoundError("No bulk counts file found")
    b_candidate = cands[0]
    bulk = pd.read_csv(b_candidate, sep=None, engine="python", index_col=0)
    print(f"   Bulk: {bulk.shape[0]} genes × {bulk.shape[1]} samples")
    # log1p normalise
    bulk = np.log1p(bulk.astype(float))
    shared = [g for g in adata_hvg.var_names if g in bulk.index]
    print(f"   Shared HVGs: {len(shared)}")
    H_shared = gene_df[shared].values   # (k, gs)
    X_bulk = bulk.loc[shared].values    # (gs, samples)
    # Solve W_bulk >=0 minimising ||H^T W_bulk - X_bulk||^2 via nnls per sample
    from scipy.optimize import nnls
    W_bulk = np.zeros((X_bulk.shape[1], H_shared.shape[0]))
    for j in range(X_bulk.shape[1]):
        w, _ = nnls(H_shared.T, X_bulk[:, j])
        W_bulk[j] = w
    bulk_df = pd.DataFrame(W_bulk, index=bulk.columns, columns=program_names)
    bulk_df.to_csv(TAB / "nmf_programs_bulk_activity.csv")
    print(f"   Saved bulk activity: {bulk_df.shape}")
except Exception as e:
    print(f"   [SKIP] Bulk projection failed: {e}")

# ── Spatial projection (Visium GSE238264) ──────────────────────────────
print("\n[SPATIAL] Projecting NMF programs onto Visium …")
sp_dir = PROJECT / "data" / "spatial_data" / "GSE238264_HCC_IO_Visium"
try:
    sp_files = sorted(sp_dir.glob("*.h5ad"))
    print(f"   Found {len(sp_files)} Visium h5ad files")
    rows = []
    for f in sp_files:
        try:
            sp = sc.read_h5ad(f)
            sp.var_names_make_unique()
            # normalise if necessary
            if sp.X.max() > 100:
                sc.pp.normalize_total(sp, target_sum=1e4)
                sc.pp.log1p(sp)
            shared = [g for g in adata_hvg.var_names if g in sp.var_names]
            if len(shared) < 50:
                continue
            H_shared = gene_df[shared].values
            Xs = sp[:, shared].X
            if sparse.issparse(Xs):
                Xs = Xs.toarray()
            Xs = np.clip(Xs, 0, None).astype(np.float32)
            # project each spot: W_spot = X_spot @ pinv(H)^T = X @ H_shared.T @ inv(HH^T)
            HH = H_shared @ H_shared.T
            Ws = Xs @ H_shared.T @ np.linalg.pinv(HH)
            Ws = np.clip(Ws, 0, None)
            # section-level mean
            tag = f.stem.replace(".h5ad", "")
            for k_idx, p in enumerate(program_names):
                rows.append(dict(section=tag, program=p,
                                   mean=float(Ws[:, k_idx].mean()),
                                   n_spots=int(Xs.shape[0])))
        except Exception as e:
            print(f"   [warn] {f.name}: {e}")
    sp_df = pd.DataFrame(rows)
    if len(sp_df):
        sp_df.to_csv(TAB / "nmf_programs_spatial_activity.csv", index=False)
        # Response mapping
        RESP_MAP = {"HCC1R": "Responder", "HCC2R": "Responder", "HCC3R": "Responder",
                     "HCC4R": "Responder", "HCC5NR": "Non-Responder",
                     "HCC6NR": "Non-Responder", "HCC7NR": "Non-Responder"}
        sp_df["response"] = sp_df["section"].str.extract(r"(HCC\d+(?:NR|R))")[0].map(RESP_MAP)
        # heatmap: programs × sections
        mat = sp_df.pivot_table(index="program", columns="section",
                                  values="mean", observed=True)
        sec_order = sorted(mat.columns, key=lambda s: (RESP_MAP.get(s, "Z"), s))
        mat = mat[sec_order]
        fig, ax = plt.subplots(figsize=(max(7, len(sec_order) * 0.8),
                                          max(5, k_chosen * 0.35)),
                                 constrained_layout=True)
        sns.heatmap(mat, cmap=CMAP["sequential"], ax=ax,
                     cbar_kws=dict(label="Mean program activity",
                                     shrink=0.6),
                     xticklabels=True, yticklabels=True, linewidths=0.3)
        # annotate R / NR under x-axis
        for i, s in enumerate(sec_order):
            resp = RESP_MAP.get(s, "?")
            col = PAL["Responder"] if resp == "Responder" else PAL["Non_Responder"]
            ax.get_xticklabels()[i].set_color(col)
        ax.set_title("NMF program activity across Visium sections "
                     "(blue=R, red=NR)", fontsize=11, fontweight="bold")
        for fmt in ("png", "pdf"):
            fig.savefig(FIG / f"05b_nmf_spatial_projection.{fmt}",
                         dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"   Saved spatial projection: {len(sp_df)} rows")
except Exception as e:
    print(f"   [SKIP] Spatial projection failed: {e}")

print("\n[DONE] NMF programs ready.")
