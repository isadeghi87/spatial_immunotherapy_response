#!/usr/bin/env python3
"""02g_infercnv_spatial.py
InferCNV analysis on HCC Visium spatial transcriptomics samples using infercnvpy.

For each GSE238264 (pre/post-treatment IO) and GSE281759 (multi-region) Visium
sample, this script:
  1. Loads the spot-level h5ad.
  2. Assigns reference (immune/stromal) vs observation (hepatocyte/tumor) per
     spot using cell2location-derived annotations from prior analysis.
  3. Runs infercnvpy.tl.infercnv() with appropriate HMM smoothing.
  4. Writes:
       - Per-spot CNV summary (cnv_score, cnv_burd)
       - Per-arm CNV summary (mean correction factor per chr arm)
       - Genome-wide heatmap per sample
       - Updated integrated spatial AnnData with cnv keys
  5. Appends an entry to the run log CSV.

Outputs:
  out/02b_spatial_biomarkers/cnv/tables/spatial_infercnvpy_spot.csv
  out/02b_spatial_biomarkers/cnv/tables/spatial_infercnvpy_arm.csv
  out/02b_spatial_biomarkers/cnv/tables/spatial_infercnvpy_run_log.csv
  out/02b_spatial_biomarkers/cnv/figures/spatial_cnv_{sample_id}.png
"""

import warnings
warnings.filterwarnings("ignore")

import os
import sys
import re
import gc
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse
import scanpy as sc
import infercnvpy as cnv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT    = Path(__file__).resolve().parents[2]
DATA       = PROJECT / "data"
ANN_DIR    = PROJECT / "out" / "02b_spatial_biomarkers" / "tables"
CNV_DIR    = PROJECT / "out" / "02b_spatial_biomarkers" / "cnv"
TABLE_DIR  = CNV_DIR / "tables"
FIG_DIR    = CNV_DIR / "figures"
NPZ_DIR    = CNV_DIR / "npz"
GENOME_TSV = DATA / "references" / "gencode.v44.basic.annotation.gtf.gz"

GSE238264_REAL = Path("/data/projects/HE_ST_machine_learning/data/GEO_HCC_processed/GSE238264")
GSE281759_REAL = Path("/data/projects/HE_ST_machine_learning/data/GEO_HCC_processed/GSE281759")

for d in (TABLE_DIR, FIG_DIR, NPZ_DIR):
    d.mkdir(parents=True, exist_ok=True)

# ── Cell-type to reference/observation mapping ────────────────────────────────
GSE238264_REF_TYPES = {
    "CD8+ T cell", "CD4+ T cell", "Treg", "NK cell", "B cell",
    "Macrophage", "DC", "Endothelial", "Fibroblast",
}
GSE281759_REF_TYPES = {
    "CD8_T", "CD4_T", "Treg", "NK", "B_cell",
    "Macrophage", "DC", "Endothelial", "Fibroblast",
}

MIN_REF = 20
MIN_OBS = 20
LOW_CONF_REF = 50
LOW_CONF_OBS = 100
KNOWN_QC_OUTLIERS = {
    "HCC2R": "low Visium spot-QC pass fraction in upstream QC",
    "HCC7NR": "low reference-spot support expected for inferCNV",
}

# ── Sample manifests ──────────────────────────────────────────────────────────
GSE238264_SAMPLES = {
    "HCC1R":  "Responder",
    "HCC2R":  "Responder",
    "HCC3R":  "Responder",
    "HCC4R":  "Responder",
    "HCC5NR": "Non_Responder",
    "HCC6NR": "Non_Responder",
    "HCC7NR": "Non_Responder",
}


def _h5ad_path(dataset: str, sample_id: str) -> Path:
    if dataset == "GSE238264":
        return GSE238264_REAL / sample_id / "st" / f"{sample_id}.h5ad"
    return Path("/data/projects/HE_ST_machine_learning/data/GEO_HCC_processed") \
           / dataset / sample_id / "st" / f"{sample_id}.h5ad"


# ── Annotation loaders ────────────────────────────────────────────────────────
def load_gse238264_annotations() -> pd.DataFrame:
    path = ANN_DIR / "gse238264_celltype_annotations.csv"
    df = pd.read_csv(path, index_col=0)
    ct_col = "celltype" if "celltype" in df.columns else "cell_type"
    # Index like HCC1R_AAACAAGTATCTCCCA-1-0 → extract sample + barcode
    idx = pd.Series(df.index.astype(str), name="idx")
    sample_part = idx.str.split("_", n=1, expand=True)
    df["sample_id"] = sample_part[0].values
    raw_bc = sample_part[1].values  # e.g. AAACAAGTATCTCCCA-1-0
    # Strip trailing -N (concat index) if present: keep e.g. AAACAAGTATCTCCCA-1
    bc_clean = np.where(
        pd.Series(raw_bc).str.count("-") > 1,
        pd.Series(raw_bc).str.rsplit("-", n=1, expand=True)[0].values,
        raw_bc,
    )
    df["barcode"] = bc_clean
    df = df.rename(columns={ct_col: "cell_type"})
    return df[["sample_id", "barcode", "cell_type"]]


def load_gse281759_annotations() -> pd.DataFrame:
    path = ANN_DIR / "gse281759_celltype_annotations.csv"
    if not path.exists():
        return pd.DataFrame(columns=["sample_id", "barcode", "cell_type"])
    df = pd.read_csv(path)
    if "dominant_cell_type" not in df.columns:
        return pd.DataFrame(columns=["sample_id", "barcode", "cell_type"])
    return pd.DataFrame({
        "sample_id": df["sample_id"].astype(str),
        "barcode":   df["barcodes"].astype(str),
        "cell_type": df["dominant_cell_type"].astype(str),
    })


def discover_gse281759_samples() -> list:
    symlink_dir = DATA / "spatial_data" / "GSE281759_HCC_MultiRegion"
    if not symlink_dir.exists():
        return []
    samples = []
    for p in sorted(symlink_dir.iterdir()):
        if p.suffix == ".h5ad":
            sid = p.stem
            real = _h5ad_path("GSE281759", sid)
            if real.exists():
                samples.append(sid)
    return samples


# ── Gene-order helper ─────────────────────────────────────────────────────────
def build_gene_order(gtf_path: Path) -> pd.DataFrame:
    """Parse GENCODE GTF → gene_symbol, chromosome, start, end."""
    import gzip
    cache_path = TABLE_DIR / "gene_order_infercnvpy.tsv"
    if cache_path.exists():
        print("  Gene order: loading from cache")
        return pd.read_csv(cache_path, sep="\t")

    print(f"  Parsing GENCODE GTF: {gtf_path}")
    AUTOSOMES = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    records = []
    open_fn = gzip.open if str(gtf_path).endswith(".gz") else open
    with open_fn(str(gtf_path), "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0]
            if chrom not in AUTOSOMES:
                continue
            attrs = parts[8]
            gene_name = None
            for field in attrs.split(";"):
                field = field.strip()
                if field.startswith("gene_name"):
                    gene_name = field.split('"')[1]
                    break
            if gene_name is None:
                continue
            records.append((gene_name, chrom, int(parts[3]), int(parts[4])))

    df = (pd.DataFrame(records, columns=["gene_symbol", "chromosome", "start", "end"])
            .drop_duplicates("gene_symbol", keep="first")
            .sort_values(["chromosome", "start"])
            .reset_index(drop=True))
    df.to_csv(cache_path, sep="\t", index=False)
    print(f"  Gene order: {len(df)} genes cached to {cache_path}")
    return df


def infer_confidence(sample_id: str, n_ref: int, n_obs: int, status: str) -> tuple[str, str]:
    """Return a conservative confidence class and explanation for spatial CNV."""
    if status != "completed":
        return "skipped", "sample did not complete inferCNV"
    reasons = []
    if n_ref < LOW_CONF_REF:
        reasons.append(f"low reference spots ({n_ref} < {LOW_CONF_REF})")
    if n_obs < LOW_CONF_OBS:
        reasons.append(f"low observation spots ({n_obs} < {LOW_CONF_OBS})")
    if sample_id in KNOWN_QC_OUTLIERS:
        reasons.append(KNOWN_QC_OUTLIERS[sample_id])
    if reasons:
        return "low_confidence", "; ".join(reasons)
    return "pass", "reference/observation support adequate"


def render_spatial_cnv_figure(adata, sample_id: str, fig_path: Path) -> str:
    """Render a robust per-sample spatial CNV diagnostic without infercnvpy ax conflicts."""
    if "X_cnv" not in adata.obsm:
        return "skipped_no_X_cnv"
    Xcnv = adata.obsm["X_cnv"]
    X = Xcnv.toarray() if scipy.sparse.issparse(Xcnv) else np.asarray(Xcnv)
    if X.size == 0:
        return "skipped_empty_X_cnv"

    rng = np.random.default_rng(0)
    score = (
        adata.obs["cnv_score"].to_numpy(dtype=float)
        if "cnv_score" in adata.obs.columns
        else np.nanmean(np.abs(X), axis=1)
    )
    order_df = pd.DataFrame({
        "idx": np.arange(adata.n_obs),
        "group": adata.obs["cnv_group"].astype(str).values,
        "score": score,
    }).sort_values(["group", "score"], ascending=[True, True])
    max_rows = 1600
    if len(order_df) > max_rows:
        keep = np.sort(rng.choice(order_df["idx"].values, max_rows, replace=False))
        order = order_df[order_df["idx"].isin(keep)]["idx"].to_numpy()
    else:
        order = order_df["idx"].to_numpy()

    X_plot = X[order, :]
    vmax = float(np.nanpercentile(np.abs(X_plot), 99)) if np.isfinite(X_plot).any() else 1.0
    vmax = max(vmax, 1e-3)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), constrained_layout=True)
    ax = axes[0]
    im = ax.imshow(X_plot, aspect="auto", interpolation="nearest", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_title(f"{sample_id} CNV window heatmap", fontsize=11, fontweight="bold")
    ax.set_xlabel("Genomic windows")
    ax.set_ylabel(f"Spots (ordered; n={len(order):,})")
    ax.set_yticks([])
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("Smoothed CNV deviation", fontsize=8)

    ax = axes[1]
    if "spatial" in adata.obsm:
        xy = adata.obsm["spatial"]
        sca = ax.scatter(
            xy[:, 0], xy[:, 1], c=score, s=8, cmap="magma", edgecolors="none", rasterized=True
        )
        ax.invert_yaxis()
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"{sample_id} spatial CNV burden", fontsize=11, fontweight="bold")
        ax.set_xticks([])
        ax.set_yticks([])
        cb2 = fig.colorbar(sca, ax=ax, fraction=0.035, pad=0.02)
        cb2.set_label("CNV score", fontsize=8)
    else:
        ax.hist(score[np.isfinite(score)], bins=40, color="#6b7280", edgecolor="white", linewidth=0.4)
        ax.set_title(f"{sample_id} CNV score distribution", fontsize=11, fontweight="bold")
        ax.set_xlabel("CNV score")
        ax.set_ylabel("Spots")

    fig.suptitle("Spatial inferCNV diagnostic (exploratory)", fontsize=12, fontweight="bold")
    fig.savefig(str(fig_path.with_suffix(".png")), dpi=600, facecolor="white")
    fig.savefig(str(fig_path.with_suffix(".pdf")), facecolor="white")
    plt.close(fig)
    return "rendered_fallback_heatmap_spatial"


# ── Per-sample runner ─────────────────────────────────────────────────────────
def run_sample(
    sample_id: str,
    dataset: str,
    response: str,
    h5ad_path: Path,
    ann_df: pd.DataFrame,          # filtered annotations for this sample
    ref_types: set,
    gene_order: pd.DataFrame,
    overwrite: bool = False,
) -> dict:
    done_flag = TABLE_DIR / f".done_spatial_cnv_{sample_id}"
    if done_flag.exists() and not overwrite:
        npz_path = NPZ_DIR / f"{sample_id}_cnv.npz"
        if npz_path.exists():
            return {"sample_id": sample_id, "status": "already_done", "reason": "cached_completed",
                    "figure_status": "cached", "confidence_class": "cached"}
        done_flag.unlink(missing_ok=True)

    # --- Load AnnData ------
    try:
        adata = sc.read_h5ad(str(h5ad_path))
    except Exception as e:
        return {"sample_id": sample_id, "status": "failed", "reason": f"read h5ad: {e}"}

    # gene symbols
    if adata.var.index[0].startswith("ENS"):
        if "gene_name" in adata.var.columns:
            adata.var_names = adata.var["gene_name"].astype(str).values
        elif "SYMBOL" in adata.var.columns:
            adata.var_names = adata.var["SYMBOL"].astype(str).values
    adata.var_names_make_unique()

    # --- Assign cell types -------------------------------------------------
    bc2ct = dict(zip(ann_df["barcode"].values, ann_df["cell_type"].values))
    barcodes = adata.obs_names.astype(str)

    def lookup_bc(bc):
        ct = bc2ct.get(bc)
        if ct is None:
            bc2 = bc.rsplit("-", 1)[0] if bc.count("-") > 1 else bc
            ct = bc2ct.get(bc2)
        return ct or "Unassigned"

    celltypes = np.array([lookup_bc(bc) for bc in barcodes])
    groups = np.where(np.isin(celltypes, list(ref_types)), "reference", "observation")

    n_ref = (groups == "reference").sum()
    n_obs = (groups == "observation").sum()
    if n_ref < MIN_REF or n_obs < MIN_OBS:
        done_flag.unlink(missing_ok=True)
        return {
            "sample_id": sample_id,
            "status": "skipped",
            "reason": f"insufficient: ref={n_ref} obs={n_obs}",
            "n_ref": int(n_ref),
            "n_obs": int(n_obs),
            "figure_status": "not_run_skipped",
            "confidence_class": "skipped",
            "confidence_reason": "insufficient reference or observation support",
        }

    adata.obs["cell_type"]   = celltypes
    adata.obs["cnv_group"]   = groups

    # --- Annotate genes with genomic position ------------------------------
    common_genes = adata.var_names[adata.var_names.isin(gene_order["gene_symbol"].values)]
    if len(common_genes) < 500:
        return {
            "sample_id": sample_id,
            "status": "skipped",
            "reason": f"too_few_common_genes: {len(common_genes)}",
        }
    adata = adata[:, common_genes].copy()

    go_sub = gene_order.set_index("gene_symbol").loc[common_genes]
    adata.var["chromosome"] = go_sub["chromosome"].values
    adata.var["start"]      = go_sub["start"].values
    adata.var["end"]        = go_sub["end"].values

    # Ensure raw counts
    if scipy.sparse.issparse(adata.X):
        maxval = adata.X.max()
    else:
        maxval = adata.X.max()
    if maxval <= 20:                  # likely normalised
        print(f"    [{sample_id}] WARNING: X may be normalised (max={maxval:.1f}), continuing")

    # --- Run infercnvpy -----------------------------------------------------
    print(f"    [{sample_id}] Running infercnvpy: ref={n_ref}, obs={n_obs}, genes={adata.n_vars}")
    try:
        cnv.tl.infercnv(
            adata,
            reference_key="cnv_group",
            reference_cat="reference",
            lfc_clip=3.0,
            window_size=100,
            step=10,
            exclude_chromosomes=None,
            n_jobs=4,
        )
        cnv.pp.neighbors(adata)
        cnv.tl.pca(adata)
        cnv.tl.leiden(adata)
        cnv.tl.cnv_score(adata)
    except Exception as e:
        return {"sample_id": sample_id, "status": "failed", "reason": f"infercnvpy: {e}"}

    # --- Save figure --------------------------------------------------------
    fig_stem = FIG_DIR / f"spatial_cnv_{sample_id}"
    figure_status = "not_run"
    try:
        figure_status = render_spatial_cnv_figure(adata, sample_id, fig_stem)
    except Exception as e:
        figure_status = f"failed: {e}"
        print(f"    [{sample_id}] Figure error (non-fatal): {e}")

    # --- Extract spot-level summary -----------------------------------------
    spot_cols = {}
    spot_cols["sample_id"] = sample_id
    spot_cols["dataset"]   = dataset
    spot_cols["response"]  = response
    spot_cols["barcode"]   = barcodes
    spot_cols["cell_type"] = celltypes
    spot_cols["group"]     = groups
    if "cnv_score" in adata.obs.columns:
        spot_cols["cnv_score"] = adata.obs["cnv_score"].values
    if "cnv_leiden" in adata.obs.columns:
        spot_cols["cnv_cluster"] = adata.obs["cnv_leiden"].values
    spot_df = pd.DataFrame(spot_cols)
    spot_csv = TABLE_DIR / f"{sample_id}_spatial_cnv_spot.csv"
    spot_df.to_csv(spot_csv, index=False)

    # --- Save per-spot CNV matrix (.npz) for universal heatmap (P1.2) -------
    if "X_cnv" in adata.obsm:
        Xcnv = adata.obsm["X_cnv"]
        Xcnv_dense = Xcnv.toarray() if scipy.sparse.issparse(Xcnv) else np.asarray(Xcnv)
        cnv_matrix_mode = "gene" if Xcnv_dense.shape[1] == adata.n_vars else "window"
        # chr_pos = mapping window-index -> chromosome (per infercnvpy uns)
        chr_pos = adata.uns.get("cnv", {}).get("chr_pos", {})
        # Persist gene-level positions as fallback
        gene_chr   = adata.var["chromosome"].astype(str).values
        gene_start = adata.var["start"].astype(np.int64).values
        gene_end   = adata.var["end"].astype(np.int64).values
        gene_names = adata.var_names.astype(str).values
        npz_path = NPZ_DIR / f"{sample_id}_cnv.npz"
        np.savez_compressed(
            str(npz_path),
            X_cnv       = Xcnv_dense.astype(np.float32),
            barcodes    = barcodes.astype("U64"),
            cell_types  = celltypes.astype("U32"),
            groups      = groups.astype("U16"),
            cnv_score   = adata.obs["cnv_score"].values.astype(np.float32) if "cnv_score" in adata.obs.columns else np.zeros(adata.n_obs, dtype=np.float32),
            cnv_leiden  = adata.obs["cnv_leiden"].astype(str).values.astype("U16") if "cnv_leiden" in adata.obs.columns else np.array(["0"] * adata.n_obs, dtype="U16"),
            chr_pos_keys   = np.array(list(chr_pos.keys()), dtype="U16"),
            chr_pos_values = np.array(list(chr_pos.values()), dtype=np.int64),
            gene_names  = gene_names.astype("U64"),
            gene_chr    = gene_chr,
            gene_start  = gene_start,
            gene_end    = gene_end,
            sample_id   = np.array(sample_id, dtype="U64"),
            dataset     = np.array(dataset, dtype="U16"),
            response    = np.array(str(response), dtype="U32"),
            cnv_matrix_mode = np.array(cnv_matrix_mode, dtype="U16"),
        )
        print(f"    [{sample_id}] Saved CNV npz: shape={Xcnv_dense.shape} → {npz_path.name}")
    else:
        cnv_matrix_mode = "missing"

    # --- Extract arm-level summary ------------------------------------------
    # Use the smoothed CNV matrix stored in adata.obsm["X_cnv"]
    if "X_cnv" in adata.obsm:
        Xcnv = adata.obsm["X_cnv"]  # spots × gene-windows
        if Xcnv.shape[1] != adata.n_vars:
            # infercnvpy v0.6 stores window-level (not gene-level) values; arm
            # slicing via var indices is not valid → skip arm CSV gracefully.
            print(f"    [{sample_id}] Skipping arm-level CSV: X_cnv has "
                  f"{Xcnv.shape[1]} windows vs {adata.n_vars} genes (window-based mode)")
        else:
            chr_arr = adata.var["chromosome"].values
            start_arr = adata.var["start"].values
            arm_records = []
            for chrom in np.unique(chr_arr):
                mask = chr_arr == chrom
                # Approximate arm: split at median start
                starts = start_arr[mask]
                if len(starts) == 0:
                    continue
                mid = np.median(starts)
                for arm_suffix, arm_mask in [("p", starts <= mid), ("q", starts > mid)]:
                    idx_g  = np.where(mask)[0][arm_mask]
                    if len(idx_g) < 5:
                        continue
                    arm_label = f"{chrom}{arm_suffix}"
                    arm_mean  = Xcnv[:, idx_g].mean(axis=1) if not scipy.sparse.issparse(Xcnv) \
                                else Xcnv[:, idx_g].toarray().mean(axis=1)
                    for bc_i, bv in enumerate(np.asarray(arm_mean).ravel()):
                        arm_records.append({
                            "sample_id": sample_id,
                            "barcode":   barcodes[bc_i],
                            "group":     groups[bc_i],
                            "arm":       arm_label,
                            "cnv_mean":  float(bv),
                        })
            if arm_records:
                arm_df = pd.DataFrame(arm_records)
                arm_csv = TABLE_DIR / f"{sample_id}_spatial_cnv_arm.csv"
                arm_df.to_csv(arm_csv, index=False)

    confidence_class, confidence_reason = infer_confidence(sample_id, int(n_ref), int(n_obs), "completed")
    done_flag.touch()
    return {
        "sample_id": sample_id,
        "status": "completed",
        "reason": "",
        "n_ref": int(n_ref),
        "n_obs": int(n_obs),
        "n_genes": int(adata.n_vars),
        "cnv_windows": int(adata.obsm["X_cnv"].shape[1]) if "X_cnv" in adata.obsm else 0,
        "cnv_matrix_mode": cnv_matrix_mode,
        "figure_status": figure_status,
        "confidence_class": confidence_class,
        "confidence_reason": confidence_reason,
    }


# ── Main ─────────────────────────────────────────────────────────────────────
def main(overwrite: bool = False):
    print("=== Spatial inferCNVpy pipeline ===")
    print(f"  infercnvpy {cnv.__version__}")

    gene_order = build_gene_order(GENOME_TSV)

    # Load annotation tables
    print("Loading cell-type annotations …")
    gse238264_ann = load_gse238264_annotations()
    gse281759_ann = load_gse281759_annotations()

    run_log = []

    # ── GSE238264 ───────────────────────────────────────────────────────────
    print("\n== GSE238264 (HCC IO Visium) ==")
    for sid, resp in GSE238264_SAMPLES.items():
        h5ad = _h5ad_path("GSE238264", sid)
        print(f"  [{sid}] {resp}")
        if not h5ad.exists():
            print(f"    MISSING: {h5ad}")
            run_log.append({"sample_id": sid, "dataset": "GSE238264", "response": resp,
                            "status": "missing", "reason": str(h5ad)})
            continue
        ann_sub = gse238264_ann[gse238264_ann["sample_id"] == sid]
        result = run_sample(
            sample_id=sid, dataset="GSE238264", response=resp,
            h5ad_path=h5ad, ann_df=ann_sub,
            ref_types=GSE238264_REF_TYPES,
            gene_order=gene_order, overwrite=overwrite,
        )
        result["dataset"] = "GSE238264"
        result["response"] = resp
        run_log.append(result)
        print(f"    → {result['status']}: {result['reason']}")
        gc.collect()

    # ── GSE281759 ───────────────────────────────────────────────────────────
    print("\n== GSE281759 (HCC MultiRegion) ==")
    gse281759_samples = discover_gse281759_samples()
    print(f"  Found {len(gse281759_samples)} samples")
    for sid in gse281759_samples:
        h5ad = _h5ad_path("GSE281759", sid)
        print(f"  [{sid}]")
        ann_sub = gse281759_ann[gse281759_ann["sample_id"] == sid] \
            if "sample_id" in gse281759_ann.columns else gse281759_ann
        result = run_sample(
            sample_id=sid, dataset="GSE281759", response="",
            h5ad_path=h5ad, ann_df=ann_sub,
            ref_types=GSE281759_REF_TYPES,
            gene_order=gene_order, overwrite=overwrite,
        )
        result["dataset"] = "GSE281759"
        result["response"] = ""
        run_log.append(result)
        print(f"    → {result['status']}: {result['reason']}")
        gc.collect()

    # ── Merge per-sample tables ─────────────────────────────────────────────
    spot_files = sorted(TABLE_DIR.glob("*_spatial_cnv_spot.csv"))
    arm_files  = sorted(TABLE_DIR.glob("*_spatial_cnv_arm.csv"))
    if spot_files:
        spot_all = pd.concat([pd.read_csv(f) for f in spot_files], ignore_index=True)
        spot_all.to_csv(TABLE_DIR / "spatial_infercnvpy_spot.csv", index=False)
        print(f"\nSpot table → {len(spot_all)} rows")
    if arm_files:
        arm_all = pd.concat([pd.read_csv(f) for f in arm_files], ignore_index=True)
        arm_all.to_csv(TABLE_DIR / "spatial_infercnvpy_arm.csv", index=False)
        print(f"Arm table  → {len(arm_all)} rows")

    # ── Run log ─────────────────────────────────────────────────────────────
    log_df = pd.DataFrame(run_log)
    log_df.to_csv(TABLE_DIR / "spatial_infercnvpy_run_log.csv", index=False)
    print("\n=== Summary ===")
    print(log_df["status"].value_counts().to_string())

    n_ok = log_df["status"].isin(["completed", "already_done"]).sum()
    if n_ok == 0:
        print("WARNING: No samples completed inferCNV.", file=sys.stderr)
        sys.exit(1)
    print(f"\n{n_ok} samples completed spatial inferCNVpy.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    main(overwrite=args.overwrite)
