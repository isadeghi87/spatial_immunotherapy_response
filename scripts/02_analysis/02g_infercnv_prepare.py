#!/usr/bin/env python3
"""02g_infercnv_prepare.py
Pre-process Visium h5ad files for the R inferCNV pipeline.

For each sample in GSE238264 (HCC1R–HCC7NR) and GSE281759, this script:
  1. Resolves spot cell-type assignments from cell2location annotation tables.
  2. Reads raw count matrices from the h5ad files.
  3. Writes per-sample inferCNV input files:
       matrix.mtx.gz   – raw counts (genes × spots, Market Exchange format)
       barcodes.tsv     – spot barcodes
       features.tsv     – gene symbols
       annotations.tsv  – barcode  cell_group (reference | observation)
  4. Writes a genome-wide gene-order file from the GENCODE GTF once.

Outputs land in:
  out/02b_spatial_biomarkers/cnv/infercnv_inputs/{sample_id}/...
  out/02b_spatial_biomarkers/cnv/infercnv_inputs/gene_order.tsv
"""

import warnings
warnings.filterwarnings("ignore")

import gzip, io, os, sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
import scanpy as sc

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT = Path(__file__).resolve().parents[2]
DATA    = PROJECT / "data"
OUT_ROOT = PROJECT / "out" / "02b_spatial_biomarkers" / "cnv" / "infercnv_inputs"
ANN_DIR  = PROJECT / "out" / "02b_spatial_biomarkers" / "tables"

GSE238264_REAL  = Path("/data/projects/HE_ST_machine_learning/data/GEO_HCC_processed/GSE238264")
GSE281759_REAL  = Path("/data/projects/HE_ST_machine_learning/data/GEO_HCC_processed/GSE281759")
GTF_PATH        = DATA / "references" / "gencode.v44.basic.annotation.gtf.gz"

# Cell types treated as normal/immune reference for inferCNV
# (everything else is the "observation" ≈ potentially malignant group)
GSE238264_REF_TYPES = {
    "CD8+ T cell", "CD4+ T cell", "Treg", "NK cell", "B cell",
    "Macrophage", "DC", "Endothelial", "Fibroblast",
}
GSE281759_REF_TYPES = {
    "CD8_T", "CD4_T", "Treg", "NK", "B_cell",
    "Macrophage", "DC", "Endothelial", "Fibroblast",
}

# Minimum spots in each group for inferCNV to be valid
MIN_REF = 20
MIN_OBS = 20

# ── Sample manifests ──────────────────────────────────────────────────────────
GSE238264_SAMPLES = {
    "HCC1R": "Responder",
    "HCC2R": "Responder",
    "HCC3R": "Responder",
    "HCC4R": "Responder",
    "HCC5NR": "Non_Responder",
    "HCC6NR": "Non_Responder",
    "HCC7NR": "Non_Responder",
}


def _gse238264_h5ad(sample_id: str) -> Path:
    return GSE238264_REAL / sample_id / "st" / f"{sample_id}.h5ad"


def _gse281759_h5ad(sample_id: str) -> Path:
    return GSE281759_REAL / sample_id / "st" / f"{sample_id}.h5ad"


# ── Gene-order file from GENCODE GTF ─────────────────────────────────────────
def build_gene_order(gtf_path: Path, out_path: Path) -> pd.DataFrame:
    if out_path.exists():
        print(f"  Gene-order file already exists: {out_path}")
        return pd.read_csv(out_path, sep="\t")

    print(f"  Parsing GENCODE GTF: {gtf_path}")
    records = []
    open_fn = gzip.open if str(gtf_path).endswith(".gz") else open
    with open_fn(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0]
            if not (chrom.startswith("chr") and chrom[3:].split("_")[0].lstrip("chrXYM").isdigit()
                    or chrom in ("chrX", "chrY")):
                # keep only chr1–chr22, chrX, chrY
                if chrom not in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]:
                    continue
            start = int(parts[3])
            attrs = parts[8]
            gene_name = None
            for field in attrs.split(";"):
                field = field.strip()
                if field.startswith("gene_name"):
                    gene_name = field.split('"')[1]
                    break
            if gene_name is None:
                continue
            records.append((gene_name, chrom, start, int(parts[4])))

    df = (pd.DataFrame(records, columns=["gene_symbol", "chr", "start", "end"])
            .drop_duplicates("gene_symbol", keep="first")
            .sort_values(["chr", "start"])
            .reset_index(drop=True))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  Gene-order: {len(df)} genes → {out_path}")
    return df


# ── Write 10X-style MTX input for inferCNV ───────────────────────────────────
def write_sample_inputs(
    sample_id: str,
    h5ad_path: Path,
    celltype_map: dict,      # barcode → cell_type string
    ref_types: set,
    out_dir: Path,
    gene_order_df: pd.DataFrame,
    overwrite: bool = False,
) -> dict:
    """
    Returns a status dict:  status ∈ {"completed","skipped","failed"},  reason
    """
    done_flag = out_dir / ".done"
    if done_flag.exists() and not overwrite:
        return {"sample_id": sample_id, "status": "already_done", "reason": ""}

    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Load h5ad --------------------------------------------------------
    try:
        adata = sc.read_h5ad(str(h5ad_path))
    except Exception as e:
        return {"sample_id": sample_id, "status": "failed", "reason": f"read h5ad: {e}"}

    # gene names
    if adata.var.index.name == "gene_ids" or adata.var.index[0].startswith("ENS"):
        # index is Ensembl IDs; try to find a gene_name column
        if "gene_name" in adata.var.columns:
            gene_names = adata.var["gene_name"].values
        elif "SYMBOL" in adata.var.columns:
            gene_names = adata.var["SYMBOL"].values
        else:
            gene_names = adata.var.index.values.astype(str)
    else:
        gene_names = adata.var.index.values.astype(str)

    barcodes = adata.obs.index.values.astype(str)

    # --- Match barcodes to cell types ------------------------------------
    # celltype_map keys might have sample prefix e.g. HCC1R_BARCODE-1-0
    # Try both with and without prefix
    def lookup(bc):
        ct = celltype_map.get(bc)
        if ct is None:  # try with numeric suffix stripped
            bc_short = bc.rsplit("-", 1)[0] if bc.count("-") > 1 else bc
            ct = celltype_map.get(bc_short)
        return ct or "Unassigned"

    celltypes = np.array([lookup(bc) for bc in barcodes])

    # --- Assign reference / observation ----------------------------------
    groups = np.where(np.isin(celltypes, list(ref_types)), "reference", "observation")
    n_ref = (groups == "reference").sum()
    n_obs = (groups == "observation").sum()

    if n_ref < MIN_REF or n_obs < MIN_OBS:
        return {
            "sample_id": sample_id,
            "status": "skipped",
            "reason": f"insufficient_spots: ref={n_ref}, obs={n_obs}",
        }

    print(f"    {sample_id}: {n_ref} reference, {n_obs} observation spots")

    # --- Keep only genes present in gene-order file ----------------------
    valid_genes_set = set(gene_order_df["gene_symbol"].values)
    gene_mask = np.isin(gene_names, list(valid_genes_set))
    if gene_mask.sum() < 500:
        return {
            "sample_id": sample_id,
            "status": "skipped",
            "reason": f"too_few_genes_in_gene_order: {gene_mask.sum()}",
        }

    X = adata.X
    if scipy.sparse.issparse(X):
        X = X.tocsr()
        X_filt = X[:, gene_mask]
    else:
        X_filt = scipy.sparse.csr_matrix(X[:, gene_mask])

    genes_filt = gene_names[gene_mask]

    # inferCNV expects genes × cells (transposed from AnnData convention)
    X_t = X_filt.T.tocsr()

    # --- Write files -----------------------------------------------------
    # matrix.mtx.gz
    mtx_path = out_dir / "matrix.mtx.gz"
    with gzip.open(str(mtx_path), "wb") as fh:
        buf = io.BytesIO()
        scipy.io.mmwrite(buf, X_t.astype(np.float32))
        fh.write(buf.getvalue())

    # barcodes.tsv
    with open(out_dir / "barcodes.tsv", "w") as fh:
        fh.write("\n".join(barcodes) + "\n")

    # features.tsv
    with open(out_dir / "features.tsv", "w") as fh:
        fh.write("\n".join(genes_filt) + "\n")

    # annotations.tsv  (barcode  group — two tab-separated columns, no header)
    ann_df = pd.DataFrame({"barcode": barcodes, "group": groups})
    ann_df.to_csv(str(out_dir / "annotations.tsv"), sep="\t", index=False, header=False)

    # celltype.tsv for diagnostics
    pd.DataFrame({"barcode": barcodes, "cell_type": celltypes, "group": groups}).to_csv(
        str(out_dir / "celltype_debug.tsv"), sep="\t", index=False
    )

    done_flag.touch()
    return {"sample_id": sample_id, "status": "completed", "reason": ""}


# ── Load annotation tables ───────────────────────────────────────────────────
def load_gse238264_celltype_map() -> dict:
    ann_path = ANN_DIR / "gse238264_celltype_annotations.csv"
    if not ann_path.exists():
        raise FileNotFoundError(f"Missing: {ann_path}")
    df = pd.read_csv(ann_path, index_col=0)
    # Index is like HCC1R_AAACAAGTATCTCCCA-1-0
    # We need barcode (after removing sample prefix) → celltype
    # but also keep full key for lookup
    celltype_col = "celltype" if "celltype" in df.columns else "cell_type"
    result = {}
    for full_key, row in df.iterrows():
        result[str(full_key)] = str(row[celltype_col])
        # also store barcode-only key (strip sample prefix + trailing -N index)
        parts = str(full_key).split("_", 1)
        if len(parts) == 2:
            bc_raw = parts[1]
            # barcode might end in -0, -1, etc from concat – strip last -N
            bc_clean = bc_raw.rsplit("-", 1)[0] if bc_raw.count("-") > 1 else bc_raw
            result[bc_clean] = str(row[celltype_col])
    return result


def load_gse281759_celltype_map() -> dict:
    ann_path = ANN_DIR / "gse281759_celltype_annotations.csv"
    if not ann_path.exists():
        return {}
    df = pd.read_csv(ann_path)
    if "dominant_cell_type" not in df.columns:
        return {}
    result = {}
    for _, row in df.iterrows():
        bc = str(row["barcodes"])
        ct = str(row["dominant_cell_type"])
        result[bc] = ct
    return result


def discover_gse281759_samples() -> list:
    symlink_dir = DATA / "spatial_data" / "GSE281759_HCC_MultiRegion"
    samples = []
    if not symlink_dir.exists():
        return samples
    for p in sorted(symlink_dir.iterdir()):
        if p.suffix == ".h5ad":
            sid = p.stem
            real = _gse281759_h5ad(sid)
            if real.exists():
                samples.append(sid)
    return samples


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    gene_order_df = build_gene_order(GTF_PATH, OUT_ROOT / "gene_order.tsv")

    run_log = []

    # ── GSE238264 ────────────────────────────────────────────────────────
    print("\n=== GSE238264 (HCC IO Visium) ===")
    gse238264_map = load_gse238264_celltype_map()
    print(f"  Loaded {len(gse238264_map)} barcode annotations")

    for sample_id, response in GSE238264_SAMPLES.items():
        h5ad = _gse238264_h5ad(sample_id)
        print(f"  [{sample_id}] {response}")
        if not h5ad.exists():
            run_log.append({"sample_id": sample_id, "status": "missing", "reason": str(h5ad)})
            print(f"    MISSING: {h5ad}")
            continue
        result = write_sample_inputs(
            sample_id=sample_id,
            h5ad_path=h5ad,
            celltype_map=gse238264_map,
            ref_types=GSE238264_REF_TYPES,
            out_dir=OUT_ROOT / sample_id,
            gene_order_df=gene_order_df,
        )
        result["dataset"] = "GSE238264"
        result["response"] = response
        run_log.append(result)
        print(f"    → {result['status']}: {result['reason']}")

    # ── GSE281759 ────────────────────────────────────────────────────────
    print("\n=== GSE281759 (HCC MultiRegion) ===")
    gse281759_map = load_gse281759_celltype_map()
    print(f"  Loaded {len(gse281759_map)} barcode annotations")

    gse281759_samples = discover_gse281759_samples()
    print(f"  Found {len(gse281759_samples)} samples")

    for sample_id in gse281759_samples:
        h5ad = _gse281759_h5ad(sample_id)
        print(f"  [{sample_id}]")
        result = write_sample_inputs(
            sample_id=sample_id,
            h5ad_path=h5ad,
            celltype_map=gse281759_map,
            ref_types=GSE281759_REF_TYPES,
            out_dir=OUT_ROOT / sample_id,
            gene_order_df=gene_order_df,
        )
        result["dataset"] = "GSE281759"
        result["response"] = ""
        run_log.append(result)
        print(f"    → {result['status']}: {result['reason']}")

    # ── Summary ──────────────────────────────────────────────────────────
    log_df = pd.DataFrame(run_log)
    log_path = OUT_ROOT / "prepare_log.csv"
    log_df.to_csv(log_path, index=False)

    print(f"\n=== Preparation summary ===")
    print(log_df["status"].value_counts().to_string())
    print(f"Log → {log_path}")

    n_ready = (log_df["status"].isin(["completed", "already_done"])).sum()
    if n_ready == 0:
        print("ERROR: No samples are ready for inferCNV.", file=sys.stderr)
        sys.exit(1)
    print(f"\n{n_ready} samples ready for inferCNV.")


if __name__ == "__main__":
    main()
