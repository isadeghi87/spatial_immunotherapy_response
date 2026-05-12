#!/usr/bin/env python3
"""
05h_geneformer_patient_features.py - Extract patient-level Geneformer features.

Question
--------
Does a frozen Geneformer representation add a transfer-capable patient-level
signal beyond the existing L1-L5 feature stack?

Method
------
1. Stream labelled cells from the integrated atlas counts layer.
2. Convert each cell to Geneformer V2 rank-encoded input using the local
   Gencode symbol-to-Ensembl map and the Geneformer token / median dictionaries.
3. Run a frozen Geneformer model and collect CLS embeddings.
4. Aggregate embeddings to the patient level as mean CLS vectors.

Outputs
-------
  out/05_multilayer/tables/geneformer_patient_features_<model>.csv
  out/05_multilayer/tables/geneformer_patient_qc_<model>.csv
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from collections import defaultdict
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))

import h5py
import numpy as np
import pandas as pd
import torch

from geneformer_utils import (
    DEFAULT_DATASETS,
    DEFAULT_MODEL_VARIANT,
    build_feature_map,
    build_symbol_to_ensembl_map,
    default_output_tag,
    get_model_input_size,
    infer_cls_embeddings,
    load_model,
    load_v2_gene_dictionaries,
    load_var_names,
    resolve_device,
    sanitize_model_name,
    select_cells,
    seed_everything,
    tokenize_count_row,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract patient-level Geneformer features")
    parser.add_argument(
        "--atlas",
        type=Path,
        default=Path("data/processed/integrated_io_atlas.h5ad"),
        help="Integrated atlas h5ad file with the counts layer.",
    )
    parser.add_argument(
        "--asset-dir",
        type=Path,
        default=Path("data/references/geneformer"),
        help="Directory created by 05g_download_geneformer_assets.py.",
    )
    parser.add_argument(
        "--model-variant",
        type=str,
        default=DEFAULT_MODEL_VARIANT,
        help="Geneformer model subfolder to use.",
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        default=Path("data/references/gencode.v44.basic.annotation.gtf.gz"),
        help="Gencode GTF used to map gene symbols to Ensembl IDs.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=8,
        help="Number of tokenized cells per forward pass.",
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=list(DEFAULT_DATASETS),
        help="Datasets to include when building patient-level Geneformer features.",
    )
    parser.add_argument(
        "--max-cells-per-patient",
        type=int,
        default=None,
        help="Optional cap for exploratory runs and smoke tests.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed used when max-cells-per-patient is set.",
    )
    parser.add_argument(
        "--force-cpu",
        action="store_true",
        help="Force CPU inference even if CUDA is available.",
    )
    parser.add_argument(
        "--output-tag",
        type=str,
        default=None,
        help="Optional tag appended to output filenames. Defaults to the model tag, or an automatic subsample tag.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    seed_everything(args.seed)

    out_dir = PROJECT / "out" / "05_multilayer" / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    asset_dir = (PROJECT / args.asset_dir).resolve() if not args.asset_dir.is_absolute() else args.asset_dir.resolve()
    atlas_path = (PROJECT / args.atlas).resolve() if not args.atlas.is_absolute() else args.atlas.resolve()
    gtf_path = (PROJECT / args.gtf).resolve() if not args.gtf.is_absolute() else args.gtf.resolve()
    model_dir = asset_dir / args.model_variant
    model_tag = sanitize_model_name(args.model_variant)
    if args.output_tag is not None:
        output_tag = args.output_tag
    else:
        output_tag = default_output_tag(args.model_variant, args.max_cells_per_patient)
    feature_path = out_dir / f"geneformer_patient_features_{output_tag}.csv"
    qc_path = out_dir / f"geneformer_patient_qc_{output_tag}.csv"

    if not atlas_path.exists():
        raise FileNotFoundError(f"Atlas file not found: {atlas_path}")
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")
    if not model_dir.exists():
        raise FileNotFoundError(f"Model directory not found: {model_dir}")

    print("=" * 70)
    print("[05h] Geneformer patient feature extraction")
    print("=" * 70)
    print(f"Atlas          : {atlas_path}")
    print(f"Assets         : {asset_dir}")
    print(f"Model variant  : {args.model_variant}")
    print(f"Datasets       : {', '.join(args.datasets)}")
    print(f"Batch size     : {args.batch_size}")
    print(f"Output tag     : {output_tag}")

    token_dict, median_dict = load_v2_gene_dictionaries(asset_dir)
    model_input_size = get_model_input_size(args.model_variant)
    cls_token_id = int(token_dict["<cls>"])
    eos_token_id = int(token_dict["<eos>"])
    pad_token_id = int(token_dict["<pad>"])

    cache_path = asset_dir / "geneformer" / "symbol_to_ensembl_gencode_v44.pkl"
    symbol_to_ensembl = build_symbol_to_ensembl_map(gtf_path, cache_path=cache_path)
    var_names = load_var_names(atlas_path)
    feature_map = build_feature_map(var_names, symbol_to_ensembl, token_dict, median_dict)
    print(
        f"Mapped genes   : {feature_map.n_mapped_genes}/{feature_map.n_total_genes} "
        f"({feature_map.n_mapped_genes / feature_map.n_total_genes:.1%})"
    )

    selected = select_cells(
        atlas_path,
        datasets=tuple(args.datasets),
        max_cells_per_patient=args.max_cells_per_patient,
        seed=args.seed,
    )
    if len(selected.cell_indices) == 0:
        raise RuntimeError("No cells matched the requested Geneformer extraction filters.")
    print(f"Selected cells : {len(selected.cell_indices):,}")
    print(f"Patients       : {len(np.unique(selected.patient_ids))}")

    device = resolve_device(force_cpu=args.force_cpu)
    print(f"Device         : {device}")
    model = load_model(model_dir, device)
    emb_dim = int(model.config.hidden_size)
    print(f"Embedding dim  : {emb_dim}")

    patient_sums: dict[str, np.ndarray] = defaultdict(lambda: np.zeros(emb_dim, dtype=np.float64))
    patient_cell_counts: dict[str, int] = defaultdict(int)
    patient_response: dict[str, str] = {}
    patient_dataset: dict[str, str] = {}

    tokenized_cells = 0
    skipped_cells = 0
    batch_sequences: list[np.ndarray] = []
    batch_patients: list[str] = []
    batch_responses: list[str] = []
    batch_datasets: list[str] = []

    with h5py.File(atlas_path, "r") as handle:
        counts = handle["layers"]["counts"]
        indptr = counts["indptr"]
        indices = counts["indices"]
        data = counts["data"]

        for counter, (cell_index, patient_id, response, dataset) in enumerate(
            zip(selected.cell_indices, selected.patient_ids, selected.responses, selected.datasets),
            start=1,
        ):
            start = int(indptr[cell_index])
            end = int(indptr[cell_index + 1])
            row_feature_indices = indices[start:end]
            row_counts = data[start:end]
            sequence = tokenize_count_row(
                row_counts=row_counts,
                row_feature_indices=row_feature_indices,
                feature_map=feature_map,
                cls_token_id=cls_token_id,
                eos_token_id=eos_token_id,
                model_input_size=model_input_size,
            )
            if sequence is None:
                skipped_cells += 1
                continue

            batch_sequences.append(sequence)
            batch_patients.append(str(patient_id))
            batch_responses.append(str(response))
            batch_datasets.append(str(dataset))

            if len(batch_sequences) == args.batch_size:
                embeddings = infer_cls_embeddings(model, batch_sequences, pad_token_id, device)
                for embedding, batch_patient, batch_response, batch_dataset in zip(
                    embeddings, batch_patients, batch_responses, batch_datasets
                ):
                    patient_sums[batch_patient] += embedding
                    patient_cell_counts[batch_patient] += 1
                    patient_response[batch_patient] = batch_response
                    patient_dataset[batch_patient] = batch_dataset
                    tokenized_cells += 1
                batch_sequences.clear()
                batch_patients.clear()
                batch_responses.clear()
                batch_datasets.clear()

            if counter % 2000 == 0:
                print(
                    f"Processed cells : {counter:,} / {len(selected.cell_indices):,} "
                    f"(tokenized={tokenized_cells:,}, skipped={skipped_cells:,})",
                    flush=True,
                )

    if batch_sequences:
        embeddings = infer_cls_embeddings(model, batch_sequences, pad_token_id, device)
        for embedding, batch_patient, batch_response, batch_dataset in zip(
            embeddings, batch_patients, batch_responses, batch_datasets
        ):
            patient_sums[batch_patient] += embedding
            patient_cell_counts[batch_patient] += 1
            patient_response[batch_patient] = batch_response
            patient_dataset[batch_patient] = batch_dataset
            tokenized_cells += 1

    rows = []
    for patient_id in sorted(patient_sums):
        mean_embedding = patient_sums[patient_id] / patient_cell_counts[patient_id]
        row = {
            "patient_id": patient_id,
            "response": patient_response[patient_id],
            "dataset": patient_dataset[patient_id],
            "n_cells_geneformer": patient_cell_counts[patient_id],
            "geneformer_model": args.model_variant,
        }
        row.update({f"gf_cls_{index:04d}": float(value) for index, value in enumerate(mean_embedding)})
        rows.append(row)

    feature_df = pd.DataFrame(rows)
    if feature_df.empty:
        raise RuntimeError("Geneformer extraction produced no patient-level features.")
    feature_df.to_csv(feature_path, index=False)

    patient_counts = np.array(list(patient_cell_counts.values()), dtype=float)
    run_timestamp = datetime.now(timezone.utc).isoformat()

    qc_df = pd.DataFrame(
        [
            {
                "run_timestamp_utc": run_timestamp,
                "model_variant": args.model_variant,
                "model_dir": str(model_dir),
                "atlas": str(atlas_path),
                "asset_dir": str(asset_dir),
                "gtf": str(gtf_path),
                "device": str(device),
                "seed": args.seed,
                "batch_size": args.batch_size,
                "datasets": "|".join(args.datasets),
                "n_total_genes": feature_map.n_total_genes,
                "n_mapped_genes": feature_map.n_mapped_genes,
                "mapped_gene_fraction": feature_map.n_mapped_genes / feature_map.n_total_genes,
                "n_selected_cells": len(selected.cell_indices),
                "n_tokenized_cells": tokenized_cells,
                "n_skipped_cells": skipped_cells,
                "n_patients": len(feature_df),
                "embedding_dim": emb_dim,
                "max_cells_per_patient": args.max_cells_per_patient,
                "mean_cells_per_patient": float(patient_counts.mean()),
                "median_cells_per_patient": float(np.median(patient_counts)),
                "min_cells_per_patient": int(patient_counts.min()),
                "max_cells_per_patient_observed": int(patient_counts.max()),
                "output_tag": output_tag,
            }
        ]
    )
    qc_df.to_csv(qc_path, index=False)

    summary = {
        "feature_file": str(feature_path),
        "qc_file": str(qc_path),
        "output_tag": output_tag,
        "n_patients": int(len(feature_df)),
        "n_tokenized_cells": int(tokenized_cells),
        "n_skipped_cells": int(skipped_cells),
    }
    print(json.dumps(summary, indent=2))
    print("[DONE]")


if __name__ == "__main__":
    main()