#!/usr/bin/env python3
"""
05j_geneformer_finetune.py - Fine-tune Geneformer V2-104M_CLcancer for IO response.

Question
--------
Can cell-level fine-tuning of the top Geneformer transformer layers produce a
patient-level IO-response signal that transfers to an independent cohort, beyond
what frozen CLS embeddings provide?

Method
------
1. Load CLcancer pre-trained weights (BertModel encoder, 12 layers, hidden 768).
2. Freeze layers 0-8; unfreeze layers 9-11 + pooler.
3. Add a 2-class linear head on the [CLS] token.
4. LOOCV on discovery patients (GSE206325, n=24):  each patient is held out in
   turn; all cells from the remaining patients train the model; the held-out
   patient's cells are scored and averaged to a patient-level logit.
5. Retrain on all discovery patients (final model). Evaluate on Cappuyns2023.
6. Save: fine-tuned checkpoint, LOOCV predictions, Cappuyns predictions, and
   patient-level CLS embeddings from the fine-tuned encoder (for downstream
   L6 benchmark integration).

Outputs
-------
  out/05_multilayer/models/geneformer_finetuned_final/
  out/05_multilayer/tables/geneformer_finetune_loocv_results.csv
  out/05_multilayer/tables/geneformer_finetune_cappuyns_results.csv
  out/05_multilayer/tables/geneformer_patient_features_finetuned.csv
  out/05_multilayer/tables/geneformer_finetune_qc.csv
  out/05_multilayer/figures/05j_finetune_loocv_summary.png

QC
--
- Class balance checked; loss weighted for R/NR imbalance.
- Gradient norms monitored per epoch.
- LOOCV AUC reported with 95% CI via bootstrap.
- Memory usage logged at start.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "scripts" / "05_multilayer"))

import h5py
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.amp import GradScaler, autocast
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, average_precision_score

from geneformer_utils import (
    DEFAULT_DATASETS,
    DEFAULT_MODEL_VARIANT,
    build_feature_map,
    build_symbol_to_ensembl_map,
    default_output_tag,
    get_model_input_size,
    load_v2_gene_dictionaries,
    load_var_names,
    pad_token_batch,
    resolve_device,
    sanitize_model_name,
    seed_everything,
    select_cells,
    tokenize_count_row,
)


# ---------------------------------------------------------------------------
# Classifier wrapper
# ---------------------------------------------------------------------------

class GeneformerClassifier(nn.Module):
    """BERT encoder + linear classification head on [CLS] token."""

    def __init__(self, encoder: nn.Module, hidden_size: int = 768, n_classes: int = 2, dropout: float = 0.1) -> None:
        super().__init__()
        self.encoder = encoder
        self.dropout = nn.Dropout(dropout)
        self.classifier = nn.Linear(hidden_size, n_classes)

    def forward(self, input_ids: torch.Tensor, attention_mask: torch.Tensor) -> torch.Tensor:
        outputs = self.encoder(input_ids=input_ids, attention_mask=attention_mask)
        cls = outputs.last_hidden_state[:, 0, :]
        return self.classifier(self.dropout(cls))

    def cls_embeddings(self, input_ids: torch.Tensor, attention_mask: torch.Tensor) -> torch.Tensor:
        with torch.no_grad():
            outputs = self.encoder(input_ids=input_ids, attention_mask=attention_mask)
            return outputs.last_hidden_state[:, 0, :]


def freeze_layers(encoder: nn.Module, n_freeze: int = 9) -> None:
    """Freeze the first n_freeze transformer layers (0-indexed)."""
    # Freeze embeddings
    for param in encoder.embeddings.parameters():
        param.requires_grad = False
    # Freeze layers 0..(n_freeze-1)
    for i in range(n_freeze):
        for param in encoder.encoder.layer[i].parameters():
            param.requires_grad = False
    # Layers n_freeze.. remain trainable (set explicitly for clarity)
    for i in range(n_freeze, len(encoder.encoder.layer)):
        for param in encoder.encoder.layer[i].parameters():
            param.requires_grad = True
    # Pooler trainable
    if hasattr(encoder, "pooler") and encoder.pooler is not None:
        for param in encoder.pooler.parameters():
            param.requires_grad = True


def count_trainable(model: nn.Module) -> tuple[int, int]:
    total = sum(p.numel() for p in model.parameters())
    trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    return trainable, total


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------

def build_cell_dataset(
    atlas_path: Path,
    patient_ids_filter: list[str],
    selected_cells,  # SelectedCells
    feature_map,
    cls_token_id: int,
    eos_token_id: int,
    model_input_size: int,
    label_map: dict[str, int],
    max_cells_per_patient: int | None,
    rng: np.random.Generator,
) -> tuple[list[np.ndarray], list[int], list[str]]:
    """Tokenize cells for a given set of patients; return (sequences, labels, patient_ids)."""
    patient_set = set(patient_ids_filter)
    # Build per-patient lists first for optional subsampling
    per_patient: dict[str, list[tuple[np.ndarray, int]]] = defaultdict(list)

    with h5py.File(atlas_path, "r") as handle:
        counts = handle["layers"]["counts"]
        indptr = counts["indptr"][:]
        indices = counts["indices"]
        data_arr = counts["data"]

        for cell_index, patient_id, response, *_ in zip(
            selected_cells.cell_indices,
            selected_cells.patient_ids,
            selected_cells.responses,
            *([selected_cells.datasets] if hasattr(selected_cells, "datasets") else [[]]),
        ):
            pid = str(patient_id)
            if pid not in patient_set:
                continue
            label = label_map.get(str(response))
            if label is None:
                continue
            start = int(indptr[cell_index])
            end = int(indptr[cell_index + 1])
            seq = tokenize_count_row(
                row_counts=data_arr[start:end],
                row_feature_indices=indices[start:end],
                feature_map=feature_map,
                cls_token_id=cls_token_id,
                eos_token_id=eos_token_id,
                model_input_size=model_input_size,
            )
            if seq is None:
                continue
            per_patient[pid].append((seq, label))

    sequences: list[np.ndarray] = []
    labels: list[int] = []
    patient_out: list[str] = []

    for pid, items in per_patient.items():
        if max_cells_per_patient is not None and len(items) > max_cells_per_patient:
            chosen_idx = rng.choice(len(items), size=max_cells_per_patient, replace=False)
            items = [items[i] for i in sorted(chosen_idx)]
        for seq, lbl in items:
            sequences.append(seq)
            labels.append(lbl)
            patient_out.append(pid)

    return sequences, labels, patient_out


# ---------------------------------------------------------------------------
# Training helpers
# ---------------------------------------------------------------------------

def train_epoch(
    model: GeneformerClassifier,
    sequences: list[np.ndarray],
    labels: list[int],
    pad_token_id: int,
    device: torch.device,
    optimizer: torch.optim.Optimizer,
    scaler: GradScaler,
    batch_size: int,
    class_weights: torch.Tensor,
    rng: np.random.Generator,
) -> tuple[float, float]:
    model.train()
    criterion = nn.CrossEntropyLoss(weight=class_weights.to(device))
    total_loss = 0.0
    total_norm = 0.0
    n_batches = 0

    perm = rng.permutation(len(sequences))
    for start in range(0, len(sequences), batch_size):
        batch_idx = perm[start : start + batch_size]
        batch_seqs = [sequences[i] for i in batch_idx]
        batch_labels = torch.tensor([labels[i] for i in batch_idx], dtype=torch.long, device=device)
        input_ids, attention_mask = pad_token_batch(batch_seqs, pad_token_id, device)

        optimizer.zero_grad(set_to_none=True)
        with autocast("cuda", enabled=device.type == "cuda"):
            logits = model(input_ids, attention_mask)
            loss = criterion(logits, batch_labels)

        scaler.scale(loss).backward()
        scaler.unscale_(optimizer)
        norm = nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        scaler.step(optimizer)
        scaler.update()

        total_loss += loss.item()
        total_norm += float(norm)
        n_batches += 1

    return total_loss / max(n_batches, 1), total_norm / max(n_batches, 1)


@torch.no_grad()
def predict_patient_scores(
    model: GeneformerClassifier,
    sequences: list[np.ndarray],
    patient_ids: list[str],
    pad_token_id: int,
    device: torch.device,
    batch_size: int,
) -> dict[str, float]:
    """Return per-patient mean Responder class probability."""
    model.eval()
    patient_logit_sums: dict[str, float] = defaultdict(float)
    patient_counts: dict[str, int] = defaultdict(int)

    for start in range(0, len(sequences), batch_size):
        batch_seqs = sequences[start : start + batch_size]
        batch_pids = patient_ids[start : start + batch_size]
        input_ids, attention_mask = pad_token_batch(batch_seqs, pad_token_id, device)
        with autocast("cuda", enabled=device.type == "cuda"):
            logits = model(input_ids, attention_mask)
        probs = torch.softmax(logits.float(), dim=-1)[:, 1].cpu().numpy()
        for pid, prob in zip(batch_pids, probs):
            patient_logit_sums[pid] += float(prob)
            patient_counts[pid] += 1

    return {pid: patient_logit_sums[pid] / patient_counts[pid] for pid in patient_logit_sums}


@torch.no_grad()
def extract_patient_embeddings(
    model: GeneformerClassifier,
    sequences: list[np.ndarray],
    patient_ids: list[str],
    pad_token_id: int,
    device: torch.device,
    batch_size: int,
    emb_dim: int,
) -> dict[str, np.ndarray]:
    """Return per-patient mean CLS embeddings from fine-tuned encoder."""
    model.eval()
    patient_sums: dict[str, np.ndarray] = defaultdict(lambda: np.zeros(emb_dim, dtype=np.float64))
    patient_counts: dict[str, int] = defaultdict(int)

    for start in range(0, len(sequences), batch_size):
        batch_seqs = sequences[start : start + batch_size]
        batch_pids = patient_ids[start : start + batch_size]
        input_ids, attention_mask = pad_token_batch(batch_seqs, pad_token_id, device)
        embs = model.cls_embeddings(input_ids, attention_mask).cpu().numpy()
        for pid, emb in zip(batch_pids, embs):
            patient_sums[pid] += emb.astype(np.float64)
            patient_counts[pid] += 1

    return {pid: patient_sums[pid] / patient_counts[pid] for pid in patient_sums}


# ---------------------------------------------------------------------------
# Bootstrap CI helper
# ---------------------------------------------------------------------------

def bootstrap_auc_ci(y_true: np.ndarray, y_score: np.ndarray, n: int = 1000, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    aucs = []
    for _ in range(n):
        idx = rng.integers(0, len(y_true), size=len(y_true))
        if len(np.unique(y_true[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y_true[idx], y_score[idx]))
    aucs = np.array(aucs)
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


# ---------------------------------------------------------------------------
# Model factory (re-instantiates fresh model for each LOOCV fold)
# ---------------------------------------------------------------------------

def build_fresh_model(model_dir: Path, device: torch.device, n_freeze: int, hidden_size: int, enable_gradient_checkpointing: bool) -> GeneformerClassifier:
    from transformers import AutoModel

    encoder = AutoModel.from_pretrained(
        str(model_dir),
        trust_remote_code=False,
        attn_implementation="sdpa",
    )
    if enable_gradient_checkpointing:
        encoder.gradient_checkpointing_enable()
    freeze_layers(encoder, n_freeze=n_freeze)
    clf = GeneformerClassifier(encoder, hidden_size=hidden_size)
    clf.to(device)
    return clf


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fine-tune Geneformer V2 for IO response classification")
    parser.add_argument("--atlas", type=Path, default=Path("data/processed/integrated_io_atlas.h5ad"))
    parser.add_argument("--asset-dir", type=Path, default=Path("data/references/geneformer"))
    parser.add_argument("--model-variant", type=str, default=DEFAULT_MODEL_VARIANT)
    parser.add_argument("--gtf", type=Path, default=Path("data/references/gencode.v44.basic.annotation.gtf.gz"))
    parser.add_argument("--max-cells-per-patient", type=int, default=200)
    parser.add_argument("--batch-size", type=int, default=4, help="Cells per forward/backward pass (keep low for 8GB GPU)")
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--n-freeze", type=int, default=9, help="Number of BERT layers to freeze (0-indexed; default: 9/12)")
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--force-cpu", action="store_true")
    parser.add_argument("--no-gradient-checkpointing", action="store_true")
    parser.add_argument("--smoke-test", action="store_true", help="Quick validation: 4 patients, 30 cells, 1 epoch")
    parser.add_argument("--max-patients", type=int, default=None, help="Cap on number of discovery patients (for testing)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    seed_everything(args.seed)

    if args.smoke_test:
        args.max_cells_per_patient = min(args.max_cells_per_patient or 9999, 30)
        args.max_patients = args.max_patients or 4
        args.epochs = 1
        print("[SMOKE TEST] max_cells_per_patient=30, max_patients=4, epochs=1")

    device = resolve_device(force_cpu=args.force_cpu)

    PROJECT_ROOT = Path(__file__).resolve().parents[2]
    asset_dir = (PROJECT_ROOT / args.asset_dir).resolve() if not args.asset_dir.is_absolute() else args.asset_dir.resolve()
    atlas_path = (PROJECT_ROOT / args.atlas).resolve() if not args.atlas.is_absolute() else args.atlas.resolve()
    gtf_path = (PROJECT_ROOT / args.gtf).resolve() if not args.gtf.is_absolute() else args.gtf.resolve()
    model_dir = asset_dir / args.model_variant

    out_tables = PROJECT_ROOT / "out" / "05_multilayer" / "tables"
    out_models = PROJECT_ROOT / "out" / "05_multilayer" / "models"
    out_figs = PROJECT_ROOT / "out" / "05_multilayer" / "figures"
    for d in (out_tables, out_models, out_figs):
        d.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("[05j] Geneformer fine-tuning for IO response")
    print("=" * 70)
    print(f"Device         : {device}")
    if device.type == "cuda":
        props = torch.cuda.get_device_properties(device)
        print(f"GPU            : {props.name} ({props.total_memory / 1e9:.1f} GB)")
    print(f"Model          : {args.model_variant}")
    print(f"Freeze layers  : 0-{args.n_freeze - 1} (fine-tune {args.n_freeze}-11)")
    print(f"Epochs         : {args.epochs}  LR: {args.lr}  Batch: {args.batch_size}")

    # --- Load dictionaries and build feature map ---
    token_dict, median_dict = load_v2_gene_dictionaries(asset_dir)
    model_input_size = get_model_input_size(args.model_variant)
    cls_token_id = int(token_dict["<cls>"])
    eos_token_id = int(token_dict["<eos>"])
    pad_token_id = int(token_dict["<pad>"])

    cache_path = asset_dir / "geneformer" / "symbol_to_ensembl_gencode_v44.pkl"
    symbol_to_ensembl = build_symbol_to_ensembl_map(gtf_path, cache_path=cache_path)
    var_names = load_var_names(atlas_path)
    feature_map = build_feature_map(var_names, symbol_to_ensembl, token_dict, median_dict)
    print(f"Mapped genes   : {feature_map.n_mapped_genes}/{feature_map.n_total_genes}")

    # --- Load architecture config for hidden_size ---
    config_path = model_dir / "config.json"
    with config_path.open() as handle:
        model_config = json.load(handle)
    hidden_size = int(model_config["hidden_size"])
    n_total_layers = int(model_config["num_hidden_layers"])
    print(f"Architecture   : hidden={hidden_size}, layers={n_total_layers}")

    # --- Select all labelled cells (discovery + Cappuyns) ---
    all_selected = select_cells(atlas_path, datasets=tuple(DEFAULT_DATASETS))
    disc_patients = sorted(set(
        str(p) for p, d in zip(all_selected.patient_ids, all_selected.datasets)
        if str(d) != "Cappuyns2023"
    ))
    capp_patients = sorted(set(
        str(p) for p, d in zip(all_selected.patient_ids, all_selected.datasets)
        if str(d) == "Cappuyns2023"
    ))

    if args.max_patients is not None:
        disc_patients = disc_patients[: args.max_patients]

    # Build metadata map
    pid_response: dict[str, str] = {}
    pid_dataset: dict[str, str] = {}
    for pid, resp, dset in zip(all_selected.patient_ids, all_selected.responses, all_selected.datasets):
        pid_response[str(pid)] = str(resp)
        pid_dataset[str(pid)] = str(dset)

    label_map = {"Responder": 1, "Non-Responder": 0}
    disc_labels_by_patient = {p: label_map[pid_response[p]] for p in disc_patients if pid_response.get(p) in label_map}
    disc_patients = [p for p in disc_patients if p in disc_labels_by_patient]

    n_R = sum(disc_labels_by_patient[p] == 1 for p in disc_patients)
    n_NR = sum(disc_labels_by_patient[p] == 0 for p in disc_patients)
    print(f"Discovery      : n={len(disc_patients)} (R={n_R}, NR={n_NR})")
    print(f"Cappuyns       : n={len(capp_patients)}")

    # Class weights inversely proportional to frequency
    n_total_cls = n_R + n_NR
    class_weights = torch.tensor(
        [n_total_cls / (2.0 * max(n_NR, 1)), n_total_cls / (2.0 * max(n_R, 1))],
        dtype=torch.float32,
    )
    print(f"Class weights  : NR={class_weights[0]:.3f}, R={class_weights[1]:.3f}")

    rng = np.random.default_rng(args.seed)
    enable_gc = not args.no_gradient_checkpointing

    # --- LOOCV on discovery ---
    print("\n--- Leave-one-out cross-validation on discovery ---")
    loocv_rows = []
    t0_loocv = time.time()

    for fold_i, held_out in enumerate(disc_patients):
        train_patients = [p for p in disc_patients if p != held_out]
        print(f"  Fold {fold_i + 1}/{len(disc_patients)}: held-out={held_out} "
              f"({pid_response[held_out]})", flush=True)

        # Build training data
        fold_rng = np.random.default_rng(args.seed + fold_i)
        train_seqs, train_labels, train_pids = build_cell_dataset(
            atlas_path=atlas_path,
            patient_ids_filter=train_patients,
            selected_cells=all_selected,
            feature_map=feature_map,
            cls_token_id=cls_token_id,
            eos_token_id=eos_token_id,
            model_input_size=model_input_size,
            label_map=label_map,
            max_cells_per_patient=args.max_cells_per_patient,
            rng=fold_rng,
        )
        test_seqs, _, test_pids = build_cell_dataset(
            atlas_path=atlas_path,
            patient_ids_filter=[held_out],
            selected_cells=all_selected,
            feature_map=feature_map,
            cls_token_id=cls_token_id,
            eos_token_id=eos_token_id,
            model_input_size=model_input_size,
            label_map=label_map,
            max_cells_per_patient=args.max_cells_per_patient,
            rng=fold_rng,
        )

        if len(train_seqs) == 0 or len(test_seqs) == 0:
            print(f"    [WARN] Insufficient data for fold {fold_i + 1}, skipping")
            continue

        # Build fresh model for each fold
        clf = build_fresh_model(model_dir, device, args.n_freeze, hidden_size, enable_gc)
        trainable, total = count_trainable(clf)
        if fold_i == 0:
            print(f"  Trainable params: {trainable:,} / {total:,} ({trainable / total:.1%})")

        optimizer = torch.optim.AdamW(
            [p for p in clf.parameters() if p.requires_grad],
            lr=args.lr, weight_decay=0.01
        )
        scaler = GradScaler("cuda", enabled=device.type == "cuda")

        for epoch in range(args.epochs):
            loss, grad_norm = train_epoch(
                clf, train_seqs, train_labels, pad_token_id, device,
                optimizer, scaler, args.batch_size, class_weights, fold_rng
            )
            if fold_i == 0 or (fold_i == len(disc_patients) - 1):
                print(f"    Epoch {epoch + 1}/{args.epochs}: loss={loss:.4f} grad_norm={grad_norm:.3f}")

        patient_scores = predict_patient_scores(clf, test_seqs, test_pids, pad_token_id, device, args.batch_size)
        pred_score = patient_scores.get(held_out, 0.5)
        true_label = disc_labels_by_patient[held_out]
        print(f"    score={pred_score:.4f}  true={pid_response[held_out]}")
        loocv_rows.append({
            "patient_id": held_out,
            "response": pid_response[held_out],
            "dataset": pid_dataset.get(held_out, ""),
            "pred_score": pred_score,
            "true_label": true_label,
            "n_train_cells": len(train_seqs),
            "n_test_cells": len(test_seqs),
        })

        # Free fold model
        del clf
        if device.type == "cuda":
            torch.cuda.empty_cache()

    loocv_df = pd.DataFrame(loocv_rows)
    loocv_path = out_tables / "geneformer_finetune_loocv_results.csv"
    loocv_df.to_csv(loocv_path, index=False)

    if len(loocv_df) >= 2 and len(loocv_df["true_label"].unique()) == 2:
        disc_auc = roc_auc_score(loocv_df["true_label"], loocv_df["pred_score"])
        ci_lo, ci_hi = bootstrap_auc_ci(loocv_df["true_label"].values, loocv_df["pred_score"].values)
        print(f"\nLOOCV AUC      : {disc_auc:.3f}  95% CI: [{ci_lo:.3f}, {ci_hi:.3f}]")
    else:
        disc_auc = float("nan")
        ci_lo, ci_hi = float("nan"), float("nan")
        print("[WARN] Not enough patients/classes for LOOCV AUC (smoke-test mode?)")

    elapsed_loocv = time.time() - t0_loocv
    print(f"LOOCV time     : {elapsed_loocv / 60:.1f} min")

    # --- Final model: retrain on all discovery ---
    print("\n--- Retraining final model on all discovery patients ---")
    final_rng = np.random.default_rng(args.seed + 999)
    all_disc_seqs, all_disc_labels, all_disc_pids = build_cell_dataset(
        atlas_path=atlas_path,
        patient_ids_filter=disc_patients,
        selected_cells=all_selected,
        feature_map=feature_map,
        cls_token_id=cls_token_id,
        eos_token_id=eos_token_id,
        model_input_size=model_input_size,
        label_map=label_map,
        max_cells_per_patient=args.max_cells_per_patient,
        rng=final_rng,
    )
    print(f"Final training : {len(all_disc_seqs):,} cells from {len(disc_patients)} patients")

    final_clf = build_fresh_model(model_dir, device, args.n_freeze, hidden_size, enable_gc)
    final_optimizer = torch.optim.AdamW(
        [p for p in final_clf.parameters() if p.requires_grad],
        lr=args.lr, weight_decay=0.01
    )
    final_scaler = GradScaler("cuda", enabled=device.type == "cuda")

    for epoch in range(args.epochs):
        loss, grad_norm = train_epoch(
            final_clf, all_disc_seqs, all_disc_labels, pad_token_id, device,
            final_optimizer, final_scaler, args.batch_size, class_weights, final_rng
        )
        print(f"  Epoch {epoch + 1}/{args.epochs}: loss={loss:.4f} grad_norm={grad_norm:.3f}")

    # Save final model
    final_model_path = out_models / "geneformer_finetuned_final"
    final_model_path.mkdir(parents=True, exist_ok=True)
    torch.save(final_clf.state_dict(), final_model_path / "pytorch_model.bin")
    # Copy config for reference
    import shutil
    shutil.copy(model_dir / "config.json", final_model_path / "base_config.json")
    finetune_config = {
        "base_model": str(model_dir),
        "n_freeze": args.n_freeze,
        "lr": args.lr,
        "epochs": args.epochs,
        "max_cells_per_patient": args.max_cells_per_patient,
        "seed": args.seed,
        "hidden_size": hidden_size,
        "n_classes": 2,
        "dropout": args.dropout,
        "trained_on": "GSE206325_discovery",
    }
    with (final_model_path / "finetune_config.json").open("w") as handle:
        json.dump(finetune_config, handle, indent=2)
    print(f"Saved final model: {final_model_path}")

    # --- Evaluate on Cappuyns ---
    if capp_patients and not args.smoke_test:
        print("\n--- Evaluating on Cappuyns2023 ---")
        capp_seqs, _, capp_pids = build_cell_dataset(
            atlas_path=atlas_path,
            patient_ids_filter=capp_patients,
            selected_cells=all_selected,
            feature_map=feature_map,
            cls_token_id=cls_token_id,
            eos_token_id=eos_token_id,
            model_input_size=model_input_size,
            label_map=label_map,
            max_cells_per_patient=args.max_cells_per_patient,
            rng=np.random.default_rng(args.seed),
        )
        capp_scores = predict_patient_scores(final_clf, capp_seqs, capp_pids, pad_token_id, device, args.batch_size)
        capp_rows = []
        for pid in capp_patients:
            if pid not in capp_scores:
                continue
            label = label_map.get(pid_response.get(pid, ""), -1)
            capp_rows.append({
                "patient_id": pid,
                "response": pid_response.get(pid, ""),
                "dataset": "Cappuyns2023",
                "pred_score": capp_scores[pid],
                "true_label": label,
            })
        capp_df = pd.DataFrame(capp_rows)
        capp_path = out_tables / "geneformer_finetune_cappuyns_results.csv"
        capp_df.to_csv(capp_path, index=False)
        if len(capp_df) >= 2 and len(capp_df["true_label"].unique()) == 2:
            capp_auc = roc_auc_score(capp_df["true_label"], capp_df["pred_score"])
            capp_prauc = average_precision_score(capp_df["true_label"], capp_df["pred_score"])
            capp_bal_acc = balanced_accuracy_score(
                capp_df["true_label"], (capp_df["pred_score"] >= 0.5).astype(int)
            )
            print(f"Cappuyns AUC   : {capp_auc:.3f}")
            print(f"Cappuyns PR-AUC: {capp_prauc:.3f}")
            print(f"Cappuyns BalAcc: {capp_bal_acc:.3f}")
        else:
            capp_auc = capp_prauc = capp_bal_acc = float("nan")
    else:
        capp_auc = capp_prauc = capp_bal_acc = float("nan")
        capp_df = pd.DataFrame()
        if args.smoke_test:
            print("[SMOKE TEST] Skipping Cappuyns evaluation")

    # --- Extract fine-tuned patient CLS embeddings (discovery + Cappuyns) ---
    if not args.smoke_test:
        print("\n--- Extracting fine-tuned patient embeddings ---")
        all_seqs, _, all_pids = build_cell_dataset(
            atlas_path=atlas_path,
            patient_ids_filter=disc_patients + capp_patients,
            selected_cells=all_selected,
            feature_map=feature_map,
            cls_token_id=cls_token_id,
            eos_token_id=eos_token_id,
            model_input_size=model_input_size,
            label_map=label_map,
            max_cells_per_patient=args.max_cells_per_patient,
            rng=np.random.default_rng(args.seed),
        )
        emb_dict = extract_patient_embeddings(
            final_clf, all_seqs, all_pids, pad_token_id, device, args.batch_size, hidden_size
        )
        emb_rows = []
        for pid in sorted(emb_dict):
            row = {
                "patient_id": pid,
                "response": pid_response.get(pid, ""),
                "dataset": pid_dataset.get(pid, ""),
                "geneformer_model": f"finetuned_{args.model_variant}",
            }
            row.update({f"gf_cls_{i:04d}": float(v) for i, v in enumerate(emb_dict[pid])})
            emb_rows.append(row)
        emb_df = pd.DataFrame(emb_rows)
        emb_path = out_tables / "geneformer_patient_features_finetuned.csv"
        emb_df.to_csv(emb_path, index=False)
        print(f"Saved embeddings: {emb_path}")
    else:
        emb_path = None

    # --- QC record ---
    qc = {
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "model_variant": args.model_variant,
        "n_freeze": args.n_freeze,
        "epochs": args.epochs,
        "lr": args.lr,
        "batch_size": args.batch_size,
        "max_cells_per_patient": args.max_cells_per_patient,
        "seed": args.seed,
        "device": str(device),
        "n_disc_patients": len(disc_patients),
        "n_disc_R": n_R,
        "n_disc_NR": n_NR,
        "disc_loocv_auc": round(disc_auc, 4) if not np.isnan(disc_auc) else None,
        "disc_loocv_auc_ci_lo": round(ci_lo, 4) if not np.isnan(ci_lo) else None,
        "disc_loocv_auc_ci_hi": round(ci_hi, 4) if not np.isnan(ci_hi) else None,
        "capp_auc": round(capp_auc, 4) if not np.isnan(capp_auc) else None,
        "capp_prauc": round(capp_prauc, 4) if not np.isnan(capp_prauc) else None,
        "capp_balanced_acc": round(capp_bal_acc, 4) if not np.isnan(capp_bal_acc) else None,
        "smoke_test": args.smoke_test,
        "loocv_time_min": round(elapsed_loocv / 60, 2),
    }
    qc_df = pd.DataFrame([qc])
    qc_path = out_tables / "geneformer_finetune_qc.csv"
    qc_df.to_csv(qc_path, index=False)

    # --- Summary figure ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        fig, axes = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)

        # Left: LOOCV score distribution
        if not loocv_df.empty:
            for resp, grp in loocv_df.groupby("response"):
                color = "#CC3333" if resp == "Responder" else "#3366CC"
                axes[0].scatter(
                    [resp] * len(grp), grp["pred_score"], color=color, alpha=0.8, s=50, zorder=3
                )
            axes[0].axhline(0.5, color="gray", lw=1, ls="--")
            axes[0].set_ylabel("Fine-tuned Responder score")
            axes[0].set_title(f"LOOCV discovery (AUC={disc_auc:.3f})", fontsize=10, fontweight="bold")
            axes[0].set_ylim(0, 1)

        # Right: Cappuyns score distribution
        if not capp_df.empty and not args.smoke_test:
            for resp, grp in capp_df.groupby("response"):
                color = "#CC3333" if resp == "Responder" else "#3366CC"
                axes[1].scatter(
                    [resp] * len(grp), grp["pred_score"], color=color, alpha=0.8, s=50, zorder=3
                )
            axes[1].axhline(0.5, color="gray", lw=1, ls="--")
            axes[1].set_ylabel("")
            axes[1].set_title(f"Cappuyns validation (AUC={capp_auc:.3f})", fontsize=10, fontweight="bold")
            axes[1].set_ylim(0, 1)

        legend_handles = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#CC3333", markersize=8, label="Responder"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#3366CC", markersize=8, label="Non-Responder"),
        ]
        fig.legend(handles=legend_handles, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.05))
        fig.suptitle("Geneformer fine-tuned IO response scores", fontsize=12, fontweight="bold")

        fig_path = out_figs / "05j_finetune_loocv_summary.png"
        fig.savefig(fig_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved figure   : {fig_path}")
    except Exception as exc:
        print(f"[WARN] Figure generation failed: {exc}")

    print("\n" + "=" * 70)
    print("[05j] COMPLETE")
    print(f"  LOOCV AUC: {disc_auc:.3f}" if not np.isnan(disc_auc) else "  LOOCV AUC: N/A (smoke test)")
    if not np.isnan(capp_auc):
        print(f"  Cappuyns AUC: {capp_auc:.3f}")
    print(f"  Model saved: {final_model_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
