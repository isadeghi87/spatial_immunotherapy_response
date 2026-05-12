#!/usr/bin/env python3
"""
05k_geneformer_insilico_perturbation.py - In silico gene perturbation for IO response.

Question
--------
Which gene perturbations (knockouts / activations) in specific immune/tumour
populations most strongly shift the predicted IO-response probability of
individual cells, as learned by the fine-tuned Geneformer encoder?

Method
------
1. Load fine-tuned Geneformer classifier from the checkpoint produced by 05j.
   Falls back to frozen CLcancer if no fine-tuned checkpoint exists.
2. For each target cell population (CD8_Dysfunctional-effector,
   CD8_Dysfunctional-progenitor, CD8_Dysfunctional-terminal, CD4_Tfh-like,
   B_Plasma): select cells from the atlas.
3. For each cell: tokenize → get baseline Responder logit.
4. Perturbation loop: for each gene token that appears in at least FREQ_THR of
   cells in this population, zero the token (replace with pad) → recompute
   logit → record Δlogit = perturbed - baseline (positive = gene suppresses R).
5. Aggregate: mean Δlogit, std Δlogit, fraction of cells carrying gene.
6. Rank by |mean Δlogit|; save top-100 per cell type.
7. Generate dot-plot figures for top-20 genes per target.

Outputs
-------
  out/05_multilayer/tables/geneformer_perturbation_{celltype}.csv
  out/05_multilayer/figures/05k_perturbation_top20_{celltype}.png
  out/05_multilayer/tables/geneformer_perturbation_qc.csv

QC
--
- Reports number of cells per target, unique tokens perturbed, runtime.
- Bootstrapped CI on top-5 genes per target.
- Warns if fewer than 20 cells available for a target.
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
from torch.cuda.amp import autocast

from geneformer_utils import (
    DEFAULT_MODEL_VARIANT,
    build_feature_map,
    build_symbol_to_ensembl_map,
    get_model_input_size,
    load_v2_gene_dictionaries,
    load_var_names,
    pad_token_batch,
    resolve_device,
    seed_everything,
    tokenize_count_row,
)

# Target cell populations (cell_type obs field)
DEFAULT_TARGETS = [
    "CD8_Dysfunctional-effector",
    "CD8_Dysfunctional-progenitor",
    "CD8_Dysfunctional-terminal",
    "CD4_Tfh-like",
    "B_Plasma",
]


# ---------------------------------------------------------------------------
# Classifier (same architecture as 05j — must match checkpoint)
# ---------------------------------------------------------------------------

class GeneformerClassifier(nn.Module):
    def __init__(self, encoder: nn.Module, hidden_size: int = 768, n_classes: int = 2, dropout: float = 0.1) -> None:
        super().__init__()
        self.encoder = encoder
        self.dropout = nn.Dropout(dropout)
        self.classifier = nn.Linear(hidden_size, n_classes)

    def forward(self, input_ids: torch.Tensor, attention_mask: torch.Tensor) -> torch.Tensor:
        outputs = self.encoder(input_ids=input_ids, attention_mask=attention_mask)
        cls = outputs.last_hidden_state[:, 0, :]
        return self.classifier(self.dropout(cls))


def load_classifier(model_dir: Path, finetune_checkpoint: Path | None, device: torch.device, hidden_size: int) -> GeneformerClassifier:
    from transformers import AutoModel

    encoder = AutoModel.from_pretrained(
        str(model_dir),
        trust_remote_code=False,
        attn_implementation="sdpa",
    )
    clf = GeneformerClassifier(encoder, hidden_size=hidden_size)

    if finetune_checkpoint is not None and finetune_checkpoint.exists():
        state = torch.load(finetune_checkpoint, map_location=device, weights_only=True)
        clf.load_state_dict(state)
        print(f"Loaded fine-tuned weights: {finetune_checkpoint}")
    else:
        print(f"[WARN] No fine-tuned checkpoint found at {finetune_checkpoint}. Using frozen CLcancer weights.")
        print("       Results reflect frozen model perturbation sensitivities.")

    clf.eval()
    clf.to(device)
    return clf


# ---------------------------------------------------------------------------
# Cell selection by cell type
# ---------------------------------------------------------------------------

def select_cells_by_type(
    atlas_path: Path,
    target_celltype: str,
    datasets: tuple[str, ...],
    max_cells: int | None,
    rng: np.random.Generator,
) -> tuple[np.ndarray, list[str], list[str]]:
    """Return (cell_indices, patient_ids, responses) for a target cell type."""
    with h5py.File(atlas_path, "r") as handle:
        obs = handle["obs"]

        # cell_type categorical
        ct_cats = [x.decode() if isinstance(x, bytes) else x for x in obs["cell_type"]["categories"][:]]
        ct_codes = obs["cell_type"]["codes"][:]
        target_code_set = {i for i, cat in enumerate(ct_cats) if cat == target_celltype}
        if not target_code_set:
            return np.array([], dtype=int), [], []

        # dataset categorical
        ds_cats = [x.decode() if isinstance(x, bytes) else x for x in obs["dataset"]["categories"][:]]
        ds_codes = obs["dataset"]["codes"][:]
        ds_code_set = {i for i, cat in enumerate(ds_cats) if cat in datasets}

        # response categorical
        resp_cats = [x.decode() if isinstance(x, bytes) else x for x in obs["response"]["categories"][:]]
        resp_codes = obs["response"]["codes"][:]
        valid_resp_set = {i for i, cat in enumerate(resp_cats) if cat in ("Responder", "Non-Responder")}

        mask = (
            np.isin(ct_codes, list(target_code_set))
            & np.isin(ds_codes, list(ds_code_set))
            & np.isin(resp_codes, list(valid_resp_set))
        )
        cell_indices = np.flatnonzero(mask)
        if len(cell_indices) == 0:
            return np.array([], dtype=int), [], []

        # Decode patient and response for selected cells
        patient_key = "patient_ID" if "patient_ID" in obs else "Patient_ID"
        patient_cats = [x.decode() if isinstance(x, bytes) else x for x in obs[patient_key]["categories"][:]]
        patient_codes = obs[patient_key]["codes"][cell_indices]
        patient_ids_out = [patient_cats[c] for c in patient_codes]
        resp_out = [resp_cats[resp_codes[i]] for i in cell_indices]

    if max_cells is not None and len(cell_indices) > max_cells:
        chosen = rng.choice(len(cell_indices), size=max_cells, replace=False)
        cell_indices = cell_indices[np.sort(chosen)]
        patient_ids_out = [patient_ids_out[i] for i in np.sort(chosen)]
        resp_out = [resp_out[i] for i in np.sort(chosen)]

    return cell_indices, patient_ids_out, resp_out


# ---------------------------------------------------------------------------
# Tokenize cells from atlas
# ---------------------------------------------------------------------------

def tokenize_cells(
    atlas_path: Path,
    cell_indices: np.ndarray,
    feature_map,
    cls_token_id: int,
    eos_token_id: int,
    model_input_size: int,
) -> tuple[list[np.ndarray], list[int]]:
    """Return (sequences, original_positions) for cells that tokenize successfully."""
    sequences = []
    valid_positions = []
    with h5py.File(atlas_path, "r") as handle:
        counts = handle["layers"]["counts"]
        indptr = counts["indptr"][:]
        indices_ds = counts["indices"]
        data_ds = counts["data"]
        for pos, cell_index in enumerate(cell_indices):
            start = int(indptr[cell_index])
            end = int(indptr[cell_index + 1])
            seq = tokenize_count_row(
                row_counts=data_ds[start:end],
                row_feature_indices=indices_ds[start:end],
                feature_map=feature_map,
                cls_token_id=cls_token_id,
                eos_token_id=eos_token_id,
                model_input_size=model_input_size,
            )
            if seq is not None:
                sequences.append(seq)
                valid_positions.append(pos)
    return sequences, valid_positions


# ---------------------------------------------------------------------------
# Baseline logit computation
# ---------------------------------------------------------------------------

@torch.no_grad()
def compute_baseline_logits(
    clf: GeneformerClassifier,
    sequences: list[np.ndarray],
    pad_token_id: int,
    device: torch.device,
    batch_size: int,
) -> np.ndarray:
    """Return Responder class probability for each cell (index 1 of softmax)."""
    all_probs = []
    clf.eval()
    for start in range(0, len(sequences), batch_size):
        batch = sequences[start : start + batch_size]
        input_ids, attention_mask = pad_token_batch(batch, pad_token_id, device)
        with autocast(enabled=device.type == "cuda"):
            logits = clf(input_ids, attention_mask)
        probs = torch.softmax(logits.float(), dim=-1)[:, 1].cpu().numpy()
        all_probs.append(probs)
    return np.concatenate(all_probs) if all_probs else np.array([])


# ---------------------------------------------------------------------------
# In silico perturbation kernel
# ---------------------------------------------------------------------------

@torch.no_grad()
def perturb_single_token(
    clf: GeneformerClassifier,
    sequences: list[np.ndarray],
    baseline_probs: np.ndarray,
    target_token: int,
    pad_token_id: int,
    device: torch.device,
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Zero out `target_token` in all sequences where it appears.
    Return (delta_logits, cell_mask) where cell_mask is True for cells carrying the gene.
    delta_logit = perturbed_prob - baseline_prob  (positive = gene suppresses R score)
    """
    clf.eval()
    perturbed_seqs = []
    cell_has_token = []

    for seq in sequences:
        if target_token in seq:
            new_seq = seq.copy()
            new_seq[new_seq == target_token] = pad_token_id
            perturbed_seqs.append(new_seq)
            cell_has_token.append(True)
        else:
            perturbed_seqs.append(None)
            cell_has_token.append(False)

    cell_mask = np.array(cell_has_token)
    carrying_seqs = [s for s, has in zip(perturbed_seqs, cell_has_token) if has]
    if not carrying_seqs:
        return np.array([]), cell_mask

    perturbed_probs = []
    for start in range(0, len(carrying_seqs), batch_size):
        batch = carrying_seqs[start : start + batch_size]
        input_ids, attention_mask = pad_token_batch(batch, pad_token_id, device)
        with autocast(enabled=device.type == "cuda"):
            logits = clf(input_ids, attention_mask)
        probs = torch.softmax(logits.float(), dim=-1)[:, 1].cpu().numpy()
        perturbed_probs.append(probs)
    perturbed_probs_arr = np.concatenate(perturbed_probs)
    baseline_carrying = baseline_probs[cell_mask]
    delta = perturbed_probs_arr - baseline_carrying
    return delta, cell_mask


# ---------------------------------------------------------------------------
# Bootstrap CI for top genes
# ---------------------------------------------------------------------------

def bootstrap_mean_ci(values: np.ndarray, n: int = 500, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    means = [rng.choice(values, size=len(values), replace=True).mean() for _ in range(n)]
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


# ---------------------------------------------------------------------------
# Token -> gene symbol reverse mapping
# ---------------------------------------------------------------------------

def build_token_to_symbol(token_dict: dict[str, int], ensembl_map: dict[str, str]) -> dict[int, str]:
    """Build reverse map: token_id -> gene symbol (via ensembl_id)."""
    ensembl_to_symbol = {v: k for k, v in ensembl_map.items()}
    token_to_symbol: dict[int, str] = {}
    for key, token_id in token_dict.items():
        if key.startswith("ENSG"):
            symbol = ensembl_to_symbol.get(key, key)
            token_to_symbol[int(token_id)] = symbol
    return token_to_symbol


def load_ensembl_map(asset_dir: Path) -> dict[str, str]:
    import pickle
    map_path = asset_dir / "geneformer" / "ensembl_mapping_dict_gc104M.pkl"
    if not map_path.exists():
        return {}
    with map_path.open("rb") as handle:
        return pickle.load(handle)


# ---------------------------------------------------------------------------
# Per-target perturbation analysis
# ---------------------------------------------------------------------------

def run_perturbation(
    clf: GeneformerClassifier,
    atlas_path: Path,
    target: str,
    datasets: tuple[str, ...],
    feature_map,
    token_dict: dict[str, int],
    cls_token_id: int,
    eos_token_id: int,
    pad_token_id: int,
    model_input_size: int,
    token_to_symbol: dict[int, str],
    device: torch.device,
    batch_size: int,
    max_cells: int,
    freq_threshold: float,
    rng: np.random.Generator,
    smoke_test: bool,
) -> pd.DataFrame:
    print(f"\n  Target: {target}")
    cell_indices, patient_ids, responses = select_cells_by_type(atlas_path, target, datasets, max_cells, rng)
    n_cells = len(cell_indices)
    print(f"  Cells selected: {n_cells}")

    if n_cells < 5:
        print(f"  [WARN] Too few cells for {target}, skipping.")
        return pd.DataFrame()

    if n_cells < 20:
        print(f"  [WARN] Only {n_cells} cells — results may be unstable.")

    sequences, valid_positions = tokenize_cells(
        atlas_path, cell_indices, feature_map, cls_token_id, eos_token_id, model_input_size
    )
    n_tokenized = len(sequences)
    print(f"  Tokenized: {n_tokenized}/{n_cells}")
    if n_tokenized == 0:
        return pd.DataFrame()

    # Expand patient/response to match tokenized cells
    valid_responses = [responses[p] for p in valid_positions]

    baseline_probs = compute_baseline_logits(clf, sequences, pad_token_id, device, batch_size)
    print(f"  Baseline Responder prob: {baseline_probs.mean():.3f} ± {baseline_probs.std():.3f}")

    # Collect unique tokens, filter by frequency threshold
    token_counts: dict[int, int] = defaultdict(int)
    for seq in sequences:
        unique_tokens = set(seq.tolist())
        for t in unique_tokens:
            if t not in (cls_token_id, eos_token_id, pad_token_id):
                token_counts[t] += 1

    min_cells_with_gene = max(1, int(freq_threshold * n_tokenized))
    candidate_tokens = [t for t, cnt in token_counts.items() if cnt >= min_cells_with_gene]
    print(f"  Unique gene tokens: {len(token_counts):,}  "
          f"(freq>={freq_threshold:.0%}: {len(candidate_tokens):,})")

    if smoke_test:
        # Only perturb first 50 genes in smoke test
        candidate_tokens = candidate_tokens[:50]
        print(f"  [SMOKE TEST] Restricting to {len(candidate_tokens)} genes")

    # Perturbation loop
    rows = []
    t0 = time.time()
    for token_i, token_id in enumerate(candidate_tokens):
        delta, cell_mask = perturb_single_token(
            clf, sequences, baseline_probs, token_id, pad_token_id, device, batch_size
        )
        if len(delta) == 0:
            continue
        mean_delta = float(delta.mean())
        std_delta = float(delta.std())
        frac = float(cell_mask.sum()) / n_tokenized
        gene_symbol = token_to_symbol.get(token_id, f"token_{token_id}")

        ci_lo, ci_hi = float("nan"), float("nan")
        if len(delta) >= 5:
            ci_lo, ci_hi = bootstrap_mean_ci(delta, n=200, seed=token_i)

        rows.append({
            "gene_symbol": gene_symbol,
            "token_id": token_id,
            "mean_delta_prob": mean_delta,
            "std_delta_prob": std_delta,
            "abs_mean_delta": abs(mean_delta),
            "frac_cells_with_gene": frac,
            "n_cells_with_gene": int(cell_mask.sum()),
            "ci_lo": ci_lo,
            "ci_hi": ci_hi,
        })

        if (token_i + 1) % 200 == 0:
            elapsed = time.time() - t0
            print(f"  Perturbed {token_i + 1}/{len(candidate_tokens)} genes ({elapsed:.0f}s)", flush=True)

    df = pd.DataFrame(rows).sort_values("abs_mean_delta", ascending=False).reset_index(drop=True)
    df.insert(0, "rank", range(1, len(df) + 1))
    df.insert(1, "cell_type", target)
    df.insert(2, "n_cells_total", n_tokenized)
    print(f"  Perturbation complete: {len(df):,} genes, top gene: "
          f"{df.iloc[0]['gene_symbol'] if len(df) > 0 else 'N/A'} "
          f"(Δprob={df.iloc[0]['mean_delta_prob']:.4f})" if len(df) > 0 else "  No genes perturbed.")
    return df


# ---------------------------------------------------------------------------
# Figure: dot plot of top-20 genes
# ---------------------------------------------------------------------------

def plot_top20(df: pd.DataFrame, target: str, out_path: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        top = df.head(20).copy()
        if top.empty:
            return

        top = top.iloc[::-1]  # reverse for horizontal barh
        colors = ["#CC3333" if v < 0 else "#3366CC" for v in top["mean_delta_prob"]]

        fig, ax = plt.subplots(figsize=(7, max(4, len(top) * 0.35 + 1)), constrained_layout=True)
        bars = ax.barh(top["gene_symbol"], top["mean_delta_prob"], color=colors, alpha=0.85)

        # Error bars (CI)
        has_ci = top["ci_lo"].notna() & top["ci_hi"].notna()
        if has_ci.any():
            errs_lo = (top["mean_delta_prob"] - top["ci_lo"]).clip(lower=0).values
            errs_hi = (top["ci_hi"] - top["mean_delta_prob"]).clip(lower=0).values
            ax.errorbar(
                top["mean_delta_prob"][has_ci],
                top["gene_symbol"][has_ci],
                xerr=[errs_lo[has_ci], errs_hi[has_ci]],
                fmt="none", color="black", capsize=3, linewidth=0.8,
            )

        ax.axvline(0, color="gray", lw=0.8, ls="--")
        ax.set_xlabel("Mean Δ Responder probability (perturbation vs baseline)")
        ax.set_title(f"{target}\nTop-20 genes by |Δprob|", fontsize=10, fontweight="bold")

        # Fraction annotation
        for i, (_, row) in enumerate(top.iterrows()):
            ax.text(
                ax.get_xlim()[1] * 0.01 + (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.98,
                i,
                f"f={row['frac_cells_with_gene']:.2f}",
                va="center", ha="right", fontsize=6.5, color="gray",
            )

        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
    except Exception as exc:
        print(f"  [WARN] Figure for {target} failed: {exc}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="In silico gene perturbation for IO response")
    parser.add_argument("--atlas", type=Path, default=Path("data/processed/integrated_io_atlas.h5ad"))
    parser.add_argument("--asset-dir", type=Path, default=Path("data/references/geneformer"))
    parser.add_argument("--model-variant", type=str, default=DEFAULT_MODEL_VARIANT)
    parser.add_argument("--gtf", type=Path, default=Path("data/references/gencode.v44.basic.annotation.gtf.gz"))
    parser.add_argument(
        "--finetune-checkpoint",
        type=Path,
        default=Path("out/05_multilayer/models/geneformer_finetuned_final/pytorch_model.bin"),
        help="Fine-tuned model checkpoint (from 05j). Falls back to frozen if missing.",
    )
    parser.add_argument("--targets", nargs="+", default=DEFAULT_TARGETS)
    parser.add_argument("--max-cells", type=int, default=150, help="Max cells per target (GPU memory)")
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument(
        "--freq-threshold", type=float, default=0.05,
        help="Min fraction of cells carrying a gene for it to be tested (0.05 = 5%%)",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--force-cpu", action="store_true")
    parser.add_argument("--smoke-test", action="store_true", help="Test with 20 cells and 50 genes per target")
    parser.add_argument(
        "--datasets", nargs="+", default=["GSE206325", "Cappuyns2023"],
        help="Datasets to draw target cells from",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    seed_everything(args.seed)

    if args.smoke_test:
        args.max_cells = min(args.max_cells, 20)
        print("[SMOKE TEST] max_cells=20, 50 genes per target")

    device = resolve_device(force_cpu=args.force_cpu)

    PROJECT_ROOT = Path(__file__).resolve().parents[2]
    asset_dir = (PROJECT_ROOT / args.asset_dir).resolve() if not args.asset_dir.is_absolute() else args.asset_dir.resolve()
    atlas_path = (PROJECT_ROOT / args.atlas).resolve() if not args.atlas.is_absolute() else args.atlas.resolve()
    gtf_path = (PROJECT_ROOT / args.gtf).resolve() if not args.gtf.is_absolute() else args.gtf.resolve()
    model_dir = asset_dir / args.model_variant
    finetune_path = (PROJECT_ROOT / args.finetune_checkpoint).resolve() if not args.finetune_checkpoint.is_absolute() else args.finetune_checkpoint.resolve()

    out_tables = PROJECT_ROOT / "out" / "05_multilayer" / "tables"
    out_figs = PROJECT_ROOT / "out" / "05_multilayer" / "figures"
    for d in (out_tables, out_figs):
        d.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("[05k] Geneformer in silico perturbation")
    print("=" * 70)
    print(f"Device         : {device}")
    if device.type == "cuda":
        props = torch.cuda.get_device_properties(device)
        print(f"GPU            : {props.name} ({props.total_memory / 1e9:.1f} GB)")
    print(f"Model          : {args.model_variant}")
    print(f"Checkpoint     : {finetune_path}")
    print(f"Targets        : {', '.join(args.targets)}")
    print(f"Max cells      : {args.max_cells}  Freq threshold: {args.freq_threshold:.0%}")

    # --- Dictionaries & feature map ---
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

    # Reverse token map for symbol annotation
    ensembl_map = load_ensembl_map(asset_dir)
    token_to_symbol = build_token_to_symbol(token_dict, ensembl_map)

    # --- Load model config ---
    config_path = model_dir / "config.json"
    with config_path.open() as handle:
        model_config = json.load(handle)
    hidden_size = int(model_config["hidden_size"])

    # --- Load classifier ---
    clf = load_classifier(model_dir, finetune_path, device, hidden_size)

    rng = np.random.default_rng(args.seed)
    datasets_tuple = tuple(args.datasets)

    # --- Run per-target perturbation ---
    all_results = []
    qc_rows = []
    t_start = time.time()

    for target in args.targets:
        t0_target = time.time()
        df_target = run_perturbation(
            clf=clf,
            atlas_path=atlas_path,
            target=target,
            datasets=datasets_tuple,
            feature_map=feature_map,
            token_dict=token_dict,
            cls_token_id=cls_token_id,
            eos_token_id=eos_token_id,
            pad_token_id=pad_token_id,
            model_input_size=model_input_size,
            token_to_symbol=token_to_symbol,
            device=device,
            batch_size=args.batch_size,
            max_cells=args.max_cells,
            freq_threshold=args.freq_threshold,
            rng=rng,
            smoke_test=args.smoke_test,
        )

        target_safe = target.replace("/", "_").replace(" ", "_")
        elapsed_target = time.time() - t0_target

        if not df_target.empty:
            out_csv = out_tables / f"geneformer_perturbation_{target_safe}.csv"
            df_target.to_csv(out_csv, index=False)
            print(f"  Saved: {out_csv}")

            fig_path = out_figs / f"05k_perturbation_top20_{target_safe}.png"
            plot_top20(df_target, target, fig_path)
            print(f"  Figure: {fig_path}")
            all_results.append(df_target)

        qc_rows.append({
            "target": target,
            "n_genes_tested": len(df_target) if not df_target.empty else 0,
            "n_cells": df_target["n_cells_total"].iloc[0] if not df_target.empty else 0,
            "top_gene": df_target["gene_symbol"].iloc[0] if not df_target.empty else "N/A",
            "top_gene_delta": df_target["mean_delta_prob"].iloc[0] if not df_target.empty else float("nan"),
            "elapsed_min": round(elapsed_target / 60, 2),
        })

    # --- Combined summary ---
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined_path = out_tables / "geneformer_perturbation_all_targets.csv"
        combined.to_csv(combined_path, index=False)
        print(f"\nCombined table: {combined_path}")

    # --- QC record ---
    qc_df = pd.DataFrame(qc_rows)
    qc_df["run_timestamp_utc"] = datetime.now(timezone.utc).isoformat()
    qc_df["model_variant"] = args.model_variant
    qc_df["finetune_checkpoint"] = str(finetune_path)
    qc_df["device"] = str(device)
    qc_df["smoke_test"] = args.smoke_test
    qc_df["total_elapsed_min"] = round((time.time() - t_start) / 60, 2)
    qc_path = out_tables / "geneformer_perturbation_qc.csv"
    qc_df.to_csv(qc_path, index=False)
    print(f"QC table       : {qc_path}")

    print("\n" + "=" * 70)
    print("[05k] COMPLETE")
    print(f"  Targets analysed: {len(qc_rows)}")
    print(f"  Total time: {(time.time() - t_start) / 60:.1f} min")
    print("=" * 70)


if __name__ == "__main__":
    main()
