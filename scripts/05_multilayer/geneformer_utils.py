#!/usr/bin/env python3
"""Utilities for local Geneformer benchmarking on the integrated atlas."""

from __future__ import annotations

import gzip
import pickle
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable

import h5py
import numpy as np
import torch

if TYPE_CHECKING:
    from transformers import AutoModel


DEFAULT_MODEL_VARIANT = "Geneformer-V2-104M_CLcancer"
DEFAULT_DATASETS = ("GSE206325", "Cappuyns2023")
MODEL_INPUT_SIZES = {
    "Geneformer-V1-10M": 2048,
    "Geneformer-V2-104M": 4096,
    "Geneformer-V2-316M": 4096,
    "Geneformer-V2-104M_CLcancer": 4096,
}
V2_COMMON_FILES = (
    "geneformer/token_dictionary_gc104M.pkl",
    "geneformer/gene_median_dictionary_gc104M.pkl",
    "geneformer/ensembl_mapping_dict_gc104M.pkl",
)
GTF_ATTRIBUTE_RE = re.compile(r'(\S+) "([^"]+)";')


@dataclass
class FeatureMap:
    token_by_feature: np.ndarray
    median_by_feature: np.ndarray
    n_mapped_genes: int
    n_total_genes: int


@dataclass
class SelectedCells:
    cell_indices: np.ndarray
    patient_ids: np.ndarray
    responses: np.ndarray
    datasets: np.ndarray


def sanitize_model_name(model_variant: str) -> str:
    return model_variant.lower().replace("-", "_")


def default_output_tag(model_variant: str, max_cells_per_patient: int | None) -> str:
    model_tag = sanitize_model_name(model_variant)
    if max_cells_per_patient is None:
        return model_tag
    return f"{model_tag}_subsample_{max_cells_per_patient}cpp"


def decode_string_array(dataset: h5py.Dataset) -> list[str]:
    values = dataset[:]
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def decode_categorical_values(group: h5py.Group, indices: np.ndarray | None = None) -> np.ndarray:
    categories = decode_string_array(group["categories"])
    codes = group["codes"][:] if indices is None else group["codes"][indices]
    return np.array([categories[code] if code >= 0 else None for code in codes], dtype=object)


def get_model_input_size(model_variant: str) -> int:
    if model_variant not in MODEL_INPUT_SIZES:
        valid = ", ".join(sorted(MODEL_INPUT_SIZES))
        raise ValueError(f"Unsupported Geneformer model '{model_variant}'. Valid options: {valid}")
    return MODEL_INPUT_SIZES[model_variant]


def load_v2_gene_dictionaries(asset_dir: Path) -> tuple[dict[str, int], dict[str, float]]:
    token_path = asset_dir / "geneformer" / "token_dictionary_gc104M.pkl"
    median_path = asset_dir / "geneformer" / "gene_median_dictionary_gc104M.pkl"
    if not token_path.exists() or not median_path.exists():
        raise FileNotFoundError(
            "Geneformer V2 dictionaries are missing. Run 05g_download_geneformer_assets.py first."
        )
    with token_path.open("rb") as handle:
        token_dict = pickle.load(handle)
    with median_path.open("rb") as handle:
        median_dict = pickle.load(handle)
    return token_dict, median_dict


def build_symbol_to_ensembl_map(gtf_path: Path, cache_path: Path | None = None) -> dict[str, str]:
    if cache_path is not None and cache_path.exists():
        with cache_path.open("rb") as handle:
            return pickle.load(handle)

    symbol_to_ensembl: dict[str, str] = {}
    open_fn = gzip.open if gtf_path.suffix == ".gz" else open
    with open_fn(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = dict(GTF_ATTRIBUTE_RE.findall(fields[8]))
            gene_id = attrs.get("gene_id")
            gene_name = attrs.get("gene_name")
            if not gene_id or not gene_name:
                continue
            symbol_to_ensembl.setdefault(gene_name, gene_id.split(".")[0])

    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with cache_path.open("wb") as handle:
            pickle.dump(symbol_to_ensembl, handle)
    return symbol_to_ensembl


def load_var_names(h5_path: Path) -> np.ndarray:
    with h5py.File(h5_path, "r") as handle:
        return np.array(decode_string_array(handle["var"]["_index"]), dtype=object)


def build_feature_map(
    var_names: Iterable[str],
    symbol_to_ensembl: dict[str, str],
    token_dict: dict[str, int],
    median_dict: dict[str, float],
) -> FeatureMap:
    var_names = list(var_names)
    token_by_feature = np.full(len(var_names), -1, dtype=np.int32)
    median_by_feature = np.zeros(len(var_names), dtype=np.float32)

    mapped = 0
    for index, symbol in enumerate(var_names):
        ensembl_id = symbol_to_ensembl.get(str(symbol))
        if ensembl_id is None:
            continue
        token = token_dict.get(ensembl_id)
        median = median_dict.get(ensembl_id)
        if token is None or median is None or median <= 0:
            continue
        token_by_feature[index] = int(token)
        median_by_feature[index] = float(median)
        mapped += 1

    return FeatureMap(
        token_by_feature=token_by_feature,
        median_by_feature=median_by_feature,
        n_mapped_genes=mapped,
        n_total_genes=len(var_names),
    )


def select_cells(
    h5_path: Path,
    datasets: tuple[str, ...] = DEFAULT_DATASETS,
    responses: tuple[str, ...] = ("Responder", "Non-Responder"),
    max_cells_per_patient: int | None = None,
    seed: int = 0,
) -> SelectedCells:
    with h5py.File(h5_path, "r") as handle:
        obs = handle["obs"]
        response_codes = obs["response"]["codes"][:]
        response_categories = decode_string_array(obs["response"]["categories"])
        dataset_codes = obs["dataset"]["codes"][:]
        dataset_categories = decode_string_array(obs["dataset"]["categories"])

        response_code_set = {i for i, value in enumerate(response_categories) if value in responses}
        dataset_code_set = {i for i, value in enumerate(dataset_categories) if value in datasets}
        mask = np.isin(response_codes, list(response_code_set)) & np.isin(dataset_codes, list(dataset_code_set))
        cell_indices = np.flatnonzero(mask)

        patient_key = "patient_ID" if "patient_ID" in obs else "Patient_ID"
        patient_ids = decode_categorical_values(obs[patient_key], cell_indices)
        response_values = decode_categorical_values(obs["response"], cell_indices)
        dataset_values = decode_categorical_values(obs["dataset"], cell_indices)

    if max_cells_per_patient is not None:
        rng = np.random.default_rng(seed)
        keep = []
        for patient_id in np.unique(patient_ids):
            patient_loc = np.flatnonzero(patient_ids == patient_id)
            if len(patient_loc) > max_cells_per_patient:
                patient_loc = np.sort(rng.choice(patient_loc, size=max_cells_per_patient, replace=False))
            keep.append(patient_loc)
        keep_idx = np.concatenate(keep) if keep else np.array([], dtype=int)
        cell_indices = cell_indices[keep_idx]
        patient_ids = patient_ids[keep_idx]
        response_values = response_values[keep_idx]
        dataset_values = dataset_values[keep_idx]

    return SelectedCells(
        cell_indices=cell_indices,
        patient_ids=patient_ids,
        responses=response_values,
        datasets=dataset_values,
    )


def resolve_device(force_cpu: bool = False) -> torch.device:
    if force_cpu or not torch.cuda.is_available():
        return torch.device("cpu")
    return torch.device("cuda")


def seed_everything(seed: int) -> None:
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    if hasattr(torch.backends, "cudnn"):
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


def load_model(model_dir: Path, device: torch.device) -> Any:
    from transformers import AutoModel
    import torch

    model = AutoModel.from_pretrained(
        str(model_dir),
        trust_remote_code=False,
        attn_implementation="sdpa",
        dtype=torch.float16,
    )
    model.eval()
    model.to(device)
    return model


def tokenize_count_row(
    row_counts: np.ndarray,
    row_feature_indices: np.ndarray,
    feature_map: FeatureMap,
    cls_token_id: int,
    eos_token_id: int,
    model_input_size: int,
    target_sum: float = 10000.0,
) -> np.ndarray | None:
    if row_counts.size == 0:
        return None
    mapped_tokens = feature_map.token_by_feature[row_feature_indices]
    keep = mapped_tokens >= 0
    if not np.any(keep):
        return None

    row_counts = np.asarray(row_counts[keep], dtype=np.float32)
    row_feature_indices = row_feature_indices[keep]
    total_counts = float(row_counts.sum())
    if total_counts <= 0:
        return None

    scaled = row_counts / total_counts * target_sum
    medians = feature_map.median_by_feature[row_feature_indices]
    rank_scores = scaled / medians
    order = np.argsort(-rank_scores, kind="stable")
    tokens = mapped_tokens[keep][order][: model_input_size - 2]
    if tokens.size == 0:
        return None

    return np.concatenate(
        (
            np.array([cls_token_id], dtype=np.int32),
            np.asarray(tokens, dtype=np.int32),
            np.array([eos_token_id], dtype=np.int32),
        )
    )


def pad_token_batch(sequences: list[np.ndarray], pad_token_id: int, device: torch.device) -> tuple[torch.Tensor, torch.Tensor]:
    max_len = max(len(sequence) for sequence in sequences)
    batch = np.full((len(sequences), max_len), pad_token_id, dtype=np.int64)
    mask = np.zeros((len(sequences), max_len), dtype=np.int64)
    for row_index, sequence in enumerate(sequences):
        seq_len = len(sequence)
        batch[row_index, :seq_len] = sequence
        mask[row_index, :seq_len] = 1
    return torch.from_numpy(batch).to(device), torch.from_numpy(mask).to(device)


def infer_cls_embeddings(
    model: Any,
    sequences: list[np.ndarray],
    pad_token_id: int,
    device: torch.device,
) -> np.ndarray:
    input_ids, attention_mask = pad_token_batch(sequences, pad_token_id, device)
    with torch.no_grad():
        outputs = model(input_ids=input_ids, attention_mask=attention_mask)
        cls_embeddings = outputs.last_hidden_state[:, 0, :]
    return cls_embeddings.detach().cpu().numpy()