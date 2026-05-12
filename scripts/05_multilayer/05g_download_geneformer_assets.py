#!/usr/bin/env python3
"""
05g_download_geneformer_assets.py - Download the Geneformer assets needed for
local patient-level benchmarking.

This script keeps the download intentionally narrow: one or more model
subfolders plus the V2 tokenizer dictionaries required for rank-encoded input.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from huggingface_hub import snapshot_download

from geneformer_utils import (
    DEFAULT_MODEL_VARIANT,
    V2_COMMON_FILES,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download Geneformer model and tokenizer assets")
    parser.add_argument(
        "--asset-dir",
        type=Path,
        default=Path("data/references/geneformer"),
        help="Local directory where the Geneformer files will be stored.",
    )
    parser.add_argument(
        "--model-variant",
        dest="model_variants",
        action="append",
        default=None,
        help="Geneformer model subfolder to download. Repeat to fetch multiple variants.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    model_variants = args.model_variants or [DEFAULT_MODEL_VARIANT]
    allow_patterns = list(V2_COMMON_FILES)
    for model_variant in model_variants:
        allow_patterns.extend(
            [
                f"{model_variant}/config.json",
                f"{model_variant}/generation_config.json",
                f"{model_variant}/model.safetensors",
                f"{model_variant}/training_args.bin",
            ]
        )

    asset_dir = args.asset_dir.resolve()
    asset_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("[05g] Download Geneformer assets")
    print("=" * 70)
    print(f"Target directory : {asset_dir}")
    print(f"Model variants   : {', '.join(model_variants)}")
    print(f"Files requested  : {len(allow_patterns)}")

    snapshot_download(
        repo_id="ctheodoris/Geneformer",
        repo_type="model",
        allow_patterns=allow_patterns,
        local_dir=str(asset_dir),
    )

    manifest = {
        "repo_id": "ctheodoris/Geneformer",
        "model_variants": model_variants,
        "common_files": list(V2_COMMON_FILES),
    }
    manifest_path = asset_dir / "download_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"Saved manifest   : {manifest_path}")
    for model_variant in model_variants:
        model_dir = asset_dir / model_variant
        print(f"Model directory  : {model_dir}")
    print("[DONE]")


if __name__ == "__main__":
    main()