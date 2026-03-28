#!/usr/bin/env python3
"""
Spatial Immunotherapy Response — Data Acquisition
Downloads all public HCC immunotherapy scRNA-seq datasets.
"""

import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"

IO_DATASETS = [
    {"geo": "GSE206325", "study": "Liu2022_antiPD1", "treatment": "anti-PD1",
     "n_patients": 30, "cells_approx": 120000},
    {"geo": "GSE202069", "study": "Zhu2023_antiPDL1", "treatment": "anti-PD1/PD-L1",
     "n_patients": 14, "cells_approx": 80000},
    {"geo": "GSE234241", "study": "Xue2023_combo", "treatment": "anti-PD1+TKI",
     "n_patients": 20, "cells_approx": 150000},
    {"geo": "GSE125449", "study": "Ma2019_antiPD1", "treatment": "anti-PD1",
     "n_patients": 8, "cells_approx": 25000},
]


def download_geo(geo_id: str, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        print(f"  [SKIP] {geo_id}")
        return True

    ftp = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_id[:-3]}nnn/{geo_id}/suppl/"
    result = subprocess.run([
        "wget", "-r", "-np", "-nd", "-P", str(output_dir),
        "-A", "*.tar.gz,*.h5,*.h5ad,*.csv.gz,*.tsv.gz,*.mtx.gz",
        "--timeout=300", "--tries=3", ftp
    ], capture_output=True, timeout=600)
    return result.returncode == 0


def main():
    print("=" * 70)
    print("Spatial Immunotherapy Response — Data Acquisition")
    print("=" * 70)

    for ds in IO_DATASETS:
        ds_dir = DATA_DIR / "pre_treatment" / f"{ds['geo']}_{ds['study']}"
        print(f"\n  {ds['geo']} ({ds['study']}): {ds['treatment']}, ~{ds['cells_approx']} cells")
        success = download_geo(ds["geo"], ds_dir)
        print(f"  Status: {'OK' if success else 'FAILED'}")

    print("\n[DONE]")


if __name__ == "__main__":
    main()
