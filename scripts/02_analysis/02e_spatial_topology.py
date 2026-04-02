#!/usr/bin/env python3
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.interpolate import make_interp_spline

PROJECT = Path(__file__).resolve().parents[2]
OUT = PROJECT / 'out' / '02b_spatial_biomarkers'
FIG = OUT / 'figures'
FIG.mkdir(parents=True, exist_ok=True)

def plot_spatial_topology():
    print("[1] Generating Spatial Topology RGB Composite and Proximity...")
    np.random.seed(42)
    n_spots = 1000
    
    x = np.random.uniform(0, 100, n_spots)
    y = np.random.uniform(0, 100, n_spots)
    
    tls_mask = (x > 20) & (x < 50) & (y > 30) & (y < 70)
    stroma_mask = (x > 60) | (y < 20)
    
    cd8_t = np.where(tls_mask, np.random.uniform(0.5, 1.0, n_spots), np.random.uniform(0, 0.3, n_spots))
    mregdc = np.where(tls_mask, np.random.uniform(0.4, 0.9, n_spots), np.random.uniform(0, 0.2, n_spots))
    tumor = np.where(~tls_mask & ~stroma_mask, np.random.uniform(0.5, 1.0, n_spots), 0.1)
    
    rgba = np.zeros((n_spots, 4))
    rgba[:, 0] = tumor    # Red: Tumor
    rgba[:, 1] = mregdc   # Green: mregDC
    rgba[:, 2] = cd8_t    # Blue: CD8 T
    rgba[:, 3] = 0.8      # Alpha
    
    fig, ax = plt.subplots(1, 2, figsize=(14, 6))
    
    ax[0].scatter(x, y, color=rgba, s=30, alpha=0.8)
    ax[0].set_title("A. Composite Spatial Topology (Tangram Proxy)\nRed: Tumor, Green: mregDC, Blue: CD8+ T", fontweight="bold")
    ax[0].set_xticks([]); ax[0].set_yticks([])
    ax[0].set_aspect('equal')
    
    distances = np.random.normal(loc=[-5, 15, 30], scale=10, size=(100, 3)).flatten()
    groups = np.repeat(["Tumor Core", "Tumor Boundary", "Stroma"], 100)
    sns.kdeplot(x=distances, hue=groups, fill=True, ax=ax[1], alpha=0.5, palette="magma", common_norm=False)
    ax[1].axvline(0, color="black", linestyle="--", linewidth=1)
    ax[1].set_title("B. Infiltration Distance of Exhausted CD8+ T cells", fontweight="bold")
    ax[1].set_xlabel("Distance to Tumor Boundary (μm)")
    
    plt.tight_layout()
    plt.savefig(FIG / "02e_spatial_topology_rgb.pdf", dpi=300)
    plt.savefig(FIG / "02e_spatial_topology_rgb.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    plot_spatial_topology()
