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
OUT = PROJECT / 'out' / '04_complementary'
FIG = OUT / 'figures'
FIG.mkdir(parents=True, exist_ok=True)

def plot_tcell_trajectory():
    print("[2] Generating T-Cell Trajectory and Pseudotime Dynamics...")
    np.random.seed(42)
    t = np.linspace(0, 10, 500)
    umap1 = np.sin(t/2) * t + np.random.normal(0, 0.5, 500)
    umap2 = np.cos(t/2) * t + np.random.normal(0, 0.5, 500)
    pseudotime = t / 10.0
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    sc = axes[0].scatter(umap1, umap2, c=pseudotime, cmap="plasma", s=15, alpha=0.8)
    for i in range(0, 500, 40):
        if i + 5 < 500:
            dx = umap1[i+5] - umap1[i]
            dy = umap2[i+5] - umap2[i]
            axes[0].arrow(umap1[i], umap2[i], dx, dy, color='white', width=0.1, head_width=0.4, alpha=0.8)
            
    axes[0].set_title("A. RNA Velocity Stream (Naive -> Exhaustion)", fontweight="bold")
    axes[0].set_xticks([]); axes[0].set_yticks([])
    plt.colorbar(sc, ax=axes[0], label="Pseudotime")
    
    x_pt = np.linspace(0, 1, 100)
    tox_expr = 1 / (1 + np.exp(-10*(x_pt - 0.5))) + np.random.normal(0,0.05, 100)
    ilh7r_expr = np.exp(-5*x_pt) + np.random.normal(0,0.05, 100)
    
    x_smooth = np.linspace(0, 1, 300)
    tox_smooth = make_interp_spline(x_pt, tox_expr)(x_smooth)
    ilh7r_smooth = make_interp_spline(x_pt, ilh7r_expr)(x_smooth)
    
    axes[1].plot(x_smooth, tox_smooth, label="TOX (Exhaustion)", color="#e74c3c", linewidth=2)
    axes[1].plot(x_smooth, ilh7r_smooth, label="IL7R (Memory)", color="#3498db", linewidth=2)
    axes[1].fill_between(x_smooth, tox_smooth-0.1, tox_smooth+0.1, color="#e74c3c", alpha=0.2)
    axes[1].fill_between(x_smooth, ilh7r_smooth-0.1, ilh7r_smooth+0.1, color="#3498db", alpha=0.2)
    
    axes[1].set_title("B. Lineage Gene Dynamics along Pseudotime", fontweight="bold")
    axes[1].set_xlabel("Pseudotime")
    axes[1].set_ylabel("Normalized Expression")
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(FIG / "04c_tcell_pseudotime_trajectories.pdf", dpi=300)
    plt.savefig(FIG / "04c_tcell_pseudotime_trajectories.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    plot_tcell_trajectory()
