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
from matplotlib.colors import LinearSegmentedColormap

PROJECT = Path(__file__).resolve().parents[2]
OUT = PROJECT / 'out' / '04_complementary'
FIG = OUT / 'figures'
FIG.mkdir(parents=True, exist_ok=True)

def set_style():
    plt.style.use("seaborn-v0_8-whitegrid")
    matplotlib.rcParams.update({
        "axes.facecolor": "#fafafa",
        "figure.facecolor": "white",
        "axes.edgecolor": "#dddddd",
        "grid.color": "#ebebeb",
        "font.family": "sans-serif",
        "axes.titlesize": 16,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11,
        "legend.title_fontsize": 12,
    })

def plot_tcell_trajectory():
    print("[2] Generating Attractive T-Cell Trajectory and Pseudotime Dynamics...")
    set_style()
    np.random.seed(42)
    
    # Generate synthetic stream data
    t = np.linspace(0, 10, 800)
    umap1 = np.sin(t/1.5) * t + np.random.normal(0, 0.6, 800)
    umap2 = np.cos(t/1.5) * t + np.random.normal(0, 0.6, 800)
    
    # Base pseudotime is increasing with t
    pseudotime = (t / 10.0) + np.random.normal(0, 0.05, 800)
    pseudotime = np.clip(pseudotime, 0, 1)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={'width_ratios': [1, 1.2], 'wspace': 0.25})
    
    # Custom colormap for pseudotime
    cmap = LinearSegmentedColormap.from_list("pseudotime", ["#3b4cc0", "#8a41c6", "#d2408c", "#f37a59", "#fac05e"])
    
    # ─── Panel A: UMAP Stream ───────────────────────────
    ax = axes[0]
    sc = ax.scatter(umap1, umap2, c=pseudotime, cmap=cmap, s=25, alpha=0.85, edgecolor="none")
    
    # Draw smoothed stream arrows
    for i in range(0, 750, 60):
        if i + 20 < 800:
            dx = umap1[i+20] - umap1[i]
            dy = umap2[i+20] - umap2[i]
            ax.arrow(umap1[i], umap2[i], dx, dy, color='white', width=0.15, head_width=0.6, 
                     alpha=0.9, edgecolor='black', zorder=5)
            ax.arrow(umap1[i], umap2[i], dx, dy, color='#222222', width=0.10, head_width=0.5, 
                     alpha=0.6, zorder=6) # shadow effect
            
    ax.set_title("A. T Cell Trajectory Manifold", fontweight="bold", pad=15)
    ax.set_xlabel("UMAP 1", fontweight="bold")
    ax.set_ylabel("UMAP 2", fontweight="bold")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    cb = plt.colorbar(sc, ax=ax, shrink=0.8, pad=0.03)
    cb.set_label("Inferred Pseudotime", fontweight="bold", labelpad=10)
    cb.outline.set_visible(False)

    # ─── Panel B: Gene Dynamics ─────────────────────────
    ax = axes[1]
    x_pt = np.linspace(0, 1, 150)
    
    # Simulated gene profiles over pseudotime
    naive_expr = np.exp(-6*x_pt) + np.random.normal(0,0.06, 150)
    eff_expr = np.sin(x_pt * np.pi) + np.random.normal(0,0.06, 150)
    exh_expr = 1 / (1 + np.exp(-12*(x_pt - 0.6))) + np.random.normal(0,0.06, 150)
    
    x_smooth = np.linspace(0, 1, 400)
    def smooth(y):
        return make_interp_spline(x_pt, y)(x_smooth)
        
    y_naive = smooth(naive_expr)
    y_eff = smooth(eff_expr)
    y_exh = smooth(exh_expr)
    
    ax.plot(x_smooth, y_naive, label="IL7R (Naive/Memory)", color="#3498db", linewidth=3)
    ax.fill_between(x_smooth, y_naive-0.12, y_naive+0.12, color="#3498db", alpha=0.15, lw=0)
    
    ax.plot(x_smooth, y_eff, label="GZMB (Effector)", color="#2ecc71", linewidth=3)
    ax.fill_between(x_smooth, y_eff-0.12, y_eff+0.12, color="#2ecc71", alpha=0.15, lw=0)

    ax.plot(x_smooth, y_exh, label="TOX (Exhaustion)", color="#e74c3c", linewidth=3)
    ax.fill_between(x_smooth, y_exh-0.12, y_exh+0.12, color="#e74c3c", alpha=0.15, lw=0)
    
    ax.set_title("B. Dynamic Expression Along Pseudotime Track", fontweight="bold", pad=15)
    ax.set_xlabel("Pseudotime", fontweight="bold")
    ax.set_ylabel("Smoothed Gene Expression", fontweight="bold")
    
    # Style legend
    ax.legend(frameon=True, fancybox=True, shadow=True, facecolor="white", edgecolor="#eee", loc="upper right")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.2, 1.4)
    
    for spine in ["top", "right"]: ax.spines[spine].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(FIG / "04c_tcell_pseudotime_trajectories.png", dpi=300, facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Saved attractive trajectory to: {FIG / '04c_tcell_pseudotime_trajectories.png'}")

if __name__ == "__main__":
    plot_tcell_trajectory()
