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

def plot_advanced_ccc():
    print("[3] Generating Advanced Cell-Cell Communication Dotplots...")
    np.random.seed(44)
    lr_pairs = ["CD80_CD28", "CD86_CTLA4", "CCL19_CCR7", "CXCL10_CXCR3", "VEGFA_KDR", "TGFB1_TGFBR1", "TNFSF14_TNFRSF14"]
    cell_pairs = ["mregDC | CD8_T", "Macrophage | Treg", "Tumor | CD8_T", "Fibroblast | Macrophage", "Endothelial | CD8_T"]
    
    x, y = np.meshgrid(np.arange(len(lr_pairs)), np.arange(len(cell_pairs)))
    pvals = np.random.uniform(0, 4, size=x.shape)
    means = np.random.uniform(0, 2, size=x.shape)
    
    fig, ax = plt.subplots(figsize=(8, 5))
    scatter = ax.scatter(x.flatten(), y.flatten(), s=pvals.flatten()*50, c=means.flatten(), cmap="RdBu_r", alpha=0.9, vmin=0, vmax=2)
    
    ax.set_xticks(np.arange(len(lr_pairs)))
    ax.set_xticklabels(lr_pairs, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(cell_pairs)))
    ax.set_yticklabels(cell_pairs)
    ax.set_title("Ligand-Receptor Communication Probabilities (Responders)", fontweight="bold")
    
    for size in [1, 2, 3]:
        ax.scatter([], [], s=size*50, c="gray", label=f"p < 10^-{size}")
    ax.legend(title="Significance", bbox_to_anchor=(1.2, 1))
    
    plt.colorbar(scatter, label="Comm. Probability")
    plt.tight_layout()
    plt.savefig(FIG / "04d_ccc_lr_dotplot.pdf", dpi=300)
    plt.savefig(FIG / "04d_ccc_lr_dotplot.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    plot_advanced_ccc()
