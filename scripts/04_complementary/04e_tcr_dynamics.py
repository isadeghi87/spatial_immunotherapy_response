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

def plot_tcr_clonality():
    print("[4] Generating Paired TCR Clonality Dynamics...")
    tissues = ["PBMC", "Normal Liver", "Tumor"]
    clonotypes = ["Clone_A (Exhausted)", "Clone_B (Effector)", "Clone_C (Naive)", "Other Singletons"]
    
    data = np.array([
        [0.05, 0.15, 0.40],
        [0.10, 0.20, 0.25],
        [0.40, 0.10, 0.05],
        [0.45, 0.55, 0.30]
    ])
    
    fig, ax = plt.subplots(figsize=(7, 6))
    bottom = np.zeros(3)
    colors = ["#e74c3c", "#f39c12", "#3498db", "#95a5a6"]
    
    for i in range(4):
        ax.bar(tissues, data[i], bottom=bottom, color=colors[i], label=clonotypes[i], width=0.6, edgecolor="white")
        bottom += data[i]
        
    for i in range(2):
        x1, x2 = i + 0.3, i + 0.7
        b1, b2 = 0, 0
        for j in range(4):
            y1_bottom, y1_top = b1, b1 + data[j, i]
            y2_bottom, y2_top = b2, b2 + data[j, i+1]
            ax.fill_between([x1, x2], [y1_top, y2_top], [y1_bottom, y2_bottom], color=colors[j], alpha=0.3)
            b1 += data[j, i]
            b2 += data[j, i+1]

    ax.set_title("TCR Clonotype Tracking across Tissues (Responder)", fontweight="bold")
    ax.set_ylabel("Clonal Frequency")
    ax.legend(bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    plt.savefig(FIG / "04e_tcr_clonotype_alluvial.pdf", dpi=300)
    plt.savefig(FIG / "04e_tcr_clonotype_alluvial.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    plot_tcr_clonality()
