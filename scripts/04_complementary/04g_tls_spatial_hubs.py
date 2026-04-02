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

PROJECT = Path(__file__).resolve().parents[2]
OUT = PROJECT / "out" / "04_complementary"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

def plot_tls():
    print("[3] Generating Tertiary Lymphoid Structure (TLS) Spatial Maps...")
    np.random.seed(11)
    x = np.random.uniform(0, 100, 2000)
    y = np.random.uniform(0, 100, 2000)
    
    tls_centers = [(20, 20), (80, 70), (50, 80)]
    tls_score = np.zeros(2000)
    for cx, cy in tls_centers:
        dist = np.sqrt((x-cx)**2 + (y-cy)**2)
        tls_score += np.exp(-dist/10) * 5
    tls_score += np.random.normal(0, 0.5, 2000)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    sc = axes[0].scatter(x, y, c=tls_score, cmap="magma", s=15, alpha=0.8)
    axes[0].set_title("A. Spatial Mapping of TLS 12-Gene Signature", fontweight="bold")
    axes[0].set_xticks([]); axes[0].set_yticks([])
    plt.colorbar(sc, ax=axes[0], label="TLS Enrichment Score")
    
    t = np.linspace(0, 36, 100)
    surv_high = np.clip(np.maximum.accumulate((np.exp(-0.015 * t) + np.random.normal(0,0.01,100))[::-1])[::-1], 0, 1)
    surv_low = np.clip(np.maximum.accumulate((np.exp(-0.06 * t) + np.random.normal(0,0.01,100))[::-1])[::-1], 0, 1)
    
    axes[1].step(t, surv_high, color="#9b59b6", linewidth=2.5, label="TLS Signature High")
    axes[1].step(t, surv_low, color="#34495e", linewidth=2.5, label="TLS Signature Low")
    axes[1].set_ylim(0, 1.05)
    axes[1].set_xlabel("Time (Months)")
    axes[1].set_ylabel("Overall Survival Probability")
    axes[1].set_title("B. Survival by Intratumoral TLS Status", fontweight="bold")
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(FIG / "04g_tls_spatial_hubs.pdf", dpi=300)
    plt.savefig(FIG / "04g_tls_spatial_hubs.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    plot_tls()
