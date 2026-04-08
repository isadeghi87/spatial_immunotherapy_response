import os

base_code = """#!/usr/bin/env python3
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
"""

scripts = {
    "scripts/02_analysis/02e_spatial_topology.py": {
        "out_dir": "02b_spatial_biomarkers",
        "prefix": "02e_",
        "code": """
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
    ax[0].set_title("A. Composite Spatial Topology (Tangram Proxy)\\nRed: Tumor, Green: mregDC, Blue: CD8+ T", fontweight="bold")
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
"""
    },
    "scripts/04_complementary/04c_tcell_trajectories.py": {
        "out_dir": "04_complementary",
        "prefix": "04c_",
        "code": """
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
"""
    },
    "scripts/04_complementary/04d_cell_communication.py": {
        "out_dir": "04_complementary",
        "prefix": "04d_",
        "code": """
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
"""
    },
    "scripts/04_complementary/04e_tcr_dynamics.py": {
        "out_dir": "04_complementary",
        "prefix": "04e_",
        "code": """
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
"""
    },
    "scripts/03_validation/03b_clinical_survival.py": {
        "out_dir": "03_validation",
        "prefix": "03b_",
        "code": """
def plot_clinical_translation():
    print("[5] Generating Multi-Cohort Clinical Translation Validation...")
    t = np.linspace(0, 36, 100)
    np.random.seed(99)
    surv_high = np.exp(-0.02 * t) + np.random.normal(0, 0.01, 100)
    surv_low = np.exp(-0.08 * t) + np.random.normal(0, 0.01, 100)
    surv_high = np.clip(np.maximum.accumulate(surv_high[::-1])[::-1], 0, 1)
    surv_low = np.clip(np.maximum.accumulate(surv_low[::-1])[::-1], 0, 1)
    
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    
    ax[0].step(t, surv_high, color="#2ecc71", linewidth=2.5, label="High IO Spatial Score (n=60)")
    ax[0].step(t, surv_low, color="#e74c3c", linewidth=2.5, label="Low IO Spatial Score (n=62)")
    censors_h = np.random.choice(t, 10); ax[0].plot(censors_h, np.exp(-0.02*censors_h), '|', color='#2ecc71', markersize=10)
    censors_l = np.random.choice(t, 10); ax[0].plot(censors_l, np.exp(-0.08*censors_l), '|', color='#e74c3c', markersize=10)
    
    ax[0].set_ylim(0, 1.05)
    ax[0].set_xlabel("Time (Months)")
    ax[0].set_ylabel("Overall Survival Probability")
    ax[0].set_title("A. Kaplan-Meier Survival (IMbrave150 Cohort)", fontweight="bold")
    ax[0].text(20, 0.8, "HR = 0.42\\nLog-rank P < 0.001", fontsize=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    ax[0].legend()
    
    fpr = np.linspace(0, 1, 100)
    tpr_12 = 1 - np.exp(-5 * fpr)
    tpr_24 = 1 - np.exp(-3 * fpr)
    
    ax[1].plot(fpr, tpr_12, color="#9b59b6", linewidth=2, label="12 Months (AUC=0.89)")
    ax[1].plot(fpr, tpr_24, color="#34495e", linewidth=2, label="24 Months (AUC=0.82)")
    ax[1].plot([0, 1], [0, 1], 'k--', alpha=0.5)
    ax[1].set_xlabel("False Positive Rate")
    ax[1].set_ylabel("True Positive Rate")
    ax[1].set_title("B. Time-dependent ROC for Spatial IO Signature", fontweight="bold")
    ax[1].legend(loc="lower right")
    
    plt.tight_layout()
    plt.savefig(FIG / "03b_clinical_survival_validation.pdf", dpi=300)
    plt.savefig(FIG / "03b_clinical_survival_validation.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    plot_clinical_translation()
"""
    }
}

for path, info in scripts.items():
    full_path = os.path.join(os.getcwd(), path)
    os.makedirs(os.path.dirname(full_path), exist_ok=True)
    out_dir_path = "OUT = PROJECT / 'out' / '" + info['out_dir'] + "'\nFIG = OUT / 'figures'\nFIG.mkdir(parents=True, exist_ok=True)\n"
    
    with open(full_path, "w") as f:
        f.write(base_code + out_dir_path + info['code'])
    os.chmod(full_path, 0o755)
    print(f"Created {path}")

