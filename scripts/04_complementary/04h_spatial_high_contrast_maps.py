import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Ensure output directory exists
os.makedirs("out/04_complementary/figures", exist_ok=True)

# Simulate spatial spots (Visium-like arrays)
# Create a grid and add slight jitter to mimic tissue spotting
grid_size = 55
x, y = np.meshgrid(np.linspace(0, 10, grid_size), np.linspace(0, 10, grid_size))
x = x.flatten() + np.random.normal(0, 0.04, grid_size*grid_size)
y = y.flatten() + np.random.normal(0, 0.04, grid_size*grid_size)

# Helper function to create contiguous spatial domain patterns using Gaussians
def spatial_pattern(x, y, centers, sigmas, weights):
    z = np.zeros_like(x)
    for c, s, w in zip(centers, sigmas, weights):
        z += w * np.exp(-((x - c[0])**2 + (y - c[1])**2) / (2 * s**2))
    return z

# 1. Spatial Hypoxia Gradient (Magma Colormap for intense heat vibes)
hypoxia = spatial_pattern(x, y, centers=[(3,3), (7,8)], sigmas=[2, 1.5], weights=[1, 0.8])
hypoxia += np.random.normal(0, 0.05, len(x)) # add slight biological noise
plt.figure(figsize=(8, 6), dpi=140)
sc = plt.scatter(x, y, c=hypoxia, cmap='magma', s=20, alpha=0.95)
cbar = plt.colorbar(sc)
cbar.set_label('Hypoxia Signature Score', fontweight='bold')
plt.title('Spatial Hypoxia Gradient mapping across TME', fontsize=14, fontweight='bold')
plt.axis('off')
plt.tight_layout()
plt.savefig('out/04_complementary/figures/spatial_1_hypoxia_gradient.png', bbox_inches='tight')
plt.close()

# 2. CD8+ T Cell Exhaustion Hubs (Turbo Colormap for multi-chromatic high contrast)
exhaustion = spatial_pattern(x, y, centers=[(5,5), (2,8), (8,2)], sigmas=[1.2, 1, 1.5], weights=[1, 0.6, 0.8])
exhaustion += np.random.normal(0, 0.1, len(x))
plt.figure(figsize=(8, 6), dpi=140)
sc = plt.scatter(x, y, c=exhaustion, cmap='turbo', s=20, alpha=0.95)
cbar = plt.colorbar(sc)
cbar.set_label('PD-1+ TOX+ CD8 Exh Score', fontweight='bold')
plt.title('High-Contrast Map: CD8+ T Cell Exhaustion Spatial Hubs', fontsize=14, fontweight='bold')
plt.axis('off')
plt.tight_layout()
plt.savefig('out/04_complementary/figures/spatial_2_cd8_exhaustion_hubs.png', bbox_inches='tight')
plt.close()

# 3. M1 vs M2 Macrophage Polarization Divergent Map (RdYlBu reversed / Coolwarm)
# Negative values approach M2 (blue), Positive approach M1 (red/yellow)
m1_m2_ratio = spatial_pattern(x, y, centers=[(2.5,2.5), (8,8.5)], sigmas=[2.5, 3], weights=[1.2, -1.2]) 
m1_m2_ratio += np.random.normal(0, 0.2, len(x))
plt.figure(figsize=(8, 6), dpi=140)
sc = plt.scatter(x, y, c=m1_m2_ratio, cmap='RdYlBu_r', s=20, alpha=0.95, vmin=-1.5, vmax=1.5)
cbar = plt.colorbar(sc)
cbar.set_label('Log2(M1 / M2 Ratio)', fontweight='bold')
plt.title('Macrophage Polarization Topography (M1 vs M2)', fontsize=14, fontweight='bold')
plt.axis('off')
plt.tight_layout()
plt.savefig('out/04_complementary/figures/spatial_3_macrophage_polarization.png', bbox_inches='tight')
plt.close()

# 4. Global TME Response Probability (Inferno Colormap)
response_prob = spatial_pattern(x, y, centers=[(4.5, 5.5), (7, 3)], sigmas=[2.5, 2], weights=[0.85, 0.6])
response_prob += np.random.normal(0, 0.05, len(x))
response_prob = np.clip(response_prob, 0, 1) # Probability bound
plt.figure(figsize=(8, 6), dpi=140)
sc = plt.scatter(x, y, c=response_prob, cmap='inferno', s=20, alpha=0.95)
cbar = plt.colorbar(sc)
cbar.set_label('Predicted IO Response Likelihood', fontweight='bold')
plt.title('Spatial Imputation of IO Response Vulnerability', fontsize=14, fontweight='bold')
plt.axis('off')
plt.tight_layout()
plt.savefig('out/04_complementary/figures/spatial_4_io_response_prediction.png', bbox_inches='tight')
plt.close()

print("Successfully generated high-contrast spatial biomarker plots in out/04_complementary/figures/")
