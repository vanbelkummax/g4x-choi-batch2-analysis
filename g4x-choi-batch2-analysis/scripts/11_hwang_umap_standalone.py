#!/usr/bin/env python3
"""
Create standalone Hwang Cancer Score UMAP visualization.
"""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Setup
output_dir = Path('/home/user/g4x-choi-batch2-analysis/output/figures')
output_dir.mkdir(parents=True, exist_ok=True)

# Load annotated data
print("Loading annotated data...")
adata = sc.read_h5ad('/home/user/g4x-choi-batch2-analysis/output/pilot_final_annotated.h5ad')
print(f"Loaded {adata.n_obs:,} cells")
print(f"obsm keys: {list(adata.obsm.keys())}")

# Compute UMAP if not present
if 'X_umap' not in adata.obsm:
    print("Computing UMAP (per-sample for consistency)...")

    # Process each sample separately for better visualization
    umaps = {}
    for sample in adata.obs['sample_id'].unique():
        print(f"  Processing {sample}...")
        adata_sub = adata[adata.obs['sample_id'] == sample].copy()

        # Normalize if needed
        if 'log1p' not in adata_sub.uns:
            sc.pp.normalize_total(adata_sub, target_sum=1e4)
            sc.pp.log1p(adata_sub)

        # PCA and neighbors
        sc.pp.pca(adata_sub, n_comps=min(50, adata_sub.n_vars - 1))
        sc.pp.neighbors(adata_sub, n_neighbors=30, n_pcs=20)
        sc.tl.umap(adata_sub, random_state=42)

        # Store
        umaps[sample] = adata_sub.obsm['X_umap']

    # Combine back - keeping per-sample UMAP coordinates
    umap_combined = np.zeros((adata.n_obs, 2))
    for sample in adata.obs['sample_id'].unique():
        mask = adata.obs['sample_id'] == sample
        umap_combined[mask.values] = umaps[sample]

    adata.obsm['X_umap'] = umap_combined
    print("UMAP computed")

# Data already has 'stage' and 'hwang_cancer' columns
stage_order = ['Normal', 'Metaplasia', 'Cancer']

# Create standalone Hwang UMAP - one panel per sample
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for idx, stage in enumerate(stage_order):
    ax = axes[idx]
    mask = adata.obs['stage'] == stage
    adata_sub = adata[mask].copy()

    # Get UMAP coordinates and scores
    umap_coords = adata_sub.obsm['X_umap']
    scores = adata_sub.obs['hwang_cancer'].values

    # Calculate mean for title
    mean_score = np.mean(scores)

    # Plot
    sc_plot = ax.scatter(
        umap_coords[:, 0],
        umap_coords[:, 1],
        c=scores,
        cmap='RdBu_r',
        s=0.5,
        alpha=0.7,
        vmin=-0.3,
        vmax=0.3,
        rasterized=True
    )

    # Sample info
    sample_id = adata_sub.obs['sample_id'].iloc[0]
    n_cells = mask.sum()

    ax.set_title(f'{stage} ({sample_id})\nn={n_cells:,} | Mean: {mean_score:.3f}', fontsize=12)
    ax.set_xlabel('UMAP1', fontsize=10)
    ax.set_ylabel('UMAP2' if idx == 0 else '', fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Add colorbar
cbar_ax = fig.add_axes([0.92, 0.2, 0.02, 0.6])
cbar = plt.colorbar(sc_plot, cax=cbar_ax)
cbar.set_label('Hwang Cancer Score', fontsize=11)

plt.suptitle('Hwang 32-Gene Cancer Signature Score (UMAP)', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout(rect=[0, 0, 0.9, 1])

# Save
fig.savefig(output_dir / 'fig5_hwang_umap.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(output_dir / 'fig5_hwang_umap.pdf', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {output_dir / 'fig5_hwang_umap.png'}")

# Also create a combined spatial + UMAP version for Hwang score
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Top row: Spatial
for idx, stage in enumerate(stage_order):
    ax = axes[0, idx]
    mask = adata.obs['stage'] == stage
    adata_sub = adata[mask].copy()

    spatial_coords = adata_sub.obsm['spatial']
    scores = adata_sub.obs['hwang_cancer'].values

    sc1 = ax.scatter(
        spatial_coords[:, 0],
        spatial_coords[:, 1],
        c=scores,
        cmap='RdBu_r',
        s=0.3,
        alpha=0.7,
        vmin=-0.3,
        vmax=0.3,
        rasterized=True
    )

    sample_id = adata_sub.obs['sample_id'].iloc[0]
    n_cells = mask.sum()
    mean_score = np.mean(scores)

    ax.set_title(f'{stage} ({sample_id})\nn={n_cells:,} | Mean: {mean_score:.3f}', fontsize=11)
    ax.set_xlabel('')
    ax.set_ylabel('Spatial' if idx == 0 else '')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    ax.invert_yaxis()

# Bottom row: UMAP
for idx, stage in enumerate(stage_order):
    ax = axes[1, idx]
    mask = adata.obs['stage'] == stage
    adata_sub = adata[mask].copy()

    umap_coords = adata_sub.obsm['X_umap']
    scores = adata_sub.obs['hwang_cancer'].values

    sc2 = ax.scatter(
        umap_coords[:, 0],
        umap_coords[:, 1],
        c=scores,
        cmap='RdBu_r',
        s=0.5,
        alpha=0.7,
        vmin=-0.3,
        vmax=0.3,
        rasterized=True
    )

    ax.set_xlabel('')
    ax.set_ylabel('UMAP' if idx == 0 else '')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Add colorbar
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
cbar = plt.colorbar(sc2, cax=cbar_ax)
cbar.set_label('Hwang Cancer Score', fontsize=11)

plt.suptitle('Hwang 32-Gene Cancer Signature: Spatial & UMAP', fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 0.9, 0.95])

fig.savefig(output_dir / 'fig5_hwang_spatial_umap.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(output_dir / 'fig5_hwang_spatial_umap.pdf', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {output_dir / 'fig5_hwang_spatial_umap.png'}")

# Print summary stats
print("\n=== Hwang Cancer Score Summary ===")
for stage in stage_order:
    mask = adata.obs['stage'] == stage
    scores = adata[mask].obs['hwang_cancer']
    print(f"{stage}: Mean={scores.mean():.4f}, Std={scores.std():.4f}, "
          f"Min={scores.min():.4f}, Max={scores.max():.4f}")

print("\nDone!")
