#!/usr/bin/env python
"""Generate UMAP panel with optimal clustering - GPU accelerated."""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# Try GPU acceleration
try:
    import rapids_singlecell as rsc
    USE_GPU = True
    print("Using GPU (rapids-singlecell)")
except ImportError:
    USE_GPU = False
    print("Using CPU (scanpy)")

# Settings
sc.settings.verbosity = 1
plt.rcParams['figure.dpi'] = 150

# Sample info
samples = ['E02', 'F02', 'G02']
stages = ['Normal', 'Metaplasia', 'Cancer']

# Optimal params
RESOLUTION = 0.3

# Load and process each sample
adatas = []
for sample in samples:
    print(f"\nProcessing {sample}...")
    adata = sc.read_h5ad(f'/home/user/g4x-choi-batch2-analysis/results/pilot/{sample}_annotated.h5ad')
    print(f"  {adata.n_obs} cells, {adata.n_vars} genes")

    if USE_GPU:
        # GPU pipeline
        rsc.get.anndata_to_GPU(adata)
        rsc.pp.pca(adata, n_comps=30)
        rsc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
        rsc.tl.umap(adata)
        rsc.tl.leiden(adata, resolution=RESOLUTION, key_added='leiden_optimal')
        rsc.get.anndata_to_CPU(adata)
    else:
        # CPU pipeline
        sc.tl.pca(adata, n_comps=30)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=RESOLUTION, key_added='leiden_optimal')

    n_clusters = adata.obs['leiden_optimal'].nunique()
    print(f"  {n_clusters} clusters at res={RESOLUTION}")
    adatas.append(adata)

# Create panel figure
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Top row: UMAPs colored by cluster
for idx, (adata, sample, stage) in enumerate(zip(adatas, samples, stages)):
    ax = axes[0, idx]
    sc.pl.umap(adata, color='leiden_optimal', ax=ax, show=False,
               title=f'{sample} - {stage}\n(n={adata.n_obs:,} cells)',
               legend_loc='right margin' if idx == 2 else 'none',
               frameon=True)
    n_clusters = adata.obs['leiden_optimal'].nunique()
    ax.text(0.02, 0.98, f'{n_clusters} clusters', transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Bottom row: Total counts
for idx, (adata, sample, stage) in enumerate(zip(adatas, samples, stages)):
    ax = axes[1, idx]
    if 'total_counts' not in adata.obs:
        adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()
    sc.pl.umap(adata, color='total_counts', ax=ax, show=False,
               title=f'Expression Intensity', cmap='viridis', frameon=True,
               colorbar_loc='right' if idx == 2 else None)

fig.suptitle('UMAP Embeddings with Optimal Clustering (res=0.3, all 337 genes)\n'
             'Top: Leiden clusters | Bottom: Total expression',
             fontsize=14, fontweight='bold', y=1.02)

plt.tight_layout()

# Save
out_path = '/home/user/g4x-choi-batch2-analysis/output/figures/umap_panel_optimal.png'
plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {out_path}")

desktop_path = '/mnt/c/Users/User/Desktop/G4X_Analysis/figures/umap_panel_optimal.png'
plt.savefig(desktop_path, dpi=150, bbox_inches='tight', facecolor='white')
print(f"Saved: {desktop_path}")

print("\n" + "="*60)
print("CLUSTER SUMMARY (res=0.3)")
print("="*60)
for adata, sample, stage in zip(adatas, samples, stages):
    n_clusters = adata.obs['leiden_optimal'].nunique()
    print(f"{sample} ({stage}): {n_clusters} clusters, {adata.n_obs:,} cells")
print("="*60)
