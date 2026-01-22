#!/usr/bin/env python3
"""
01_hvg_optimize_parallel.py - FULLY PARALLELIZED HVG optimization

Uses:
- GPU: PCA + silhouette (PyTorch)
- Multiprocessing: Parallel samples AND configs
- All 24 CPU cores
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import torch
from sklearn.metrics import calinski_harabasz_score
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
from functools import partial
import time

warnings.filterwarnings('ignore')

# Force spawn for CUDA compatibility
mp.set_start_method('spawn', force=True)

# Config
INPUT = Path('/home/user/g4x-choi-batch2-analysis/results/pilot/merged_pilot.h5ad')
OUTPUT = Path('/home/user/g4x-choi-batch2-analysis/results/pilot')
FIGURES = OUTPUT / 'figures'
FIGURES.mkdir(parents=True, exist_ok=True)

SAMPLES = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}
HVG_GRID = [30, 50, 100, 200, 'all']
RESOLUTION_GRID = [0.3, 0.5, 0.8, 1.0]
N_WORKERS = 20  # Use 20 of 24 cores


def gpu_pca(X, n_components=30):
    """GPU PCA using PyTorch SVD."""
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    X_tensor = torch.tensor(X, dtype=torch.float32, device=device)
    X_centered = X_tensor - X_tensor.mean(dim=0)
    U, S, V = torch.svd(X_centered)
    return (U[:, :n_components] * S[:n_components]).cpu().numpy()


def fast_silhouette(X_pca, labels, n_sample=5000):
    """Fast centroid-based silhouette approximation."""
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    n = len(labels)
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)

    if n_clusters < 2:
        return 0.0

    if n > n_sample:
        idx = np.random.choice(n, n_sample, replace=False)
        X_pca = X_pca[idx]
        labels = labels[idx]
        n = n_sample

    X_tensor = torch.tensor(X_pca, dtype=torch.float32, device=device)

    # Cluster centroids
    centroids = []
    for lab in unique_labels:
        mask = labels == lab
        centroids.append(X_tensor[torch.tensor(mask, device=device)].mean(dim=0))
    centroids = torch.stack(centroids)

    # Distances to centroids
    dists = torch.cdist(X_tensor, centroids, p=2)

    label_to_idx = {l: i for i, l in enumerate(unique_labels)}
    label_indices = torch.tensor([label_to_idx[l] for l in labels], device=device)

    # Own cluster distance
    a = dists[torch.arange(n, device=device), label_indices]

    # Nearest other cluster
    dists_masked = dists.clone()
    dists_masked[torch.arange(n, device=device), label_indices] = float('inf')
    b = dists_masked.min(dim=1).values

    silhouette = (b - a) / torch.maximum(a, b)
    silhouette[torch.isnan(silhouette)] = 0

    return silhouette.mean().item()


def process_config(args):
    """Process single HVG+resolution config. Run in subprocess."""
    sample_data, sample_id, stage, n_hvg, resolution = args

    # Reconstruct AnnData in subprocess
    ad = sc.AnnData(X=sample_data['X'], obs=sample_data['obs'], var=sample_data['var'])
    if 'counts' in sample_data:
        ad.layers['counts'] = sample_data['counts']

    # HVG selection
    if n_hvg != 'all':
        n_hvg_use = min(n_hvg, ad.n_vars - 1)
        sc.pp.highly_variable_genes(
            ad, n_top_genes=n_hvg_use, flavor='seurat_v3',
            layer='counts', span=0.3
        )
        ad = ad[:, ad.var.highly_variable].copy()

    n_genes = ad.n_vars

    # PCA (GPU)
    n_pcs = min(30, n_genes - 1)
    X = ad.X.toarray() if hasattr(ad.X, 'toarray') else ad.X
    ad.obsm['X_pca'] = gpu_pca(X, n_components=n_pcs)

    # Neighbors + Leiden (CPU - use threads)
    sc.pp.neighbors(ad, n_neighbors=15, n_pcs=n_pcs, use_rep='X_pca')
    sc.tl.leiden(ad, resolution=resolution)

    n_clusters = ad.obs['leiden'].nunique()

    if n_clusters > 1:
        sil = fast_silhouette(ad.obsm['X_pca'], ad.obs['leiden'].values)
        ch = calinski_harabasz_score(ad.obsm['X_pca'], ad.obs['leiden'])
    else:
        sil, ch = 0.0, 0.0

    return {
        'sample_id': sample_id,
        'stage': stage,
        'n_hvg': str(n_hvg),
        'n_genes_used': n_genes,
        'resolution': resolution,
        'n_clusters': n_clusters,
        'silhouette': sil,
        'calinski_harabasz': ch
    }


def prepare_sample_data(adata, sample_id):
    """Prepare sample data for multiprocessing (pickle-able)."""
    sample = adata[adata.obs['sample_id'] == sample_id].copy()
    return {
        'X': sample.X.toarray() if hasattr(sample.X, 'toarray') else sample.X,
        'obs': sample.obs.copy(),
        'var': sample.var.copy(),
        'counts': sample.layers.get('counts', sample.X).toarray()
                  if hasattr(sample.layers.get('counts', sample.X), 'toarray')
                  else sample.layers.get('counts', sample.X)
    }


def main():
    start_time = time.time()

    print("=" * 60)
    print("PARALLEL HVG Optimization (20 workers)")
    print("=" * 60)

    # Check GPU
    if torch.cuda.is_available():
        print(f"GPU: {torch.cuda.get_device_name(0)}")

    # Load data
    print(f"\nLoading: {INPUT}")
    adata = sc.read_h5ad(INPUT)
    print(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars} genes")

    # Prepare all tasks
    tasks = []
    for sample_id, stage in SAMPLES.items():
        sample_data = prepare_sample_data(adata, sample_id)
        n_cells = sample_data['X'].shape[0]
        print(f"  {sample_id} ({stage}): {n_cells:,} cells")

        for n_hvg in HVG_GRID:
            for res in RESOLUTION_GRID:
                tasks.append((sample_data, sample_id, stage, n_hvg, res))

    print(f"\nTotal configs: {len(tasks)} (running {N_WORKERS} in parallel)")

    # Run in parallel
    results = []
    completed = 0

    with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        futures = {executor.submit(process_config, task): task for task in tasks}

        for future in as_completed(futures):
            try:
                result = future.result()
                results.append(result)
                completed += 1

                print(f"  [{completed:2d}/{len(tasks)}] {result['sample_id']} "
                      f"HVG={result['n_hvg']:>4} res={result['resolution']:.1f}: "
                      f"{result['n_clusters']:2d} clusters, sil={result['silhouette']:.3f}")

            except Exception as e:
                print(f"  ERROR: {e}")

    elapsed = time.time() - start_time
    print(f"\nCompleted in {elapsed:.1f}s ({len(tasks)/elapsed:.1f} configs/sec)")

    # Create DataFrame
    results_df = pd.DataFrame(results)
    results_df.to_csv(OUTPUT / 'hvg_optimization.csv', index=False)
    print(f"Saved: {OUTPUT / 'hvg_optimization.csv'}")

    # Best configs per sample
    print("\n" + "=" * 60)
    print("BEST CONFIGURATIONS")
    print("=" * 60)

    best_configs = {}
    for sample_id in SAMPLES:
        sample_df = results_df[results_df['sample_id'] == sample_id]
        best = sample_df.loc[sample_df['silhouette'].idxmax()]
        best_configs[sample_id] = best
        print(f"{sample_id}: HVG={best['n_hvg']}, res={best['resolution']}, "
              f"sil={best['silhouette']:.3f}, clusters={best['n_clusters']}")

    # Global recommendation
    avg_by_hvg = results_df.groupby('n_hvg')['silhouette'].mean()
    best_hvg = avg_by_hvg.idxmax()
    print(f"\nRECOMMENDATION: HVG={best_hvg} (avg sil={avg_by_hvg[best_hvg]:.3f})")

    # Save recommendation
    with open(OUTPUT / 'hvg_recommendation.txt', 'w') as f:
        f.write(f"Best HVG: {best_hvg}\n")
        f.write(f"Avg silhouette: {avg_by_hvg[best_hvg]:.3f}\n\n")
        for sample_id, cfg in best_configs.items():
            f.write(f"{sample_id}: HVG={cfg['n_hvg']}, res={cfg['resolution']}, "
                    f"sil={cfg['silhouette']:.3f}\n")

    # Visualization
    print("\nCreating heatmap...")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for i, (sample_id, stage) in enumerate(SAMPLES.items()):
        sample_df = results_df[results_df['sample_id'] == sample_id]
        pivot = sample_df.pivot(index='n_hvg', columns='resolution', values='silhouette')
        row_order = ['30', '50', '100', '200', 'all']
        pivot = pivot.reindex([r for r in row_order if r in pivot.index])

        sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlGn',
                    ax=axes[i], cbar_kws={'label': 'Silhouette'})
        axes[i].set_title(f'{sample_id}: {stage}', fontweight='bold')

    plt.suptitle(f'HVG Optimization (completed in {elapsed:.0f}s)', fontweight='bold')
    plt.tight_layout()
    fig.savefig(FIGURES / 'hvg_silhouette_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {FIGURES / 'hvg_silhouette_heatmap.png'}")


if __name__ == '__main__':
    main()
