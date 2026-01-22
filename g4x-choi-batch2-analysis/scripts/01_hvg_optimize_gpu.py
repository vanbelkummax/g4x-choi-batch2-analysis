#!/usr/bin/env python3
"""
01_hvg_optimize_gpu.py - GPU-accelerated HVG optimization

Uses PyTorch GPU for PCA and silhouette computation on full samples.
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import torch
import torch.nn.functional as F
from sklearn.metrics import calinski_harabasz_score
warnings.filterwarnings('ignore')

# Check GPU
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"[{'GPU' if device.type == 'cuda' else 'CPU'}] Using {device}")
if device.type == 'cuda':
    print(f"  Device: {torch.cuda.get_device_name(0)}")
    print(f"  Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

# === CONFIG ===
INPUT = Path('/home/user/g4x-choi-batch2-analysis/results/pilot/merged_pilot.h5ad')
OUTPUT = Path('/home/user/g4x-choi-batch2-analysis/results/pilot')
FIGURES = OUTPUT / 'figures'
FIGURES.mkdir(parents=True, exist_ok=True)

SAMPLES = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}
HVG_GRID = [30, 50, 100, 200, 'all']
RESOLUTION_GRID = [0.3, 0.5, 0.8, 1.0]


def gpu_pca(X, n_components=30):
    """GPU-accelerated PCA using PyTorch."""
    X_tensor = torch.tensor(X, dtype=torch.float32, device=device)
    # Center the data
    X_centered = X_tensor - X_tensor.mean(dim=0)
    # SVD
    U, S, V = torch.svd(X_centered)
    # Project to n_components
    return (U[:, :n_components] * S[:n_components]).cpu().numpy()


def gpu_silhouette_sample(X_pca, labels, n_sample=10000):
    """
    GPU-accelerated silhouette score with sampling.
    Uses batched distance computation for memory efficiency.
    """
    n = len(labels)
    if n > n_sample:
        idx = np.random.choice(n, n_sample, replace=False)
        X_pca = X_pca[idx]
        labels = labels[idx]
        n = n_sample

    X_tensor = torch.tensor(X_pca, dtype=torch.float32, device=device)
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)

    if n_clusters < 2:
        return 0.0

    # Create cluster assignment matrix
    label_to_idx = {l: i for i, l in enumerate(unique_labels)}
    label_indices = np.array([label_to_idx[l] for l in labels])

    # Compute silhouette in batches to avoid OOM
    batch_size = 5000
    silhouettes = []

    for i in range(0, n, batch_size):
        end = min(i + batch_size, n)
        batch = X_tensor[i:end]
        batch_labels = label_indices[i:end]

        # Compute distances from batch to all points
        # Using squared Euclidean for speed (same ranking)
        dists = torch.cdist(batch, X_tensor, p=2)  # [batch, n]

        # Compute mean intra-cluster distance (a) and min inter-cluster distance (b)
        batch_sil = []
        for j, (drow, lab) in enumerate(zip(dists, batch_labels)):
            same_cluster = label_indices == lab
            other_clusters = ~same_cluster

            # a(i) = mean distance to same cluster (excluding self)
            same_mask = torch.tensor(same_cluster, device=device)
            same_mask[i + j] = False  # Exclude self
            if same_mask.sum() > 0:
                a = drow[same_mask].mean().item()
            else:
                a = 0

            # b(i) = min mean distance to other clusters
            b = float('inf')
            for k, other_lab in enumerate(unique_labels):
                if k != lab:
                    other_mask = label_indices == k
                    if other_mask.sum() > 0:
                        mean_dist = drow[torch.tensor(other_mask, device=device)].mean().item()
                        b = min(b, mean_dist)

            if b == float('inf'):
                b = 0

            # Silhouette for this point
            if max(a, b) > 0:
                s = (b - a) / max(a, b)
            else:
                s = 0
            batch_sil.append(s)

        silhouettes.extend(batch_sil)

    return np.mean(silhouettes)


def gpu_silhouette_fast(X_pca, labels, n_sample=5000):
    """
    Fast GPU silhouette using cluster centroids approximation.
    Much faster for large datasets.
    """
    n = len(labels)
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)

    if n_clusters < 2:
        return 0.0

    # Sample if needed
    if n > n_sample:
        idx = np.random.choice(n, n_sample, replace=False)
        X_pca = X_pca[idx]
        labels = labels[idx]
        n = n_sample

    X_tensor = torch.tensor(X_pca, dtype=torch.float32, device=device)

    # Compute cluster centroids
    centroids = []
    for lab in unique_labels:
        mask = labels == lab
        centroids.append(X_tensor[torch.tensor(mask, device=device)].mean(dim=0))
    centroids = torch.stack(centroids)  # [n_clusters, n_features]

    # Distance from each point to all centroids
    dists_to_centroids = torch.cdist(X_tensor, centroids, p=2)  # [n, n_clusters]

    # For each point: a = dist to own centroid, b = min dist to other centroids
    label_to_idx = {l: i for i, l in enumerate(unique_labels)}
    label_indices = torch.tensor([label_to_idx[l] for l in labels], device=device)

    # Get own cluster distance
    a = dists_to_centroids[torch.arange(n, device=device), label_indices]

    # Get nearest other cluster distance
    # Set own cluster to inf, then take min
    dists_masked = dists_to_centroids.clone()
    dists_masked[torch.arange(n, device=device), label_indices] = float('inf')
    b = dists_masked.min(dim=1).values

    # Compute silhouette
    silhouette = (b - a) / torch.maximum(a, b)
    silhouette[torch.isnan(silhouette)] = 0

    return silhouette.mean().item()


def process_hvg_config(adata: sc.AnnData, n_hvg, resolution: float) -> dict:
    """Process one HVG + resolution configuration with GPU acceleration."""
    ad = adata.copy()

    # Select HVGs or use all genes
    if n_hvg != 'all':
        n_hvg = min(n_hvg, ad.n_vars - 1)
        sc.pp.highly_variable_genes(
            ad, n_top_genes=n_hvg, flavor='seurat_v3',
            layer='counts', span=0.3
        )
        ad = ad[:, ad.var.highly_variable].copy()

    n_genes_used = ad.n_vars

    # GPU PCA
    n_pcs = min(30, n_genes_used - 1)
    X = ad.X.toarray() if hasattr(ad.X, 'toarray') else ad.X
    ad.obsm['X_pca'] = gpu_pca(X, n_components=n_pcs)

    # Neighbors + clustering (CPU - fast enough)
    sc.pp.neighbors(ad, n_neighbors=15, n_pcs=n_pcs, use_rep='X_pca')
    sc.tl.leiden(ad, resolution=resolution)

    # Metrics
    n_clusters = ad.obs['leiden'].nunique()

    if n_clusters > 1:
        # GPU silhouette (fast centroid-based approximation)
        sil = gpu_silhouette_fast(ad.obsm['X_pca'], ad.obs['leiden'].values)
        # CPU Calinski-Harabasz (fast enough)
        ch = calinski_harabasz_score(ad.obsm['X_pca'], ad.obs['leiden'])
    else:
        sil = 0.0
        ch = 0.0

    return {
        'n_hvg': str(n_hvg),
        'n_genes_used': n_genes_used,
        'resolution': resolution,
        'n_clusters': n_clusters,
        'silhouette': sil,
        'calinski_harabasz': ch,
    }


def main():
    print("=" * 60)
    print("HVG Optimization with GPU Acceleration")
    print("=" * 60)

    # Load data
    print(f"\nLoading: {INPUT}")
    adata = sc.read_h5ad(INPUT)
    print(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars} genes")

    all_results = []
    best_configs = {}

    for sample_id, stage in SAMPLES.items():
        print(f"\n{'=' * 50}")
        print(f"{sample_id} ({stage})")
        print(f"{'=' * 50}")

        sample = adata[adata.obs['sample_id'] == sample_id].copy()
        print(f"Cells: {sample.n_obs:,}")

        sample_results = []

        for n_hvg in HVG_GRID:
            for res in RESOLUTION_GRID:
                result = process_hvg_config(sample, n_hvg, res)
                result['sample_id'] = sample_id
                result['stage'] = stage
                sample_results.append(result)

                print(f"  HVG={n_hvg:>4}, res={res:.1f}: "
                      f"{result['n_clusters']:2d} clusters, "
                      f"sil={result['silhouette']:.3f}, "
                      f"CH={result['calinski_harabasz']:.0f}")

        # Find best config
        best = max(sample_results, key=lambda x: x['silhouette'])
        best_configs[sample_id] = best
        print(f"\n  BEST: HVG={best['n_hvg']}, res={best['resolution']:.1f}, "
              f"sil={best['silhouette']:.3f}")

        all_results.extend(sample_results)

    # Save results
    results_df = pd.DataFrame(all_results)
    results_path = OUTPUT / 'hvg_optimization.csv'
    results_df.to_csv(results_path, index=False)
    print(f"\nSaved: {results_path}")

    # Visualization
    print("\nCreating visualization...")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, (sample_id, stage) in enumerate(SAMPLES.items()):
        ax = axes[i]
        sample_df = results_df[results_df['sample_id'] == sample_id]
        pivot = sample_df.pivot(index='n_hvg', columns='resolution', values='silhouette')
        row_order = ['30', '50', '100', '200', 'all']
        pivot = pivot.reindex([r for r in row_order if r in pivot.index])

        sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlGn',
                    ax=ax, cbar_kws={'label': 'Silhouette'})
        ax.set_title(f'{sample_id}: {stage}', fontweight='bold')
        ax.set_xlabel('Resolution')
        ax.set_ylabel('HVG')

    plt.suptitle('HVG Optimization: Silhouette Score (GPU-Accelerated)', fontsize=14, fontweight='bold')
    plt.tight_layout()

    fig_path = FIGURES / 'hvg_silhouette_heatmap.png'
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_path}")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY: Best Configurations")
    print("=" * 60)
    print(f"{'Sample':<8} {'Stage':<12} {'HVG':>6} {'Res':>5} {'Clusters':>8} {'Silhouette':>10}")
    print("-" * 60)
    for sample_id, cfg in best_configs.items():
        print(f"{sample_id:<8} {cfg['stage']:<12} {cfg['n_hvg']:>6} "
              f"{cfg['resolution']:>5.1f} {cfg['n_clusters']:>8} {cfg['silhouette']:>10.3f}")

    # Recommendation
    print("\n" + "-" * 60)
    avg_by_hvg = results_df.groupby('n_hvg')['silhouette'].mean()
    best_hvg = avg_by_hvg.idxmax()
    print(f"RECOMMENDATION: Use HVG={best_hvg} (avg silhouette={avg_by_hvg[best_hvg]:.3f})")

    # Save recommendation
    rec_path = OUTPUT / 'hvg_recommendation.txt'
    with open(rec_path, 'w') as f:
        f.write(f"Best HVG setting: {best_hvg}\n")
        f.write(f"Average silhouette across samples: {avg_by_hvg[best_hvg]:.3f}\n")
        f.write("\nPer-sample optimal configs:\n")
        for sample_id, cfg in best_configs.items():
            f.write(f"  {sample_id}: HVG={cfg['n_hvg']}, res={cfg['resolution']}, "
                    f"sil={cfg['silhouette']:.3f}\n")
    print(f"Saved: {rec_path}")

    # GPU memory summary
    if device.type == 'cuda':
        print(f"\nGPU Memory Used: {torch.cuda.max_memory_allocated() / 1e9:.2f} GB")


if __name__ == '__main__':
    main()
