#!/usr/bin/env python3
"""
01_hvg_optimize.py - HVG optimization with grid search

Grid search over HVG counts and clustering resolutions.
Evaluates using silhouette score and Calinski-Harabasz index.
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# === CONFIG ===
INPUT = Path('/home/user/g4x-choi-batch2-analysis/results/pilot/merged_pilot.h5ad')
OUTPUT = Path('/home/user/g4x-choi-batch2-analysis/results/pilot')
FIGURES = OUTPUT / 'figures'
FIGURES.mkdir(parents=True, exist_ok=True)

SAMPLES = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}

# HVG grid: includes 'all' for using all 337 genes
HVG_GRID = [30, 50, 100, 200, 'all']

# Resolution grid for Leiden clustering
RESOLUTION_GRID = [0.3, 0.5, 0.8, 1.0]

# GPU UMAP if available
try:
    from cuml.manifold import UMAP as cumlUMAP
    import cupy as cp
    def run_umap(X):
        return cp.asnumpy(cumlUMAP(n_neighbors=15, min_dist=0.3, verbose=False)
                         .fit_transform(cp.asarray(X.astype(np.float32))))
    print("[GPU] cuML UMAP enabled")
except ImportError:
    import umap
    def run_umap(X):
        return umap.UMAP(n_neighbors=15, min_dist=0.3).fit_transform(X)
    print("[CPU] Using standard UMAP")


def process_hvg_config(adata: sc.AnnData, n_hvg, resolution: float) -> dict:
    """Process one HVG + resolution configuration."""
    ad = adata.copy()

    # Select HVGs or use all genes
    if n_hvg != 'all':
        n_hvg = min(n_hvg, ad.n_vars - 1)  # Can't select more HVGs than genes
        sc.pp.highly_variable_genes(
            ad, n_top_genes=n_hvg, flavor='seurat_v3',
            layer='counts', span=0.3
        )
        ad = ad[:, ad.var.highly_variable].copy()

    n_genes_used = ad.n_vars

    # PCA
    n_pcs = min(30, n_genes_used - 1)
    sc.pp.pca(ad, n_comps=n_pcs)

    # Neighbors + clustering
    sc.pp.neighbors(ad, n_neighbors=15, n_pcs=n_pcs)
    sc.tl.leiden(ad, resolution=resolution)

    # Calculate metrics
    n_clusters = ad.obs['leiden'].nunique()

    if n_clusters > 1:
        # Subsample for fast silhouette (O(n²) is too slow for 30K+ cells)
        n_sample = min(5000, ad.n_obs)
        idx = np.random.choice(ad.n_obs, n_sample, replace=False)
        sil = silhouette_score(ad.obsm['X_pca'][idx], ad.obs['leiden'].iloc[idx])
        ch = calinski_harabasz_score(ad.obsm['X_pca'][idx], ad.obs['leiden'].iloc[idx])
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
        'adata': ad  # Keep for best config
    }


def main():
    print("="*60)
    print("HVG Optimization with Grid Search")
    print("="*60)

    # Load data
    print(f"\nLoading: {INPUT}")
    adata = sc.read_h5ad(INPUT)
    print(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars} genes")

    all_results = []
    best_configs = {}

    for sample_id, stage in SAMPLES.items():
        print(f"\n{'='*50}")
        print(f"{sample_id} ({stage})")
        print(f"{'='*50}")

        # Subset to sample
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

        # Find best config (highest silhouette)
        best = max(sample_results, key=lambda x: x['silhouette'])
        best_configs[sample_id] = best
        print(f"\n  BEST: HVG={best['n_hvg']}, res={best['resolution']:.1f}, "
              f"sil={best['silhouette']:.3f}")

        all_results.extend(sample_results)

    # Create results DataFrame (drop adata column)
    results_df = pd.DataFrame([{k: v for k, v in r.items() if k != 'adata'}
                               for r in all_results])

    # Save results
    results_path = OUTPUT / 'hvg_optimization.csv'
    results_df.to_csv(results_path, index=False)
    print(f"\nSaved: {results_path}")

    # === Create Visualization ===
    print("\nCreating visualization...")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, (sample_id, stage) in enumerate(SAMPLES.items()):
        ax = axes[i]
        sample_df = results_df[results_df['sample_id'] == sample_id]

        # Pivot for heatmap
        pivot = sample_df.pivot(index='n_hvg', columns='resolution', values='silhouette')
        # Reorder rows
        row_order = ['30', '50', '100', '200', 'all']
        pivot = pivot.reindex([r for r in row_order if r in pivot.index])

        sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlGn',
                    ax=ax, cbar_kws={'label': 'Silhouette'})
        ax.set_title(f'{sample_id}: {stage}', fontweight='bold')
        ax.set_xlabel('Resolution')
        ax.set_ylabel('HVG')

    plt.suptitle('HVG Optimization: Silhouette Score by Config', fontsize=14, fontweight='bold')
    plt.tight_layout()

    fig_path = FIGURES / 'hvg_silhouette_heatmap.png'
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_path}")

    # === Summary Table ===
    print("\n" + "="*60)
    print("SUMMARY: Best Configurations")
    print("="*60)
    print(f"{'Sample':<8} {'Stage':<12} {'HVG':>6} {'Res':>5} {'Clusters':>8} {'Silhouette':>10}")
    print("-"*60)
    for sample_id, cfg in best_configs.items():
        print(f"{sample_id:<8} {cfg['stage']:<12} {cfg['n_hvg']:>6} "
              f"{cfg['resolution']:>5.1f} {cfg['n_clusters']:>8} {cfg['silhouette']:>10.3f}")

    # Global recommendation
    print("\n" + "-"*60)
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


if __name__ == '__main__':
    main()
