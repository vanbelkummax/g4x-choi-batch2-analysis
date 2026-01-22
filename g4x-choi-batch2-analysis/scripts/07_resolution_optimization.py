#!/usr/bin/env python3
"""
07_resolution_optimization.py - Multi-resolution clustering optimization

Tests multiple Leiden resolutions and computes biology-informed metrics:
- NMI: Normalized Mutual Information between clusters and reference labels
- ARI: Adjusted Rand Index between clusters and reference labels
- Purity: Mean proportion of dominant cell type per cluster
- Silhouette: Cluster separation in PCA space

Finds the optimal resolution that best recovers biological cell types.
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.metrics import (
    normalized_mutual_info_score,
    adjusted_rand_score,
    silhouette_score
)
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# === CONFIG ===
BASE = Path('/home/user/g4x-choi-batch2-analysis')
RESULTS = BASE / 'results/pilot'
FIGURES = RESULTS / 'figures'
FIGURES.mkdir(parents=True, exist_ok=True)

SAMPLES = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}

# Resolution range to test
RESOLUTIONS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0]

# Label column priority for computing metrics
LABEL_COLS = ['celltype_reference', 'celltype_markers', 'celltype_gating']

# Composite score weights
WEIGHTS = {
    'nmi': 0.35,      # Agreement with reference labels
    'ari': 0.25,      # Cluster-label correspondence
    'purity': 0.25,   # Cluster homogeneity
    'silhouette': 0.15  # Cluster separation
}


def get_best_label_col(adata):
    """Find the best available label column."""
    for col in LABEL_COLS:
        if col in adata.obs.columns:
            # Check that it has valid annotations
            valid = ~adata.obs[col].isin(['Unknown', 'Transfer_Failed', 'Not_Run', 'Failed'])
            if valid.sum() > adata.n_obs * 0.3:  # At least 30% valid
                return col
    return None


def compute_purity(cluster_labels, cell_type_labels):
    """Compute mean cluster purity."""
    purities = []
    for cluster in np.unique(cluster_labels):
        mask = cluster_labels == cluster
        if mask.sum() == 0:
            continue
        types = cell_type_labels[mask]
        # Filter out Unknown types
        valid_types = types[~np.isin(types, ['Unknown', 'Transfer_Failed'])]
        if len(valid_types) > 0:
            type_counts = Counter(valid_types)
            dominant_count = type_counts.most_common(1)[0][1]
            purity = dominant_count / len(valid_types)
            purities.append(purity)

    return np.mean(purities) if purities else 0.0


def compute_metrics_at_resolution(adata, resolution, label_col, hvg_setting='all'):
    """Compute all metrics at a given resolution."""
    # Make a copy
    ad = adata.copy()

    # Select HVGs if specified
    if hvg_setting != 'all':
        n_hvg = int(hvg_setting)
        try:
            sc.pp.highly_variable_genes(
                ad, n_top_genes=min(n_hvg, ad.n_vars - 1),
                flavor='seurat_v3', layer='counts', span=0.3
            )
            ad = ad[:, ad.var.highly_variable].copy()
        except:
            pass

    # PCA
    n_pcs = min(30, ad.n_vars - 1)
    sc.pp.pca(ad, n_comps=n_pcs)

    # Neighbors
    sc.pp.neighbors(ad, n_neighbors=15, n_pcs=n_pcs)

    # Cluster at this resolution
    sc.tl.leiden(ad, resolution=resolution, key_added='leiden_test')

    # Get labels
    cluster_labels = ad.obs['leiden_test'].astype(str).values

    # Get cell type labels
    cell_type_labels = ad.obs[label_col].astype(str).values

    # Filter to valid (non-Unknown) labels
    valid_mask = ~np.isin(cell_type_labels, ['Unknown', 'Transfer_Failed', 'Not_Run', 'Failed'])

    metrics = {
        'resolution': resolution,
        'n_clusters': len(np.unique(cluster_labels)),
    }

    # NMI - normalized mutual information
    if valid_mask.sum() > 100:
        metrics['nmi'] = normalized_mutual_info_score(
            cell_type_labels[valid_mask],
            cluster_labels[valid_mask]
        )

        # ARI - adjusted rand index
        metrics['ari'] = adjusted_rand_score(
            cell_type_labels[valid_mask],
            cluster_labels[valid_mask]
        )

        # Purity
        metrics['purity'] = compute_purity(cluster_labels, cell_type_labels)
    else:
        metrics['nmi'] = 0.0
        metrics['ari'] = 0.0
        metrics['purity'] = 0.0

    # Silhouette score (in PCA space)
    try:
        X_pca = ad.obsm['X_pca']
        # Sample if too many cells
        if ad.n_obs > 10000:
            idx = np.random.choice(ad.n_obs, 10000, replace=False)
            X_sub = X_pca[idx]
            labels_sub = cluster_labels[idx]
        else:
            X_sub = X_pca
            labels_sub = cluster_labels

        if len(np.unique(labels_sub)) > 1:
            metrics['silhouette'] = silhouette_score(X_sub, labels_sub)
        else:
            metrics['silhouette'] = 0.0
    except:
        metrics['silhouette'] = 0.0

    # Composite score
    metrics['composite'] = (
        WEIGHTS['nmi'] * metrics['nmi'] +
        WEIGHTS['ari'] * metrics['ari'] +
        WEIGHTS['purity'] * metrics['purity'] +
        WEIGHTS['silhouette'] * max(0, metrics['silhouette'])  # Silhouette can be negative
    )

    return metrics


def optimize_sample(sample_id: str, stage: str):
    """Run resolution optimization for one sample."""
    # Try to load reference-annotated first, fall back to regular annotated
    ref_path = RESULTS / f'{sample_id}_reference_annotated.h5ad'
    ann_path = RESULTS / f'{sample_id}_annotated.h5ad'

    if ref_path.exists():
        adata = sc.read_h5ad(ref_path)
        print(f"    Loaded reference-annotated: {ref_path}")
    elif ann_path.exists():
        adata = sc.read_h5ad(ann_path)
        print(f"    Loaded: {ann_path}")
    else:
        print(f"    WARNING: No annotated data found for {sample_id}")
        return None

    # Find best label column
    label_col = get_best_label_col(adata)
    if label_col is None:
        print(f"    WARNING: No valid label column found for {sample_id}")
        return None

    print(f"    Using labels: {label_col}")
    print(f"    Cells: {adata.n_obs:,}")

    # Test each resolution
    all_metrics = []
    for res in RESOLUTIONS:
        metrics = compute_metrics_at_resolution(adata, res, label_col)
        metrics['sample_id'] = sample_id
        metrics['stage'] = stage
        metrics['label_col'] = label_col
        all_metrics.append(metrics)

        print(f"      res={res:.1f}: NMI={metrics['nmi']:.3f}, "
              f"ARI={metrics['ari']:.3f}, Purity={metrics['purity']:.3f}, "
              f"Clusters={metrics['n_clusters']}")

    return pd.DataFrame(all_metrics)


def find_optimal_resolution(results_df: pd.DataFrame) -> dict:
    """Find optimal resolution based on composite score."""
    # Mean composite score across samples
    avg_by_res = results_df.groupby('resolution')['composite'].mean()
    optimal_res = avg_by_res.idxmax()
    optimal_score = avg_by_res.max()

    # Per-sample optimal
    per_sample = results_df.loc[
        results_df.groupby('sample_id')['composite'].idxmax(),
        ['sample_id', 'resolution', 'composite']
    ]

    return {
        'optimal_resolution': optimal_res,
        'optimal_score': optimal_score,
        'per_sample_optimal': per_sample.to_dict('records'),
        'avg_by_resolution': avg_by_res.to_dict()
    }


def recluster_at_optimal(sample_id: str, optimal_res: float):
    """Recluster sample at optimal resolution and save."""
    # Load
    ref_path = RESULTS / f'{sample_id}_reference_annotated.h5ad'
    ann_path = RESULTS / f'{sample_id}_annotated.h5ad'

    if ref_path.exists():
        adata = sc.read_h5ad(ref_path)
    elif ann_path.exists():
        adata = sc.read_h5ad(ann_path)
    else:
        return

    # Recluster
    if 'X_pca' not in adata.obsm:
        n_pcs = min(30, adata.n_vars - 1)
        sc.pp.pca(adata, n_comps=n_pcs)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)

    sc.tl.leiden(adata, resolution=optimal_res, key_added='leiden_optimal')

    # Recompute UMAP for visualization
    sc.tl.umap(adata)

    # Save
    output_path = RESULTS / f'{sample_id}_optimal.h5ad'
    adata.write_h5ad(output_path)
    print(f"  Saved: {output_path}")

    return adata


def create_visualizations(results_df: pd.DataFrame, optimal_info: dict):
    """Create optimization result visualizations."""
    optimal_res = optimal_info['optimal_resolution']

    # 1. Resolution vs Metrics line plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    metrics_to_plot = ['nmi', 'ari', 'purity', 'composite']
    titles = ['NMI (Cluster-Label Agreement)', 'ARI (Adjusted Rand Index)',
              'Cluster Purity', 'Composite Score']
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

    for ax, metric, title, color in zip(axes.flatten(), metrics_to_plot, titles, colors):
        for sample_id in SAMPLES.keys():
            sample_df = results_df[results_df['sample_id'] == sample_id]
            ax.plot(sample_df['resolution'], sample_df[metric],
                   marker='o', label=sample_id, alpha=0.7)

        # Mean line
        mean_by_res = results_df.groupby('resolution')[metric].mean()
        ax.plot(mean_by_res.index, mean_by_res.values,
               linewidth=2.5, color='black', linestyle='--', label='Mean')

        # Mark optimal
        ax.axvline(optimal_res, color='red', linestyle=':', alpha=0.7, label=f'Optimal ({optimal_res})')

        ax.set_xlabel('Leiden Resolution')
        ax.set_ylabel(metric.upper())
        ax.set_title(title, fontweight='bold')
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.suptitle('Resolution Optimization Results', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(FIGURES / 'resolution_optimization.png', dpi=150)
    plt.close()
    print(f"  Saved: {FIGURES / 'resolution_optimization.png'}")

    # 2. Clusters vs Resolution
    fig, ax = plt.subplots(figsize=(10, 6))

    for sample_id, stage in SAMPLES.items():
        sample_df = results_df[results_df['sample_id'] == sample_id]
        ax.plot(sample_df['resolution'], sample_df['n_clusters'],
               marker='o', label=f'{sample_id} ({stage})', linewidth=2)

    ax.axvline(optimal_res, color='red', linestyle=':', alpha=0.7,
              label=f'Optimal ({optimal_res})')
    ax.set_xlabel('Leiden Resolution')
    ax.set_ylabel('Number of Clusters')
    ax.set_title('Cluster Count vs Resolution', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(FIGURES / 'clusters_vs_resolution.png', dpi=150)
    plt.close()
    print(f"  Saved: {FIGURES / 'clusters_vs_resolution.png'}")

    # 3. UMAP at optimal resolution
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for i, (sample_id, stage) in enumerate(SAMPLES.items()):
        ax = axes[i]

        # Load optimal
        opt_path = RESULTS / f'{sample_id}_optimal.h5ad'
        if opt_path.exists():
            adata = sc.read_h5ad(opt_path)

            if 'X_umap' in adata.obsm and 'leiden_optimal' in adata.obs.columns:
                sc.pl.umap(adata, color='leiden_optimal', ax=ax, show=False,
                          legend_loc='on data', legend_fontsize=6,
                          title=f'{sample_id}: {stage}\n(res={optimal_res})')
            else:
                ax.text(0.5, 0.5, 'UMAP not available', ha='center', va='center')
                ax.set_title(f'{sample_id}: {stage}')

    plt.suptitle(f'Optimal Clustering (Resolution = {optimal_res})', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(FIGURES / 'optimal_clusters_umap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {FIGURES / 'optimal_clusters_umap.png'}")


def create_validation_plot(results_df: pd.DataFrame):
    """Create validation plot comparing marker vs reference proportions."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, (sample_id, stage) in enumerate(SAMPLES.items()):
        ax = axes[i]

        # Load data with both annotations
        ref_path = RESULTS / f'{sample_id}_reference_annotated.h5ad'
        if not ref_path.exists():
            ref_path = RESULTS / f'{sample_id}_annotated.h5ad'

        if ref_path.exists():
            adata = sc.read_h5ad(ref_path)

            # Get proportions
            if 'celltype_markers' in adata.obs.columns:
                marker_props = adata.obs['celltype_markers'].value_counts(normalize=True)
            else:
                marker_props = pd.Series()

            if 'celltype_reference' in adata.obs.columns:
                ref_props = adata.obs['celltype_reference'].value_counts(normalize=True)
            else:
                ref_props = pd.Series()

            # Plot comparison for key types
            key_types = ['Goblet', 'Epithelial', 'Fibroblast', 'T_Cell', 'Macrophage']
            x = np.arange(len(key_types))
            width = 0.35

            marker_vals = [marker_props.get(ct, 0) * 100 for ct in key_types]
            ref_vals = [ref_props.get(ct, 0) * 100 for ct in key_types]

            ax.bar(x - width/2, marker_vals, width, label='Marker-based', color='steelblue')
            ax.bar(x + width/2, ref_vals, width, label='Reference-based', color='coral')

            ax.set_xlabel('Cell Type')
            ax.set_ylabel('Proportion (%)')
            ax.set_title(f'{sample_id}: {stage}', fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels(key_types, rotation=45, ha='right')
            ax.legend()

    plt.suptitle('Marker vs Reference Annotation Comparison', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(FIGURES / 'validation_comparison.png', dpi=150)
    plt.close()
    print(f"  Saved: {FIGURES / 'validation_comparison.png'}")


def main():
    print("=" * 60)
    print("Multi-Resolution Clustering Optimization")
    print("=" * 60)

    # Step 1: Run optimization for each sample
    print("\n[1/4] Testing resolutions across samples...")
    all_results = []

    for sample_id, stage in SAMPLES.items():
        print(f"\n  {sample_id} ({stage}):")
        sample_results = optimize_sample(sample_id, stage)
        if sample_results is not None:
            all_results.append(sample_results)

    if not all_results:
        print("ERROR: No samples could be optimized")
        return

    results_df = pd.concat(all_results, ignore_index=True)

    # Save raw results
    results_df.to_csv(RESULTS / 'resolution_optimization.csv', index=False)
    print(f"\nSaved: {RESULTS / 'resolution_optimization.csv'}")

    # Step 2: Find optimal resolution
    print("\n[2/4] Finding optimal resolution...")
    optimal_info = find_optimal_resolution(results_df)

    print(f"\n  OPTIMAL RESOLUTION: {optimal_info['optimal_resolution']}")
    print(f"  Mean composite score: {optimal_info['optimal_score']:.4f}")

    print("\n  Per-sample optimal:")
    for item in optimal_info['per_sample_optimal']:
        print(f"    {item['sample_id']}: res={item['resolution']}, score={item['composite']:.3f}")

    # Save recommendation
    with open(RESULTS / 'optimal_resolution.txt', 'w') as f:
        f.write(f"Optimal Resolution: {optimal_info['optimal_resolution']}\n")
        f.write(f"Mean Composite Score: {optimal_info['optimal_score']:.4f}\n")
        f.write(f"\nWeight breakdown:\n")
        for metric, weight in WEIGHTS.items():
            f.write(f"  {metric}: {weight}\n")
        f.write(f"\nPer-sample optimal:\n")
        for item in optimal_info['per_sample_optimal']:
            f.write(f"  {item['sample_id']}: {item['resolution']} ({item['composite']:.3f})\n")

    print(f"Saved: {RESULTS / 'optimal_resolution.txt'}")

    # Step 3: Recluster at optimal resolution
    print("\n[3/4] Reclustering at optimal resolution...")
    optimal_res = optimal_info['optimal_resolution']

    for sample_id in SAMPLES.keys():
        recluster_at_optimal(sample_id, optimal_res)

    # Step 4: Create visualizations
    print("\n[4/4] Creating visualizations...")
    create_visualizations(results_df, optimal_info)
    create_validation_plot(results_df)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"\nOptimal Resolution: {optimal_res}")
    print(f"Composite Score: {optimal_info['optimal_score']:.4f}")

    print(f"\nMetrics at optimal resolution:")
    opt_metrics = results_df[results_df['resolution'] == optimal_res]
    for _, row in opt_metrics.iterrows():
        print(f"  {row['sample_id']}: NMI={row['nmi']:.3f}, ARI={row['ari']:.3f}, "
              f"Purity={row['purity']:.3f}, Clusters={row['n_clusters']}")

    print(f"\nOutputs saved to: {RESULTS}")
    print("\nKey files:")
    print(f"  - resolution_optimization.csv: Metrics at all resolutions")
    print(f"  - optimal_resolution.txt: Recommended resolution")
    print(f"  - *_optimal.h5ad: Data reclustered at optimal resolution")


if __name__ == '__main__':
    main()
