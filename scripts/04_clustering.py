#!/usr/bin/env python3
"""
Spatial Biology Hackathon 2026 - Clustering Pipeline
=====================================================

Leiden clustering with multi-resolution evaluation.

Platform-specific processing:
- Visium: Standard HVG selection, full preprocessing
- G4X: Skip HVG (targeted panel), scale with max_value

Author: Max Van Belkum
Date: 2026-01-20
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Optional, Dict, List, Tuple
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Configuration
# =============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
ADATA_DIR = OUTPUT_DIR / "adata"
FIG_DIR = OUTPUT_DIR / "figures" / "clustering"
TABLE_DIR = OUTPUT_DIR / "tables"

# Create directories
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)
(ADATA_DIR / "clustered").mkdir(parents=True, exist_ok=True)

# Clustering parameters
RESOLUTIONS = [0.3, 0.5, 0.8, 1.0, 1.5]
DEFAULT_RESOLUTION = 0.5  # Start conservative for spatial data

# =============================================================================
# Preprocessing Functions
# =============================================================================

def preprocess_visium(adata: ad.AnnData, n_top_genes: int = 2000) -> ad.AnnData:
    """
    Standard Visium preprocessing: normalize, HVG, PCA, neighbors, UMAP.

    Args:
        adata: Raw Visium AnnData
        n_top_genes: Number of highly variable genes to select

    Returns:
        Preprocessed AnnData
    """
    print(f"  Preprocessing Visium: {adata.n_obs} spots, {adata.n_vars} genes")

    # Store raw counts
    adata.layers['counts'] = adata.X.copy()

    # Normalize to 10,000 counts per spot
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Log transform
    sc.pp.log1p(adata)

    # Store normalized for later
    adata.layers['normalized'] = adata.X.copy()

    # Highly variable genes
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor='seurat_v3',
        layer='counts'  # Use counts for HVG selection
    )
    print(f"  Selected {adata.var['highly_variable'].sum()} HVGs")

    # Scale for PCA
    sc.pp.scale(adata, max_value=10)

    # PCA
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)

    # Neighbors and UMAP
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)

    return adata


def preprocess_g4x(adata: ad.AnnData) -> ad.AnnData:
    """
    G4X preprocessing: targeted panel, skip HVG selection.

    Args:
        adata: Raw G4X AnnData

    Returns:
        Preprocessed AnnData
    """
    print(f"  Preprocessing G4X: {adata.n_obs} cells, {adata.n_vars} features")

    # Store raw counts
    adata.layers['counts'] = adata.X.copy()

    # Normalize per cell
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Log transform
    sc.pp.log1p(adata)

    # Store normalized
    adata.layers['normalized'] = adata.X.copy()

    # Mark all genes as highly variable (targeted panel)
    adata.var['highly_variable'] = True

    # Scale
    sc.pp.scale(adata, max_value=10)

    # PCA - use all features since it's a targeted panel
    n_comps = min(50, adata.n_vars - 1, adata.n_obs - 1)
    sc.tl.pca(adata, n_comps=n_comps)

    # Neighbors and UMAP
    n_pcs = min(30, n_comps)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
    sc.tl.umap(adata)

    return adata


def preprocess_sample(adata: ad.AnnData) -> ad.AnnData:
    """Route to platform-specific preprocessing."""
    platform = adata.obs['platform'].iloc[0] if 'platform' in adata.obs else 'Visium'

    if platform == 'Visium':
        return preprocess_visium(adata)
    elif platform == 'G4X':
        return preprocess_g4x(adata)
    else:
        print(f"  Unknown platform {platform}, using Visium preprocessing")
        return preprocess_visium(adata)


# =============================================================================
# Clustering Functions
# =============================================================================

def run_leiden_multiresolution(
    adata: ad.AnnData,
    resolutions: List[float] = RESOLUTIONS
) -> ad.AnnData:
    """
    Run Leiden clustering at multiple resolutions.

    Args:
        adata: Preprocessed AnnData with neighbors computed
        resolutions: List of resolution parameters to try

    Returns:
        AnnData with leiden columns for each resolution
    """
    print(f"  Running Leiden at resolutions: {resolutions}")

    for res in resolutions:
        key = f"leiden_{str(res).replace('.', '_')}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        n_clusters = adata.obs[key].nunique()
        print(f"    Resolution {res}: {n_clusters} clusters")

    # Set default resolution
    adata.obs['leiden'] = adata.obs[f"leiden_{str(DEFAULT_RESOLUTION).replace('.', '_')}"].copy()

    return adata


def evaluate_clustering_stability(
    adata: ad.AnnData,
    resolutions: List[float] = RESOLUTIONS
) -> pd.DataFrame:
    """
    Evaluate clustering quality across resolutions.

    Returns DataFrame with cluster counts and sizes per resolution.
    """
    results = []

    for res in resolutions:
        key = f"leiden_{str(res).replace('.', '_')}"
        if key in adata.obs:
            clusters = adata.obs[key]
            sizes = clusters.value_counts()
            results.append({
                'resolution': res,
                'n_clusters': len(sizes),
                'min_size': sizes.min(),
                'max_size': sizes.max(),
                'median_size': sizes.median(),
                'mean_size': sizes.mean()
            })

    return pd.DataFrame(results)


def select_optimal_resolution(
    adata: ad.AnnData,
    resolutions: List[float] = RESOLUTIONS,
    min_cluster_size: int = 10,
    target_n_clusters: Optional[int] = None
) -> float:
    """
    Select optimal resolution based on cluster quality.

    Criteria:
    - Minimum cluster size >= min_cluster_size
    - If target_n_clusters provided, closest to that
    - Otherwise, prefer moderate number of clusters
    """
    stats = evaluate_clustering_stability(adata, resolutions)

    # Filter by minimum cluster size
    valid = stats[stats['min_size'] >= min_cluster_size]

    if len(valid) == 0:
        print(f"  Warning: No resolution has min cluster size >= {min_cluster_size}")
        return DEFAULT_RESOLUTION

    if target_n_clusters:
        # Find closest to target
        valid['dist_to_target'] = abs(valid['n_clusters'] - target_n_clusters)
        best = valid.loc[valid['dist_to_target'].idxmin()]
    else:
        # Prefer moderate number (5-15 clusters for spatial)
        valid['quality'] = abs(valid['n_clusters'] - 8)
        best = valid.loc[valid['quality'].idxmin()]

    return best['resolution']


# =============================================================================
# Visualization Functions
# =============================================================================

def plot_clustering(
    adata: ad.AnnData,
    sample_name: str,
    resolution: float = DEFAULT_RESOLUTION
) -> None:
    """Generate clustering visualization plots."""

    key = f"leiden_{str(resolution).replace('.', '_')}"
    if key not in adata.obs:
        key = 'leiden'

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # UMAP colored by clusters
    sc.pl.umap(adata, color=key, ax=axes[0], show=False, legend_loc='on data',
               title=f'{sample_name} - UMAP Clusters (res={resolution})')

    # Spatial plot if coordinates exist
    if 'spatial' in adata.obsm:
        sc.pl.embedding(adata, basis='spatial', color=key, ax=axes[1], show=False,
                       title=f'{sample_name} - Spatial Clusters')
    else:
        axes[1].text(0.5, 0.5, 'No spatial coordinates', ha='center', va='center')
        axes[1].set_title(f'{sample_name} - Spatial (N/A)')

    # Cluster sizes
    cluster_sizes = adata.obs[key].value_counts().sort_index()
    axes[2].bar(cluster_sizes.index.astype(str), cluster_sizes.values)
    axes[2].set_xlabel('Cluster')
    axes[2].set_ylabel('Count')
    axes[2].set_title(f'{sample_name} - Cluster Sizes')
    axes[2].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(FIG_DIR / f"{sample_name}_clustering.png", dpi=150, bbox_inches='tight')
    plt.close()


def plot_resolution_comparison(
    adata: ad.AnnData,
    sample_name: str,
    resolutions: List[float] = RESOLUTIONS
) -> None:
    """Plot clustering at different resolutions."""

    n_res = len(resolutions)
    fig, axes = plt.subplots(2, n_res, figsize=(4*n_res, 8))

    for i, res in enumerate(resolutions):
        key = f"leiden_{str(res).replace('.', '_')}"
        if key not in adata.obs:
            continue

        n_clusters = adata.obs[key].nunique()

        # UMAP
        sc.pl.umap(adata, color=key, ax=axes[0, i], show=False, legend_loc=None,
                   title=f'res={res}, n={n_clusters}')

        # Spatial
        if 'spatial' in adata.obsm:
            sc.pl.embedding(adata, basis='spatial', color=key, ax=axes[1, i],
                           show=False, legend_loc=None, title='')
        else:
            axes[1, i].axis('off')

    plt.suptitle(f'{sample_name} - Resolution Comparison', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / f"{sample_name}_resolution_comparison.png", dpi=150, bbox_inches='tight')
    plt.close()


# =============================================================================
# Main Pipeline
# =============================================================================

def cluster_sample(
    adata: ad.AnnData,
    sample_name: str,
    resolutions: List[float] = RESOLUTIONS,
    save: bool = True
) -> ad.AnnData:
    """
    Full clustering pipeline for one sample.

    Args:
        adata: Raw or QC'd AnnData
        sample_name: Sample identifier
        resolutions: Resolutions to test
        save: Whether to save outputs

    Returns:
        Clustered AnnData
    """
    print(f"\n{'='*60}")
    print(f"Clustering: {sample_name}")
    print(f"{'='*60}")

    # Check if already preprocessed
    if 'pca' not in adata.obsm:
        print("  Running preprocessing...")
        adata = preprocess_sample(adata)
    else:
        print("  Already preprocessed, skipping...")

    # Run clustering at multiple resolutions
    adata = run_leiden_multiresolution(adata, resolutions)

    # Evaluate and select optimal
    stats = evaluate_clustering_stability(adata, resolutions)
    print(f"\n  Resolution comparison:")
    print(stats.to_string(index=False))

    optimal_res = select_optimal_resolution(adata, resolutions)
    print(f"\n  Selected resolution: {optimal_res}")

    # Update default leiden to optimal
    opt_key = f"leiden_{str(optimal_res).replace('.', '_')}"
    adata.obs['leiden'] = adata.obs[opt_key].copy()
    adata.uns['leiden_resolution'] = optimal_res

    # Generate plots
    if save:
        print(f"  Generating plots...")
        plot_clustering(adata, sample_name, optimal_res)
        plot_resolution_comparison(adata, sample_name, resolutions)

    # Save
    if save:
        out_path = ADATA_DIR / "clustered" / f"{sample_name}_clustered.h5ad"
        adata.write_h5ad(out_path)
        print(f"  Saved to: {out_path}")

    return adata


def cluster_all_samples(
    input_dir: Optional[Path] = None,
    samples: Optional[List[str]] = None
) -> Dict[str, ad.AnnData]:
    """
    Cluster all preprocessed samples.

    Args:
        input_dir: Directory with preprocessed h5ad files
        samples: Specific samples to process (None = all)

    Returns:
        Dictionary of clustered AnnData objects
    """
    if input_dir is None:
        input_dir = ADATA_DIR / "preprocessed"

    # Find all h5ad files
    h5ad_files = list(input_dir.glob("*.h5ad"))

    if samples:
        h5ad_files = [f for f in h5ad_files if f.stem in samples]

    print(f"Found {len(h5ad_files)} samples to cluster")

    clustered = {}
    stats_list = []

    for h5ad_path in sorted(h5ad_files):
        sample_name = h5ad_path.stem

        try:
            # Load
            adata = sc.read_h5ad(h5ad_path)

            # Cluster
            adata = cluster_sample(adata, sample_name)
            clustered[sample_name] = adata

            # Collect stats
            n_clusters = adata.obs['leiden'].nunique()
            resolution = adata.uns.get('leiden_resolution', DEFAULT_RESOLUTION)
            stats_list.append({
                'sample': sample_name,
                'platform': adata.obs['platform'].iloc[0],
                'n_cells': adata.n_obs,
                'n_clusters': n_clusters,
                'resolution': resolution,
                'cluster_sizes': adata.obs['leiden'].value_counts().to_dict()
            })

        except Exception as e:
            print(f"  ERROR clustering {sample_name}: {e}")
            import traceback
            traceback.print_exc()

    # Save summary stats
    stats_df = pd.DataFrame(stats_list)
    stats_df.to_csv(TABLE_DIR / "cluster_stats.csv", index=False)
    print(f"\nCluster stats saved to: {TABLE_DIR / 'cluster_stats.csv'}")

    return clustered


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Cluster spatial transcriptomics samples")
    parser.add_argument("--sample", type=str, help="Single sample to process")
    parser.add_argument("--input-dir", type=str, help="Input directory with h5ad files")
    parser.add_argument("--resolution", type=float, default=None, help="Single resolution to use")
    args = parser.parse_args()

    print("="*60)
    print("SPATIAL BIOLOGY HACKATHON 2026 - CLUSTERING")
    print("="*60)

    if args.sample:
        # Single sample mode
        input_dir = Path(args.input_dir) if args.input_dir else ADATA_DIR / "preprocessed"
        adata = sc.read_h5ad(input_dir / f"{args.sample}.h5ad")
        cluster_sample(adata, args.sample)
    else:
        # Batch mode
        input_dir = Path(args.input_dir) if args.input_dir else None
        cluster_all_samples(input_dir)

    print("\n" + "="*60)
    print("CLUSTERING COMPLETE")
    print("="*60)
