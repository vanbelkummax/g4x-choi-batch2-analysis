#!/usr/bin/env python3
"""
Spatial Biology Hackathon 2026 - Spatial Analysis
==================================================

Squidpy-based spatial statistics:
- Spatial neighborhood graphs
- Neighborhood enrichment (co-localization)
- Spatial autocorrelation (Moran's I)
- Spatially variable genes (SVGs)

Author: Max Van Belkum
Date: 2026-01-20
"""

import scanpy as sc
import squidpy as sq
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
FIG_DIR = OUTPUT_DIR / "figures" / "spatial"
TABLE_DIR = OUTPUT_DIR / "tables"

# Create directories
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)
(ADATA_DIR / "spatial").mkdir(parents=True, exist_ok=True)

# Spatial parameters
N_NEIGHBORS_SPATIAL = 6  # Hexagonal grid for Visium
COORD_TYPE = {
    'Visium': 'grid',
    'G4X': 'generic'
}

# =============================================================================
# Spatial Graph Construction
# =============================================================================

def build_spatial_graph(
    adata: ad.AnnData,
    n_neighbors: int = N_NEIGHBORS_SPATIAL,
    platform: Optional[str] = None
) -> ad.AnnData:
    """
    Build spatial connectivity graph.

    Args:
        adata: AnnData with spatial coordinates in .obsm['spatial']
        n_neighbors: Number of spatial neighbors
        platform: Platform type for coordinate handling

    Returns:
        AnnData with spatial_connectivities and spatial_distances
    """
    if platform is None:
        platform = adata.obs['platform'].iloc[0] if 'platform' in adata.obs else 'Visium'

    coord_type = COORD_TYPE.get(platform, 'generic')

    print(f"  Building spatial graph: {coord_type} coordinates, n_neighbors={n_neighbors}")

    try:
        if coord_type == 'grid':
            # Visium hexagonal grid
            sq.gr.spatial_neighbors(
                adata,
                n_neighs=n_neighbors,
                coord_type='grid',
                spatial_key='spatial'
            )
        else:
            # Generic coordinates (G4X, etc.)
            sq.gr.spatial_neighbors(
                adata,
                n_neighs=n_neighbors,
                coord_type='generic',
                spatial_key='spatial'
            )

        n_edges = adata.obsp['spatial_connectivities'].sum() / 2
        print(f"  Created spatial graph with {int(n_edges)} edges")

    except Exception as e:
        print(f"  Warning: Failed to build spatial graph: {e}")
        print(f"  Attempting fallback with generic coordinates...")

        # Fallback: use KNN on spatial coordinates
        from sklearn.neighbors import kneighbors_graph

        if 'spatial' in adata.obsm:
            spatial_coords = adata.obsm['spatial']
            conn = kneighbors_graph(spatial_coords, n_neighbors=n_neighbors, mode='connectivity')
            adata.obsp['spatial_connectivities'] = conn
            adata.obsp['spatial_distances'] = kneighbors_graph(spatial_coords, n_neighbors=n_neighbors, mode='distance')

    return adata


# =============================================================================
# Neighborhood Enrichment
# =============================================================================

def compute_neighborhood_enrichment(
    adata: ad.AnnData,
    cluster_key: str = 'leiden',
    n_perms: int = 500
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute neighborhood enrichment (co-localization patterns).

    This identifies which cell types/clusters tend to be spatially close.

    Args:
        adata: AnnData with spatial graph and clusters
        cluster_key: Key in .obs with cluster labels
        n_perms: Number of permutations for significance testing

    Returns:
        z_scores: Enrichment z-scores matrix
        counts: Observed neighbor counts matrix
    """
    print(f"  Computing neighborhood enrichment for {cluster_key}...")

    # Ensure spatial graph exists
    if 'spatial_connectivities' not in adata.obsp:
        adata = build_spatial_graph(adata)

    try:
        sq.gr.nhood_enrichment(
            adata,
            cluster_key=cluster_key,
            n_perms=n_perms,
            seed=42
        )

        z_scores = adata.uns[f'{cluster_key}_nhood_enrichment']['zscore']
        counts = adata.uns[f'{cluster_key}_nhood_enrichment'].get('count', None)

        n_clusters = z_scores.shape[0]
        print(f"  Computed enrichment matrix: {n_clusters}x{n_clusters}")

        return z_scores, counts

    except Exception as e:
        print(f"  Warning: Neighborhood enrichment failed: {e}")
        return None, None


def plot_neighborhood_enrichment(
    adata: ad.AnnData,
    cluster_key: str,
    sample_name: str
) -> None:
    """Plot neighborhood enrichment heatmap."""

    fig, ax = plt.subplots(figsize=(8, 6))

    try:
        sq.pl.nhood_enrichment(
            adata,
            cluster_key=cluster_key,
            ax=ax,
            title=f'{sample_name} - Neighborhood Enrichment'
        )
    except:
        # Manual plotting if squidpy plotting fails
        z_scores = adata.uns.get(f'{cluster_key}_nhood_enrichment', {}).get('zscore')
        if z_scores is not None:
            clusters = adata.obs[cluster_key].cat.categories
            sns.heatmap(z_scores, xticklabels=clusters, yticklabels=clusters,
                       center=0, cmap='RdBu_r', ax=ax)
            ax.set_title(f'{sample_name} - Neighborhood Enrichment (z-scores)')

    plt.tight_layout()
    plt.savefig(FIG_DIR / f"{sample_name}_neighborhood_enrichment.png", dpi=150, bbox_inches='tight')
    plt.close()


# =============================================================================
# Spatial Autocorrelation
# =============================================================================

def compute_spatial_autocorrelation(
    adata: ad.AnnData,
    genes: Optional[List[str]] = None,
    n_genes: int = 100,
    mode: str = 'moran'
) -> pd.DataFrame:
    """
    Compute spatial autocorrelation (Moran's I or Geary's C).

    Identifies genes with non-random spatial patterns.

    Args:
        adata: AnnData with spatial graph
        genes: Specific genes to test (None = use HVGs or top expressed)
        n_genes: Number of top genes if genes not specified
        mode: 'moran' for Moran's I or 'geary' for Geary's C

    Returns:
        DataFrame with autocorrelation statistics per gene
    """
    print(f"  Computing spatial autocorrelation ({mode})...")

    # Ensure spatial graph exists
    if 'spatial_connectivities' not in adata.obsp:
        adata = build_spatial_graph(adata)

    # Select genes
    if genes is None:
        if 'highly_variable' in adata.var and adata.var['highly_variable'].sum() > 0:
            genes = adata.var_names[adata.var['highly_variable']][:n_genes].tolist()
        else:
            # Use top expressed genes
            mean_expr = np.asarray(adata.X.mean(axis=0)).flatten()
            top_idx = np.argsort(mean_expr)[-n_genes:]
            genes = adata.var_names[top_idx].tolist()

    print(f"  Testing {len(genes)} genes...")

    try:
        sq.gr.spatial_autocorr(
            adata,
            mode=mode,
            genes=genes,
            n_perms=100,
            n_jobs=4
        )

        result_key = f'moranI' if mode == 'moran' else f'gearyC'
        results = adata.uns[result_key].copy()
        results = results.sort_values('I' if mode == 'moran' else 'C', ascending=False)

        print(f"  Top 5 spatially autocorrelated genes:")
        print(results.head())

        return results

    except Exception as e:
        print(f"  Warning: Spatial autocorrelation failed: {e}")
        return pd.DataFrame()


# =============================================================================
# Spatially Variable Genes
# =============================================================================

def find_spatially_variable_genes(
    adata: ad.AnnData,
    n_genes: int = 100,
    method: str = 'moran'
) -> pd.DataFrame:
    """
    Identify spatially variable genes (SVGs).

    Args:
        adata: AnnData with spatial graph
        n_genes: Number of genes to test
        method: Method for SVG detection ('moran' or 'geary')

    Returns:
        DataFrame with SVG statistics
    """
    results = compute_spatial_autocorrelation(adata, n_genes=n_genes, mode=method)

    if len(results) > 0:
        # Filter significant SVGs
        if 'pval_norm' in results.columns:
            svg_df = results[results['pval_norm'] < 0.05].copy()
        else:
            svg_df = results.head(50).copy()  # Top 50 by score

        print(f"  Found {len(svg_df)} significant SVGs")
        return svg_df

    return pd.DataFrame()


def plot_spatially_variable_genes(
    adata: ad.AnnData,
    svg_df: pd.DataFrame,
    sample_name: str,
    n_plot: int = 6
) -> None:
    """Plot top spatially variable genes."""

    if len(svg_df) == 0:
        print(f"  No SVGs to plot for {sample_name}")
        return

    top_genes = svg_df.head(n_plot).index.tolist()
    available_genes = [g for g in top_genes if g in adata.var_names]

    if len(available_genes) == 0:
        return

    n_cols = min(3, len(available_genes))
    n_rows = (len(available_genes) + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_rows * n_cols == 1:
        axes = np.array([[axes]])
    axes = axes.flatten()

    for i, gene in enumerate(available_genes):
        ax = axes[i]
        if 'spatial' in adata.obsm:
            sc.pl.embedding(adata, basis='spatial', color=gene, ax=ax, show=False,
                           title=gene, use_raw=False, layer='normalized' if 'normalized' in adata.layers else None)
        else:
            sc.pl.umap(adata, color=gene, ax=ax, show=False, title=gene)

    # Hide unused axes
    for i in range(len(available_genes), len(axes)):
        axes[i].axis('off')

    plt.suptitle(f'{sample_name} - Spatially Variable Genes', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / f"{sample_name}_svg_genes.png", dpi=150, bbox_inches='tight')
    plt.close()


# =============================================================================
# Main Pipeline
# =============================================================================

def run_spatial_analysis(
    adata: ad.AnnData,
    sample_name: str,
    cluster_key: str = 'leiden',
    save: bool = True
) -> ad.AnnData:
    """
    Full spatial analysis pipeline for one sample.

    Args:
        adata: Clustered AnnData
        sample_name: Sample identifier
        cluster_key: Cluster annotation to use
        save: Whether to save outputs

    Returns:
        AnnData with spatial analysis results
    """
    print(f"\n{'='*60}")
    print(f"Spatial Analysis: {sample_name}")
    print(f"{'='*60}")

    # Check for spatial coordinates
    if 'spatial' not in adata.obsm:
        print(f"  WARNING: No spatial coordinates found, skipping spatial analysis")
        return adata

    # 1. Build spatial graph
    print("\n1. Building spatial graph...")
    adata = build_spatial_graph(adata)

    # 2. Neighborhood enrichment
    print("\n2. Computing neighborhood enrichment...")
    if cluster_key in adata.obs:
        z_scores, counts = compute_neighborhood_enrichment(adata, cluster_key)
        if z_scores is not None and save:
            plot_neighborhood_enrichment(adata, cluster_key, sample_name)
    else:
        print(f"  Skipping: {cluster_key} not found")

    # 3. Spatial autocorrelation
    print("\n3. Computing spatial autocorrelation...")
    autocorr_df = compute_spatial_autocorrelation(adata, n_genes=100)
    if len(autocorr_df) > 0:
        adata.uns['spatial_autocorr'] = autocorr_df

    # 4. Spatially variable genes
    print("\n4. Identifying spatially variable genes...")
    svg_df = find_spatially_variable_genes(adata, n_genes=100)
    if len(svg_df) > 0:
        adata.uns['svg'] = svg_df
        if save:
            plot_spatially_variable_genes(adata, svg_df, sample_name)

    # Save
    if save:
        out_path = ADATA_DIR / "spatial" / f"{sample_name}_spatial.h5ad"
        adata.write_h5ad(out_path)
        print(f"\n  Saved to: {out_path}")

    return adata


def run_spatial_analysis_all(
    input_dir: Optional[Path] = None,
    samples: Optional[List[str]] = None
) -> Dict[str, ad.AnnData]:
    """
    Run spatial analysis on all clustered samples.

    Args:
        input_dir: Directory with clustered h5ad files
        samples: Specific samples to process (None = all)

    Returns:
        Dictionary of AnnData objects with spatial analysis
    """
    if input_dir is None:
        input_dir = ADATA_DIR / "clustered"

    # Find all clustered h5ad files
    h5ad_files = list(input_dir.glob("*_clustered.h5ad"))

    if samples:
        h5ad_files = [f for f in h5ad_files if f.stem.replace('_clustered', '') in samples]

    print(f"Found {len(h5ad_files)} samples for spatial analysis")

    analyzed = {}
    all_svgs = []

    for h5ad_path in sorted(h5ad_files):
        sample_name = h5ad_path.stem.replace('_clustered', '')

        try:
            # Load
            adata = sc.read_h5ad(h5ad_path)

            # Run spatial analysis
            adata = run_spatial_analysis(adata, sample_name)
            analyzed[sample_name] = adata

            # Collect SVGs
            if 'svg' in adata.uns:
                svg_df = adata.uns['svg'].copy()
                svg_df['sample'] = sample_name
                all_svgs.append(svg_df)

        except Exception as e:
            print(f"  ERROR in spatial analysis for {sample_name}: {e}")
            import traceback
            traceback.print_exc()

    # Save combined SVG table
    if all_svgs:
        combined_svgs = pd.concat(all_svgs, ignore_index=False)
        combined_svgs.to_csv(TABLE_DIR / "spatially_variable_genes.csv")
        print(f"\nSVG table saved to: {TABLE_DIR / 'spatially_variable_genes.csv'}")

    return analyzed


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Spatial analysis of clustered samples")
    parser.add_argument("--sample", type=str, help="Single sample to process")
    parser.add_argument("--input-dir", type=str, help="Input directory with clustered h5ad files")
    args = parser.parse_args()

    print("="*60)
    print("SPATIAL BIOLOGY HACKATHON 2026 - SPATIAL ANALYSIS")
    print("="*60)

    if args.sample:
        # Single sample mode
        input_dir = Path(args.input_dir) if args.input_dir else ADATA_DIR / "clustered"
        adata = sc.read_h5ad(input_dir / f"{args.sample}_clustered.h5ad")
        run_spatial_analysis(adata, args.sample)
    else:
        # Batch mode
        input_dir = Path(args.input_dir) if args.input_dir else None
        run_spatial_analysis_all(input_dir)

    print("\n" + "="*60)
    print("SPATIAL ANALYSIS COMPLETE")
    print("="*60)
