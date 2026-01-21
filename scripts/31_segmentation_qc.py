#!/usr/bin/env python3
"""
G4X Segmentation QC and Sensitivity Analysis
=============================================
Spatial Hackathon 2026 - Max Van Belkum

Per Mitchel et al. 2026: Run ALL downstream results three ways:
1. Baseline: Original segmentation
2. QC-strict: Remove top 10% most suspicious cells
3. Corrected: After cellAdmix correction (future)

Report stability: If claims flip across (1-3), it's pipeline sensitivity, not biology.

Usage:
    python 31_segmentation_qc.py [--sample SAMPLE_ID] [--all]
"""

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import sparse, stats
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import argparse
import warnings
import logging
import gc

warnings.filterwarnings('ignore')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# =============================================================================
# Configuration
# =============================================================================

INPUT_DIR = Path("/home/user/spatial-hackathon-2026/results/g4x_choi_batch2/wnn_integrated")
OUTPUT_DIR = Path("/home/user/spatial-hackathon-2026/results/g4x_choi_batch2/segmentation_qc")
FIGURE_DIR = Path("/home/user/spatial-hackathon-2026/figures/g4x/segmentation_qc")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# QC thresholds
QC_CONFIG = {
    # Cell area (pixels^2) - needs calibration for G4X
    'min_cell_area': 50,
    'max_cell_area': 2000,

    # Gene/count filters
    'min_genes': 5,
    'max_genes': 300,
    'min_counts': 10,
    'max_counts': 2000,

    # Admixture filter (from WNN script)
    'admix_strict_percentile': 90,  # Remove top 10% by admixture score

    # Boundary cell detection
    'boundary_distance_threshold': 50,  # pixels from ROI edge

    # Density filter
    'density_percentile': 95,  # Flag high-density regions
}


# =============================================================================
# QC Functions
# =============================================================================

def compute_cell_density(adata, radius=50):
    """
    Compute local cell density (cells within radius of each cell).
    High density regions are prone to segmentation errors.
    """
    logger.info("Computing local cell density...")

    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    elif 'cell_x' in adata.obs and 'cell_y' in adata.obs:
        coords = adata.obs[['cell_x', 'cell_y']].values
    else:
        return np.zeros(adata.n_obs)

    # Use ball tree for efficient radius query
    from sklearn.neighbors import BallTree
    tree = BallTree(coords)
    counts = tree.query_radius(coords, r=radius, count_only=True)

    return counts.astype(float)


def compute_boundary_distance(adata):
    """
    Compute distance of each cell to the ROI boundary.
    Cells near boundaries are more prone to segmentation artifacts.
    """
    logger.info("Computing boundary distances...")

    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    elif 'cell_x' in adata.obs and 'cell_y' in adata.obs:
        coords = adata.obs[['cell_x', 'cell_y']].values
    else:
        return np.ones(adata.n_obs) * np.inf

    # Compute distance to convex hull boundary
    from scipy.spatial import ConvexHull, Delaunay

    try:
        hull = ConvexHull(coords)
        # For each point, find minimum distance to hull edges
        hull_points = coords[hull.vertices]

        distances = np.zeros(adata.n_obs)
        for i in range(adata.n_obs):
            point = coords[i]
            # Distance to each hull edge
            edge_distances = []
            for j in range(len(hull_points)):
                p1 = hull_points[j]
                p2 = hull_points[(j + 1) % len(hull_points)]
                d = point_to_line_distance(point, p1, p2)
                edge_distances.append(d)
            distances[i] = min(edge_distances)

        return distances
    except Exception as e:
        logger.warning(f"Could not compute boundary distances: {e}")
        return np.ones(adata.n_obs) * np.inf


def point_to_line_distance(point, line_p1, line_p2):
    """Calculate perpendicular distance from point to line segment."""
    x0, y0 = point
    x1, y1 = line_p1
    x2, y2 = line_p2

    num = abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1)
    den = np.sqrt((y2 - y1) ** 2 + (x2 - x1) ** 2)

    return num / (den + 1e-10)


def apply_hard_qc_filters(adata, config=QC_CONFIG):
    """
    Apply hard QC filters based on cell/gene metrics.

    Returns mask of cells passing QC.
    """
    logger.info("Applying hard QC filters...")

    n_cells = adata.n_obs
    pass_qc = np.ones(n_cells, dtype=bool)

    # Gene count filter
    if 'n_genes_by_counts' in adata.obs:
        gene_counts = adata.obs['n_genes_by_counts']
    else:
        gene_counts = np.array((adata.X > 0).sum(axis=1)).flatten()

    min_genes_mask = gene_counts >= config['min_genes']
    max_genes_mask = gene_counts <= config['max_genes']
    pass_qc &= min_genes_mask & max_genes_mask

    logger.info(f"  Gene count filter: {np.sum(~(min_genes_mask & max_genes_mask))} cells removed")

    # Total count filter
    if 'total_counts' in adata.obs:
        total_counts = adata.obs['total_counts']
    else:
        total_counts = np.array(adata.X.sum(axis=1)).flatten()

    min_counts_mask = total_counts >= config['min_counts']
    max_counts_mask = total_counts <= config['max_counts']
    pass_qc &= min_counts_mask & max_counts_mask

    logger.info(f"  Total count filter: {np.sum(~(min_counts_mask & max_counts_mask))} cells removed")

    # Cell area filter (if available)
    if 'nuclei_area' in adata.obs:
        area = adata.obs['nuclei_area']
        min_area_mask = area >= config['min_cell_area']
        max_area_mask = area <= config['max_cell_area']
        pass_qc &= min_area_mask & max_area_mask
        logger.info(f"  Cell area filter: {np.sum(~(min_area_mask & max_area_mask))} cells removed")

    logger.info(f"  Total passing hard QC: {np.sum(pass_qc)} / {n_cells} ({100*np.mean(pass_qc):.1f}%)")

    return pass_qc


def compute_qc_scores(adata, config=QC_CONFIG):
    """
    Compute composite QC scores for each cell.

    Combines:
    - Admixture score (from WNN script)
    - Cell density
    - Boundary distance
    - Gene/count metrics
    """
    logger.info("Computing composite QC scores...")

    qc_df = pd.DataFrame(index=adata.obs_names)

    # 1. Admixture score (from WNN script)
    if 'admixture_score' in adata.obs:
        qc_df['admix_score'] = adata.obs['admixture_score']
    else:
        qc_df['admix_score'] = 0

    # 2. Cell density
    density = compute_cell_density(adata)
    density_norm = (density - np.percentile(density, 5)) / (np.percentile(density, 95) - np.percentile(density, 5) + 1e-10)
    qc_df['density_score'] = np.clip(density_norm, 0, 1)
    qc_df['cell_density'] = density

    # 3. Boundary distance (inverted - closer to edge = higher score)
    boundary_dist = compute_boundary_distance(adata)
    # Normalize: max score at boundary, 0 at center
    max_dist = np.percentile(boundary_dist[np.isfinite(boundary_dist)], 95)
    boundary_score = 1 - np.clip(boundary_dist / max_dist, 0, 1)
    qc_df['boundary_score'] = boundary_score
    qc_df['boundary_distance'] = boundary_dist

    # 4. Gene detection rate (low = suspicious)
    if 'n_genes_by_counts' in adata.obs:
        gene_counts = adata.obs['n_genes_by_counts']
    else:
        gene_counts = np.array((adata.X > 0).sum(axis=1)).flatten()

    gene_pct = np.percentile(gene_counts, [5, 95])
    gene_norm = (gene_counts - gene_pct[0]) / (gene_pct[1] - gene_pct[0] + 1e-10)
    qc_df['gene_score'] = 1 - np.clip(gene_norm, 0, 1)  # Inverted: low gene count = high score

    # Composite QC score (weighted)
    qc_df['composite_qc_score'] = (
        0.4 * qc_df['admix_score'] +
        0.3 * qc_df['density_score'] +
        0.2 * qc_df['boundary_score'] +
        0.1 * qc_df['gene_score']
    )

    return qc_df


def create_qc_subsets(adata, qc_df, config=QC_CONFIG):
    """
    Create three analysis subsets:
    1. Baseline: All cells (after hard QC only)
    2. QC-strict: Remove top 10% by composite QC score
    3. Clean: Admixture-filtered (admix_flag = False)

    Returns dict of masks.
    """
    logger.info("Creating QC subsets...")

    n_cells = adata.n_obs

    # Baseline: hard QC only
    baseline_mask = apply_hard_qc_filters(adata, config)

    # QC-strict: remove top X% by composite score
    strict_threshold = np.percentile(qc_df['composite_qc_score'], config['admix_strict_percentile'])
    strict_mask = baseline_mask & (qc_df['composite_qc_score'] < strict_threshold)

    # Clean: admixture-filtered
    if 'admixture_flag' in adata.obs:
        clean_mask = baseline_mask & (~adata.obs['admixture_flag'])
    else:
        # Fall back to admix score threshold
        admix_threshold = np.percentile(qc_df['admix_score'], 70)
        clean_mask = baseline_mask & (qc_df['admix_score'] < admix_threshold)

    subsets = {
        'baseline': baseline_mask,
        'qc_strict': strict_mask,
        'clean': clean_mask,
    }

    logger.info(f"  Baseline: {np.sum(baseline_mask)} cells ({100*np.mean(baseline_mask):.1f}%)")
    logger.info(f"  QC-strict: {np.sum(strict_mask)} cells ({100*np.mean(strict_mask):.1f}%)")
    logger.info(f"  Clean: {np.sum(clean_mask)} cells ({100*np.mean(clean_mask):.1f}%)")

    return subsets


# =============================================================================
# Robustness Analysis
# =============================================================================

def compare_cluster_proportions(adata, subsets, cluster_key='leiden_wnn_0_5'):
    """
    Compare cluster proportions across QC subsets.
    """
    if cluster_key not in adata.obs:
        return None

    results = {}
    for subset_name, mask in subsets.items():
        if np.sum(mask) == 0:
            continue
        proportions = adata.obs.loc[mask, cluster_key].value_counts(normalize=True)
        results[subset_name] = proportions

    # Combine into DataFrame
    df = pd.DataFrame(results).fillna(0)

    # Compute stability metric (correlation between subsets)
    if len(df.columns) >= 2:
        stability = df.corr()
    else:
        stability = None

    return df, stability


def compare_spatial_metrics(adata, subsets, cluster_key='leiden_wnn_0_5'):
    """
    Compare neighborhood enrichment across QC subsets.
    """
    results = {}

    for subset_name, mask in subsets.items():
        if np.sum(mask) < 100:
            continue

        adata_sub = adata[mask].copy()

        if cluster_key not in adata_sub.obs or 'spatial' not in adata_sub.obsm:
            continue

        try:
            # Build spatial graph
            sq.gr.spatial_neighbors(adata_sub, coord_type='generic', n_neighs=6)

            # Neighborhood enrichment
            sq.gr.nhood_enrichment(adata_sub, cluster_key=cluster_key)

            zscore = adata_sub.uns[f'{cluster_key}_nhood_enrichment']['zscore']
            results[subset_name] = zscore
        except Exception as e:
            logger.warning(f"Could not compute spatial metrics for {subset_name}: {e}")

    return results


# =============================================================================
# Visualization
# =============================================================================

def plot_qc_distributions(qc_df, output_path):
    """Plot QC score distributions."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Admixture score
    axes[0, 0].hist(qc_df['admix_score'], bins=50, edgecolor='black', alpha=0.7)
    axes[0, 0].set_xlabel('Admixture Score')
    axes[0, 0].set_title('Admixture Score Distribution')

    # Density score
    axes[0, 1].hist(qc_df['density_score'], bins=50, edgecolor='black', alpha=0.7, color='orange')
    axes[0, 1].set_xlabel('Density Score')
    axes[0, 1].set_title('Cell Density Score')

    # Boundary score
    axes[0, 2].hist(qc_df['boundary_score'], bins=50, edgecolor='black', alpha=0.7, color='green')
    axes[0, 2].set_xlabel('Boundary Score')
    axes[0, 2].set_title('Boundary Distance Score')

    # Composite score
    axes[1, 0].hist(qc_df['composite_qc_score'], bins=50, edgecolor='black', alpha=0.7, color='red')
    axes[1, 0].set_xlabel('Composite QC Score')
    axes[1, 0].set_title('Composite QC Score')

    # Score correlations
    scatter = axes[1, 1].scatter(
        qc_df['admix_score'],
        qc_df['density_score'],
        c=qc_df['composite_qc_score'],
        cmap='RdYlBu_r',
        alpha=0.3,
        s=1
    )
    axes[1, 1].set_xlabel('Admixture Score')
    axes[1, 1].set_ylabel('Density Score')
    axes[1, 1].set_title('Admixture vs Density')
    plt.colorbar(scatter, ax=axes[1, 1])

    # Cumulative distribution
    sorted_scores = np.sort(qc_df['composite_qc_score'])
    cdf = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores)
    axes[1, 2].plot(sorted_scores, cdf)
    axes[1, 2].axhline(0.9, color='red', linestyle='--', label='90th percentile')
    axes[1, 2].set_xlabel('Composite QC Score')
    axes[1, 2].set_ylabel('Cumulative Proportion')
    axes[1, 2].set_title('CDF of QC Scores')
    axes[1, 2].legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_spatial_qc(adata, qc_df, output_path):
    """Plot spatial distribution of QC scores."""
    if 'spatial' not in adata.obsm:
        return

    coords = adata.obsm['spatial']

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Subsample for large datasets
    n = min(30000, adata.n_obs)
    idx = np.random.choice(adata.n_obs, n, replace=False)

    # Composite QC score
    scatter = axes[0, 0].scatter(
        coords[idx, 0], coords[idx, 1],
        c=qc_df['composite_qc_score'].values[idx],
        cmap='RdYlBu_r', s=1, alpha=0.5
    )
    axes[0, 0].set_title('Composite QC Score')
    axes[0, 0].set_aspect('equal')
    plt.colorbar(scatter, ax=axes[0, 0])

    # Density
    scatter = axes[0, 1].scatter(
        coords[idx, 0], coords[idx, 1],
        c=qc_df['cell_density'].values[idx],
        cmap='YlOrRd', s=1, alpha=0.5
    )
    axes[0, 1].set_title('Cell Density')
    axes[0, 1].set_aspect('equal')
    plt.colorbar(scatter, ax=axes[0, 1])

    # Admixture
    scatter = axes[1, 0].scatter(
        coords[idx, 0], coords[idx, 1],
        c=qc_df['admix_score'].values[idx],
        cmap='RdYlBu_r', s=1, alpha=0.5
    )
    axes[1, 0].set_title('Admixture Score')
    axes[1, 0].set_aspect('equal')
    plt.colorbar(scatter, ax=axes[1, 0])

    # Boundary distance
    scatter = axes[1, 1].scatter(
        coords[idx, 0], coords[idx, 1],
        c=np.clip(qc_df['boundary_distance'].values[idx], 0, 500),
        cmap='viridis', s=1, alpha=0.5
    )
    axes[1, 1].set_title('Boundary Distance')
    axes[1, 1].set_aspect('equal')
    plt.colorbar(scatter, ax=axes[1, 1])

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_subset_comparison(proportions_df, output_path):
    """Plot cluster proportion comparison across subsets."""
    if proportions_df is None:
        return

    fig, ax = plt.subplots(figsize=(12, 6))
    proportions_df.plot(kind='bar', ax=ax)
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Proportion')
    ax.set_title('Cluster Proportions by QC Subset')
    ax.legend(title='Subset')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


# =============================================================================
# Main Pipeline
# =============================================================================

def process_sample(sample_path):
    """Process single sample through segmentation QC."""
    sample_id = sample_path.stem.replace('_wnn', '')

    try:
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing: {sample_id}")
        logger.info(f"{'='*60}")

        # Load WNN-integrated data
        adata = sc.read_h5ad(sample_path)
        logger.info(f"  Loaded: {adata.n_obs} cells")

        # Compute QC scores
        qc_df = compute_qc_scores(adata)

        # Add QC scores to adata
        for col in qc_df.columns:
            adata.obs[col] = qc_df[col].values

        # Create QC subsets
        subsets = create_qc_subsets(adata, qc_df)

        # Store subset masks
        for name, mask in subsets.items():
            adata.obs[f'qc_pass_{name}'] = mask

        # Compare cluster proportions
        proportions, stability = compare_cluster_proportions(adata, subsets)

        # Compare spatial metrics
        spatial_metrics = compare_spatial_metrics(adata, subsets)

        # Save results
        out_path = OUTPUT_DIR / f"{sample_id}_qc.h5ad"
        adata.write(out_path)
        logger.info(f"  Saved: {out_path}")

        # Generate figures
        plot_qc_distributions(qc_df, FIGURE_DIR / f"{sample_id}_qc_distributions.png")
        plot_spatial_qc(adata, qc_df, FIGURE_DIR / f"{sample_id}_spatial_qc.png")
        if proportions is not None:
            plot_subset_comparison(proportions, FIGURE_DIR / f"{sample_id}_subset_comparison.png")

        # Summary
        result = {
            'sample_id': sample_id,
            'n_cells': adata.n_obs,
            'n_baseline': int(subsets['baseline'].sum()),
            'n_strict': int(subsets['qc_strict'].sum()),
            'n_clean': int(subsets['clean'].sum()),
            'pct_baseline': float(100 * subsets['baseline'].mean()),
            'pct_strict': float(100 * subsets['qc_strict'].mean()),
            'pct_clean': float(100 * subsets['clean'].mean()),
            'mean_composite_qc': float(qc_df['composite_qc_score'].mean()),
            'median_density': float(qc_df['cell_density'].median()),
        }

        # Add stability metric
        if stability is not None and 'baseline' in stability and 'qc_strict' in stability:
            result['stability_baseline_strict'] = float(stability.loc['baseline', 'qc_strict'])

        del adata
        gc.collect()

        return result

    except Exception as e:
        logger.error(f"  Error processing {sample_id}: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    parser = argparse.ArgumentParser(description='G4X Segmentation QC')
    parser.add_argument('--sample', type=str, help='Process specific sample')
    parser.add_argument('--all', action='store_true', help='Process all')
    args = parser.parse_args()

    print("=" * 70)
    print("G4X Segmentation QC + Sensitivity Analysis")
    print("=" * 70)

    sample_files = sorted(INPUT_DIR.glob("*_wnn.h5ad"))
    logger.info(f"Found {len(sample_files)} WNN-integrated samples")

    if args.sample:
        sample_files = [f for f in sample_files if args.sample in f.stem]

    results = []
    for f in tqdm(sample_files, desc="Processing"):
        result = process_sample(f)
        if result:
            results.append(result)

    # Summary
    if results:
        df = pd.DataFrame(results)
        df.to_csv(OUTPUT_DIR / "segmentation_qc_summary.csv", index=False)
        print("\n" + "=" * 70)
        print("Summary")
        print("=" * 70)
        print(df.to_string(index=False))

    print("\nDone!")


if __name__ == "__main__":
    main()
