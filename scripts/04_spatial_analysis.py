#!/usr/bin/env python3
"""
G4X Choi Batch 2 - Spatial Analysis
====================================

Analyze spatial patterns and cell-cell interactions.
- Neighborhood enrichment
- Spatial niche identification
- Cell-cell interaction analysis
"""

import os
import sys
import gzip
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from scipy.spatial import cKDTree
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Configuration
DATA_ROOT = Path('/mnt/x/Choi_Batch_2_Tuesday')
OUTPUT_DIR = Path('/home/user/g4x-choi-batch2-analysis/results')
FIG_DIR = OUTPUT_DIR / 'figures'
TABLE_DIR = OUTPUT_DIR / 'tables'

# Cell type colors
CELL_TYPE_COLORS = {
    'Epithelial': '#1f77b4',
    'CD8_T': '#ff7f0e',
    'CD4_T': '#2ca02c',
    'Treg': '#d62728',
    'T_cell': '#9467bd',
    'B_cell': '#8c564b',
    'Macrophage': '#e377c2',
    'DC': '#7f7f7f',
    'Stromal': '#bcbd22',
    'Endothelial': '#17becf',
    'Immune_other': '#ff9896',
    'Other': '#c7c7c7'
}

def load_sample_data(lane_idx, sample_prefix):
    """Load cell metadata with coordinates and cell types."""
    lane = f'L00{lane_idx + 1}'
    sample = f'{sample_prefix}{lane_idx + 1}'
    lane_dir = list(DATA_ROOT.glob(f'g4-028-083-FC1-{lane}_*'))[0]
    sample_path = lane_dir / sample / 'single_cell_data' / 'cell_metadata.csv.gz'

    if not sample_path.exists():
        return None

    with gzip.open(sample_path, 'rt') as f:
        df = pd.read_csv(f)

    return df

def gate_cell_types_quick(df, quantile_thresh=0.8):
    """Quick cell type gating (same as annotation script)."""

    thresholds = {}
    protein_cols = ['PD1', 'CD8', 'PDL1', 'FOXP3', 'KI67', 'CD45', 'PanCK',
                    'CD20', 'CD4', 'CD11c', 'HLA-DR', 'CD68', 'CD31', 'CD3', 'aSMA']

    for marker in protein_cols:
        col = f'{marker}_intensity_mean'
        if col in df.columns:
            thresholds[marker] = df[col].quantile(quantile_thresh)

    def get_val(marker):
        col = f'{marker}_intensity_mean'
        return df[col].values if col in df.columns else np.zeros(len(df))

    def get_t(marker):
        return thresholds.get(marker, 0)

    cell_types = np.full(len(df), 'Other', dtype=object)

    # Gating
    cd45 = get_val('CD45') > get_t('CD45')
    cd3 = get_val('CD3') > get_t('CD3')
    cd4 = get_val('CD4') > get_t('CD4')
    cd8 = get_val('CD8') > get_t('CD8')
    foxp3 = get_val('FOXP3') > get_t('FOXP3')
    cd20 = get_val('CD20') > get_t('CD20')
    cd68 = get_val('CD68') > get_t('CD68')
    cd11c = get_val('CD11c') > get_t('CD11c')
    hladr = get_val('HLA-DR') > get_t('HLA-DR')
    panck = get_val('PanCK') > get_t('PanCK')
    asma = get_val('aSMA') > get_t('aSMA')
    cd31 = get_val('CD31') > get_t('CD31')

    cell_types[panck & ~cd45] = 'Epithelial'
    cell_types[asma & ~panck & ~cd45] = 'Stromal'
    cell_types[cd31 & ~panck & ~cd45] = 'Endothelial'
    cell_types[cd45 & cd3 & cd8 & ~cd4] = 'CD8_T'
    cell_types[cd45 & cd3 & cd4 & ~cd8 & foxp3] = 'Treg'
    cell_types[cd45 & cd3 & cd4 & ~cd8 & ~foxp3] = 'CD4_T'
    cell_types[cd45 & cd3 & ~cd4 & ~cd8] = 'T_cell'
    cell_types[cd45 & cd20 & ~cd3] = 'B_cell'
    cell_types[cd45 & cd68 & ~cd3 & ~cd20] = 'Macrophage'
    cell_types[cd45 & cd11c & hladr & ~cd3 & ~cd20 & ~cd68] = 'DC'
    cell_types[cd45 & (cell_types == 'Other')] = 'Immune_other'

    return cell_types

def compute_neighborhood_enrichment(coords, cell_types, radius=50, n_permutations=100):
    """
    Compute neighborhood enrichment scores.

    For each cell type pair, calculate if they co-occur more or less
    than expected by chance.
    """

    # Build KD-tree for fast neighbor lookup
    tree = cKDTree(coords)

    # Get unique cell types (excluding rare ones)
    type_counts = pd.Series(cell_types).value_counts()
    valid_types = type_counts[type_counts > 100].index.tolist()

    n_types = len(valid_types)
    type_to_idx = {t: i for i, t in enumerate(valid_types)}

    # Count observed co-occurrences
    observed = np.zeros((n_types, n_types))

    for i, (coord, ct) in enumerate(zip(coords, cell_types)):
        if ct not in type_to_idx:
            continue

        # Find neighbors within radius
        neighbors = tree.query_ball_point(coord, radius)

        for j in neighbors:
            if i == j:
                continue
            neighbor_type = cell_types[j]
            if neighbor_type not in type_to_idx:
                continue

            observed[type_to_idx[ct], type_to_idx[neighbor_type]] += 1

    # Permutation test
    expected = np.zeros((n_types, n_types))
    expected_std = np.zeros((n_types, n_types))

    perm_counts = []
    for _ in range(n_permutations):
        # Shuffle cell type labels
        shuffled_types = np.random.permutation(cell_types)
        perm_observed = np.zeros((n_types, n_types))

        for i, (coord, ct) in enumerate(zip(coords, shuffled_types)):
            if ct not in type_to_idx:
                continue

            neighbors = tree.query_ball_point(coord, radius)
            for j in neighbors:
                if i == j:
                    continue
                neighbor_type = shuffled_types[j]
                if neighbor_type not in type_to_idx:
                    continue
                perm_observed[type_to_idx[ct], type_to_idx[neighbor_type]] += 1

        perm_counts.append(perm_observed)

    perm_counts = np.array(perm_counts)
    expected = perm_counts.mean(axis=0)
    expected_std = perm_counts.std(axis=0)

    # Compute z-scores
    with np.errstate(divide='ignore', invalid='ignore'):
        zscore = (observed - expected) / (expected_std + 1e-10)
        zscore = np.nan_to_num(zscore, nan=0, posinf=5, neginf=-5)

    return zscore, valid_types, observed, expected

def analyze_sample_spatial(lane_idx, sample_prefix, subsample=20000):
    """Analyze spatial patterns for a single sample."""

    df = load_sample_data(lane_idx, sample_prefix)
    if df is None:
        return None

    sample_id = f'{sample_prefix}{lane_idx + 1}'

    # Add cell types
    cell_types = gate_cell_types_quick(df)
    df['cell_type'] = cell_types

    # Subsample for speed
    if len(df) > subsample:
        df = df.sample(subsample, random_state=42)

    # Get coordinates
    coords = df[['cell_x', 'cell_y']].values
    cell_types = df['cell_type'].values

    # Compute neighborhood enrichment
    try:
        zscore, valid_types, observed, expected = compute_neighborhood_enrichment(
            coords, cell_types, radius=30, n_permutations=50
        )
    except Exception as e:
        print(f"    Warning: {e}")
        return None

    return {
        'sample_id': sample_id,
        'zscore': zscore,
        'cell_types': valid_types,
        'observed': observed,
        'expected': expected,
        'n_cells': len(df)
    }

def plot_neighborhood_heatmap(zscore_mean, cell_types, output_dir):
    """Plot neighborhood enrichment heatmap."""

    fig, ax = plt.subplots(figsize=(10, 8))

    im = ax.imshow(zscore_mean, cmap='RdBu_r', vmin=-3, vmax=3)
    ax.set_xticks(range(len(cell_types)))
    ax.set_xticklabels(cell_types, rotation=45, ha='right')
    ax.set_yticks(range(len(cell_types)))
    ax.set_yticklabels(cell_types)
    ax.set_xlabel('Neighbor Cell Type')
    ax.set_ylabel('Central Cell Type')
    ax.set_title('Neighborhood Enrichment (Mean Z-score)')
    plt.colorbar(im, ax=ax, label='Z-score')

    # Add text annotations
    for i in range(len(cell_types)):
        for j in range(len(cell_types)):
            text = f'{zscore_mean[i,j]:.1f}'
            color = 'white' if abs(zscore_mean[i,j]) > 1.5 else 'black'
            ax.text(j, i, text, ha='center', va='center', color=color, fontsize=8)

    plt.tight_layout()
    plt.savefig(output_dir / 'neighborhood_enrichment.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'neighborhood_enrichment.pdf', bbox_inches='tight')
    plt.close()

    print(f"Saved neighborhood enrichment to {output_dir / 'neighborhood_enrichment.png'}")

def plot_spatial_maps(sample_results, output_dir):
    """Plot spatial distribution for selected samples."""

    # Select diverse samples
    n_show = min(4, len(sample_results))
    selected = sample_results[:n_show]

    fig, axes = plt.subplots(2, 2, figsize=(14, 14))
    axes = axes.flatten()

    for idx, result in enumerate(selected):
        ax = axes[idx]

        sample_id = result['sample_id']

        # Load full data for visualization
        lane_idx = int(sample_id[-1]) - 1
        sample_prefix = sample_id[:-1]
        df = load_sample_data(lane_idx, sample_prefix)

        if df is None:
            continue

        # Add cell types
        cell_types = gate_cell_types_quick(df)
        df['cell_type'] = cell_types

        # Subsample for plotting
        if len(df) > 10000:
            df = df.sample(10000, random_state=42)

        # Plot
        for ct in CELL_TYPE_COLORS.keys():
            ct_data = df[df['cell_type'] == ct]
            if len(ct_data) > 0:
                ax.scatter(ct_data['cell_x'], ct_data['cell_y'],
                          c=CELL_TYPE_COLORS[ct], s=1, alpha=0.5, label=ct)

        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_title(f'Sample {sample_id} (n={result["n_cells"]:,})')
        ax.set_aspect('equal')

    # Add legend to last plot
    axes[-1].legend(bbox_to_anchor=(1.02, 1), loc='upper left', markerscale=5, fontsize=8)

    plt.tight_layout()
    plt.savefig(output_dir / 'spatial_maps.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'spatial_maps.pdf', bbox_inches='tight')
    plt.close()

    print(f"Saved spatial maps to {output_dir / 'spatial_maps.png'}")

def main():
    """Run spatial analysis."""
    print("=" * 60)
    print("G4X Choi Batch 2 - Spatial Analysis")
    print("=" * 60)

    # Process samples
    results = []

    print("\nComputing neighborhood enrichment...")
    for lane_idx in range(4):
        for sample_prefix in ['A0', 'B0', 'C0', 'D0', 'E0', 'F0', 'G0', 'H0']:
            sample_id = f'{sample_prefix}{lane_idx + 1}'
            print(f"  Processing {sample_id}...", end=' ')

            try:
                result = analyze_sample_spatial(lane_idx, sample_prefix)
                if result is not None:
                    results.append(result)
                    print(f"done ({result['n_cells']:,} cells)")
                else:
                    print("skipped")
            except Exception as e:
                print(f"error: {e}")

    print(f"\nProcessed {len(results)} samples")

    # Average neighborhood enrichment across samples
    if len(results) > 0:
        # Find common cell types
        all_types = set()
        for r in results:
            all_types.update(r['cell_types'])

        # Create aligned z-score matrices
        common_types = sorted(all_types)
        n_types = len(common_types)
        type_to_idx = {t: i for i, t in enumerate(common_types)}

        aligned_zscores = []
        for r in results:
            aligned = np.zeros((n_types, n_types))
            for i, ct1 in enumerate(r['cell_types']):
                for j, ct2 in enumerate(r['cell_types']):
                    if ct1 in type_to_idx and ct2 in type_to_idx:
                        aligned[type_to_idx[ct1], type_to_idx[ct2]] = r['zscore'][i, j]
            aligned_zscores.append(aligned)

        zscore_mean = np.mean(aligned_zscores, axis=0)
        zscore_std = np.std(aligned_zscores, axis=0)

        # Save results
        zscore_df = pd.DataFrame(zscore_mean, index=common_types, columns=common_types)
        zscore_df.to_csv(TABLE_DIR / 'neighborhood_enrichment_mean.csv')

        # Generate plots
        print("\n" + "-" * 40)
        print("Generating plots...")

        plot_neighborhood_heatmap(zscore_mean, common_types, FIG_DIR)
        plot_spatial_maps(results, FIG_DIR)

        # Summary statistics
        print("\n" + "=" * 60)
        print("SPATIAL ANALYSIS SUMMARY")
        print("=" * 60)

        print("\nSTRONGEST POSITIVE INTERACTIONS (co-localization):")
        for i, ct1 in enumerate(common_types):
            for j, ct2 in enumerate(common_types):
                if i < j and zscore_mean[i, j] > 2:
                    print(f"  {ct1} <-> {ct2}: z = {zscore_mean[i,j]:.2f}")

        print("\nSTRONGEST NEGATIVE INTERACTIONS (avoidance):")
        for i, ct1 in enumerate(common_types):
            for j, ct2 in enumerate(common_types):
                if i < j and zscore_mean[i, j] < -2:
                    print(f"  {ct1} <-> {ct2}: z = {zscore_mean[i,j]:.2f}")

        print("\nSELF-AGGREGATION (homotypic clustering):")
        for i, ct in enumerate(common_types):
            if zscore_mean[i, i] > 2:
                print(f"  {ct}: z = {zscore_mean[i,i]:.2f} (clusters together)")
            elif zscore_mean[i, i] < -1:
                print(f"  {ct}: z = {zscore_mean[i,i]:.2f} (dispersed)")

    print("\n" + "=" * 60)
    print("Spatial Analysis Complete!")
    print(f"Results saved to: {TABLE_DIR}")
    print("=" * 60)

if __name__ == '__main__':
    main()
