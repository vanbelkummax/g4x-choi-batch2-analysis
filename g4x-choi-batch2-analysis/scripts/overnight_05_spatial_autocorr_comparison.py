#!/usr/bin/env python
"""
Spatial Autocorrelation Methods Comparison for G4X Cancer Composite

Compares 6 methods for detecting spatial clustering of cancer cells:
1. Raw Expression - cancer_composite score
2. Moran's I - Global spatial autocorrelation statistic
3. Local Moran's I (LISA) - Local clusters (HH, LL, HL, LH)
4. Getis-Ord Gi* - Hot/cold spot analysis
5. Ripley's K/L - Point pattern at multiple scales
6. Neighborhood Enrichment - Cell type co-localization z-scores

Output: 6-row × 3-column figure comparing methods across E02/F02/G02
"""

import warnings
warnings.filterwarnings('ignore')

import os
import sys
# Disable CUDA for this script (avoid squidpy GPU issues)
os.environ['CUDA_VISIBLE_DEVICES'] = ''
# Force unbuffered output
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import zscore as scipy_zscore
from libpysal.weights import KNN
from esda.moran import Moran, Moran_Local
from esda.getisord import G_Local
import seaborn as sns
from joblib import Parallel, delayed
import multiprocessing

# Use all available cores
N_JOBS = multiprocessing.cpu_count()
print(f"Using {N_JOBS} CPU cores for parallel processing")

# Output paths
OUTPUT_DIR = '/home/user/g4x-choi-batch2-analysis/output/figures'
WINDOWS_DIR = '/mnt/c/Users/User/Desktop/G4X_Overnight_Results'

def build_spatial_weights(coords, k=15):
    """Build k-nearest neighbor spatial weights matrix using libpysal."""
    w = KNN.from_array(coords, k=k)
    w.transform = 'r'  # Row standardize
    return w

def compute_global_morans_i(adata, feature='cancer_composite', k=15):
    """Compute global Moran's I for each sample."""
    results = {}
    for sample in adata.obs['sample_id'].unique():
        mask = adata.obs['sample_id'] == sample
        coords = adata.obsm['spatial'][mask]
        values = adata.obs.loc[mask, feature].values

        w = build_spatial_weights(coords, k=k)
        mi = Moran(values, w)
        results[sample] = {
            'I': mi.I,
            'p_value': mi.p_sim,
            'z_score': mi.z_sim,
            'expected_I': mi.EI
        }
    return results

def compute_local_morans_i(coords, values, k=15):
    """Compute Local Moran's I (LISA) for identifying spatial clusters."""
    w = build_spatial_weights(coords, k=k)
    lisa = Moran_Local(values, w, permutations=99)  # Reduced for speed

    # Classify into HH, LL, HL, LH quadrants
    # q: 1=HH, 2=LH, 3=LL, 4=HL
    quadrant = lisa.q
    sig = lisa.p_sim < 0.05

    # Create labels
    labels = np.full(len(values), 'NS', dtype=object)  # Not significant
    labels[(quadrant == 1) & sig] = 'HH'  # High-High cluster
    labels[(quadrant == 2) & sig] = 'LH'  # Low-High outlier
    labels[(quadrant == 3) & sig] = 'LL'  # Low-Low cluster
    labels[(quadrant == 4) & sig] = 'HL'  # High-Low outlier

    return {
        'Is': lisa.Is,           # Local I values
        'p_sim': lisa.p_sim,     # p-values
        'quadrant': quadrant,
        'labels': labels,
        'z_sim': lisa.z_sim
    }

def compute_getis_ord(coords, values, k=15):
    """Compute Getis-Ord Gi* statistic for hot/cold spots."""
    w = build_spatial_weights(coords, k=k)
    # Use star=True for Gi* (includes self in neighborhood)
    g = G_Local(values, w, star=True, permutations=99)  # Reduced for speed

    # Classify hot/cold spots
    labels = np.full(len(values), 'NS', dtype=object)
    sig = g.p_sim < 0.05
    labels[(g.Zs > 0) & sig] = 'Hot'   # Hot spot
    labels[(g.Zs < 0) & sig] = 'Cold'  # Cold spot

    return {
        'Gs': g.Gs,
        'Zs': g.Zs,
        'p_sim': g.p_sim,
        'labels': labels
    }

def _single_permutation(celltypes_perm, neighbors, type_to_idx, n_types):
    """Single permutation for parallel execution."""
    perm_counts = np.zeros((n_types, n_types))
    for i, cell_type in enumerate(celltypes_perm):
        type_i = type_to_idx[cell_type]
        for j in neighbors[i]:
            type_j = type_to_idx[celltypes_perm[j]]
            perm_counts[type_i, type_j] += 1
    return perm_counts

def compute_nhood_enrichment(coords, celltypes, k=15, n_perms=50):
    """
    Compute neighborhood enrichment z-scores (parallelized implementation).
    Tests whether cell types are spatially enriched or depleted as neighbors.
    """
    # Get unique cell types
    unique_types = np.unique(celltypes)
    n_types = len(unique_types)
    type_to_idx = {t: i for i, t in enumerate(unique_types)}

    # Build KNN graph
    tree = cKDTree(coords)
    _, neighbors = tree.query(coords, k=k+1)  # +1 because first neighbor is self
    neighbors = neighbors[:, 1:]  # Remove self

    # Count observed neighbor pairs (vectorized)
    celltypes_idx = np.array([type_to_idx[t] for t in celltypes])
    observed = np.zeros((n_types, n_types))
    for i in range(len(celltypes)):
        type_i = celltypes_idx[i]
        neighbor_types = celltypes_idx[neighbors[i]]
        for type_j in neighbor_types:
            observed[type_i, type_j] += 1

    # Parallel permutation test
    permutations = [np.random.permutation(celltypes) for _ in range(n_perms)]
    expected_list = Parallel(n_jobs=N_JOBS)(
        delayed(_single_permutation)(perm, neighbors, type_to_idx, n_types)
        for perm in permutations
    )

    expected_list = np.array(expected_list)
    expected_mean = np.mean(expected_list, axis=0)
    expected_std = np.std(expected_list, axis=0)

    # Compute z-scores
    with np.errstate(divide='ignore', invalid='ignore'):
        zscore_matrix = (observed - expected_mean) / expected_std
        zscore_matrix = np.nan_to_num(zscore_matrix, nan=0, posinf=0, neginf=0)

    return zscore_matrix, list(unique_types)

def compute_ripleys_k(coords, values, radii=None, threshold_percentile=75, max_points=10000):
    """
    Compute Ripley's K/L function for high-value cells.
    Tests clustering of cells with high cancer_composite scores.
    Subsamples if > max_points to avoid O(n²) blowup.
    """
    if radii is None:
        # Auto-determine radii based on data extent
        extent = np.max(coords, axis=0) - np.min(coords, axis=0)
        max_r = np.min(extent) / 4
        radii = np.linspace(0, max_r, 30)  # Reduced from 50

    # Get high-value cells (top 25%)
    threshold = np.percentile(values, threshold_percentile)
    high_mask = values >= threshold
    high_coords = coords[high_mask]

    n = len(high_coords)
    if n < 10:
        return {'radii': radii, 'K': np.zeros_like(radii), 'L': np.zeros_like(radii),
                'n_high': n, 'threshold': threshold}

    # Subsample if too many points (O(n²) is expensive)
    if n > max_points:
        idx = np.random.choice(n, max_points, replace=False)
        high_coords = high_coords[idx]
        n = max_points

    # Build KD-tree for efficient neighbor counting
    tree = cKDTree(high_coords)

    # Compute K(r) = (A/n^2) * sum of pairs within r
    area = np.prod(np.max(coords, axis=0) - np.min(coords, axis=0))

    K = np.zeros(len(radii))
    for i, r in enumerate(radii):
        if r == 0:
            K[i] = 0
        else:
            # Count pairs within radius r
            pairs = tree.query_pairs(r)
            K[i] = (area / (n * (n-1))) * 2 * len(pairs) if n > 1 else 0

    # L(r) = sqrt(K(r)/pi) - r (Besag's transformation)
    # Under CSR, L(r) ≈ 0
    L = np.sqrt(K / np.pi) - radii

    return {
        'radii': radii,
        'K': K,
        'L': L,
        'n_high': n,
        'threshold': threshold
    }

def create_comparison_figure(adata, k=15):
    """Create 6-row × 3-column comparison figure."""

    samples = ['E02', 'F02', 'G02']
    stages = ['Normal', 'Metaplasia', 'Cancer']

    fig, axes = plt.subplots(6, 3, figsize=(18, 24))
    fig.suptitle('Spatial Autocorrelation Methods Comparison: cancer_composite',
                 fontsize=16, fontweight='bold', y=0.995)

    # Color maps
    cmap_raw = 'RdBu_r'
    cmap_lisa = {'HH': '#d7191c', 'LL': '#2c7bb6', 'HL': '#fdae61',
                 'LH': '#abd9e9', 'NS': '#eeeeee'}
    cmap_hotcold = {'Hot': '#d7191c', 'Cold': '#2c7bb6', 'NS': '#eeeeee'}

    # Storage for global statistics
    global_stats = []

    for col, (sample, stage) in enumerate(zip(samples, stages)):
        print(f"\nProcessing {sample} ({stage})...")
        mask = adata.obs['sample_id'] == sample
        coords = adata.obsm['spatial'][mask]
        values = adata.obs.loc[mask, 'cancer_composite'].values

        # Subsample for speed if needed (keep all for accuracy)
        n_cells = np.sum(mask)
        print(f"  {n_cells:,} cells")

        # === Row 0: Raw Expression ===
        ax = axes[0, col]
        sc = ax.scatter(coords[:, 0], coords[:, 1],
                       c=values, cmap=cmap_raw, s=0.5, alpha=0.7,
                       vmin=-1, vmax=2)
        ax.set_title(f'{sample}: {stage}\n(n={n_cells:,})', fontsize=12)
        ax.set_aspect('equal')
        ax.axis('off')
        if col == 2:
            plt.colorbar(sc, ax=ax, label='cancer_composite', shrink=0.7)
        if col == 0:
            ax.set_ylabel('Raw Expression', fontsize=12, fontweight='bold')

        # === Row 1: Global Moran's I (just display statistic) ===
        ax = axes[1, col]
        w = build_spatial_weights(coords, k=k)
        mi = Moran(values, w, permutations=99)  # Reduced for speed

        # Show the values as a bar/text
        ax.text(0.5, 0.6, f"Moran's I = {mi.I:.4f}",
                ha='center', va='center', fontsize=16, fontweight='bold',
                transform=ax.transAxes)
        ax.text(0.5, 0.4, f"p-value = {mi.p_sim:.4f}",
                ha='center', va='center', fontsize=14,
                transform=ax.transAxes)
        ax.text(0.5, 0.25, f"z-score = {mi.z_sim:.2f}",
                ha='center', va='center', fontsize=12,
                transform=ax.transAxes)

        # Color code significance
        bg_color = '#ffcccc' if mi.p_sim < 0.05 else '#cccccc'
        ax.set_facecolor(bg_color)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        if col == 0:
            ax.set_ylabel("Global Moran's I", fontsize=12, fontweight='bold')

        global_stats.append({
            'sample': sample, 'stage': stage,
            'morans_I': mi.I, 'morans_p': mi.p_sim, 'morans_z': mi.z_sim
        })
        print(f"  Moran's I: {mi.I:.4f} (p={mi.p_sim:.4f})")

        # === Row 2: Local Moran's I (LISA) ===
        ax = axes[2, col]
        print("  Computing Local Moran's I...")
        lisa_result = compute_local_morans_i(coords, values, k=k)

        # Plot by quadrant
        for label, color in cmap_lisa.items():
            mask_label = lisa_result['labels'] == label
            if np.any(mask_label):
                ax.scatter(coords[mask_label, 0], coords[mask_label, 1],
                          c=color, s=0.5, alpha=0.8, label=label)

        # Count significant
        n_hh = np.sum(lisa_result['labels'] == 'HH')
        n_ll = np.sum(lisa_result['labels'] == 'LL')
        ax.set_title(f'HH={n_hh:,}, LL={n_ll:,}', fontsize=10)
        ax.set_aspect('equal')
        ax.axis('off')
        if col == 2:
            ax.legend(loc='upper right', markerscale=5, fontsize=8)
        if col == 0:
            ax.set_ylabel('Local Moran (LISA)', fontsize=12, fontweight='bold')

        print(f"  LISA: HH={n_hh}, LL={n_ll}")

        # === Row 3: Getis-Ord Gi* ===
        ax = axes[3, col]
        print("  Computing Getis-Ord Gi*...")
        gi_result = compute_getis_ord(coords, values, k=k)

        # Plot hot/cold spots
        for label, color in cmap_hotcold.items():
            mask_label = gi_result['labels'] == label
            if np.any(mask_label):
                ax.scatter(coords[mask_label, 0], coords[mask_label, 1],
                          c=color, s=0.5, alpha=0.8, label=label)

        n_hot = np.sum(gi_result['labels'] == 'Hot')
        n_cold = np.sum(gi_result['labels'] == 'Cold')
        ax.set_title(f'Hot={n_hot:,}, Cold={n_cold:,}', fontsize=10)
        ax.set_aspect('equal')
        ax.axis('off')
        if col == 2:
            ax.legend(loc='upper right', markerscale=5, fontsize=8)
        if col == 0:
            ax.set_ylabel('Getis-Ord Gi*', fontsize=12, fontweight='bold')

        print(f"  Gi*: Hot={n_hot}, Cold={n_cold}")

        # === Row 4: Ripley's L function ===
        ax = axes[4, col]
        print("  Computing Ripley's K/L...")
        ripley_result = compute_ripleys_k(coords, values)

        ax.plot(ripley_result['radii'], ripley_result['L'], 'b-', linewidth=2)
        ax.axhline(y=0, color='k', linestyle='--', alpha=0.5, label='CSR')
        ax.fill_between(ripley_result['radii'], ripley_result['L'], 0,
                       where=ripley_result['L'] > 0, alpha=0.3, color='red',
                       label='Clustered')
        ax.fill_between(ripley_result['radii'], ripley_result['L'], 0,
                       where=ripley_result['L'] < 0, alpha=0.3, color='blue',
                       label='Dispersed')
        ax.set_xlabel('Distance (r)', fontsize=10)
        ax.set_ylabel('L(r) - r', fontsize=10)
        ax.set_title(f"High cells (top 25%): n={ripley_result['n_high']:,}", fontsize=10)
        if col == 2:
            ax.legend(loc='upper right', fontsize=8)
        if col == 0:
            ax.set_ylabel("Ripley's L Function", fontsize=12, fontweight='bold')

        # Store max L for comparison
        max_L = np.max(ripley_result['L'])
        print(f"  Ripley L max: {max_L:.2f}")

    # === Row 5: Neighborhood Enrichment (pure Python implementation) ===
    print("\nComputing neighborhood enrichment (pure Python)...")

    for col, sample in enumerate(samples):
        ax = axes[5, col]
        mask = adata.obs['sample_id'] == sample
        coords = adata.obsm['spatial'][mask]
        celltypes = adata.obs.loc[mask, 'celltype_final'].values

        # Compute neighborhood enrichment manually
        zscore_matrix, celltype_labels = compute_nhood_enrichment(coords, celltypes, k=k)

        # Plot heatmap
        im = ax.imshow(zscore_matrix, cmap='RdBu_r', vmin=-10, vmax=10, aspect='auto')

        n_types = len(celltype_labels)

        # Only label if not too many
        if n_types <= 15:
            ax.set_xticks(range(n_types))
            ax.set_yticks(range(n_types))
            ax.set_xticklabels(celltype_labels, rotation=90, fontsize=6)
            ax.set_yticklabels(celltype_labels, fontsize=6)
        else:
            ax.set_xticks([])
            ax.set_yticks([])

        ax.set_title(f'{n_types} cell types', fontsize=10)
        if col == 2:
            plt.colorbar(im, ax=ax, label='z-score', shrink=0.7)
        if col == 0:
            ax.set_ylabel('Nhood Enrichment', fontsize=12, fontweight='bold')

    plt.tight_layout()

    # Save
    import os
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fig.savefig(f'{OUTPUT_DIR}/fig_spatial_autocorr_comparison.png', dpi=200, bbox_inches='tight')
    fig.savefig(f'{WINDOWS_DIR}/fig_spatial_autocorr_comparison.png', dpi=200, bbox_inches='tight')
    print(f"\nSaved to {OUTPUT_DIR}/fig_spatial_autocorr_comparison.png")
    print(f"Saved to {WINDOWS_DIR}/fig_spatial_autocorr_comparison.png")

    plt.close()

    # Save global statistics
    stats_df = pd.DataFrame(global_stats)
    stats_df.to_csv(f'{OUTPUT_DIR}/spatial_autocorr_stats.csv', index=False)
    print(f"\nGlobal Moran's I statistics:")
    print(stats_df.to_string(index=False))

    return stats_df

def main():
    print("Loading data...")
    adata = sc.read_h5ad('output/overnight_annotated_v2.h5ad')
    print(f"Loaded {adata.n_obs:,} cells, {adata.n_vars} genes")

    # Run comparison
    stats = create_comparison_figure(adata, k=15)

    print("\n" + "="*60)
    print("SPATIAL AUTOCORRELATION COMPARISON COMPLETE")
    print("="*60)
    print(f"\nKey findings:")
    for _, row in stats.iterrows():
        sig = "***" if row['morans_p'] < 0.001 else "**" if row['morans_p'] < 0.01 else "*" if row['morans_p'] < 0.05 else ""
        print(f"  {row['sample']} ({row['stage']}): I={row['morans_I']:.4f} {sig}")

    print("\nExpected: Higher Moran's I in Cancer (G02) = stronger spatial clustering")

if __name__ == '__main__':
    main()
