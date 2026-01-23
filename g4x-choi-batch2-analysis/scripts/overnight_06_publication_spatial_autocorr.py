#!/usr/bin/env python
"""
Publication-Quality Spatial Autocorrelation Analysis
G4X Gastric Cancer Progression: Normal → Metaplasia → Cancer

Generates:
1. 5-row × 3-column comparison figure (excludes Ripley's L)
2. Publication-quality results table
3. Summary statistics CSV

Methods: Raw Expression, Moran's I, LISA, Getis-Ord Gi*, Neighborhood Enrichment
"""

import warnings
warnings.filterwarnings('ignore')

import os
import sys
os.environ['CUDA_VISIBLE_DEVICES'] = ''
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import cKDTree
from libpysal.weights import KNN
from esda.moran import Moran, Moran_Local
from esda.getisord import G_Local
from joblib import Parallel, delayed
import multiprocessing

N_JOBS = multiprocessing.cpu_count()
print(f"Using {N_JOBS} CPU cores")

# Output
OUTPUT_DIR = '/mnt/c/Users/User/Desktop/G4X_Spatial_Autocorrelation'
os.makedirs(OUTPUT_DIR, exist_ok=True)

def build_spatial_weights(coords, k=15):
    w = KNN.from_array(coords, k=k)
    w.transform = 'r'
    return w

def compute_local_morans_i(coords, values, k=15):
    w = build_spatial_weights(coords, k=k)
    lisa = Moran_Local(values, w, permutations=99)
    quadrant = lisa.q
    sig = lisa.p_sim < 0.05
    labels = np.full(len(values), 'NS', dtype=object)
    labels[(quadrant == 1) & sig] = 'HH'
    labels[(quadrant == 2) & sig] = 'LH'
    labels[(quadrant == 3) & sig] = 'LL'
    labels[(quadrant == 4) & sig] = 'HL'
    return {'labels': labels, 'Is': lisa.Is, 'p_sim': lisa.p_sim}

def compute_getis_ord(coords, values, k=15):
    w = build_spatial_weights(coords, k=k)
    g = G_Local(values, w, star=True, permutations=99)
    labels = np.full(len(values), 'NS', dtype=object)
    sig = g.p_sim < 0.05
    labels[(g.Zs > 0) & sig] = 'Hot'
    labels[(g.Zs < 0) & sig] = 'Cold'
    return {'labels': labels, 'Zs': g.Zs, 'p_sim': g.p_sim}

def _single_permutation(celltypes_perm, neighbors, type_to_idx, n_types):
    perm_counts = np.zeros((n_types, n_types))
    for i, cell_type in enumerate(celltypes_perm):
        type_i = type_to_idx[cell_type]
        for j in neighbors[i]:
            type_j = type_to_idx[celltypes_perm[j]]
            perm_counts[type_i, type_j] += 1
    return perm_counts

def compute_nhood_enrichment(coords, celltypes, k=15, n_perms=50):
    unique_types = np.unique(celltypes)
    n_types = len(unique_types)
    type_to_idx = {t: i for i, t in enumerate(unique_types)}

    tree = cKDTree(coords)
    _, neighbors = tree.query(coords, k=k+1)
    neighbors = neighbors[:, 1:]

    celltypes_idx = np.array([type_to_idx[t] for t in celltypes])
    observed = np.zeros((n_types, n_types))
    for i in range(len(celltypes)):
        type_i = celltypes_idx[i]
        for type_j in celltypes_idx[neighbors[i]]:
            observed[type_i, type_j] += 1

    permutations = [np.random.permutation(celltypes) for _ in range(n_perms)]
    expected_list = Parallel(n_jobs=N_JOBS)(
        delayed(_single_permutation)(perm, neighbors, type_to_idx, n_types)
        for perm in permutations
    )

    expected_list = np.array(expected_list)
    expected_mean = np.mean(expected_list, axis=0)
    expected_std = np.std(expected_list, axis=0)

    with np.errstate(divide='ignore', invalid='ignore'):
        zscore_matrix = (observed - expected_mean) / expected_std
        zscore_matrix = np.nan_to_num(zscore_matrix, nan=0, posinf=0, neginf=0)

    return zscore_matrix, list(unique_types)

def create_publication_figure(adata, k=15):
    """Create 5-row × 3-column publication figure (no Ripley's L)."""

    samples = ['E02', 'F02', 'G02']
    stages = ['Normal', 'Metaplasia', 'Cancer']

    fig, axes = plt.subplots(5, 3, figsize=(14, 18))

    cmap_raw = 'RdBu_r'
    cmap_lisa = {'HH': '#d7191c', 'LL': '#2c7bb6', 'HL': '#fdae61', 'LH': '#abd9e9', 'NS': '#eeeeee'}
    cmap_hotcold = {'Hot': '#d7191c', 'Cold': '#2c7bb6', 'NS': '#eeeeee'}

    all_stats = []

    for col, (sample, stage) in enumerate(zip(samples, stages)):
        print(f"\nProcessing {sample} ({stage})...")
        mask = adata.obs['sample_id'] == sample
        coords = adata.obsm['spatial'][mask]
        values = adata.obs.loc[mask, 'cancer_composite'].values
        celltypes = adata.obs.loc[mask, 'celltype_final'].values
        n_cells = np.sum(mask)

        # Row 0: Raw Expression
        ax = axes[0, col]
        sc_plot = ax.scatter(coords[:, 0], coords[:, 1], c=values, cmap=cmap_raw,
                            s=0.3, alpha=0.7, vmin=-1, vmax=2, rasterized=True)
        ax.set_title(f'{sample}: {stage}\n(n={n_cells:,})', fontsize=11, fontweight='bold')
        ax.set_aspect('equal')
        ax.axis('off')
        if col == 2:
            cbar = plt.colorbar(sc_plot, ax=ax, shrink=0.6, pad=0.02)
            cbar.set_label('Cancer Composite', fontsize=9)
        if col == 0:
            ax.text(-0.15, 0.5, 'A. Raw\nExpression', transform=ax.transAxes,
                   fontsize=11, fontweight='bold', va='center', ha='right')

        # Row 1: Global Moran's I
        ax = axes[1, col]
        w = build_spatial_weights(coords, k=k)
        mi = Moran(values, w, permutations=99)

        ax.text(0.5, 0.65, f"I = {mi.I:.3f}", ha='center', va='center',
               fontsize=20, fontweight='bold', transform=ax.transAxes)
        ax.text(0.5, 0.40, f"z = {mi.z_sim:.1f}", ha='center', va='center',
               fontsize=14, transform=ax.transAxes)
        sig_text = "p < 0.01 ***" if mi.p_sim < 0.01 else f"p = {mi.p_sim:.3f}"
        ax.text(0.5, 0.20, sig_text, ha='center', va='center',
               fontsize=12, transform=ax.transAxes)

        bg_color = '#ffe6e6' if mi.p_sim < 0.05 else '#f0f0f0'
        ax.set_facecolor(bg_color)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        for spine in ax.spines.values(): spine.set_visible(False)
        ax.set_xticks([]); ax.set_yticks([])
        if col == 0:
            ax.text(-0.15, 0.5, "B. Global\nMoran's I", transform=ax.transAxes,
                   fontsize=11, fontweight='bold', va='center', ha='right')

        all_stats.append({
            'Sample': sample, 'Stage': stage, 'n_cells': n_cells,
            'Morans_I': mi.I, 'Morans_z': mi.z_sim, 'Morans_p': mi.p_sim
        })

        # Row 2: Local Moran's I (LISA)
        ax = axes[2, col]
        print("  Computing LISA...")
        lisa = compute_local_morans_i(coords, values, k=k)

        for label in ['NS', 'LL', 'LH', 'HL', 'HH']:  # Plot NS first
            mask_label = lisa['labels'] == label
            if np.any(mask_label):
                ax.scatter(coords[mask_label, 0], coords[mask_label, 1],
                          c=cmap_lisa[label], s=0.3, alpha=0.8, rasterized=True)

        n_hh = np.sum(lisa['labels'] == 'HH')
        n_ll = np.sum(lisa['labels'] == 'LL')
        ax.set_title(f'HH: {n_hh:,} | LL: {n_ll:,}', fontsize=9)
        ax.set_aspect('equal')
        ax.axis('off')
        if col == 0:
            ax.text(-0.15, 0.5, 'C. Local\nMoran (LISA)', transform=ax.transAxes,
                   fontsize=11, fontweight='bold', va='center', ha='right')

        all_stats[-1]['LISA_HH'] = n_hh
        all_stats[-1]['LISA_LL'] = n_ll

        # Row 3: Getis-Ord Gi*
        ax = axes[3, col]
        print("  Computing Gi*...")
        gi = compute_getis_ord(coords, values, k=k)

        for label in ['NS', 'Cold', 'Hot']:
            mask_label = gi['labels'] == label
            if np.any(mask_label):
                ax.scatter(coords[mask_label, 0], coords[mask_label, 1],
                          c=cmap_hotcold[label], s=0.3, alpha=0.8, rasterized=True)

        n_hot = np.sum(gi['labels'] == 'Hot')
        n_cold = np.sum(gi['labels'] == 'Cold')
        ax.set_title(f'Hot: {n_hot:,} | Cold: {n_cold:,}', fontsize=9)
        ax.set_aspect('equal')
        ax.axis('off')
        if col == 0:
            ax.text(-0.15, 0.5, 'D. Getis-Ord\nGi* Hot/Cold', transform=ax.transAxes,
                   fontsize=11, fontweight='bold', va='center', ha='right')

        all_stats[-1]['Gi_Hot'] = n_hot
        all_stats[-1]['Gi_Cold'] = n_cold

    # Row 4: Neighborhood Enrichment
    print("\nComputing neighborhood enrichment...")
    for col, sample in enumerate(samples):
        ax = axes[4, col]
        mask = adata.obs['sample_id'] == sample
        coords = adata.obsm['spatial'][mask]
        celltypes = adata.obs.loc[mask, 'celltype_final'].values

        zscore_matrix, celltype_labels = compute_nhood_enrichment(coords, celltypes, k=k)

        im = ax.imshow(zscore_matrix, cmap='RdBu_r', vmin=-10, vmax=10, aspect='auto')
        n_types = len(celltype_labels)

        if n_types <= 12:
            ax.set_xticks(range(n_types))
            ax.set_yticks(range(n_types))
            ax.set_xticklabels(celltype_labels, rotation=90, fontsize=5)
            ax.set_yticklabels(celltype_labels, fontsize=5)
        else:
            ax.set_xticks([]); ax.set_yticks([])

        ax.set_title(f'{n_types} cell types', fontsize=9)
        if col == 2:
            cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
            cbar.set_label('z-score', fontsize=9)
        if col == 0:
            ax.text(-0.15, 0.5, 'E. Neighborhood\nEnrichment', transform=ax.transAxes,
                   fontsize=11, fontweight='bold', va='center', ha='right')

    plt.tight_layout(rect=[0.05, 0.02, 1, 0.96])
    fig.suptitle('Spatial Autocorrelation of Cancer Composite Score\nGastric Cancer Progression: Normal → Metaplasia → Cancer',
                fontsize=14, fontweight='bold', y=0.99)

    # Save figure
    fig.savefig(f'{OUTPUT_DIR}/Figure_Spatial_Autocorrelation.png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(f'{OUTPUT_DIR}/Figure_Spatial_Autocorrelation.pdf', dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\nSaved figure to {OUTPUT_DIR}/")
    plt.close()

    return pd.DataFrame(all_stats)

def create_results_table(stats_df):
    """Create publication-quality results table."""

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.axis('off')

    # Format data for table
    table_data = []
    for _, row in stats_df.iterrows():
        sig = '***' if row['Morans_p'] < 0.001 else '**' if row['Morans_p'] < 0.01 else '*' if row['Morans_p'] < 0.05 else ''
        table_data.append([
            row['Sample'],
            row['Stage'],
            f"{row['n_cells']:,}",
            f"{row['Morans_I']:.3f}{sig}",
            f"{row['Morans_z']:.1f}",
            f"{row['LISA_HH']:,}",
            f"{row['LISA_LL']:,}",
            f"{row['Gi_Hot']:,}",
            f"{row['Gi_Cold']:,}"
        ])

    columns = ['Sample', 'Stage', 'Cells', "Moran's I", 'z-score',
               'LISA HH', 'LISA LL', 'Hot Spots', 'Cold Spots']

    table = ax.table(cellText=table_data, colLabels=columns, loc='center',
                    cellLoc='center', colColours=['#4472C4']*len(columns))

    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.8)

    # Style header
    for i in range(len(columns)):
        table[(0, i)].set_text_props(color='white', fontweight='bold')

    # Highlight cancer row
    for i in range(len(columns)):
        table[(3, i)].set_facecolor('#FFE6E6')

    plt.title('Table 1. Spatial Autocorrelation Statistics\n', fontsize=12, fontweight='bold')

    fig.savefig(f'{OUTPUT_DIR}/Table_Spatial_Statistics.png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(f'{OUTPUT_DIR}/Table_Spatial_Statistics.pdf', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    # Save CSV
    stats_df.to_csv(f'{OUTPUT_DIR}/spatial_autocorrelation_results.csv', index=False)
    print(f"Saved table to {OUTPUT_DIR}/")

def main():
    print("Loading data...")
    adata = sc.read_h5ad('output/overnight_annotated_v2.h5ad')
    print(f"Loaded {adata.n_obs:,} cells")

    # Generate figure and stats
    stats_df = create_publication_figure(adata, k=15)

    # Create table
    create_results_table(stats_df)

    print("\n" + "="*60)
    print("PUBLICATION OUTPUTS COMPLETE")
    print("="*60)
    print(f"\nFiles saved to: {OUTPUT_DIR}")
    print("\nKey Finding:")
    cancer_I = stats_df[stats_df['Stage'] == 'Cancer']['Morans_I'].values[0]
    normal_I = stats_df[stats_df['Stage'] == 'Normal']['Morans_I'].values[0]
    print(f"  Cancer Moran's I ({cancer_I:.3f}) is {cancer_I/normal_I:.1f}× higher than Normal ({normal_I:.3f})")

if __name__ == '__main__':
    main()
