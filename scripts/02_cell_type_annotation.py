#!/usr/bin/env python3
"""
G4X Choi Batch 2 - Cell Type Annotation
=======================================

Cell type annotation using protein marker expression.
Hierarchical gating strategy based on established immune panels.
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
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, MiniBatchKMeans
from scipy.stats import zscore
import warnings
warnings.filterwarnings('ignore')

# Configuration
DATA_ROOT = Path('/mnt/x/Choi_Batch_2_Tuesday')
OUTPUT_DIR = Path('/home/user/g4x-choi-batch2-analysis/results')
FIG_DIR = OUTPUT_DIR / 'figures'
TABLE_DIR = OUTPUT_DIR / 'tables'

# Protein markers available
PROTEIN_COLS = [
    'PD1_intensity_mean', 'CD8_intensity_mean', 'PDL1_intensity_mean',
    'FOXP3_intensity_mean', 'KI67_intensity_mean', 'CD45_intensity_mean',
    'PanCK_intensity_mean', 'Isotype_intensity_mean', 'CD20_intensity_mean',
    'CD4_intensity_mean', 'CD11c_intensity_mean', 'HLA-DR_intensity_mean',
    'CD68_intensity_mean', 'CD31_intensity_mean', 'ATPase_intensity_mean',
    'CD3_intensity_mean', 'aSMA_intensity_mean'
]

# Clean marker names
MARKER_MAP = {col: col.replace('_intensity_mean', '') for col in PROTEIN_COLS}

def load_sample_metadata(lane_idx, sample_prefix):
    """Load cell metadata for a sample."""
    lane = f'L00{lane_idx + 1}'
    sample = f'{sample_prefix}{lane_idx + 1}'
    lane_dir = list(DATA_ROOT.glob(f'g4-028-083-FC1-{lane}_*'))[0]
    sample_path = lane_dir / sample / 'single_cell_data' / 'cell_metadata.csv.gz'

    if not sample_path.exists():
        return None

    with gzip.open(sample_path, 'rt') as f:
        df = pd.read_csv(f)
    return df

def gate_cell_types(df, quantile_thresh=0.8):
    """
    Gate cell types based on protein marker expression.

    Hierarchical gating strategy:
    1. CD45+ → Immune cells
       - CD3+ → T cells
         - CD4+/CD8- → CD4+ T cells
           - FOXP3+ → Tregs
         - CD8+/CD4- → CD8+ T cells
       - CD20+ → B cells
       - CD68+ → Macrophages
       - CD11c+/HLA-DR+ → Dendritic cells
    2. PanCK+ → Epithelial cells
    3. aSMA+ → Stromal/Fibroblasts
    4. CD31+ → Endothelial cells
    """

    # Calculate dynamic thresholds based on quantiles
    thresholds = {}
    for col in PROTEIN_COLS:
        if col in df.columns:
            marker = MARKER_MAP[col]
            # Use 80th percentile as positive threshold
            thresholds[marker] = df[col].quantile(quantile_thresh)

    # Initialize cell type column
    cell_types = np.full(len(df), 'Other', dtype=object)

    # Get marker intensities
    def get_marker(marker):
        col = f'{marker}_intensity_mean'
        return df[col].values if col in df.columns else np.zeros(len(df))

    # Get thresholds
    def get_thresh(marker, default=0):
        return thresholds.get(marker, default)

    # Step 1: CD45+ Immune cells
    cd45 = get_marker('CD45')
    cd45_pos = cd45 > get_thresh('CD45')

    # Step 2: T cells (CD3+)
    cd3 = get_marker('CD3')
    cd3_pos = cd3 > get_thresh('CD3')

    # Step 3: CD4 and CD8
    cd4 = get_marker('CD4')
    cd8 = get_marker('CD8')
    cd4_pos = cd4 > get_thresh('CD4')
    cd8_pos = cd8 > get_thresh('CD8')

    # Step 4: FOXP3 for Tregs
    foxp3 = get_marker('FOXP3')
    foxp3_pos = foxp3 > get_thresh('FOXP3')

    # Step 5: B cells (CD20+)
    cd20 = get_marker('CD20')
    cd20_pos = cd20 > get_thresh('CD20')

    # Step 6: Macrophages (CD68+)
    cd68 = get_marker('CD68')
    cd68_pos = cd68 > get_thresh('CD68')

    # Step 7: DCs (CD11c+, HLA-DR+)
    cd11c = get_marker('CD11c')
    hladr = get_marker('HLA-DR')
    dc_pos = (cd11c > get_thresh('CD11c')) & (hladr > get_thresh('HLA-DR'))

    # Step 8: Epithelial (PanCK+)
    panck = get_marker('PanCK')
    panck_pos = panck > get_thresh('PanCK')

    # Step 9: Stromal (aSMA+)
    asma = get_marker('aSMA')
    asma_pos = asma > get_thresh('aSMA')

    # Step 10: Endothelial (CD31+)
    cd31 = get_marker('CD31')
    cd31_pos = cd31 > get_thresh('CD31')

    # Apply gating hierarchy
    # Epithelial first (usually most distinctive)
    cell_types[panck_pos & ~cd45_pos] = 'Epithelial'

    # Stromal
    cell_types[(asma_pos & ~panck_pos & ~cd45_pos)] = 'Stromal'

    # Endothelial
    cell_types[(cd31_pos & ~panck_pos & ~cd45_pos)] = 'Endothelial'

    # Immune cells (CD45+)
    # T cells
    cell_types[(cd45_pos & cd3_pos & cd8_pos & ~cd4_pos)] = 'CD8_T'
    cell_types[(cd45_pos & cd3_pos & cd4_pos & ~cd8_pos & foxp3_pos)] = 'Treg'
    cell_types[(cd45_pos & cd3_pos & cd4_pos & ~cd8_pos & ~foxp3_pos)] = 'CD4_T'
    cell_types[(cd45_pos & cd3_pos & ~cd4_pos & ~cd8_pos)] = 'T_cell'

    # B cells
    cell_types[(cd45_pos & cd20_pos & ~cd3_pos)] = 'B_cell'

    # Macrophages
    cell_types[(cd45_pos & cd68_pos & ~cd3_pos & ~cd20_pos)] = 'Macrophage'

    # Dendritic cells
    cell_types[(cd45_pos & dc_pos & ~cd3_pos & ~cd20_pos & ~cd68_pos)] = 'DC'

    # Generic immune
    cell_types[(cd45_pos & (cell_types == 'Other'))] = 'Immune_other'

    return cell_types

def analyze_sample(lane_idx, sample_prefix):
    """Analyze single sample cell type composition."""
    df = load_sample_metadata(lane_idx, sample_prefix)
    if df is None:
        return None

    sample_id = f'{sample_prefix}{lane_idx + 1}'

    # Gate cell types
    cell_types = gate_cell_types(df)
    df['cell_type'] = cell_types

    # Get counts
    type_counts = df['cell_type'].value_counts()

    # Summary
    summary = {
        'sample_id': sample_id,
        'total_cells': len(df),
    }

    # Add cell type percentages
    for ct in ['Epithelial', 'CD8_T', 'CD4_T', 'Treg', 'T_cell',
               'B_cell', 'Macrophage', 'DC', 'Stromal', 'Endothelial',
               'Immune_other', 'Other']:
        summary[f'{ct}_pct'] = type_counts.get(ct, 0) / len(df) * 100

    return summary, df

def run_clustering(protein_data, n_clusters=12):
    """Run unsupervised clustering on protein data."""

    # Standardize
    scaler = StandardScaler()
    scaled = scaler.fit_transform(protein_data)

    # MiniBatch KMeans for large data
    kmeans = MiniBatchKMeans(n_clusters=n_clusters, random_state=42,
                             batch_size=10000, n_init=3)
    clusters = kmeans.fit_predict(scaled)

    return clusters, kmeans, scaler

def plot_cell_type_composition(results_df, output_dir):
    """Plot cell type composition across samples."""

    # Cell type columns
    ct_cols = [c for c in results_df.columns if c.endswith('_pct')]
    ct_names = [c.replace('_pct', '') for c in ct_cols]

    # Sort by total immune content
    immune_cols = ['CD8_T_pct', 'CD4_T_pct', 'Treg_pct', 'T_cell_pct',
                   'B_cell_pct', 'Macrophage_pct', 'DC_pct', 'Immune_other_pct']
    results_df['total_immune'] = results_df[[c for c in immune_cols if c in results_df.columns]].sum(axis=1)
    results_df = results_df.sort_values('total_immune', ascending=False)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Stacked bar chart
    ax = axes[0, 0]
    ct_data = results_df[ct_cols].values
    bottom = np.zeros(len(results_df))

    colors = plt.cm.tab20(np.linspace(0, 1, len(ct_cols)))
    for i, (ct, color) in enumerate(zip(ct_cols, colors)):
        ax.bar(range(len(results_df)), results_df[ct].values, bottom=bottom,
               label=ct.replace('_pct', ''), color=color)
        bottom += results_df[ct].values

    ax.set_xlabel('Sample (sorted by immune content)')
    ax.set_ylabel('Cell Type (%)')
    ax.set_title('Cell Type Composition')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    ax.set_xticks(range(len(results_df)))
    ax.set_xticklabels(results_df['sample_id'], rotation=45, ha='right', fontsize=7)

    # 2. Heatmap of cell type proportions
    ax = axes[0, 1]
    ct_matrix = results_df[ct_cols].values
    im = ax.imshow(ct_matrix, aspect='auto', cmap='viridis')
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Sample')
    ax.set_xticks(range(len(ct_cols)))
    ax.set_xticklabels([c.replace('_pct', '') for c in ct_cols], rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(len(results_df)))
    ax.set_yticklabels(results_df['sample_id'].values, fontsize=7)
    ax.set_title('Cell Type Proportions Heatmap')
    plt.colorbar(im, ax=ax, label='%')

    # 3. Box plot of main cell types
    ax = axes[1, 0]
    main_cts = ['Epithelial', 'CD8_T', 'CD4_T', 'Macrophage', 'B_cell', 'Stromal']
    box_data = [results_df[f'{ct}_pct'].values for ct in main_cts if f'{ct}_pct' in results_df.columns]
    ax.boxplot(box_data, labels=[ct for ct in main_cts if f'{ct}_pct' in results_df.columns])
    ax.set_ylabel('Cell Type (%)')
    ax.set_title('Distribution of Major Cell Types')
    ax.tick_params(axis='x', rotation=45)

    # 4. Immune vs non-immune
    ax = axes[1, 1]
    non_immune = results_df['Epithelial_pct'] + results_df['Stromal_pct'] + results_df['Endothelial_pct'] + results_df['Other_pct']
    ax.scatter(non_immune, results_df['total_immune'], alpha=0.6)
    ax.set_xlabel('Non-Immune Cells (%)')
    ax.set_ylabel('Immune Cells (%)')
    ax.set_title('Immune vs Non-Immune Content')

    # Add sample labels for outliers
    for i, row in results_df.iterrows():
        ni = row['Epithelial_pct'] + row['Stromal_pct'] + row['Endothelial_pct'] + row['Other_pct']
        if row['total_immune'] > 40 or ni > 80:
            ax.annotate(row['sample_id'], (ni, row['total_immune']), fontsize=8)

    plt.tight_layout()
    plt.savefig(output_dir / 'cell_type_composition.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'cell_type_composition.pdf', bbox_inches='tight')
    plt.close()

    print(f"Saved cell type composition to {output_dir / 'cell_type_composition.png'}")

def plot_marker_expression(marker_stats, output_dir):
    """Plot marker expression across cell types."""

    fig, ax = plt.subplots(figsize=(14, 10))

    # Create heatmap of mean expression by cell type
    cell_types = marker_stats.index.tolist()
    markers = marker_stats.columns.tolist()

    # Z-score normalize
    zscore_data = marker_stats.apply(zscore, axis=0)

    im = ax.imshow(zscore_data.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_xticks(range(len(markers)))
    ax.set_xticklabels(markers, rotation=45, ha='right')
    ax.set_yticks(range(len(cell_types)))
    ax.set_yticklabels(cell_types)
    ax.set_xlabel('Protein Marker')
    ax.set_ylabel('Cell Type')
    ax.set_title('Marker Expression by Cell Type (Z-score)')
    plt.colorbar(im, ax=ax, label='Z-score')

    plt.tight_layout()
    plt.savefig(output_dir / 'marker_expression_by_celltype.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'marker_expression_by_celltype.pdf', bbox_inches='tight')
    plt.close()

    print(f"Saved marker expression heatmap to {output_dir / 'marker_expression_by_celltype.png'}")

def main():
    """Run cell type annotation."""
    print("=" * 60)
    print("G4X Choi Batch 2 - Cell Type Annotation")
    print("=" * 60)

    # Process all samples
    results = []
    all_marker_data = []

    print("\nProcessing samples...")
    for lane_idx in range(4):
        for sample_prefix in ['A0', 'B0', 'C0', 'D0', 'E0', 'F0', 'G0', 'H0']:
            sample_id = f'{sample_prefix}{lane_idx + 1}'
            print(f"  Processing {sample_id}...", end=' ')

            try:
                result = analyze_sample(lane_idx, sample_prefix)
                if result is not None:
                    summary, df = result
                    results.append(summary)

                    # Collect marker data for aggregate analysis
                    if len(all_marker_data) < 200000:  # Limit for memory
                        sample_size = min(5000, len(df))
                        sampled = df.sample(sample_size, random_state=42)
                        sampled['sample_id'] = sample_id
                        all_marker_data.append(sampled)

                    print(f"done ({summary['total_cells']:,} cells)")
                else:
                    print("skipped (no data)")
            except Exception as e:
                print(f"error: {e}")

    # Create results DataFrame
    results_df = pd.DataFrame(results)
    print(f"\nProcessed {len(results_df)} samples")

    # Save results
    results_df.to_csv(TABLE_DIR / 'cell_type_proportions.csv', index=False)

    # Combine marker data
    print("\nAggregating marker expression data...")
    combined_df = pd.concat(all_marker_data, ignore_index=True)
    print(f"  Collected {len(combined_df):,} cells for analysis")

    # Calculate marker expression by cell type
    marker_cols = [c for c in combined_df.columns if '_intensity_mean' in c and c in MARKER_MAP]
    marker_names = [MARKER_MAP[c] for c in marker_cols]

    marker_stats = combined_df.groupby('cell_type')[marker_cols].mean()
    marker_stats.columns = marker_names
    marker_stats.to_csv(TABLE_DIR / 'marker_expression_by_celltype.csv')

    # Generate plots
    print("\n" + "-" * 40)
    print("Generating plots...")

    plot_cell_type_composition(results_df, FIG_DIR)
    plot_marker_expression(marker_stats, FIG_DIR)

    # Summary statistics
    print("\n" + "=" * 60)
    print("CELL TYPE ANNOTATION SUMMARY")
    print("=" * 60)

    ct_cols = [c for c in results_df.columns if c.endswith('_pct')]
    for ct_col in ct_cols:
        ct = ct_col.replace('_pct', '')
        mean_pct = results_df[ct_col].mean()
        if mean_pct > 0.5:  # Only show cell types with >0.5% mean
            print(f"  {ct}: {mean_pct:.1f}% (range: {results_df[ct_col].min():.1f}-{results_df[ct_col].max():.1f}%)")

    print("\n" + "-" * 40)
    print("TOP 5 SAMPLES BY IMMUNE INFILTRATION:")
    top_immune = results_df.nlargest(5, 'total_immune')
    for _, row in top_immune.iterrows():
        print(f"  {row['sample_id']}: {row['total_immune']:.1f}% immune cells")

    print("\n" + "-" * 40)
    print("TOP 5 SAMPLES BY EPITHELIAL CONTENT:")
    top_epi = results_df.nlargest(5, 'Epithelial_pct')
    for _, row in top_epi.iterrows():
        print(f"  {row['sample_id']}: {row['Epithelial_pct']:.1f}% epithelial cells")

    print("\n" + "=" * 60)
    print("Cell Type Annotation Complete!")
    print(f"Results saved to: {TABLE_DIR}")
    print("=" * 60)

if __name__ == '__main__':
    main()
