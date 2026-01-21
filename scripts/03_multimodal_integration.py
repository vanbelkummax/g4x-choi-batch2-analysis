#!/usr/bin/env python3
"""
G4X Choi Batch 2 - Multimodal Integration Analysis
==================================================

Analyze RNA-protein correlation and cross-modality relationships.
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
from scipy.stats import spearmanr, pearsonr
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# Configuration
DATA_ROOT = Path('/mnt/x/Choi_Batch_2_Tuesday')
OUTPUT_DIR = Path('/home/user/g4x-choi-batch2-analysis/results')
FIG_DIR = OUTPUT_DIR / 'figures'
TABLE_DIR = OUTPUT_DIR / 'tables'

# Protein to gene mapping
PROTEIN_TO_GENE = {
    'CD3': ['CD3D', 'CD3E', 'CD3G'],
    'CD4': ['CD4'],
    'CD8': ['CD8A', 'CD8B'],
    'CD20': ['MS4A1'],
    'CD45': ['PTPRC'],
    'CD68': ['CD68'],
    'CD31': ['PECAM1'],
    'CD11c': ['ITGAX'],
    'HLA-DR': ['HLA-DRA', 'HLA-DRB1'],
    'PanCK': ['KRT8', 'KRT18', 'KRT19', 'KRT7'],  # Pan-cytokeratin
    'aSMA': ['ACTA2'],
    'KI67': ['MKI67'],
    'PD1': ['PDCD1'],
    'PDL1': ['CD274'],
    'FOXP3': ['FOXP3'],
    'ATPase': ['ATP1A1', 'ATP1B1'],  # Na/K ATPase
}

def load_sample_data(lane_idx, sample_prefix):
    """Load protein and RNA data for a sample."""
    lane = f'L00{lane_idx + 1}'
    sample = f'{sample_prefix}{lane_idx + 1}'
    lane_dir = list(DATA_ROOT.glob(f'g4-028-083-FC1-{lane}_*'))[0]
    sample_path = lane_dir / sample

    # Load cell metadata (has protein intensities)
    metadata_file = sample_path / 'single_cell_data' / 'cell_metadata.csv.gz'
    if not metadata_file.exists():
        return None, None, None

    with gzip.open(metadata_file, 'rt') as f:
        protein_df = pd.read_csv(f)

    # Load H5 file for RNA
    h5_file = sample_path / 'single_cell_data' / 'feature_matrix.h5'
    if not h5_file.exists():
        return None, None, None

    with h5py.File(h5_file, 'r') as f:
        X = f['X'][:]
        cell_ids = [c.decode() if isinstance(c, bytes) else c for c in f['obs/cell_id'][:]]
        gene_ids = [g.decode() if isinstance(g, bytes) else g for g in f['var/gene_id'][:]]

    # Create RNA dataframe
    rna_df = pd.DataFrame(X, columns=gene_ids)
    rna_df['cell_id'] = cell_ids

    return protein_df, rna_df, gene_ids

def calculate_protein_rna_correlation(protein_df, rna_df, gene_list):
    """Calculate correlation between protein and RNA for matching genes."""

    correlations = []

    for protein, genes in PROTEIN_TO_GENE.items():
        protein_col = f'{protein}_intensity_mean'
        if protein_col not in protein_df.columns:
            continue

        protein_values = protein_df[protein_col].values

        for gene in genes:
            if gene not in rna_df.columns:
                continue

            rna_values = rna_df[gene].values

            # Remove zeros for more meaningful correlation
            mask = (protein_values > 0) | (rna_values > 0)
            if mask.sum() < 100:
                continue

            # Spearman correlation (rank-based, more robust)
            spearman_r, spearman_p = spearmanr(protein_values[mask], rna_values[mask])

            # Pearson correlation
            pearson_r, pearson_p = pearsonr(protein_values[mask], rna_values[mask])

            correlations.append({
                'protein': protein,
                'gene': gene,
                'spearman_r': spearman_r,
                'spearman_p': spearman_p,
                'pearson_r': pearson_r,
                'pearson_p': pearson_p,
                'n_cells': mask.sum()
            })

    return pd.DataFrame(correlations)

def analyze_sample_correlations(lane_idx, sample_prefix):
    """Analyze correlations for a single sample."""
    protein_df, rna_df, gene_list = load_sample_data(lane_idx, sample_prefix)
    if protein_df is None:
        return None

    corr_df = calculate_protein_rna_correlation(protein_df, rna_df, gene_list)
    if len(corr_df) == 0:
        return None

    sample_id = f'{sample_prefix}{lane_idx + 1}'
    corr_df['sample_id'] = sample_id

    return corr_df

def plot_correlation_summary(all_corr_df, output_dir):
    """Plot summary of protein-RNA correlations."""

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Heatmap of correlations by protein-gene pair
    ax = axes[0, 0]
    pivot_df = all_corr_df.pivot_table(
        index='sample_id',
        columns=['protein', 'gene'],
        values='spearman_r',
        aggfunc='mean'
    )

    # Take mean across samples
    mean_corr = all_corr_df.groupby(['protein', 'gene'])['spearman_r'].mean()
    mean_corr_df = mean_corr.reset_index()
    mean_corr_df.columns = ['protein', 'gene', 'mean_spearman_r']

    # Create pivot for heatmap
    unique_pairs = mean_corr_df[['protein', 'gene']].drop_duplicates()
    labels = [f"{row['protein']}-{row['gene']}" for _, row in unique_pairs.iterrows()]

    ax.barh(range(len(mean_corr_df)), mean_corr_df['mean_spearman_r'].values)
    ax.set_yticks(range(len(mean_corr_df)))
    ax.set_yticklabels([f"{row['protein']}-{row['gene']}" for _, row in mean_corr_df.iterrows()], fontsize=9)
    ax.set_xlabel('Mean Spearman Correlation')
    ax.set_title('Protein-RNA Correlation (Mean Across Samples)')
    ax.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    ax.axvline(x=0.3, color='g', linestyle='--', alpha=0.5, label='r=0.3')
    ax.axvline(x=-0.3, color='g', linestyle='--', alpha=0.5)
    ax.legend()

    # 2. Distribution of correlations
    ax = axes[0, 1]
    ax.hist(all_corr_df['spearman_r'], bins=30, edgecolor='black', alpha=0.7)
    ax.axvline(x=all_corr_df['spearman_r'].mean(), color='r', linestyle='--',
               label=f'Mean: {all_corr_df["spearman_r"].mean():.3f}')
    ax.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    ax.set_xlabel('Spearman Correlation')
    ax.set_ylabel('Frequency')
    ax.set_title('Distribution of Protein-RNA Correlations')
    ax.legend()

    # 3. Sample-level variation
    ax = axes[1, 0]
    sample_means = all_corr_df.groupby('sample_id')['spearman_r'].mean().sort_values()
    ax.barh(range(len(sample_means)), sample_means.values)
    ax.set_yticks(range(len(sample_means)))
    ax.set_yticklabels(sample_means.index, fontsize=7)
    ax.set_xlabel('Mean Spearman Correlation')
    ax.set_title('Mean Protein-RNA Correlation by Sample')
    ax.axvline(x=sample_means.mean(), color='r', linestyle='--')

    # 4. Correlation vs significance
    ax = axes[1, 1]
    ax.scatter(all_corr_df['spearman_r'], -np.log10(all_corr_df['spearman_p']+1e-100),
               alpha=0.3, s=10)
    ax.axhline(y=-np.log10(0.05), color='r', linestyle='--', label='p=0.05')
    ax.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    ax.set_xlabel('Spearman Correlation')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title('Correlation vs Significance')
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_dir / 'protein_rna_correlation.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'protein_rna_correlation.pdf', bbox_inches='tight')
    plt.close()

    print(f"Saved correlation plots to {output_dir / 'protein_rna_correlation.png'}")

def plot_detailed_correlations(all_corr_df, output_dir):
    """Plot detailed per-protein correlations."""

    # Get unique proteins with data
    proteins = all_corr_df['protein'].unique()
    n_proteins = len(proteins)

    fig, axes = plt.subplots(4, 4, figsize=(16, 16))
    axes = axes.flatten()

    for i, protein in enumerate(proteins[:16]):
        ax = axes[i]
        protein_data = all_corr_df[all_corr_df['protein'] == protein]

        if len(protein_data) == 0:
            ax.set_visible(False)
            continue

        # Box plot of correlations across samples
        genes = protein_data['gene'].unique()
        gene_data = [protein_data[protein_data['gene'] == g]['spearman_r'].values for g in genes]

        bp = ax.boxplot(gene_data, labels=genes)
        ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        ax.set_xlabel('Gene')
        ax.set_ylabel('Spearman r')
        ax.set_title(f'{protein} Protein')
        ax.tick_params(axis='x', rotation=45)

    # Hide unused axes
    for j in range(len(proteins), 16):
        axes[j].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_dir / 'protein_rna_by_marker.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'protein_rna_by_marker.pdf', bbox_inches='tight')
    plt.close()

    print(f"Saved per-marker plots to {output_dir / 'protein_rna_by_marker.png'}")

def main():
    """Run multimodal integration analysis."""
    print("=" * 60)
    print("G4X Choi Batch 2 - Multimodal Integration Analysis")
    print("=" * 60)

    # Process all samples
    all_correlations = []

    print("\nCalculating protein-RNA correlations...")
    for lane_idx in range(4):
        for sample_prefix in ['A0', 'B0', 'C0', 'D0', 'E0', 'F0', 'G0', 'H0']:
            sample_id = f'{sample_prefix}{lane_idx + 1}'
            print(f"  Processing {sample_id}...", end=' ')

            try:
                corr_df = analyze_sample_correlations(lane_idx, sample_prefix)
                if corr_df is not None and len(corr_df) > 0:
                    all_correlations.append(corr_df)
                    print(f"done ({len(corr_df)} pairs)")
                else:
                    print("no matching genes")
            except Exception as e:
                print(f"error: {e}")

    # Combine all correlations
    all_corr_df = pd.concat(all_correlations, ignore_index=True)
    print(f"\nTotal correlation measurements: {len(all_corr_df)}")

    # Save results
    all_corr_df.to_csv(TABLE_DIR / 'protein_rna_correlations.csv', index=False)

    # Generate summary statistics
    summary = all_corr_df.groupby(['protein', 'gene']).agg({
        'spearman_r': ['mean', 'std', 'min', 'max'],
        'pearson_r': ['mean', 'std'],
        'n_cells': 'mean'
    }).round(3)
    summary.columns = ['_'.join(col).strip() for col in summary.columns.values]
    summary = summary.reset_index()
    summary.to_csv(TABLE_DIR / 'protein_rna_correlation_summary.csv', index=False)

    # Generate plots
    print("\n" + "-" * 40)
    print("Generating plots...")

    plot_correlation_summary(all_corr_df, FIG_DIR)
    plot_detailed_correlations(all_corr_df, FIG_DIR)

    # Print summary
    print("\n" + "=" * 60)
    print("MULTIMODAL INTEGRATION SUMMARY")
    print("=" * 60)

    print("\nProtein-RNA Correlation Statistics:")
    print(f"  Mean Spearman r: {all_corr_df['spearman_r'].mean():.3f}")
    print(f"  Median Spearman r: {all_corr_df['spearman_r'].median():.3f}")
    print(f"  Std Spearman r: {all_corr_df['spearman_r'].std():.3f}")

    # Best correlations
    print("\n" + "-" * 40)
    print("TOP PROTEIN-RNA CORRELATIONS (Mean across samples):")
    top_corr = summary.nlargest(10, 'spearman_r_mean')
    for _, row in top_corr.iterrows():
        print(f"  {row['protein']}-{row['gene']}: r = {row['spearman_r_mean']:.3f} (±{row['spearman_r_std']:.3f})")

    # Weakest correlations
    print("\n" + "-" * 40)
    print("WEAKEST PROTEIN-RNA CORRELATIONS (potential discordance):")
    weak_corr = summary.nsmallest(5, 'spearman_r_mean')
    for _, row in weak_corr.iterrows():
        print(f"  {row['protein']}-{row['gene']}: r = {row['spearman_r_mean']:.3f} (±{row['spearman_r_std']:.3f})")

    print("\n" + "-" * 40)
    print("INTERPRETATION:")
    mean_r = all_corr_df['spearman_r'].mean()
    if mean_r > 0.3:
        print("  Strong overall protein-RNA correlation suggests good concordance.")
    elif mean_r > 0.1:
        print("  Moderate protein-RNA correlation - typical for imaging-based platforms.")
    else:
        print("  Weak protein-RNA correlation - may indicate biological post-transcriptional regulation")
        print("  or technical factors (different sensitivity, dynamic range).")

    print("\n" + "=" * 60)
    print("Multimodal Integration Analysis Complete!")
    print(f"Results saved to: {TABLE_DIR}")
    print("=" * 60)

if __name__ == '__main__':
    main()
