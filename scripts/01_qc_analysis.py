#!/usr/bin/env python3
"""
G4X Choi Batch 2 - Quality Control Analysis
============================================

Comprehensive QC of multimodal G4X spatial transcriptomics data.
- Cell-level QC metrics
- Protein signal-to-noise analysis
- Transcript quality assessment
- Sample-level comparisons
"""

import os
import sys
import h5py
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Configuration
DATA_ROOT = Path('/mnt/x/Choi_Batch_2_Tuesday')
OUTPUT_DIR = Path('/home/user/g4x-choi-batch2-analysis/results')
FIG_DIR = OUTPUT_DIR / 'figures'
TABLE_DIR = OUTPUT_DIR / 'tables'

# Create output directories
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

# Sample layout
LANES = ['L001', 'L002', 'L003', 'L004']
SAMPLES_PER_LANE = ['A0', 'B0', 'C0', 'D0', 'E0', 'F0', 'G0', 'H0']

def get_sample_path(lane_idx, sample_prefix):
    """Get path to sample directory."""
    lane = f'L00{lane_idx + 1}'
    sample = f'{sample_prefix}{lane_idx + 1}'
    lane_dir = list(DATA_ROOT.glob(f'g4-028-083-FC1-{lane}_*'))[0]
    return lane_dir / sample

def load_core_metrics():
    """Load pre-computed core metrics."""
    metrics = pd.read_csv(DATA_ROOT / 'choi_preGC_b2_core_metrics.csv')
    return metrics

def load_protein_metrics():
    """Load protein core metrics."""
    protein_metrics = pd.read_csv(DATA_ROOT / 'choi_preGC_b2_protein_core_metrics.csv')
    return protein_metrics

def load_sample_metadata(sample_path):
    """Load cell metadata for a sample."""
    metadata_file = sample_path / 'single_cell_data' / 'cell_metadata.csv.gz'
    if not metadata_file.exists():
        return None
    with gzip.open(metadata_file, 'rt') as f:
        df = pd.read_csv(f)
    return df

def load_h5_data(sample_path):
    """Load expression matrix from H5 file."""
    h5_file = sample_path / 'single_cell_data' / 'feature_matrix.h5'
    if not h5_file.exists():
        return None, None, None

    with h5py.File(h5_file, 'r') as f:
        X = f['X'][:]
        cell_ids = [c.decode() if isinstance(c, bytes) else c for c in f['obs/cell_id'][:]]
        gene_ids = [g.decode() if isinstance(g, bytes) else g for g in f['var/gene_id'][:]]

    return X, cell_ids, gene_ids

# Protein markers
PROTEIN_MARKERS = [
    'ATPase', 'CD11c', 'CD20', 'CD3', 'CD31', 'CD4', 'CD45',
    'CD68', 'CD8', 'FOXP3', 'HLA-DR', 'Isotype', 'KI67',
    'PD1', 'PDL1', 'PanCK', 'aSMA'
]

# Cell type markers for gating
CELL_TYPE_MARKERS = {
    'T_cells': ['CD3'],
    'CD4_T': ['CD3', 'CD4'],
    'CD8_T': ['CD3', 'CD8'],
    'Tregs': ['CD3', 'CD4', 'FOXP3'],
    'B_cells': ['CD20'],
    'Macrophages': ['CD68'],
    'DCs': ['CD11c', 'HLA-DR'],
    'Endothelial': ['CD31'],
    'Epithelial': ['PanCK'],
    'Stromal': ['aSMA']
}

def analyze_protein_qc(core_metrics, protein_metrics):
    """Analyze protein signal quality across samples."""

    # Extract SNR columns
    snr_cols = [c for c in protein_metrics.columns if '_snr' in c]

    # Create SNR summary
    snr_data = []
    for col in snr_cols:
        marker = col.replace('_snr', '')
        snr_data.append({
            'Marker': marker,
            'Mean_SNR': protein_metrics[col].mean(),
            'Median_SNR': protein_metrics[col].median(),
            'Min_SNR': protein_metrics[col].min(),
            'Max_SNR': protein_metrics[col].max(),
            'Std_SNR': protein_metrics[col].std()
        })

    snr_df = pd.DataFrame(snr_data)
    snr_df = snr_df.sort_values('Median_SNR', ascending=False)

    return snr_df

def plot_qc_overview(core_metrics, output_dir):
    """Generate QC overview plots."""

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # 1. Cell counts per sample
    ax = axes[0, 0]
    ax.bar(range(len(core_metrics)), core_metrics['number_cells'].values / 1000)
    ax.set_xlabel('Sample Index')
    ax.set_ylabel('Cells (thousands)')
    ax.set_title('Cell Counts per Sample')
    ax.axhline(y=core_metrics['number_cells'].mean()/1000, color='r', linestyle='--',
               label=f'Mean: {core_metrics["number_cells"].mean()/1000:.1f}K')
    ax.legend()

    # 2. Transcripts per cell
    ax = axes[0, 1]
    ax.hist(core_metrics['median_transcripts_per_cell'], bins=15, edgecolor='black')
    ax.set_xlabel('Median Transcripts per Cell')
    ax.set_ylabel('Number of Samples')
    ax.set_title('Transcript Detection Distribution')
    ax.axvline(x=50, color='r', linestyle='--', label='QC Threshold (50)')
    ax.legend()

    # 3. Genes per cell
    ax = axes[0, 2]
    ax.hist(core_metrics['median_unique_genes_per_cell'], bins=15, edgecolor='black')
    ax.set_xlabel('Median Genes per Cell')
    ax.set_ylabel('Number of Samples')
    ax.set_title('Gene Detection Distribution')
    ax.axvline(x=25, color='r', linestyle='--', label='QC Threshold (25)')
    ax.legend()

    # 4. Empty cell percentage
    ax = axes[1, 0]
    ax.bar(range(len(core_metrics)), core_metrics['pct_empty_cells'].values)
    ax.set_xlabel('Sample Index')
    ax.set_ylabel('% Empty Cells')
    ax.set_title('Empty Cell Percentage')
    ax.axhline(y=5, color='r', linestyle='--', label='5% threshold')
    ax.legend()

    # 5. Transcript capture efficiency
    ax = axes[1, 1]
    ax.scatter(core_metrics['total_tissue_area_mm^2'],
               core_metrics['total_transcripts']/1e6,
               c=core_metrics['number_cells']/1000, cmap='viridis')
    ax.set_xlabel('Tissue Area (mm²)')
    ax.set_ylabel('Total Transcripts (millions)')
    ax.set_title('Transcripts vs Tissue Area')
    plt.colorbar(ax.collections[0], ax=ax, label='Cells (K)')

    # 6. Quality score (Q20)
    ax = axes[1, 2]
    ax.hist(core_metrics['pct_q20'], bins=15, edgecolor='black')
    ax.set_xlabel('% Q20')
    ax.set_ylabel('Number of Samples')
    ax.set_title('Sequencing Quality (Q20)')

    plt.tight_layout()
    plt.savefig(output_dir / 'qc_overview.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'qc_overview.pdf', bbox_inches='tight')
    plt.close()

    print(f"Saved QC overview to {output_dir / 'qc_overview.png'}")

def plot_protein_snr(snr_df, protein_metrics, output_dir):
    """Plot protein signal-to-noise ratios."""

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 1. SNR by marker (summary)
    ax = axes[0]
    colors = ['green' if x > 2 else 'orange' if x > 1 else 'red'
              for x in snr_df['Median_SNR']]
    ax.barh(range(len(snr_df)), snr_df['Median_SNR'].values, color=colors)
    ax.set_yticks(range(len(snr_df)))
    ax.set_yticklabels(snr_df['Marker'].values)
    ax.set_xlabel('Median SNR')
    ax.set_title('Protein Signal-to-Noise Ratio')
    ax.axvline(x=2, color='green', linestyle='--', alpha=0.5, label='Good (>2)')
    ax.axvline(x=1, color='orange', linestyle='--', alpha=0.5, label='Marginal (>1)')
    ax.legend(loc='lower right')

    # 2. SNR heatmap across samples
    ax = axes[1]
    snr_cols = [c for c in protein_metrics.columns if '_snr' in c]
    snr_matrix = protein_metrics[snr_cols].values
    markers = [c.replace('_snr', '') for c in snr_cols]

    im = ax.imshow(snr_matrix.T, aspect='auto', cmap='RdYlGn', vmin=0, vmax=5)
    ax.set_xlabel('Sample Index')
    ax.set_ylabel('Protein Marker')
    ax.set_yticks(range(len(markers)))
    ax.set_yticklabels(markers)
    ax.set_title('SNR Across Samples')
    plt.colorbar(im, ax=ax, label='SNR')

    plt.tight_layout()
    plt.savefig(output_dir / 'protein_snr.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'protein_snr.pdf', bbox_inches='tight')
    plt.close()

    print(f"Saved protein SNR plots to {output_dir / 'protein_snr.png'}")

def plot_sample_variability(core_metrics, output_dir):
    """Plot sample-to-sample variability metrics."""

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Lane comparison - cell counts
    ax = axes[0, 0]
    lane_data = []
    for i, lane in enumerate(['L001', 'L002', 'L003', 'L004']):
        lane_metrics = core_metrics.iloc[i*8:(i+1)*8]
        lane_data.extend([(lane, c/1000) for c in lane_metrics['number_cells']])

    lane_df = pd.DataFrame(lane_data, columns=['Lane', 'Cells_K'])
    sns.boxplot(data=lane_df, x='Lane', y='Cells_K', ax=ax)
    ax.set_ylabel('Cells (thousands)')
    ax.set_title('Cell Counts by Lane')

    # 2. Lane comparison - transcripts
    ax = axes[0, 1]
    lane_data = []
    for i, lane in enumerate(['L001', 'L002', 'L003', 'L004']):
        lane_metrics = core_metrics.iloc[i*8:(i+1)*8]
        lane_data.extend([(lane, t) for t in lane_metrics['median_transcripts_per_cell']])

    lane_df = pd.DataFrame(lane_data, columns=['Lane', 'Transcripts'])
    sns.boxplot(data=lane_df, x='Lane', y='Transcripts', ax=ax)
    ax.set_ylabel('Median Transcripts/Cell')
    ax.set_title('Transcript Detection by Lane')

    # 3. Cell area distribution
    ax = axes[1, 0]
    ax.scatter(core_metrics['median_cell_area_um^2'],
               core_metrics['median_transcripts_per_cell'],
               s=core_metrics['number_cells']/2000, alpha=0.6)
    ax.set_xlabel('Median Cell Area (µm²)')
    ax.set_ylabel('Median Transcripts/Cell')
    ax.set_title('Cell Size vs Transcript Detection')

    # 4. Correlation: transcripts vs genes
    ax = axes[1, 1]
    ax.scatter(core_metrics['median_transcripts_per_cell'],
               core_metrics['median_unique_genes_per_cell'])

    # Add correlation
    corr = np.corrcoef(core_metrics['median_transcripts_per_cell'],
                       core_metrics['median_unique_genes_per_cell'])[0,1]
    ax.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax.transAxes,
            fontsize=12, verticalalignment='top')
    ax.set_xlabel('Median Transcripts/Cell')
    ax.set_ylabel('Median Genes/Cell')
    ax.set_title('Transcript vs Gene Detection')

    plt.tight_layout()
    plt.savefig(output_dir / 'sample_variability.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'sample_variability.pdf', bbox_inches='tight')
    plt.close()

    print(f"Saved sample variability plots to {output_dir / 'sample_variability.png'}")

def generate_qc_summary(core_metrics, snr_df, output_dir):
    """Generate summary statistics table."""

    summary = {
        'Metric': [
            'Total Samples',
            'Total Cells',
            'Total Transcripts',
            'Mean Cells/Sample',
            'Median Cells/Sample',
            'Mean Transcripts/Cell',
            'Median Genes/Cell (mean)',
            'Mean Empty Cell %',
            'Mean Q20 %',
            'Samples with >50K cells',
            'Samples with low trans (<40)',
            'High SNR proteins (>2)',
            'Low SNR proteins (<1)'
        ],
        'Value': [
            len(core_metrics),
            f"{core_metrics['number_cells'].sum():,}",
            f"{core_metrics['total_transcripts'].sum():,}",
            f"{core_metrics['number_cells'].mean():,.0f}",
            f"{core_metrics['number_cells'].median():,.0f}",
            f"{core_metrics['median_transcripts_per_cell'].mean():.1f}",
            f"{core_metrics['median_unique_genes_per_cell'].mean():.1f}",
            f"{core_metrics['pct_empty_cells'].mean():.2f}%",
            f"{core_metrics['pct_q20'].mean():.1f}%",
            (core_metrics['number_cells'] > 50000).sum(),
            (core_metrics['median_transcripts_per_cell'] < 40).sum(),
            (snr_df['Median_SNR'] > 2).sum(),
            (snr_df['Median_SNR'] < 1).sum()
        ]
    }

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output_dir / 'qc_summary.csv', index=False)

    # Also save detailed metrics
    core_metrics.to_csv(output_dir / 'sample_metrics_detailed.csv', index=False)
    snr_df.to_csv(output_dir / 'protein_snr_summary.csv', index=False)

    print(f"Saved QC summary to {output_dir / 'qc_summary.csv'}")

    return summary_df

def main():
    """Run QC analysis."""
    print("=" * 60)
    print("G4X Choi Batch 2 - Quality Control Analysis")
    print("=" * 60)

    # Load metrics
    print("\nLoading core metrics...")
    core_metrics = load_core_metrics()
    print(f"  Loaded {len(core_metrics)} samples")

    print("\nLoading protein metrics...")
    protein_metrics = load_protein_metrics()
    print(f"  Loaded protein data for {len(protein_metrics)} samples")

    # Analyze protein SNR
    print("\nAnalyzing protein signal-to-noise ratios...")
    snr_df = analyze_protein_qc(core_metrics, protein_metrics)
    print("\nProtein SNR Summary (sorted by quality):")
    print(snr_df[['Marker', 'Median_SNR', 'Min_SNR', 'Max_SNR']].to_string(index=False))

    # Generate plots
    print("\n" + "-" * 40)
    print("Generating QC plots...")

    plot_qc_overview(core_metrics, FIG_DIR)
    plot_protein_snr(snr_df, protein_metrics, FIG_DIR)
    plot_sample_variability(core_metrics, FIG_DIR)

    # Generate summary
    print("\n" + "-" * 40)
    print("Generating QC summary...")
    summary_df = generate_qc_summary(core_metrics, snr_df, TABLE_DIR)

    print("\n" + "=" * 60)
    print("QC SUMMARY")
    print("=" * 60)
    for _, row in summary_df.iterrows():
        print(f"  {row['Metric']}: {row['Value']}")

    # Flag samples with potential issues
    print("\n" + "-" * 40)
    print("SAMPLES REQUIRING ATTENTION:")

    low_trans = core_metrics[core_metrics['median_transcripts_per_cell'] < 40]
    if len(low_trans) > 0:
        print(f"\n  Low transcript detection (<40/cell):")
        for _, row in low_trans.iterrows():
            print(f"    - {row['sample_id']}: {row['median_transcripts_per_cell']:.0f} trans/cell")

    high_empty = core_metrics[core_metrics['pct_empty_cells'] > 5]
    if len(high_empty) > 0:
        print(f"\n  High empty cell % (>5%):")
        for _, row in high_empty.iterrows():
            print(f"    - {row['sample_id']}: {row['pct_empty_cells']:.1f}%")

    print("\n" + "=" * 60)
    print("QC Analysis Complete!")
    print(f"Results saved to: {OUTPUT_DIR}")
    print("=" * 60)

if __name__ == '__main__':
    main()
