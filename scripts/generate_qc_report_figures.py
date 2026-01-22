#!/usr/bin/env python3
"""
Generate comprehensive QC report figures for G4X Choi Batch 2 Analysis.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import scanpy as sc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Paths
RESULTS_DIR = Path("/home/user/g4x-choi-batch2-analysis/results/qc_all_samples")
FIGURES_DIR = Path("/home/user/g4x-choi-batch2-analysis/reports/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Color schemes
STAGE_COLORS = {
    'control': '#2ecc71',
    'normal': '#3498db',
    'metaplasia': '#f39c12',
    'cancer': '#e74c3c'
}

LANE_COLORS = {
    'L001': '#1abc9c',
    'L002': '#9b59b6',
    'L003': '#e67e22',
    'L004': '#34495e'
}

QC_COLORS = {
    'PASS': '#27ae60',
    'WARN': '#f1c40f',
    'FAIL': '#c0392b'
}


def load_data():
    """Load QC summary and merged data."""
    print("Loading data...")

    # Sample QC summary
    qc_df = pd.read_csv(RESULTS_DIR / "sample_qc_summary.csv")

    # Merged corrected data (for UMAP plots)
    adata = None
    corrected_path = RESULTS_DIR / "merged" / "merged_corrected.h5ad"
    if corrected_path.exists():
        print("Loading merged corrected data...")
        adata = sc.read_h5ad(corrected_path)

    return qc_df, adata


def fig1_sample_overview(qc_df):
    """Figure 1: Sample overview - cells per sample by stage."""
    print("Generating Figure 1: Sample Overview...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1a: Cells per sample by stage
    ax = axes[0, 0]
    stages = ['control', 'normal', 'metaplasia', 'cancer']
    for stage in stages:
        subset = qc_df[qc_df['stage'] == stage]
        ax.bar(subset['sample_id'], subset['n_cells'],
               color=STAGE_COLORS[stage], label=stage, alpha=0.8)
    ax.set_xlabel('Sample ID')
    ax.set_ylabel('Number of Cells')
    ax.set_title('A) Cells per Sample by Disease Stage')
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title='Stage')
    ax.set_ylim(0, qc_df['n_cells'].max() * 1.1)

    # 1b: Cells by lane
    ax = axes[0, 1]
    lane_totals = qc_df.groupby('lane')['n_cells'].sum()
    colors = [LANE_COLORS[lane] for lane in lane_totals.index]
    bars = ax.bar(lane_totals.index, lane_totals.values, color=colors, alpha=0.8)
    ax.set_xlabel('Lane')
    ax.set_ylabel('Total Cells')
    ax.set_title('B) Total Cells per Lane')
    for bar, val in zip(bars, lane_totals.values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5000,
                f'{val:,}', ha='center', va='bottom', fontsize=10)

    # 1c: Median transcripts per cell by sample
    ax = axes[1, 0]
    for stage in stages:
        subset = qc_df[qc_df['stage'] == stage]
        ax.scatter(subset['sample_id'], subset['median_counts'],
                   c=STAGE_COLORS[stage], s=100, label=stage, alpha=0.8)
    ax.axhline(y=30, color='red', linestyle='--', alpha=0.5, label='Min threshold (30)')
    ax.set_xlabel('Sample ID')
    ax.set_ylabel('Median Transcripts per Cell')
    ax.set_title('C) Median Transcripts per Cell by Sample')
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title='Stage', loc='upper right')

    # 1d: QC Status distribution
    ax = axes[1, 1]
    qc_counts = qc_df['qc_status'].value_counts()
    colors = [QC_COLORS[status] for status in qc_counts.index]
    wedges, texts, autotexts = ax.pie(qc_counts.values, labels=qc_counts.index,
                                       colors=colors, autopct='%1.0f%%',
                                       startangle=90, explode=[0.02]*len(qc_counts))
    ax.set_title('D) Sample QC Status Distribution')

    # Add counts to legend
    legend_labels = [f'{status}: {count}' for status, count in qc_counts.items()]
    ax.legend(wedges, legend_labels, loc='lower right')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig1_sample_overview.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig1_sample_overview.png")


def fig2_qc_metrics(qc_df):
    """Figure 2: QC metrics distributions."""
    print("Generating Figure 2: QC Metrics...")

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # 2a: Median genes per cell
    ax = axes[0, 0]
    for stage in ['control', 'normal', 'metaplasia', 'cancer']:
        subset = qc_df[qc_df['stage'] == stage]
        ax.scatter(range(len(subset)), subset['median_genes'].values,
                   c=STAGE_COLORS[stage], s=80, label=stage, alpha=0.7)
    ax.axhline(y=20, color='red', linestyle='--', alpha=0.5, label='Min threshold')
    ax.set_xlabel('Sample Index')
    ax.set_ylabel('Median Genes per Cell')
    ax.set_title('A) Median Genes Detected per Cell')
    ax.legend()

    # 2b: Empty cell percentage
    ax = axes[0, 1]
    bars = ax.bar(qc_df['sample_id'], qc_df['pct_zero_count_cells'],
                  color=[STAGE_COLORS[s] for s in qc_df['stage']], alpha=0.8)
    ax.axhline(y=5, color='red', linestyle='--', alpha=0.5, label='Warning threshold (5%)')
    ax.set_xlabel('Sample ID')
    ax.set_ylabel('% Empty Cells')
    ax.set_title('B) Percentage of Empty Cells per Sample')
    ax.tick_params(axis='x', rotation=45)
    ax.legend()

    # 2c: Mean vs Std counts (quality indicator)
    ax = axes[0, 2]
    for stage in ['control', 'normal', 'metaplasia', 'cancer']:
        subset = qc_df[qc_df['stage'] == stage]
        ax.scatter(subset['mean_counts'], subset['std_counts'],
                   c=STAGE_COLORS[stage], s=100, label=stage, alpha=0.7)
    ax.set_xlabel('Mean Counts per Cell')
    ax.set_ylabel('Std Dev Counts')
    ax.set_title('C) Count Distribution Quality')
    ax.legend()

    # 2d: Cells by stage (boxplot)
    ax = axes[1, 0]
    stage_order = ['control', 'normal', 'metaplasia', 'cancer']
    bp = ax.boxplot([qc_df[qc_df['stage'] == s]['n_cells'].values for s in stage_order],
                    labels=stage_order, patch_artist=True)
    for patch, stage in zip(bp['boxes'], stage_order):
        patch.set_facecolor(STAGE_COLORS[stage])
        patch.set_alpha(0.7)
    ax.set_xlabel('Disease Stage')
    ax.set_ylabel('Number of Cells')
    ax.set_title('D) Cell Count Distribution by Stage')

    # 2e: Median counts by stage
    ax = axes[1, 1]
    bp = ax.boxplot([qc_df[qc_df['stage'] == s]['median_counts'].values for s in stage_order],
                    labels=stage_order, patch_artist=True)
    for patch, stage in zip(bp['boxes'], stage_order):
        patch.set_facecolor(STAGE_COLORS[stage])
        patch.set_alpha(0.7)
    ax.set_xlabel('Disease Stage')
    ax.set_ylabel('Median Transcripts per Cell')
    ax.set_title('E) Transcript Distribution by Stage')

    # 2f: Protein positive percentage
    ax = axes[1, 2]
    prot_data = qc_df.dropna(subset=['pct_protein_positive'])
    ax.bar(prot_data['sample_id'], prot_data['pct_protein_positive'],
           color=[STAGE_COLORS[s] for s in prot_data['stage']], alpha=0.8)
    ax.set_xlabel('Sample ID')
    ax.set_ylabel('% Protein Positive Cells')
    ax.set_title('F) Protein Detection Rate')
    ax.tick_params(axis='x', rotation=45)
    ax.set_ylim(99.9, 100.1)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig2_qc_metrics.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig2_qc_metrics.png")


def fig3_failed_samples(qc_df):
    """Figure 3: Failed/warned samples analysis."""
    print("Generating Figure 3: Failed Samples Analysis...")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Get problem samples
    problem_samples = qc_df[qc_df['qc_status'].isin(['WARN', 'FAIL'])]
    passing_samples = qc_df[qc_df['qc_status'] == 'PASS']

    # 3a: Comparison of key metrics
    ax = axes[0]
    metrics = ['median_counts', 'median_genes', 'pct_zero_count_cells']
    x = np.arange(len(metrics))
    width = 0.35

    pass_means = [passing_samples[m].mean() for m in metrics]
    prob_means = [problem_samples[m].mean() for m in metrics]

    # Normalize for visualization
    pass_norm = [v / max(pass_means[i], prob_means[i]) for i, v in enumerate(pass_means)]
    prob_norm = [v / max(pass_means[i], prob_means[i]) for i, v in enumerate(prob_means)]

    bars1 = ax.bar(x - width/2, pass_norm, width, label='PASS', color=QC_COLORS['PASS'], alpha=0.8)
    bars2 = ax.bar(x + width/2, prob_norm, width, label='WARN/FAIL', color=QC_COLORS['FAIL'], alpha=0.8)

    ax.set_ylabel('Normalized Value')
    ax.set_title('A) Metric Comparison: PASS vs Problem Samples')
    ax.set_xticks(x)
    ax.set_xticklabels(['Median\nTranscripts', 'Median\nGenes', 'Empty\nCells %'])
    ax.legend()
    ax.set_ylim(0, 1.2)

    # 3b: H04 (failed) vs other cancer samples
    ax = axes[1]
    cancer_samples = qc_df[qc_df['stage'] == 'cancer']
    h04 = cancer_samples[cancer_samples['sample_id'] == 'H04']
    other_cancer = cancer_samples[cancer_samples['sample_id'] != 'H04']

    metrics_compare = ['median_counts', 'median_genes']
    x = np.arange(len(metrics_compare))

    other_means = [other_cancer[m].mean() for m in metrics_compare]
    h04_vals = [h04[m].values[0] for m in metrics_compare]

    bars1 = ax.bar(x - width/2, other_means, width, label='Other Cancer', color='#3498db', alpha=0.8)
    bars2 = ax.bar(x + width/2, h04_vals, width, label='H04 (Failed)', color=QC_COLORS['FAIL'], alpha=0.8)

    ax.set_ylabel('Value')
    ax.set_title('B) H04 vs Other Cancer Samples')
    ax.set_xticks(x)
    ax.set_xticklabels(['Median Transcripts', 'Median Genes'])
    ax.legend()

    # Add threshold lines
    ax.axhline(y=30, color='red', linestyle='--', alpha=0.3)
    ax.axhline(y=20, color='red', linestyle='--', alpha=0.3)

    # 3c: D01 warning analysis
    ax = axes[2]
    d01 = qc_df[qc_df['sample_id'] == 'D01']
    other_pass = qc_df[(qc_df['qc_status'] == 'PASS')]

    # Show empty cell % comparison
    categories = ['D01 (WARN)', 'PASS Mean', 'Threshold']
    values = [
        d01['pct_zero_count_cells'].values[0],
        other_pass['pct_zero_count_cells'].mean(),
        5.0
    ]
    colors_bar = [QC_COLORS['WARN'], QC_COLORS['PASS'], 'gray']

    bars = ax.bar(categories, values, color=colors_bar, alpha=0.8)
    ax.set_ylabel('Empty Cell %')
    ax.set_title('C) D01 Empty Cell Rate Analysis')

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
                f'{val:.1f}%', ha='center', va='bottom')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig3_failed_samples.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig3_failed_samples.png")


def fig4_batch_effects(adata):
    """Figure 4: Batch effect visualization - UMAP by lane and stage."""
    print("Generating Figure 4: Batch Effects...")

    if adata is None:
        print("  Skipping - no merged data available")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Subsample for plotting (too many cells)
    np.random.seed(42)
    n_plot = min(50000, adata.n_obs)
    idx = np.random.choice(adata.n_obs, n_plot, replace=False)
    adata_sub = adata[idx].copy()

    # 4a: UMAP colored by lane
    ax = axes[0]
    for lane in sorted(adata_sub.obs['lane'].unique()):
        mask = adata_sub.obs['lane'] == lane
        ax.scatter(adata_sub.obsm['X_umap'][mask, 0],
                   adata_sub.obsm['X_umap'][mask, 1],
                   c=LANE_COLORS.get(lane, 'gray'), s=1, alpha=0.3, label=lane)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title('A) UMAP by Lane (Post-Harmony)')
    ax.legend(markerscale=10, loc='upper right')

    # 4b: UMAP colored by stage
    ax = axes[1]
    for stage in ['control', 'normal', 'metaplasia', 'cancer']:
        mask = adata_sub.obs['stage'] == stage
        ax.scatter(adata_sub.obsm['X_umap'][mask, 0],
                   adata_sub.obsm['X_umap'][mask, 1],
                   c=STAGE_COLORS[stage], s=1, alpha=0.3, label=stage)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title('B) UMAP by Disease Stage (Post-Harmony)')
    ax.legend(markerscale=10, loc='upper right')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig4_batch_effects.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig4_batch_effects.png")


def fig5_lane_comparison(qc_df):
    """Figure 5: Lane-to-lane comparison."""
    print("Generating Figure 5: Lane Comparison...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    lanes = ['L001', 'L002', 'L003', 'L004']

    # 5a: Cells per lane by stage
    ax = axes[0, 0]
    stage_order = ['control', 'normal', 'metaplasia', 'cancer']
    x = np.arange(len(lanes))
    width = 0.2

    for i, stage in enumerate(stage_order):
        vals = [qc_df[(qc_df['lane'] == lane) & (qc_df['stage'] == stage)]['n_cells'].sum()
                for lane in lanes]
        ax.bar(x + i*width, vals, width, label=stage, color=STAGE_COLORS[stage], alpha=0.8)

    ax.set_xlabel('Lane')
    ax.set_ylabel('Total Cells')
    ax.set_title('A) Cells per Lane by Stage')
    ax.set_xticks(x + width*1.5)
    ax.set_xticklabels(lanes)
    ax.legend()

    # 5b: Median transcripts by lane
    ax = axes[0, 1]
    bp = ax.boxplot([qc_df[qc_df['lane'] == lane]['median_counts'].values for lane in lanes],
                    labels=lanes, patch_artist=True)
    for patch, lane in zip(bp['boxes'], lanes):
        patch.set_facecolor(LANE_COLORS[lane])
        patch.set_alpha(0.7)
    ax.set_xlabel('Lane')
    ax.set_ylabel('Median Transcripts per Cell')
    ax.set_title('B) Transcript Distribution by Lane')

    # 5c: QC status by lane
    ax = axes[1, 0]
    qc_by_lane = qc_df.groupby(['lane', 'qc_status']).size().unstack(fill_value=0)
    qc_by_lane = qc_by_lane.reindex(columns=['PASS', 'WARN', 'FAIL'], fill_value=0)

    bottom = np.zeros(len(lanes))
    for status in ['PASS', 'WARN', 'FAIL']:
        if status in qc_by_lane.columns:
            vals = qc_by_lane[status].values
            ax.bar(lanes, vals, bottom=bottom, label=status, color=QC_COLORS[status], alpha=0.8)
            bottom += vals

    ax.set_xlabel('Lane')
    ax.set_ylabel('Number of Samples')
    ax.set_title('C) QC Status by Lane')
    ax.legend()

    # 5d: Patient distribution by lane
    ax = axes[1, 1]
    patient_by_lane = qc_df.groupby(['lane', 'patient']).size().unstack(fill_value=0)
    patient_by_lane.plot(kind='bar', ax=ax, alpha=0.8)
    ax.set_xlabel('Lane')
    ax.set_ylabel('Number of Samples')
    ax.set_title('D) Patient Distribution by Lane')
    ax.tick_params(axis='x', rotation=0)
    ax.legend(title='Patient', bbox_to_anchor=(1.02, 1))

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig5_lane_comparison.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig5_lane_comparison.png")


def fig6_summary_stats(qc_df):
    """Figure 6: Summary statistics table as figure."""
    print("Generating Figure 6: Summary Statistics...")

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')

    # Calculate summary stats
    stats = {
        'Metric': [
            'Total Samples',
            'Samples Passed QC',
            'Samples Warned',
            'Samples Failed',
            'Total Cells (Raw)',
            'Cells After QC',
            'Genes in Panel',
            'Proteins in Panel',
            'Median Transcripts/Cell',
            'Median Genes/Cell',
            'LISI Score (Pre)',
            'LISI Score (Post)',
            'Batch Correction'
        ],
        'Value': [
            '32',
            '30',
            '1 (D01)',
            '1 (H04)',
            '2,308,968',
            '1,835,026',
            '341',
            '17',
            f"{qc_df['median_counts'].median():.0f}",
            f"{qc_df['median_genes'].median():.0f}",
            '2.46',
            '2.66',
            'Harmony (GPU)'
        ]
    }

    df_stats = pd.DataFrame(stats)

    # Create table
    table = ax.table(cellText=df_stats.values,
                     colLabels=df_stats.columns,
                     cellLoc='center',
                     loc='center',
                     colWidths=[0.5, 0.3])

    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.2, 1.8)

    # Style header
    for i in range(len(df_stats.columns)):
        table[(0, i)].set_facecolor('#3498db')
        table[(0, i)].set_text_props(color='white', fontweight='bold')

    # Alternate row colors
    for i in range(1, len(df_stats) + 1):
        for j in range(len(df_stats.columns)):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#ecf0f1')

    ax.set_title('G4X Choi Batch 2 - QC Summary Statistics', fontsize=14, fontweight='bold', pad=20)

    plt.savefig(FIGURES_DIR / "fig6_summary_stats.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig6_summary_stats.png")


def main():
    """Generate all figures."""
    print("=" * 60)
    print("G4X QC Report Figure Generation")
    print("=" * 60)

    # Load data
    qc_df, adata = load_data()

    # Generate figures
    fig1_sample_overview(qc_df)
    fig2_qc_metrics(qc_df)
    fig3_failed_samples(qc_df)
    fig4_batch_effects(adata)
    fig5_lane_comparison(qc_df)
    fig6_summary_stats(qc_df)

    print("=" * 60)
    print(f"All figures saved to: {FIGURES_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
