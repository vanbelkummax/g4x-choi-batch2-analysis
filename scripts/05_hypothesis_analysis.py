#!/usr/bin/env python3
"""
G4X Choi Batch 2 - Hypothesis-Driven Analysis
==============================================

Targeted analysis of biological hypotheses:
1. Immune infiltration patterns
2. Checkpoint expression (PD1/PDL1)
3. T cell exhaustion signatures
4. TME heterogeneity across samples
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
from scipy.stats import spearmanr, mannwhitneyu, kruskal
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# Configuration
DATA_ROOT = Path('/mnt/x/Choi_Batch_2_Tuesday')
OUTPUT_DIR = Path('/home/user/g4x-choi-batch2-analysis/results')
FIG_DIR = OUTPUT_DIR / 'figures'
TABLE_DIR = OUTPUT_DIR / 'tables'

def load_sample_data(lane_idx, sample_prefix):
    """Load cell metadata with protein intensities."""
    lane = f'L00{lane_idx + 1}'
    sample = f'{sample_prefix}{lane_idx + 1}'
    lane_dir = list(DATA_ROOT.glob(f'g4-028-083-FC1-{lane}_*'))[0]
    sample_path = lane_dir / sample / 'single_cell_data' / 'cell_metadata.csv.gz'

    if not sample_path.exists():
        return None

    with gzip.open(sample_path, 'rt') as f:
        df = pd.read_csv(f)

    return df

def gate_cell_types(df):
    """Gate cell types (same logic as other scripts)."""
    quantile_thresh = 0.8
    protein_markers = ['PD1', 'CD8', 'PDL1', 'FOXP3', 'KI67', 'CD45', 'PanCK',
                       'CD20', 'CD4', 'CD11c', 'HLA-DR', 'CD68', 'CD31', 'CD3', 'aSMA']

    thresholds = {}
    for marker in protein_markers:
        col = f'{marker}_intensity_mean'
        if col in df.columns:
            thresholds[marker] = df[col].quantile(quantile_thresh)

    def get_val(marker):
        col = f'{marker}_intensity_mean'
        return df[col].values if col in df.columns else np.zeros(len(df))

    def get_t(marker):
        return thresholds.get(marker, 0)

    cell_types = np.full(len(df), 'Other', dtype=object)

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

def analyze_checkpoint_expression(df, cell_types):
    """Analyze PD1/PDL1 checkpoint expression."""

    df['cell_type'] = cell_types

    # PD1 expression by cell type
    pd1_by_type = df.groupby('cell_type')['PD1_intensity_mean'].agg(['mean', 'std', 'median'])

    # PDL1 expression by cell type
    pdl1_by_type = df.groupby('cell_type')['PDL1_intensity_mean'].agg(['mean', 'std', 'median'])

    # PD1+ T cells
    t_cells = df[df['cell_type'].isin(['CD8_T', 'CD4_T', 'T_cell', 'Treg'])]
    pd1_thresh = df['PD1_intensity_mean'].quantile(0.8)
    pd1_pos_t = (t_cells['PD1_intensity_mean'] > pd1_thresh).sum() / len(t_cells) * 100 if len(t_cells) > 0 else 0

    # PDL1+ tumor cells
    tumor_cells = df[df['cell_type'] == 'Epithelial']
    pdl1_thresh = df['PDL1_intensity_mean'].quantile(0.8)
    pdl1_pos_tumor = (tumor_cells['PDL1_intensity_mean'] > pdl1_thresh).sum() / len(tumor_cells) * 100 if len(tumor_cells) > 0 else 0

    return {
        'pd1_by_type': pd1_by_type,
        'pdl1_by_type': pdl1_by_type,
        'pct_pd1_pos_t_cells': pd1_pos_t,
        'pct_pdl1_pos_tumor': pdl1_pos_tumor
    }

def analyze_t_cell_exhaustion(df, cell_types):
    """Analyze T cell exhaustion markers."""

    df['cell_type'] = cell_types

    # CD8+ T cells
    cd8_t = df[df['cell_type'] == 'CD8_T']

    if len(cd8_t) < 10:
        return None

    # Calculate exhaustion score based on:
    # High PD1, low KI67 (proliferation)
    pd1_norm = (cd8_t['PD1_intensity_mean'] - cd8_t['PD1_intensity_mean'].min()) / \
               (cd8_t['PD1_intensity_mean'].max() - cd8_t['PD1_intensity_mean'].min() + 1e-10)
    ki67_norm = (cd8_t['KI67_intensity_mean'] - cd8_t['KI67_intensity_mean'].min()) / \
                (cd8_t['KI67_intensity_mean'].max() - cd8_t['KI67_intensity_mean'].min() + 1e-10)

    # Exhaustion score: high PD1, low proliferation
    exhaustion_score = pd1_norm * (1 - ki67_norm)

    # Classify into exhausted vs active
    exhausted = exhaustion_score > exhaustion_score.quantile(0.75)

    return {
        'mean_exhaustion_score': exhaustion_score.mean(),
        'pct_exhausted': exhausted.sum() / len(cd8_t) * 100,
        'mean_pd1': cd8_t['PD1_intensity_mean'].mean(),
        'mean_ki67': cd8_t['KI67_intensity_mean'].mean(),
        'n_cd8_t': len(cd8_t)
    }

def analyze_tme_composition(df, cell_types):
    """Analyze tumor microenvironment composition."""

    df['cell_type'] = cell_types
    type_counts = df['cell_type'].value_counts()
    type_pcts = type_counts / len(df) * 100

    # Calculate TME metrics
    immune_types = ['CD8_T', 'CD4_T', 'Treg', 'T_cell', 'B_cell', 'Macrophage', 'DC', 'Immune_other']
    total_immune = type_pcts[[t for t in immune_types if t in type_pcts.index]].sum()

    # CD8/Treg ratio (important for immunotherapy)
    cd8_pct = type_pcts.get('CD8_T', 0)
    treg_pct = type_pcts.get('Treg', 0)
    cd8_treg_ratio = cd8_pct / (treg_pct + 0.01)

    # Macrophage infiltration
    mac_pct = type_pcts.get('Macrophage', 0)

    # M1/M2 proxy using HLA-DR expression in macrophages
    macs = df[df['cell_type'] == 'Macrophage']
    if len(macs) > 10:
        hladr_high = (macs['HLA-DR_intensity_mean'] > macs['HLA-DR_intensity_mean'].median()).sum()
        m1_proxy = hladr_high / len(macs) * 100  # High HLA-DR = M1-like
    else:
        m1_proxy = None

    return {
        'total_immune_pct': total_immune,
        'cd8_pct': cd8_pct,
        'treg_pct': treg_pct,
        'cd8_treg_ratio': cd8_treg_ratio,
        'macrophage_pct': mac_pct,
        'm1_like_pct': m1_proxy,
        'epithelial_pct': type_pcts.get('Epithelial', 0),
        'stromal_pct': type_pcts.get('Stromal', 0)
    }

def analyze_sample(lane_idx, sample_prefix):
    """Run all hypothesis analyses for a sample."""

    df = load_sample_data(lane_idx, sample_prefix)
    if df is None:
        return None

    sample_id = f'{sample_prefix}{lane_idx + 1}'
    cell_types = gate_cell_types(df)

    # Run analyses
    checkpoint = analyze_checkpoint_expression(df, cell_types)
    exhaustion = analyze_t_cell_exhaustion(df, cell_types)
    tme = analyze_tme_composition(df, cell_types)

    result = {
        'sample_id': sample_id,
        'n_cells': len(df),
        **tme,
        'pct_pd1_pos_t_cells': checkpoint['pct_pd1_pos_t_cells'],
        'pct_pdl1_pos_tumor': checkpoint['pct_pdl1_pos_tumor']
    }

    if exhaustion is not None:
        result.update({
            'cd8_exhaustion_score': exhaustion['mean_exhaustion_score'],
            'pct_cd8_exhausted': exhaustion['pct_exhausted'],
            'n_cd8_t_cells': exhaustion['n_cd8_t']
        })

    return result

def plot_immune_infiltration(results_df, output_dir):
    """Plot immune infiltration patterns."""

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Immune composition bar chart
    ax = axes[0, 0]
    immune_cols = ['cd8_pct', 'cd8_treg_ratio', 'macrophage_pct', 'total_immune_pct']
    results_sorted = results_df.sort_values('total_immune_pct', ascending=False)

    x = range(len(results_sorted))
    ax.bar(x, results_sorted['total_immune_pct'], label='Total Immune', alpha=0.7)
    ax.bar(x, results_sorted['cd8_pct'], label='CD8+ T', alpha=0.7)
    ax.bar(x, results_sorted['treg_pct'], label='Treg', alpha=0.7)
    ax.set_xlabel('Sample (sorted by immune %)')
    ax.set_ylabel('Cell Type (%)')
    ax.set_title('Immune Cell Infiltration')
    ax.legend()
    ax.set_xticks(x[::4])
    ax.set_xticklabels(results_sorted['sample_id'].values[::4], rotation=45, ha='right')

    # 2. CD8/Treg ratio
    ax = axes[0, 1]
    ax.bar(range(len(results_df)), results_df['cd8_treg_ratio'].values)
    ax.axhline(y=1, color='r', linestyle='--', label='Ratio = 1')
    ax.set_xlabel('Sample')
    ax.set_ylabel('CD8/Treg Ratio')
    ax.set_title('CD8+ T / Treg Ratio')
    ax.legend()

    # 3. Checkpoint expression
    ax = axes[1, 0]
    ax.scatter(results_df['pct_pd1_pos_t_cells'], results_df['pct_pdl1_pos_tumor'],
               c=results_df['total_immune_pct'], cmap='viridis', s=50)
    ax.set_xlabel('PD1+ T cells (%)')
    ax.set_ylabel('PDL1+ Tumor cells (%)')
    ax.set_title('Checkpoint Expression')
    plt.colorbar(ax.collections[0], ax=ax, label='Total Immune %')

    # 4. TME composition heatmap
    ax = axes[1, 1]
    tme_cols = ['epithelial_pct', 'stromal_pct', 'cd8_pct', 'cd8_treg_ratio',
                'treg_pct', 'macrophage_pct', 'total_immune_pct']
    tme_cols = [c for c in tme_cols if c in results_df.columns]

    tme_data = results_df[tme_cols].values
    # Normalize for heatmap
    tme_norm = (tme_data - tme_data.mean(axis=0)) / (tme_data.std(axis=0) + 1e-10)

    im = ax.imshow(tme_norm.T, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_xlabel('Sample')
    ax.set_ylabel('TME Feature')
    ax.set_yticks(range(len(tme_cols)))
    ax.set_yticklabels([c.replace('_pct', '').replace('_', ' ') for c in tme_cols])
    ax.set_title('TME Heterogeneity (Z-score)')
    plt.colorbar(im, ax=ax, label='Z-score')

    plt.tight_layout()
    plt.savefig(output_dir / 'immune_infiltration.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'immune_infiltration.pdf', bbox_inches='tight')
    plt.close()

    print(f"Saved immune infiltration plots to {output_dir / 'immune_infiltration.png'}")

def plot_exhaustion_analysis(results_df, output_dir):
    """Plot T cell exhaustion analysis."""

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Filter samples with CD8 data
    cd8_df = results_df.dropna(subset=['cd8_exhaustion_score'])

    if len(cd8_df) == 0:
        print("No CD8+ T cell data available for exhaustion analysis")
        return

    # 1. Exhaustion score distribution
    ax = axes[0]
    ax.hist(cd8_df['cd8_exhaustion_score'], bins=15, edgecolor='black')
    ax.axvline(x=cd8_df['cd8_exhaustion_score'].median(), color='r', linestyle='--',
               label=f'Median: {cd8_df["cd8_exhaustion_score"].median():.3f}')
    ax.set_xlabel('Exhaustion Score')
    ax.set_ylabel('Number of Samples')
    ax.set_title('CD8+ T Cell Exhaustion Score')
    ax.legend()

    # 2. Exhaustion vs immune infiltration
    ax = axes[1]
    ax.scatter(cd8_df['total_immune_pct'], cd8_df['pct_cd8_exhausted'])
    corr, pval = spearmanr(cd8_df['total_immune_pct'], cd8_df['pct_cd8_exhausted'])
    ax.text(0.05, 0.95, f'r = {corr:.2f}\np = {pval:.3f}', transform=ax.transAxes,
            fontsize=10, verticalalignment='top')
    ax.set_xlabel('Total Immune Cells (%)')
    ax.set_ylabel('Exhausted CD8+ T Cells (%)')
    ax.set_title('Exhaustion vs Immune Infiltration')

    # 3. Sample ranking
    ax = axes[2]
    sorted_df = cd8_df.sort_values('pct_cd8_exhausted', ascending=False)
    colors = ['red' if x > 50 else 'orange' if x > 25 else 'green' for x in sorted_df['pct_cd8_exhausted']]
    ax.barh(range(len(sorted_df)), sorted_df['pct_cd8_exhausted'], color=colors)
    ax.set_yticks(range(len(sorted_df)))
    ax.set_yticklabels(sorted_df['sample_id'], fontsize=7)
    ax.set_xlabel('Exhausted CD8+ T Cells (%)')
    ax.set_title('Sample Ranking by T Cell Exhaustion')

    plt.tight_layout()
    plt.savefig(output_dir / 'exhaustion_analysis.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'exhaustion_analysis.pdf', bbox_inches='tight')
    plt.close()

    print(f"Saved exhaustion analysis to {output_dir / 'exhaustion_analysis.png'}")

def main():
    """Run hypothesis-driven analysis."""
    print("=" * 60)
    print("G4X Choi Batch 2 - Hypothesis-Driven Analysis")
    print("=" * 60)

    # Process all samples
    results = []

    print("\nAnalyzing samples...")
    for lane_idx in range(4):
        for sample_prefix in ['A0', 'B0', 'C0', 'D0', 'E0', 'F0', 'G0', 'H0']:
            sample_id = f'{sample_prefix}{lane_idx + 1}'
            print(f"  Processing {sample_id}...", end=' ')

            try:
                result = analyze_sample(lane_idx, sample_prefix)
                if result is not None:
                    results.append(result)
                    print(f"done")
                else:
                    print("skipped")
            except Exception as e:
                print(f"error: {e}")

    # Create results DataFrame
    results_df = pd.DataFrame(results)
    print(f"\nProcessed {len(results_df)} samples")

    # Save results
    results_df.to_csv(TABLE_DIR / 'hypothesis_analysis_results.csv', index=False)

    # Generate plots
    print("\n" + "-" * 40)
    print("Generating plots...")

    plot_immune_infiltration(results_df, FIG_DIR)
    plot_exhaustion_analysis(results_df, FIG_DIR)

    # Summary statistics
    print("\n" + "=" * 60)
    print("HYPOTHESIS-DRIVEN ANALYSIS SUMMARY")
    print("=" * 60)

    print("\n--- IMMUNE INFILTRATION ---")
    print(f"  Total Immune: {results_df['total_immune_pct'].mean():.1f}% ± {results_df['total_immune_pct'].std():.1f}%")
    print(f"  CD8+ T cells: {results_df['cd8_pct'].mean():.1f}% ± {results_df['cd8_pct'].std():.1f}%")
    print(f"  Tregs: {results_df['treg_pct'].mean():.1f}% ± {results_df['treg_pct'].std():.1f}%")
    print(f"  CD8/Treg ratio: {results_df['cd8_treg_ratio'].mean():.1f} ± {results_df['cd8_treg_ratio'].std():.1f}")

    print("\n--- CHECKPOINT EXPRESSION ---")
    print(f"  PD1+ T cells: {results_df['pct_pd1_pos_t_cells'].mean():.1f}% ± {results_df['pct_pd1_pos_t_cells'].std():.1f}%")
    print(f"  PDL1+ Tumor: {results_df['pct_pdl1_pos_tumor'].mean():.1f}% ± {results_df['pct_pdl1_pos_tumor'].std():.1f}%")

    if 'cd8_exhaustion_score' in results_df.columns:
        print("\n--- T CELL EXHAUSTION ---")
        ex_df = results_df.dropna(subset=['cd8_exhaustion_score'])
        if len(ex_df) > 0:
            print(f"  Exhaustion score: {ex_df['cd8_exhaustion_score'].mean():.3f} ± {ex_df['cd8_exhaustion_score'].std():.3f}")
            print(f"  % Exhausted CD8: {ex_df['pct_cd8_exhausted'].mean():.1f}% ± {ex_df['pct_cd8_exhausted'].std():.1f}%")

    print("\n--- SAMPLE HETEROGENEITY ---")
    print(f"  Samples with high immune (>25%): {(results_df['total_immune_pct'] > 25).sum()}")
    print(f"  Samples with high CD8/Treg (>5): {(results_df['cd8_treg_ratio'] > 5).sum()}")
    print(f"  Samples with high PDL1+ tumor (>30%): {(results_df['pct_pdl1_pos_tumor'] > 30).sum()}")

    # Identify interesting samples
    print("\n" + "-" * 40)
    print("NOTABLE SAMPLES:")

    # High immune
    high_immune = results_df.nlargest(3, 'total_immune_pct')
    print("\n  Highest immune infiltration:")
    for _, row in high_immune.iterrows():
        print(f"    {row['sample_id']}: {row['total_immune_pct']:.1f}%")

    # High CD8/Treg ratio (favorable)
    high_ratio = results_df.nlargest(3, 'cd8_treg_ratio')
    print("\n  Highest CD8/Treg ratio (favorable):")
    for _, row in high_ratio.iterrows():
        print(f"    {row['sample_id']}: {row['cd8_treg_ratio']:.1f}")

    # Low CD8/Treg ratio (unfavorable)
    low_ratio = results_df.nsmallest(3, 'cd8_treg_ratio')
    print("\n  Lowest CD8/Treg ratio (unfavorable):")
    for _, row in low_ratio.iterrows():
        print(f"    {row['sample_id']}: {row['cd8_treg_ratio']:.1f}")

    print("\n" + "=" * 60)
    print("Hypothesis-Driven Analysis Complete!")
    print(f"Results saved to: {TABLE_DIR}")
    print("=" * 60)

if __name__ == '__main__':
    main()
