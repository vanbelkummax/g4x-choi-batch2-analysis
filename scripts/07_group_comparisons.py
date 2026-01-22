#!/usr/bin/env python3
"""
G4X Choi Batch 2 - Between-Group Statistical Comparisons
=========================================================

Compares samples across defined groups with statistical testing.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configuration
OUTPUT_DIR = Path('/home/user/g4x-choi-batch2-analysis/results')
FIG_DIR = OUTPUT_DIR / 'figures'
TABLE_DIR = OUTPUT_DIR / 'tables'

def load_analysis_data():
    """Load hypothesis analysis results."""
    df = pd.read_csv(TABLE_DIR / 'hypothesis_analysis_results.csv')

    # Extract lane from sample_id (last digit)
    df['lane'] = df['sample_id'].str[-1].astype(int)

    return df

def define_groups(df):
    """Define comparison groups based on biological thresholds."""

    # CD8/Treg ratio groups (median split)
    median_ratio = df['cd8_treg_ratio'].median()
    df['cd8_treg_group'] = np.where(df['cd8_treg_ratio'] >= median_ratio, 'High_CD8/Treg', 'Low_CD8/Treg')

    # PDL1+ tumor groups (median split)
    median_pdl1 = df['pct_pdl1_pos_tumor'].median()
    df['pdl1_group'] = np.where(df['pct_pdl1_pos_tumor'] >= median_pdl1, 'High_PDL1', 'Low_PDL1')

    # Immune infiltration groups
    df['immune_group'] = pd.cut(
        df['total_immune_pct'],
        bins=[0, 18, 22, 100],
        labels=['Cold', 'Intermediate', 'Hot']
    )

    # PD1+ T cell exhaustion groups
    median_pd1 = df['pct_pd1_pos_t_cells'].median()
    df['exhaustion_group'] = np.where(df['pct_pd1_pos_t_cells'] >= median_pd1, 'High_Exhaustion', 'Low_Exhaustion')

    return df

def mann_whitney_comparison(group1_values, group2_values):
    """Perform Mann-Whitney U test with effect size."""
    # Remove NaN
    g1 = group1_values.dropna()
    g2 = group2_values.dropna()

    if len(g1) < 3 or len(g2) < 3:
        return np.nan, np.nan, np.nan

    stat, p = stats.mannwhitneyu(g1, g2, alternative='two-sided')

    # Rank-biserial correlation as effect size
    n1, n2 = len(g1), len(g2)
    r = 1 - (2 * stat) / (n1 * n2)

    return stat, p, r

def kruskal_wallis_comparison(groups_dict):
    """Perform Kruskal-Wallis test across multiple groups."""
    group_values = [v.dropna().values for v in groups_dict.values() if len(v.dropna()) >= 3]

    if len(group_values) < 2:
        return np.nan, np.nan

    stat, p = stats.kruskal(*group_values)
    return stat, p

def compare_by_group(df, group_var, metrics):
    """Compare metrics between groups."""
    results = []

    groups = df[group_var].unique()
    groups = [g for g in groups if pd.notna(g)]

    if len(groups) == 2:
        # Two-group comparison (Mann-Whitney)
        g1_name, g2_name = sorted(groups)
        g1_data = df[df[group_var] == g1_name]
        g2_data = df[df[group_var] == g2_name]

        for metric in metrics:
            stat, p, effect = mann_whitney_comparison(g1_data[metric], g2_data[metric])
            results.append({
                'group_var': group_var,
                'group1': g1_name,
                'group2': g2_name,
                'metric': metric,
                'group1_mean': g1_data[metric].mean(),
                'group1_std': g1_data[metric].std(),
                'group2_mean': g2_data[metric].mean(),
                'group2_std': g2_data[metric].std(),
                'statistic': stat,
                'p_value': p,
                'effect_size': effect,
                'test': 'Mann-Whitney U'
            })
    else:
        # Multi-group comparison (Kruskal-Wallis)
        for metric in metrics:
            groups_dict = {g: df[df[group_var] == g][metric] for g in groups}
            stat, p = kruskal_wallis_comparison(groups_dict)

            # Calculate means for each group
            means = {g: df[df[group_var] == g][metric].mean() for g in groups}

            results.append({
                'group_var': group_var,
                'groups': ', '.join(map(str, sorted(groups))),
                'metric': metric,
                'group_means': str(means),
                'statistic': stat,
                'p_value': p,
                'test': 'Kruskal-Wallis'
            })

    return pd.DataFrame(results)

def plot_group_comparisons(df, output_path):
    """Generate comprehensive comparison figure."""

    fig, axes = plt.subplots(3, 3, figsize=(18, 18))

    # Key metrics to compare
    metrics = ['cd8_treg_ratio', 'pct_pdl1_pos_tumor', 'pct_pd1_pos_t_cells',
               'macrophage_pct', 'stromal_pct', 'epithelial_pct']

    # Row 1: Lane comparisons (batch effects)
    for i, metric in enumerate(['total_immune_pct', 'cd8_treg_ratio', 'pct_pdl1_pos_tumor']):
        ax = axes[0, i]
        sns.boxplot(data=df, x='lane', y=metric, ax=ax, palette='Set2')
        ax.set_title(f'{metric} by Lane', fontweight='bold')
        ax.set_xlabel('Lane')

        # Add Kruskal-Wallis p-value
        lane_groups = {l: df[df['lane'] == l][metric] for l in df['lane'].unique()}
        _, p = kruskal_wallis_comparison(lane_groups)
        ax.text(0.02, 0.98, f'KW p = {p:.3f}', transform=ax.transAxes, va='top', fontsize=10,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    axes[0, 0].set_ylabel('Percentage')

    # Row 2: CD8/Treg group comparisons
    for i, metric in enumerate(['pct_pdl1_pos_tumor', 'pct_pd1_pos_t_cells', 'macrophage_pct']):
        ax = axes[1, i]
        sns.boxplot(data=df, x='cd8_treg_group', y=metric, ax=ax,
                   order=['Low_CD8/Treg', 'High_CD8/Treg'], palette='coolwarm')
        ax.set_title(f'{metric} by CD8/Treg ratio', fontweight='bold')
        ax.set_xlabel('')

        # Add Mann-Whitney p-value
        low = df[df['cd8_treg_group'] == 'Low_CD8/Treg'][metric]
        high = df[df['cd8_treg_group'] == 'High_CD8/Treg'][metric]
        _, p, _ = mann_whitney_comparison(low, high)
        ax.text(0.02, 0.98, f'MW p = {p:.3f}', transform=ax.transAxes, va='top', fontsize=10,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Row 3: PDL1 group comparisons
    for i, metric in enumerate(['cd8_treg_ratio', 'pct_pd1_pos_t_cells', 'total_immune_pct']):
        ax = axes[2, i]
        sns.boxplot(data=df, x='pdl1_group', y=metric, ax=ax,
                   order=['Low_PDL1', 'High_PDL1'], palette='RdYlGn_r')
        ax.set_title(f'{metric} by PDL1 status', fontweight='bold')
        ax.set_xlabel('')

        # Add Mann-Whitney p-value
        low = df[df['pdl1_group'] == 'Low_PDL1'][metric]
        high = df[df['pdl1_group'] == 'High_PDL1'][metric]
        _, p, _ = mann_whitney_comparison(low, high)
        ax.text(0.02, 0.98, f'MW p = {p:.3f}', transform=ax.transAxes, va='top', fontsize=10,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.suptitle('Between-Group Comparisons (N=32 samples)', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_path}")

def plot_correlation_matrix(df, output_path):
    """Plot correlation matrix of key metrics."""

    numeric_cols = ['n_cells', 'total_immune_pct', 'cd8_pct', 'treg_pct', 'cd8_treg_ratio',
                   'macrophage_pct', 'epithelial_pct', 'stromal_pct',
                   'pct_pd1_pos_t_cells', 'pct_pdl1_pos_tumor']

    corr_matrix = df[numeric_cols].corr()

    fig, ax = plt.subplots(figsize=(12, 10))
    mask = np.triu(np.ones_like(corr_matrix), k=1)

    sns.heatmap(corr_matrix, mask=mask, annot=True, fmt='.2f',
                cmap='coolwarm', center=0, ax=ax,
                cbar_kws={'label': 'Correlation'})

    ax.set_title('Correlation Matrix of Key Metrics (N=32)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_path}")

def main():
    """Run all group comparisons."""
    print("="*60)
    print("G4X Choi Batch 2 - Between-Group Comparisons")
    print("="*60)

    # Load data
    df = load_analysis_data()
    print(f"Loaded {len(df)} samples")

    # Define groups
    df = define_groups(df)
    print("\nGroup definitions:")
    print(f"  CD8/Treg groups: {df['cd8_treg_group'].value_counts().to_dict()}")
    print(f"  PDL1 groups: {df['pdl1_group'].value_counts().to_dict()}")
    print(f"  Immune groups: {df['immune_group'].value_counts().to_dict()}")
    print(f"  Exhaustion groups: {df['exhaustion_group'].value_counts().to_dict()}")

    # Metrics to compare
    metrics = ['n_cells', 'total_immune_pct', 'cd8_pct', 'treg_pct', 'cd8_treg_ratio',
               'macrophage_pct', 'epithelial_pct', 'stromal_pct',
               'pct_pd1_pos_t_cells', 'pct_pdl1_pos_tumor', 'cd8_exhaustion_score']

    # Run comparisons
    all_results = []

    print("\n--- Lane Comparisons (Batch Effects) ---")
    lane_results = compare_by_group(df, 'lane', metrics)
    all_results.append(lane_results)
    sig_lanes = lane_results[lane_results['p_value'] < 0.05]
    if len(sig_lanes) > 0:
        print(f"  Significant lane effects (p<0.05): {len(sig_lanes)} metrics")
        for _, row in sig_lanes.iterrows():
            print(f"    - {row['metric']}: p={row['p_value']:.4f}")
    else:
        print("  No significant lane effects (good!)")

    print("\n--- CD8/Treg Ratio Comparisons ---")
    cd8_results = compare_by_group(df, 'cd8_treg_group', metrics)
    all_results.append(cd8_results)
    sig_cd8 = cd8_results[cd8_results['p_value'] < 0.05]
    if len(sig_cd8) > 0:
        print(f"  Significant differences: {len(sig_cd8)} metrics")
        for _, row in sig_cd8.iterrows():
            print(f"    - {row['metric']}: p={row['p_value']:.4f}, effect={row['effect_size']:.2f}")

    print("\n--- PDL1 Status Comparisons ---")
    pdl1_results = compare_by_group(df, 'pdl1_group', metrics)
    all_results.append(pdl1_results)
    sig_pdl1 = pdl1_results[pdl1_results['p_value'] < 0.05]
    if len(sig_pdl1) > 0:
        print(f"  Significant differences: {len(sig_pdl1)} metrics")
        for _, row in sig_pdl1.iterrows():
            print(f"    - {row['metric']}: p={row['p_value']:.4f}, effect={row['effect_size']:.2f}")

    print("\n--- Exhaustion Group Comparisons ---")
    exh_results = compare_by_group(df, 'exhaustion_group', metrics)
    all_results.append(exh_results)
    sig_exh = exh_results[exh_results['p_value'] < 0.05]
    if len(sig_exh) > 0:
        print(f"  Significant differences: {len(sig_exh)} metrics")
        for _, row in sig_exh.iterrows():
            print(f"    - {row['metric']}: p={row['p_value']:.4f}, effect={row['effect_size']:.2f}")

    # Save results
    combined_results = pd.concat(all_results, ignore_index=True)
    output_csv = TABLE_DIR / 'group_comparison_results.csv'
    combined_results.to_csv(output_csv, index=False)
    print(f"\nSaved results to: {output_csv}")

    # Generate figures
    print("\nGenerating figures...")
    plot_group_comparisons(df, FIG_DIR / 'group_comparisons.png')
    plot_correlation_matrix(df, FIG_DIR / 'correlation_matrix.png')

    # Save grouped data
    df.to_csv(TABLE_DIR / 'samples_with_groups.csv', index=False)
    print(f"Saved grouped data to: {TABLE_DIR / 'samples_with_groups.csv'}")

    print("\n" + "="*60)
    print("COMPLETE")
    print("="*60)

if __name__ == '__main__':
    main()
