#!/usr/bin/env python3
"""
04_exploratory_deg.py - Exploratory DEG analysis

IMPORTANT CAVEATS:
- N = 1 patient (all 3 samples from same individual)
- P-values are NOT meaningful for group comparisons
- Rankings by fold-change and effect size ONLY
- Findings require validation in independent cohorts

Analysis:
1. Stage comparisons (all cells): Normal vs Metaplasia vs Cancer
2. Cell-type specific: Within each major cell type, compare stages
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# === CONFIG ===
INPUT = Path('/home/user/g4x-choi-batch2-analysis/results/pilot')
OUTPUT = INPUT
FIGURES = OUTPUT / 'figures'
FIGURES.mkdir(exist_ok=True)

SAMPLES = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}
STAGE_ORDER = ['Normal', 'Metaplasia', 'Cancer']

# Annotation method to use
ANNOTATION_METHODS = ['celltype_markers', 'celltype_celltypist', 'celltype_gating']

# DEG parameters
MIN_CELLS_PER_GROUP = 50  # Minimum cells to include comparison
TOP_DEGS = 50  # Top genes to report per comparison


def get_annotation_col(adata):
    """Get best available annotation column."""
    for method in ANNOTATION_METHODS:
        if method in adata.obs.columns:
            return method
    return None


def calculate_fold_change(adata, group_col, group1, group2):
    """
    Calculate log2 fold change between two groups.

    Returns DataFrame with gene, log2FC, mean1, mean2, pct1, pct2.
    NOTE: P-values from Wilcoxon are included but should NOT be interpreted
    as meaningful given n=1 patient.
    """
    mask1 = adata.obs[group_col] == group1
    mask2 = adata.obs[group_col] == group2

    n1, n2 = mask1.sum(), mask2.sum()
    if n1 < MIN_CELLS_PER_GROUP or n2 < MIN_CELLS_PER_GROUP:
        return None

    # Get expression values
    X = adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()

    X1 = X[mask1]
    X2 = X[mask2]

    results = []
    for i, gene in enumerate(adata.var_names):
        expr1 = X1[:, i]
        expr2 = X2[:, i]

        mean1 = np.mean(expr1)
        mean2 = np.mean(expr2)

        # Percent expressed
        pct1 = (expr1 > 0).mean() * 100
        pct2 = (expr2 > 0).mean() * 100

        # Log2 fold change (add pseudocount)
        log2fc = np.log2((mean2 + 0.01) / (mean1 + 0.01))

        # Wilcoxon test (for ranking, not significance!)
        try:
            _, pval = stats.mannwhitneyu(expr1, expr2, alternative='two-sided')
        except:
            pval = 1.0

        results.append({
            'gene': gene,
            'log2FC': log2fc,
            f'mean_{group1}': mean1,
            f'mean_{group2}': mean2,
            f'pct_{group1}': pct1,
            f'pct_{group2}': pct2,
            'pval_wilcoxon': pval,  # For ranking only!
            'n_cells_1': n1,
            'n_cells_2': n2,
        })

    return pd.DataFrame(results)


def create_volcano_plot(df, title, output_path, group1, group2):
    """
    Create volcano-style plot.
    NOTE: Using -log10(pval) for y-axis but p-values are NOT meaningful!
    This is purely for visualization.
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Color by fold change direction
    colors = np.where(df['log2FC'] > 1, '#E74C3C',  # Up in group2
              np.where(df['log2FC'] < -1, '#3498DB',  # Up in group1
              '#95A5A6'))  # Not changed

    ax.scatter(df['log2FC'], -np.log10(df['pval_wilcoxon'] + 1e-300),
               c=colors, alpha=0.6, s=20, rasterized=True)

    # Label top genes by |log2FC|
    top_genes = pd.concat([df.nlargest(15, 'log2FC'), df.nsmallest(10, 'log2FC')])

    for _, row in top_genes.iterrows():
        ax.annotate(row['gene'],
                    (row['log2FC'], -np.log10(row['pval_wilcoxon'] + 1e-300)),
                    fontsize=8, alpha=0.8)

    ax.axvline(1, ls='--', color='gray', alpha=0.5)
    ax.axvline(-1, ls='--', color='gray', alpha=0.5)

    ax.set_xlabel(f'log2 Fold Change ({group2} vs {group1})')
    ax.set_ylabel('-log10(p-value) [FOR RANKING ONLY]')
    ax.set_title(f'{title}\nN=1 patient - p-values not meaningful', fontweight='bold')

    # Add warning box
    ax.text(0.02, 0.98, 'EXPLORATORY ONLY\nn=1 patient',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))

    plt.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close()


def main():
    print("="*60)
    print("Exploratory DEG Analysis")
    print("="*60)
    print("\nCAVEAT: N=1 patient - p-values are NOT meaningful!")
    print("Rankings by fold-change only.\n")

    # Load and merge all samples
    print("Loading annotated samples...")
    adatas = []
    for sample_id, stage in SAMPLES.items():
        path = INPUT / f'{sample_id}_annotated.h5ad'
        if path.exists():
            adata = sc.read_h5ad(path)
            adata.obs['stage'] = stage
            adatas.append(adata)
            print(f"  {sample_id} ({stage}): {adata.n_obs:,} cells")
        else:
            print(f"  {sample_id}: NOT FOUND")

    if len(adatas) != 3:
        print("ERROR: Need all 3 samples for comparison")
        return

    # Merge
    merged = sc.concat(adatas, join='inner')
    merged.obs_names_make_unique()
    print(f"\nMerged: {merged.n_obs:,} cells x {merged.n_vars} genes")

    annotation_col = get_annotation_col(merged)
    print(f"Using annotation: {annotation_col}")

    all_degs = []

    # === Global Stage Comparisons ===
    print("\n" + "="*50)
    print("GLOBAL STAGE COMPARISONS (all cells)")
    print("="*50)

    stage_comparisons = [
        ('Normal', 'Metaplasia'),
        ('Normal', 'Cancer'),
        ('Metaplasia', 'Cancer'),
    ]

    for stage1, stage2 in stage_comparisons:
        print(f"\n{stage1} vs {stage2}:")
        df = calculate_fold_change(merged, 'stage', stage1, stage2)

        if df is None:
            print("  Insufficient cells")
            continue

        # Sort by absolute fold change
        df['abs_log2FC'] = df['log2FC'].abs()
        df = df.sort_values('abs_log2FC', ascending=False)

        # Top upregulated in stage2
        up = df[df['log2FC'] > 0].head(10)
        print(f"  Top up in {stage2}:")
        for _, row in up.head(5).iterrows():
            print(f"    {row['gene']}: log2FC={row['log2FC']:.2f}")

        # Top upregulated in stage1 (negative FC)
        down = df[df['log2FC'] < 0].head(10)
        print(f"  Top up in {stage1}:")
        for _, row in down.head(5).iterrows():
            print(f"    {row['gene']}: log2FC={row['log2FC']:.2f}")

        # Save
        df['comparison'] = f'{stage2}_vs_{stage1}'
        df['comparison_type'] = 'global'
        all_degs.append(df)

        # Volcano plot
        create_volcano_plot(
            df, f'Global: {stage2} vs {stage1}',
            FIGURES / f'deg_volcano_global_{stage2}_vs_{stage1}.png',
            stage1, stage2
        )

    # === Cell-Type Specific Comparisons ===
    if annotation_col:
        print("\n" + "="*50)
        print("CELL-TYPE SPECIFIC COMPARISONS")
        print("="*50)

        major_types = merged.obs[annotation_col].value_counts()
        major_types = major_types[major_types >= MIN_CELLS_PER_GROUP * 2].index
        major_types = [t for t in major_types if t not in ['Unknown', 'Not_Run', 'Failed']]

        for cell_type in major_types[:6]:  # Top 6 cell types
            print(f"\n{cell_type}:")
            subset = merged[merged.obs[annotation_col] == cell_type].copy()

            if subset.n_obs < MIN_CELLS_PER_GROUP * 2:
                print("  Insufficient cells")
                continue

            for stage1, stage2 in [('Normal', 'Cancer')]:  # Focus comparison
                df = calculate_fold_change(subset, 'stage', stage1, stage2)

                if df is None:
                    print(f"  {stage1} vs {stage2}: Insufficient cells")
                    continue

                df['abs_log2FC'] = df['log2FC'].abs()
                df = df.sort_values('abs_log2FC', ascending=False)

                print(f"  {stage1} vs {stage2} (n={df['n_cells_1'].iloc[0]:,} vs {df['n_cells_2'].iloc[0]:,}):")
                for _, row in df.head(3).iterrows():
                    direction = '↑' if row['log2FC'] > 0 else '↓'
                    print(f"    {row['gene']}: {direction} log2FC={row['log2FC']:.2f}")

                df['comparison'] = f'{cell_type}_{stage2}_vs_{stage1}'
                df['comparison_type'] = 'cell_type_specific'
                df['cell_type'] = cell_type
                all_degs.append(df)

    # === Save Results ===
    if all_degs:
        deg_df = pd.concat(all_degs, ignore_index=True)

        # Save full results
        deg_df.to_csv(OUTPUT / 'deg_exploratory_full.csv', index=False)
        print(f"\nSaved: deg_exploratory_full.csv")

        # Save top DEGs summary
        top_degs = deg_df.nlargest(TOP_DEGS * len(all_degs), 'abs_log2FC')
        top_degs.to_csv(OUTPUT / 'deg_exploratory_top.csv', index=False)
        print(f"Saved: deg_exploratory_top.csv")

        # === Create Summary Heatmap ===
        print("\nCreating summary heatmap...")
        create_deg_heatmap(deg_df, FIGURES / 'deg_heatmap_top.png')

    print(f"\nAll results saved to: {OUTPUT}")


def create_deg_heatmap(deg_df, output_path):
    """Create heatmap of top DEGs across comparisons."""
    # Get global comparisons only for clarity
    global_degs = deg_df[deg_df['comparison_type'] == 'global'].copy()

    if global_degs.empty:
        return

    # Get top genes by max |log2FC| across any comparison
    gene_max_fc = global_degs.groupby('gene')['abs_log2FC'].max().nlargest(30).index

    # Pivot
    pivot_data = global_degs[global_degs['gene'].isin(gene_max_fc)]
    pivot = pivot_data.pivot(index='gene', columns='comparison', values='log2FC')

    # Reorder
    pivot = pivot.reindex(index=gene_max_fc)

    fig, ax = plt.subplots(figsize=(10, 12))
    sns.heatmap(pivot, cmap='RdBu_r', center=0, ax=ax,
                cbar_kws={'label': 'log2 Fold Change'},
                annot=True, fmt='.1f', annot_kws={'size': 8})

    ax.set_title('Top DEGs Across Stage Comparisons\n(EXPLORATORY - N=1 patient)',
                 fontweight='bold')
    ax.set_xlabel('Comparison')
    ax.set_ylabel('Gene')

    plt.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close()


if __name__ == '__main__':
    main()
