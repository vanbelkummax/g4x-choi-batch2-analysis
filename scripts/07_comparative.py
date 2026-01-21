#!/usr/bin/env python3
"""
Spatial Biology Hackathon 2026 - Comparative Analysis
======================================================

Compare treatment responders (R) vs non-responders (NR) in PDAC:
- Cell type proportions
- Spatial organization (neighborhood enrichment)
- Differential gene expression
- Response signature identification

Key comparisons:
1. R vs NR at baseline (Pre)
2. R vs NR post-treatment (Post)
3. Pre vs Post in Responders
4. Pre vs Post in Non-Responders

Author: Max Van Belkum
Date: 2026-01-20
"""

import scanpy as sc
import squidpy as sq
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Optional, Dict, List, Tuple
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Configuration
# =============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
ADATA_DIR = OUTPUT_DIR / "adata"
FIG_DIR = OUTPUT_DIR / "figures" / "comparative"
TABLE_DIR = OUTPUT_DIR / "tables"

# Create directories
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Sample groupings
# NOTE: YP04A excluded (134 cells after QC - too sparse for reliable analysis)
SAMPLE_GROUPS = {
    'R_Pre': ['YP12A', 'YP15A'],
    'R_Post': ['YP12C', 'YP15C'],
    'NR_Pre': ['YP03A'],  # YP04A excluded
    'NR_Post': ['YP03C', 'YP04C'],
}

# Clinical metadata
# YP04A excluded from analysis (134 cells after QC - too sparse)
PDAC_METADATA = {
    "YP03A": {"patient": "YP03", "timepoint": "Pre", "response": "NR"},
    "YP03C": {"patient": "YP03", "timepoint": "Post", "response": "NR"},
    # "YP04A": excluded - 134 cells after QC
    "YP04C": {"patient": "YP04", "timepoint": "Post", "response": "NR"},
    "YP12A": {"patient": "YP12", "timepoint": "Pre", "response": "R"},
    "YP12C": {"patient": "YP12", "timepoint": "Post", "response": "R"},
    "YP15A": {"patient": "YP15", "timepoint": "Pre", "response": "R"},
    "YP15C": {"patient": "YP15", "timepoint": "Post", "response": "R"},
}

# =============================================================================
# Data Loading
# =============================================================================

def load_annotated_samples(
    input_dir: Optional[Path] = None,
    samples: Optional[List[str]] = None
) -> Dict[str, ad.AnnData]:
    """Load annotated PDAC samples."""

    if input_dir is None:
        input_dir = ADATA_DIR / "annotated"

    if samples is None:
        samples = list(PDAC_METADATA.keys())

    adatas = {}
    for sample in samples:
        h5ad_path = input_dir / f"{sample}_annotated.h5ad"
        if h5ad_path.exists():
            adata = sc.read_h5ad(h5ad_path)
            # Ensure metadata
            if sample in PDAC_METADATA:
                for k, v in PDAC_METADATA[sample].items():
                    adata.obs[k] = v
            adatas[sample] = adata
            print(f"  Loaded {sample}: {adata.n_obs} cells")
        else:
            print(f"  Warning: {sample} not found at {h5ad_path}")

    return adatas


def merge_samples(
    adatas: Dict[str, ad.AnnData],
    samples: Optional[List[str]] = None
) -> ad.AnnData:
    """Merge multiple samples into one AnnData."""

    if samples:
        adatas = {k: v for k, v in adatas.items() if k in samples}

    if len(adatas) == 0:
        return None

    # Find common genes
    common_genes = set(list(adatas.values())[0].var_names)
    for adata in adatas.values():
        common_genes &= set(adata.var_names)
    common_genes = list(common_genes)

    print(f"  Merging {len(adatas)} samples with {len(common_genes)} common genes")

    # Subset and concatenate
    subsets = []
    for name, adata in adatas.items():
        subset = adata[:, common_genes].copy()
        subset.obs['sample'] = name
        subsets.append(subset)

    merged = ad.concat(subsets, join='inner')
    return merged


# =============================================================================
# Cell Type Proportion Analysis
# =============================================================================

def compare_cell_type_proportions(
    adatas: Dict[str, ad.AnnData],
    group1_samples: List[str],
    group2_samples: List[str],
    group1_name: str = 'Group1',
    group2_name: str = 'Group2'
) -> pd.DataFrame:
    """
    Compare cell type proportions between two groups.

    Returns DataFrame with proportions and statistical tests.
    """
    print(f"\n  Comparing {group1_name} vs {group2_name}")

    # Get proportions for each sample
    all_props = []

    for sample in group1_samples + group2_samples:
        if sample not in adatas:
            continue

        adata = adatas[sample]
        if 'cell_type' not in adata.obs:
            continue

        props = adata.obs['cell_type'].value_counts(normalize=True)
        prop_dict = props.to_dict()
        prop_dict['sample'] = sample
        prop_dict['group'] = group1_name if sample in group1_samples else group2_name
        all_props.append(prop_dict)

    if len(all_props) == 0:
        return pd.DataFrame()

    props_df = pd.DataFrame(all_props).fillna(0)

    # Get cell types
    cell_types = [c for c in props_df.columns if c not in ['sample', 'group']]

    # Statistical comparison (Welch's t-test for better power with small n)
    results = []
    for ct in cell_types:
        g1_vals = props_df[props_df['group'] == group1_name][ct].values
        g2_vals = props_df[props_df['group'] == group2_name][ct].values

        if len(g1_vals) > 1 and len(g2_vals) > 1:
            stat, pval = stats.ttest_ind(g1_vals, g2_vals, equal_var=False)  # Welch's t-test
        else:
            pval = 1.0

        results.append({
            'cell_type': ct,
            f'{group1_name}_mean': np.mean(g1_vals),
            f'{group2_name}_mean': np.mean(g2_vals),
            'log2FC': np.log2((np.mean(g2_vals) + 0.001) / (np.mean(g1_vals) + 0.001)),
            'pvalue': pval
        })

    result_df = pd.DataFrame(results).sort_values('pvalue')
    return result_df


def plot_cell_type_comparison(
    adatas: Dict[str, ad.AnnData],
    comparisons: Dict[str, Tuple[List[str], List[str]]]
) -> None:
    """Plot cell type proportions for all comparisons."""

    # Gather all proportions
    all_data = []
    for sample, adata in adatas.items():
        if 'cell_type' not in adata.obs:
            continue

        meta = PDAC_METADATA.get(sample, {})
        props = adata.obs['cell_type'].value_counts(normalize=True)

        for ct, prop in props.items():
            all_data.append({
                'sample': sample,
                'response': meta.get('response', 'Unknown'),
                'timepoint': meta.get('timepoint', 'Unknown'),
                'cell_type': ct,
                'proportion': prop
            })

    if len(all_data) == 0:
        print("  No data for cell type comparison plot")
        return

    df = pd.DataFrame(all_data)

    # Create plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Plot 1: By response group
    cell_types = df['cell_type'].unique()
    ct_order = df.groupby('cell_type')['proportion'].mean().sort_values(ascending=False).index

    ax = axes[0]
    df_pivot = df.pivot_table(index='cell_type', columns=['response', 'timepoint'],
                              values='proportion', aggfunc='mean').fillna(0)
    df_pivot = df_pivot.reindex(ct_order)
    df_pivot.plot(kind='barh', ax=ax, width=0.8)
    ax.set_xlabel('Proportion')
    ax.set_title('Cell Type Proportions by Response & Timepoint')
    ax.legend(title='Group', bbox_to_anchor=(1.02, 1), loc='upper left')

    # Plot 2: R vs NR boxplot for top cell types
    ax = axes[1]
    top_cts = ct_order[:6]  # Top 6 cell types
    df_top = df[df['cell_type'].isin(top_cts)]
    sns.boxplot(data=df_top, x='cell_type', y='proportion', hue='response', ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_title('Cell Type Proportions: R vs NR')
    ax.legend(title='Response')

    plt.tight_layout()
    plt.savefig(FIG_DIR / "cell_type_proportions.png", dpi=150, bbox_inches='tight')
    plt.close()


# =============================================================================
# Spatial Organization Analysis
# =============================================================================

def compare_neighborhood_enrichment(
    adatas: Dict[str, ad.AnnData],
    group1_samples: List[str],
    group2_samples: List[str],
    cluster_key: str = 'cell_type'
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compare neighborhood enrichment patterns between groups.

    Returns mean z-scores for each group.
    """
    print(f"\n  Comparing spatial organization...")

    def get_enrichment(samples):
        matrices = []
        for sample in samples:
            if sample not in adatas:
                continue
            adata = adatas[sample]

            # Check if enrichment already computed
            key = f'{cluster_key}_nhood_enrichment'
            if key not in adata.uns:
                # Try to compute
                try:
                    if 'spatial_connectivities' not in adata.obsp:
                        sq.gr.spatial_neighbors(adata, n_neighs=6)
                    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key, n_perms=100)
                except:
                    continue

            if key in adata.uns:
                matrices.append(adata.uns[key]['zscore'])

        if len(matrices) > 0:
            # Ensure same shape (may differ if cell types differ)
            return np.nanmean(matrices, axis=0)
        return None

    z1 = get_enrichment(group1_samples)
    z2 = get_enrichment(group2_samples)

    return z1, z2


def plot_spatial_comparison(
    adatas: Dict[str, ad.AnnData],
    cluster_key: str = 'cell_type'
) -> None:
    """Plot spatial organization comparison between R and NR."""

    # Get enrichment for each group
    groups = {
        'R_Pre': SAMPLE_GROUPS['R_Pre'],
        'R_Post': SAMPLE_GROUPS['R_Post'],
        'NR_Pre': SAMPLE_GROUPS['NR_Pre'],
        'NR_Post': SAMPLE_GROUPS['NR_Post'],
    }

    enrichments = {}
    for name, samples in groups.items():
        z, _ = compare_neighborhood_enrichment(adatas, samples, [], cluster_key)
        if z is not None:
            enrichments[name] = z

    if len(enrichments) < 2:
        print("  Not enough data for spatial comparison plot")
        return

    # Plot
    n_groups = len(enrichments)
    fig, axes = plt.subplots(1, n_groups, figsize=(5*n_groups, 4))
    if n_groups == 1:
        axes = [axes]

    for i, (name, z) in enumerate(enrichments.items()):
        ax = axes[i]
        sns.heatmap(z, cmap='RdBu_r', center=0, ax=ax,
                   xticklabels=False, yticklabels=False)
        ax.set_title(f'{name}\nNeighborhood Enrichment')

    plt.suptitle('Spatial Organization by Response Group', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "spatial_organization_comparison.png", dpi=150, bbox_inches='tight')
    plt.close()


# =============================================================================
# Differential Expression
# =============================================================================

def run_differential_expression(
    adatas: Dict[str, ad.AnnData],
    group1_samples: List[str],
    group2_samples: List[str],
    group1_name: str = 'Group1',
    group2_name: str = 'Group2',
    cell_type: Optional[str] = None
) -> pd.DataFrame:
    """
    Run differential expression between groups.

    Args:
        adatas: Dictionary of AnnData objects
        group1_samples, group2_samples: Sample lists for each group
        cell_type: Restrict to specific cell type (None = all cells)

    Returns:
        DataFrame with DE results
    """
    print(f"\n  Running DE: {group1_name} vs {group2_name}")
    if cell_type:
        print(f"    Restricted to: {cell_type}")

    # Merge samples
    merged = merge_samples(adatas, group1_samples + group2_samples)
    if merged is None:
        return pd.DataFrame()

    # Filter to cell type if specified
    if cell_type and 'cell_type' in merged.obs:
        merged = merged[merged.obs['cell_type'] == cell_type].copy()
        if merged.n_obs < 10:
            print(f"    Too few cells ({merged.n_obs}), skipping")
            return pd.DataFrame()

    # Add group labels
    merged.obs['de_group'] = merged.obs['sample'].apply(
        lambda x: group1_name if x in group1_samples else group2_name
    )

    # Run DE
    try:
        sc.tl.rank_genes_groups(
            merged,
            groupby='de_group',
            method='wilcoxon',
            reference=group1_name
        )

        # Extract results for group2 vs group1
        result = sc.get.rank_genes_groups_df(merged, group=group2_name)
        result['comparison'] = f'{group2_name}_vs_{group1_name}'
        if cell_type:
            result['cell_type'] = cell_type

        print(f"    Found {(result['pvals_adj'] < 0.05).sum()} significant genes (FDR < 0.05)")

        return result

    except Exception as e:
        print(f"    DE failed: {e}")
        return pd.DataFrame()


def plot_volcano(
    de_results: pd.DataFrame,
    title: str = "Differential Expression"
) -> None:
    """Plot volcano plot from DE results."""

    if len(de_results) == 0:
        return

    fig, ax = plt.subplots(figsize=(10, 8))

    # Add -log10 pvalue
    de_results['-log10_padj'] = -np.log10(de_results['pvals_adj'] + 1e-100)

    # Color by significance
    colors = np.where(
        (de_results['pvals_adj'] < 0.05) & (abs(de_results['logfoldchanges']) > 0.5),
        np.where(de_results['logfoldchanges'] > 0, 'red', 'blue'),
        'gray'
    )

    ax.scatter(
        de_results['logfoldchanges'],
        de_results['-log10_padj'],
        c=colors,
        alpha=0.5,
        s=5
    )

    # Add lines
    ax.axhline(-np.log10(0.05), ls='--', c='black', alpha=0.5)
    ax.axvline(-0.5, ls='--', c='black', alpha=0.5)
    ax.axvline(0.5, ls='--', c='black', alpha=0.5)

    # Label top genes
    top_genes = de_results.nlargest(10, '-log10_padj')
    for _, row in top_genes.iterrows():
        ax.annotate(row['names'], (row['logfoldchanges'], row['-log10_padj']),
                   fontsize=8, alpha=0.8)

    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-Log10 Adjusted P-value')
    ax.set_title(title)

    plt.tight_layout()
    fname = title.replace(' ', '_').lower() + '_volcano.png'
    plt.savefig(FIG_DIR / fname, dpi=150, bbox_inches='tight')
    plt.close()


# =============================================================================
# Main Pipeline
# =============================================================================

def run_comparative_analysis(
    adatas: Optional[Dict[str, ad.AnnData]] = None,
    input_dir: Optional[Path] = None
) -> Dict[str, pd.DataFrame]:
    """
    Run full comparative analysis pipeline.

    Returns dictionary of result DataFrames.
    """
    print("="*60)
    print("COMPARATIVE ANALYSIS: Responders vs Non-Responders")
    print("="*60)

    # Load data if not provided
    if adatas is None:
        print("\n1. Loading annotated samples...")
        adatas = load_annotated_samples(input_dir)

    if len(adatas) < 2:
        print("  ERROR: Not enough samples loaded for comparison")
        return {}

    results = {}

    # 2. Cell type proportion comparison
    print("\n2. Cell Type Proportion Analysis")
    print("-"*40)

    # R vs NR at Pre
    ct_pre = compare_cell_type_proportions(
        adatas, SAMPLE_GROUPS['R_Pre'], SAMPLE_GROUPS['NR_Pre'],
        'R_Pre', 'NR_Pre'
    )
    if len(ct_pre) > 0:
        results['celltype_R_vs_NR_Pre'] = ct_pre
        print("\n  Pre-treatment comparison (R vs NR):")
        print(ct_pre[['cell_type', 'R_Pre_mean', 'NR_Pre_mean', 'log2FC', 'pvalue']].head(10).to_string())

    # R vs NR at Post
    ct_post = compare_cell_type_proportions(
        adatas, SAMPLE_GROUPS['R_Post'], SAMPLE_GROUPS['NR_Post'],
        'R_Post', 'NR_Post'
    )
    if len(ct_post) > 0:
        results['celltype_R_vs_NR_Post'] = ct_post

    # Plot
    plot_cell_type_comparison(adatas, {
        'R_vs_NR_Pre': (SAMPLE_GROUPS['R_Pre'], SAMPLE_GROUPS['NR_Pre']),
        'R_vs_NR_Post': (SAMPLE_GROUPS['R_Post'], SAMPLE_GROUPS['NR_Post']),
    })

    # 3. Spatial organization comparison
    print("\n3. Spatial Organization Analysis")
    print("-"*40)
    plot_spatial_comparison(adatas, cluster_key='cell_type')

    # 4. Differential expression
    print("\n4. Differential Expression Analysis")
    print("-"*40)

    # R vs NR at baseline (key comparison)
    de_pre = run_differential_expression(
        adatas, SAMPLE_GROUPS['NR_Pre'], SAMPLE_GROUPS['R_Pre'],
        'NR_Pre', 'R_Pre'
    )
    if len(de_pre) > 0:
        results['de_R_vs_NR_Pre'] = de_pre
        plot_volcano(de_pre, "R vs NR Pre-treatment")

    # R vs NR post-treatment
    de_post = run_differential_expression(
        adatas, SAMPLE_GROUPS['NR_Post'], SAMPLE_GROUPS['R_Post'],
        'NR_Post', 'R_Post'
    )
    if len(de_post) > 0:
        results['de_R_vs_NR_Post'] = de_post
        plot_volcano(de_post, "R vs NR Post-treatment")

    # 5. Save results
    print("\n5. Saving Results")
    print("-"*40)

    for name, df in results.items():
        out_path = TABLE_DIR / f"{name}.csv"
        df.to_csv(out_path, index=False)
        print(f"  Saved: {out_path}")

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    if 'celltype_R_vs_NR_Pre' in results:
        sig_ct = results['celltype_R_vs_NR_Pre'][results['celltype_R_vs_NR_Pre']['pvalue'] < 0.1]
        print(f"\nCell types with p<0.1 at baseline:")
        for _, row in sig_ct.iterrows():
            direction = "higher in R" if row['log2FC'] < 0 else "higher in NR"
            print(f"  - {row['cell_type']}: {direction} (p={row['pvalue']:.3f})")

    if 'de_R_vs_NR_Pre' in results:
        de = results['de_R_vs_NR_Pre']
        up_in_r = de[(de['pvals_adj'] < 0.05) & (de['logfoldchanges'] > 0.5)]
        down_in_r = de[(de['pvals_adj'] < 0.05) & (de['logfoldchanges'] < -0.5)]
        print(f"\nDE genes at baseline (FDR < 0.05, |log2FC| > 0.5):")
        print(f"  - Up in Responders: {len(up_in_r)} genes")
        print(f"  - Down in Responders: {len(down_in_r)} genes")

        if len(up_in_r) > 0:
            print(f"  - Top up in R: {', '.join(up_in_r.head(5)['names'].tolist())}")

    return results


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Comparative analysis: R vs NR")
    parser.add_argument("--input-dir", type=str, help="Input directory with annotated h5ad files")
    args = parser.parse_args()

    print("="*60)
    print("SPATIAL BIOLOGY HACKATHON 2026 - COMPARATIVE ANALYSIS")
    print("="*60)

    input_dir = Path(args.input_dir) if args.input_dir else None
    results = run_comparative_analysis(input_dir=input_dir)

    print("\n" + "="*60)
    print("COMPARATIVE ANALYSIS COMPLETE")
    print("="*60)
