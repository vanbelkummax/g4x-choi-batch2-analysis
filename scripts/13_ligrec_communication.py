#!/usr/bin/env python3
"""
Day 4: Ligand-Receptor Communication Analysis
==============================================

Implements squidpy ligrec analysis for cell-cell communication.
Compares L-R interaction patterns between Responders vs Non-Responders.

Author: Max Van Belkum
Date: 2026-01-20
"""

import matplotlib
matplotlib.use('Agg')

import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from scipy.stats import mannwhitneyu
import warnings
import logging

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

# =============================================================================
# Configuration
# =============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
ADATA_DIR = OUTPUT_DIR / "adata" / "annotated"
FIG_DIR = OUTPUT_DIR / "figures" / "ligrec"
TABLE_DIR = OUTPUT_DIR / "tables"

FIG_DIR.mkdir(parents=True, exist_ok=True)

# Sample metadata
PDAC_METADATA = {
    "YP03A": {"patient": "YP03", "timepoint": "Pre", "response": "NR"},
    "YP03C": {"patient": "YP03", "timepoint": "Post", "response": "NR"},
    "YP04C": {"patient": "YP04", "timepoint": "Post", "response": "NR"},
    "YP12A": {"patient": "YP12", "timepoint": "Pre", "response": "R"},
    "YP12C": {"patient": "YP12", "timepoint": "Post", "response": "R"},
    "YP15A": {"patient": "YP15", "timepoint": "Pre", "response": "R"},
    "YP15C": {"patient": "YP15", "timepoint": "Post", "response": "R"},
}

RESPONSE_COLORS = {'R': '#2ecc71', 'NR': '#e74c3c'}


def run_ligrec_analysis(adata, sample_name: str, n_perms: int = 100) -> Optional[Dict]:
    """
    Run ligand-receptor analysis using squidpy.

    Uses CellPhoneDB database via squidpy for permutation-based
    testing of L-R interactions.
    """
    logger.info(f"  Running L-R analysis for {sample_name} (n_perms={n_perms})...")

    if 'cell_type' not in adata.obs.columns:
        logger.warning("  No cell_type column, skipping L-R analysis")
        return None

    try:
        # Ensure we have enough cells per type
        cell_counts = adata.obs['cell_type'].value_counts()
        valid_types = cell_counts[cell_counts >= 10].index.tolist()

        if len(valid_types) < 2:
            logger.warning(f"  Only {len(valid_types)} cell types with >=10 cells, skipping")
            return None

        # Filter to valid cell types
        adata_filt = adata[adata.obs['cell_type'].isin(valid_types)].copy()

        # Run permutation-based L-R analysis
        sq.gr.ligrec(
            adata_filt,
            n_perms=n_perms,
            cluster_key='cell_type',
            copy=False,
            use_raw=False,
            transmitter_params={'categories': 'ligand'},
            receiver_params={'categories': 'receptor'},
            seed=42,
        )

        if 'cell_type_ligrec' in adata_filt.uns:
            pvals = adata_filt.uns['cell_type_ligrec']['pvalues']
            means = adata_filt.uns['cell_type_ligrec']['means']

            # Count significant interactions
            n_sig = (pvals < 0.05).sum().sum()
            n_total = pvals.notna().sum().sum()

            logger.info(f"  Found {n_sig}/{n_total} significant L-R interactions (p<0.05)")

            return {
                'pvalues': pvals,
                'means': means,
                'n_significant': n_sig,
                'n_total': n_total,
                'cell_types': valid_types,
                'adata': adata_filt,
            }
        else:
            logger.warning("  No ligrec results found")
            return None

    except Exception as e:
        logger.error(f"  L-R analysis failed: {e}")
        return None


def get_top_interactions(pvals: pd.DataFrame, means: pd.DataFrame,
                         n_top: int = 20) -> pd.DataFrame:
    """Get top significant L-R interactions."""
    # Flatten the multi-index dataframes
    sig_interactions = []

    for (source, target), pval_series in pvals.items():
        for lr_pair, pval in pval_series.items():
            if pd.notna(pval) and pval < 0.05:
                mean_val = means.loc[lr_pair, (source, target)]
                sig_interactions.append({
                    'source': source,
                    'target': target,
                    'lr_pair': lr_pair,
                    'pvalue': pval,
                    'mean': mean_val,
                })

    if not sig_interactions:
        return pd.DataFrame()

    df = pd.DataFrame(sig_interactions)
    df = df.sort_values('mean', ascending=False).head(n_top)
    return df


def compare_ligrec_responses(results_dict: Dict) -> pd.DataFrame:
    """Compare L-R interaction counts between R and NR."""
    comparison = []

    for sample, res in results_dict.items():
        if res is None:
            continue
        meta = PDAC_METADATA[sample]
        comparison.append({
            'sample': sample,
            'response': meta['response'],
            'timepoint': meta['timepoint'],
            'n_significant': res['n_significant'],
            'n_total': res['n_total'],
            'pct_significant': 100 * res['n_significant'] / res['n_total'] if res['n_total'] > 0 else 0,
            'n_cell_types': len(res['cell_types']),
        })

    return pd.DataFrame(comparison)


def plot_ligrec_summary(comparison_df: pd.DataFrame, fig_dir: Path):
    """Plot L-R interaction summary comparing R vs NR."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    # 1. Total significant interactions
    ax = axes[0]
    r_data = comparison_df[comparison_df['response'] == 'R']['n_significant']
    nr_data = comparison_df[comparison_df['response'] == 'NR']['n_significant']

    x = [0, 1]
    colors = [RESPONSE_COLORS['R'], RESPONSE_COLORS['NR']]

    ax.bar(x, [r_data.mean(), nr_data.mean()], color=colors, alpha=0.7, width=0.6)
    ax.errorbar(x, [r_data.mean(), nr_data.mean()],
                yerr=[r_data.std(), nr_data.std()],
                fmt='none', color='black', capsize=5)

    # Add individual points
    for val in r_data:
        ax.scatter(0 + np.random.uniform(-0.1, 0.1), val, color='black', s=50, zorder=5)
    for val in nr_data:
        ax.scatter(1 + np.random.uniform(-0.1, 0.1), val, color='black', s=50, zorder=5)

    # Statistical test
    if len(r_data) > 1 and len(nr_data) > 1:
        stat, pval = mannwhitneyu(r_data, nr_data, alternative='two-sided')
        ax.set_title(f'Significant L-R Interactions\np={pval:.3f}', fontsize=12)
    else:
        ax.set_title('Significant L-R Interactions', fontsize=12)

    ax.set_xticks(x)
    ax.set_xticklabels(['Responders', 'Non-Responders'])
    ax.set_ylabel('Count')

    # 2. Percent significant
    ax = axes[1]
    r_pct = comparison_df[comparison_df['response'] == 'R']['pct_significant']
    nr_pct = comparison_df[comparison_df['response'] == 'NR']['pct_significant']

    ax.bar(x, [r_pct.mean(), nr_pct.mean()], color=colors, alpha=0.7, width=0.6)
    ax.errorbar(x, [r_pct.mean(), nr_pct.mean()],
                yerr=[r_pct.std(), nr_pct.std()],
                fmt='none', color='black', capsize=5)

    for val in r_pct:
        ax.scatter(0 + np.random.uniform(-0.1, 0.1), val, color='black', s=50, zorder=5)
    for val in nr_pct:
        ax.scatter(1 + np.random.uniform(-0.1, 0.1), val, color='black', s=50, zorder=5)

    if len(r_pct) > 1 and len(nr_pct) > 1:
        stat, pval = mannwhitneyu(r_pct, nr_pct, alternative='two-sided')
        ax.set_title(f'% Significant Interactions\np={pval:.3f}', fontsize=12)
    else:
        ax.set_title('% Significant Interactions', fontsize=12)

    ax.set_xticks(x)
    ax.set_xticklabels(['Responders', 'Non-Responders'])
    ax.set_ylabel('Percentage')

    # 3. Per-sample breakdown
    ax = axes[2]
    sample_order = comparison_df.sort_values(['response', 'timepoint'])['sample'].tolist()
    colors_per_sample = [RESPONSE_COLORS[PDAC_METADATA[s]['response']] for s in sample_order]

    y_pos = range(len(sample_order))
    vals = [comparison_df[comparison_df['sample'] == s]['n_significant'].values[0] for s in sample_order]

    ax.barh(y_pos, vals, color=colors_per_sample, alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sample_order)
    ax.set_xlabel('Significant L-R Interactions')
    ax.set_title('Per-Sample Breakdown', fontsize=12)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=RESPONSE_COLORS['R'], label='Responder'),
                       Patch(facecolor=RESPONSE_COLORS['NR'], label='Non-Responder')]
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()

    for ext in ['png', 'pdf']:
        fig.savefig(fig_dir / f'fig18_ligrec_summary.{ext}', dpi=150, bbox_inches='tight')
    plt.close()
    logger.info("Saved fig18_ligrec_summary")


def plot_top_interactions_heatmap(all_top_interactions: Dict, fig_dir: Path):
    """Plot heatmap of top L-R interactions across samples."""
    # Collect all interactions (convert tuples to strings)
    all_pairs = set()
    for sample, df in all_top_interactions.items():
        if df is not None and not df.empty:
            for lr in df['lr_pair'].tolist():
                # Convert tuple to string if needed
                if isinstance(lr, tuple):
                    lr = f"{lr[0]}->{lr[1]}"
                all_pairs.add(lr)

    if not all_pairs:
        logger.warning("No significant interactions to plot")
        return

    # Build matrix
    samples = list(PDAC_METADATA.keys())
    matrix = pd.DataFrame(0.0, index=list(all_pairs), columns=samples)

    for sample, df in all_top_interactions.items():
        if df is not None and not df.empty:
            for _, row in df.iterrows():
                lr = row['lr_pair']
                if isinstance(lr, tuple):
                    lr = f"{lr[0]}->{lr[1]}"
                matrix.loc[lr, sample] = row['mean']

    # Filter to interactions present in at least 2 samples
    matrix = matrix[matrix.sum(axis=1) > 0]
    if len(matrix) > 30:
        # Take top 30 by total signal
        matrix = matrix.loc[matrix.sum(axis=1).nlargest(30).index]

    if matrix.empty:
        logger.warning("No interactions to plot after filtering")
        return

    # Create heatmap
    fig, ax = plt.subplots(figsize=(10, max(8, len(matrix) * 0.3)))

    # Sort samples by response
    sample_order = ['YP12A', 'YP12C', 'YP15A', 'YP15C', 'YP03A', 'YP03C', 'YP04C']
    sample_order = [s for s in sample_order if s in matrix.columns]
    matrix = matrix[sample_order]

    # Cluster rows
    from scipy.cluster.hierarchy import linkage, leaves_list
    if len(matrix) > 2:
        Z = linkage(matrix.values, method='average')
        order = leaves_list(Z)
        matrix = matrix.iloc[order]

    sns.heatmap(matrix, cmap='YlOrRd', ax=ax, cbar_kws={'label': 'Mean Expression'},
                xticklabels=True, yticklabels=True)

    ax.set_xlabel('Sample')
    ax.set_ylabel('Ligand-Receptor Pair')
    ax.set_title('Top L-R Interactions Across PDAC Samples', fontsize=12)

    # Color x-axis labels by response
    for i, label in enumerate(ax.get_xticklabels()):
        sample = label.get_text()
        color = RESPONSE_COLORS[PDAC_METADATA[sample]['response']]
        label.set_color(color)
        label.set_fontweight('bold')

    plt.tight_layout()

    for ext in ['png', 'pdf']:
        fig.savefig(fig_dir / f'fig19_ligrec_heatmap.{ext}', dpi=150, bbox_inches='tight')
    plt.close()
    logger.info("Saved fig19_ligrec_heatmap")


def plot_ligrec_dotplot(results_dict: Dict, fig_dir: Path, sample: str = None):
    """Plot squidpy's ligrec dotplot for a representative sample."""
    # Find a sample with results
    if sample is None:
        for s, res in results_dict.items():
            if res is not None and res['n_significant'] > 0:
                sample = s
                break

    if sample is None or results_dict.get(sample) is None:
        logger.warning("No sample with L-R results for dotplot")
        return

    res = results_dict[sample]
    adata = res['adata']

    try:
        fig, ax = plt.subplots(figsize=(12, 10))
        sq.pl.ligrec(
            adata,
            cluster_key='cell_type',
            source_groups=res['cell_types'][:5],  # Limit for readability
            target_groups=res['cell_types'][:5],
            means_range=(0.1, np.inf),
            pvalue_threshold=0.05,
            ax=ax,
        )
        ax.set_title(f'L-R Interactions: {sample} ({PDAC_METADATA[sample]["response"]})', fontsize=12)

        plt.tight_layout()
        for ext in ['png', 'pdf']:
            fig.savefig(fig_dir / f'fig20_ligrec_dotplot_{sample}.{ext}', dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved fig20_ligrec_dotplot_{sample}")
    except Exception as e:
        logger.warning(f"Could not generate dotplot: {e}")


def main():
    """Main analysis function."""
    logger.info("=" * 60)
    logger.info("LIGAND-RECEPTOR COMMUNICATION ANALYSIS")
    logger.info("=" * 60)

    # Load samples and run L-R analysis
    results_dict = {}
    top_interactions = {}

    for sample_name in PDAC_METADATA.keys():
        adata_path = ADATA_DIR / f"{sample_name}_annotated.h5ad"

        if not adata_path.exists():
            logger.warning(f"Missing: {sample_name}")
            continue

        logger.info(f"\nProcessing {sample_name}...")
        adata = sc.read_h5ad(adata_path)

        # Run L-R analysis
        result = run_ligrec_analysis(adata, sample_name, n_perms=100)
        results_dict[sample_name] = result

        if result is not None:
            # Get top interactions
            top_df = get_top_interactions(result['pvalues'], result['means'], n_top=30)
            top_interactions[sample_name] = top_df

            # Save per-sample results
            if not top_df.empty:
                top_df.to_csv(TABLE_DIR / f'ligrec_{sample_name}_top.csv', index=False)

    # Compare R vs NR
    logger.info("\n" + "=" * 60)
    logger.info("COMPARING RESPONDERS vs NON-RESPONDERS")
    logger.info("=" * 60)

    comparison_df = compare_ligrec_responses(results_dict)
    if not comparison_df.empty:
        comparison_df.to_csv(TABLE_DIR / 'ligrec_comparison.csv', index=False)

        # Print summary
        r_sig = comparison_df[comparison_df['response'] == 'R']['n_significant'].mean()
        nr_sig = comparison_df[comparison_df['response'] == 'NR']['n_significant'].mean()
        logger.info(f"  Mean significant interactions: R={r_sig:.1f}, NR={nr_sig:.1f}")

        # Generate figures
        plot_ligrec_summary(comparison_df, FIG_DIR)
        plot_top_interactions_heatmap(top_interactions, FIG_DIR)

        # Dotplot for one R and one NR sample
        for sample in ['YP12A', 'YP03A']:
            if sample in results_dict and results_dict[sample] is not None:
                plot_ligrec_dotplot(results_dict, FIG_DIR, sample)

    logger.info("\n" + "=" * 60)
    logger.info("L-R ANALYSIS COMPLETE")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
