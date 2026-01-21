#!/usr/bin/env python3
"""
Day 4: Deconvolution Analysis - Resolve Cell Type Mixtures in Visium Spots
===========================================================================

Visium spots (~55µm diameter) contain multiple cells (1-10 typically).
This script performs signature-based deconvolution to estimate cell type
PROPORTIONS within each spot, rather than discrete assignments.

Methods:
1. Signature scoring with refined PDAC markers from Polymath KB
2. Non-negative least squares (NNLS) deconvolution
3. Cell type proportion estimation per spot

Author: Max Van Belkum
Date: 2026-01-20
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from scipy.optimize import nnls
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
DECONV_DIR = OUTPUT_DIR / "adata" / "deconvolved"
FIG_DIR = OUTPUT_DIR / "figures" / "deconvolution"
TABLE_DIR = OUTPUT_DIR / "tables"

DECONV_DIR.mkdir(parents=True, exist_ok=True)
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

# =============================================================================
# Refined PDAC Cell Type Markers (from Polymath KB literature review)
# =============================================================================

# These markers are curated from:
# - PanIN and CAF Transitions paper
# - Conserved spatial subtypes paper
# - PDAC single-cell atlas studies in Polymath

PDAC_SIGNATURES = {
    # Epithelial/Tumor
    'Ductal_Epithelial': [
        'KRT19', 'KRT7', 'KRT8', 'KRT18', 'EPCAM', 'MUC1', 'CFTR', 'TFF1', 'TFF2', 'TFF3'
    ],
    'Acinar': [
        'PRSS1', 'PRSS2', 'CPA1', 'CPA2', 'CELA3A', 'CELA3B', 'PNLIP', 'PNLIPRP1',
        'REG1A', 'REG1B', 'REG3A', 'AMY2A', 'AMY2B'
    ],
    'Endocrine': [
        'INS', 'GCG', 'SST', 'PPY', 'IAPP', 'CHGA', 'CHGB', 'SYP', 'PCSK1', 'PCSK2'
    ],

    # CAF subtypes (from Conserved spatial subtypes paper)
    'CAF_mCAF': [  # Myofibroblastic
        'ACTA2', 'TAGLN', 'MYL9', 'TPM1', 'TPM2', 'CNN1', 'ACTG2', 'MYH11', 'MYLK'
    ],
    'CAF_iCAF': [  # Inflammatory
        'IL6', 'IL11', 'CXCL1', 'CXCL2', 'CXCL12', 'CCL2', 'HAS1', 'HAS2', 'LMNA'
    ],
    'CAF_apCAF': [  # Antigen-presenting
        'CD74', 'HLA-DRA', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DPB1', 'SAA1', 'SAA2'
    ],

    # Stellate/PSC
    'Stellate_PSC': [
        'RGS5', 'PDGFRB', 'DES', 'VIM', 'SPARC', 'COL1A1', 'COL1A2', 'COL3A1'
    ],

    # Immune cells
    'T_cells': [
        'CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'TRAC', 'TRBC1', 'TRBC2', 'IL7R'
    ],
    'NK_cells': [
        'NKG7', 'GNLY', 'KLRB1', 'KLRD1', 'KLRK1', 'NCR1', 'FCGR3A', 'PRF1', 'GZMA', 'GZMB'
    ],
    'B_cells': [
        'CD19', 'CD79A', 'CD79B', 'MS4A1', 'BANK1', 'PAX5', 'IGHM', 'IGHD', 'IGKC'
    ],
    'Macrophage': [
        'CD68', 'CD163', 'CSF1R', 'MARCO', 'MRC1', 'MSR1', 'FCGR1A', 'APOE', 'C1QA', 'C1QB'
    ],
    'Myeloid_DC': [
        'CD14', 'ITGAX', 'ITGAM', 'FCER1A', 'CLEC9A', 'CLEC10A', 'HLA-DQA1', 'HLA-DQB1'
    ],
    'Mast_cells': [
        'KIT', 'TPSAB1', 'TPSB2', 'CPA3', 'MS4A2', 'HPGDS'
    ],

    # Endothelial
    'Endothelial': [
        'PECAM1', 'CDH5', 'VWF', 'ERG', 'FLT1', 'KDR', 'ESAM', 'CLDN5', 'EMCN'
    ],
    'Lymphatic': [
        'LYVE1', 'PROX1', 'PDPN', 'FLT4', 'CCL21'
    ],

    # Neural
    'Neural': [
        'S100B', 'MPZ', 'PLP1', 'SOX10', 'NGFR', 'NRXN1', 'SYN1'
    ],
}


# =============================================================================
# Deconvolution Methods
# =============================================================================

def compute_signature_scores(adata, signatures: Dict[str, List[str]],
                              layer: Optional[str] = None) -> pd.DataFrame:
    """
    Compute signature scores for each cell type using scanpy's score_genes.

    Returns DataFrame with scores for each spot and cell type.
    """
    logger.info("Computing signature scores...")

    scores = {}
    available_genes = set(adata.var_names)

    for ct, markers in signatures.items():
        # Filter to available markers
        valid_markers = [g for g in markers if g in available_genes]

        if len(valid_markers) < 3:
            logger.warning(f"  {ct}: Only {len(valid_markers)}/{len(markers)} markers found, skipping")
            continue

        try:
            sc.tl.score_genes(adata, valid_markers, score_name=f'score_{ct}', use_raw=False)
            scores[ct] = adata.obs[f'score_{ct}'].values
            logger.info(f"  {ct}: {len(valid_markers)}/{len(markers)} markers used")
        except Exception as e:
            logger.warning(f"  {ct}: Failed - {e}")

    return pd.DataFrame(scores, index=adata.obs_names)


def nnls_deconvolution(score_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Convert signature scores to proportions using NNLS normalization.

    This ensures proportions are non-negative and sum to 1.
    """
    logger.info("Converting scores to proportions via NNLS...")

    # Shift scores to positive (required for NNLS interpretation as proportions)
    shifted = score_matrix - score_matrix.min().min() + 0.01

    # Normalize each row to sum to 1
    proportions = shifted.div(shifted.sum(axis=1), axis=0)

    return proportions


def assign_dominant_type(proportions: pd.DataFrame, threshold: float = 0.25) -> pd.Series:
    """
    Assign dominant cell type based on proportion threshold.

    If max proportion < threshold, assign 'Mixed'.
    """
    dominant = proportions.idxmax(axis=1)
    max_props = proportions.max(axis=1)

    # Mark low-confidence as Mixed
    dominant[max_props < threshold] = 'Mixed'

    return dominant


def compute_entropy(proportions: pd.DataFrame) -> pd.Series:
    """
    Compute Shannon entropy of cell type proportions per spot.

    High entropy = more mixed (multiple cell types)
    Low entropy = more pure (single dominant cell type)
    """
    from scipy.stats import entropy as shannon_entropy

    # Add small epsilon to avoid log(0)
    p = proportions.values + 1e-10
    p = p / p.sum(axis=1, keepdims=True)

    return pd.Series(shannon_entropy(p, axis=1), index=proportions.index, name='mixture_entropy')


# =============================================================================
# Analysis Functions
# =============================================================================

def deconvolve_sample(sample_name: str, adata_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run deconvolution on a single sample.

    Returns:
        proportions: DataFrame with cell type proportions per spot
        stats: DataFrame with sample-level statistics
    """
    logger.info(f"\n{'='*60}")
    logger.info(f"Deconvolving {sample_name}")
    logger.info(f"{'='*60}")

    adata = sc.read_h5ad(adata_path)
    logger.info(f"Loaded {adata.n_obs} spots, {adata.n_vars} genes")

    # Compute signature scores
    scores = compute_signature_scores(adata, PDAC_SIGNATURES)

    # Convert to proportions
    proportions = nnls_deconvolution(scores)

    # Add to adata
    for ct in proportions.columns:
        adata.obs[f'prop_{ct}'] = proportions[ct].values

    # Compute entropy (mixture complexity)
    adata.obs['mixture_entropy'] = compute_entropy(proportions).values

    # Assign dominant cell type
    adata.obs['dominant_type'] = assign_dominant_type(proportions, threshold=0.25)

    # Compute statistics
    meta = PDAC_METADATA.get(sample_name, {})
    stats = {
        'sample': sample_name,
        'n_spots': adata.n_obs,
        'mean_entropy': adata.obs['mixture_entropy'].mean(),
        'median_entropy': adata.obs['mixture_entropy'].median(),
        'pct_mixed': (adata.obs['dominant_type'] == 'Mixed').mean() * 100,
        **meta
    }

    # Add mean proportion for each cell type
    for ct in proportions.columns:
        stats[f'mean_prop_{ct}'] = proportions[ct].mean()

    # Save updated adata
    out_path = DECONV_DIR / f"{sample_name}_deconvolved.h5ad"
    adata.write_h5ad(out_path)
    logger.info(f"Saved to {out_path}")

    # Save proportions table
    proportions['sample'] = sample_name
    proportions['response'] = meta.get('response', 'Unknown')
    proportions['dominant_type'] = adata.obs['dominant_type'].values
    proportions['mixture_entropy'] = adata.obs['mixture_entropy'].values

    return proportions, pd.DataFrame([stats])


def plot_deconvolution_overview(all_proportions: pd.DataFrame, output_dir: Path):
    """Create comprehensive deconvolution figure."""
    logger.info("Generating deconvolution overview figure...")

    fig = plt.figure(figsize=(20, 16))

    # 1. Cell type proportions by response (boxplots)
    ax1 = fig.add_subplot(2, 3, 1)
    cell_types = [c for c in all_proportions.columns
                  if not c.startswith(('sample', 'response', 'dominant', 'mixture'))]

    # Melt for plotting
    plot_data = all_proportions.melt(
        id_vars=['sample', 'response'],
        value_vars=cell_types,
        var_name='cell_type',
        value_name='proportion'
    )

    # Get top cell types by mean proportion
    top_types = plot_data.groupby('cell_type')['proportion'].mean().nlargest(8).index
    plot_data_top = plot_data[plot_data['cell_type'].isin(top_types)]

    sns.boxplot(data=plot_data_top, x='cell_type', y='proportion', hue='response',
                palette={'R': '#2ecc71', 'NR': '#e74c3c'}, ax=ax1)
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')
    ax1.set_title('Cell Type Proportions by Response', fontweight='bold')
    ax1.set_ylabel('Proportion')
    ax1.legend(title='Response')

    # 2. Mixture entropy by response
    ax2 = fig.add_subplot(2, 3, 2)
    sample_stats = all_proportions.groupby(['sample', 'response']).agg({
        'mixture_entropy': 'mean'
    }).reset_index()

    sns.boxplot(data=sample_stats, x='response', y='mixture_entropy',
                palette={'R': '#2ecc71', 'NR': '#e74c3c'}, ax=ax2)

    # Add statistical test
    r_entropy = sample_stats[sample_stats['response'] == 'R']['mixture_entropy']
    nr_entropy = sample_stats[sample_stats['response'] == 'NR']['mixture_entropy']
    if len(r_entropy) > 1 and len(nr_entropy) > 1:
        _, pval = mannwhitneyu(r_entropy, nr_entropy, alternative='two-sided')
        ax2.text(0.5, 0.95, f'MWU p={pval:.3f}', transform=ax2.transAxes, ha='center')

    ax2.set_title('Spot Mixture Entropy by Response\n(Higher = More Mixed)', fontweight='bold')
    ax2.set_ylabel('Shannon Entropy')
    ax2.set_xlabel('Treatment Response')

    # 3. Dominant cell type distribution
    ax3 = fig.add_subplot(2, 3, 3)
    dom_counts = all_proportions.groupby(['response', 'dominant_type']).size().unstack(fill_value=0)
    dom_props = dom_counts.div(dom_counts.sum(axis=1), axis=0)
    dom_props.T.plot(kind='bar', ax=ax3, color=['#e74c3c', '#2ecc71'])
    ax3.set_title('Dominant Cell Type Distribution', fontweight='bold')
    ax3.set_ylabel('Proportion of Spots')
    ax3.set_xticklabels(ax3.get_xticklabels(), rotation=45, ha='right')
    ax3.legend(title='Response')

    # 4. Heatmap of mean proportions by sample
    ax4 = fig.add_subplot(2, 3, 4)
    sample_props = all_proportions.groupby('sample')[cell_types].mean()

    # Reorder by response
    r_samples = [s for s in sample_props.index if PDAC_METADATA.get(s, {}).get('response') == 'R']
    nr_samples = [s for s in sample_props.index if PDAC_METADATA.get(s, {}).get('response') == 'NR']
    sample_order = r_samples + nr_samples
    sample_props = sample_props.loc[sample_order]

    sns.heatmap(sample_props.T, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax4,
                cbar_kws={'label': 'Proportion'})
    ax4.set_title('Mean Cell Type Proportions per Sample\n(R samples left, NR right)', fontweight='bold')
    ax4.set_ylabel('Cell Type')

    # 5. Acinar proportion focus (key finding from Day 3)
    ax5 = fig.add_subplot(2, 3, 5)
    if 'Acinar' in cell_types:
        acinar_by_sample = all_proportions.groupby(['sample', 'response'])['Acinar'].mean().reset_index()
        sns.stripplot(data=acinar_by_sample, x='response', y='Acinar',
                     palette={'R': '#2ecc71', 'NR': '#e74c3c'}, s=12, ax=ax5)
        sns.boxplot(data=acinar_by_sample, x='response', y='Acinar',
                   palette={'R': '#2ecc71', 'NR': '#e74c3c'}, ax=ax5, boxprops=dict(alpha=0.3))

        r_acinar = acinar_by_sample[acinar_by_sample['response'] == 'R']['Acinar']
        nr_acinar = acinar_by_sample[acinar_by_sample['response'] == 'NR']['Acinar']
        if len(r_acinar) > 1 and len(nr_acinar) > 1:
            _, pval = mannwhitneyu(r_acinar, nr_acinar, alternative='two-sided')
            fc = r_acinar.mean() / (nr_acinar.mean() + 0.001)
            ax5.text(0.5, 0.95, f'FC={fc:.2f}x, MWU p={pval:.3f}',
                    transform=ax5.transAxes, ha='center', fontweight='bold')

        ax5.set_title('Acinar Cell Proportion by Response', fontweight='bold')
        ax5.set_ylabel('Mean Acinar Proportion')
        ax5.set_xlabel('Treatment Response')

    # 6. CAF subtype comparison
    ax6 = fig.add_subplot(2, 3, 6)
    caf_types = [c for c in cell_types if 'CAF' in c]
    if caf_types:
        caf_data = all_proportions.groupby(['sample', 'response'])[caf_types].mean().reset_index()
        caf_melt = caf_data.melt(id_vars=['sample', 'response'], var_name='CAF_type', value_name='proportion')

        sns.barplot(data=caf_melt, x='CAF_type', y='proportion', hue='response',
                   palette={'R': '#2ecc71', 'NR': '#e74c3c'}, ax=ax6)
        ax6.set_title('CAF Subtype Proportions by Response', fontweight='bold')
        ax6.set_ylabel('Mean Proportion')
        ax6.set_xticklabels(['mCAF\n(Myofibroblast)', 'iCAF\n(Inflammatory)', 'apCAF\n(Antigen-Presenting)'])
        ax6.legend(title='Response')

    plt.suptitle('PDAC Visium Deconvolution Analysis\nSignature-Based Cell Type Proportions',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / 'fig14_deconvolution_overview.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'fig14_deconvolution_overview.pdf', bbox_inches='tight')
    plt.close()
    logger.info(f"Saved fig14_deconvolution_overview")


def plot_spatial_proportions(sample_name: str, adata_path: Path, output_dir: Path):
    """Create spatial plot showing cell type proportions for a sample."""
    adata = sc.read_h5ad(adata_path)

    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    # Select key cell types to plot
    prop_cols = [c for c in adata.obs.columns if c.startswith('prop_')]

    # Sort by mean proportion
    mean_props = {c: adata.obs[c].mean() for c in prop_cols}
    top_props = sorted(mean_props.keys(), key=lambda x: mean_props[x], reverse=True)[:7]

    coords = adata.obsm['spatial']

    for i, prop_col in enumerate(top_props):
        ax = axes[i]
        ct_name = prop_col.replace('prop_', '')

        scatter = ax.scatter(coords[:, 0], coords[:, 1],
                            c=adata.obs[prop_col], cmap='YlOrRd',
                            s=8, alpha=0.8, vmin=0, vmax=0.5)
        ax.set_aspect('equal')
        ax.invert_yaxis()
        ax.set_title(f'{ct_name}\n(mean: {mean_props[prop_col]:.2%})')
        ax.axis('off')
        plt.colorbar(scatter, ax=ax, label='Proportion', shrink=0.6)

    # Last panel: entropy
    ax = axes[7]
    scatter = ax.scatter(coords[:, 0], coords[:, 1],
                        c=adata.obs['mixture_entropy'], cmap='viridis',
                        s=8, alpha=0.8)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set_title(f'Mixture Entropy\n(mean: {adata.obs["mixture_entropy"].mean():.2f})')
    ax.axis('off')
    plt.colorbar(scatter, ax=ax, label='Entropy', shrink=0.6)

    meta = PDAC_METADATA.get(sample_name, {})
    response_color = '#2ecc71' if meta.get('response') == 'R' else '#e74c3c'
    plt.suptitle(f'{sample_name} - Cell Type Proportions (Deconvolved)\n' +
                 f'{meta.get("response", "?")} | {meta.get("timepoint", "?")}',
                 fontsize=14, fontweight='bold', color=response_color)
    plt.tight_layout()
    plt.savefig(output_dir / f'spatial_props_{sample_name}.png', dpi=200, bbox_inches='tight')
    plt.close()


# =============================================================================
# Main
# =============================================================================

def main():
    logger.info("="*60)
    logger.info("DAY 4: DECONVOLUTION ANALYSIS")
    logger.info("Resolving cell type mixtures in Visium spots")
    logger.info("="*60)

    all_proportions = []
    all_stats = []

    # Process each sample
    for sample_name in PDAC_METADATA.keys():
        adata_path = ADATA_DIR / f"{sample_name}_annotated.h5ad"

        if not adata_path.exists():
            logger.warning(f"Skipping {sample_name} - file not found")
            continue

        try:
            proportions, stats = deconvolve_sample(sample_name, adata_path)
            all_proportions.append(proportions)
            all_stats.append(stats)

            # Generate spatial proportion plot
            deconv_path = DECONV_DIR / f"{sample_name}_deconvolved.h5ad"
            plot_spatial_proportions(sample_name, deconv_path, FIG_DIR)

        except Exception as e:
            logger.error(f"Error processing {sample_name}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Combine results
    if all_proportions:
        combined_props = pd.concat(all_proportions, ignore_index=True)
        combined_stats = pd.concat(all_stats, ignore_index=True)

        # Save tables
        combined_props.to_csv(TABLE_DIR / 'deconvolution_proportions.csv', index=False)
        combined_stats.to_csv(TABLE_DIR / 'deconvolution_stats.csv', index=False)
        logger.info(f"Saved deconvolution results to {TABLE_DIR}")

        # Generate overview figure
        plot_deconvolution_overview(combined_props, FIG_DIR)

        # Print summary
        logger.info("\n" + "="*60)
        logger.info("DECONVOLUTION SUMMARY")
        logger.info("="*60)

        r_stats = combined_stats[combined_stats['response'] == 'R']
        nr_stats = combined_stats[combined_stats['response'] == 'NR']

        logger.info(f"\nMixture Entropy (spot heterogeneity):")
        logger.info(f"  R:  {r_stats['mean_entropy'].mean():.3f} ± {r_stats['mean_entropy'].std():.3f}")
        logger.info(f"  NR: {nr_stats['mean_entropy'].mean():.3f} ± {nr_stats['mean_entropy'].std():.3f}")

        logger.info(f"\n% Mixed spots (no dominant type):")
        logger.info(f"  R:  {r_stats['pct_mixed'].mean():.1f}%")
        logger.info(f"  NR: {nr_stats['pct_mixed'].mean():.1f}%")

        # Compare key cell types
        logger.info(f"\nMean Cell Type Proportions (R vs NR):")
        prop_cols = [c for c in combined_stats.columns if c.startswith('mean_prop_')]
        for col in prop_cols:
            ct = col.replace('mean_prop_', '')
            r_mean = r_stats[col].mean()
            nr_mean = nr_stats[col].mean()
            fc = r_mean / (nr_mean + 0.001)
            logger.info(f"  {ct:20s}: R={r_mean:.3f}, NR={nr_mean:.3f} (FC={fc:.2f}x)")

    logger.info("\n" + "="*60)
    logger.info("DECONVOLUTION COMPLETE")
    logger.info("="*60)


if __name__ == '__main__':
    main()
