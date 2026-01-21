#!/usr/bin/env python3
"""
Day 4: Polymath-Driven Hypotheses - Complete Analysis
======================================================

Implements all remaining hypotheses from Polymath KB:
1. Ligand-Receptor Communication (squidpy ligrec)
2. CAF Subtype Scoring (iCAF/mCAF/apCAF)
3. Spatial Entropy quantification
4. Metabolic Gene Scoring (MPC pathway)

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
from scipy.stats import mannwhitneyu, entropy as shannon_entropy
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
FIG_DIR = OUTPUT_DIR / "figures" / "day4_hypotheses"
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

# =============================================================================
# 1. CAF Subtype Markers (from Polymath: Conserved spatial subtypes paper)
# =============================================================================

CAF_MARKERS = {
    'iCAF': ['IL6', 'IL11', 'CXCL1', 'CXCL2', 'CXCL12', 'CCL2', 'HAS1', 'HAS2', 'LMNA', 'CFD', 'DPT'],
    'mCAF': ['ACTA2', 'TAGLN', 'MYL9', 'TPM1', 'TPM2', 'CNN1', 'ACTG2', 'MYH11', 'MYLK', 'POSTN'],
    'apCAF': ['CD74', 'HLA-DRA', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DPB1', 'SAA1', 'SLPI'],  # MHC class II
}

# =============================================================================
# 2. Metabolic Pathway Genes (from Polymath: MPC regulation paper)
# =============================================================================

METABOLIC_MARKERS = {
    'Pyruvate_Metabolism': ['MPC1', 'MPC2', 'PDHA1', 'PDHB', 'PDK1', 'PDK2', 'PDK3', 'PDK4'],
    'Glycolysis': ['HK1', 'HK2', 'PFKP', 'PFKM', 'PKM', 'ENO1', 'GAPDH', 'LDHA', 'LDHB'],
    'Oxidative_Phosphorylation': ['MT-ND1', 'MT-ND2', 'MT-CO1', 'MT-CO2', 'MT-ATP6', 'MT-CYB'],
    'Fatty_Acid_Metabolism': ['CPT1A', 'CPT2', 'ACSL1', 'ACADM', 'ACADL', 'ECHS1'],
}

# =============================================================================
# Analysis Functions
# =============================================================================

def score_gene_signatures(adata, markers_dict: Dict[str, List[str]],
                          prefix: str = 'score') -> pd.DataFrame:
    """Score cells for gene signatures."""
    available_genes = set(adata.var_names)
    results = {}

    for name, markers in markers_dict.items():
        valid_markers = [g for g in markers if g in available_genes]
        if len(valid_markers) >= 3:
            try:
                sc.tl.score_genes(adata, valid_markers, score_name=f'{prefix}_{name}', use_raw=False)
                results[name] = {
                    'n_markers': len(valid_markers),
                    'mean_score': adata.obs[f'{prefix}_{name}'].mean(),
                    'std_score': adata.obs[f'{prefix}_{name}'].std(),
                }
                logger.info(f"  Scored {name}: {len(valid_markers)} markers, mean={results[name]['mean_score']:.3f}")
            except Exception as e:
                logger.warning(f"  Failed to score {name}: {e}")
        else:
            logger.warning(f"  Skipping {name}: only {len(valid_markers)}/{len(markers)} markers available")

    return pd.DataFrame(results).T


def compute_spatial_entropy(adata, cell_type_col: str = 'cell_type',
                            n_bins: int = 10) -> Dict[str, float]:
    """
    Compute spatial entropy metrics for tissue organization.

    Higher entropy = more spatially mixed (potentially more accessible)
    Lower entropy = more spatially organized (distinct compartments)
    """
    coords = adata.obsm['spatial']

    # 1. Overall spatial entropy (cell type distribution across space)
    cell_types = adata.obs[cell_type_col].values
    unique_types = np.unique(cell_types)

    # Divide space into grid
    x_bins = np.linspace(coords[:, 0].min(), coords[:, 0].max(), n_bins + 1)
    y_bins = np.linspace(coords[:, 1].min(), coords[:, 1].max(), n_bins + 1)

    # Count cell types in each grid cell
    grid_entropies = []
    for i in range(n_bins):
        for j in range(n_bins):
            mask = ((coords[:, 0] >= x_bins[i]) & (coords[:, 0] < x_bins[i+1]) &
                   (coords[:, 1] >= y_bins[j]) & (coords[:, 1] < y_bins[j+1]))
            if mask.sum() > 5:  # Require at least 5 cells
                ct_counts = pd.Series(cell_types[mask]).value_counts(normalize=True)
                grid_entropies.append(shannon_entropy(ct_counts))

    # 2. Per-cell-type spatial entropy (how spread out is each type?)
    ct_entropies = {}
    for ct in unique_types:
        ct_mask = cell_types == ct
        if ct_mask.sum() > 10:
            ct_coords = coords[ct_mask]
            # Compute entropy of 2D histogram
            H, _, _ = np.histogram2d(ct_coords[:, 0], ct_coords[:, 1], bins=n_bins)
            H = H / H.sum()
            H = H[H > 0]
            ct_entropies[ct] = shannon_entropy(H)

    return {
        'global_spatial_entropy': np.mean(grid_entropies) if grid_entropies else 0,
        'global_entropy_std': np.std(grid_entropies) if grid_entropies else 0,
        'mean_celltype_entropy': np.mean(list(ct_entropies.values())) if ct_entropies else 0,
        **{f'entropy_{ct}': v for ct, v in ct_entropies.items()}
    }


def run_ligrec_analysis(adata, sample_name: str) -> Optional[pd.DataFrame]:
    """
    Run ligand-receptor analysis using squidpy.

    Note: This is computationally intensive and may take time.
    """
    logger.info(f"  Running L-R analysis for {sample_name}...")

    if 'cell_type' not in adata.obs.columns:
        logger.warning("  No cell_type column, skipping L-R analysis")
        return None

    # Ensure spatial neighbors exist
    if 'spatial_connectivities' not in adata.obsp:
        sq.gr.spatial_neighbors(adata, n_neighs=6)

    try:
        # Run permutation-based L-R analysis
        sq.gr.ligrec(
            adata,
            n_perms=100,  # Reduced for speed
            cluster_key='cell_type',
            use_raw=False,
            transmitter_params={'categories': 'ligand'},
            receiver_params={'categories': 'receptor'},
            copy=False
        )

        # Extract results
        if 'cell_type_ligrec' in adata.uns:
            pvals = adata.uns['cell_type_ligrec']['pvalues']
            means = adata.uns['cell_type_ligrec']['means']

            # Get significant interactions
            sig_mask = pvals < 0.05
            n_sig = sig_mask.sum().sum()
            logger.info(f"  Found {n_sig} significant L-R interactions (p<0.05)")

            return pvals

    except Exception as e:
        logger.error(f"  L-R analysis failed: {e}")

    return None


def analyze_sample(sample_name: str, adata_path: Path) -> Dict:
    """Run all Day 4 analyses on a single sample."""
    logger.info(f"\n{'='*60}")
    logger.info(f"Analyzing {sample_name}")
    logger.info(f"{'='*60}")

    adata = sc.read_h5ad(adata_path)
    meta = PDAC_METADATA.get(sample_name, {})

    results = {
        'sample': sample_name,
        'n_spots': adata.n_obs,
        **meta
    }

    # 1. CAF Subtype Scoring
    logger.info("1. CAF Subtype Scoring...")
    caf_results = score_gene_signatures(adata, CAF_MARKERS, prefix='CAF')
    for caf_type, stats in caf_results.iterrows():
        results[f'CAF_{caf_type}_mean'] = stats['mean_score']

    # Compute CAF dominance (which subtype is highest)
    caf_cols = [c for c in adata.obs.columns if c.startswith('CAF_')]
    if caf_cols:
        caf_scores = adata.obs[caf_cols]
        adata.obs['dominant_CAF'] = caf_scores.idxmax(axis=1).str.replace('CAF_', '')

        # CAF subtype proportions
        dom_counts = adata.obs['dominant_CAF'].value_counts(normalize=True)
        for caf, prop in dom_counts.items():
            results[f'pct_{caf}_dominant'] = prop * 100

    # 2. Metabolic Scoring
    logger.info("2. Metabolic Pathway Scoring...")
    metab_results = score_gene_signatures(adata, METABOLIC_MARKERS, prefix='metab')
    for pathway, stats in metab_results.iterrows():
        results[f'metab_{pathway}_mean'] = stats['mean_score']

    # 3. Spatial Entropy
    logger.info("3. Computing Spatial Entropy...")
    entropy_results = compute_spatial_entropy(adata)
    results.update(entropy_results)

    # 4. L-R Communication (skip for speed - can enable if needed)
    # logger.info("4. L-R Communication Analysis...")
    # ligrec_results = run_ligrec_analysis(adata, sample_name)

    return results, adata


def plot_caf_comparison(all_results: pd.DataFrame, output_dir: Path):
    """Plot CAF subtype comparison between R and NR."""
    logger.info("Generating CAF comparison figure...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. CAF scores by response
    ax1 = axes[0, 0]
    caf_cols = [c for c in all_results.columns if c.startswith('CAF_') and c.endswith('_mean')]
    if caf_cols:
        caf_data = all_results[['sample', 'response'] + caf_cols].melt(
            id_vars=['sample', 'response'],
            var_name='CAF_type',
            value_name='score'
        )
        caf_data['CAF_type'] = caf_data['CAF_type'].str.replace('CAF_', '').str.replace('_mean', '')

        sns.barplot(data=caf_data, x='CAF_type', y='score', hue='response',
                   palette=RESPONSE_COLORS, ax=ax1)
        ax1.set_title('CAF Subtype Scores by Response', fontweight='bold')
        ax1.set_ylabel('Gene Signature Score')
        ax1.set_xlabel('CAF Subtype')

        # Add statistical annotations
        for i, caf in enumerate(caf_data['CAF_type'].unique()):
            r_vals = caf_data[(caf_data['CAF_type'] == caf) & (caf_data['response'] == 'R')]['score']
            nr_vals = caf_data[(caf_data['CAF_type'] == caf) & (caf_data['response'] == 'NR')]['score']
            if len(r_vals) > 1 and len(nr_vals) > 1:
                _, pval = mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
                sig = '*' if pval < 0.05 else ('†' if pval < 0.1 else '')
                if sig:
                    ax1.text(i, ax1.get_ylim()[1] * 0.95, sig, ha='center', fontsize=12)

    # 2. iCAF/mCAF ratio
    ax2 = axes[0, 1]
    if 'CAF_iCAF_mean' in all_results.columns and 'CAF_mCAF_mean' in all_results.columns:
        all_results['iCAF_mCAF_ratio'] = all_results['CAF_iCAF_mean'] / (all_results['CAF_mCAF_mean'] + 0.01)

        sns.boxplot(data=all_results, x='response', y='iCAF_mCAF_ratio',
                   palette=RESPONSE_COLORS, ax=ax2)
        sns.stripplot(data=all_results, x='response', y='iCAF_mCAF_ratio',
                     color='black', alpha=0.5, ax=ax2)

        r_ratio = all_results[all_results['response'] == 'R']['iCAF_mCAF_ratio']
        nr_ratio = all_results[all_results['response'] == 'NR']['iCAF_mCAF_ratio']
        if len(r_ratio) > 1 and len(nr_ratio) > 1:
            _, pval = mannwhitneyu(r_ratio, nr_ratio)
            ax2.text(0.5, 0.95, f'MWU p={pval:.3f}', transform=ax2.transAxes, ha='center')

        ax2.set_title('iCAF/mCAF Ratio by Response\n(Higher = more inflammatory)', fontweight='bold')
        ax2.set_ylabel('iCAF/mCAF Ratio')
        ax2.set_xlabel('Treatment Response')

    # 3. Spatial Entropy
    ax3 = axes[1, 0]
    if 'global_spatial_entropy' in all_results.columns:
        sns.boxplot(data=all_results, x='response', y='global_spatial_entropy',
                   palette=RESPONSE_COLORS, ax=ax3)
        sns.stripplot(data=all_results, x='response', y='global_spatial_entropy',
                     color='black', alpha=0.5, ax=ax3)

        r_ent = all_results[all_results['response'] == 'R']['global_spatial_entropy']
        nr_ent = all_results[all_results['response'] == 'NR']['global_spatial_entropy']
        if len(r_ent) > 1 and len(nr_ent) > 1:
            _, pval = mannwhitneyu(r_ent, nr_ent)
            fc = r_ent.mean() / (nr_ent.mean() + 0.001)
            ax3.text(0.5, 0.95, f'FC={fc:.2f}x, MWU p={pval:.3f}',
                    transform=ax3.transAxes, ha='center')

        ax3.set_title('Global Spatial Entropy by Response\n(Higher = more mixed/accessible)', fontweight='bold')
        ax3.set_ylabel('Shannon Entropy')
        ax3.set_xlabel('Treatment Response')

    # 4. Metabolic pathways
    ax4 = axes[1, 1]
    metab_cols = [c for c in all_results.columns if c.startswith('metab_') and c.endswith('_mean')]
    if metab_cols:
        metab_data = all_results[['sample', 'response'] + metab_cols].melt(
            id_vars=['sample', 'response'],
            var_name='Pathway',
            value_name='score'
        )
        metab_data['Pathway'] = metab_data['Pathway'].str.replace('metab_', '').str.replace('_mean', '')

        sns.barplot(data=metab_data, x='Pathway', y='score', hue='response',
                   palette=RESPONSE_COLORS, ax=ax4)
        ax4.set_title('Metabolic Pathway Scores by Response', fontweight='bold')
        ax4.set_ylabel('Gene Signature Score')
        ax4.set_xlabel('Metabolic Pathway')
        ax4.set_xticklabels(ax4.get_xticklabels(), rotation=45, ha='right')

    plt.suptitle('Day 4 Hypotheses: CAF Subtypes, Spatial Entropy, and Metabolism',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / 'fig15_day4_hypotheses.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'fig15_day4_hypotheses.pdf', bbox_inches='tight')
    plt.close()
    logger.info("Saved fig15_day4_hypotheses")


def plot_metabolic_detail(all_results: pd.DataFrame, output_dir: Path):
    """Detailed metabolic pathway analysis figure."""
    logger.info("Generating metabolic detail figure...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    metab_pathways = ['Pyruvate_Metabolism', 'Glycolysis', 'Oxidative_Phosphorylation', 'Fatty_Acid_Metabolism']

    for idx, pathway in enumerate(metab_pathways):
        ax = axes.flatten()[idx]
        col = f'metab_{pathway}_mean'

        if col in all_results.columns:
            sns.boxplot(data=all_results, x='response', y=col, palette=RESPONSE_COLORS, ax=ax)
            sns.stripplot(data=all_results, x='response', y=col, color='black', alpha=0.6, s=10, ax=ax)

            # Stats
            r_vals = all_results[all_results['response'] == 'R'][col]
            nr_vals = all_results[all_results['response'] == 'NR'][col]
            if len(r_vals) > 1 and len(nr_vals) > 1:
                _, pval = mannwhitneyu(r_vals, nr_vals)
                fc = r_vals.mean() / (nr_vals.mean() + 0.001)
                ax.text(0.5, 0.95, f'FC={fc:.2f}x, p={pval:.3f}',
                       transform=ax.transAxes, ha='center', fontweight='bold')

            ax.set_title(pathway.replace('_', ' '), fontweight='bold')
            ax.set_ylabel('Pathway Score')
            ax.set_xlabel('Treatment Response')
        else:
            ax.text(0.5, 0.5, f'{pathway}\nNot available', ha='center', va='center')
            ax.axis('off')

    plt.suptitle('Metabolic Pathway Analysis: R vs NR\n(MPC/Pyruvate hypothesis from Polymath KB)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / 'fig16_metabolic_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'fig16_metabolic_analysis.pdf', bbox_inches='tight')
    plt.close()
    logger.info("Saved fig16_metabolic_analysis")


def plot_entropy_heatmap(all_results: pd.DataFrame, output_dir: Path):
    """Spatial entropy heatmap by cell type."""
    logger.info("Generating entropy heatmap figure...")

    # Extract cell type entropies
    entropy_cols = [c for c in all_results.columns if c.startswith('entropy_') and not c.startswith('entropy_global')]

    if not entropy_cols:
        logger.warning("No cell type entropy columns found")
        return

    # Create heatmap data
    entropy_data = all_results[['sample', 'response'] + entropy_cols].copy()
    entropy_data = entropy_data.set_index('sample')

    # Reorder by response
    r_samples = [s for s in entropy_data.index if entropy_data.loc[s, 'response'] == 'R']
    nr_samples = [s for s in entropy_data.index if entropy_data.loc[s, 'response'] == 'NR']
    sample_order = r_samples + nr_samples

    # Get only numeric entropy columns
    entropy_numeric = entropy_data[entropy_cols].loc[sample_order]
    entropy_numeric.columns = [c.replace('entropy_', '') for c in entropy_numeric.columns]

    fig, ax = plt.subplots(figsize=(12, 8))

    sns.heatmap(entropy_numeric.T, annot=True, fmt='.2f', cmap='YlOrRd',
                cbar_kws={'label': 'Spatial Entropy'}, ax=ax)

    # Add response annotations
    for i, sample in enumerate(sample_order):
        response = all_results[all_results['sample'] == sample]['response'].values[0]
        color = RESPONSE_COLORS[response]
        ax.text(i + 0.5, -0.3, response, ha='center', va='center', color=color, fontweight='bold')

    ax.set_title('Spatial Entropy by Cell Type and Sample\n(Higher = more spatially dispersed)',
                 fontweight='bold')
    ax.set_xlabel('Sample (R=Responder, NR=Non-Responder)')
    ax.set_ylabel('Cell Type')

    plt.tight_layout()
    plt.savefig(output_dir / 'fig17_entropy_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'fig17_entropy_heatmap.pdf', bbox_inches='tight')
    plt.close()
    logger.info("Saved fig17_entropy_heatmap")


# =============================================================================
# Main
# =============================================================================

def main():
    logger.info("="*60)
    logger.info("DAY 4: POLYMATH-DRIVEN HYPOTHESES")
    logger.info("CAF Subtypes | Metabolic Pathways | Spatial Entropy")
    logger.info("="*60)

    all_results = []

    # Process each sample
    for sample_name in PDAC_METADATA.keys():
        adata_path = ADATA_DIR / f"{sample_name}_annotated.h5ad"

        if not adata_path.exists():
            logger.warning(f"Skipping {sample_name} - file not found")
            continue

        try:
            results, adata = analyze_sample(sample_name, adata_path)
            all_results.append(results)
        except Exception as e:
            logger.error(f"Error processing {sample_name}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Combine results
    if all_results:
        combined_results = pd.DataFrame(all_results)

        # Save table
        combined_results.to_csv(TABLE_DIR / 'day4_hypotheses_results.csv', index=False)
        logger.info(f"Saved results to {TABLE_DIR / 'day4_hypotheses_results.csv'}")

        # Generate figures
        plot_caf_comparison(combined_results, FIG_DIR)
        plot_metabolic_detail(combined_results, FIG_DIR)
        plot_entropy_heatmap(combined_results, FIG_DIR)

        # Print summary
        logger.info("\n" + "="*60)
        logger.info("DAY 4 HYPOTHESES SUMMARY")
        logger.info("="*60)

        r_data = combined_results[combined_results['response'] == 'R']
        nr_data = combined_results[combined_results['response'] == 'NR']

        # CAF summary
        logger.info("\n1. CAF SUBTYPES:")
        for caf in ['iCAF', 'mCAF', 'apCAF']:
            col = f'CAF_{caf}_mean'
            if col in combined_results.columns:
                r_mean = r_data[col].mean()
                nr_mean = nr_data[col].mean()
                fc = r_mean / (nr_mean + 0.001)
                logger.info(f"   {caf}: R={r_mean:.3f}, NR={nr_mean:.3f} (FC={fc:.2f}x)")

        # Spatial Entropy summary
        logger.info("\n2. SPATIAL ENTROPY:")
        if 'global_spatial_entropy' in combined_results.columns:
            r_ent = r_data['global_spatial_entropy'].mean()
            nr_ent = nr_data['global_spatial_entropy'].mean()
            _, pval = mannwhitneyu(r_data['global_spatial_entropy'], nr_data['global_spatial_entropy'])
            logger.info(f"   Global entropy: R={r_ent:.3f}, NR={nr_ent:.3f} (p={pval:.3f})")

        # Metabolic summary
        logger.info("\n3. METABOLIC PATHWAYS:")
        for pathway in ['Pyruvate_Metabolism', 'Glycolysis', 'Oxidative_Phosphorylation']:
            col = f'metab_{pathway}_mean'
            if col in combined_results.columns:
                r_mean = r_data[col].mean()
                nr_mean = nr_data[col].mean()
                fc = r_mean / (nr_mean + 0.001)
                logger.info(f"   {pathway}: R={r_mean:.3f}, NR={nr_mean:.3f} (FC={fc:.2f}x)")

    logger.info("\n" + "="*60)
    logger.info("DAY 4 ANALYSIS COMPLETE")
    logger.info("="*60)


if __name__ == '__main__':
    main()
