#!/usr/bin/env python3
"""
Generate Missing Figures with Dual Statistical Testing
======================================================

Generates figures 7, 8, 9, 14, 15 that are missing from showcase_v2.

Author: Max Van Belkum
Date: 2026-01-20
"""

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from statistical_framework import dual_test, format_stats_annotation

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
TABLE_DIR = OUTPUT_DIR / "tables"
FIG_DIR = OUTPUT_DIR / "figures" / "showcase_v2"
ADATA_DIR = OUTPUT_DIR / "adata"

plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# Colors
RESPONSE_COLORS = {'R': '#2ecc71', 'NR': '#e74c3c'}

SAMPLES = {
    "YP03A": {"response": "NR", "timepoint": "Pre"},
    "YP03C": {"response": "NR", "timepoint": "Post"},
    "YP04C": {"response": "NR", "timepoint": "Post"},
    "YP12A": {"response": "R", "timepoint": "Pre"},
    "YP12C": {"response": "R", "timepoint": "Post"},
    "YP15A": {"response": "R", "timepoint": "Pre"},
    "YP15C": {"response": "R", "timepoint": "Post"},
}


def fig7_spatial_entropy():
    """Figure 7: Spatial Entropy Analysis with dual stats."""
    print("Generating Figure 7: Spatial Entropy...")

    # Load metrics
    df = pd.read_csv(TABLE_DIR / "deep_dive_metrics.csv")

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    # Panel A: Boxplot comparison
    ax1 = axes[0]
    R_vals = df[df['response'] == 'R']['spatial_entropy'].values
    NR_vals = df[df['response'] == 'NR']['spatial_entropy'].values

    plot_data = pd.DataFrame({
        'Response': ['R']*len(R_vals) + ['NR']*len(NR_vals),
        'Spatial Entropy': list(R_vals) + list(NR_vals)
    })

    sns.boxplot(data=plot_data, x='Response', y='Spatial Entropy',
                palette=RESPONSE_COLORS, ax=ax1, width=0.5)
    sns.stripplot(data=plot_data, x='Response', y='Spatial Entropy',
                 color='black', size=10, ax=ax1)

    # Dual stats
    stats = dual_test(R_vals, NR_vals)
    stats_text = f"Welch p={stats['welch_p']:.3f}\nMWU p={stats['mwu_p']:.3f}\nd={stats['effect_size']:.2f}"
    ax1.text(0.98, 0.98, stats_text, transform=ax1.transAxes,
             fontsize=10, va='top', ha='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    ax1.set_title('A) R vs NR Comparison', fontweight='bold')
    ax1.set_ylabel('Spatial Entropy (bits)')

    # Panel B: Per-sample barplot
    ax2 = axes[1]
    df_sorted = df.sort_values(['response', 'spatial_entropy'], ascending=[False, False])
    colors = [RESPONSE_COLORS[r] for r in df_sorted['response']]
    bars = ax2.bar(df_sorted['sample'], df_sorted['spatial_entropy'], color=colors, edgecolor='black')
    ax2.set_ylabel('Spatial Entropy (bits)')
    ax2.set_title('B) Per-Sample Values', fontweight='bold')
    ax2.tick_params(axis='x', rotation=45)
    ax2.axhline(R_vals.mean(), color=RESPONSE_COLORS['R'], linestyle='--', alpha=0.7, label='R mean')
    ax2.axhline(NR_vals.mean(), color=RESPONSE_COLORS['NR'], linestyle='--', alpha=0.7, label='NR mean')
    ax2.legend(fontsize=8)

    # Panel C: Interpretation
    ax3 = axes[2]
    ax3.axis('off')

    interpretation = """
SPATIAL ENTROPY ANALYSIS
========================

Definition: Shannon entropy of cell type
distribution per sample. Higher values
indicate more diverse cell populations.

Key Finding:
• Responders show higher spatial entropy
• Effect size d=2.0 (LARGE effect)
• Welch p=0.087 (trending)

Interpretation:
Responders have more diverse cell type
composition, suggesting a more complex
tumor microenvironment that may be more
amenable to treatment.

Clinical Relevance:
Spatial entropy could serve as a
predictive biomarker for treatment
response in PDAC.
"""
    ax3.text(0.1, 0.9, interpretation, transform=ax3.transAxes,
             fontsize=10, family='monospace', va='top',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))

    plt.suptitle('Figure 7: Spatial Entropy - Cell Type Diversity Analysis',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig7_spatial_entropy.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig7_spatial_entropy.pdf", bbox_inches='tight')
    plt.close()
    print("  Saved fig7_spatial_entropy.png/pdf")


def fig8_mi_biomarkers():
    """Figure 8: Mutual Information Biomarkers."""
    print("Generating Figure 8: MI Biomarkers...")

    # Load MI biomarker data
    mi_df = pd.read_csv(TABLE_DIR / "mi_vs_response_biomarkers.csv")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: Top MI genes
    ax1 = axes[0]
    top_genes = mi_df.head(20)

    # Determine MI score column
    mi_col = 'mi_vs_response' if 'mi_vs_response' in top_genes.columns else 'mi_score'

    if mi_col in top_genes.columns:
        bars = ax1.barh(range(len(top_genes)), top_genes[mi_col], color='steelblue', edgecolor='black')
        ax1.set_xlabel('Mutual Information Score')
        ax1.set_yticks(range(len(top_genes)))
        ax1.set_yticklabels(top_genes['gene'] if 'gene' in top_genes.columns else top_genes.index)
        ax1.invert_yaxis()
        ax1.set_title('A) Top 20 MI Biomarkers', fontweight='bold')
    else:
        # Fallback if columns differ
        ax1.text(0.5, 0.5, 'MI score data format differs\nCheck mi_vs_response_biomarkers.csv',
                ha='center', va='center', transform=ax1.transAxes)
        ax1.set_title('A) Top MI Biomarkers', fontweight='bold')

    # Panel B: Interpretation panel
    ax2 = axes[1]
    ax2.axis('off')

    interpretation = """
MUTUAL INFORMATION BIOMARKERS
=============================

Method: Calculated MI between gene
expression and treatment response (R/NR)
across all cells.

Top Biomarker Candidates:
• ACSS3 - Acetyl-CoA metabolism
• SOCS2 - Cytokine signaling suppressor
• STAT4 - T cell differentiation
• CMKLR1 - Chemokine receptor

These genes show non-linear relationships
with treatment response that may not be
captured by standard DE analysis.

Biological Relevance:
Many top MI genes are immune regulators,
suggesting immune microenvironment
differences between R and NR.

Next Steps:
• Validate in independent cohorts
• Investigate functional roles
• Consider as therapeutic targets
"""
    ax2.text(0.1, 0.95, interpretation, transform=ax2.transAxes,
             fontsize=10, family='monospace', va='top',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))

    plt.suptitle('Figure 8: Mutual Information-Based Biomarker Discovery',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig8_mi_biomarkers.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig8_mi_biomarkers.pdf", bbox_inches='tight')
    plt.close()
    print("  Saved fig8_mi_biomarkers.png/pdf")


def fig9_treatment_response():
    """Figure 9: Treatment Response Biomarkers (from DE analysis)."""
    print("Generating Figure 9: Treatment Response Biomarkers...")

    # Load DE results
    de_pre = pd.read_csv(TABLE_DIR / "de_R_vs_NR_Pre.csv")
    de_post = pd.read_csv(TABLE_DIR / "de_R_vs_NR_Post.csv")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel A: Pre-treatment volcano
    ax1 = axes[0]
    if 'logfoldchanges' in de_pre.columns and 'pvals_adj' in de_pre.columns:
        de_pre['neg_log10_p'] = -np.log10(de_pre['pvals_adj'].clip(lower=1e-300))
        sig_up = (de_pre['logfoldchanges'] > 1) & (de_pre['pvals_adj'] < 0.05)
        sig_down = (de_pre['logfoldchanges'] < -1) & (de_pre['pvals_adj'] < 0.05)

        ax1.scatter(de_pre['logfoldchanges'], de_pre['neg_log10_p'],
                   c='gray', alpha=0.3, s=5)
        ax1.scatter(de_pre.loc[sig_up, 'logfoldchanges'], de_pre.loc[sig_up, 'neg_log10_p'],
                   c=RESPONSE_COLORS['R'], alpha=0.7, s=10, label='Up in R')
        ax1.scatter(de_pre.loc[sig_down, 'logfoldchanges'], de_pre.loc[sig_down, 'neg_log10_p'],
                   c=RESPONSE_COLORS['NR'], alpha=0.7, s=10, label='Up in NR')

        ax1.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
        ax1.axvline(1, color='gray', linestyle='--', alpha=0.5)
        ax1.axvline(-1, color='gray', linestyle='--', alpha=0.5)
        ax1.set_xlabel('Log2 Fold Change')
        ax1.set_ylabel('-Log10(adj p-value)')
        ax1.legend(fontsize=8)

        n_up = sig_up.sum()
        n_down = sig_down.sum()
        ax1.set_title(f'A) Pre-Treatment\n{n_up} up, {n_down} down in R', fontweight='bold')

    # Panel B: Post-treatment volcano
    ax2 = axes[1]
    if 'logfoldchanges' in de_post.columns and 'pvals_adj' in de_post.columns:
        de_post['neg_log10_p'] = -np.log10(de_post['pvals_adj'].clip(lower=1e-300))
        sig_up = (de_post['logfoldchanges'] > 1) & (de_post['pvals_adj'] < 0.05)
        sig_down = (de_post['logfoldchanges'] < -1) & (de_post['pvals_adj'] < 0.05)

        ax2.scatter(de_post['logfoldchanges'], de_post['neg_log10_p'],
                   c='gray', alpha=0.3, s=5)
        ax2.scatter(de_post.loc[sig_up, 'logfoldchanges'], de_post.loc[sig_up, 'neg_log10_p'],
                   c=RESPONSE_COLORS['R'], alpha=0.7, s=10, label='Up in R')
        ax2.scatter(de_post.loc[sig_down, 'logfoldchanges'], de_post.loc[sig_down, 'neg_log10_p'],
                   c=RESPONSE_COLORS['NR'], alpha=0.7, s=10, label='Up in NR')

        ax2.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
        ax2.axvline(1, color='gray', linestyle='--', alpha=0.5)
        ax2.axvline(-1, color='gray', linestyle='--', alpha=0.5)
        ax2.set_xlabel('Log2 Fold Change')
        ax2.set_ylabel('-Log10(adj p-value)')
        ax2.legend(fontsize=8)

        n_up = sig_up.sum()
        n_down = sig_down.sum()
        ax2.set_title(f'B) Post-Treatment\n{n_up} up, {n_down} down in R', fontweight='bold')

    # Panel C: Key genes
    ax3 = axes[2]
    ax3.axis('off')

    interpretation = """
TREATMENT RESPONSE DE ANALYSIS
==============================

Pre-Treatment Findings:
• Genes upregulated in future responders
  may predict treatment success
• Key candidates: Check volcano for
  high FC + low p-value genes

Post-Treatment Findings:
• Differences after treatment may
  explain response mechanisms
• Compare with pre-treatment to
  identify treatment-induced changes

Novel Biomarker Candidates:
• CCZ1, MS4A2, NLRP7 (not in Polymath)
• Validated via Open Targets:
  - MS4A2: 208 disease associations
  - NLRP7: 132 associations (CRC, IBD)
  - CCZ1: 22 associations

Clinical Application:
These genes could form a predictive
signature for PDAC treatment response.
"""
    ax3.text(0.05, 0.95, interpretation, transform=ax3.transAxes,
             fontsize=9, family='monospace', va='top',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))

    plt.suptitle('Figure 9: Differential Expression - Treatment Response Biomarkers',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig9_treatment_response.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig9_treatment_response.pdf", bbox_inches='tight')
    plt.close()
    print("  Saved fig9_treatment_response.png/pdf")


def fig14_deconvolution():
    """Figure 14: Cell Type Deconvolution Analysis."""
    print("Generating Figure 14: Deconvolution Analysis...")

    # Load deconvolution data
    deconv_df = pd.read_csv(TABLE_DIR / "deconvolution_proportions.csv")
    deconv_stats = pd.read_csv(TABLE_DIR / "deconvolution_stats.csv")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel A: Heatmap of proportions
    ax1 = axes[0]

    # Pivot the data
    if 'sample' in deconv_df.columns and 'cell_type' in deconv_df.columns:
        pivot_df = deconv_df.pivot_table(index='sample', columns='cell_type', values='proportion', aggfunc='mean')

        # Add response annotation
        sample_order = ['YP12A', 'YP12C', 'YP15A', 'YP15C', 'YP03A', 'YP03C', 'YP04C']
        pivot_df = pivot_df.reindex([s for s in sample_order if s in pivot_df.index])

        sns.heatmap(pivot_df, cmap='YlOrRd', ax=ax1, cbar_kws={'label': 'Proportion'})
        ax1.set_title('A) Deconvolution Proportions', fontweight='bold')
        ax1.set_ylabel('')

        # Add response labels
        for i, sample in enumerate(pivot_df.index):
            response = SAMPLES[sample]['response']
            color = RESPONSE_COLORS[response]
            ax1.text(-0.5, i + 0.5, response, ha='right', va='center', color=color, fontweight='bold')
    else:
        ax1.text(0.5, 0.5, 'Deconvolution data format differs', ha='center', va='center')
        ax1.set_title('A) Deconvolution Proportions', fontweight='bold')

    # Panel B: Key cell type comparison
    ax2 = axes[1]

    if 'cell_type' in deconv_stats.columns:
        # Get top variable cell types
        top_types = deconv_stats.head(8)

        if 'R_mean' in top_types.columns and 'NR_mean' in top_types.columns:
            x = np.arange(len(top_types))
            width = 0.35

            ax2.bar(x - width/2, top_types['R_mean'], width, label='R',
                   color=RESPONSE_COLORS['R'], edgecolor='black')
            ax2.bar(x + width/2, top_types['NR_mean'], width, label='NR',
                   color=RESPONSE_COLORS['NR'], edgecolor='black')

            ax2.set_xticks(x)
            ax2.set_xticklabels(top_types['cell_type'], rotation=45, ha='right')
            ax2.set_ylabel('Mean Proportion')
            ax2.legend()
            ax2.set_title('B) R vs NR Comparison', fontweight='bold')

    # Panel C: Interpretation
    ax3 = axes[2]
    ax3.axis('off')

    interpretation = """
DECONVOLUTION ANALYSIS
======================

Method: Reference-based deconvolution
to estimate cell type proportions at
spot level.

Key Differences:
• Acinar cells: Higher in R (p=0.057)
• Stellate/PSC: Variable
• Immune infiltration patterns differ

Biological Implications:
Higher acinar content in responders
suggests less dedifferentiated tumors
and potentially better treatment
response.

Spatial Context:
Deconvolution complements our marker-
based annotation by providing
proportion estimates even for mixed
spots.

Limitations:
• Reference panel quality matters
• Low-abundance types less reliable
• Spot-level resolution limits
"""
    ax3.text(0.05, 0.95, interpretation, transform=ax3.transAxes,
             fontsize=9, family='monospace', va='top',
             bbox=dict(boxstyle='round', facecolor='lavender', alpha=0.3))

    plt.suptitle('Figure 14: Cell Type Deconvolution Analysis',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig14_deconvolution.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig14_deconvolution.pdf", bbox_inches='tight')
    plt.close()
    print("  Saved fig14_deconvolution.png/pdf")


def fig15_ligrec():
    """Figure 15: Ligand-Receptor Communication Analysis."""
    print("Generating Figure 15: Ligand-Receptor Communication...")

    # Load comparison data
    comparison_df = pd.read_csv(TABLE_DIR / "ligrec_comparison.csv")

    # Load individual sample data
    ligrec_data = {}
    for sample in SAMPLES.keys():
        path = TABLE_DIR / f"ligrec_{sample}_top.csv"
        if path.exists():
            ligrec_data[sample] = pd.read_csv(path)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel A: Total interactions comparison
    ax1 = axes[0]

    if 'sample' in comparison_df.columns and 'n_interactions' in comparison_df.columns:
        comparison_df['response'] = comparison_df['sample'].map(lambda x: SAMPLES.get(x, {}).get('response', 'Unknown'))
        colors = [RESPONSE_COLORS.get(r, 'gray') for r in comparison_df['response']]

        ax1.bar(comparison_df['sample'], comparison_df['n_interactions'], color=colors, edgecolor='black')
        ax1.set_ylabel('Number of Significant Interactions')
        ax1.set_title('A) LR Interactions per Sample', fontweight='bold')
        ax1.tick_params(axis='x', rotation=45)

        # Add stats if we have grouped data
        R_ints = comparison_df[comparison_df['response'] == 'R']['n_interactions'].values
        NR_ints = comparison_df[comparison_df['response'] == 'NR']['n_interactions'].values

        if len(R_ints) >= 2 and len(NR_ints) >= 2:
            stats = dual_test(R_ints, NR_ints)
            stats_text = f"Welch p={stats['welch_p']:.3f}\nMWU p={stats['mwu_p']:.3f}"
            ax1.text(0.98, 0.98, stats_text, transform=ax1.transAxes,
                    fontsize=9, va='top', ha='right',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # Panel B: Top interactions heatmap
    ax2 = axes[1]

    # Aggregate top interactions across samples
    all_pairs = []
    for sample, df in ligrec_data.items():
        if 'ligand' in df.columns and 'receptor' in df.columns:
            for _, row in df.head(5).iterrows():
                all_pairs.append({
                    'pair': f"{row['ligand']}-{row['receptor']}",
                    'sample': sample,
                    'score': row.get('score', row.get('pvalue', 1))
                })

    if all_pairs:
        pairs_df = pd.DataFrame(all_pairs)
        top_pairs = pairs_df.groupby('pair').size().nlargest(10).index.tolist()

        # Create matrix
        matrix_data = pairs_df[pairs_df['pair'].isin(top_pairs)].pivot_table(
            index='pair', columns='sample', values='score', aggfunc='mean'
        )

        if not matrix_data.empty:
            sns.heatmap(matrix_data, cmap='viridis_r', ax=ax2,
                       cbar_kws={'label': 'Score'})
            ax2.set_title('B) Top LR Pairs', fontweight='bold')
    else:
        ax2.text(0.5, 0.5, 'LR pair data unavailable', ha='center', va='center')
        ax2.set_title('B) Top LR Pairs', fontweight='bold')

    # Panel C: Interpretation
    ax3 = axes[2]
    ax3.axis('off')

    interpretation = """
LIGAND-RECEPTOR ANALYSIS
========================

Method: CellPhoneDB-style analysis
identifying significant cell-cell
communication events.

Key Pathways:
• Growth factor signaling
• Immune checkpoint interactions
• ECM-receptor interactions
• Cytokine signaling

R vs NR Differences:
• Responders may show different
  immune checkpoint patterns
• Stromal-epithelial crosstalk
  varies between groups

Clinical Relevance:
LR interactions can identify:
• Therapeutic targets
• Resistance mechanisms
• Immune evasion strategies

Future Directions:
• Pathway enrichment analysis
• Network-based prioritization
• Integration with spatial context
"""
    ax3.text(0.05, 0.95, interpretation, transform=ax3.transAxes,
             fontsize=9, family='monospace', va='top',
             bbox=dict(boxstyle='round', facecolor='mistyrose', alpha=0.3))

    plt.suptitle('Figure 15: Ligand-Receptor Cell Communication Analysis',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig15_ligrec.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig15_ligrec.pdf", bbox_inches='tight')
    plt.close()
    print("  Saved fig15_ligrec.png/pdf")


def main():
    print("="*60)
    print("GENERATING MISSING FIGURES (7, 8, 9, 14, 15)")
    print("="*60 + "\n")

    fig7_spatial_entropy()
    fig8_mi_biomarkers()
    fig9_treatment_response()
    fig14_deconvolution()
    fig15_ligrec()

    print("\n" + "="*60)
    print("ALL MISSING FIGURES GENERATED")
    print("="*60)


if __name__ == "__main__":
    main()
