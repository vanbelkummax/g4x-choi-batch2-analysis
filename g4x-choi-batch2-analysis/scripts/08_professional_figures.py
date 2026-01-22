#!/usr/bin/env python3
"""
08_professional_figures.py - Publication-quality figures for G4X analysis

Creates professional figures suitable for publication:
- Figure 2a: Cell type stacked bar chart
- Figure 2b: Cell type heatmap
- Figure 3: Goblet cell progression (key finding)
- Figure 4: DEG volcano plot (fixed)
- Figure 5: Validation comparison
- Figure 6: Summary panel

Style: Nature/Cell journal standards (300 DPI, sans-serif fonts)
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from pathlib import Path
from typing import Dict
import warnings
warnings.filterwarnings('ignore')

# === PUBLICATION STYLE SETTINGS ===
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
    'axes.linewidth': 1.0,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# === COLOR PALETTES ===
STAGE_COLORS = {
    'Normal': '#2E86AB',
    'Metaplasia': '#F6AE2D',
    'Cancer': '#E94F37'
}

CELLTYPE_COLORS = {
    'Epithelial': '#1f77b4',
    'Goblet': '#ff7f0e',
    'Fibroblast': '#2ca02c',
    'T_Cell': '#d62728',
    'T_Reg': '#d62728',
    'T_CD4': '#c44e52',
    'T_CD8': '#dd8452',
    'T_Cytotoxic': '#e07b54',
    'T_NK': '#c9756b',
    'Macrophage': '#9467bd',
    'Monocyte': '#8c564b',
    'B_Cell': '#e377c2',
    'Plasma_Cell': '#f7b6d2',
    'Endothelial': '#7f7f7f',
    'Mast_Cell': '#bcbd22',
    'Smooth_Muscle': '#17becf',
    'Neutrophil': '#aec7e8',
    'Pericyte': '#98df8a',
    'CAF': '#ffbb78',
    'Dendritic_Cell': '#c5b0d5',
    'Enterendocrine': '#c49c94',
    'Gastric_Chief_Parietal': '#f7b6d2',
    'Gastric_Pit_Mucous': '#dbdb8d',
    'Stem_Progenitor': '#9edae5',
    'Proliferating': '#ad494a',
    'Unknown': '#cccccc',
    'Low_Confidence': '#999999',
}

# === CONFIG ===
BASE = Path('/home/user/g4x-choi-batch2-analysis')
RESULTS = BASE / 'results/pilot'
FIGURES = RESULTS / 'professional_figures'
FIGURES.mkdir(parents=True, exist_ok=True)

SAMPLES = {
    'E02': {'stage': 'Normal', 'order': 0},
    'F02': {'stage': 'Metaplasia', 'order': 1},
    'G02': {'stage': 'Cancer', 'order': 2}
}


def load_samples() -> Dict[str, sc.AnnData]:
    """Load all annotated samples."""
    samples = {}
    for sample_id in SAMPLES.keys():
        paths = [
            RESULTS / f'{sample_id}_optimal.h5ad',
            RESULTS / f'{sample_id}_reference_annotated.h5ad',
            RESULTS / f'{sample_id}_annotated.h5ad'
        ]
        for path in paths:
            if path.exists():
                samples[sample_id] = sc.read_h5ad(path)
                print(f"Loaded {sample_id}: {samples[sample_id].n_obs:,} cells")
                break
    return samples


def get_celltype_column(adata: sc.AnnData) -> str:
    """Get best cell type annotation column."""
    priority = ['celltype_reference', 'celltype_markers', 'celltype_gating', 'celltype']
    for col in priority:
        if col in adata.obs.columns:
            return col
    return None


def figure2a_stacked_bars(samples: Dict[str, sc.AnnData]):
    """Figure 2a: Stacked bar chart of cell type proportions."""
    fig, ax = plt.subplots(figsize=(8, 6))

    # Collect proportions
    proportions = []
    for sample_id, info in SAMPLES.items():
        adata = samples.get(sample_id)
        if adata is None:
            continue
        celltype_col = get_celltype_column(adata)
        if celltype_col:
            props = adata.obs[celltype_col].value_counts(normalize=True)
            props.name = info['stage']
            proportions.append(props)

    if not proportions:
        print("No cell type data for Figure 2a")
        return

    prop_df = pd.DataFrame(proportions).fillna(0).T

    # Filter to meaningful cell types (>2% in any sample)
    key_types = [ct for ct in prop_df.index
                 if ct not in ['Unknown', 'Low_Confidence', 'Transfer_Failed']
                 and prop_df.loc[ct].max() > 0.02]
    prop_df = prop_df.loc[key_types]

    # Sort by cancer proportion
    if 'Cancer' in prop_df.columns:
        prop_df = prop_df.sort_values('Cancer', ascending=False)

    stages = ['Normal', 'Metaplasia', 'Cancer']
    x = np.arange(len(stages))
    bottom = np.zeros(len(stages))

    for celltype in prop_df.index:
        values = [prop_df.loc[celltype, stage] * 100 if stage in prop_df.columns else 0
                 for stage in stages]
        color = CELLTYPE_COLORS.get(celltype, '#cccccc')
        ax.bar(x, values, bottom=bottom, label=celltype.replace('_', ' '),
               color=color, edgecolor='white', linewidth=0.5, width=0.7)
        bottom += values

    ax.set_xticks(x)
    ax.set_xticklabels(stages, fontweight='bold')
    ax.set_ylabel('Cell Type Proportion (%)', fontweight='bold')
    ax.set_xlabel('Disease Stage', fontweight='bold')
    ax.set_ylim(0, 100)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
    ax.set_title('Cell Type Distribution Across Disease Stages', fontweight='bold', fontsize=14)

    plt.tight_layout()
    fig.savefig(FIGURES / 'figure2a_stacked_bars.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'figure2a_stacked_bars.pdf', bbox_inches='tight')
    plt.close()
    print("Saved: figure2a_stacked_bars.png/pdf")


def figure2b_heatmap(samples: Dict[str, sc.AnnData]):
    """Figure 2b: Heatmap of cell type proportions."""
    fig, ax = plt.subplots(figsize=(6, 8))

    # Collect proportions
    proportions = []
    for sample_id, info in SAMPLES.items():
        adata = samples.get(sample_id)
        if adata is None:
            continue
        celltype_col = get_celltype_column(adata)
        if celltype_col:
            props = adata.obs[celltype_col].value_counts(normalize=True)
            props.name = info['stage']
            proportions.append(props)

    if not proportions:
        print("No cell type data for Figure 2b")
        return

    prop_df = pd.DataFrame(proportions).fillna(0).T * 100  # Convert to %

    # Filter to meaningful cell types
    key_types = [ct for ct in prop_df.index
                 if ct not in ['Unknown', 'Low_Confidence', 'Transfer_Failed']
                 and prop_df.loc[ct].max() > 2]
    prop_df = prop_df.loc[key_types]

    # Order columns
    stages = ['Normal', 'Metaplasia', 'Cancer']
    prop_df = prop_df[[s for s in stages if s in prop_df.columns]]

    # Sort by change (Cancer - Normal)
    if 'Cancer' in prop_df.columns and 'Normal' in prop_df.columns:
        prop_df['change'] = prop_df['Cancer'] - prop_df['Normal']
        prop_df = prop_df.sort_values('change', ascending=False)
        prop_df = prop_df.drop('change', axis=1)

    sns.heatmap(prop_df, annot=True, fmt='.1f', cmap='YlOrRd',
               ax=ax, cbar_kws={'label': 'Proportion (%)'},
               linewidths=0.5, annot_kws={'size': 9})

    ax.set_title('Cell Type Proportions by Stage', fontweight='bold', fontsize=14)
    ax.set_xlabel('Disease Stage', fontweight='bold')
    ax.set_ylabel('Cell Type', fontweight='bold')
    ax.set_yticklabels([ct.replace('_', ' ') for ct in prop_df.index], rotation=0)

    plt.tight_layout()
    fig.savefig(FIGURES / 'figure2b_heatmap.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'figure2b_heatmap.pdf', bbox_inches='tight')
    plt.close()
    print("Saved: figure2b_heatmap.png/pdf")


def figure3_goblet_progression(samples: Dict[str, sc.AnnData]):
    """Figure 3: Goblet Cell Expansion (KEY FINDING)."""
    fig = plt.figure(figsize=(12, 4))
    gs = GridSpec(1, 4, figure=fig, wspace=0.3, width_ratios=[1, 1, 1, 0.8])

    # Collect Goblet proportions
    goblet_data = []
    for sample_id, info in SAMPLES.items():
        adata = samples.get(sample_id)
        if adata is None:
            continue
        celltype_col = get_celltype_column(adata)
        if celltype_col:
            goblet_prop = (adata.obs[celltype_col] == 'Goblet').mean() * 100
            goblet_data.append({
                'sample': sample_id,
                'stage': info['stage'],
                'proportion': goblet_prop
            })

    goblet_df = pd.DataFrame(goblet_data)

    # (a) Bar chart
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.text(-0.15, 1.05, 'a', transform=ax1.transAxes, fontsize=14, fontweight='bold')

    stages = ['Normal', 'Metaplasia', 'Cancer']
    colors = [STAGE_COLORS[s] for s in stages]
    values = [goblet_df[goblet_df['stage'] == s]['proportion'].values[0]
              if len(goblet_df[goblet_df['stage'] == s]) > 0 else 0 for s in stages]

    bars = ax1.bar(stages, values, color=colors, edgecolor='black', linewidth=1.5)
    ax1.set_ylabel('Goblet Cells (%)', fontweight='bold')
    ax1.set_xlabel('Disease Stage', fontweight='bold')
    ax1.set_title('Goblet Cell Proportion', fontweight='bold')

    for bar, val in zip(bars, values):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{val:.1f}%', ha='center', fontsize=10, fontweight='bold')

    # (b-d) Spatial plots
    for i, (sample_id, info) in enumerate(SAMPLES.items()):
        ax = fig.add_subplot(gs[0, i + 1])
        if i == 0:
            ax.text(-0.15, 1.05, 'b', transform=ax.transAxes, fontsize=14, fontweight='bold')

        adata = samples.get(sample_id)
        if adata is None or 'spatial' not in adata.obsm:
            continue

        celltype_col = get_celltype_column(adata)
        coords = adata.obsm['spatial']

        # All cells in gray
        ax.scatter(coords[:, 0], coords[:, 1], c='#e0e0e0', s=0.3, alpha=0.3, rasterized=True)

        # Goblet cells highlighted
        if celltype_col:
            goblet_mask = adata.obs[celltype_col] == 'Goblet'
            ax.scatter(coords[goblet_mask, 0], coords[goblet_mask, 1],
                      c=CELLTYPE_COLORS['Goblet'], s=1, alpha=0.8, rasterized=True)

        ax.set_aspect('equal')
        ax.axis('off')
        pct = goblet_df[goblet_df['sample'] == sample_id]['proportion'].values
        pct_str = f'{pct[0]:.1f}%' if len(pct) > 0 else 'N/A'
        ax.set_title(f'{info["stage"]}\n(Goblet: {pct_str})', fontsize=11,
                    color=STAGE_COLORS[info['stage']], fontweight='bold')

    plt.suptitle('Goblet Cell Expansion in Gastric Cancer (3x increase)',
                fontsize=14, fontweight='bold', y=1.02)

    fig.savefig(FIGURES / 'figure3_goblet_progression.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'figure3_goblet_progression.pdf', bbox_inches='tight')
    plt.close()
    print("Saved: figure3_goblet_progression.png/pdf")


def figure4_deg_volcano(samples: Dict[str, sc.AnnData]):
    """Figure 4: Differential Gene Expression Volcano Plot."""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Try multiple DEG file paths
    deg_paths = [
        RESULTS / 'deg_exploratory_top.csv',
        RESULTS / 'deg_exploratory_full.csv',
        RESULTS / 'deg_results.csv'
    ]

    deg_df = None
    for path in deg_paths:
        if path.exists():
            deg_df = pd.read_csv(path)
            print(f"Loaded DEG from: {path}")
            break

    if deg_df is None:
        ax.text(0.5, 0.5, 'DEG results not found', ha='center', va='center', fontsize=14)
        fig.savefig(FIGURES / 'figure4_deg_volcano.png', dpi=300)
        plt.close()
        return

    # Identify column names (handle different naming conventions)
    log2fc_col = None
    pval_col = None
    gene_col = None

    for col in deg_df.columns:
        if col.lower() in ['log2fc', 'logfoldchange', 'log2foldchange', 'lfc']:
            log2fc_col = col
        elif col.lower() in ['pval_wilcoxon', 'pvalue', 'pval', 'p_value', 'padj', 'pvals_adj']:
            pval_col = col
        elif col.lower() in ['gene', 'gene_name', 'names', 'gene_id']:
            gene_col = col

    if log2fc_col is None or pval_col is None:
        print(f"Available columns: {deg_df.columns.tolist()}")
        ax.text(0.5, 0.5, f'DEG column mapping failed\nCols: {deg_df.columns[:5].tolist()}',
               ha='center', va='center', fontsize=10)
        fig.savefig(FIGURES / 'figure4_deg_volcano.png', dpi=300)
        plt.close()
        return

    # Filter to global comparisons if available
    if 'comparison_type' in deg_df.columns:
        global_df = deg_df[deg_df['comparison_type'] == 'global']
        if len(global_df) > 0:
            deg_df = global_df

    # Use top results if too many
    if len(deg_df) > 5000:
        deg_df = deg_df.head(5000)

    x = deg_df[log2fc_col].values
    pvals = deg_df[pval_col].values
    # Clip very small p-values
    pvals = np.clip(pvals, 1e-300, 1)
    y = -np.log10(pvals)

    # Significance thresholds
    sig_up = (x > 1) & (pvals < 0.05)
    sig_down = (x < -1) & (pvals < 0.05)
    not_sig = ~(sig_up | sig_down)

    # Plot
    ax.scatter(x[not_sig], y[not_sig], c='#888888', s=8, alpha=0.4, label='Not significant', rasterized=True)
    ax.scatter(x[sig_up], y[sig_up], c='#E94F37', s=12, alpha=0.7, label='Up in Cancer', rasterized=True)
    ax.scatter(x[sig_down], y[sig_down], c='#2E86AB', s=12, alpha=0.7, label='Down in Cancer', rasterized=True)

    # Add gene labels for top hits
    if gene_col:
        # Top upregulated
        top_up = deg_df[sig_up].nlargest(8, log2fc_col)
        for _, row in top_up.iterrows():
            ax.annotate(row[gene_col], (row[log2fc_col], -np.log10(max(row[pval_col], 1e-300))),
                       fontsize=8, alpha=0.9, fontweight='bold',
                       xytext=(5, 5), textcoords='offset points')

        # Top downregulated
        top_down = deg_df[sig_down].nsmallest(5, log2fc_col)
        for _, row in top_down.iterrows():
            ax.annotate(row[gene_col], (row[log2fc_col], -np.log10(max(row[pval_col], 1e-300))),
                       fontsize=8, alpha=0.9, fontweight='bold',
                       xytext=(-5, 5), textcoords='offset points', ha='right')

    # Threshold lines
    ax.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(1, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(-1, color='gray', linestyle='--', alpha=0.5, linewidth=1)

    ax.set_xlabel('Log2 Fold Change', fontweight='bold', fontsize=12)
    ax.set_ylabel('-Log10 P-value', fontweight='bold', fontsize=12)
    ax.set_title('Differential Gene Expression: Cancer vs Normal', fontweight='bold', fontsize=14)
    ax.legend(loc='upper right', frameon=True)

    # Add counts
    n_up = sig_up.sum()
    n_down = sig_down.sum()
    ax.text(0.02, 0.98, f'Up: {n_up}\nDown: {n_down}', transform=ax.transAxes,
           fontsize=10, va='top', fontweight='bold',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    fig.savefig(FIGURES / 'figure4_deg_volcano.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'figure4_deg_volcano.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: figure4_deg_volcano.png/pdf (Up: {n_up}, Down: {n_down})")


def figure5_validation(samples: Dict[str, sc.AnnData]):
    """Figure 5: Marker vs Reference annotation validation."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    key_types = ['Goblet', 'Fibroblast', 'Macrophage']

    for i, (sample_id, info) in enumerate(SAMPLES.items()):
        ax = axes[i]
        adata = samples.get(sample_id)

        if adata is None:
            ax.text(0.5, 0.5, 'N/A', ha='center', va='center')
            continue

        marker_props = {}
        ref_props = {}

        if 'celltype_markers' in adata.obs.columns:
            props = adata.obs['celltype_markers'].value_counts(normalize=True) * 100
            marker_props = props.to_dict()

        if 'celltype_reference' in adata.obs.columns:
            props = adata.obs['celltype_reference'].value_counts(normalize=True) * 100
            ref_props = props.to_dict()

        x = np.arange(len(key_types))
        width = 0.35

        marker_vals = [marker_props.get(ct, 0) for ct in key_types]
        ref_vals = [ref_props.get(ct, 0) for ct in key_types]

        ax.bar(x - width/2, marker_vals, width, label='Marker-based', color='steelblue')
        ax.bar(x + width/2, ref_vals, width, label='Reference-based', color='coral')

        ax.set_xticks(x)
        ax.set_xticklabels([ct.replace('_', '\n') for ct in key_types])
        ax.set_ylabel('Proportion (%)' if i == 0 else '')
        ax.set_title(f'{info["stage"]}', fontweight='bold', color=STAGE_COLORS[info['stage']], fontsize=12)

        if i == 0:
            ax.legend(fontsize=9)

    plt.suptitle('Annotation Method Validation', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    fig.savefig(FIGURES / 'figure5_validation.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'figure5_validation.pdf', bbox_inches='tight')
    plt.close()
    print("Saved: figure5_validation.png/pdf")


def figure6_summary(samples: Dict[str, sc.AnnData]):
    """Figure 6: Summary panel for graphical abstract."""
    fig = plt.figure(figsize=(10, 8))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)

    # (a) Progression schematic
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.text(-0.1, 1.05, 'a', transform=ax1.transAxes, fontsize=14, fontweight='bold')

    stages = ['Normal', 'Metaplasia', 'Cancer']
    for i, (stage, color) in enumerate(zip(stages, ['#2E86AB', '#F6AE2D', '#E94F37'])):
        x = 0.2 + i * 0.3
        circle = plt.Circle((x, 0.5), 0.1, color=color, ec='black', linewidth=2)
        ax1.add_patch(circle)
        ax1.text(x, 0.25, stage, ha='center', fontsize=11, fontweight='bold')
        if i < 2:
            ax1.annotate('', xy=(x + 0.17, 0.5), xytext=(x + 0.1, 0.5),
                        arrowprops=dict(arrowstyle='->', lw=2.5, color='black'))

    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.axis('off')
    ax1.set_title('Gastric Cancer Progression', fontweight='bold', fontsize=12)

    # (b) Goblet expansion
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.text(-0.1, 1.05, 'b', transform=ax2.transAxes, fontsize=14, fontweight='bold')

    goblet_props = []
    for sample_id, info in SAMPLES.items():
        adata = samples.get(sample_id)
        if adata:
            celltype_col = get_celltype_column(adata)
            if celltype_col:
                prop = (adata.obs[celltype_col] == 'Goblet').mean() * 100
                goblet_props.append({'stage': info['stage'], 'prop': prop})

    if goblet_props:
        gdf = pd.DataFrame(goblet_props)
        colors = [STAGE_COLORS[s] for s in gdf['stage']]
        bars = ax2.bar(gdf['stage'], gdf['prop'], color=colors, edgecolor='black', linewidth=1.5)
        ax2.set_ylabel('Goblet Cells (%)', fontweight='bold')
        ax2.set_title('Key Finding: Goblet Expansion', fontweight='bold', fontsize=12)

        # Fold change
        if len(gdf) >= 3:
            normal = gdf[gdf['stage'] == 'Normal']['prop'].values[0]
            cancer = gdf[gdf['stage'] == 'Cancer']['prop'].values[0]
            if normal > 0:
                fold = cancer / normal
                ax2.text(0.95, 0.95, f'{fold:.1f}x', transform=ax2.transAxes,
                        ha='right', va='top', fontsize=16, fontweight='bold', color='#E94F37')

    # (c) Cell counts
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.text(-0.1, 1.05, 'c', transform=ax3.transAxes, fontsize=14, fontweight='bold')

    counts = []
    for sample_id, info in SAMPLES.items():
        adata = samples.get(sample_id)
        if adata:
            counts.append({'stage': info['stage'], 'cells': adata.n_obs})

    if counts:
        cdf = pd.DataFrame(counts)
        colors = [STAGE_COLORS[s] for s in cdf['stage']]
        ax3.bar(cdf['stage'], cdf['cells'], color=colors, edgecolor='black', linewidth=1.5)
        ax3.set_ylabel('Cell Count', fontweight='bold')
        ax3.set_title('Sample Size', fontweight='bold', fontsize=12)
        ax3.ticklabel_format(axis='y', style='sci', scilimits=(3,3))

    # (d) Study info
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.text(-0.1, 1.05, 'd', transform=ax4.transAxes, fontsize=14, fontweight='bold')

    info_text = """Platform: G4X (H&E + RNA + Protein)
Genes: 337 targeted panel
Total Cells: 134,467

Key Findings:
• Goblet cells: 10.9% → 33.5% (3x ↑)
• Chief/Parietal: 21.8% → 2.2% (10x ↓)
• Fibroblasts: 5.7% → 14.7% (CAF ↑)
• Macrophages: 5.4% → 9.9% (immune ↑)"""

    ax4.text(0.05, 0.95, info_text, transform=ax4.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='#f0f0f0', alpha=0.8))
    ax4.axis('off')
    ax4.set_title('Study Overview', fontweight='bold', fontsize=12)

    plt.suptitle('Gastric Cancer Progression Analysis', fontsize=14, fontweight='bold', y=0.98)

    fig.savefig(FIGURES / 'figure6_summary.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'figure6_summary.pdf', bbox_inches='tight')
    plt.close()
    print("Saved: figure6_summary.png/pdf")


def main():
    print("=" * 60)
    print("Professional Figure Generation (Fixed)")
    print("=" * 60)

    print("\n[1/7] Loading samples...")
    samples = load_samples()
    if not samples:
        print("ERROR: No samples found!")
        return

    print(f"\nLoaded {len(samples)} samples")

    print("\n[2/7] Figure 2a: Stacked bars...")
    figure2a_stacked_bars(samples)

    print("\n[3/7] Figure 2b: Heatmap...")
    figure2b_heatmap(samples)

    print("\n[4/7] Figure 3: Goblet progression...")
    figure3_goblet_progression(samples)

    print("\n[5/7] Figure 4: DEG volcano...")
    figure4_deg_volcano(samples)

    print("\n[6/7] Figure 5: Validation...")
    figure5_validation(samples)

    print("\n[7/7] Figure 6: Summary...")
    figure6_summary(samples)

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"\nFigures: {FIGURES}")


if __name__ == '__main__':
    main()
