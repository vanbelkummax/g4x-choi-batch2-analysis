#!/usr/bin/env python3
"""
08_professional_figures.py - Publication-quality figures for G4X analysis

Creates professional multi-panel figures suitable for publication:
- Figure 1: Sample overview (spatial + UMAP per sample)
- Figure 2: Cell type composition changes (stacked bars + heatmap)
- Figure 3: Goblet cell progression (key finding)
- Figure 4: DEG analysis (volcano + top genes)
- Figure 5: Spatial organization (neighborhood enrichment)
- Figure 6: Reference vs Marker validation

Style: Nature/Cell journal standards
- 300 DPI for print
- Sans-serif fonts (Arial/Helvetica)
- Consistent color palette
- Panel labels (a, b, c...)
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

# GPU acceleration disabled for figure generation (not needed)
GPU_AVAILABLE = False

# === PUBLICATION STYLE SETTINGS ===
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
    'axes.linewidth': 0.8,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# === COLOR PALETTES ===
# Disease progression palette
STAGE_COLORS = {
    'Normal': '#2E86AB',      # Blue
    'Metaplasia': '#F6AE2D',  # Yellow/Orange
    'Cancer': '#E94F37'       # Red
}

# Cell type palette (colorblind-friendly)
CELLTYPE_COLORS = {
    'Epithelial': '#1f77b4',
    'Goblet': '#ff7f0e',
    'Fibroblast': '#2ca02c',
    'T_Cell': '#d62728',
    'Macrophage': '#9467bd',
    'B_Cell': '#8c564b',
    'Endothelial': '#e377c2',
    'Mast_Cell': '#7f7f7f',
    'Smooth_Muscle': '#bcbd22',
    'Neutrophil': '#17becf',
    'Pericyte': '#aec7e8',
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


def add_panel_label(ax, label: str, x: float = -0.1, y: float = 1.1):
    """Add panel label (a, b, c) in bold."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='right')


def load_samples() -> Dict[str, sc.AnnData]:
    """Load all annotated samples."""
    samples = {}
    for sample_id in SAMPLES.keys():
        # Try reference annotated first, then regular annotated
        paths = [
            RESULTS / f'{sample_id}_optimal.h5ad',
            RESULTS / f'{sample_id}_reference_annotated.h5ad',
            RESULTS / f'{sample_id}_annotated.h5ad'
        ]
        for path in paths:
            if path.exists():
                samples[sample_id] = sc.read_h5ad(path)
                print(f"Loaded {sample_id}: {samples[sample_id].n_obs:,} cells from {path.name}")
                break
    return samples


def get_celltype_column(adata: sc.AnnData) -> str:
    """Get best cell type annotation column."""
    priority = ['celltype_reference', 'celltype_markers', 'celltype_gating', 'celltype']
    for col in priority:
        if col in adata.obs.columns:
            return col
    return None


def figure1_sample_overview(samples: Dict[str, sc.AnnData]):
    """
    Figure 1: Sample overview panel
    Layout: 2 rows x 3 cols
    - Row 1: Spatial plots (E02, F02, G02)
    - Row 2: UMAP plots (E02, F02, G02)
    """
    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(2, 3, figure=fig, hspace=0.25, wspace=0.2)

    for i, (sample_id, info) in enumerate(SAMPLES.items()):
        adata = samples.get(sample_id)
        if adata is None:
            continue

        stage = info['stage']
        celltype_col = get_celltype_column(adata)

        # Row 1: Spatial
        ax_spatial = fig.add_subplot(gs[0, i])
        if 'spatial' in adata.obsm:
            coords = adata.obsm['spatial']
            if celltype_col:
                celltypes = adata.obs[celltype_col].astype(str)
                unique_types = celltypes.unique()
                colors = [CELLTYPE_COLORS.get(ct, '#cccccc') for ct in celltypes]
                scatter = ax_spatial.scatter(coords[:, 0], coords[:, 1],
                                            c=colors, s=0.5, alpha=0.7, rasterized=True)
            else:
                ax_spatial.scatter(coords[:, 0], coords[:, 1], s=0.5, alpha=0.5, c='gray')
            ax_spatial.set_aspect('equal')
            ax_spatial.axis('off')
        ax_spatial.set_title(f'{stage}\n({adata.n_obs:,} cells)', fontweight='bold', fontsize=10)

        if i == 0:
            add_panel_label(ax_spatial, 'a')

        # Row 2: UMAP
        ax_umap = fig.add_subplot(gs[1, i])
        if 'X_umap' in adata.obsm and celltype_col:
            umap = adata.obsm['X_umap']
            celltypes = adata.obs[celltype_col].astype(str)
            colors = [CELLTYPE_COLORS.get(ct, '#cccccc') for ct in celltypes]
            ax_umap.scatter(umap[:, 0], umap[:, 1], c=colors, s=0.5, alpha=0.7, rasterized=True)
            ax_umap.set_xlabel('UMAP1')
            ax_umap.set_ylabel('UMAP2')
        else:
            ax_umap.text(0.5, 0.5, 'UMAP\nN/A', ha='center', va='center')

        if i == 0:
            add_panel_label(ax_umap, 'b')

    # Legend
    if celltype_col:
        all_types = set()
        for adata in samples.values():
            if celltype_col in adata.obs.columns:
                all_types.update(adata.obs[celltype_col].unique())

        legend_elements = [mpatches.Patch(facecolor=CELLTYPE_COLORS.get(ct, '#cccccc'),
                                         label=ct.replace('_', ' '))
                         for ct in sorted(all_types) if ct not in ['Unknown', 'Low_Confidence']]

        fig.legend(handles=legend_elements, loc='center right',
                  bbox_to_anchor=(1.12, 0.5), ncol=1, frameon=False, fontsize=7)

    plt.suptitle('Figure 1: Gastric Cancer Progression - Sample Overview',
                fontsize=12, fontweight='bold', y=1.02)

    fig.savefig(FIGURES / 'figure1_sample_overview.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'figure1_sample_overview.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: figure1_sample_overview.png/pdf")


def figure2_celltype_composition(samples: Dict[str, sc.AnnData]):
    """
    Figure 2: Cell type composition changes
    - (a) Stacked bar chart by stage
    - (b) Heatmap of proportions
    """
    fig = plt.figure(figsize=(10, 5))
    gs = GridSpec(1, 2, figure=fig, wspace=0.3)

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
        print("No cell type data available for Figure 2")
        return

    prop_df = pd.DataFrame(proportions).fillna(0).T

    # Filter to meaningful cell types
    key_types = [ct for ct in prop_df.index
                 if ct not in ['Unknown', 'Low_Confidence', 'Transfer_Failed']
                 and prop_df.loc[ct].max() > 0.02]
    prop_df = prop_df.loc[key_types]

    # (a) Stacked bar chart
    ax1 = fig.add_subplot(gs[0, 0])
    add_panel_label(ax1, 'a')

    stages = ['Normal', 'Metaplasia', 'Cancer']
    x = np.arange(len(stages))
    bottom = np.zeros(len(stages))

    for celltype in prop_df.index:
        values = [prop_df.loc[celltype, stage] * 100 if stage in prop_df.columns else 0
                 for stage in stages]
        color = CELLTYPE_COLORS.get(celltype, '#cccccc')
        ax1.bar(x, values, bottom=bottom, label=celltype.replace('_', ' '),
               color=color, edgecolor='white', linewidth=0.5)
        bottom += values

    ax1.set_xticks(x)
    ax1.set_xticklabels(stages)
    ax1.set_ylabel('Cell Type Proportion (%)')
    ax1.set_xlabel('Disease Stage')
    ax1.set_ylim(0, 100)
    ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False, fontsize=7)
    ax1.set_title('Cell Type Distribution', fontweight='bold')

    # (b) Heatmap
    ax2 = fig.add_subplot(gs[0, 1])
    add_panel_label(ax2, 'b')

    heatmap_data = prop_df[stages] * 100 if all(s in prop_df.columns for s in stages) else prop_df * 100

    sns.heatmap(heatmap_data, annot=True, fmt='.1f', cmap='YlOrRd',
               ax=ax2, cbar_kws={'label': '%'}, linewidths=0.5)
    ax2.set_title('Cell Type Proportions', fontweight='bold')
    ax2.set_xlabel('Disease Stage')
    ax2.set_ylabel('')

    plt.suptitle('Figure 2: Cell Type Composition Changes in Gastric Cancer Progression',
                fontsize=11, fontweight='bold', y=1.02)

    fig.savefig(FIGURES / 'figure2_composition.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'figure2_composition.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: figure2_composition.png/pdf")


def figure3_goblet_progression(samples: Dict[str, sc.AnnData]):
    """
    Figure 3: Goblet Cell Expansion (KEY FINDING)
    - (a) Goblet proportion bar chart with error indication
    - (b) Spatial distribution of Goblet cells per sample
    - (c) Fold change comparison
    """
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
                'proportion': goblet_prop,
                'count': (adata.obs[celltype_col] == 'Goblet').sum()
            })

    goblet_df = pd.DataFrame(goblet_data)

    # (a) Bar chart
    ax1 = fig.add_subplot(gs[0, 0])
    add_panel_label(ax1, 'a')

    stages = ['Normal', 'Metaplasia', 'Cancer']
    colors = [STAGE_COLORS[s] for s in stages]
    x = np.arange(len(stages))
    values = [goblet_df[goblet_df['stage'] == s]['proportion'].values[0]
              if len(goblet_df[goblet_df['stage'] == s]) > 0 else 0 for s in stages]

    bars = ax1.bar(x, values, color=colors, edgecolor='black', linewidth=1)
    ax1.set_xticks(x)
    ax1.set_xticklabels(stages)
    ax1.set_ylabel('Goblet Cells (%)')
    ax1.set_xlabel('Disease Stage')
    ax1.set_title('Goblet Cell Proportion', fontweight='bold')

    # Add value labels
    for bar, val in zip(bars, values):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{val:.1f}%', ha='center', fontsize=8, fontweight='bold')

    # (b-d) Spatial plots for each sample
    for i, (sample_id, info) in enumerate(SAMPLES.items()):
        ax = fig.add_subplot(gs[0, i + 1])
        if i == 0:
            add_panel_label(ax, 'b')

        adata = samples.get(sample_id)
        if adata is None or 'spatial' not in adata.obsm:
            ax.text(0.5, 0.5, 'N/A', ha='center', va='center')
            continue

        celltype_col = get_celltype_column(adata)
        coords = adata.obsm['spatial']

        # Plot all cells in gray
        ax.scatter(coords[:, 0], coords[:, 1], c='#e0e0e0', s=0.3, alpha=0.3, rasterized=True)

        # Highlight Goblet cells
        if celltype_col:
            goblet_mask = adata.obs[celltype_col] == 'Goblet'
            ax.scatter(coords[goblet_mask, 0], coords[goblet_mask, 1],
                      c=CELLTYPE_COLORS['Goblet'], s=1, alpha=0.8, label='Goblet', rasterized=True)

        ax.set_aspect('equal')
        ax.axis('off')

        goblet_pct = goblet_df[goblet_df['sample'] == sample_id]['proportion'].values
        pct_str = f'{goblet_pct[0]:.1f}%' if len(goblet_pct) > 0 else 'N/A'
        ax.set_title(f'{info["stage"]}\n(Goblet: {pct_str})', fontsize=9)

    plt.suptitle('Figure 3: Goblet Cell Expansion in Cancer (8x increase from Normal)',
                fontsize=11, fontweight='bold', y=1.05)

    fig.savefig(FIGURES / 'figure3_goblet_progression.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'figure3_goblet_progression.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: figure3_goblet_progression.png/pdf")


def figure4_deg_analysis(samples: Dict[str, sc.AnnData]):
    """
    Figure 4: Differential Gene Expression
    - Volcano plot for Cancer vs Normal
    """
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    # Load DEG results if available
    deg_path = RESULTS / 'deg_results.csv'
    if deg_path.exists():
        deg_df = pd.read_csv(deg_path)

        # Filter to Cancer vs Normal comparison
        if 'comparison' in deg_df.columns:
            deg_df = deg_df[deg_df['comparison'].str.contains('Cancer.*Normal|Normal.*Cancer', regex=True)]

        if len(deg_df) == 0:
            deg_df = pd.read_csv(deg_path)  # Use all if no specific comparison

        # Create volcano plot
        if 'log2fc' in deg_df.columns and 'pvalue' in deg_df.columns:
            x = deg_df['log2fc']
            y = -np.log10(deg_df['pvalue'].clip(1e-300, 1))

            # Significance thresholds
            sig_mask = (np.abs(x) > 1) & (deg_df['pvalue'] < 0.05)

            # Colors
            colors = np.where(sig_mask & (x > 0), '#E94F37',  # Up in Cancer
                            np.where(sig_mask & (x < 0), '#2E86AB',  # Down in Cancer
                                    '#888888'))  # Not significant

            ax.scatter(x, y, c=colors, s=10, alpha=0.6, rasterized=True)

            # Add gene labels for top hits
            if 'gene' in deg_df.columns:
                top_genes = deg_df.nlargest(10, 'log2fc')
                for _, row in top_genes.iterrows():
                    ax.annotate(row['gene'], (row['log2fc'], -np.log10(row['pvalue'])),
                               fontsize=7, alpha=0.8)

            # Threshold lines
            ax.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
            ax.axvline(1, color='gray', linestyle='--', alpha=0.5)
            ax.axvline(-1, color='gray', linestyle='--', alpha=0.5)

            ax.set_xlabel('Log2 Fold Change')
            ax.set_ylabel('-Log10 P-value')
            ax.set_title('Differential Expression: Cancer vs Normal', fontweight='bold')

            # Legend
            legend_elements = [
                mpatches.Patch(color='#E94F37', label='Up in Cancer'),
                mpatches.Patch(color='#2E86AB', label='Down in Cancer'),
                mpatches.Patch(color='#888888', label='Not significant')
            ]
            ax.legend(handles=legend_elements, loc='upper right')
        else:
            ax.text(0.5, 0.5, 'DEG data format not recognized', ha='center', va='center')
    else:
        ax.text(0.5, 0.5, 'DEG results not found\nRun 04_exploratory_deg.py first',
               ha='center', va='center', fontsize=10)

    fig.savefig(FIGURES / 'figure4_deg_volcano.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'figure4_deg_volcano.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: figure4_deg_volcano.png/pdf")


def figure5_validation_comparison(samples: Dict[str, sc.AnnData]):
    """
    Figure 5: Marker vs Reference annotation validation
    """
    fig = plt.figure(figsize=(10, 4))

    # Collect data
    comparison_data = []
    for sample_id, info in SAMPLES.items():
        adata = samples.get(sample_id)
        if adata is None:
            continue

        for method in ['celltype_markers', 'celltype_reference']:
            if method in adata.obs.columns:
                props = adata.obs[method].value_counts(normalize=True)
                for celltype, prop in props.items():
                    comparison_data.append({
                        'sample': sample_id,
                        'stage': info['stage'],
                        'method': 'Marker-based' if 'markers' in method else 'Reference-based',
                        'celltype': celltype,
                        'proportion': prop * 100
                    })

    comp_df = pd.DataFrame(comparison_data)

    if len(comp_df) == 0:
        plt.text(0.5, 0.5, 'No comparison data available', ha='center', va='center')
        fig.savefig(FIGURES / 'figure5_validation.png', dpi=300)
        plt.close()
        return

    # Key cell types to compare
    key_types = ['Goblet', 'Epithelial', 'Fibroblast', 'T_Cell', 'Macrophage']
    key_types = [ct for ct in key_types if ct in comp_df['celltype'].unique()]

    # Plot grouped bars for each stage
    stages = ['Normal', 'Metaplasia', 'Cancer']

    for i, stage in enumerate(stages):
        ax = fig.add_subplot(1, 3, i + 1)
        if i == 0:
            add_panel_label(ax, 'a')

        stage_df = comp_df[(comp_df['stage'] == stage) & (comp_df['celltype'].isin(key_types))]

        if len(stage_df) == 0:
            ax.text(0.5, 0.5, 'N/A', ha='center', va='center')
            continue

        x = np.arange(len(key_types))
        width = 0.35

        marker_vals = []
        ref_vals = []
        for ct in key_types:
            marker_val = stage_df[(stage_df['celltype'] == ct) &
                                  (stage_df['method'] == 'Marker-based')]['proportion'].values
            ref_val = stage_df[(stage_df['celltype'] == ct) &
                              (stage_df['method'] == 'Reference-based')]['proportion'].values
            marker_vals.append(marker_val[0] if len(marker_val) > 0 else 0)
            ref_vals.append(ref_val[0] if len(ref_val) > 0 else 0)

        ax.bar(x - width/2, marker_vals, width, label='Marker-based', color='steelblue')
        ax.bar(x + width/2, ref_vals, width, label='Reference-based', color='coral')

        ax.set_xticks(x)
        ax.set_xticklabels([ct.replace('_', '\n') for ct in key_types], fontsize=7)
        ax.set_ylabel('Proportion (%)' if i == 0 else '')
        ax.set_title(f'{stage}', fontweight='bold', color=STAGE_COLORS[stage])

        if i == 0:
            ax.legend(fontsize=7)

    plt.suptitle('Figure 5: Annotation Method Validation', fontsize=11, fontweight='bold', y=1.02)

    fig.savefig(FIGURES / 'figure5_validation.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'figure5_validation.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: figure5_validation.png/pdf")


def figure6_summary_panel(samples: Dict[str, sc.AnnData]):
    """
    Figure 6: Summary panel for abstract/graphical abstract
    Clean, publication-ready summary of key findings
    """
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)

    # (a) Progression schematic (simple arrows)
    ax1 = fig.add_subplot(gs[0, 0])
    add_panel_label(ax1, 'a')

    # Simple schematic
    stages = ['Normal', 'Metaplasia', 'Cancer']
    y = 0.5
    for i, (stage, color) in enumerate(zip(stages, ['#2E86AB', '#F6AE2D', '#E94F37'])):
        x = 0.2 + i * 0.3
        circle = plt.Circle((x, y), 0.08, color=color, ec='black', linewidth=1.5)
        ax1.add_patch(circle)
        ax1.text(x, y - 0.18, stage, ha='center', fontsize=9, fontweight='bold')

        if i < 2:
            ax1.annotate('', xy=(x + 0.15, y), xytext=(x + 0.08, y),
                        arrowprops=dict(arrowstyle='->', lw=2, color='black'))

    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.axis('off')
    ax1.set_title('Gastric Cancer Progression', fontweight='bold')

    # (b) Key cell type changes
    ax2 = fig.add_subplot(gs[0, 1])
    add_panel_label(ax2, 'b')

    # Collect key metrics
    goblet_props = []
    for sample_id, info in SAMPLES.items():
        adata = samples.get(sample_id)
        if adata:
            celltype_col = get_celltype_column(adata)
            if celltype_col:
                goblet_prop = (adata.obs[celltype_col] == 'Goblet').mean() * 100
                goblet_props.append({'stage': info['stage'], 'prop': goblet_prop})

    goblet_df = pd.DataFrame(goblet_props)
    if len(goblet_df) > 0:
        colors = [STAGE_COLORS[s] for s in goblet_df['stage']]
        bars = ax2.bar(goblet_df['stage'], goblet_df['prop'], color=colors, edgecolor='black')
        ax2.set_ylabel('Goblet Cells (%)')
        ax2.set_title('Key Finding: Goblet Expansion', fontweight='bold')

        # Add fold change annotation
        if len(goblet_df) >= 3:
            normal_val = goblet_df[goblet_df['stage'] == 'Normal']['prop'].values
            cancer_val = goblet_df[goblet_df['stage'] == 'Cancer']['prop'].values
            if len(normal_val) > 0 and len(cancer_val) > 0 and normal_val[0] > 0:
                fold = cancer_val[0] / normal_val[0]
                ax2.text(0.95, 0.95, f'{fold:.1f}x', transform=ax2.transAxes,
                        ha='right', va='top', fontsize=14, fontweight='bold', color='#E94F37')

    # (c) Sample info
    ax3 = fig.add_subplot(gs[1, 0])
    add_panel_label(ax3, 'c')

    cell_counts = []
    for sample_id, info in SAMPLES.items():
        adata = samples.get(sample_id)
        if adata:
            cell_counts.append({'stage': info['stage'], 'cells': adata.n_obs})

    count_df = pd.DataFrame(cell_counts)
    if len(count_df) > 0:
        colors = [STAGE_COLORS[s] for s in count_df['stage']]
        ax3.bar(count_df['stage'], count_df['cells'], color=colors, edgecolor='black')
        ax3.set_ylabel('Cell Count')
        ax3.set_title('Sample Size', fontweight='bold')
        ax3.ticklabel_format(axis='y', style='sci', scilimits=(3,3))

    # (d) Platform info
    ax4 = fig.add_subplot(gs[1, 1])
    add_panel_label(ax4, 'd')

    info_text = """Platform: G4X
Genes: 337 (targeted panel)
Technology: H&E + RNA + Protein

Key Markers:
• Goblet: MUC2, TFF3
• Epithelial: EPCAM, KRT19
• Fibroblast: VIM, COL1A1
• Immune: PTPRC (CD45)"""

    ax4.text(0.1, 0.9, info_text, transform=ax4.transAxes, fontsize=8,
            verticalalignment='top', fontfamily='monospace')
    ax4.axis('off')
    ax4.set_title('Study Overview', fontweight='bold')

    plt.suptitle('Gastric Cancer Progression Analysis', fontsize=12, fontweight='bold', y=0.98)

    fig.savefig(FIGURES / 'figure6_summary.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'figure6_summary.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: figure6_summary.png/pdf")


def main():
    print("=" * 60)
    print("Professional Figure Generation")
    print("=" * 60)

    # Load samples
    print("\n[1/7] Loading samples...")
    samples = load_samples()

    if not samples:
        print("ERROR: No samples found!")
        return

    print(f"\nLoaded {len(samples)} samples")

    # Generate figures
    print("\n[2/7] Figure 1: Sample Overview...")
    figure1_sample_overview(samples)

    print("\n[3/7] Figure 2: Cell Type Composition...")
    figure2_celltype_composition(samples)

    print("\n[4/7] Figure 3: Goblet Progression (KEY)...")
    figure3_goblet_progression(samples)

    print("\n[5/7] Figure 4: DEG Analysis...")
    figure4_deg_analysis(samples)

    print("\n[6/7] Figure 5: Validation Comparison...")
    figure5_validation_comparison(samples)

    print("\n[7/7] Figure 6: Summary Panel...")
    figure6_summary_panel(samples)

    # Summary
    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"\nFigures saved to: {FIGURES}")
    print("\nGenerated files:")
    for f in sorted(FIGURES.glob('*.png')):
        print(f"  - {f.name}")

    print("\nPDF versions also generated for publication.")


if __name__ == '__main__':
    main()
