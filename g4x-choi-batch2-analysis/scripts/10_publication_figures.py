#!/usr/bin/env python3
"""
Publication Figures for G4X RCTD + IM Analysis

Generates integrated figures showing:
1. RCTD cell type spatial distribution
2. IM score spatial patterns
3. Hwang cancer score progression
4. QC metrics panel

Critical fixes:
- NaN coordinate masking
- Reproducible random seed
- Proper weights alignment validation
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

# CRITICAL: Reproducible random seed
RNG = np.random.default_rng(42)

# ============================================================================
# Configuration
# ============================================================================

CONFIG = {
    'data_path': 'output/pilot_final_annotated.h5ad',
    'output_dir': 'output/figures',
    'dpi': 300,
    'figsize_panel': (16, 12),
    'figsize_single': (8, 6),
}

# Color palettes
CELL_TYPE_COLORS = {
    'Cancer': '#e41a1c',
    'PMC': '#377eb8',
    'Neck_like': '#4daf4a',
    'Enterocyte': '#984ea3',
    'Goblet': '#ff7f00',
    'B_cell': '#ffff33',
    'T_cell': '#a65628',
    'Fibroblast': '#f781bf',
    'Endothelial': '#999999',
    'MSC': '#66c2a5',
    'Macrophage': '#fc8d62',
    'Mast': '#8da0cb',
    'Enteroendocrine': '#e78ac3',
    'Chief': '#a6d854',
    'GMC': '#ffd92f',
    'PC': '#e5c494',
    'Smooth_Muscle': '#b3b3b3',
    'Unknown': '#cccccc',
}

STAGE_COLORS = {
    'Normal': '#2166ac',
    'Metaplasia': '#92c5de',
    'Cancer': '#b2182b',
}

# ============================================================================
# Helper Functions
# ============================================================================

def get_valid_coords(sample_data):
    """
    Extract coordinates with NaN masking.

    Returns:
        coords: (N, 2) array of coordinates
        valid_mask: boolean array indicating valid (non-NaN) coordinates
    """
    # Try obsm['spatial'] first
    if 'spatial' in sample_data.obsm:
        coords = sample_data.obsm['spatial']
    # Fall back to obs columns
    elif 'cell_x' in sample_data.obs and 'cell_y' in sample_data.obs:
        coords = np.column_stack([
            sample_data.obs['cell_x'].values,
            sample_data.obs['cell_y'].values
        ])
    else:
        return None, None

    # CRITICAL: Mask NaN coordinates
    valid_mask = ~np.isnan(coords).any(axis=1)

    return coords, valid_mask


def plot_spatial_category(ax, coords, categories, valid_mask, color_map,
                         title='', point_size=1, alpha=0.7):
    """Plot spatial scatter with categorical coloring."""
    # Apply valid mask
    coords_valid = coords[valid_mask]
    cats_valid = categories[valid_mask]

    # Get unique categories
    unique_cats = [c for c in cats_valid.unique() if pd.notna(c)]

    # Shuffle plotting order for better visibility
    shuffle_idx = RNG.permutation(len(coords_valid))

    for cat in unique_cats:
        cat_mask = cats_valid.values == cat
        cat_idx = shuffle_idx[cat_mask[shuffle_idx]]

        color = color_map.get(cat, '#cccccc')
        ax.scatter(
            coords_valid[cat_idx, 0],
            coords_valid[cat_idx, 1],
            c=color,
            s=point_size,
            alpha=alpha,
            label=cat,
            rasterized=True
        )

    ax.set_title(title)
    ax.set_aspect('equal')
    ax.axis('off')


def plot_spatial_continuous(ax, coords, values, valid_mask, title='',
                           cmap='viridis', point_size=1, alpha=0.7,
                           vmin=None, vmax=None):
    """Plot spatial scatter with continuous coloring."""
    # Apply valid mask and remove NaN values
    coords_valid = coords[valid_mask]
    vals_valid = values[valid_mask]

    # Further mask for NaN values
    val_mask = ~np.isnan(vals_valid)
    coords_plot = coords_valid[val_mask]
    vals_plot = vals_valid[val_mask]

    # Shuffle for better overlap handling
    shuffle_idx = RNG.permutation(len(coords_plot))

    scatter = ax.scatter(
        coords_plot[shuffle_idx, 0],
        coords_plot[shuffle_idx, 1],
        c=vals_plot[shuffle_idx],
        cmap=cmap,
        s=point_size,
        alpha=alpha,
        vmin=vmin,
        vmax=vmax,
        rasterized=True
    )

    ax.set_title(title)
    ax.set_aspect('equal')
    ax.axis('off')

    return scatter


# ============================================================================
# Figure Generation Functions
# ============================================================================

def generate_integrated_panel(adata, output_dir):
    """
    Generate main integrated figure panel.

    Layout:
    Row 1: RCTD cell types per sample (E02, F02, G02)
    Row 2: Hwang cancer score per sample
    Row 3: IM ratio per sample
    """
    print("\n--- Generating Integrated Panel ---")

    samples = ['E02', 'F02', 'G02']
    stages = ['Normal', 'Metaplasia', 'Cancer']

    fig, axes = plt.subplots(3, 3, figsize=(18, 16))

    for col_idx, (sample, stage) in enumerate(zip(samples, stages)):
        # Get sample data
        sample_mask = adata.obs['sample_id'] == sample
        sample_data = adata[sample_mask].copy()

        coords, valid_mask = get_valid_coords(sample_data)

        if coords is None:
            print(f"  Warning: No coordinates for {sample}")
            continue

        n_valid = valid_mask.sum()
        n_total = len(valid_mask)
        print(f"  {sample}: {n_valid}/{n_total} valid coordinates")

        # Row 1: RCTD cell types
        ax = axes[0, col_idx]
        if 'cell_type_rctd' in sample_data.obs:
            plot_spatial_category(
                ax, coords, sample_data.obs['cell_type_rctd'], valid_mask,
                CELL_TYPE_COLORS,
                title=f'{sample} ({stage}) - Cell Types',
                point_size=0.5
            )
        else:
            ax.text(0.5, 0.5, 'No RCTD data', ha='center', va='center',
                   transform=ax.transAxes)
            ax.axis('off')

        # Row 2: Hwang cancer score
        ax = axes[1, col_idx]
        if 'hwang_cancer' in sample_data.obs:
            scatter = plot_spatial_continuous(
                ax, coords, sample_data.obs['hwang_cancer'].values, valid_mask,
                title=f'{sample} - Hwang Cancer Score',
                cmap='RdBu_r',
                point_size=0.5,
                vmin=-0.5, vmax=0.5
            )
            plt.colorbar(scatter, ax=ax, shrink=0.5, label='Score')
        else:
            ax.text(0.5, 0.5, 'No Hwang data', ha='center', va='center',
                   transform=ax.transAxes)
            ax.axis('off')

        # Row 3: IM ratio
        ax = axes[2, col_idx]
        if 'im_ratio' in sample_data.obs:
            scatter = plot_spatial_continuous(
                ax, coords, sample_data.obs['im_ratio'].values, valid_mask,
                title=f'{sample} - IM Ratio (Intestinal/Total)',
                cmap='PiYG',
                point_size=0.5,
                vmin=0.3, vmax=0.7
            )
            plt.colorbar(scatter, ax=ax, shrink=0.5, label='Ratio')
        else:
            ax.text(0.5, 0.5, 'No IM data', ha='center', va='center',
                   transform=ax.transAxes)
            ax.axis('off')

    # Add legend for cell types (top right)
    if 'cell_type_rctd' in adata.obs:
        unique_types = adata.obs['cell_type_rctd'].dropna().unique()
        legend_patches = [
            mpatches.Patch(color=CELL_TYPE_COLORS.get(ct, '#cccccc'), label=ct)
            for ct in sorted(unique_types) if ct in CELL_TYPE_COLORS
        ]
        fig.legend(handles=legend_patches, loc='upper right', ncol=2,
                  fontsize=8, title='Cell Types')

    plt.suptitle('G4X Gastric Cancer Progression: RCTD + IM Analysis',
                fontsize=14, y=0.98)
    plt.tight_layout(rect=[0, 0, 0.9, 0.96])

    # Save
    for fmt in ['png', 'pdf']:
        path = output_dir / f'fig_integrated_panel.{fmt}'
        fig.savefig(path, dpi=CONFIG['dpi'], bbox_inches='tight')
        print(f"  Saved: {path}")

    plt.close(fig)


def generate_qc_panel(adata, output_dir):
    """
    Generate QC metrics panel.

    Shows:
    - Cell type distribution per stage
    - IM score distributions
    - Hwang score distributions
    - RCTD rejection rates
    """
    print("\n--- Generating QC Panel ---")

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Panel 1: Cell type distribution stacked bar
    ax = axes[0, 0]
    if 'cell_type_rctd' in adata.obs:
        ct_counts = pd.crosstab(
            adata.obs['stage'],
            adata.obs['cell_type_rctd'],
            normalize='index'
        )
        # Keep top 8 cell types, group rest as "Other"
        top_types = ct_counts.sum().nlargest(8).index
        ct_plot = ct_counts[top_types].copy()
        ct_plot['Other'] = ct_counts.drop(columns=top_types, errors='ignore').sum(axis=1)

        ct_plot.plot(kind='bar', stacked=True, ax=ax,
                    color=[CELL_TYPE_COLORS.get(ct, '#cccccc') for ct in ct_plot.columns])
        ax.set_title('Cell Type Distribution by Stage')
        ax.set_ylabel('Proportion')
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    else:
        ax.text(0.5, 0.5, 'No RCTD data', ha='center', va='center')

    # Panel 2: IM intestinal score violin
    ax = axes[0, 1]
    if 'im_intestinal' in adata.obs:
        data = [adata.obs[adata.obs['stage'] == s]['im_intestinal'].dropna()
               for s in ['Normal', 'Metaplasia', 'Cancer']]
        parts = ax.violinplot(data, positions=[1, 2, 3], showmeans=True)
        for pc, color in zip(parts['bodies'], [STAGE_COLORS['Normal'],
                                                STAGE_COLORS['Metaplasia'],
                                                STAGE_COLORS['Cancer']]):
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(['Normal', 'Metaplasia', 'Cancer'])
        ax.set_title('IM Intestinal Score')
        ax.set_ylabel('Score')
    else:
        ax.text(0.5, 0.5, 'No IM data', ha='center', va='center')

    # Panel 3: Hwang score violin
    ax = axes[0, 2]
    if 'hwang_cancer' in adata.obs:
        data = [adata.obs[adata.obs['stage'] == s]['hwang_cancer'].dropna()
               for s in ['Normal', 'Metaplasia', 'Cancer']]
        parts = ax.violinplot(data, positions=[1, 2, 3], showmeans=True)
        for pc, color in zip(parts['bodies'], [STAGE_COLORS['Normal'],
                                                STAGE_COLORS['Metaplasia'],
                                                STAGE_COLORS['Cancer']]):
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(['Normal', 'Metaplasia', 'Cancer'])
        ax.set_title('Hwang Cancer Score')
        ax.set_ylabel('Score')
    else:
        ax.text(0.5, 0.5, 'No Hwang data', ha='center', va='center')

    # Panel 4: RCTD spot class distribution
    ax = axes[1, 0]
    if 'spot_class' in adata.obs:
        spot_counts = pd.crosstab(
            adata.obs['stage'],
            adata.obs['spot_class'],
            normalize='index'
        ) * 100
        spot_counts.plot(kind='bar', ax=ax, width=0.8)
        ax.set_title('RCTD Spot Classification')
        ax.set_ylabel('Percentage')
        ax.legend(title='Spot Class', fontsize=8)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    else:
        ax.text(0.5, 0.5, 'No spot class data', ha='center', va='center')

    # Panel 5: IM ratio histogram by stage
    ax = axes[1, 1]
    if 'im_ratio' in adata.obs:
        for stage, color in STAGE_COLORS.items():
            data = adata.obs[adata.obs['stage'] == stage]['im_ratio'].dropna()
            ax.hist(data, bins=30, alpha=0.5, label=stage, color=color, density=True)
        ax.set_xlabel('IM Ratio')
        ax.set_ylabel('Density')
        ax.set_title('IM Ratio Distribution')
        ax.legend()
    else:
        ax.text(0.5, 0.5, 'No IM ratio data', ha='center', va='center')

    # Panel 6: Summary stats table
    ax = axes[1, 2]
    ax.axis('off')

    stats_text = "Summary Statistics\n" + "=" * 30 + "\n\n"
    stats_text += f"Total cells: {adata.n_obs:,}\n"
    stats_text += f"Genes: {adata.n_vars}\n\n"

    for sample in ['E02', 'F02', 'G02']:
        n = (adata.obs['sample_id'] == sample).sum()
        stage = adata.obs[adata.obs['sample_id'] == sample]['stage'].iloc[0]
        stats_text += f"{sample} ({stage}): {n:,} cells\n"

    if 'cell_type_rctd' in adata.obs:
        stats_text += f"\nRCTD assigned: {adata.obs['cell_type_rctd'].notna().sum():,}\n"

    if 'spot_class' in adata.obs:
        reject_rate = (adata.obs['spot_class'] == 'reject').mean() * 100
        stats_text += f"RCTD reject rate: {reject_rate:.1f}%\n"

    ax.text(0.1, 0.9, stats_text, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', fontfamily='monospace')

    plt.suptitle('G4X Analysis QC Panel', fontsize=14, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save
    for fmt in ['png', 'pdf']:
        path = output_dir / f'fig_qc_panel.{fmt}'
        fig.savefig(path, dpi=CONFIG['dpi'], bbox_inches='tight')
        print(f"  Saved: {path}")

    plt.close(fig)


def generate_cell_type_focus(adata, output_dir):
    """
    Generate focused cell type figures.

    Shows individual cell types across all samples.
    """
    print("\n--- Generating Cell Type Focus Figures ---")

    if 'cell_type_rctd' not in adata.obs:
        print("  No RCTD data, skipping")
        return

    # Get top cell types
    top_types = adata.obs['cell_type_rctd'].value_counts().nlargest(6).index.tolist()

    samples = ['E02', 'F02', 'G02']

    fig, axes = plt.subplots(len(top_types), 3, figsize=(15, 3 * len(top_types)))

    for row_idx, ct in enumerate(top_types):
        for col_idx, sample in enumerate(samples):
            ax = axes[row_idx, col_idx]

            # Get sample data
            sample_mask = adata.obs['sample_id'] == sample
            sample_data = adata[sample_mask].copy()

            coords, valid_mask = get_valid_coords(sample_data)

            if coords is None:
                ax.axis('off')
                continue

            # Highlight this cell type
            ct_mask = (sample_data.obs['cell_type_rctd'] == ct).values & valid_mask

            # Plot background (gray)
            bg_mask = valid_mask & ~ct_mask
            if bg_mask.sum() > 0:
                ax.scatter(
                    coords[bg_mask, 0], coords[bg_mask, 1],
                    c='#eeeeee', s=0.3, alpha=0.3, rasterized=True
                )

            # Plot highlighted cells
            if ct_mask.sum() > 0:
                ax.scatter(
                    coords[ct_mask, 0], coords[ct_mask, 1],
                    c=CELL_TYPE_COLORS.get(ct, 'red'), s=1, alpha=0.8,
                    rasterized=True
                )

            n_cells = ct_mask.sum()
            stage = adata.obs[adata.obs['sample_id'] == sample]['stage'].iloc[0]
            ax.set_title(f'{sample} ({stage})\nn={n_cells:,}', fontsize=10)
            ax.set_aspect('equal')
            ax.axis('off')

        # Row label
        axes[row_idx, 0].set_ylabel(ct, fontsize=12, rotation=0, ha='right',
                                    labelpad=50)

    plt.suptitle('Cell Type Spatial Distribution', fontsize=14, y=0.99)
    plt.tight_layout(rect=[0.05, 0, 1, 0.97])

    # Save
    for fmt in ['png', 'pdf']:
        path = output_dir / f'fig_celltype_focus.{fmt}'
        fig.savefig(path, dpi=CONFIG['dpi'], bbox_inches='tight')
        print(f"  Saved: {path}")

    plt.close(fig)


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 60)
    print("Publication Figures Generation")
    print("=" * 60)

    # Setup output directory
    output_dir = Path(CONFIG['output_dir'])
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print("\n--- Loading Data ---")
    data_path = Path(CONFIG['data_path'])

    if not data_path.exists():
        print(f"Error: Data not found at {data_path}")
        print("Run 09_merge_results.py first!")
        return

    adata = sc.read_h5ad(data_path)
    print(f"Shape: {adata.shape}")
    print(f"Obs columns: {list(adata.obs.columns)}")

    # Generate figures
    generate_integrated_panel(adata, output_dir)
    generate_qc_panel(adata, output_dir)
    generate_cell_type_focus(adata, output_dir)

    print("\n" + "=" * 60)
    print("Figure Generation Complete")
    print("=" * 60)
    print(f"\nOutput directory: {output_dir}")
    print(f"Files generated:")
    for f in sorted(output_dir.glob('fig_*.png')):
        print(f"  - {f.name}")


if __name__ == '__main__':
    main()
