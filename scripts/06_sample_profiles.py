#!/usr/bin/env python3
"""
G4X Choi Batch 2 - Sample Profile Dashboards
=============================================

Generates comprehensive 8-panel visual profiles for each of 32 samples.
"""

import os
import sys
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from pathlib import Path
from PIL import Image
from scipy import stats
from scipy.ndimage import gaussian_filter
import warnings
warnings.filterwarnings('ignore')

# Configuration
DATA_ROOT = Path('/mnt/x/Choi_Batch_2_Tuesday')
OUTPUT_DIR = Path('/home/user/g4x-choi-batch2-analysis/results')
FIG_DIR = OUTPUT_DIR / 'figures' / 'sample_profiles'
TABLE_DIR = OUTPUT_DIR / 'tables'

# Create output directories
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Sample layout
LANES = ['L001', 'L002', 'L003', 'L004']
SAMPLES_PER_LANE = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Protein markers
PROTEIN_MARKERS = [
    'CD3', 'CD4', 'CD8', 'FOXP3', 'CD20', 'CD68', 'CD11c', 'HLA-DR',
    'PD1', 'PDL1', 'PanCK', 'aSMA', 'CD31', 'CD45', 'KI67'
]

# Cell type colors
CELL_TYPE_COLORS = {
    'Epithelial': '#E41A1C',
    'CD8_T': '#377EB8',
    'CD4_T': '#4DAF4A',
    'T_cell': '#8DD3C7',  # Generic T cells
    'Treg': '#984EA3',
    'B_cell': '#FF7F00',
    'Macrophage': '#A65628',
    'DC': '#F781BF',
    'Stromal': '#999999',
    'Endothelial': '#66C2A5',
    'Immune_other': '#FFED6F',
    'Other': '#CCCCCC'
}

def get_sample_path(lane_idx, sample_letter):
    """Get path to sample directory."""
    lane = f'L00{lane_idx + 1}'
    sample = f'{sample_letter}0{lane_idx + 1}'
    lane_dirs = list(DATA_ROOT.glob(f'g4-028-083-FC1-{lane}_*'))
    if not lane_dirs:
        return None
    return lane_dirs[0] / sample

def load_cell_metadata(sample_path):
    """Load cell metadata with positions."""
    metadata_file = sample_path / 'single_cell_data' / 'cell_metadata.csv.gz'
    if not metadata_file.exists():
        return None
    with gzip.open(metadata_file, 'rt') as f:
        df = pd.read_csv(f)
    return df

def load_protein_data(sample_path):
    """Load protein expression data."""
    protein_file = sample_path / 'single_cell_data' / 'cell_by_protein.csv.gz'
    if not protein_file.exists():
        return None
    with gzip.open(protein_file, 'rt') as f:
        df = pd.read_csv(f)
    # Replace Inf with NaN, then fill NaN with 0 for intensity columns
    intensity_cols = [c for c in df.columns if '_intensity_mean' in c]
    for col in intensity_cols:
        df[col] = df[col].replace([np.inf, -np.inf], np.nan).fillna(0)
    return df

def load_he_thumbnail(sample_path):
    """Load H&E thumbnail image."""
    he_thumb = sample_path / 'h_and_e' / 'h_and_e_thumbnail.jpg'
    if not he_thumb.exists():
        return None
    return Image.open(he_thumb)

def assign_cell_types(protein_df):
    """Simple cell type assignment based on protein expression."""
    cell_types = pd.Series('Other', index=protein_df.index)

    # Get column names - handle both formats (with and without _intensity_mean suffix)
    cols = protein_df.columns.tolist()

    # Map marker names to actual column names
    def get_col(marker):
        if marker in cols:
            return marker
        elif f'{marker}_intensity_mean' in cols:
            return f'{marker}_intensity_mean'
        return None

    # Helper to check marker positivity (top 30%)
    def is_positive(marker):
        col = get_col(marker)
        if col is None:
            return pd.Series(False, index=protein_df.index)
        threshold = protein_df[col].quantile(0.7)
        return protein_df[col] > threshold

    # Assign in order of specificity (most specific first)
    # Epithelial
    if get_col('PanCK'):
        cell_types[is_positive('PanCK')] = 'Epithelial'

    # Stromal
    if get_col('aSMA'):
        cell_types[is_positive('aSMA') & ~is_positive('PanCK')] = 'Stromal'

    # Endothelial
    if get_col('CD31'):
        cell_types[is_positive('CD31')] = 'Endothelial'

    # Immune cells (CD45+)
    if get_col('CD45'):
        immune_mask = is_positive('CD45')

        # Macrophages
        if get_col('CD68'):
            cell_types[immune_mask & is_positive('CD68')] = 'Macrophage'

        # B cells
        if get_col('CD20'):
            cell_types[immune_mask & is_positive('CD20')] = 'B_cell'

        # T cells (CD3+)
        if get_col('CD3'):
            t_cell_mask = immune_mask & is_positive('CD3')
            cell_types[t_cell_mask] = 'T_cell'

            # CD8 T cells
            if get_col('CD8'):
                cell_types[t_cell_mask & is_positive('CD8')] = 'CD8_T'

            # CD4 T cells
            if get_col('CD4'):
                cd4_mask = t_cell_mask & is_positive('CD4')
                cell_types[cd4_mask] = 'CD4_T'

                # Tregs
                if get_col('FOXP3'):
                    cell_types[cd4_mask & is_positive('FOXP3')] = 'Treg'

        # DCs
        if get_col('CD11c') and get_col('HLA-DR'):
            dc_mask = immune_mask & is_positive('CD11c') & is_positive('HLA-DR')
            cell_types[dc_mask] = 'DC'

    return cell_types

def generate_sample_profile(sample_id, lane_idx, sample_letter):
    """Generate comprehensive 8-panel profile for one sample."""

    sample_path = get_sample_path(lane_idx, sample_letter)
    if sample_path is None or not sample_path.exists():
        print(f"  Sample path not found for {sample_id}")
        return False

    # Load data
    metadata = load_cell_metadata(sample_path)
    protein_df = load_protein_data(sample_path)
    he_image = load_he_thumbnail(sample_path)

    if metadata is None or protein_df is None:
        print(f"  Data not found for {sample_id}")
        return False

    # Assign cell types
    cell_types = assign_cell_types(protein_df)
    metadata['cell_type'] = cell_types.values

    # Create figure with 8 panels (4x2)
    fig = plt.figure(figsize=(24, 28))
    gs = GridSpec(4, 2, figure=fig, hspace=0.3, wspace=0.25)

    # --- Panel 1: Spatial Map by Cell Type ---
    ax1 = fig.add_subplot(gs[0, 0])

    # Check for position columns
    if 'cell_x' in metadata.columns and 'cell_y' in metadata.columns:
        x_col, y_col = 'cell_x', 'cell_y'
    elif 'center_x_um' in metadata.columns and 'center_y_um' in metadata.columns:
        x_col, y_col = 'center_x_um', 'center_y_um'
    elif 'CenterX_global_px' in metadata.columns:
        x_col, y_col = 'CenterX_global_px', 'CenterY_global_px'
    else:
        x_col, y_col = metadata.columns[2], metadata.columns[3]  # fallback

    # Subsample for plotting if too many cells
    n_cells = len(metadata)
    if n_cells > 50000:
        plot_idx = np.random.choice(n_cells, 50000, replace=False)
        plot_df = metadata.iloc[plot_idx]
    else:
        plot_df = metadata

    # Plot by cell type
    for ct, color in CELL_TYPE_COLORS.items():
        mask = plot_df['cell_type'] == ct
        if mask.sum() > 0:
            ax1.scatter(
                plot_df.loc[mask, x_col], plot_df.loc[mask, y_col],
                c=color, s=0.5, alpha=0.6, label=ct, rasterized=True
            )

    ax1.set_title(f'Spatial Distribution - {sample_id} (n={n_cells:,})', fontsize=14, fontweight='bold')
    ax1.set_xlabel('X Position (μm)')
    ax1.set_ylabel('Y Position (μm)')
    ax1.legend(loc='upper right', fontsize=8, markerscale=10)
    ax1.set_aspect('equal')

    # --- Panel 2: H&E Thumbnail ---
    ax2 = fig.add_subplot(gs[0, 1])
    if he_image is not None:
        ax2.imshow(he_image)
        ax2.set_title('H&E Image', fontsize=14, fontweight='bold')
    else:
        ax2.text(0.5, 0.5, 'H&E not available', ha='center', va='center', transform=ax2.transAxes)
    ax2.axis('off')

    # --- Panel 3: Protein Expression Heatmap ---
    ax3 = fig.add_subplot(gs[1, 0])

    # Get mean expression by cell type - handle column name formats
    def get_marker_col(m):
        if m in protein_df.columns:
            return m
        elif f'{m}_intensity_mean' in protein_df.columns:
            return f'{m}_intensity_mean'
        return None

    available_markers = []
    marker_cols = []
    for m in PROTEIN_MARKERS:
        col = get_marker_col(m)
        if col:
            available_markers.append(m)
            marker_cols.append(col)

    protein_df['cell_type'] = cell_types.values

    if len(marker_cols) > 0:
        expr_by_type = protein_df.groupby('cell_type')[marker_cols].mean()
        # Rename columns to marker names
        expr_by_type.columns = available_markers
        expr_by_type = expr_by_type.loc[expr_by_type.index.isin(CELL_TYPE_COLORS.keys())]

        if len(expr_by_type) > 0 and len(expr_by_type.columns) > 0:
            # Z-score normalize, handling NaN/Inf
            expr_by_type = expr_by_type.replace([np.inf, -np.inf], np.nan)
            expr_z = (expr_by_type - expr_by_type.mean()) / (expr_by_type.std() + 1e-10)
            expr_z = expr_z.fillna(0)  # Fill remaining NaN with 0
            sns.heatmap(
                expr_z.T, ax=ax3, cmap='RdBu_r', center=0,
                cbar_kws={'label': 'Z-score'},
                xticklabels=True, yticklabels=True
            )
            ax3.set_title('Protein Expression by Cell Type', fontsize=14, fontweight='bold')
            ax3.set_xlabel('Cell Type')
            ax3.set_ylabel('Marker')
        else:
            ax3.text(0.5, 0.5, 'No cell types detected', ha='center', va='center', transform=ax3.transAxes)
            ax3.set_title('Protein Expression by Cell Type', fontsize=14, fontweight='bold')
    else:
        ax3.text(0.5, 0.5, 'No markers available', ha='center', va='center', transform=ax3.transAxes)
        ax3.set_title('Protein Expression by Cell Type', fontsize=14, fontweight='bold')

    # --- Panel 4: Cell Type Composition ---
    ax4 = fig.add_subplot(gs[1, 1])

    type_counts = metadata['cell_type'].value_counts()
    # Sort by standard order
    ordered_types = [t for t in CELL_TYPE_COLORS.keys() if t in type_counts.index]
    type_counts = type_counts[ordered_types]
    colors = [CELL_TYPE_COLORS[t] for t in ordered_types]

    wedges, texts, autotexts = ax4.pie(
        type_counts.values,
        labels=ordered_types,
        colors=colors,
        autopct=lambda p: f'{p:.1f}%' if p > 3 else '',
        startangle=90
    )
    ax4.set_title('Cell Type Composition', fontsize=14, fontweight='bold')

    # --- Panel 5: Marker Violin Plots ---
    ax5 = fig.add_subplot(gs[2, 0])

    # Select key markers for visualization
    key_marker_names = ['CD3', 'CD8', 'CD4', 'FOXP3', 'PD1', 'PDL1', 'PanCK', 'CD68']
    key_markers = []
    key_cols = []
    for m in key_marker_names:
        col = get_marker_col(m)
        if col:
            key_markers.append(m)
            key_cols.append(col)

    if len(key_cols) > 0:
        # Subsample for violin plot
        if len(protein_df) > 10000:
            violin_df = protein_df.sample(10000)
        else:
            violin_df = protein_df

        # Create violin data with proper column names
        violin_subset = violin_df[key_cols].copy()
        violin_subset.columns = key_markers
        # Remove NaN/Inf values
        violin_subset = violin_subset.replace([np.inf, -np.inf], np.nan).dropna()
        if len(violin_subset) > 0:
            violin_data = violin_subset.melt(var_name='Marker', value_name='Expression')
            sns.violinplot(data=violin_data, x='Marker', y='Expression', ax=ax5,
                          palette='Set2', inner='quartile', cut=0)
            ax5.set_title('Key Marker Distributions', fontsize=14, fontweight='bold')
            ax5.set_xticklabels(ax5.get_xticklabels(), rotation=45, ha='right')
            ax5.set_ylabel('Expression')
        else:
            ax5.text(0.5, 0.5, 'Data contains NaN values', ha='center', va='center', transform=ax5.transAxes)
            ax5.set_title('Key Marker Distributions', fontsize=14, fontweight='bold')

    # --- Panel 6: Spatial Density (KDE) ---
    ax6 = fig.add_subplot(gs[2, 1])

    # Focus on immune cells
    immune_types = ['CD8_T', 'CD4_T', 'Treg', 'Macrophage', 'B_cell', 'DC']
    immune_mask = metadata['cell_type'].isin(immune_types)

    if immune_mask.sum() > 100:
        x = metadata.loc[immune_mask, x_col].values
        y = metadata.loc[immune_mask, y_col].values

        # Create 2D histogram for density
        x_range = (x.min(), x.max())
        y_range = (y.min(), y.max())

        h, xedges, yedges = np.histogram2d(x, y, bins=50, range=[x_range, y_range])
        h_smooth = gaussian_filter(h.T, sigma=1.5)

        im = ax6.imshow(
            h_smooth, origin='lower',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
            cmap='YlOrRd', aspect='auto'
        )
        plt.colorbar(im, ax=ax6, label='Immune Cell Density')
        ax6.set_title('Immune Cell Spatial Density', fontsize=14, fontweight='bold')
        ax6.set_xlabel('X Position (μm)')
        ax6.set_ylabel('Y Position (μm)')
    else:
        ax6.text(0.5, 0.5, 'Insufficient immune cells', ha='center', va='center', transform=ax6.transAxes)

    # --- Panel 7: Cell Type Proportions Bar Chart ---
    ax7 = fig.add_subplot(gs[3, 0])

    # Calculate proportions for specific types
    immune_types_detail = ['CD8_T', 'CD4_T', 'Treg', 'B_cell', 'Macrophage', 'DC']
    proportions = []
    labels = []
    bar_colors = []

    for ct in immune_types_detail:
        count = (metadata['cell_type'] == ct).sum()
        prop = count / len(metadata) * 100
        proportions.append(prop)
        labels.append(ct)
        bar_colors.append(CELL_TYPE_COLORS.get(ct, '#CCCCCC'))

    bars = ax7.barh(labels, proportions, color=bar_colors)
    ax7.set_xlabel('Percentage of Total Cells')
    ax7.set_title('Immune Cell Proportions', fontsize=14, fontweight='bold')

    # Add value labels
    for bar, prop in zip(bars, proportions):
        ax7.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
                f'{prop:.1f}%', va='center', fontsize=10)

    ax7.set_xlim(0, max(proportions) * 1.3 if proportions else 10)

    # --- Panel 8: Checkpoint Co-expression ---
    ax8 = fig.add_subplot(gs[3, 1])

    pd1_col = get_marker_col('PD1')
    pdl1_col = get_marker_col('PDL1')

    if pd1_col and pdl1_col:
        # Subsample for scatter
        if len(protein_df) > 5000:
            scatter_df = protein_df.sample(5000)
        else:
            scatter_df = protein_df

        ax8.scatter(
            scatter_df[pd1_col], scatter_df[pdl1_col],
            c=scatter_df['cell_type'].map(CELL_TYPE_COLORS).fillna('#CCCCCC'),
            s=3, alpha=0.5, rasterized=True
        )

        # Add correlation
        r, p = stats.pearsonr(protein_df[pd1_col], protein_df[pdl1_col])
        ax8.text(0.05, 0.95, f'r = {r:.3f}\np = {p:.2e}', transform=ax8.transAxes,
                fontsize=10, va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax8.set_xlabel('PD1 Expression')
        ax8.set_ylabel('PDL1 Expression')
        ax8.set_title('Checkpoint Co-expression', fontsize=14, fontweight='bold')
    else:
        ax8.text(0.5, 0.5, 'PD1/PDL1 not available', ha='center', va='center', transform=ax8.transAxes)

    # Overall title
    fig.suptitle(f'Sample Profile: {sample_id}', fontsize=18, fontweight='bold', y=0.98)

    # Save
    output_path = FIG_DIR / f'{sample_id}_profile.png'
    fig.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    print(f"  ✓ Generated profile for {sample_id}")
    return True

def main():
    """Generate profiles for all 32 samples."""
    print("="*60)
    print("G4X Choi Batch 2 - Sample Profile Generation")
    print("="*60)

    successful = 0
    failed = 0

    for lane_idx in range(4):
        lane_num = lane_idx + 1
        print(f"\nProcessing Lane {lane_num}...")

        for sample_letter in SAMPLES_PER_LANE:
            sample_id = f'{sample_letter}0{lane_num}'
            print(f"  Processing {sample_id}...")

            try:
                if generate_sample_profile(sample_id, lane_idx, sample_letter):
                    successful += 1
                else:
                    failed += 1
            except Exception as e:
                print(f"  ERROR: {e}")
                failed += 1

    print("\n" + "="*60)
    print(f"COMPLETE: {successful} successful, {failed} failed")
    print(f"Output: {FIG_DIR}")
    print("="*60)

if __name__ == '__main__':
    main()
