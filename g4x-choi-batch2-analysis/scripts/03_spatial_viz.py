#!/usr/bin/env python3
"""
03_spatial_viz.py - Spatial visualization of annotated cell types

Creates:
1. UMAP colored by cell type (per sample)
2. Spatial scatter plots colored by cell type
3. Cell type density heatmaps
4. 3-panel progression figure (E02 → F02 → G02)
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# === CONFIG ===
INPUT = Path('/home/user/g4x-choi-batch2-analysis/results/pilot')
OUTPUT = INPUT / 'figures'
OUTPUT.mkdir(exist_ok=True)

SAMPLES = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}

# Color palette for cell types
COLORS = {
    # Epithelial (greens)
    'Epithelial': '#2ECC71',
    'Epithelial_General': '#27AE60',
    'Gastric_Pit': '#1ABC9C',
    'Gastric_Chief': '#16A085',
    'Goblet': '#117A65',
    'Intestinal_Meta': '#0E6655',
    'Enteroendocrine': '#0B5345',
    'Stem_Progenitor': '#82E0AA',

    # T cells (blues)
    'T_Cell': '#3498DB',
    'T_CD4': '#5DADE2',
    'T_CD8_Cytotoxic': '#2980B9',
    'T_Reg': '#1F618D',
    'T_Exhausted': '#154360',

    # B/Plasma (purples)
    'B_Cell': '#9B59B6',
    'Plasma': '#7D3C98',

    # Myeloid (oranges)
    'Macrophage': '#E67E22',
    'Monocyte': '#D35400',

    # Stromal (reds)
    'Fibroblast': '#E74C3C',
    'CAF': '#C0392B',
    'Endothelial': '#EC7063',

    # Unknown
    'Unknown': '#95A5A6',
}

# Annotation method to use (priority order)
ANNOTATION_METHODS = ['celltype_markers', 'celltype_celltypist', 'celltype_gating']


def get_annotation_col(adata):
    """Get best available annotation column."""
    for method in ANNOTATION_METHODS:
        if method in adata.obs.columns:
            return method
    return None


def plot_umap(adata, title, path, annotation_col):
    """Plot UMAP colored by cell type."""
    fig, ax = plt.subplots(figsize=(10, 8), facecolor='white')

    cell_types = adata.obs[annotation_col].unique()

    for ct in sorted(cell_types):
        mask = adata.obs[annotation_col] == ct
        if mask.sum() > 0:
            ax.scatter(
                adata.obsm['X_umap'][mask, 0],
                adata.obsm['X_umap'][mask, 1],
                c=COLORS.get(ct, '#95A5A6'),
                s=3,
                alpha=0.7,
                label=f'{ct} ({mask.sum():,})',
                rasterized=True
            )

    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.set_title(f'{title}\n{adata.n_obs:,} cells', fontsize=14, fontweight='bold')
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=8, markerscale=3)
    ax.set_aspect('equal')

    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_spatial(adata, title, path, annotation_col):
    """Plot cells at spatial coordinates, colored by type."""
    fig, ax = plt.subplots(figsize=(12, 10), facecolor='white')

    # Get spatial coordinates
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    elif 'cell_x' in adata.obs.columns and 'cell_y' in adata.obs.columns:
        coords = adata.obs[['cell_x', 'cell_y']].values
    else:
        print(f"    WARNING: No spatial coordinates found for {title}")
        plt.close()
        return

    cell_types = adata.obs[annotation_col].unique()

    for ct in sorted(cell_types):
        mask = adata.obs[annotation_col] == ct
        if mask.sum() > 0:
            ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                c=COLORS.get(ct, '#95A5A6'),
                s=0.5,
                alpha=0.7,
                label=ct,
                rasterized=True
            )

    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(f'{title}\n{adata.n_obs:,} cells', fontsize=14, fontweight='bold')

    # Legend
    handles = [mpatches.Patch(color=COLORS.get(ct, '#95A5A6'), label=ct)
               for ct in adata.obs[annotation_col].value_counts().index[:15]]  # Top 15
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1, 1), fontsize=9)

    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_progression(samples_dict, output_path):
    """Create 3-panel progression figure."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 7), facecolor='white')

    all_cell_types = set()

    for i, (sample_id, stage) in enumerate(SAMPLES.items()):
        path = INPUT / f'{sample_id}_annotated.h5ad'
        if not path.exists():
            axes[i].text(0.5, 0.5, 'Data not found', ha='center', va='center')
            continue

        adata = sc.read_h5ad(path)
        annotation_col = get_annotation_col(adata)
        if annotation_col is None:
            axes[i].text(0.5, 0.5, 'No annotations', ha='center', va='center')
            continue

        ax = axes[i]

        # Get spatial coordinates
        if 'spatial' in adata.obsm:
            coords = adata.obsm['spatial']
        elif 'cell_x' in adata.obs.columns and 'cell_y' in adata.obs.columns:
            coords = adata.obs[['cell_x', 'cell_y']].values
        else:
            continue

        cell_types = adata.obs[annotation_col].unique()
        all_cell_types.update(cell_types)

        for ct in sorted(cell_types):
            mask = adata.obs[annotation_col] == ct
            if mask.sum() > 0:
                ax.scatter(
                    coords[mask, 0],
                    coords[mask, 1],
                    c=COLORS.get(ct, '#95A5A6'),
                    s=0.3,
                    alpha=0.6,
                    rasterized=True
                )

        ax.invert_yaxis()
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(f'{sample_id}: {stage}\n({adata.n_obs:,} cells)', fontsize=12, fontweight='bold')

    # Shared legend
    handles = [mpatches.Patch(color=COLORS.get(ct, '#95A5A6'), label=ct)
               for ct in sorted(all_cell_types) if ct != 'Unknown']
    fig.legend(handles=handles, loc='center right', bbox_to_anchor=(1.15, 0.5), fontsize=9)

    plt.suptitle('Gastric Cancer Progression: Normal → Metaplasia → Cancer',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_cell_type_proportions(output_path):
    """Bar plot of cell type proportions across stages."""
    data = []

    for sample_id, stage in SAMPLES.items():
        path = INPUT / f'{sample_id}_annotated.h5ad'
        if not path.exists():
            continue

        adata = sc.read_h5ad(path)
        annotation_col = get_annotation_col(adata)
        if annotation_col is None:
            continue

        counts = adata.obs[annotation_col].value_counts()
        proportions = counts / counts.sum()

        for ct, prop in proportions.items():
            data.append({
                'sample_id': sample_id,
                'stage': stage,
                'cell_type': ct,
                'proportion': prop,
                'count': counts[ct]
            })

    if not data:
        return

    df = pd.DataFrame(data)

    # Get top cell types (excluding Unknown)
    top_types = df[df['cell_type'] != 'Unknown'].groupby('cell_type')['proportion'].mean().nlargest(10).index

    fig, ax = plt.subplots(figsize=(12, 6))

    df_plot = df[df['cell_type'].isin(top_types)]
    pivot = df_plot.pivot(index='cell_type', columns='stage', values='proportion')

    # Reorder columns
    col_order = ['Normal', 'Metaplasia', 'Cancer']
    pivot = pivot[[c for c in col_order if c in pivot.columns]]

    pivot.plot(kind='bar', ax=ax, color=['#3498DB', '#F39C12', '#E74C3C'])
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Proportion')
    ax.set_title('Cell Type Proportions by Stage', fontweight='bold')
    ax.legend(title='Stage')
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close()


def main():
    print("="*60)
    print("Spatial Visualization")
    print("="*60)

    # Individual sample plots
    for sample_id, stage in SAMPLES.items():
        path = INPUT / f'{sample_id}_annotated.h5ad'
        if not path.exists():
            print(f"  {sample_id}: NOT FOUND - run 02_annotate_multimethod.py first")
            continue

        print(f"\n{sample_id} ({stage})")
        adata = sc.read_h5ad(path)

        annotation_col = get_annotation_col(adata)
        if annotation_col is None:
            print(f"  WARNING: No annotation columns found")
            continue

        print(f"  Using annotation: {annotation_col}")

        # UMAP plot
        if 'X_umap' in adata.obsm:
            plot_umap(adata, f'{sample_id}: {stage}', OUTPUT / f'{sample_id}_umap.png', annotation_col)
            print(f"  Saved: {sample_id}_umap.png")
        else:
            print(f"  WARNING: No UMAP found")

        # Spatial plot
        plot_spatial(adata, f'{sample_id}: {stage}', OUTPUT / f'{sample_id}_spatial.png', annotation_col)
        print(f"  Saved: {sample_id}_spatial.png")

    # Combined progression plot
    print("\nCreating progression figure...")
    plot_progression(SAMPLES, OUTPUT / 'progression_3panel.png')
    print(f"Saved: progression_3panel.png")

    # Cell type proportions
    print("Creating proportion plot...")
    plot_cell_type_proportions(OUTPUT / 'celltype_proportions.png')
    print(f"Saved: celltype_proportions.png")

    print(f"\nAll figures saved to: {OUTPUT}")


if __name__ == '__main__':
    main()
