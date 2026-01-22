#!/usr/bin/env python3
"""Spatial visualization of annotated cell types."""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

SAMPLES = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}
INPUT = Path('/home/user/g4x-choi-batch2-analysis/results/pilot')
OUTPUT = INPUT / 'figures'
OUTPUT.mkdir(exist_ok=True)

COLORS = {
    'Epithelial': '#2ECC71', 'Pit_Mucous': '#1ABC9C', 'Goblet': '#117A65',
    'Intestinal_Meta': '#0E6655', 'T_Cell': '#3498DB', 'Macrophage': '#E67E22',
    'B_Cell': '#9B59B6', 'Fibroblast': '#E74C3C', 'CAF': '#C0392B',
    'Endothelial': '#EC7063'
}

def plot_spatial(adata, title, path):
    """Plot cells at spatial coordinates, colored by type."""
    fig, ax = plt.subplots(figsize=(10, 10), facecolor='white')

    for ct in adata.obs['cell_type'].unique():
        mask = adata.obs['cell_type'] == ct
        ax.scatter(adata.obs.loc[mask, 'cell_x'], adata.obs.loc[mask, 'cell_y'],
                  c=COLORS.get(ct, '#95A5A6'), s=0.5, alpha=0.7, label=ct, rasterized=True)

    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(f'{title}\n{adata.n_obs:,} cells', fontsize=14, fontweight='bold')

    handles = [mpatches.Patch(color=COLORS.get(ct, '#95A5A6'), label=ct)
               for ct in adata.obs['cell_type'].value_counts().index]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1, 1), fontsize=9)

    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()

def main():
    print("Creating spatial plots...")

    # Individual plots
    for sid, stage in SAMPLES.items():
        path = INPUT / f'{sid}_annotated.h5ad'
        if not path.exists():
            print(f"  {sid}: NOT FOUND - run 02_annotate.py first")
            continue

        adata = sc.read_h5ad(path)
        plot_spatial(adata, f'{sid}: {stage}', OUTPUT / f'{sid}_spatial.png')
        print(f"  {sid}: saved")

    # Combined progression plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), facecolor='white')
    for i, (sid, stage) in enumerate(SAMPLES.items()):
        path = INPUT / f'{sid}_annotated.h5ad'
        if not path.exists():
            continue

        adata = sc.read_h5ad(path)
        ax = axes[i]

        for ct in adata.obs['cell_type'].unique():
            mask = adata.obs['cell_type'] == ct
            ax.scatter(adata.obs.loc[mask, 'cell_x'], adata.obs.loc[mask, 'cell_y'],
                      c=COLORS.get(ct, '#95A5A6'), s=0.3, alpha=0.6, rasterized=True)

        ax.invert_yaxis()
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(f'{sid}: {stage}\n({adata.n_obs:,} cells)', fontsize=12, fontweight='bold')

    plt.suptitle('Gastric Cancer Progression: Normal → Metaplasia → Cancer', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT / 'progression.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {OUTPUT / 'progression.png'}")

if __name__ == '__main__':
    main()
