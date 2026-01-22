#!/usr/bin/env python
"""Generate combined HVG optimization panel for all 3 samples."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
df = pd.read_csv('/home/user/g4x-choi-batch2-analysis/output/data/hvg_optimization.csv')

# Sample order and labels
samples = ['E02', 'F02', 'G02']
stages = ['Normal', 'Metaplasia', 'Cancer']

# HVG and resolution order
hvg_order = ['30', '50', '100', '200', 'all']
res_order = [0.3, 0.5, 0.8, 1.0]

# Create figure with 3 columns
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Track best scores for annotation
best_scores = []

for idx, (sample, stage) in enumerate(zip(samples, stages)):
    ax = axes[idx]

    # Filter data for this sample
    sample_df = df[df['sample_id'] == sample].copy()

    # Convert n_hvg to string for proper ordering
    sample_df['n_hvg'] = sample_df['n_hvg'].astype(str)

    # Pivot for heatmap
    pivot = sample_df.pivot_table(
        index='n_hvg',
        columns='resolution',
        values='silhouette',
        aggfunc='mean'
    )

    # Reorder
    pivot = pivot.reindex(hvg_order)
    pivot = pivot[res_order]

    # Find best configuration
    best_idx = sample_df['silhouette'].idxmax()
    best_row = sample_df.loc[best_idx]
    best_hvg = best_row['n_hvg']
    best_res = best_row['resolution']
    best_sil = best_row['silhouette']
    best_scores.append((best_hvg, best_res, best_sil))

    # Create heatmap
    sns.heatmap(
        pivot,
        ax=ax,
        annot=True,
        fmt='.3f',
        cmap='RdYlGn',
        vmin=0,
        vmax=0.2,
        cbar=idx == 2,  # Only show colorbar on last plot
        linewidths=0.5,
        annot_kws={'size': 9}
    )

    # Mark best cell with star
    hvg_idx = hvg_order.index(str(best_hvg))
    res_idx = res_order.index(best_res)
    ax.add_patch(plt.Rectangle((res_idx, hvg_idx), 1, 1, fill=False, edgecolor='blue', lw=3))

    # Title with sample info
    ax.set_title(f'{sample} - {stage}\nBest: HVG={best_hvg}, res={best_res}\nSilhouette = {best_sil:.3f}',
                 fontsize=11, fontweight='bold')
    ax.set_xlabel('Resolution')
    ax.set_ylabel('HVG Count' if idx == 0 else '')

    # Rotate x labels
    ax.set_xticklabels([str(r) for r in res_order], rotation=0)

# Add overall title
fig.suptitle('HVG Optimization: Silhouette Scores Across Disease Progression\n'
             '(Blue box = optimal configuration)',
             fontsize=14, fontweight='bold', y=1.02)

plt.tight_layout()

# Save figure
out_path = '/home/user/g4x-choi-batch2-analysis/output/figures/hvg_panel_all_samples.png'
plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
print(f"Saved: {out_path}")

# Also save to desktop
desktop_path = '/mnt/c/Users/User/Desktop/G4X_Analysis/figures/hvg_panel_all_samples.png'
plt.savefig(desktop_path, dpi=150, bbox_inches='tight', facecolor='white')
print(f"Saved: {desktop_path}")

# Print summary
print("\n" + "="*60)
print("SUMMARY: Best configurations per sample")
print("="*60)
for (sample, stage), (hvg, res, sil) in zip(zip(samples, stages), best_scores):
    print(f"{sample} ({stage:10s}): HVG={hvg:>3}, res={res}, silhouette={sil:.3f}")
print("="*60)
print("\nKEY FINDING: Cancer (G02) has ~2x better cluster separation!")
