#!/usr/bin/env python
"""
Comprehensive QC with per-sample figure generation
G4X Gastric Cancer Hackathon 2026
"""

import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# === CONFIG ===
DATA_DIR = Path('/home/user/spatial-hackathon-2026/results/g4x_choi_batch2/annotated_v2')
OUT_DIR = Path('/home/user/g4x-choi-batch2-analysis/results/qc_figures')
OUT_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES = ['A01', 'B01', 'C01', 'D01', 'E01', 'F01', 'G01', 'H01']

# RNA-Protein pairs for concordance
RNA_PROTEIN_PAIRS = [
    ('CD8A', 'CD8'), ('PDCD1', 'PD1'), ('MS4A1', 'CD20'),
    ('MKI67', 'KI67'), ('CD4', 'CD4'), ('CD3E', 'CD3'),
    ('ACTA2', 'aSMA'), ('PTPRC', 'CD45'), ('PECAM1', 'CD31')
]

print("="*60)
print("G4X COMPREHENSIVE QC - FIGURE GENERATION")
print("="*60)

# === PHASE 1: PER-SAMPLE QC ===
print("\n[PHASE 1] Per-Sample QC Figures...")

summary_stats = []

for sample in SAMPLES:
    print(f"\n  Processing {sample}...")
    fpath = DATA_DIR / f'{sample}_annotated.h5ad'
    if not fpath.exists():
        print(f"    WARNING: {fpath} not found, skipping")
        continue

    adata = sc.read(fpath)

    # Calculate QC metrics if not present
    if 'total_counts' not in adata.obs:
        sc.pp.calculate_qc_metrics(adata, inplace=True)

    # Store summary stats
    summary_stats.append({
        'sample': sample,
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'median_counts': adata.obs['total_counts'].median() if 'total_counts' in adata.obs else np.nan,
        'median_genes': adata.obs['n_genes_by_counts'].median() if 'n_genes_by_counts' in adata.obs else np.nan,
    })

    # --- Figure 1: Violin QC Metrics ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle(f'{sample} - Cell-Level QC Metrics', fontsize=14, fontweight='bold')

    metrics = ['total_counts', 'n_genes_by_counts']
    for i, metric in enumerate(metrics):
        if metric in adata.obs:
            axes[i].violinplot(adata.obs[metric].values, positions=[0], showmedians=True)
            axes[i].set_title(metric)
            axes[i].set_ylabel(metric)
            median_val = adata.obs[metric].median()
            axes[i].axhline(median_val, color='red', linestyle='--', alpha=0.5)
            axes[i].text(0.02, 0.98, f'Median: {median_val:.0f}', transform=axes[i].transAxes, va='top')

    # Histogram: genes per cell
    if 'n_genes_by_counts' in adata.obs:
        axes[2].hist(adata.obs['n_genes_by_counts'], bins=50, edgecolor='black', alpha=0.7)
        axes[2].set_xlabel('Genes per cell')
        axes[2].set_ylabel('Count')
        axes[2].set_title('Gene Detection Distribution')
        axes[2].axvline(adata.obs['n_genes_by_counts'].median(), color='red', linestyle='--')

    plt.tight_layout()
    plt.savefig(OUT_DIR / f'{sample}_violin_qc.png', dpi=150, bbox_inches='tight')
    plt.close()

    # --- Figure 2: Spatial QC ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'{sample} - Spatial QC', fontsize=14, fontweight='bold')

    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']

        # Total counts spatial
        if 'total_counts' in adata.obs:
            scatter = axes[0].scatter(coords[:, 0], coords[:, 1],
                                      c=adata.obs['total_counts'], cmap='viridis', s=1, alpha=0.5)
            plt.colorbar(scatter, ax=axes[0], label='Total Counts')
            axes[0].set_title('Spatial: Total Counts')
            axes[0].set_aspect('equal')

        # Genes detected spatial
        if 'n_genes_by_counts' in adata.obs:
            scatter = axes[1].scatter(coords[:, 0], coords[:, 1],
                                      c=adata.obs['n_genes_by_counts'], cmap='plasma', s=1, alpha=0.5)
            plt.colorbar(scatter, ax=axes[1], label='Genes Detected')
            axes[1].set_title('Spatial: Genes per Cell')
            axes[1].set_aspect('equal')

        # Cell type spatial
        if 'cell_type' in adata.obs:
            cell_types = adata.obs['cell_type'].astype('category')
            colors = plt.cm.tab20(np.linspace(0, 1, len(cell_types.cat.categories)))
            color_map = dict(zip(cell_types.cat.categories, colors))
            cell_colors = [color_map[ct] for ct in adata.obs['cell_type']]
            axes[2].scatter(coords[:, 0], coords[:, 1], c=cell_colors, s=1, alpha=0.5)
            axes[2].set_title('Spatial: Cell Types')
            axes[2].set_aspect('equal')
            # Legend
            handles = [plt.scatter([], [], c=[color_map[ct]], label=ct, s=20)
                       for ct in list(cell_types.cat.categories)[:10]]
            axes[2].legend(handles=handles, loc='upper right', fontsize=6, ncol=2)

    for ax in axes:
        ax.set_xlabel('X')
        ax.set_ylabel('Y')

    plt.tight_layout()
    plt.savefig(OUT_DIR / f'{sample}_spatial_qc.png', dpi=150, bbox_inches='tight')
    plt.close()

    # --- Figure 3: Protein QC ---
    if 'protein' in adata.obsm:
        protein_data = pd.DataFrame(adata.obsm['protein'])
        if 'protein_names' in adata.uns:
            protein_data.columns = adata.uns['protein_names']

        fig, axes = plt.subplots(2, 1, figsize=(14, 8))
        fig.suptitle(f'{sample} - Protein QC', fontsize=14, fontweight='bold')

        # Violin plot of all proteins
        protein_melted = protein_data.melt(var_name='Protein', value_name='Expression')
        sns.violinplot(data=protein_melted, x='Protein', y='Expression', ax=axes[0], scale='width')
        axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=45, ha='right')
        axes[0].set_title('Protein Expression Distribution')

        # Protein correlation heatmap
        corr = protein_data.corr()
        sns.heatmap(corr, annot=False, cmap='coolwarm', center=0, ax=axes[1],
                    xticklabels=True, yticklabels=True)
        axes[1].set_title('Protein-Protein Correlations')

        plt.tight_layout()
        plt.savefig(OUT_DIR / f'{sample}_protein_qc.png', dpi=150, bbox_inches='tight')
        plt.close()

    print(f"    {sample}: {adata.n_obs:,} cells, {adata.n_vars} genes")

# === PHASE 2: CROSS-SAMPLE COMPARISON ===
print("\n[PHASE 2] Cross-Sample Comparison...")

# Load all samples for comparison
all_obs = []
for sample in SAMPLES:
    fpath = DATA_DIR / f'{sample}_annotated.h5ad'
    if fpath.exists():
        adata = sc.read(fpath)
        if 'total_counts' not in adata.obs:
            sc.pp.calculate_qc_metrics(adata, inplace=True)
        obs = adata.obs.copy()
        obs['sample'] = sample
        all_obs.append(obs)

combined_obs = pd.concat(all_obs, ignore_index=True)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Cross-Sample QC Comparison', fontsize=14, fontweight='bold')

# Ridge plot: total counts
if 'total_counts' in combined_obs:
    for i, sample in enumerate(SAMPLES):
        subset = combined_obs[combined_obs['sample'] == sample]['total_counts']
        if len(subset) > 0:
            axes[0, 0].hist(subset, bins=50, alpha=0.5, label=sample, density=True)
    axes[0, 0].set_xlabel('Total Counts')
    axes[0, 0].set_ylabel('Density')
    axes[0, 0].set_title('Total Counts by Sample')
    axes[0, 0].legend(fontsize=8)

# Box plot: cells per sample
cells_per_sample = combined_obs.groupby('sample').size()
axes[0, 1].bar(cells_per_sample.index, cells_per_sample.values, color='steelblue', edgecolor='black')
axes[0, 1].set_xlabel('Sample')
axes[0, 1].set_ylabel('Number of Cells')
axes[0, 1].set_title('Cells per Sample')
for i, v in enumerate(cells_per_sample.values):
    axes[0, 1].text(i, v + 1000, f'{v:,}', ha='center', fontsize=8)

# Stacked bar: cell type proportions
if 'cell_type' in combined_obs:
    ct_props = pd.crosstab(combined_obs['sample'], combined_obs['cell_type'], normalize='index')
    ct_props.plot(kind='bar', stacked=True, ax=axes[1, 0], colormap='tab20')
    axes[1, 0].set_xlabel('Sample')
    axes[1, 0].set_ylabel('Proportion')
    axes[1, 0].set_title('Cell Type Proportions by Sample')
    axes[1, 0].legend(loc='upper right', fontsize=6, ncol=2)
    axes[1, 0].tick_params(axis='x', rotation=45)

# Summary stats table
summary_df = pd.DataFrame(summary_stats)
axes[1, 1].axis('off')
table = axes[1, 1].table(cellText=summary_df.values,
                          colLabels=summary_df.columns,
                          loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1.2, 1.5)
axes[1, 1].set_title('Summary Statistics', pad=20)

plt.tight_layout()
plt.savefig(OUT_DIR / 'cross_sample_comparison.png', dpi=150, bbox_inches='tight')
plt.close()

# === PHASE 3: RNA-PROTEIN CONCORDANCE ===
print("\n[PHASE 3] RNA-Protein Concordance...")

# Use first sample with protein data
sample = 'A01'
adata = sc.read(DATA_DIR / f'{sample}_annotated.h5ad')

if 'protein' in adata.obsm and 'protein_names' in adata.uns:
    protein_names = list(adata.uns['protein_names'])
    protein_df = pd.DataFrame(adata.obsm['protein'], columns=protein_names)

    fig, axes = plt.subplots(3, 3, figsize=(12, 12))
    fig.suptitle('RNA-Protein Concordance', fontsize=14, fontweight='bold')

    concordance_results = []
    for idx, (rna, prot) in enumerate(RNA_PROTEIN_PAIRS[:9]):
        ax = axes[idx // 3, idx % 3]

        if rna in adata.var_names and prot in protein_names:
            rna_expr = adata[:, rna].X.toarray().flatten() if hasattr(adata[:, rna].X, 'toarray') else adata[:, rna].X.flatten()
            prot_expr = protein_df[prot].values

            # Subsample for plotting
            n_plot = min(5000, len(rna_expr))
            idx_plot = np.random.choice(len(rna_expr), n_plot, replace=False)

            ax.scatter(rna_expr[idx_plot], prot_expr[idx_plot], s=1, alpha=0.3)

            # Correlation
            r, p = stats.spearmanr(rna_expr, prot_expr)
            concordance_results.append({'RNA': rna, 'Protein': prot, 'Spearman_r': r, 'p_value': p})

            ax.set_xlabel(f'{rna} (RNA)')
            ax.set_ylabel(f'{prot} (Protein)')
            ax.set_title(f'r={r:.3f}')
        else:
            ax.text(0.5, 0.5, f'{rna}/{prot}\nNot found', ha='center', va='center')
            ax.set_title('Missing')

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'rna_protein_concordance.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Save concordance table
    conc_df = pd.DataFrame(concordance_results)
    conc_df.to_csv(OUT_DIR / 'rna_protein_concordance.csv', index=False)
    print(f"  RNA-Protein correlations saved")

# === PHASE 4: SAVE SUMMARY ===
print("\n[PHASE 4] Saving Summary...")

summary_df = pd.DataFrame(summary_stats)
summary_df.to_csv(OUT_DIR / 'qc_summary_stats.csv', index=False)

print("\n" + "="*60)
print("QC COMPLETE")
print("="*60)
print(f"Figures saved to: {OUT_DIR}")
print(f"Files generated:")
for f in sorted(OUT_DIR.glob('*.png')):
    print(f"  - {f.name}")
print(f"\nSummary: {OUT_DIR / 'qc_summary_stats.csv'}")
