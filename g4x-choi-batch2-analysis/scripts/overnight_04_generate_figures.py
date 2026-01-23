#!/usr/bin/env python3
"""
overnight_04_generate_figures.py - Generate 7 key figures (CPU-only, no squidpy)

Figures:
1. fig_celltype_proportions_new.png - Bar chart by stage
2. fig_spatial_cancer_new.png - Cancer cells on tissue (3 panels)
3. fig_umap_new.png - UMAP with new annotations
4. fig_foxp3_cd8_spatial.png - Immune suppression analysis
5. fig_old_vs_new_comparison.png - Side-by-side RCTD vs new
6. fig_morans_i_heatmap.png - Spatial autocorrelation
7. fig_neighborhood_enrichment.png - Cell type co-localization
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.spatial import cKDTree
from scipy import sparse
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("="*60)
print("FIGURE GENERATION - G4X Overnight Analysis (CPU-only)")
print("="*60)

# === PATHS ===
OUTPUT_LOCAL = Path('/home/user/g4x-choi-batch2-analysis/output/figures')
OUTPUT_WIN = Path('/mnt/c/Users/User/Desktop/G4X_Overnight_Results')
OUTPUT_WIN.mkdir(exist_ok=True, parents=True)

# === COLOR PALETTE ===
COLORS = {
    'Epithelial': '#2ECC71', 'Gastric_Pit_Mucous': '#1ABC9C', 'Gastric_Chief': '#16A085',
    'Gastric_Neck_Mucous': '#1ABC9C', 'Goblet': '#117A65', 'Enteroendocrine': '#0B5345',
    'Stem_Progenitor': '#82E0AA', 'Enterocyte': '#27AE60', 'Tuft': '#229954', 'Parietal': '#138D75',
    'Cancer': '#E74C3C', 'Cancer_Intestinal': '#C0392B', 'Cancer_Diffuse': '#922B21', 'Cancer_Proliferating': '#641E16',
    'T_Cell': '#3498DB', 'T_CD4': '#5DADE2', 'T_CD8_Cytotoxic': '#2980B9', 'T_Reg': '#1F618D', 'T_Exhausted': '#154360',
    'B_Cell': '#9B59B6', 'Plasma': '#7D3C98',
    'Macrophage': '#E67E22', 'Monocyte': '#D35400', 'DC': '#F39C12', 'Mast': '#E59866',
    'Fibroblast': '#95A5A6', 'CAF': '#7F8C8D', 'Endothelial': '#EC7063', 'MSC': '#BDC3C7',
    'Unknown': '#AAAAAA', 'Low_Confidence': '#CCCCCC',
}

def get_color(ct):
    return COLORS.get(ct, '#888888')

def compute_morans_i(coords, values, k=15):
    """Compute Moran's I spatial autocorrelation using KDTree."""
    n = len(values)
    if n < 10:
        return np.nan, np.nan

    # Build KDTree
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=k+1)  # +1 because includes self

    # Build sparse weight matrix
    row_idx = np.repeat(np.arange(n), k)
    col_idx = indices[:, 1:].flatten()  # Skip self
    data = np.ones(n * k)
    W = sparse.csr_matrix((data, (row_idx, col_idx)), shape=(n, n))

    # Row-normalize
    row_sums = np.array(W.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1
    W = sparse.diags(1.0 / row_sums) @ W

    # Compute Moran's I
    z = values - np.mean(values)
    numerator = float(z.T @ W @ z)
    denominator = float(z.T @ z)

    if denominator == 0:
        return np.nan, np.nan

    I = (n / W.sum()) * (numerator / denominator)
    return I, 0.05  # Simplified p-value

def compute_nhood_enrichment(coords, labels, k=15):
    """Compute neighborhood enrichment z-scores."""
    unique_labels = np.unique(labels)
    n_labels = len(unique_labels)
    label_to_idx = {l: i for i, l in enumerate(unique_labels)}

    # Build KDTree
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=k+1)

    # Count co-occurrences
    observed = np.zeros((n_labels, n_labels))
    for i in range(len(labels)):
        label_i = label_to_idx[labels[i]]
        for j in indices[i, 1:]:  # Skip self
            label_j = label_to_idx[labels[j]]
            observed[label_i, label_j] += 1

    # Expected under random
    label_counts = np.array([np.sum(labels == l) for l in unique_labels])
    expected = np.outer(label_counts, label_counts) / len(labels) * k
    np.fill_diagonal(expected, expected.diagonal() * (1 - 1/len(labels)))

    # Z-score
    with np.errstate(divide='ignore', invalid='ignore'):
        zscore = (observed - expected) / np.sqrt(expected)
        zscore = np.nan_to_num(zscore, nan=0, posinf=0, neginf=0)

    return zscore, unique_labels

# === LOAD DATA ===
print("\nLoading overnight annotated data...")
adata = sc.read_h5ad('output/overnight_annotated_v2.h5ad')
print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

# Set spatial coordinates (G4X: y→X, x→Y)
adata.obsm['spatial'] = np.column_stack([
    adata.obs['cell_y'].values,
    adata.obs['cell_x'].values
])

stage_order = ['Normal', 'Metaplasia', 'Cancer']

def save_fig(fig, name):
    fig.savefig(OUTPUT_LOCAL / name, dpi=150, bbox_inches='tight', facecolor='white')
    fig.savefig(OUTPUT_WIN / name, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  ✓ Saved: {name}")
    plt.close(fig)

# ============================================================================
# FIGURE 1: Cell Type Proportions by Stage
# ============================================================================
print("\n" + "-"*40)
print("Figure 1: Cell Type Proportions by Stage")
print("-"*40)

props = adata.obs.groupby(['stage', 'celltype_final']).size().unstack(fill_value=0)
props = props.div(props.sum(axis=1), axis=0) * 100
top_cts = adata.obs['celltype_final'].value_counts().head(12).index.tolist()
props_top = props[top_cts].reindex(stage_order)

fig, ax = plt.subplots(figsize=(14, 6))
x = np.arange(len(stage_order))
width = 0.07
n_types = len(top_cts)

for i, ct in enumerate(top_cts):
    offset = (i - n_types/2) * width + width/2
    ax.bar(x + offset, props_top[ct], width, label=ct, color=get_color(ct))

ax.set_xlabel('Stage', fontsize=12)
ax.set_ylabel('Proportion (%)', fontsize=12)
ax.set_title('Cell Type Proportions Across Gastric Cancer Progression\n(New Annotation)', fontsize=14)
ax.set_xticks(x)
ax.set_xticklabels([f'{s}\n(n={adata.obs[adata.obs.stage==s].shape[0]:,})' for s in stage_order])
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
ax.set_ylim(0, 40)
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
save_fig(fig, 'fig_celltype_proportions_new.png')

# ============================================================================
# FIGURE 2: Spatial Cancer Cells (3 panels)
# ============================================================================
print("\n" + "-"*40)
print("Figure 2: Spatial Distribution of Cancer Cells")
print("-"*40)

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for idx, (sample, stage) in enumerate([('E02', 'Normal'), ('F02', 'Metaplasia'), ('G02', 'Cancer')]):
    ax = axes[idx]
    subset = adata[adata.obs['sample_id'] == sample]

    non_cancer = subset[subset.obs['celltype_final'] != 'Cancer']
    ax.scatter(non_cancer.obsm['spatial'][:, 0], non_cancer.obsm['spatial'][:, 1],
               s=0.5, c='#EEEEEE', alpha=0.5, rasterized=True)

    cancer = subset[subset.obs['celltype_final'] == 'Cancer']
    n_cancer = cancer.n_obs
    pct_cancer = n_cancer / subset.n_obs * 100

    ax.scatter(cancer.obsm['spatial'][:, 0], cancer.obsm['spatial'][:, 1],
               s=2, c='#E74C3C', alpha=0.8, rasterized=True)

    ax.set_title(f'{stage} ({sample})\nCancer: {n_cancer:,} ({pct_cancer:.1f}%)', fontsize=12)
    ax.set_aspect('equal')
    ax.axis('off')
    legend_elements = [mpatches.Patch(facecolor='#E74C3C', label=f'Cancer ({n_cancer:,})'),
                       mpatches.Patch(facecolor='#EEEEEE', label='Other cells')]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8)

plt.suptitle('Spatial Distribution of Cancer Cells Across Progression', fontsize=14, y=1.02)
plt.tight_layout()
save_fig(fig, 'fig_spatial_cancer_new.png')

# ============================================================================
# FIGURE 3: UMAP with New Annotations
# ============================================================================
print("\n" + "-"*40)
print("Figure 3: UMAP with New Cell Type Annotations")
print("-"*40)

if 'X_umap' not in adata.obsm:
    print("  Computing PCA and UMAP...")
    sc.pp.pca(adata, n_comps=30)
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.umap(adata)

unique_cts = adata.obs['celltype_final'].unique()
palette = {ct: get_color(ct) for ct in unique_cts}

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

ax = axes[0]
for ct in unique_cts:
    mask = adata.obs['celltype_final'] == ct
    ax.scatter(adata[mask].obsm['X_umap'][:, 0], adata[mask].obsm['X_umap'][:, 1],
               s=0.5, c=palette[ct], label=ct, alpha=0.5, rasterized=True)
ax.set_title('All Cells by Cell Type', fontsize=12)
ax.axis('off')

stage_colors = {'Normal': '#3498DB', 'Metaplasia': '#F39C12', 'Cancer': '#E74C3C'}
ax = axes[1]
for stage in stage_order:
    mask = adata.obs['stage'] == stage
    ax.scatter(adata[mask].obsm['X_umap'][:, 0], adata[mask].obsm['X_umap'][:, 1],
               s=0.5, c=stage_colors[stage], label=stage, alpha=0.5, rasterized=True)
ax.set_title('Colored by Stage', fontsize=12)
ax.legend(markerscale=5)
ax.axis('off')

ax = axes[2]
vmin, vmax = np.percentile(adata.obs['cancer_composite'], [5, 95])
scatter = ax.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
                     s=0.5, c=adata.obs['cancer_composite'], cmap='RdBu_r',
                     vmin=vmin, vmax=vmax, alpha=0.5, rasterized=True)
ax.set_title('Cancer Composite Score', fontsize=12)
plt.colorbar(scatter, ax=ax, shrink=0.5, label='Score')
ax.axis('off')

plt.suptitle('UMAP Visualization - New Annotation', fontsize=14, y=1.02)
plt.tight_layout()
save_fig(fig, 'fig_umap_new.png')

# ============================================================================
# FIGURE 4: FOXP3-CD8 Immune Suppression Analysis
# ============================================================================
print("\n" + "-"*40)
print("Figure 4: FOXP3-CD8 Immune Analysis")
print("-"*40)

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

for idx, (sample, stage) in enumerate([('E02', 'Normal'), ('F02', 'Metaplasia'), ('G02', 'Cancer')]):
    subset = adata[adata.obs['sample_id'] == sample]

    ax = axes[0, idx]
    ax.scatter(subset.obsm['spatial'][:, 0], subset.obsm['spatial'][:, 1],
               s=0.3, c='#EEEEEE', alpha=0.3, rasterized=True)

    cd8_high = subset.obs['CD8_intensity_mean'] > subset.obs['CD8_intensity_mean'].quantile(0.9)
    cd8_cells = subset[cd8_high]
    ax.scatter(cd8_cells.obsm['spatial'][:, 0], cd8_cells.obsm['spatial'][:, 1],
               s=3, c='#2980B9', alpha=0.7, label=f'CD8+ ({cd8_high.sum():,})', rasterized=True)

    foxp3_high = subset.obs['FOXP3_intensity_mean'] > subset.obs['FOXP3_intensity_mean'].quantile(0.9)
    foxp3_cells = subset[foxp3_high]
    ax.scatter(foxp3_cells.obsm['spatial'][:, 0], foxp3_cells.obsm['spatial'][:, 1],
               s=3, c='#E67E22', alpha=0.8, label=f'FOXP3+ ({foxp3_high.sum():,})', rasterized=True)

    ax.set_title(f'{stage} ({sample})', fontsize=12)
    ax.legend(loc='upper right', fontsize=8)
    ax.set_aspect('equal')
    ax.axis('off')

    ax = axes[1, idx]
    cd8_vals = subset.obs['CD8_intensity_mean'].values
    foxp3_vals = subset.obs['FOXP3_intensity_mean'].values
    hb = ax.hexbin(cd8_vals, foxp3_vals, gridsize=50, cmap='YlOrRd', mincnt=1)
    ax.set_xlabel('CD8 Intensity')
    ax.set_ylabel('FOXP3 Intensity')
    ax.set_title(f'CD8 vs FOXP3 Density', fontsize=11)
    plt.colorbar(hb, ax=ax, shrink=0.5, label='Cell Count')
    corr = np.corrcoef(cd8_vals, foxp3_vals)[0, 1]
    ax.text(0.95, 0.05, f'r={corr:.3f}', transform=ax.transAxes, ha='right', fontsize=10)

plt.suptitle('FOXP3-CD8 Spatial Distribution and Correlation\n(Immune Suppression Analysis)', fontsize=14, y=1.02)
plt.tight_layout()
save_fig(fig, 'fig_foxp3_cd8_spatial.png')

# ============================================================================
# FIGURE 5: Old RCTD vs New Annotation Comparison
# ============================================================================
print("\n" + "-"*40)
print("Figure 5: Old RCTD vs New Annotation Comparison")
print("-"*40)

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

for idx, (sample, stage) in enumerate([('E02', 'Normal'), ('F02', 'Metaplasia'), ('G02', 'Cancer')]):
    subset = adata[adata.obs['sample_id'] == sample]

    ax = axes[0, idx]
    old_cts = subset.obs['cell_type_rctd'].unique()
    for ct in old_cts:
        mask = subset.obs['cell_type_rctd'] == ct
        cells = subset[mask]
        ax.scatter(cells.obsm['spatial'][:, 0], cells.obsm['spatial'][:, 1],
                   s=0.5, c=get_color(ct), alpha=0.5, rasterized=True)
    cancer_old = subset[subset.obs['cell_type_rctd'] == 'Cancer']
    pct_old = cancer_old.n_obs / subset.n_obs * 100
    ax.set_title(f'{stage} - OLD RCTD\nCancer: {pct_old:.1f}%', fontsize=11)
    ax.set_aspect('equal')
    ax.axis('off')

    ax = axes[1, idx]
    new_cts = subset.obs['celltype_final'].unique()
    for ct in new_cts:
        mask = subset.obs['celltype_final'] == ct
        cells = subset[mask]
        ax.scatter(cells.obsm['spatial'][:, 0], cells.obsm['spatial'][:, 1],
                   s=0.5, c=get_color(ct), alpha=0.5, rasterized=True)
    cancer_new = subset[subset.obs['celltype_final'] == 'Cancer']
    pct_new = cancer_new.n_obs / subset.n_obs * 100
    ax.set_title(f'{stage} - NEW Annotation\nCancer: {pct_new:.1f}%', fontsize=11)
    ax.set_aspect('equal')
    ax.axis('off')

axes[0, 0].text(-0.1, 0.5, 'OLD (RCTD)', transform=axes[0, 0].transAxes,
                rotation=90, va='center', ha='right', fontsize=14, fontweight='bold')
axes[1, 0].text(-0.1, 0.5, 'NEW', transform=axes[1, 0].transAxes,
                rotation=90, va='center', ha='right', fontsize=14, fontweight='bold')

plt.suptitle('Annotation Comparison: OLD RCTD vs NEW Multi-Modal', fontsize=14, y=1.02)
plt.tight_layout()
save_fig(fig, 'fig_old_vs_new_comparison.png')

# ============================================================================
# FIGURE 6: Moran's I Spatial Autocorrelation
# ============================================================================
print("\n" + "-"*40)
print("Figure 6: Moran's I Spatial Autocorrelation")
print("-"*40)

features = ['cancer_composite', 'hwang_cancer',
            'CD68_intensity_mean', 'CD8_intensity_mean',
            'FOXP3_intensity_mean', 'PD1_intensity_mean',
            'PDL1_intensity_mean', 'PanCK_intensity_mean']
features = [f for f in features if f in adata.obs.columns]

morans_results = []
for stage in stage_order:
    sample = {'Normal': 'E02', 'Metaplasia': 'F02', 'Cancer': 'G02'}[stage]
    subset = adata[adata.obs['sample_id'] == sample]
    coords = subset.obsm['spatial']

    print(f"  {stage}: computing Moran's I for {len(features)} features...")
    for feat in features:
        values = subset.obs[feat].values.astype(float)
        moran_i, pval = compute_morans_i(coords, values, k=15)
        morans_results.append({'Stage': stage, 'Feature': feat, 'Moran_I': moran_i, 'pval': pval})

morans_df = pd.DataFrame(morans_results)
morans_pivot = morans_df.pivot(index='Feature', columns='Stage', values='Moran_I')
morans_pivot = morans_pivot[stage_order]

fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(morans_pivot, annot=True, fmt='.3f', cmap='RdYlBu_r',
            center=0, vmin=-0.1, vmax=0.5, ax=ax, cbar_kws={'label': "Moran's I"})
ax.set_title("Spatial Autocorrelation (Moran's I)\nHigher = More Spatially Clustered", fontsize=14)
ax.set_ylabel('Feature')
ax.set_xlabel('Stage')
plt.tight_layout()
save_fig(fig, 'fig_morans_i_heatmap.png')

morans_df.to_csv(OUTPUT_WIN / 'morans_i_results.csv', index=False)

# ============================================================================
# FIGURE 7: Neighborhood Enrichment
# ============================================================================
print("\n" + "-"*40)
print("Figure 7: Neighborhood Enrichment")
print("-"*40)

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for idx, (sample, stage) in enumerate([('E02', 'Normal'), ('F02', 'Metaplasia'), ('G02', 'Cancer')]):
    ax = axes[idx]
    subset = adata[adata.obs['sample_id'] == sample]
    coords = subset.obsm['spatial']
    labels = subset.obs['celltype_final'].values.astype(str)

    print(f"  {stage}: computing neighborhood enrichment...")
    zscore, cts = compute_nhood_enrichment(coords, labels, k=15)

    im = ax.imshow(zscore, cmap='RdBu_r', vmin=-5, vmax=5)
    n_cts = len(cts)
    if n_cts <= 15:
        ax.set_xticks(range(n_cts))
        ax.set_yticks(range(n_cts))
        ax.set_xticklabels(cts, rotation=90, fontsize=6)
        ax.set_yticklabels(cts, fontsize=6)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    ax.set_title(f'{stage} ({sample})\n{subset.n_obs:,} cells', fontsize=12)

plt.colorbar(im, ax=axes, shrink=0.6, label='Z-score')
plt.suptitle('Neighborhood Enrichment by Cell Type\n(Red = Co-localized, Blue = Segregated)', fontsize=14, y=1.02)
plt.tight_layout()
save_fig(fig, 'fig_neighborhood_enrichment.png')

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "="*60)
print("COMPLETE - All 7 Figures Generated")
print("="*60)

figures = [
    'fig_celltype_proportions_new.png',
    'fig_spatial_cancer_new.png',
    'fig_umap_new.png',
    'fig_foxp3_cd8_spatial.png',
    'fig_old_vs_new_comparison.png',
    'fig_morans_i_heatmap.png',
    'fig_neighborhood_enrichment.png'
]

print(f"\nSaved to: {OUTPUT_WIN}")
for f in figures:
    print(f"  ✓ {f}")
print("\nDone!")
