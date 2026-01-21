#!/usr/bin/env python3
"""
PCA Analyses for Spatial Hackathon 2026
========================================

Generates multiple PCA-based figures:
- Fig 16: Bulk RNA PCA (pseudobulk)
- Fig 17: Spatial PCA (kernel-weighted)
- Fig 18: ncRNA PCA (types only)
- Fig 19: ncRNA PCA (integrated)
- Fig 20: Multi-modal PCA (all features)
- Fig 21: Separation methods comparison

Author: Max Van Belkum
Date: 2026-01-20
"""

import matplotlib
matplotlib.use('Agg')

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import silhouette_score, accuracy_score
from sklearn.model_selection import LeaveOneOut
import warnings

warnings.filterwarnings('ignore')


def safe_standardize(X):
    """Standardize data, handling NaN values safely."""
    scaler = StandardScaler()
    scaled = scaler.fit_transform(X)
    # Handle NaN from constant columns (zero variance)
    scaled = np.nan_to_num(scaled, nan=0.0, posinf=0.0, neginf=0.0)
    return scaled, scaler

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
ADATA_DIR = OUTPUT_DIR / "adata"
TABLE_DIR = OUTPUT_DIR / "tables"
FIG_DIR = OUTPUT_DIR / "figures" / "showcase_v2"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Visual settings
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 10

# Sample metadata
PDAC_METADATA = {
    "YP03A": {"patient": "YP03", "timepoint": "Pre", "response": "NR"},
    "YP03C": {"patient": "YP03", "timepoint": "Post", "response": "NR"},
    "YP04C": {"patient": "YP04", "timepoint": "Post", "response": "NR"},
    "YP12A": {"patient": "YP12", "timepoint": "Pre", "response": "R"},
    "YP12C": {"patient": "YP12", "timepoint": "Post", "response": "R"},
    "YP15A": {"patient": "YP15", "timepoint": "Pre", "response": "R"},
    "YP15C": {"patient": "YP15", "timepoint": "Post", "response": "R"},
}

RESPONSE_COLORS = {'R': '#2ecc71', 'NR': '#e74c3c'}
TIMEPOINT_MARKERS = {'Pre': 'o', 'Post': 's'}


def load_all_samples() -> Dict[str, sc.AnnData]:
    """Load all PDAC samples."""
    samples = {}
    for sample in PDAC_METADATA.keys():
        adata_path = ADATA_DIR / "polymathic" / f"{sample}_polymathic.h5ad"
        if not adata_path.exists():
            adata_path = ADATA_DIR / "annotated" / f"{sample}_annotated.h5ad"
        if adata_path.exists():
            samples[sample] = sc.read_h5ad(adata_path)
    return samples


def create_pseudobulk(samples: Dict[str, sc.AnnData], n_hvg: int = 500) -> Tuple[pd.DataFrame, List[str]]:
    """Create pseudobulk expression matrix by summing counts per sample."""
    # Find common genes
    common_genes = None
    for adata in samples.values():
        genes = set(adata.var_names)
        if common_genes is None:
            common_genes = genes
        else:
            common_genes = common_genes & genes

    common_genes = sorted(common_genes)
    print(f"  Common genes across samples: {len(common_genes)}")

    # Create pseudobulk matrix
    pseudobulk = {}
    for sample, adata in samples.items():
        # Sum expression across cells
        if hasattr(adata.X, 'toarray'):
            expr = np.array(adata[:, common_genes].X.toarray().sum(axis=0)).flatten()
        else:
            expr = np.array(adata[:, common_genes].X.sum(axis=0)).flatten()
        pseudobulk[sample] = expr

    df = pd.DataFrame(pseudobulk, index=common_genes).T

    # Normalize (CPM-like)
    df = df.div(df.sum(axis=1), axis=0) * 1e6

    # Log transform
    df = np.log1p(df)

    # Select highly variable genes
    variances = df.var()
    hvg = variances.nlargest(n_hvg).index.tolist()

    return df[hvg], hvg


def get_ncrna_genes(gene_names: List[str]) -> List[str]:
    """Identify ncRNA genes by naming patterns."""
    ncrna_patterns = ['LINC', 'MIR', 'SNORD', 'SNORA', '-AS', 'AC0', 'AL0', 'lncRNA']
    ncrna = []
    for gene in gene_names:
        for pattern in ncrna_patterns:
            if pattern in gene.upper() or gene.upper().startswith(pattern):
                ncrna.append(gene)
                break
    return list(set(ncrna))


def load_polymathic_features() -> pd.DataFrame:
    """Load polymathic analysis features for each sample."""
    df = pd.read_csv(TABLE_DIR / "polymathic_analysis_results.csv")
    df = df.set_index('sample')
    # Keep only numeric columns
    df = df.select_dtypes(include=[np.number])
    return df


def load_cell_type_proportions() -> pd.DataFrame:
    """Load cell type proportions for each sample."""
    df = pd.read_csv(TABLE_DIR / "cell_type_proportions.csv")
    df = df.set_index('sample')
    # Keep only numeric columns
    df = df.select_dtypes(include=[np.number])
    return df


def fig16_bulk_rna_pca():
    """Figure 16: Bulk RNA PCA (pseudobulk from single-cell data)."""
    print("Generating Figure 16: Bulk RNA PCA")

    samples = load_all_samples()
    if len(samples) < 2:
        print("  ERROR: Not enough samples for PCA")
        return

    # Create pseudobulk
    pseudobulk, hvg = create_pseudobulk(samples, n_hvg=500)
    print(f"  Pseudobulk matrix: {pseudobulk.shape}")

    # Standardize and PCA
    scaler = StandardScaler()
    scaled = scaler.fit_transform(pseudobulk)
    pca = PCA(n_components=min(len(samples) - 1, 5))
    pcs = pca.fit_transform(scaled)

    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel A: PC1 vs PC2
    ax = axes[0]
    for i, sample in enumerate(pseudobulk.index):
        meta = PDAC_METADATA[sample]
        ax.scatter(pcs[i, 0], pcs[i, 1],
                   c=RESPONSE_COLORS[meta['response']],
                   marker=TIMEPOINT_MARKERS[meta['timepoint']],
                   s=150, edgecolors='black', linewidth=1.5, zorder=3)
        ax.annotate(sample, (pcs[i, 0], pcs[i, 1]),
                    textcoords='offset points', xytext=(5, 5), fontsize=9)

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title("A) PC1 vs PC2")
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

    # Panel B: Variance explained
    ax = axes[1]
    n_pcs = len(pca.explained_variance_ratio_)
    ax.bar(range(1, n_pcs + 1), pca.explained_variance_ratio_ * 100, color='steelblue')
    ax.set_xlabel("Principal Component")
    ax.set_ylabel("Variance Explained (%)")
    ax.set_title("B) Variance Explained")
    ax.set_xticks(range(1, n_pcs + 1))

    cumsum = np.cumsum(pca.explained_variance_ratio_) * 100
    ax2 = ax.twinx()
    ax2.plot(range(1, n_pcs + 1), cumsum, 'r-o', label='Cumulative')
    ax2.set_ylabel("Cumulative (%)", color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.set_ylim(0, 105)

    # Panel C: Top gene loadings for PC1
    ax = axes[2]
    loadings = pd.Series(pca.components_[0], index=hvg)
    top_pos = loadings.nlargest(10)
    top_neg = loadings.nsmallest(10)
    top_genes = pd.concat([top_pos, top_neg]).sort_values()

    colors = ['#e74c3c' if v < 0 else '#2ecc71' for v in top_genes.values]
    ax.barh(range(len(top_genes)), top_genes.values, color=colors)
    ax.set_yticks(range(len(top_genes)))
    ax.set_yticklabels(top_genes.index, fontsize=8)
    ax.set_xlabel("PC1 Loading")
    ax.set_title("C) Top PC1 Gene Loadings")
    ax.axvline(x=0, color='black', linewidth=0.5)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=RESPONSE_COLORS['R'],
               markersize=10, label='Responder'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=RESPONSE_COLORS['NR'],
               markersize=10, label='Non-Responder'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
               markersize=10, label='Pre-treatment'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='gray',
               markersize=10, label='Post-treatment'),
    ]
    fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.99, 0.98))

    plt.suptitle("Figure 16: Bulk RNA PCA (Pseudobulk from Single-Cell Data, Top 500 HVGs)",
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(FIG_DIR / "fig16_bulk_rna_pca.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig16_bulk_rna_pca.pdf", bbox_inches='tight')
    plt.close()
    print(f"  Saved fig16_bulk_rna_pca.png/pdf")


def fig17_spatial_pca():
    """Figure 17: Spatial PCA incorporating spatial coordinates."""
    print("Generating Figure 17: Spatial PCA")

    samples = load_all_samples()
    if len(samples) < 2:
        print("  ERROR: Not enough samples")
        return

    # Compute spatial autocorrelation features per sample
    spatial_features = []
    sample_names = []

    for sample, adata in samples.items():
        if 'spatial' not in adata.obsm:
            print(f"  WARNING: {sample} has no spatial coordinates")
            continue

        coords = adata.obsm['spatial']

        # Compute spatial statistics
        features = {
            'sample': sample,
            'n_cells': len(coords),
            'spatial_extent_x': coords[:, 0].max() - coords[:, 0].min(),
            'spatial_extent_y': coords[:, 1].max() - coords[:, 1].min(),
            'cell_density': len(coords) / ((coords[:, 0].max() - coords[:, 0].min()) *
                                           (coords[:, 1].max() - coords[:, 1].min()) + 1),
            'coord_var_x': np.var(coords[:, 0]),
            'coord_var_y': np.var(coords[:, 1]),
            'coord_cov': np.cov(coords.T)[0, 1] if coords.shape[1] >= 2 else 0,
        }

        # Add mean expression of top genes if available
        if hasattr(adata.X, 'toarray'):
            mean_expr = np.array(adata.X.toarray().mean(axis=0)).flatten()
        else:
            mean_expr = np.array(adata.X.mean(axis=0)).flatten()
        features['mean_total_expr'] = mean_expr.mean()
        features['var_total_expr'] = mean_expr.var()

        spatial_features.append(features)
        sample_names.append(sample)

    df = pd.DataFrame(spatial_features)
    df = df.set_index('sample')

    # Add polymathic features
    poly_df = load_polymathic_features()
    spatial_cols = [c for c in poly_df.columns if any(x in c for x in ['centrality', 'topology', 'betti'])]
    # Exclude non-numeric columns
    if spatial_cols:
        poly_numeric = poly_df[spatial_cols].select_dtypes(include=[np.number])
        if not poly_numeric.empty:
            df = df.join(poly_numeric, how='left')

    # Fill NaN with column means (only for numeric columns)
    df = df.select_dtypes(include=[np.number])
    df = df.fillna(df.mean())

    # Standardize and PCA
    scaler = StandardScaler()
    scaled = scaler.fit_transform(df)
    n_components = min(len(df) - 1, 5, len(df.columns))
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(scaled)

    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel A: PC1 vs PC2
    ax = axes[0]
    for i, sample in enumerate(df.index):
        meta = PDAC_METADATA[sample]
        ax.scatter(pcs[i, 0], pcs[i, 1],
                   c=RESPONSE_COLORS[meta['response']],
                   marker=TIMEPOINT_MARKERS[meta['timepoint']],
                   s=150, edgecolors='black', linewidth=1.5, zorder=3)
        ax.annotate(sample, (pcs[i, 0], pcs[i, 1]),
                    textcoords='offset points', xytext=(5, 5), fontsize=9)

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title("A) Spatial Features PC1 vs PC2")
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

    # Panel B: Feature loadings heatmap
    ax = axes[1]
    loadings_df = pd.DataFrame(pca.components_[:3].T, index=df.columns, columns=['PC1', 'PC2', 'PC3'])
    sns.heatmap(loadings_df, cmap='RdBu_r', center=0, annot=True, fmt='.2f',
                ax=ax, cbar_kws={'shrink': 0.8})
    ax.set_title("B) Feature Loadings")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)

    # Panel C: PC1 by response
    ax = axes[2]
    pc1_r = [pcs[i, 0] for i, s in enumerate(df.index) if PDAC_METADATA[s]['response'] == 'R']
    pc1_nr = [pcs[i, 0] for i, s in enumerate(df.index) if PDAC_METADATA[s]['response'] == 'NR']

    positions = [0.8, 1.2]
    bp = ax.boxplot([pc1_r, pc1_nr], positions=positions, widths=0.3, patch_artist=True)
    bp['boxes'][0].set_facecolor(RESPONSE_COLORS['R'])
    bp['boxes'][1].set_facecolor(RESPONSE_COLORS['NR'])

    # Add individual points
    for i, (pc1, label) in enumerate([(pc1_r, 'R'), (pc1_nr, 'NR')]):
        x = np.random.normal(positions[i], 0.04, len(pc1))
        ax.scatter(x, pc1, c=RESPONSE_COLORS[label], s=80, edgecolors='black', zorder=3)

    ax.set_xticks(positions)
    ax.set_xticklabels(['Responder', 'Non-Responder'])
    ax.set_ylabel("PC1 Score")
    ax.set_title("C) PC1 by Response")

    plt.suptitle("Figure 17: Spatial PCA (Spatial + Topology Features)",
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(FIG_DIR / "fig17_spatial_pca.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig17_spatial_pca.pdf", bbox_inches='tight')
    plt.close()
    print(f"  Saved fig17_spatial_pca.png/pdf")


def fig18_ncrna_pca():
    """Figure 18: ncRNA PCA (non-coding RNA genes only)."""
    print("Generating Figure 18: ncRNA PCA")

    samples = load_all_samples()
    if len(samples) < 2:
        print("  ERROR: Not enough samples")
        return

    # Find common genes and identify ncRNAs
    common_genes = None
    for adata in samples.values():
        genes = set(adata.var_names)
        if common_genes is None:
            common_genes = genes
        else:
            common_genes = common_genes & genes

    common_genes = sorted(common_genes)
    ncrna_genes = get_ncrna_genes(common_genes)
    print(f"  Found {len(ncrna_genes)} ncRNA genes out of {len(common_genes)} common genes")

    if len(ncrna_genes) < 10:
        print("  WARNING: Very few ncRNA genes found. Figure may be limited.")

    # Create pseudobulk for ncRNAs
    pseudobulk = {}
    for sample, adata in samples.items():
        genes_to_use = [g for g in ncrna_genes if g in adata.var_names]
        if hasattr(adata.X, 'toarray'):
            expr = np.array(adata[:, genes_to_use].X.toarray().sum(axis=0)).flatten()
        else:
            expr = np.array(adata[:, genes_to_use].X.sum(axis=0)).flatten()
        pseudobulk[sample] = expr

    df = pd.DataFrame(pseudobulk, index=[g for g in ncrna_genes if g in adata.var_names]).T

    # Normalize and log
    df = df.div(df.sum(axis=1) + 1, axis=0) * 1e6
    df = np.log1p(df)

    # Remove zero-variance columns
    df = df.loc[:, df.var() > 0]
    print(f"  ncRNA matrix after filtering: {df.shape}")

    if df.shape[1] < 3:
        print("  ERROR: Not enough variable ncRNAs for PCA")
        return

    # Select top variable ncRNAs
    n_top = min(50, df.shape[1])
    top_var = df.var().nlargest(n_top).index
    df_top = df[top_var]

    # Fill any NaN values with column means
    df_top = df_top.fillna(df_top.mean())
    # If still NaN (all zeros), fill with 0
    df_top = df_top.fillna(0)

    # PCA
    scaler = StandardScaler()
    scaled = scaler.fit_transform(df_top)
    # Handle any NaN from scaling (constant columns)
    scaled = np.nan_to_num(scaled, nan=0.0)
    n_components = min(len(df_top) - 1, 3)
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(scaled)

    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel A: PC1 vs PC2
    ax = axes[0]
    for i, sample in enumerate(df.index):
        meta = PDAC_METADATA[sample]
        ax.scatter(pcs[i, 0], pcs[i, 1] if pcs.shape[1] > 1 else 0,
                   c=RESPONSE_COLORS[meta['response']],
                   marker=TIMEPOINT_MARKERS[meta['timepoint']],
                   s=150, edgecolors='black', linewidth=1.5, zorder=3)
        ax.annotate(sample, (pcs[i, 0], pcs[i, 1] if pcs.shape[1] > 1 else 0),
                    textcoords='offset points', xytext=(5, 5), fontsize=9)

    var_exp = pca.explained_variance_ratio_
    ax.set_xlabel(f"PC1 ({var_exp[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_exp[1]*100:.1f}%)" if len(var_exp) > 1 else "PC2")
    ax.set_title("A) ncRNA Expression PC1 vs PC2")

    # Panel B: Top ncRNA loadings
    ax = axes[1]
    loadings = pd.Series(pca.components_[0], index=top_var)
    top_genes = loadings.abs().nlargest(15).index
    top_loadings = loadings[top_genes].sort_values()

    colors = ['#e74c3c' if v < 0 else '#2ecc71' for v in top_loadings.values]
    ax.barh(range(len(top_loadings)), top_loadings.values, color=colors)
    ax.set_yticks(range(len(top_loadings)))
    ax.set_yticklabels(top_loadings.index, fontsize=8)
    ax.set_xlabel("PC1 Loading")
    ax.set_title("B) Top ncRNA Loadings (PC1)")
    ax.axvline(x=0, color='black', linewidth=0.5)

    # Panel C: ncRNA type breakdown
    ax = axes[2]
    type_counts = {'LINC': 0, 'MIR': 0, 'SNORD/A': 0, 'AS': 0, 'Other': 0}
    for gene in ncrna_genes:
        if 'LINC' in gene.upper():
            type_counts['LINC'] += 1
        elif 'MIR' in gene.upper():
            type_counts['MIR'] += 1
        elif 'SNOR' in gene.upper():
            type_counts['SNORD/A'] += 1
        elif '-AS' in gene.upper():
            type_counts['AS'] += 1
        else:
            type_counts['Other'] += 1

    ax.pie(type_counts.values(), labels=type_counts.keys(), autopct='%1.1f%%',
           colors=sns.color_palette('Set2', len(type_counts)))
    ax.set_title(f"C) ncRNA Types (n={len(ncrna_genes)})")

    plt.suptitle("Figure 18: ncRNA PCA (Non-coding RNA Expression)",
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(FIG_DIR / "fig18_ncrna_pca.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig18_ncrna_pca.pdf", bbox_inches='tight')
    plt.close()
    print(f"  Saved fig18_ncrna_pca.png/pdf")


def fig19_ncrna_integrated_pca():
    """Figure 19: ncRNA PCA integrated with cell type proportions."""
    print("Generating Figure 19: ncRNA Integrated PCA")

    samples = load_all_samples()
    if len(samples) < 2:
        print("  ERROR: Not enough samples")
        return

    # Get ncRNA expression (simplified)
    common_genes = None
    for adata in samples.values():
        genes = set(adata.var_names)
        if common_genes is None:
            common_genes = genes
        else:
            common_genes = common_genes & genes

    ncrna_genes = get_ncrna_genes(list(common_genes))

    # Create ncRNA summary stats per sample
    ncrna_features = {}
    for sample, adata in samples.items():
        genes_to_use = [g for g in ncrna_genes if g in adata.var_names]
        if genes_to_use:
            if hasattr(adata.X, 'toarray'):
                expr = adata[:, genes_to_use].X.toarray()
            else:
                expr = adata[:, genes_to_use].X
            ncrna_features[sample] = {
                'ncrna_total': np.sum(expr),
                'ncrna_mean': np.mean(expr),
                'ncrna_var': np.var(expr.sum(axis=1)),
                'ncrna_pct_expressing': (expr.sum(axis=1) > 0).mean(),
            }

    ncrna_df = pd.DataFrame(ncrna_features).T

    # Load cell type proportions
    ct_df = load_cell_type_proportions()

    # Combine features
    combined = ncrna_df.join(ct_df, how='inner')
    combined = combined.fillna(0)
    # Keep only numeric columns
    combined = combined.select_dtypes(include=[np.number])
    print(f"  Combined features: {combined.shape}")

    if combined.shape[0] < 2:
        print("  ERROR: Not enough samples with all features")
        return

    # PCA
    scaled, _ = safe_standardize(combined)
    n_components = min(len(combined) - 1, 3)
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(scaled)

    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel A: PC1 vs PC2
    ax = axes[0]
    for i, sample in enumerate(combined.index):
        meta = PDAC_METADATA[sample]
        ax.scatter(pcs[i, 0], pcs[i, 1] if pcs.shape[1] > 1 else 0,
                   c=RESPONSE_COLORS[meta['response']],
                   marker=TIMEPOINT_MARKERS[meta['timepoint']],
                   s=150, edgecolors='black', linewidth=1.5, zorder=3)
        ax.annotate(sample, (pcs[i, 0], pcs[i, 1] if pcs.shape[1] > 1 else 0),
                    textcoords='offset points', xytext=(5, 5), fontsize=9)

    var_exp = pca.explained_variance_ratio_
    ax.set_xlabel(f"PC1 ({var_exp[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_exp[1]*100:.1f}%)" if len(var_exp) > 1 else "PC2")
    ax.set_title("A) Integrated PC1 vs PC2")

    # Panel B: Feature contributions
    ax = axes[1]
    loadings = pd.DataFrame(pca.components_.T, index=combined.columns,
                            columns=[f'PC{i+1}' for i in range(n_components)])
    top_features = loadings['PC1'].abs().nlargest(12).index
    loadings_top = loadings.loc[top_features, ['PC1', 'PC2']] if n_components > 1 else loadings.loc[top_features, ['PC1']]

    sns.heatmap(loadings_top, cmap='RdBu_r', center=0, annot=True, fmt='.2f', ax=ax)
    ax.set_title("B) Feature Loadings (Top 12)")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=9)

    # Panel C: Correlation with response
    ax = axes[2]
    pc1_values = pd.Series(pcs[:, 0], index=combined.index)
    response_binary = pd.Series([1 if PDAC_METADATA[s]['response'] == 'R' else 0
                                 for s in combined.index], index=combined.index)

    ax.scatter(pc1_values, response_binary + np.random.normal(0, 0.05, len(response_binary)),
               c=[RESPONSE_COLORS[PDAC_METADATA[s]['response']] for s in combined.index],
               s=100, edgecolors='black')

    for s in combined.index:
        ax.annotate(s, (pc1_values[s], 1 if PDAC_METADATA[s]['response'] == 'R' else 0),
                    textcoords='offset points', xytext=(5, 5), fontsize=8)

    ax.set_xlabel("PC1 Score")
    ax.set_ylabel("Response (1=R, 0=NR)")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(['NR', 'R'])
    ax.set_title("C) PC1 vs Response")

    plt.suptitle("Figure 19: ncRNA Integrated PCA (ncRNA + Cell Types)",
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(FIG_DIR / "fig19_ncrna_integrated_pca.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig19_ncrna_integrated_pca.pdf", bbox_inches='tight')
    plt.close()
    print(f"  Saved fig19_ncrna_integrated_pca.png/pdf")


def fig20_multimodal_pca():
    """Figure 20: Integrated multi-modal PCA with all features."""
    print("Generating Figure 20: Multi-modal Integrated PCA")

    samples = load_all_samples()
    if len(samples) < 2:
        print("  ERROR: Not enough samples")
        return

    # Create pseudobulk for top HVGs
    pseudobulk, hvg = create_pseudobulk(samples, n_hvg=100)

    # Load other features
    ct_df = load_cell_type_proportions()
    poly_df = load_polymathic_features()

    # Select relevant polymathic columns
    poly_cols = [c for c in poly_df.columns if any(x in c for x in
                 ['centrality', 'topology', 'betti', 'mi'])]
    poly_subset = poly_df[poly_cols] if poly_cols else pd.DataFrame(index=poly_df.index)

    # Combine all features
    combined = pseudobulk.join(ct_df, how='inner')
    combined = combined.join(poly_subset, how='inner')
    combined = combined.fillna(combined.mean())

    # Keep only numeric and remove zero-variance columns
    combined = combined.select_dtypes(include=[np.number])
    combined = combined.loc[:, combined.var() > 0]
    print(f"  Multi-modal features: {combined.shape}")

    # PCA
    scaled, _ = safe_standardize(combined)
    n_components = min(len(combined) - 1, 5)
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(scaled)

    # Create figure
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

    # Panel A: PC1 vs PC2
    ax = fig.add_subplot(gs[0, 0])
    for i, sample in enumerate(combined.index):
        meta = PDAC_METADATA[sample]
        ax.scatter(pcs[i, 0], pcs[i, 1],
                   c=RESPONSE_COLORS[meta['response']],
                   marker=TIMEPOINT_MARKERS[meta['timepoint']],
                   s=150, edgecolors='black', linewidth=1.5, zorder=3)
        ax.annotate(sample, (pcs[i, 0], pcs[i, 1]),
                    textcoords='offset points', xytext=(5, 5), fontsize=9)

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title("A) Multi-modal PC1 vs PC2")

    # Panel B: PC2 vs PC3
    ax = fig.add_subplot(gs[0, 1])
    for i, sample in enumerate(combined.index):
        meta = PDAC_METADATA[sample]
        ax.scatter(pcs[i, 1], pcs[i, 2] if n_components > 2 else 0,
                   c=RESPONSE_COLORS[meta['response']],
                   marker=TIMEPOINT_MARKERS[meta['timepoint']],
                   s=150, edgecolors='black', linewidth=1.5, zorder=3)

    ax.set_xlabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_ylabel(f"PC3 ({pca.explained_variance_ratio_[2]*100:.1f}%)" if n_components > 2 else "PC3")
    ax.set_title("B) PC2 vs PC3")

    # Panel C: Variance explained
    ax = fig.add_subplot(gs[0, 2])
    var_exp = pca.explained_variance_ratio_
    ax.bar(range(1, len(var_exp) + 1), var_exp * 100, color='steelblue')
    ax.set_xlabel("Principal Component")
    ax.set_ylabel("Variance Explained (%)")
    ax.set_title("C) Variance Explained")

    # Panel D: Feature type contributions
    ax = fig.add_subplot(gs[1, 0])

    # Categorize features
    gene_loadings = np.abs(pca.components_[0, :len(hvg)]).sum()
    ct_loadings = np.abs(pca.components_[0, len(hvg):len(hvg)+len(ct_df.columns)]).sum()
    poly_loadings = np.abs(pca.components_[0, len(hvg)+len(ct_df.columns):]).sum()

    contributions = [gene_loadings, ct_loadings, poly_loadings]
    labels = [f'Expression\n({len(hvg)} genes)',
              f'Cell Types\n({len(ct_df.columns)})',
              f'Spatial/Topology\n({len(poly_cols)})']
    colors = ['#3498db', '#2ecc71', '#9b59b6']

    ax.pie(contributions, labels=labels, autopct='%1.1f%%', colors=colors)
    ax.set_title("D) Feature Type Contribution to PC1")

    # Panel E: Top features heatmap
    ax = fig.add_subplot(gs[1, 1:])
    loadings = pd.Series(pca.components_[0], index=combined.columns)
    top_features = loadings.abs().nlargest(20).index

    # Create sample x feature heatmap
    top_df = combined[top_features]
    top_scaled = (top_df - top_df.mean()) / top_df.std()

    sns.heatmap(top_scaled, cmap='RdBu_r', center=0, ax=ax,
                yticklabels=[f"{s} ({PDAC_METADATA[s]['response']})" for s in top_scaled.index],
                xticklabels=True)
    ax.set_title("E) Top 20 Features (Standardized)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)

    plt.suptitle("Figure 20: Multi-modal Integrated PCA\n(Expression + Cell Types + Spatial/Topology)",
                 fontsize=14, fontweight='bold')
    plt.savefig(FIG_DIR / "fig20_multimodal_pca.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig20_multimodal_pca.pdf", bbox_inches='tight')
    plt.close()
    print(f"  Saved fig20_multimodal_pca.png/pdf")


def fig21_separation_methods():
    """Figure 21: Comparison of different separation methods."""
    print("Generating Figure 21: Separation Methods Comparison")

    samples = load_all_samples()
    if len(samples) < 2:
        print("  ERROR: Not enough samples")
        return

    # Prepare combined feature matrix
    pseudobulk, hvg = create_pseudobulk(samples, n_hvg=100)
    ct_df = load_cell_type_proportions()
    poly_df = load_polymathic_features()

    poly_cols = [c for c in poly_df.columns if any(x in c for x in
                 ['centrality', 'topology', 'betti'])]

    combined = pseudobulk.join(ct_df, how='inner')
    if poly_cols:
        combined = combined.join(poly_df[poly_cols], how='inner')
    combined = combined.fillna(combined.mean())
    combined = combined.select_dtypes(include=[np.number])
    combined = combined.loc[:, combined.var() > 0]

    # Labels
    y = np.array([1 if PDAC_METADATA[s]['response'] == 'R' else 0 for s in combined.index])

    # Standardize
    X, _ = safe_standardize(combined)

    # Methods to compare
    results = []

    # 1. PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)
    sil_pca = silhouette_score(X_pca, y) if len(np.unique(y)) > 1 else 0
    results.append({'method': 'PCA', 'silhouette': sil_pca, 'embedding': X_pca})

    # 2. LDA
    try:
        lda = LinearDiscriminantAnalysis(n_components=1)
        X_lda = lda.fit_transform(X, y)
        # For 1D, pad with zeros for plotting
        X_lda_2d = np.column_stack([X_lda, np.zeros(len(X_lda))])
        sil_lda = silhouette_score(X_lda, y) if len(np.unique(y)) > 1 else 0
        results.append({'method': 'LDA', 'silhouette': sil_lda, 'embedding': X_lda_2d})
    except Exception as e:
        print(f"  LDA failed: {e}")

    # 3. Leave-one-out classification accuracy
    methods_acc = {}

    # Logistic Regression
    loo = LeaveOneOut()
    preds_lr = []
    for train_idx, test_idx in loo.split(X):
        lr = LogisticRegression(max_iter=1000, random_state=42)
        lr.fit(X[train_idx], y[train_idx])
        preds_lr.append(lr.predict(X[test_idx])[0])
    methods_acc['Logistic Regression'] = accuracy_score(y, preds_lr)

    # Random Forest
    preds_rf = []
    for train_idx, test_idx in loo.split(X):
        rf = RandomForestClassifier(n_estimators=100, random_state=42)
        rf.fit(X[train_idx], y[train_idx])
        preds_rf.append(rf.predict(X[test_idx])[0])
    methods_acc['Random Forest'] = accuracy_score(y, preds_rf)

    # SVM
    preds_svm = []
    for train_idx, test_idx in loo.split(X):
        svm = SVC(kernel='linear', random_state=42)
        svm.fit(X[train_idx], y[train_idx])
        preds_svm.append(svm.predict(X[test_idx])[0])
    methods_acc['SVM (Linear)'] = accuracy_score(y, preds_svm)

    # Create figure
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Panel A: PCA embedding
    ax = axes[0, 0]
    for i, sample in enumerate(combined.index):
        meta = PDAC_METADATA[sample]
        ax.scatter(X_pca[i, 0], X_pca[i, 1],
                   c=RESPONSE_COLORS[meta['response']],
                   s=150, edgecolors='black', linewidth=1.5)
        ax.annotate(sample, (X_pca[i, 0], X_pca[i, 1]),
                    textcoords='offset points', xytext=(3, 3), fontsize=9)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"A) PCA (silhouette={sil_pca:.2f})")

    # Panel B: LDA embedding
    ax = axes[0, 1]
    if len(results) > 1:
        X_lda_2d = results[1]['embedding']
        for i, sample in enumerate(combined.index):
            meta = PDAC_METADATA[sample]
            jitter = np.random.normal(0, 0.1)
            ax.scatter(X_lda_2d[i, 0], jitter,
                       c=RESPONSE_COLORS[meta['response']],
                       s=150, edgecolors='black', linewidth=1.5)
            ax.annotate(sample, (X_lda_2d[i, 0], jitter),
                        textcoords='offset points', xytext=(3, 3), fontsize=9)
        ax.set_xlabel("LD1")
        ax.set_ylabel("(jittered)")
        ax.set_title(f"B) LDA (silhouette={results[1]['silhouette']:.2f})")
    else:
        ax.text(0.5, 0.5, "LDA not available", ha='center', va='center')
        ax.set_title("B) LDA")

    # Panel C: Classification accuracy comparison
    ax = axes[0, 2]
    methods = list(methods_acc.keys())
    accs = list(methods_acc.values())
    colors = ['#3498db', '#2ecc71', '#e74c3c']
    bars = ax.bar(methods, accs, color=colors)
    ax.set_ylabel("Leave-One-Out Accuracy")
    ax.set_title("C) Classification Accuracy")
    ax.set_ylim(0, 1.1)
    for bar, acc in zip(bars, accs):
        ax.text(bar.get_x() + bar.get_width()/2, acc + 0.02, f'{acc:.2f}',
                ha='center', fontsize=10)
    ax.axhline(y=0.5, color='gray', linestyle='--', label='Random chance')
    plt.setp(ax.get_xticklabels(), rotation=20, ha='right')

    # Panel D: Feature importance (Random Forest)
    ax = axes[1, 0]
    rf_full = RandomForestClassifier(n_estimators=100, random_state=42)
    rf_full.fit(X, y)
    importances = pd.Series(rf_full.feature_importances_, index=combined.columns)
    top_imp = importances.nlargest(15)

    ax.barh(range(len(top_imp)), top_imp.values, color='#2ecc71')
    ax.set_yticks(range(len(top_imp)))
    ax.set_yticklabels(top_imp.index, fontsize=8)
    ax.set_xlabel("Feature Importance")
    ax.set_title("D) RF Feature Importance (Top 15)")

    # Panel E: Confusion matrix for best method
    ax = axes[1, 1]
    best_method = max(methods_acc, key=methods_acc.get)
    if best_method == 'Logistic Regression':
        preds = preds_lr
    elif best_method == 'Random Forest':
        preds = preds_rf
    else:
        preds = preds_svm

    from sklearn.metrics import confusion_matrix
    cm = confusion_matrix(y, preds)
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax,
                xticklabels=['NR', 'R'], yticklabels=['NR', 'R'])
    ax.set_xlabel("Predicted")
    ax.set_ylabel("Actual")
    ax.set_title(f"E) Confusion Matrix ({best_method})")

    # Panel F: Summary statistics
    ax = axes[1, 2]
    ax.axis('off')

    summary_text = f"""
    SEPARATION METHODS SUMMARY
    {'='*35}

    Dataset: n={len(y)} samples
    Classes: {sum(y)} Responders, {len(y)-sum(y)} Non-Responders
    Features: {X.shape[1]}

    EMBEDDING METHODS:
    - PCA Silhouette: {sil_pca:.3f}
    {'- LDA Silhouette: ' + f"{results[1]['silhouette']:.3f}" if len(results) > 1 else '- LDA: N/A'}

    CLASSIFICATION (LOO-CV):
    - Logistic Regression: {methods_acc['Logistic Regression']:.1%}
    - Random Forest: {methods_acc['Random Forest']:.1%}
    - SVM (Linear): {methods_acc['SVM (Linear)']:.1%}

    Best method: {best_method}
    Accuracy: {max(methods_acc.values()):.1%}

    Note: Small sample size (n={len(y)}) limits
    statistical power. Results are exploratory.
    """

    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='#f8f9fa', edgecolor='#dee2e6'))

    plt.suptitle("Figure 21: Separation Methods Comparison\n(R vs NR Classification)",
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(FIG_DIR / "fig21_separation_methods.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig21_separation_methods.pdf", bbox_inches='tight')
    plt.close()
    print(f"  Saved fig21_separation_methods.png/pdf")


def main():
    """Run all PCA analysis figures."""
    print("=" * 60)
    print("PCA ANALYSES FOR SPATIAL HACKATHON 2026")
    print("=" * 60)

    fig16_bulk_rna_pca()
    fig17_spatial_pca()
    fig18_ncrna_pca()
    fig19_ncrna_integrated_pca()
    fig20_multimodal_pca()
    fig21_separation_methods()

    print("\n" + "=" * 60)
    print("All PCA figures generated successfully!")
    print(f"Output directory: {FIG_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
