#!/usr/bin/env python3
"""
G4X Multimodal Integration - Spatial Hackathon 2026
====================================================
Integrates RNA (387 genes) + Protein (17 markers) data

Methods:
1. RNA-Protein correlation analysis
2. Joint embedding (CCA-style or concatenated PCA)
3. Protein-guided cell type refinement
"""

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import CCA
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Configuration
ANNOTATED_DIR = Path("/home/user/spatial-hackathon-2026/results/g4x_choi_batch2/annotated")
OUTPUT_DIR = Path("/home/user/spatial-hackathon-2026/results/g4x_choi_batch2/multimodal")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR = Path("/home/user/spatial-hackathon-2026/figures/g4x")
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# Key RNA-protein pairs for correlation analysis
RNA_PROTEIN_PAIRS = [
    ('CD3D', 'CD3'),  # T cell
    ('CD3E', 'CD3'),
    ('CD4', 'CD4'),   # Helper T
    ('CD8A', 'CD8'),  # Cytotoxic T
    ('MS4A1', 'CD20'),  # B cell (CD20)
    ('CD68', 'CD68'),   # Macrophage
    ('ITGAX', 'CD11c'), # Dendritic cell
    ('PTPRC', 'CD45'),  # Pan-immune
    ('FOXP3', 'FOXP3'), # Treg
    ('PDCD1', 'PD1'),   # Checkpoint
    ('CD274', 'PDL1'),  # Checkpoint
    ('PECAM1', 'CD31'), # Endothelial
    ('EPCAM', 'PanCK'), # Epithelial (proxy)
    ('KRT19', 'PanCK'), # Epithelial
    ('ACTA2', 'aSMA'),  # Fibroblast/myofibroblast
    ('MKI67', 'KI67'),  # Proliferation
]

def get_protein_names(adata):
    """Get normalized protein names."""
    if 'protein_names' not in adata.uns:
        return None
    names = adata.uns['protein_names']
    # Convert to list if numpy array
    if hasattr(names, 'tolist'):
        return names.tolist()
    return list(names)

def compute_rna_protein_correlation(adata, pairs=RNA_PROTEIN_PAIRS):
    """Compute correlation between RNA and protein expression."""
    protein_names = get_protein_names(adata)
    if protein_names is None:
        return None

    protein_data = adata.obsm['protein']
    protein_name_to_idx = {p: i for i, p in enumerate(protein_names)}

    # Use log-normalized expression
    if 'counts' in adata.layers:
        # Re-normalize for correlation
        X = adata.X.copy() if not hasattr(adata.X, 'toarray') else adata.X.toarray()
    else:
        X = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()

    gene_names = list(adata.var_names)
    gene_name_to_idx = {g: i for i, g in enumerate(gene_names)}

    correlations = []
    for rna_gene, protein_marker in pairs:
        # Check if both exist
        rna_gene_upper = rna_gene.upper()
        matches = [g for g in gene_names if g.upper() == rna_gene_upper]

        if not matches:
            continue
        if protein_marker not in protein_name_to_idx:
            continue

        rna_idx = gene_name_to_idx[matches[0]]
        prot_idx = protein_name_to_idx[protein_marker]

        rna_expr = X[:, rna_idx]
        prot_expr = protein_data[:, prot_idx]

        # Spearman correlation
        r, p = stats.spearmanr(rna_expr, prot_expr)

        correlations.append({
            'RNA_gene': matches[0],
            'Protein': protein_marker,
            'spearman_r': r,
            'p_value': p,
            'n_cells': len(rna_expr)
        })

    return pd.DataFrame(correlations)

def joint_embedding(adata, n_components=20):
    """Create joint RNA + Protein embedding."""
    protein_names = get_protein_names(adata)
    if protein_names is None:
        return adata

    # Get normalized protein data
    protein_data = StandardScaler().fit_transform(adata.obsm['protein'])

    # Get RNA PCA
    if 'X_pca' not in adata.obsm:
        return adata

    rna_pca = adata.obsm['X_pca'][:, :min(20, adata.obsm['X_pca'].shape[1])]

    # Scale RNA PCA
    rna_scaled = StandardScaler().fit_transform(rna_pca)

    # Option 1: Concatenated PCA (simpler, more robust)
    combined = np.hstack([rna_scaled, protein_data])
    joint_pca = PCA(n_components=min(n_components, combined.shape[1] - 1))
    adata.obsm['X_joint_pca'] = joint_pca.fit_transform(combined)

    # Store explained variance
    adata.uns['joint_pca_variance_ratio'] = joint_pca.explained_variance_ratio_

    # Create joint UMAP
    temp_adata = sc.AnnData(np.zeros((adata.n_obs, 1)))
    temp_adata.obsm['X_pca'] = adata.obsm['X_joint_pca']
    sc.pp.neighbors(temp_adata, use_rep='X_pca', n_neighbors=15)
    sc.tl.umap(temp_adata)
    adata.obsm['X_joint_umap'] = temp_adata.obsm['X_umap']

    return adata

def plot_rna_protein_scatter(adata, rna_gene, protein_marker, output_path, sample_frac=0.1):
    """Plot RNA vs Protein scatter for a gene pair."""
    protein_names = get_protein_names(adata)
    if protein_names is None:
        return

    gene_names = list(adata.var_names)
    rna_matches = [g for g in gene_names if g.upper() == rna_gene.upper()]

    if not rna_matches or protein_marker not in protein_names:
        return

    rna_idx = gene_names.index(rna_matches[0])
    prot_idx = protein_names.index(protein_marker)

    X = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
    rna_expr = X[:, rna_idx]
    prot_expr = adata.obsm['protein'][:, prot_idx]

    # Subsample for plotting
    n_sample = min(int(len(rna_expr) * sample_frac), 10000)
    idx = np.random.choice(len(rna_expr), n_sample, replace=False)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(rna_expr[idx], prot_expr[idx], alpha=0.3, s=5, c='steelblue')
    ax.set_xlabel(f'{rna_matches[0]} (RNA)')
    ax.set_ylabel(f'{protein_marker} (Protein)')

    r, p = stats.spearmanr(rna_expr, prot_expr)
    ax.set_title(f'RNA-Protein Correlation\nρ = {r:.3f}, p = {p:.2e}')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()

def plot_joint_umap(adata, output_path, color_by='cell_type'):
    """Plot joint UMAP colored by cell type."""
    if 'X_joint_umap' not in adata.obsm:
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Original UMAP
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color=color_by, ax=axes[0], show=False, legend_loc='right margin')
        axes[0].set_title('RNA-only UMAP')

    # Joint UMAP
    adata_temp = adata.copy()
    adata_temp.obsm['X_umap'] = adata.obsm['X_joint_umap']
    sc.pl.umap(adata_temp, color=color_by, ax=axes[1], show=False, legend_loc='right margin')
    axes[1].set_title('Joint RNA+Protein UMAP')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

def process_sample(sample_path):
    """Process a single sample for multimodal integration."""
    sample_id = sample_path.stem.replace('_annotated', '')
    print(f"\n--- Processing {sample_id} ---")

    # Load annotated data
    adata = sc.read_h5ad(sample_path)
    print(f"  Loaded: {adata.n_obs} cells")

    # Compute RNA-Protein correlations
    print("  Computing RNA-Protein correlations...")
    corr_df = compute_rna_protein_correlation(adata)
    if corr_df is not None and len(corr_df) > 0:
        print(f"  Found {len(corr_df)} correlatable pairs")
        corr_df.to_csv(OUTPUT_DIR / f"{sample_id}_rna_protein_corr.csv", index=False)

        # Top correlations
        top_corr = corr_df.nlargest(5, 'spearman_r')
        for _, row in top_corr.iterrows():
            print(f"    {row['RNA_gene']}-{row['Protein']}: ρ={row['spearman_r']:.3f}")

    # Joint embedding
    print("  Creating joint embedding...")
    adata = joint_embedding(adata)

    if 'X_joint_pca' in adata.obsm:
        print(f"  Joint PCA: {adata.obsm['X_joint_pca'].shape[1]} components")

    # Save with multimodal features
    out_path = OUTPUT_DIR / f"{sample_id}_multimodal.h5ad"
    adata.write(out_path)
    print(f"  Saved: {out_path}")

    # Generate figures for first sample
    if sample_id == 'A01':
        print("  Generating example figures...")

        # RNA-Protein scatter for CD3
        plot_rna_protein_scatter(adata, 'CD3E', 'CD3',
                                FIGURE_DIR / f"{sample_id}_CD3_correlation.png")

        # Joint UMAP comparison
        if 'cell_type' in adata.obs.columns:
            plot_joint_umap(adata, FIGURE_DIR / f"{sample_id}_joint_umap.png")

    return {
        'sample_id': sample_id,
        'n_cells': adata.n_obs,
        'n_rna_prot_pairs': len(corr_df) if corr_df is not None else 0,
        'has_joint_pca': 'X_joint_pca' in adata.obsm,
        'has_joint_umap': 'X_joint_umap' in adata.obsm
    }

def summarize_correlations(output_dir):
    """Aggregate RNA-Protein correlations across all samples."""
    corr_files = list(output_dir.glob("*_rna_protein_corr.csv"))
    if not corr_files:
        return None

    all_corrs = []
    for f in corr_files:
        df = pd.read_csv(f)
        df['sample_id'] = f.stem.replace('_rna_protein_corr', '')
        all_corrs.append(df)

    combined = pd.concat(all_corrs, ignore_index=True)

    # Average correlation per pair
    avg_corr = combined.groupby(['RNA_gene', 'Protein']).agg({
        'spearman_r': ['mean', 'std'],
        'n_cells': 'sum'
    }).round(3)
    avg_corr.columns = ['mean_r', 'std_r', 'total_cells']
    avg_corr = avg_corr.sort_values('mean_r', ascending=False)

    avg_corr.to_csv(output_dir / "aggregated_rna_protein_correlations.csv")
    return avg_corr

def main():
    """Main multimodal integration pipeline."""
    print("=" * 60)
    print("G4X Multimodal Integration")
    print("=" * 60)

    # Get annotated samples
    sample_files = sorted(ANNOTATED_DIR.glob("*_annotated.h5ad"))
    print(f"Found {len(sample_files)} annotated samples")

    # Process each sample
    results = []
    for sample_path in tqdm(sample_files, desc="Processing samples"):
        result = process_sample(sample_path)
        results.append(result)

    # Aggregate correlations
    print("\n--- Aggregating RNA-Protein Correlations ---")
    avg_corr = summarize_correlations(OUTPUT_DIR)
    if avg_corr is not None:
        print("Top RNA-Protein Correlations (averaged across samples):")
        print(avg_corr.head(10))

    # Save summary
    summary_df = pd.DataFrame(results)
    summary_df.to_csv(OUTPUT_DIR / "multimodal_summary.csv", index=False)
    print(f"\nSummary saved to {OUTPUT_DIR / 'multimodal_summary.csv'}")

    print("\nMultimodal integration complete!")
    return results

if __name__ == "__main__":
    main()
