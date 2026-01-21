#!/usr/bin/env python3
"""
G4X MULTIMODAL ANALYSIS: RNA + Protein Spatial Transcriptomics
===============================================================

Analyzes Singular Genomics G4X data with both RNA and protein measurements.
Data: Gastric Cancer samples from Choi lab

Features:
1. Data loading and QC
2. RNA-Protein correlation analysis
3. Cell type annotation using protein markers
4. Polymathic analysis (graph centrality, TDA)
5. Cross-modal integration

Author: Max Van Belkum
Date: 2026-01-19
"""

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import anndata as ad
import h5py
import gzip
import logging
from scipy.stats import spearmanr, pearsonr
from sklearn.preprocessing import StandardScaler

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

# Paths
DATA_DIR = Path('/home/user/data/hackathon/multimodal')
OUTPUT_DIR = Path('/home/user/spatial-hackathon-2026/outputs')
ADATA_DIR = OUTPUT_DIR / 'adata' / 'g4x'
FIG_DIR = OUTPUT_DIR / 'figures' / 'g4x'

ADATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
(OUTPUT_DIR / 'tables').mkdir(parents=True, exist_ok=True)

# G4X Samples
G4X_SAMPLES = ['A01', 'B01']

# Protein panel markers and cell type associations
PROTEIN_MARKERS = {
    'CD3': 'T cells',
    'CD4': 'CD4+ T cells',
    'CD8': 'CD8+ T cells',
    'CD20': 'B cells',
    'CD45': 'Immune cells',
    'CD68': 'Macrophages',
    'FOXP3': 'Regulatory T cells',
    'PD1': 'Exhausted T cells',
    'PDL1': 'Tumor/APC',
    'PanCK': 'Epithelial/Tumor',
    'aSMA': 'CAFs/Smooth muscle',
    'CD31': 'Endothelial',
    'HLA-DR': 'Antigen presenting cells',
    'KI67': 'Proliferating cells',
    'CD11c': 'Dendritic cells',
}


def load_g4x_sample(sample: str) -> Optional[ad.AnnData]:
    """
    Load G4X multimodal data for a sample.

    Returns AnnData with:
    - .X: RNA expression
    - .obsm['protein']: Protein expression
    - .obsm['spatial']: Spatial coordinates
    """
    sample_dir = DATA_DIR / sample / 'single_cell_data'

    if not sample_dir.exists():
        logger.warning(f"Sample directory not found: {sample_dir}")
        return None

    logger.info(f"Loading {sample}...")

    # Load cell metadata
    meta_path = sample_dir / 'cell_metadata.csv.gz'
    with gzip.open(meta_path, 'rt') as f:
        cell_meta = pd.read_csv(f, index_col=0)
    logger.info(f"  Loaded {len(cell_meta)} cells metadata")

    # Load RNA expression
    rna_path = sample_dir / 'cell_by_transcript.csv.gz'
    with gzip.open(rna_path, 'rt') as f:
        rna_df = pd.read_csv(f, index_col=0)
    logger.info(f"  Loaded RNA: {rna_df.shape[0]} cells x {rna_df.shape[1]} genes")

    # Load protein expression
    protein_path = sample_dir / 'cell_by_protein.csv.gz'
    with gzip.open(protein_path, 'rt') as f:
        protein_df = pd.read_csv(f, index_col=0)
    logger.info(f"  Loaded Protein: {protein_df.shape[0]} cells x {protein_df.shape[1]} markers")

    # Align indices
    common_cells = cell_meta.index.intersection(rna_df.index).intersection(protein_df.index)
    logger.info(f"  Common cells: {len(common_cells)}")

    cell_meta = cell_meta.loc[common_cells]
    rna_df = rna_df.loc[common_cells]
    protein_df = protein_df.loc[common_cells]

    # Create AnnData
    adata = ad.AnnData(X=rna_df.values, obs=cell_meta, var=pd.DataFrame(index=rna_df.columns))
    adata.var_names = rna_df.columns
    adata.obs_names = common_cells

    # Add protein data
    adata.obsm['protein'] = protein_df.values
    adata.uns['protein_names'] = list(protein_df.columns)

    # Add spatial coordinates
    if 'x_centroid' in cell_meta.columns and 'y_centroid' in cell_meta.columns:
        adata.obsm['spatial'] = cell_meta[['x_centroid', 'y_centroid']].values
    elif 'centroid_x' in cell_meta.columns and 'centroid_y' in cell_meta.columns:
        adata.obsm['spatial'] = cell_meta[['centroid_x', 'centroid_y']].values
    elif 'cell_x' in cell_meta.columns and 'cell_y' in cell_meta.columns:
        adata.obsm['spatial'] = cell_meta[['cell_x', 'cell_y']].values
        logger.info(f"  Using cell_x/cell_y for spatial coordinates")
    else:
        logger.warning(f"  No spatial coordinates found in {list(cell_meta.columns)}")

    adata.obs['sample'] = sample

    return adata


def annotate_cells_by_protein(adata: ad.AnnData) -> ad.AnnData:
    """
    Annotate cell types based on protein marker expression.
    Uses threshold-based gating approach.
    """
    protein_df = pd.DataFrame(
        adata.obsm['protein'],
        index=adata.obs_names,
        columns=adata.uns['protein_names']
    )

    # Handle *_intensity_mean naming convention
    # Create a mapping from simple names to actual column names
    col_map = {}
    for col in protein_df.columns:
        simple_name = col.replace('_intensity_mean', '').upper()
        col_map[simple_name] = col

    # Standardize protein expression per marker
    protein_scaled = pd.DataFrame(
        StandardScaler().fit_transform(protein_df),
        index=protein_df.index,
        columns=protein_df.columns
    )

    # Helper function to get marker value
    def get_marker(row, marker_name):
        col = col_map.get(marker_name.upper())
        if col and col in row.index:
            return row[col]
        return -999  # Not found

    # Simple gating logic (vectorized for speed)
    cell_types = []
    for idx in protein_scaled.index:
        row = protein_scaled.loc[idx]

        # Priority-based annotation
        panck = get_marker(row, 'PanCK')
        cd8 = get_marker(row, 'CD8')
        cd4 = get_marker(row, 'CD4')
        foxp3 = get_marker(row, 'FOXP3')
        cd3 = get_marker(row, 'CD3')
        cd20 = get_marker(row, 'CD20')
        cd68 = get_marker(row, 'CD68')
        cd11c = get_marker(row, 'CD11c')
        asma = get_marker(row, 'aSMA')
        cd31 = get_marker(row, 'CD31')
        cd45 = get_marker(row, 'CD45')

        if panck > 1.0:
            cell_type = 'Epithelial'
        elif cd8 > 0.5:
            cell_type = 'CD8+ T cell'
        elif cd4 > 0.5:
            if foxp3 > 0.5:
                cell_type = 'Treg'
            else:
                cell_type = 'CD4+ T cell'
        elif cd3 > 0.5:
            cell_type = 'T cell'
        elif cd20 > 0.5:
            cell_type = 'B cell'
        elif cd68 > 0.5:
            cell_type = 'Macrophage'
        elif cd11c > 0.5:
            cell_type = 'Dendritic cell'
        elif asma > 0.5:
            cell_type = 'Fibroblast'
        elif cd31 > 0.5:
            cell_type = 'Endothelial'
        elif cd45 > 0.5:
            cell_type = 'Other immune'
        else:
            cell_type = 'Other'

        cell_types.append(cell_type)

    adata.obs['protein_cell_type'] = pd.Categorical(cell_types)

    # Count cell types
    counts = adata.obs['protein_cell_type'].value_counts()
    logger.info(f"  Cell type annotation:")
    for ct, count in counts.items():
        logger.info(f"    {ct}: {count} ({100*count/len(adata):.1f}%)")

    return adata


def compute_rna_protein_correlation(adata: ad.AnnData) -> pd.DataFrame:
    """
    Compute correlation between RNA and protein expression.
    Looks for matching gene names.
    """
    protein_names = adata.uns['protein_names']
    gene_names = list(adata.var_names)

    correlations = []

    protein_df = pd.DataFrame(
        adata.obsm['protein'],
        index=adata.obs_names,
        columns=protein_names
    )

    for protein in protein_names:
        # Skip non-gene markers
        if protein in ['Isotype', 'ATPase']:
            continue

        # Try to find matching gene
        gene_matches = [g for g in gene_names if protein.upper() in g.upper()]

        if gene_matches:
            gene = gene_matches[0]
            protein_expr = protein_df[protein].values
            rna_expr = adata[:, gene].X.flatten()

            if hasattr(rna_expr, 'toarray'):
                rna_expr = rna_expr.toarray().flatten()

            # Compute correlation
            r_spearman, p_spearman = spearmanr(rna_expr, protein_expr)
            r_pearson, p_pearson = pearsonr(rna_expr, protein_expr)

            correlations.append({
                'protein': protein,
                'gene': gene,
                'spearman_r': r_spearman,
                'spearman_p': p_spearman,
                'pearson_r': r_pearson,
                'pearson_p': p_pearson
            })

    return pd.DataFrame(correlations)


def analyze_multimodal_spatial(adata: ad.AnnData) -> Dict:
    """
    Perform spatial analysis on multimodal data.
    """
    results = {}

    # Build spatial graph
    if 'spatial' in adata.obsm:
        sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

        # Neighborhood enrichment by protein cell type
        if 'protein_cell_type' in adata.obs.columns:
            sq.gr.nhood_enrichment(adata, cluster_key='protein_cell_type')
            results['nhood_enrichment'] = adata.uns.get('protein_cell_type_nhood_enrichment', None)

    return results


def create_g4x_figures(adata: ad.AnnData, sample: str):
    """Create visualization figures for G4X sample."""

    # Figure 1: Spatial overview with protein cell types
    if 'spatial' in adata.obsm and 'protein_cell_type' in adata.obs.columns:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Panel 1: Spatial by cell type
        ax = axes[0]
        coords = adata.obsm['spatial']
        cell_types = adata.obs['protein_cell_type']

        for ct in cell_types.cat.categories:
            mask = cell_types == ct
            ax.scatter(coords[mask, 0], coords[mask, 1], s=3, alpha=0.5, label=ct)

        ax.set_title(f'{sample}: Cell Types (Protein-based)', fontweight='bold')
        ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1, 1))
        ax.set_aspect('equal')
        ax.axis('off')

        # Panel 2: Protein expression heatmap
        ax = axes[1]
        protein_df = pd.DataFrame(
            adata.obsm['protein'],
            index=adata.obs_names,
            columns=adata.uns['protein_names']
        )

        # Mean expression by cell type
        mean_expr = protein_df.groupby(adata.obs['protein_cell_type']).mean()
        sns.heatmap(mean_expr.T, ax=ax, cmap='viridis', xticklabels=True, yticklabels=True)
        ax.set_title('Protein Expression by Cell Type', fontweight='bold')

        plt.tight_layout()
        plt.savefig(FIG_DIR / f'{sample}_overview.png', dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"  Saved {sample}_overview.png")


def main():
    """Run G4X multimodal analysis."""
    logger.info("="*60)
    logger.info("G4X MULTIMODAL ANALYSIS")
    logger.info("="*60)

    all_adatas = {}
    all_correlations = []

    for sample in G4X_SAMPLES:
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing {sample}")
        logger.info(f"{'='*60}")

        # Load data
        adata = load_g4x_sample(sample)
        if adata is None:
            continue

        # Annotate cells by protein
        adata = annotate_cells_by_protein(adata)

        # RNA-Protein correlation
        corr_df = compute_rna_protein_correlation(adata)
        if not corr_df.empty:
            corr_df['sample'] = sample
            all_correlations.append(corr_df)
            logger.info(f"  RNA-Protein correlations computed for {len(corr_df)} markers")

        # Spatial analysis
        results = analyze_multimodal_spatial(adata)

        # Create figures
        create_g4x_figures(adata, sample)

        # Save processed data
        adata.write(ADATA_DIR / f'{sample}_g4x.h5ad')
        logger.info(f"  Saved {sample}_g4x.h5ad")

        all_adatas[sample] = adata

    # Combined analysis
    if all_correlations:
        combined_corr = pd.concat(all_correlations, ignore_index=True)
        combined_corr.to_csv(OUTPUT_DIR / 'tables' / 'g4x_rna_protein_correlations.csv', index=False)
        logger.info(f"\nSaved RNA-Protein correlations to tables/g4x_rna_protein_correlations.csv")

        # Summary figure
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot correlation by protein marker
        mean_corr = combined_corr.groupby('protein')['spearman_r'].mean().sort_values(ascending=False)
        colors = ['green' if r > 0.3 else 'orange' if r > 0 else 'red' for r in mean_corr.values]
        mean_corr.plot(kind='barh', ax=ax, color=colors)
        ax.axvline(x=0, color='gray', linestyle='--')
        ax.set_xlabel('Spearman Correlation (RNA vs Protein)')
        ax.set_title('RNA-Protein Correlation by Marker\n(Higher = better central dogma compliance)', fontweight='bold')
        ax.invert_yaxis()

        plt.tight_layout()
        plt.savefig(FIG_DIR / 'rna_protein_correlation_summary.png', dpi=300, bbox_inches='tight')
        plt.close()
        logger.info("Saved rna_protein_correlation_summary.png")

    logger.info("\n" + "="*60)
    logger.info("G4X ANALYSIS COMPLETE")
    logger.info("="*60)


if __name__ == '__main__':
    main()
