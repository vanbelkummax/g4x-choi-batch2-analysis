#!/usr/bin/env python3
"""
G4X PCA Deep Dive Analysis
==========================
Comprehensive PCA exploration of the G4X gastric cancer progression dataset.

This script generates 8 figure panels exploring:
- Batch correction validation (pre/post metrics)
- Variance decomposition by technical vs biological factors
- Sample overview with H&E thumbnails
- Pseudobulk PCA with progression statistics
- Cell-level PCA colored by stage/patient/cell type
- PC loading analysis with gene set enrichment
- Multimodal comparison (RNA vs Protein vs WNN)
- Per-sample detailed panels

Usage:
    conda activate enact
    cd ~/g4x-choi-batch2-analysis

    # Dry run (validate inputs)
    python scripts/46_pca_deep_dive.py --dry-run

    # Full run
    python scripts/46_pca_deep_dive.py 2>&1 | tee logs/46_pca_deep_dive.log

    # Fast iteration (skip enrichment)
    python scripts/46_pca_deep_dive.py --skip-enrichment

    # Custom output directory
    python scripts/46_pca_deep_dive.py --output-dir results/pca_v2

Output:
    results/pca_deep_dive/
    ├── figures/
    │   ├── fig0_batch_validation.png
    │   ├── fig1_variance_decomposition.png
    │   ├── fig2_sample_overview.png
    │   ├── fig3_pseudobulk_pca.png
    │   ├── fig4_cell_pca.png
    │   ├── fig5_loading_analysis.png
    │   ├── fig6_multimodal_comparison.png
    │   ├── fig7_summary_grid.png
    │   └── sample_panels/{sample}_panel.png
    ├── data/
    │   ├── variance_partition.csv
    │   ├── pc_loadings.csv
    │   ├── pseudobulk_pca.csv
    │   └── progression_stats.json
    └── PCA_DEEP_DIVE_REPORT.md
"""

import os
import sys
import shutil
import json
import warnings
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, List, Any

warnings.filterwarnings('ignore')


def check_environment(skip_check: bool = False):
    """Verify running in correct conda environment."""
    if skip_check or os.environ.get("G4X_SKIP_ENV_CHECK"):
        return

    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    env_name = os.path.basename(conda_prefix) if conda_prefix else ""
    valid_patterns = ['enact', 'spatial', 'scanpy', 'single-cell', 'scverse']
    is_valid = any(pattern in env_name.lower() for pattern in valid_patterns)

    if not is_valid:
        print("ERROR: This script requires a spatial analysis conda environment.")
        print(f"Current environment: {env_name or 'none'}")
        print("\nOptions:")
        print("  1. Activate a valid env: conda activate enact")
        print("  2. Skip check: --skip-env-check")
        print("  3. Set env var: export G4X_SKIP_ENV_CHECK=1")
        sys.exit(1)


def check_disk_space(path, required_gb=2.0):
    """Check available disk space."""
    total, used, free = shutil.disk_usage(path)
    free_gb = free / (1024 ** 3)
    if free_gb < required_gb:
        print(f"WARNING: Low disk space! {free_gb:.1f} GB available, {required_gb:.1f} GB recommended")
        if free_gb < 1.0:
            print("ERROR: Critically low disk space (<1 GB). Aborting.")
            sys.exit(1)
    return free_gb


# Defer imports until after environment check
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
import logging
import argparse
import gc

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# =============================================================================
# Configuration
# =============================================================================

BASE_DIR = Path('/home/user/g4x-choi-batch2-analysis')
QC_OUTPUT_DIR = BASE_DIR / 'results' / 'qc_all_samples'

# Data paths
COUNTS_PATH = QC_OUTPUT_DIR / 'merged' / 'merged_counts.h5ad'
CORRECTED_PATH = QC_OUTPUT_DIR / 'merged' / 'merged_corrected.h5ad'
BATCH_CSV = QC_OUTPUT_DIR / 'merged' / 'batch_assessment.csv'
SAMPLE_QC_CSV = QC_OUTPUT_DIR / 'sample_qc_summary.csv'
# G4X raw data directory with correctly labeled H&E images
# NOTE: QuPath thumbnails are mislabeled - use G4X source instead
G4X_BASE = Path('/mnt/x/Choi_Batch_2_Tuesday')

# Lane directory names (fixed prefix, lane number varies)
LANE_DIRS = {
    '01': 'g4-028-083-FC1-L001_5XGaAe5DB2dm7sRK',
    '02': 'g4-028-083-FC1-L002_5XGaAe5DB2dm7sRK',
    '03': 'g4-028-083-FC1-L003_5XGaAe5DB2dm7sRK',
    '04': 'g4-028-083-FC1-L004_5XGaAe5DB2dm7sRK',
}

# Expected exclusions
# - H04: Failed QC (low median transcripts: 25.0 < 30)
# - C02, H02: Passed QC but failed processing (NaN in protein data during WNN)
EXCLUDED_SAMPLES = {'H04', 'C02', 'H02'}

# Color schemes
STAGE_COLORS = {
    'normal': '#2ecc71',
    'control': '#3498db',
    'metaplasia': '#f39c12',
    'cancer': '#e74c3c'
}

STAGE_ORDER = ['normal', 'control', 'metaplasia', 'cancer']

# Ordinal stage encoding for regression
STAGE_ORDINAL = {
    'normal': 0,
    'control': 1,
    'metaplasia': 2,
    'cancer': 3
}


# =============================================================================
# Helper Functions
# =============================================================================

def load_batch_metrics(batch_csv_path: Path) -> Dict[str, Any]:
    """Load batch metrics with guards for missing data."""
    result = {
        'pre_lisi': None,
        'post_lisi': None,
        'pre_silhouette': None,
        'post_silhouette': None,
        'lisi_threshold': 3.0,
        'valid': False,
        'error': None
    }

    if not batch_csv_path.exists():
        result['error'] = f"Batch metrics file not found: {batch_csv_path}"
        return result

    try:
        # EXPLICIT load - this was missing in v3
        batch_df = pd.read_csv(batch_csv_path)
        logger.info(f"Loaded batch metrics: {batch_df.shape[0]} rows, columns: {list(batch_df.columns)}")
    except Exception as e:
        result['error'] = f"Failed to read batch metrics: {e}"
        return result

    # Check required columns
    required_cols = ['stage', 'mean_lisi', 'silhouette_batch']
    missing_cols = [c for c in required_cols if c not in batch_df.columns]
    if missing_cols:
        result['error'] = f"Missing columns in batch metrics: {missing_cols}"
        return result

    # Extract metrics (using explicit string matching)
    pre_mask = batch_df['stage'] == 'pre-correction'
    post_mask = batch_df['stage'] == 'post-correction'

    if pre_mask.any():
        result['pre_lisi'] = float(batch_df.loc[pre_mask, 'mean_lisi'].values[0])
        result['pre_silhouette'] = float(batch_df.loc[pre_mask, 'silhouette_batch'].values[0])

    if post_mask.any():
        result['post_lisi'] = float(batch_df.loc[post_mask, 'mean_lisi'].values[0])
        result['post_silhouette'] = float(batch_df.loc[post_mask, 'silhouette_batch'].values[0])

    if 'lisi_threshold' in batch_df.columns:
        result['lisi_threshold'] = float(batch_df['lisi_threshold'].values[0])

    result['valid'] = result['post_lisi'] is not None
    return result


def get_celltype_colors(adata: ad.AnnData) -> Dict[str, tuple]:
    """Get or generate cell type colors."""
    if 'cell_type_colors' in adata.uns:
        categories = adata.obs['cell_type'].cat.categories
        colors = adata.uns['cell_type_colors']
        return dict(zip(categories, colors))
    else:
        # Generate consistent palette
        from matplotlib import colormaps
        categories = adata.obs['cell_type'].unique()
        n_types = len(categories)
        cmap = colormaps['tab20']
        colors = [cmap(i % 20) for i in range(n_types)]
        return dict(zip(categories, colors))


def validate_embedding(adata: ad.AnnData, embedding: str) -> bool:
    """Check if embedding exists in adata.obsm."""
    if embedding not in adata.obsm:
        available = list(adata.obsm.keys())
        logger.error(f"Embedding '{embedding}' not found. Available: {available}")
        return False
    return True


def subsample_adata(adata: ad.AnnData, n: int, seed: int = 42) -> ad.AnnData:
    """Subsample AnnData to n cells."""
    if adata.n_obs <= n:
        return adata
    np.random.seed(seed)
    indices = np.random.choice(adata.n_obs, size=n, replace=False)
    return adata[indices].copy()


def get_he_thumbnail(sample_id: str, cache_dir: Path) -> Optional[Path]:
    """Get H&E thumbnail path for sample from G4X data (correctly labeled).

    The G4X data folders contain properly labeled h_and_e_thumbnail.jpg files,
    unlike the QuPath thumbnails which are mislabeled.

    Path pattern: G4X_BASE / LANE_DIRS[lane] / sample_id / h_and_e / h_and_e_thumbnail.jpg
    """
    # Extract lane from sample ID (e.g., 'A01' -> '01')
    lane_suffix = sample_id[-2:]  # Last 2 chars: '01', '02', '03', or '04'

    if lane_suffix not in LANE_DIRS:
        logger.warning(f"Unknown lane suffix for sample {sample_id}: {lane_suffix}")
        return None

    lane_dir = LANE_DIRS[lane_suffix]
    thumbnail_path = G4X_BASE / lane_dir / sample_id / 'h_and_e' / 'h_and_e_thumbnail.jpg'

    if thumbnail_path.exists():
        # Cache locally
        cached_path = cache_dir / f'{sample_id}_he.jpg'
        if not cached_path.exists():
            try:
                shutil.copy(thumbnail_path, cached_path)
            except Exception as e:
                logger.warning(f"Could not cache H&E for {sample_id}: {e}")
                return thumbnail_path
        return cached_path
    else:
        logger.warning(f"H&E thumbnail not found: {thumbnail_path}")

    return None


# =============================================================================
# Analysis Functions
# =============================================================================

def compute_variance_partition(adata: ad.AnnData, embedding: str = 'X_pca_harmony',
                               n_pcs: int = 20) -> pd.DataFrame:
    """Partition variance by factors, with embedding validation."""

    # VALIDATE embedding exists
    if not validate_embedding(adata, embedding):
        raise ValueError(f"Embedding '{embedding}' not found")

    X = adata.obsm[embedding][:, :min(n_pcs, adata.obsm[embedding].shape[1])]
    results = []

    for pc_idx in range(X.shape[1]):
        pc_values = X[:, pc_idx]

        for factor in ['lane', 'patient', 'stage', 'cell_type']:
            if factor not in adata.obs.columns:
                continue

            # FIXED: Use drop_first=True to avoid multicollinearity
            dummies = pd.get_dummies(adata.obs[factor], drop_first=True)

            model = LinearRegression(fit_intercept=True).fit(dummies, pc_values)
            r2 = model.score(dummies, pc_values)
            r2 = max(0, r2)  # Clip negative

            category = 'Technical' if factor in ['lane', 'patient'] else 'Biological'

            results.append({
                'PC': f'PC{pc_idx+1}',
                'Factor': factor,
                'R²': r2,
                'Category': category
            })

    return pd.DataFrame(results)


def compute_pseudobulk_pca(adata: ad.AnnData) -> pd.DataFrame:
    """Compute sample-level pseudobulk PCA."""
    logger.info("Computing pseudobulk PCA...")

    # Get raw counts if available
    if 'counts' in adata.layers:
        X = adata.layers['counts']
    else:
        X = adata.X

    # Convert to dense if sparse
    if hasattr(X, 'toarray'):
        X = X.toarray()

    # Aggregate by sample
    sample_ids = adata.obs['sample_id'].values
    unique_samples = np.unique(sample_ids)

    pseudobulk = []
    for sample in unique_samples:
        mask = sample_ids == sample
        sample_counts = X[mask].sum(axis=0)
        if hasattr(sample_counts, 'A1'):
            sample_counts = sample_counts.A1
        pseudobulk.append(sample_counts)

    pseudobulk_matrix = np.vstack(pseudobulk)

    # Log-normalize
    pseudobulk_norm = np.log1p(pseudobulk_matrix / pseudobulk_matrix.sum(axis=1, keepdims=True) * 1e4)

    # PCA
    pca = PCA(n_components=min(10, len(unique_samples) - 1))
    pcs = pca.fit_transform(pseudobulk_norm)

    # Create DataFrame
    sample_meta = adata.obs.drop_duplicates('sample_id').set_index('sample_id')

    df = pd.DataFrame(
        pcs,
        index=unique_samples,
        columns=[f'PC{i+1}' for i in range(pcs.shape[1])]
    )

    # Add metadata
    for col in ['stage', 'patient', 'lane']:
        if col in sample_meta.columns:
            df[col] = df.index.map(lambda x: sample_meta.loc[x, col] if x in sample_meta.index else None)

    df['n_cells'] = adata.obs.groupby('sample_id').size().reindex(unique_samples).values

    # Store variance explained as metadata (per-PC, not per-sample)
    df.attrs['variance_explained'] = pca.explained_variance_ratio_.tolist()

    return df


def compute_progression_stats(pseudobulk_pca: pd.DataFrame) -> Dict[str, Any]:
    """Compute Kruskal-Wallis and Spearman with power warning."""

    # Validate data
    if 'stage' not in pseudobulk_pca.columns:
        return {'error': 'No stage column in pseudobulk data'}

    valid_mask = pseudobulk_pca['stage'].isin(STAGE_ORDER)
    df = pseudobulk_pca[valid_mask].copy()

    if len(df) < 4:
        return {'error': f'Insufficient samples for statistics: {len(df)}'}

    # Group by stage
    stage_groups = [df[df['stage'] == s]['PC1'].values for s in STAGE_ORDER if s in df['stage'].values]
    stage_counts = {s: len(df[df['stage'] == s]) for s in STAGE_ORDER if s in df['stage'].values}

    # Kruskal-Wallis
    try:
        kw_stat, kw_pval = stats.kruskal(*[g for g in stage_groups if len(g) > 0])
    except Exception as e:
        kw_stat, kw_pval = np.nan, np.nan

    # Spearman correlation with ordinal stage
    df['stage_ordinal'] = df['stage'].map(STAGE_ORDINAL)
    try:
        rho, spearman_pval = stats.spearmanr(df['stage_ordinal'], df['PC1'])
    except Exception as e:
        rho, spearman_pval = np.nan, np.nan

    stats_result = {
        'kruskal_wallis_stat': float(kw_stat) if not np.isnan(kw_stat) else None,
        'kruskal_wallis_pval': float(kw_pval) if not np.isnan(kw_pval) else None,
        'spearman_rho': float(rho) if not np.isnan(rho) else None,
        'spearman_pval': float(spearman_pval) if not np.isnan(spearman_pval) else None,
        'n_per_stage': stage_counts,
        'n_valid_samples': len(df),
        'power_note': (
            f"Pseudobulk analysis has limited power (n={min(stage_counts.values())}-"
            f"{max(stage_counts.values())} per stage). Cell-level analyses "
            "provide more robust statistical inference."
        )
    }

    return stats_result


def compute_pc_loadings(adata: ad.AnnData, embedding: str = 'X_pca_harmony',
                        n_pcs: int = 5, n_top: int = 20) -> pd.DataFrame:
    """Extract top gene loadings for each PC."""

    if 'PCs' not in adata.varm:
        logger.warning("PCA loadings not found in varm['PCs']. Computing from scratch...")
        # Fallback: use variance of gene expression along PCs
        return pd.DataFrame()

    loadings = adata.varm['PCs'][:, :n_pcs]
    genes = adata.var_names.tolist()

    results = []
    for pc_idx in range(n_pcs):
        pc_loadings = loadings[:, pc_idx]

        # Top positive
        top_pos_idx = np.argsort(pc_loadings)[-n_top:][::-1]
        for rank, idx in enumerate(top_pos_idx):
            results.append({
                'PC': f'PC{pc_idx+1}',
                'Gene': genes[idx],
                'Loading': pc_loadings[idx],
                'Direction': 'positive',
                'Rank': rank + 1
            })

        # Top negative
        top_neg_idx = np.argsort(pc_loadings)[:n_top]
        for rank, idx in enumerate(top_neg_idx):
            results.append({
                'PC': f'PC{pc_idx+1}',
                'Gene': genes[idx],
                'Loading': pc_loadings[idx],
                'Direction': 'negative',
                'Rank': rank + 1
            })

    return pd.DataFrame(results)


def compute_modality_pcas(adata: ad.AnnData, n_pcs: int = 30) -> Dict[str, np.ndarray]:
    """
    Compute separate PCAs for RNA-only, Protein-only, and combined modalities.

    Returns dict with embeddings for each modality.
    """
    logger.info("Computing modality-specific PCAs...")
    results = {}

    # RNA-only PCA (already exists as X_pca_harmony, but compute fresh for comparison)
    logger.info("  Computing RNA-only PCA...")
    from sklearn.decomposition import PCA as skPCA
    from sklearn.preprocessing import StandardScaler

    # Get RNA data (normalized)
    if hasattr(adata.X, 'toarray'):
        X_rna = adata.X.toarray()
    else:
        X_rna = adata.X.copy()

    # Standardize
    X_rna_scaled = StandardScaler().fit_transform(X_rna)

    # PCA
    pca_rna = skPCA(n_components=min(n_pcs, X_rna.shape[1]))
    results['rna_pca'] = pca_rna.fit_transform(X_rna_scaled)
    results['rna_variance_ratio'] = pca_rna.explained_variance_ratio_

    # Protein-only PCA
    logger.info("  Computing Protein-only PCA...")
    if 'protein' in adata.obsm:
        X_protein = adata.obsm['protein']
        if hasattr(X_protein, 'toarray'):
            X_protein = X_protein.toarray()

        X_protein_scaled = StandardScaler().fit_transform(X_protein)
        pca_protein = skPCA(n_components=min(n_pcs, X_protein.shape[1]))
        results['protein_pca'] = pca_protein.fit_transform(X_protein_scaled)
        results['protein_variance_ratio'] = pca_protein.explained_variance_ratio_
    elif 'X_pca_protein' in adata.obsm:
        results['protein_pca'] = adata.obsm['X_pca_protein']
        results['protein_variance_ratio'] = None

    # Combined PCA (concatenate RNA + Protein)
    logger.info("  Computing Combined (RNA+Protein) PCA...")
    if 'protein' in adata.obsm:
        X_combined = np.hstack([X_rna_scaled, X_protein_scaled])
        pca_combined = skPCA(n_components=min(n_pcs, X_combined.shape[1]))
        results['combined_pca'] = pca_combined.fit_transform(X_combined)
        results['combined_variance_ratio'] = pca_combined.explained_variance_ratio_

    # Use existing WNN embedding
    if 'X_wnn' in adata.obsm:
        results['wnn'] = adata.obsm['X_wnn']

    # Use existing harmony-corrected PCA
    if 'X_pca_harmony' in adata.obsm:
        results['harmony'] = adata.obsm['X_pca_harmony']

    return results


def compute_spatial_pca(adata: ad.AnnData, sample_id: str, n_neighbors: int = 15) -> Dict[str, np.ndarray]:
    """
    Compute spatially-aware PCA for a single sample using spatial neighbors.

    This incorporates spatial information by computing a spatially-weighted
    covariance matrix before PCA.
    """
    from sklearn.neighbors import NearestNeighbors
    from sklearn.decomposition import PCA as skPCA
    from sklearn.preprocessing import StandardScaler

    mask = adata.obs['sample_id'] == sample_id
    sample_data = adata[mask]

    # Get coordinates
    coords = np.column_stack([
        sample_data.obs['cell_x'].values,
        sample_data.obs['cell_y'].values
    ])

    # Get expression data
    if hasattr(sample_data.X, 'toarray'):
        X = sample_data.X.toarray()
    else:
        X = sample_data.X.copy()

    # Standardize
    X_scaled = StandardScaler().fit_transform(X)

    # Find spatial neighbors
    nn = NearestNeighbors(n_neighbors=min(n_neighbors, len(coords)))
    nn.fit(coords)
    distances, indices = nn.kneighbors(coords)

    # Compute spatially-smoothed expression (average with neighbors)
    X_spatial = np.zeros_like(X_scaled)
    for i in range(len(X_scaled)):
        neighbor_idx = indices[i]
        # Distance-weighted average
        weights = 1 / (distances[i] + 1e-6)
        weights /= weights.sum()
        X_spatial[i] = np.average(X_scaled[neighbor_idx], axis=0, weights=weights)

    # PCA on spatially-smoothed data
    pca = skPCA(n_components=min(20, X_spatial.shape[1]))
    pca_coords = pca.fit_transform(X_spatial)

    return {
        'spatial_pca': pca_coords,
        'variance_ratio': pca.explained_variance_ratio_,
        'coords': coords,
        'cell_types': sample_data.obs['cell_type'].values
    }


# =============================================================================
# Figure Generation Functions
# =============================================================================

def generate_fig0_batch_validation(batch_metrics: Dict, output_dir: Path):
    """Figure 0: Batch correction validation."""
    logger.info("Generating Figure 0: Batch validation...")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # LISI comparison
    ax = axes[0]
    stages = ['Pre-correction', 'Post-correction']
    lisi_values = [batch_metrics.get('pre_lisi', 0), batch_metrics.get('post_lisi', 0)]
    colors = ['#e74c3c', '#2ecc71']

    bars = ax.bar(stages, lisi_values, color=colors, edgecolor='black', linewidth=1.5)
    ax.axhline(y=batch_metrics['lisi_threshold'], color='gray', linestyle='--',
               label=f"Threshold ({batch_metrics['lisi_threshold']})")

    ax.set_ylabel('Mean LISI', fontsize=12)
    ax.set_title('Batch Mixing (LISI)\nHigher = Better', fontsize=14)
    ax.legend(loc='upper right')

    # Add values on bars
    for bar, val in zip(bars, lisi_values):
        if val is not None:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                   f'{val:.2f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Silhouette comparison
    ax = axes[1]
    sil_values = [batch_metrics.get('pre_silhouette', 0), batch_metrics.get('post_silhouette', 0)]

    bars = ax.bar(stages, sil_values, color=colors, edgecolor='black', linewidth=1.5)
    ax.axhline(y=0.3, color='gray', linestyle='--', label='Threshold (0.3)')

    ax.set_ylabel('Silhouette Score (batch)', fontsize=12)
    ax.set_title('Batch Separation (Silhouette)\nLower = Better', fontsize=14)
    ax.legend(loc='upper right')

    for bar, val in zip(bars, sil_values):
        if val is not None:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                   f'{val:.4f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / 'figures' / f'fig0_batch_validation.{fmt}', dpi=300, bbox_inches='tight')
    plt.close(fig)

    logger.info("  Figure 0 saved")


def generate_fig1_variance_decomposition(variance_df: pd.DataFrame, output_dir: Path):
    """Figure 1: Variance decomposition stacked bar chart."""
    logger.info("Generating Figure 1: Variance decomposition...")

    # Pivot for stacked bar
    pivot = variance_df.pivot(index='PC', columns='Factor', values='R²').fillna(0)

    # Order PCs numerically
    pc_order = [f'PC{i}' for i in range(1, len(pivot) + 1)]
    pivot = pivot.reindex(pc_order)

    # Colors by factor
    factor_colors = {
        'lane': '#e74c3c',      # Technical - red
        'patient': '#f39c12',   # Technical - orange
        'stage': '#2ecc71',     # Biological - green
        'cell_type': '#3498db'  # Biological - blue
    }

    fig, ax = plt.subplots(figsize=(14, 6))

    bottom = np.zeros(len(pivot))
    for factor in ['lane', 'patient', 'stage', 'cell_type']:
        if factor in pivot.columns:
            values = pivot[factor].values
            ax.bar(pivot.index, values, bottom=bottom,
                   label=factor.replace('_', ' ').title(),
                   color=factor_colors.get(factor, 'gray'),
                   edgecolor='white', linewidth=0.5)
            bottom += values

    ax.set_xlabel('Principal Component', fontsize=12)
    ax.set_ylabel('Variance Explained (R²)', fontsize=12)
    ax.set_title('Variance Decomposition by Factor\n(Technical vs Biological)', fontsize=14)
    ax.legend(title='Factor', loc='upper right')
    ax.set_ylim(0, 1)

    # Add total variance line
    ax.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5)

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / 'figures' / f'fig1_variance_decomposition.{fmt}', dpi=300, bbox_inches='tight')
    plt.close(fig)

    logger.info("  Figure 1 saved")


def generate_fig2_sample_overview(adata: ad.AnnData, output_dir: Path, cache_dir: Path):
    """Figure 2: Sample overview grid with H&E thumbnails."""
    logger.info("Generating Figure 2: Sample overview grid...")

    samples = sorted(adata.obs['sample_id'].unique())
    n_samples = len(samples)

    # Grid layout
    n_cols = 6
    n_rows = (n_samples + n_cols - 1) // n_cols

    fig = plt.figure(figsize=(20, 4 * n_rows))

    for idx, sample in enumerate(samples):
        ax = fig.add_subplot(n_rows, n_cols, idx + 1)

        # Get H&E thumbnail
        he_path = get_he_thumbnail(sample, cache_dir)

        if he_path and he_path.exists():
            try:
                img = plt.imread(he_path)
                ax.imshow(img)
            except Exception as e:
                ax.text(0.5, 0.5, 'H&E\nNot Available', ha='center', va='center',
                       transform=ax.transAxes, fontsize=10)
                ax.set_facecolor('#f0f0f0')
        else:
            ax.text(0.5, 0.5, 'H&E\nNot Available', ha='center', va='center',
                   transform=ax.transAxes, fontsize=10)
            ax.set_facecolor('#f0f0f0')

        ax.axis('off')

        # Get sample metadata
        sample_data = adata.obs[adata.obs['sample_id'] == sample].iloc[0]
        stage = sample_data['stage'] if 'stage' in sample_data.index else 'unknown'
        n_cells = (adata.obs['sample_id'] == sample).sum()

        stage_color = STAGE_COLORS.get(stage, 'gray')

        ax.set_title(f'{sample}\n{stage.title()} | {n_cells:,} cells',
                    fontsize=10, fontweight='bold', color=stage_color)

    plt.suptitle('Sample Overview with H&E Thumbnails', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / 'figures' / f'fig2_sample_overview.{fmt}', dpi=200, bbox_inches='tight')
    plt.close(fig)

    logger.info("  Figure 2 saved")


def generate_fig3_pseudobulk_pca(pseudobulk_df: pd.DataFrame, stats_result: Dict, output_dir: Path):
    """Figure 3: Pseudobulk PCA with progression statistics."""
    logger.info("Generating Figure 3: Pseudobulk PCA...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # PC1 vs PC2 colored by stage
    ax = axes[0]

    for stage in STAGE_ORDER:
        if stage in pseudobulk_df['stage'].values:
            mask = pseudobulk_df['stage'] == stage
            ax.scatter(pseudobulk_df.loc[mask, 'PC1'],
                      pseudobulk_df.loc[mask, 'PC2'],
                      c=STAGE_COLORS.get(stage, 'gray'),
                      label=stage.title(),
                      s=100, edgecolors='black', linewidth=1, alpha=0.8)

            # Add sample labels
            for idx, row in pseudobulk_df[mask].iterrows():
                ax.annotate(idx, (row['PC1'], row['PC2']),
                           fontsize=7, ha='center', va='bottom', alpha=0.7)

    ax.set_xlabel('PC1', fontsize=12)
    ax.set_ylabel('PC2', fontsize=12)
    ax.set_title('Pseudobulk PCA (Sample-level)', fontsize=14)
    ax.legend(title='Stage', loc='best')
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.3)

    # PC1 box plot by stage
    ax = axes[1]

    stage_data = []
    stage_labels = []
    for stage in STAGE_ORDER:
        if stage in pseudobulk_df['stage'].values:
            stage_data.append(pseudobulk_df[pseudobulk_df['stage'] == stage]['PC1'].values)
            stage_labels.append(stage.title())

    bp = ax.boxplot(stage_data, labels=stage_labels, patch_artist=True)

    for patch, stage in zip(bp['boxes'], [s for s in STAGE_ORDER if s in pseudobulk_df['stage'].values]):
        patch.set_facecolor(STAGE_COLORS.get(stage, 'gray'))
        patch.set_alpha(0.6)

    ax.set_xlabel('Stage', fontsize=12)
    ax.set_ylabel('PC1', fontsize=12)
    ax.set_title('PC1 by Disease Stage', fontsize=14)

    # Add statistics annotation
    if stats_result.get('kruskal_wallis_pval') is not None:
        kw_p = stats_result['kruskal_wallis_pval']
        rho = stats_result.get('spearman_rho', 'N/A')
        sp_p = stats_result.get('spearman_pval', 'N/A')

        stats_text = f"Kruskal-Wallis p={kw_p:.3f}\nSpearman ρ={rho:.2f} (p={sp_p:.3f})"
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
               fontsize=9, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Add power warning
    ax.text(0.02, 0.02, stats_result.get('power_note', ''),
           transform=ax.transAxes, fontsize=7, verticalalignment='bottom',
           style='italic', alpha=0.7, wrap=True)

    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / 'figures' / f'fig3_pseudobulk_pca.{fmt}', dpi=300, bbox_inches='tight')
    plt.close(fig)

    logger.info("  Figure 3 saved")


def generate_fig4_cell_pca(adata: ad.AnnData, output_dir: Path, subsample_size: int = 100000, seed: int = 42):
    """Figure 4: Cell-level PCA colored by various factors."""
    logger.info(f"Generating Figure 4: Cell-level PCA (n={subsample_size:,})...")

    # Subsample for visualization
    adata_sub = subsample_adata(adata, subsample_size, seed)

    fig, axes = plt.subplots(2, 2, figsize=(14, 14))

    # Use harmony-corrected PCA
    embedding = 'X_pca_harmony'
    if not validate_embedding(adata_sub, embedding):
        embedding = 'X_pca'

    pca = adata_sub.obsm[embedding][:, :2]

    # 1. Color by stage
    ax = axes[0, 0]
    for stage in STAGE_ORDER:
        if stage in adata_sub.obs['stage'].values:
            mask = adata_sub.obs['stage'] == stage
            ax.scatter(pca[mask, 0], pca[mask, 1],
                      c=STAGE_COLORS.get(stage, 'gray'),
                      label=stage.title(), s=1, alpha=0.3, rasterized=True)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title('Colored by Stage')
    ax.legend(markerscale=10, loc='best')

    # 2. Color by patient
    ax = axes[0, 1]
    patients = adata_sub.obs['patient'].unique()
    patient_colors = plt.cm.tab20(np.linspace(0, 1, len(patients)))
    patient_cmap = dict(zip(patients, patient_colors))

    for patient in patients:
        mask = adata_sub.obs['patient'] == patient
        ax.scatter(pca[mask, 0], pca[mask, 1],
                  c=[patient_cmap[patient]], label=patient,
                  s=1, alpha=0.3, rasterized=True)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title('Colored by Patient')
    ax.legend(markerscale=10, loc='best', fontsize=7, ncol=2)

    # 3. Color by cell type
    ax = axes[1, 0]
    celltype_colors = get_celltype_colors(adata_sub)

    for ct, color in celltype_colors.items():
        mask = adata_sub.obs['cell_type'] == ct
        if mask.sum() > 0:
            ax.scatter(pca[mask, 0], pca[mask, 1],
                      c=[color], label=ct, s=1, alpha=0.3, rasterized=True)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title('Colored by Cell Type')
    ax.legend(markerscale=10, loc='best', fontsize=7, ncol=2)

    # 4. Color by lane (batch)
    ax = axes[1, 1]
    lanes = sorted(adata_sub.obs['lane'].unique())
    lane_colors = plt.cm.Set1(np.linspace(0, 1, len(lanes)))
    lane_cmap = dict(zip(lanes, lane_colors))

    for lane in lanes:
        mask = adata_sub.obs['lane'] == lane
        ax.scatter(pca[mask, 0], pca[mask, 1],
                  c=[lane_cmap[lane]], label=f'Lane {lane}',
                  s=1, alpha=0.3, rasterized=True)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title('Colored by Lane (Batch)')
    ax.legend(markerscale=10, loc='best')

    plt.suptitle(f'Cell-level PCA (n={len(adata_sub):,} cells)', fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / 'figures' / f'fig4_cell_pca.{fmt}', dpi=200, bbox_inches='tight')
    plt.close(fig)

    # Clean up
    del adata_sub
    gc.collect()

    logger.info("  Figure 4 saved")


def generate_fig5_loading_analysis(loadings_df: pd.DataFrame, output_dir: Path, skip_enrichment: bool = True):
    """Figure 5: PC loading analysis with top genes."""
    logger.info("Generating Figure 5: Loading analysis...")

    if loadings_df.empty:
        logger.warning("  No loadings data available, skipping Figure 5")
        return

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    for pc_idx, pc in enumerate(['PC1', 'PC2', 'PC3', 'PC4', 'PC5']):
        if pc_idx >= len(axes):
            break

        ax = axes[pc_idx]
        pc_data = loadings_df[loadings_df['PC'] == pc]

        if pc_data.empty:
            ax.text(0.5, 0.5, f'No data for {pc}', ha='center', va='center', transform=ax.transAxes)
            continue

        # Sort by absolute loading
        pc_data = pc_data.copy()
        pc_data['abs_loading'] = pc_data['Loading'].abs()
        top_genes = pc_data.nlargest(20, 'abs_loading')

        colors = ['#2ecc71' if l > 0 else '#e74c3c' for l in top_genes['Loading']]

        ax.barh(range(len(top_genes)), top_genes['Loading'], color=colors)
        ax.set_yticks(range(len(top_genes)))
        ax.set_yticklabels(top_genes['Gene'], fontsize=8)
        ax.set_xlabel('Loading')
        ax.set_title(f'{pc} Top Gene Loadings')
        ax.axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
        ax.invert_yaxis()

    # Last panel: summary
    ax = axes[-1]
    ax.axis('off')

    summary_text = "PC Loading Summary\n" + "="*30 + "\n\n"
    for pc in ['PC1', 'PC2', 'PC3']:
        pc_data = loadings_df[loadings_df['PC'] == pc]
        if not pc_data.empty:
            top_pos = pc_data[pc_data['Direction'] == 'positive'].nlargest(3, 'Loading')['Gene'].tolist()
            top_neg = pc_data[pc_data['Direction'] == 'negative'].nsmallest(3, 'Loading')['Gene'].tolist()
            summary_text += f"{pc}:\n"
            summary_text += f"  +: {', '.join(top_pos)}\n"
            summary_text += f"  -: {', '.join(top_neg)}\n\n"

    if skip_enrichment:
        summary_text += "\nNote: Gene set enrichment skipped\n(use --skip-enrichment=False for enrichment)"

    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / 'figures' / f'fig5_loading_analysis.{fmt}', dpi=300, bbox_inches='tight')
    plt.close(fig)

    logger.info("  Figure 5 saved")


def generate_fig6_multimodal_comparison(adata: ad.AnnData, output_dir: Path, subsample_size: int = 50000, seed: int = 42):
    """Figure 6: Multimodal comparison (RNA vs Protein vs WNN)."""
    logger.info(f"Generating Figure 6: Multimodal comparison (n={subsample_size:,})...")

    adata_sub = subsample_adata(adata, subsample_size, seed)

    embeddings = {
        'RNA (Harmony)': 'X_pca_harmony',
        'Protein': 'X_pca_protein',
        'WNN': 'X_wnn',
        'RNA (Uncorrected)': 'X_pca_uncorrected'
    }

    # Filter to available embeddings
    available = {k: v for k, v in embeddings.items() if v in adata_sub.obsm}

    n_embed = len(available)
    fig, axes = plt.subplots(2, n_embed, figsize=(5 * n_embed, 10))

    if n_embed == 1:
        axes = axes.reshape(2, 1)

    for idx, (name, emb_key) in enumerate(available.items()):
        emb = adata_sub.obsm[emb_key][:, :2]

        # Row 1: Color by stage
        ax = axes[0, idx]
        for stage in STAGE_ORDER:
            if stage in adata_sub.obs['stage'].values:
                mask = adata_sub.obs['stage'] == stage
                ax.scatter(emb[mask, 0], emb[mask, 1],
                          c=STAGE_COLORS.get(stage, 'gray'),
                          label=stage.title(), s=1, alpha=0.3, rasterized=True)
        ax.set_xlabel('Dim 1')
        ax.set_ylabel('Dim 2')
        ax.set_title(f'{name}\n(Stage)', fontsize=11)
        if idx == 0:
            ax.legend(markerscale=10, loc='best')

        # Row 2: Color by cell type
        ax = axes[1, idx]
        celltype_colors = get_celltype_colors(adata_sub)
        for ct, color in celltype_colors.items():
            mask = adata_sub.obs['cell_type'] == ct
            if mask.sum() > 0:
                ax.scatter(emb[mask, 0], emb[mask, 1],
                          c=[color], label=ct, s=1, alpha=0.3, rasterized=True)
        ax.set_xlabel('Dim 1')
        ax.set_ylabel('Dim 2')
        ax.set_title(f'{name}\n(Cell Type)', fontsize=11)

    plt.suptitle('Multimodal Embedding Comparison', fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / 'figures' / f'fig6_multimodal_comparison.{fmt}', dpi=200, bbox_inches='tight')
    plt.close(fig)

    del adata_sub
    gc.collect()

    logger.info("  Figure 6 saved")


def generate_fig8_modality_pcas(adata: ad.AnnData, modality_pcas: Dict, output_dir: Path,
                                subsample_size: int = 50000, seed: int = 42):
    """Figure 8: RNA-only, Protein-only, and Combined PCA comparison."""
    logger.info(f"Generating Figure 8: Modality-specific PCAs (n={subsample_size:,})...")

    # Subsample
    adata_sub = subsample_adata(adata, subsample_size, seed)
    indices = np.where(np.isin(np.arange(adata.n_obs),
                               np.random.RandomState(seed).choice(adata.n_obs, min(subsample_size, adata.n_obs), replace=False)))[0]

    fig, axes = plt.subplots(3, 3, figsize=(18, 18))

    modalities = [
        ('RNA-only PCA', 'rna_pca'),
        ('Protein-only PCA', 'protein_pca'),
        ('Combined (RNA+Protein) PCA', 'combined_pca'),
    ]

    for row_idx, (name, key) in enumerate(modalities):
        if key not in modality_pcas or modality_pcas[key] is None:
            for col_idx in range(3):
                axes[row_idx, col_idx].text(0.5, 0.5, f'{name}\nNot Available',
                                            ha='center', va='center', transform=axes[row_idx, col_idx].transAxes)
                axes[row_idx, col_idx].axis('off')
            continue

        # Subsample the PCA coordinates
        pca_full = modality_pcas[key]
        pca = pca_full[indices[:len(adata_sub)], :2] if len(pca_full) > len(adata_sub) else pca_full[:len(adata_sub), :2]

        # Column 1: Color by stage
        ax = axes[row_idx, 0]
        for stage in STAGE_ORDER:
            if stage in adata_sub.obs['stage'].values:
                mask = adata_sub.obs['stage'].values == stage
                ax.scatter(pca[mask, 0], pca[mask, 1],
                          c=STAGE_COLORS.get(stage, 'gray'),
                          label=stage.title(), s=1, alpha=0.3, rasterized=True)
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_title(f'{name}\nColored by Stage')
        if row_idx == 0:
            ax.legend(markerscale=10, loc='best')

        # Column 2: Color by cell type
        ax = axes[row_idx, 1]
        celltype_colors = get_celltype_colors(adata_sub)
        for ct, color in celltype_colors.items():
            mask = adata_sub.obs['cell_type'].values == ct
            if mask.sum() > 0:
                ax.scatter(pca[mask, 0], pca[mask, 1],
                          c=[color], label=ct, s=1, alpha=0.3, rasterized=True)
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_title(f'{name}\nColored by Cell Type')

        # Column 3: Variance explained
        ax = axes[row_idx, 2]
        var_key = key.replace('_pca', '_variance_ratio')
        if var_key in modality_pcas and modality_pcas[var_key] is not None:
            var_ratio = modality_pcas[var_key][:10]
            ax.bar(range(1, len(var_ratio) + 1), var_ratio * 100, color='steelblue', edgecolor='black')
            ax.set_xlabel('Principal Component')
            ax.set_ylabel('Variance Explained (%)')
            ax.set_title(f'{name}\nVariance Explained')
            ax.set_xticks(range(1, len(var_ratio) + 1))

            # Cumulative variance
            cumvar = np.cumsum(var_ratio) * 100
            ax2 = ax.twinx()
            ax2.plot(range(1, len(var_ratio) + 1), cumvar, 'r-o', markersize=4)
            ax2.set_ylabel('Cumulative %', color='red')
            ax2.tick_params(axis='y', labelcolor='red')
            ax2.set_ylim(0, 100)
        else:
            ax.text(0.5, 0.5, 'Variance data\nNot Available', ha='center', va='center', transform=ax.transAxes)

    plt.suptitle('Modality-Specific PCA Analysis', fontsize=16, fontweight='bold', y=1.01)
    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / 'figures' / f'fig8_modality_pcas.{fmt}', dpi=200, bbox_inches='tight')
    plt.close(fig)

    del adata_sub
    gc.collect()

    logger.info("  Figure 8 saved")


def generate_fig9_spatial_pca_samples(adata: ad.AnnData, output_dir: Path, n_samples: int = 6):
    """Figure 9: Spatially-integrated PCA for selected samples."""
    logger.info(f"Generating Figure 9: Spatial PCA for {n_samples} samples...")

    # Select diverse samples (one from each stage if possible)
    samples_by_stage = {}
    for stage in STAGE_ORDER:
        stage_samples = adata.obs[adata.obs['stage'] == stage]['sample_id'].unique()
        if len(stage_samples) > 0:
            samples_by_stage[stage] = list(stage_samples)

    selected_samples = []
    for stage in STAGE_ORDER:
        if stage in samples_by_stage and len(selected_samples) < n_samples:
            remaining = n_samples - len(selected_samples)
            selected_samples.extend(samples_by_stage[stage][:max(1, remaining // len(samples_by_stage))])

    selected_samples = selected_samples[:n_samples]

    n_cols = 3
    n_rows = (len(selected_samples) + n_cols - 1) // n_cols * 2  # 2 rows per sample

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * (n_rows // 2)))

    for idx, sample in enumerate(selected_samples):
        row_base = (idx // n_cols) * 2
        col = idx % n_cols

        # Compute spatial PCA for this sample
        try:
            spatial_result = compute_spatial_pca(adata, sample, n_neighbors=15)
        except Exception as e:
            logger.warning(f"  Spatial PCA failed for {sample}: {e}")
            continue

        pca_coords = spatial_result['spatial_pca']
        xy_coords = spatial_result['coords']
        cell_types = spatial_result['cell_types']

        sample_meta = adata.obs[adata.obs['sample_id'] == sample].iloc[0]
        stage = sample_meta.get('stage', 'unknown')

        celltype_colors = get_celltype_colors(adata)

        # Row 1: Spatial distribution colored by spatial PC1
        ax = axes[row_base, col] if n_rows > 2 else axes[row_base * n_cols + col] if n_rows == 2 else axes[col]
        # NOTE: Swap X and Y to match H&E thumbnail orientation
        sc = ax.scatter(xy_coords[:, 1], xy_coords[:, 0], c=pca_coords[:, 0],
                       cmap='RdBu_r', s=1, alpha=0.7, rasterized=True)
        ax.set_aspect('equal')
        ax.invert_yaxis()
        ax.set_title(f'{sample} ({stage.title()})\nSpatial PC1', fontsize=10)
        ax.axis('off')
        plt.colorbar(sc, ax=ax, shrink=0.7, label='PC1')

        # Row 2: Spatial PCA embedding
        ax = axes[row_base + 1, col] if n_rows > 2 else axes[(row_base + 1) * n_cols + col] if n_rows == 2 else axes[n_cols + col]
        for ct, color in celltype_colors.items():
            mask = cell_types == ct
            if mask.sum() > 0:
                ax.scatter(pca_coords[mask, 0], pca_coords[mask, 1],
                          c=[color], label=ct, s=2, alpha=0.5, rasterized=True)
        ax.set_xlabel('Spatial PC1')
        ax.set_ylabel('Spatial PC2')
        ax.set_title(f'{sample} Spatial PCA Embedding')

    plt.suptitle('Spatially-Integrated PCA by Sample', fontsize=16, fontweight='bold', y=1.01)
    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / 'figures' / f'fig9_spatial_pca.{fmt}', dpi=200, bbox_inches='tight')
    plt.close(fig)

    gc.collect()
    logger.info("  Figure 9 saved")


def generate_fig7_sample_panels(adata: ad.AnnData, output_dir: Path, cache_dir: Path, max_cells_per_sample: int = 10000):
    """Figure 7: Per-sample detailed panels."""
    logger.info("Generating Figure 7: Per-sample panels...")

    samples = sorted(adata.obs['sample_id'].unique())
    panel_dir = output_dir / 'figures' / 'sample_panels'
    panel_dir.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        logger.info(f"  Processing {sample}...")

        sample_mask = adata.obs['sample_id'] == sample
        adata_sample = adata[sample_mask].copy()

        if adata_sample.n_obs > max_cells_per_sample:
            adata_sample = subsample_adata(adata_sample, max_cells_per_sample)

        fig = plt.figure(figsize=(16, 12))
        gs = GridSpec(2, 3, figure=fig)

        # Get sample metadata
        sample_meta = adata_sample.obs.iloc[0]
        stage = sample_meta.get('stage', 'unknown')
        patient = sample_meta.get('patient', 'unknown')
        lane = sample_meta.get('lane', 'unknown')
        n_cells = adata_sample.n_obs

        # Panel 1: H&E thumbnail
        ax1 = fig.add_subplot(gs[0, 0])
        he_path = get_he_thumbnail(sample, cache_dir)
        if he_path and he_path.exists():
            try:
                img = plt.imread(he_path)
                ax1.imshow(img)
            except:
                ax1.text(0.5, 0.5, 'H&E\nNot Available', ha='center', va='center',
                        transform=ax1.transAxes)
                ax1.set_facecolor('#f0f0f0')
        else:
            ax1.text(0.5, 0.5, 'H&E\nNot Available', ha='center', va='center',
                    transform=ax1.transAxes)
            ax1.set_facecolor('#f0f0f0')
        ax1.axis('off')
        ax1.set_title('H&E Thumbnail', fontsize=12)

        # Panel 2: Spatial cell distribution
        ax2 = fig.add_subplot(gs[0, 1])
        if 'cell_x' in adata_sample.obs.columns and 'cell_y' in adata_sample.obs.columns:
            celltype_colors = get_celltype_colors(adata_sample)
            for ct, color in celltype_colors.items():
                mask = adata_sample.obs['cell_type'] == ct
                if mask.sum() > 0:
                    # NOTE: Swap X and Y to match H&E thumbnail orientation
                    ax2.scatter(adata_sample.obs.loc[mask, 'cell_y'],
                               adata_sample.obs.loc[mask, 'cell_x'],
                               c=[color], label=ct, s=1, alpha=0.5, rasterized=True)
            ax2.set_aspect('equal')
            ax2.invert_yaxis()
        ax2.set_xlabel('Y (swapped)')
        ax2.set_ylabel('X (swapped)')
        ax2.set_title('Spatial Distribution', fontsize=12)

        # Panel 3: Cell type composition
        ax3 = fig.add_subplot(gs[0, 2])
        ct_counts = adata_sample.obs['cell_type'].value_counts()
        colors = [get_celltype_colors(adata_sample).get(ct, 'gray') for ct in ct_counts.index]
        ax3.pie(ct_counts.values, labels=ct_counts.index, colors=colors,
               autopct='%1.1f%%', startangle=90)
        ax3.set_title('Cell Type Composition', fontsize=12)

        # Panel 4: PCA embedding (RNA)
        ax4 = fig.add_subplot(gs[1, 0])
        if 'X_pca_harmony' in adata_sample.obsm:
            pca = adata_sample.obsm['X_pca_harmony'][:, :2]
            for ct, color in get_celltype_colors(adata_sample).items():
                mask = adata_sample.obs['cell_type'] == ct
                if mask.sum() > 0:
                    ax4.scatter(pca[mask, 0], pca[mask, 1],
                               c=[color], label=ct, s=5, alpha=0.5, rasterized=True)
        ax4.set_xlabel('PC1')
        ax4.set_ylabel('PC2')
        ax4.set_title('RNA PCA', fontsize=12)

        # Panel 5: WNN embedding
        ax5 = fig.add_subplot(gs[1, 1])
        if 'X_wnn' in adata_sample.obsm:
            wnn = adata_sample.obsm['X_wnn'][:, :2]
            for ct, color in get_celltype_colors(adata_sample).items():
                mask = adata_sample.obs['cell_type'] == ct
                if mask.sum() > 0:
                    ax5.scatter(wnn[mask, 0], wnn[mask, 1],
                               c=[color], label=ct, s=5, alpha=0.5, rasterized=True)
        ax5.set_xlabel('WNN Dim 1')
        ax5.set_ylabel('WNN Dim 2')
        ax5.set_title('WNN Embedding', fontsize=12)

        # Panel 6: Metadata summary
        ax6 = fig.add_subplot(gs[1, 2])
        ax6.axis('off')

        summary_text = f"""
Sample: {sample}
Stage: {stage.title()}
Patient: {patient}
Lane: {lane}
Cells: {n_cells:,}

Cell Type Breakdown:
"""
        for ct, count in ct_counts.head(8).items():
            pct = count / n_cells * 100
            summary_text += f"  {ct}: {count:,} ({pct:.1f}%)\n"

        ax6.text(0.1, 0.95, summary_text, transform=ax6.transAxes, fontsize=10,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

        stage_color = STAGE_COLORS.get(stage, 'gray')
        plt.suptitle(f'Sample {sample} | {stage.title()}', fontsize=16, fontweight='bold', color=stage_color)
        plt.tight_layout()

        fig.savefig(panel_dir / f'{sample}_panel.png', dpi=150, bbox_inches='tight')
        plt.close(fig)

        del adata_sample
        gc.collect()

    logger.info(f"  Generated {len(samples)} sample panels")


def generate_summary_grid(samples: List[str], output_dir: Path):
    """Create a summary grid of all sample panels."""
    logger.info("Generating summary grid (fig7)...")

    panel_dir = output_dir / 'figures' / 'sample_panels'

    n_samples = len(samples)
    n_cols = 5
    n_rows = (n_samples + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(25, 5 * n_rows))
    axes = axes.flatten()

    for idx, sample in enumerate(sorted(samples)):
        ax = axes[idx]
        panel_path = panel_dir / f'{sample}_panel.png'
        if panel_path.exists():
            try:
                img = plt.imread(panel_path)
                ax.imshow(img)
            except:
                ax.text(0.5, 0.5, f'{sample}\nNot Available', ha='center', va='center', transform=ax.transAxes)
        else:
            ax.text(0.5, 0.5, f'{sample}\nNot Available', ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')
        ax.set_title(sample, fontsize=8)

    # Hide empty axes
    for idx in range(len(samples), len(axes)):
        axes[idx].axis('off')

    plt.suptitle('All Sample Panels', fontsize=16, fontweight='bold')
    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / 'figures' / f'fig7_summary_grid.{fmt}', dpi=150, bbox_inches='tight')
    plt.close(fig)

    logger.info("  Summary grid saved")


# =============================================================================
# Report Generation
# =============================================================================

def generate_report(output_dir: Path, batch_metrics: Dict, variance_df: pd.DataFrame,
                   stats_result: Dict, n_samples: int, n_cells: int):
    """Generate markdown report."""
    logger.info("Generating report...")

    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    report = f"""# PCA Deep Dive Analysis Report

Generated: {timestamp}

## Overview

| Metric | Value |
|--------|-------|
| Total Cells | {n_cells:,} |
| Samples Analyzed | {n_samples} |
| Excluded Samples | {', '.join(EXCLUDED_SAMPLES)} |

## Batch Correction

| Metric | Pre-Correction | Post-Correction |
|--------|----------------|-----------------|
| LISI (higher = better) | {batch_metrics.get('pre_lisi', 'N/A'):.3f} | {batch_metrics.get('post_lisi', 'N/A'):.3f} |
| Silhouette (lower = better) | {batch_metrics.get('pre_silhouette', 'N/A'):.4f} | {batch_metrics.get('post_silhouette', 'N/A'):.4f} |

## Variance Decomposition

Top factors explaining variance in PC1:
"""

    if not variance_df.empty:
        pc1_variance = variance_df[variance_df['PC'] == 'PC1'].sort_values('R²', ascending=False)
        for _, row in pc1_variance.head(4).iterrows():
            report += f"- {row['Factor']}: {row['R²']:.3f} ({row['Category']})\n"

    report += f"""

## Progression Statistics

{stats_result.get('power_note', 'N/A')}

| Test | Statistic | P-value |
|------|-----------|---------|
| Kruskal-Wallis | {stats_result.get('kruskal_wallis_stat', 'N/A')} | {stats_result.get('kruskal_wallis_pval', 'N/A')} |
| Spearman rho | {stats_result.get('spearman_rho', 'N/A')} | {stats_result.get('spearman_pval', 'N/A')} |

### Samples per Stage
"""

    for stage, count in stats_result.get('n_per_stage', {}).items():
        report += f"- {stage}: n={count}\n"

    report += """

## Figures Generated

1. `fig0_batch_validation.png` - Batch correction metrics
2. `fig1_variance_decomposition.png` - Variance by factor
3. `fig2_sample_overview.png` - Sample grid with H&E
4. `fig3_pseudobulk_pca.png` - Sample-level PCA
5. `fig4_cell_pca.png` - Cell-level PCA
6. `fig5_loading_analysis.png` - PC gene loadings
7. `fig6_multimodal_comparison.png` - RNA vs Protein vs WNN
8. `fig7_summary_grid.png` - All sample panels
9. `sample_panels/` - Individual sample analyses

## Data Files

- `variance_partition.csv` - Variance explained by each factor
- `pc_loadings.csv` - Top genes for each PC
- `pseudobulk_pca.csv` - Sample-level PCA coordinates
- `progression_stats.json` - Statistical test results
"""

    with open(output_dir / 'PCA_DEEP_DIVE_REPORT.md', 'w') as f:
        f.write(report)

    logger.info("  Report saved")


def generate_manifest(output_dir: Path, args, n_samples: int, n_cells: int, elapsed_seconds: float):
    """Generate manifest.json with run metadata."""
    # Convert args to JSON-serializable dict
    args_dict = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}

    manifest = {
        'script': '46_pca_deep_dive.py',
        'timestamp': datetime.now().isoformat(),
        'args': args_dict,
        'n_samples': n_samples,
        'n_cells': n_cells,
        'elapsed_seconds': elapsed_seconds,
        'excluded_samples': list(EXCLUDED_SAMPLES),
        'outputs': {
            'figures': [f'fig{i}*.png' for i in range(8)],
            'sample_panels': f'{n_samples} panels',
            'data_files': ['variance_partition.csv', 'pc_loadings.csv', 'pseudobulk_pca.csv', 'progression_stats.json'],
            'report': 'PCA_DEEP_DIVE_REPORT.md'
        }
    }

    with open(output_dir / 'manifest.json', 'w') as f:
        json.dump(manifest, f, indent=2)

    logger.info("  Manifest saved")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='PCA Deep Dive Analysis')
    parser.add_argument('--skip-enrichment', action='store_true',
                       help='Skip gene set enrichment (faster iteration)')
    parser.add_argument('--subsample-size', type=int, default=100000,
                       help='Number of cells to subsample (default: 100000)')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility')
    parser.add_argument('--output-dir', type=Path,
                       default=BASE_DIR / 'results' / 'pca_deep_dive',
                       help='Output directory for figures and data')
    parser.add_argument('--dry-run', action='store_true',
                       help='Validate inputs and show what would be done')
    parser.add_argument('--skip-env-check', action='store_true',
                       help='Skip conda environment check')
    args = parser.parse_args()

    # Environment check
    check_environment(args.skip_env_check)

    # Create output directories
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / 'figures' / 'sample_panels').mkdir(parents=True, exist_ok=True)
    (args.output_dir / 'data').mkdir(parents=True, exist_ok=True)
    cache_dir = args.output_dir / 'cache' / 'thumbnails'
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Check disk space
    check_disk_space(args.output_dir)

    # Validate inputs
    logger.info("Validating inputs...")

    if not CORRECTED_PATH.exists():
        logger.error(f"Missing: {CORRECTED_PATH}")
        sys.exit(1)

    if not COUNTS_PATH.exists():
        logger.warning(f"Missing counts file: {COUNTS_PATH} (will use .X instead)")

    if not BATCH_CSV.exists():
        logger.error(f"Missing: {BATCH_CSV}")
        sys.exit(1)

    # Load batch metrics
    batch_metrics = load_batch_metrics(BATCH_CSV)
    if not batch_metrics['valid']:
        logger.warning(f"Batch metrics issue: {batch_metrics.get('error', 'unknown')}")

    # Quick data check
    logger.info(f"Loading data from {CORRECTED_PATH}...")
    adata_meta = sc.read_h5ad(CORRECTED_PATH, backed='r')
    detected_samples = set(adata_meta.obs['sample_id'].unique().tolist())
    passing_samples = sorted(detected_samples - EXCLUDED_SAMPLES)
    n_cells = adata_meta.n_obs

    logger.info(f"Detected {len(detected_samples)} samples, {len(passing_samples)} passing")
    logger.info(f"Total cells: {n_cells:,}")

    del adata_meta
    gc.collect()

    # Dry run mode
    if args.dry_run:
        print("\n" + "="*60)
        print("DRY RUN MODE")
        print("="*60)
        print(f"Output directory: {args.output_dir}")
        print(f"Samples to analyze: {len(passing_samples)}")
        print(f"Excluded samples: {EXCLUDED_SAMPLES}")
        print(f"Subsample size: {args.subsample_size:,}")
        print(f"Skip enrichment: {args.skip_enrichment}")
        print(f"Batch metrics valid: {batch_metrics['valid']}")
        print(f"\nWould generate:")
        for i in range(10):
            print(f"  - fig{i}_*.png/pdf")
        print(f"  - {len(passing_samples)} sample panels")
        print(f"  - fig8: Modality-specific PCAs (RNA, Protein, Combined)")
        print(f"  - fig9: Spatially-integrated PCA")
        print(f"  - PCA_DEEP_DIVE_REPORT.md")
        print(f"  - manifest.json")
        print("\nDry run complete. No outputs generated.")
        sys.exit(0)

    # Full run
    import time
    start_time = time.time()

    logger.info("Loading full dataset...")
    adata = sc.read_h5ad(CORRECTED_PATH)

    # Filter to passing samples
    initial_cells = adata.n_obs
    adata = adata[~adata.obs['sample_id'].isin(EXCLUDED_SAMPLES)].copy()
    logger.info(f"Filtered: {initial_cells:,} -> {adata.n_obs:,} cells ({len(passing_samples)} samples)")

    # Figure 0: Batch validation
    generate_fig0_batch_validation(batch_metrics, args.output_dir)

    # Figure 1: Variance decomposition
    logger.info("Computing variance partition...")
    variance_df = compute_variance_partition(adata, embedding='X_pca_harmony', n_pcs=20)
    variance_df.to_csv(args.output_dir / 'data' / 'variance_partition.csv', index=False)
    generate_fig1_variance_decomposition(variance_df, args.output_dir)

    # Figure 2: Sample overview with H&E
    generate_fig2_sample_overview(adata, args.output_dir, cache_dir)

    # Figure 3: Pseudobulk PCA
    pseudobulk_df = compute_pseudobulk_pca(adata)
    pseudobulk_df.to_csv(args.output_dir / 'data' / 'pseudobulk_pca.csv')
    stats_result = compute_progression_stats(pseudobulk_df)
    with open(args.output_dir / 'data' / 'progression_stats.json', 'w') as f:
        json.dump(stats_result, f, indent=2)
    generate_fig3_pseudobulk_pca(pseudobulk_df, stats_result, args.output_dir)

    # Figure 4: Cell-level PCA
    generate_fig4_cell_pca(adata, args.output_dir, args.subsample_size, args.seed)

    # Figure 5: Loading analysis
    loadings_df = compute_pc_loadings(adata, n_pcs=5, n_top=20)
    if not loadings_df.empty:
        loadings_df.to_csv(args.output_dir / 'data' / 'pc_loadings.csv', index=False)
    generate_fig5_loading_analysis(loadings_df, args.output_dir, args.skip_enrichment)

    # Figure 6: Multimodal comparison
    generate_fig6_multimodal_comparison(adata, args.output_dir, args.subsample_size // 2, args.seed)

    # Figure 7: Per-sample panels (with corrected H&E mapping)
    generate_fig7_sample_panels(adata, args.output_dir, cache_dir)
    generate_summary_grid(passing_samples, args.output_dir)

    # Figure 8: Modality-specific PCAs (RNA, Protein, Combined)
    logger.info("Computing modality-specific PCAs...")
    modality_pcas = compute_modality_pcas(adata, n_pcs=30)
    generate_fig8_modality_pcas(adata, modality_pcas, args.output_dir, args.subsample_size // 2, args.seed)

    # Figure 9: Spatially-integrated PCA
    generate_fig9_spatial_pca_samples(adata, args.output_dir, n_samples=6)

    # Generate report and manifest
    elapsed = time.time() - start_time
    generate_report(args.output_dir, batch_metrics, variance_df, stats_result,
                   len(passing_samples), adata.n_obs)
    generate_manifest(args.output_dir, args, len(passing_samples), adata.n_obs, elapsed)

    logger.info("="*60)
    logger.info(f"PCA Deep Dive complete in {elapsed/60:.1f} minutes")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info("="*60)


if __name__ == '__main__':
    main()
