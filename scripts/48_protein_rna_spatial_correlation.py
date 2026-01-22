#!/usr/bin/env python3
"""
G4X Protein-RNA Spatial Correlation Analysis
=============================================

Side-by-side spatial visualization of protein vs RNA with consistent scaling.
Computes Pearson correlation, Spearman correlation, and SSIM for each protein-RNA pair.

Key Design Decisions:
1. Protein-RNA Mapping (15 valid pairs):
   - Exact matches: CD4, CD68, FOXP3
   - Gene name variants: CD3→CD3D, CD8→CD8A, HLA-DR→HLA-DRA, etc.
   - Proxy mapping: PanCK→EPCAM (KRT8/18/19 missing from RNA panel)

2. Scaling & Normalization:
   - Z-score normalize both modalities for fair SSIM comparison
   - Consistent bounding box across all samples

3. Zero-Variance Guards:
   - Skip correlation for proteins/genes with no expression variance
   - Document all NaN cases

Usage:
    conda activate enact
    cd ~/g4x-choi-batch2-analysis
    python scripts/48_protein_rna_spatial_correlation.py 2>&1 | tee logs/48_protein_rna.log

Output:
    results/protein_rna_correlation/
    ├── sample_panels/           # 29 samples × 15 proteins = 435 PNG files
    ├── protein_summary/         # 15 summary panels
    ├── correlation_details.csv  # All metrics (435 rows)
    ├── correlation_matrix.csv   # Samples × Proteins Pearson r
    ├── ssim_matrix.csv          # Samples × Proteins SSIM
    ├── zero_variance_cases.csv  # Document any NaN correlations
    ├── summary_heatmap.png
    └── PROTEIN_RNA_REPORT.md
"""

import os
import sys
import warnings
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, List, Tuple, Any

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
        sys.exit(1)


# Check environment before heavy imports
check_environment()

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr
from scipy.ndimage import gaussian_filter
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
INPUT_PATH = BASE_DIR / 'results' / 'qc_all_samples' / 'merged' / 'merged_corrected.h5ad'
OUTPUT_DIR = BASE_DIR / 'results' / 'protein_rna_correlation'

# Corrected protein to RNA gene mapping (15 valid pairs)
# Based on RNA panel availability check
PROTEIN_TO_GENE = {
    'CD4': 'CD4',           # Exact match
    'CD68': 'CD68',         # Exact match
    'FOXP3': 'FOXP3',       # Exact match
    'CD3': 'CD3D',          # Use CD3D (CD3E also available but D is primary)
    'CD8': 'CD8A',          # Primary chain (CD8B not in panel)
    'HLA-DR': 'HLA-DRA',    # MHC-II alpha chain
    'KI67': 'MKI67',        # Proliferation marker
    'PD1': 'PDCD1',         # Checkpoint
    'PDL1': 'CD274',        # Checkpoint ligand
    'aSMA': 'ACTA2',        # CAF marker
    'CD31': 'PECAM1',       # Endothelial
    'PanCK': 'EPCAM',       # ⚠️ PROXY - KRT8/18/19 missing, use EPCAM
    'CD20': 'MS4A1',        # B cells
    'CD45': 'PTPRC',        # Pan-immune
    'CD11c': 'ITGAX',       # DCs
}

# Proteins to exclude (no RNA equivalent or control)
EXCLUDED_PROTEINS = ['ATPase', 'Isotype', 'cytoplasmicstain', 'nuclearstain']

# Stage colors for visualization
STAGE_COLORS = {
    'normal': '#2ecc71',
    'control': '#3498db',
    'metaplasia': '#f39c12',
    'cancer': '#e74c3c'
}

# Grid parameters for spatial binning
GRID_SIZE = 50  # 50x50 grid for correlation computation


# =============================================================================
# Correlation Functions
# =============================================================================

def compute_global_bounds(adata: ad.AnnData) -> Tuple[float, float, float, float]:
    """Compute global bounding box across all samples for consistent grids."""
    x_min = adata.obs['cell_x'].min()
    x_max = adata.obs['cell_x'].max()
    y_min = adata.obs['cell_y'].min()
    y_max = adata.obs['cell_y'].max()
    return (x_min, x_max, y_min, y_max)


def safe_zscore(vals: np.ndarray) -> np.ndarray:
    """Z-score normalize with zero-variance protection."""
    std = np.std(vals)
    if std < 1e-9:
        return np.zeros_like(vals)
    return (vals - np.mean(vals)) / std


def bin_to_grid(x: np.ndarray, y: np.ndarray, vals: np.ndarray,
                bounds: Tuple[float, float, float, float],
                grid_size: int = 50) -> np.ndarray:
    """
    Bin spatial values to a regular grid.

    Parameters
    ----------
    x, y : array
        Cell coordinates
    vals : array
        Values to bin (expression, intensity)
    bounds : tuple
        (x_min, x_max, y_min, y_max) for consistent grid
    grid_size : int
        Number of bins per dimension

    Returns
    -------
    grid : array
        2D binned values (grid_size x grid_size)
    """
    from scipy.stats import binned_statistic_2d

    x_edges = np.linspace(bounds[0], bounds[1], grid_size + 1)
    y_edges = np.linspace(bounds[2], bounds[3], grid_size + 1)

    grid, _, _, _ = binned_statistic_2d(
        x, y, vals,
        statistic='mean',
        bins=[x_edges, y_edges]
    )

    # Replace NaN with 0 (empty bins)
    grid = np.nan_to_num(grid, nan=0.0)

    return grid


def compute_ssim(img1: np.ndarray, img2: np.ndarray,
                 window_size: int = 7, k1: float = 0.01, k2: float = 0.03) -> float:
    """
    Compute Structural Similarity Index (SSIM) between two grids.

    Parameters
    ----------
    img1, img2 : array
        2D arrays to compare
    window_size : int
        Size of Gaussian window
    k1, k2 : float
        Stability constants

    Returns
    -------
    ssim : float
        SSIM value in range [-1, 1]
    """
    # Normalize to [0, 1]
    img1_norm = (img1 - img1.min()) / (img1.max() - img1.min() + 1e-9)
    img2_norm = (img2 - img2.min()) / (img2.max() - img2.min() + 1e-9)

    # Constants
    C1 = k1 ** 2
    C2 = k2 ** 2

    # Gaussian window
    sigma = window_size / 6.0

    # Local means
    mu1 = gaussian_filter(img1_norm, sigma)
    mu2 = gaussian_filter(img2_norm, sigma)

    # Local variances and covariance
    mu1_sq = mu1 ** 2
    mu2_sq = mu2 ** 2
    mu1_mu2 = mu1 * mu2

    sigma1_sq = gaussian_filter(img1_norm ** 2, sigma) - mu1_sq
    sigma2_sq = gaussian_filter(img2_norm ** 2, sigma) - mu2_sq
    sigma12 = gaussian_filter(img1_norm * img2_norm, sigma) - mu1_mu2

    # SSIM
    ssim_map = ((2 * mu1_mu2 + C1) * (2 * sigma12 + C2)) / \
               ((mu1_sq + mu2_sq + C1) * (sigma1_sq + sigma2_sq + C2))

    return float(np.mean(ssim_map))


def safe_correlation(protein_vals: np.ndarray, rna_vals: np.ndarray) -> Dict[str, Any]:
    """
    Compute correlation with zero-variance guards.

    Parameters
    ----------
    protein_vals, rna_vals : array
        Values to correlate

    Returns
    -------
    dict with pearson_r, spearman_r, and notes
    """
    result = {
        'pearson_r': np.nan,
        'spearman_r': np.nan,
        'pearson_p': np.nan,
        'spearman_p': np.nan,
        'note': None,
    }

    # Check for zero variance
    if np.std(protein_vals) < 1e-9:
        result['note'] = 'zero_variance_protein'
        return result

    if np.std(rna_vals) < 1e-9:
        result['note'] = 'zero_variance_rna'
        return result

    # Remove NaN pairs
    mask = ~(np.isnan(protein_vals) | np.isnan(rna_vals))
    if mask.sum() < 10:
        result['note'] = 'insufficient_data'
        return result

    protein_clean = protein_vals[mask]
    rna_clean = rna_vals[mask]

    try:
        r_pearson, p_pearson = pearsonr(protein_clean, rna_clean)
        result['pearson_r'] = r_pearson
        result['pearson_p'] = p_pearson
    except Exception:
        pass

    try:
        r_spearman, p_spearman = spearmanr(protein_clean, rna_clean)
        result['spearman_r'] = r_spearman
        result['spearman_p'] = p_spearman
    except Exception:
        pass

    return result


def analyze_sample_protein(adata: ad.AnnData, sample_id: str, protein: str, gene: str,
                          global_bounds: Tuple[float, float, float, float],
                          grid_size: int = 50) -> Dict[str, Any]:
    """
    Analyze protein-RNA correlation for one sample and one protein.

    Returns correlation metrics and grid data for visualization.
    """
    # Get sample data
    sample_mask = adata.obs['sample_id'] == sample_id
    sample_data = adata[sample_mask]

    if sample_data.n_obs < 100:
        return {
            'sample_id': sample_id,
            'protein': protein,
            'gene': gene,
            'n_cells': sample_data.n_obs,
            'error': 'insufficient_cells',
        }

    # Get coordinates
    x = sample_data.obs['cell_x'].values
    y = sample_data.obs['cell_y'].values

    # Get protein values
    protein_col = f'{protein}_intensity_mean'
    if protein_col not in sample_data.obs.columns:
        return {
            'sample_id': sample_id,
            'protein': protein,
            'gene': gene,
            'error': 'protein_column_missing',
        }
    protein_vals = sample_data.obs[protein_col].values

    # Get RNA values
    if gene not in sample_data.var_names:
        return {
            'sample_id': sample_id,
            'protein': protein,
            'gene': gene,
            'error': 'gene_not_found',
        }
    gene_idx = list(sample_data.var_names).index(gene)

    if hasattr(sample_data.X, 'toarray'):
        rna_vals = sample_data.X[:, gene_idx].toarray().flatten()
    else:
        rna_vals = sample_data.X[:, gene_idx].flatten()

    # Cell-level correlation
    cell_corr = safe_correlation(protein_vals, rna_vals)

    # Grid-based correlation and SSIM
    # First, compute sample-specific bounds for this sample
    sample_bounds = (x.min(), x.max(), y.min(), y.max())

    # Bin to grid using sample bounds (for this sample's visualization)
    protein_grid = bin_to_grid(x, y, protein_vals, sample_bounds, grid_size)
    rna_grid = bin_to_grid(x, y, rna_vals, sample_bounds, grid_size)

    # Z-score normalize for SSIM
    protein_grid_z = safe_zscore(protein_grid)
    rna_grid_z = safe_zscore(rna_grid)

    # Compute SSIM
    try:
        ssim = compute_ssim(protein_grid_z, rna_grid_z)
    except Exception:
        ssim = np.nan

    # Grid-level correlation (on flattened grids)
    grid_corr = safe_correlation(protein_grid.flatten(), rna_grid.flatten())

    return {
        'sample_id': sample_id,
        'protein': protein,
        'gene': gene,
        'n_cells': sample_data.n_obs,
        'protein_mean': np.mean(protein_vals),
        'protein_std': np.std(protein_vals),
        'rna_mean': np.mean(rna_vals),
        'rna_std': np.std(rna_vals),
        # Cell-level correlations
        'cell_pearson_r': cell_corr['pearson_r'],
        'cell_spearman_r': cell_corr['spearman_r'],
        'cell_pearson_p': cell_corr['pearson_p'],
        'cell_spearman_p': cell_corr['spearman_p'],
        # Grid-level correlations
        'grid_pearson_r': grid_corr['pearson_r'],
        'grid_spearman_r': grid_corr['spearman_r'],
        'ssim': ssim,
        'note': cell_corr.get('note'),
        # Grids for plotting (optional, memory-heavy)
        '_protein_grid': protein_grid,
        '_rna_grid': rna_grid,
        '_coords': (x, y),
        '_bounds': sample_bounds,
    }


# =============================================================================
# Visualization Functions
# =============================================================================

def generate_sample_panel(result: Dict, output_dir: Path, stage: str):
    """Generate side-by-side protein-RNA spatial comparison for one sample."""

    sample_id = result['sample_id']
    protein = result['protein']
    gene = result['gene']

    if 'error' in result:
        logger.warning(f"  Skipping {sample_id} {protein}: {result['error']}")
        return

    protein_grid = result['_protein_grid']
    rna_grid = result['_rna_grid']
    x, y = result['_coords']
    bounds = result['_bounds']

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # Panel 1: Protein spatial
    ax = axes[0]
    vmin_p, vmax_p = np.percentile(protein_grid, [2, 98])
    im = ax.imshow(protein_grid.T, origin='lower', cmap='viridis',
                   aspect='auto', vmin=vmin_p, vmax=vmax_p)
    ax.set_title(f'Protein: {protein}', fontsize=11)
    ax.set_xlabel('X bins')
    ax.set_ylabel('Y bins')
    plt.colorbar(im, ax=ax, shrink=0.8, label='Mean Intensity')

    # Panel 2: RNA spatial
    ax = axes[1]
    vmin_r, vmax_r = np.percentile(rna_grid, [2, 98])
    im = ax.imshow(rna_grid.T, origin='lower', cmap='magma',
                   aspect='auto', vmin=vmin_r, vmax=vmax_r)

    # Label with proxy note if PanCK
    gene_label = gene
    if protein == 'PanCK':
        gene_label = f'{gene} (proxy)'

    ax.set_title(f'RNA: {gene_label}', fontsize=11)
    ax.set_xlabel('X bins')
    ax.set_ylabel('Y bins')
    plt.colorbar(im, ax=ax, shrink=0.8, label='Mean Expression')

    # Panel 3: Correlation scatter
    ax = axes[2]

    # Subsample for scatter plot if too many points
    n_bins = protein_grid.size
    if n_bins > 1000:
        idx = np.random.choice(n_bins, 1000, replace=False)
        p_flat = protein_grid.flatten()[idx]
        r_flat = rna_grid.flatten()[idx]
    else:
        p_flat = protein_grid.flatten()
        r_flat = rna_grid.flatten()

    ax.scatter(p_flat, r_flat, alpha=0.4, s=10, c='steelblue')

    # Add correlation line
    if not np.isnan(result['grid_pearson_r']):
        z = np.polyfit(p_flat, r_flat, 1)
        p = np.poly1d(z)
        x_line = np.linspace(p_flat.min(), p_flat.max(), 100)
        ax.plot(x_line, p(x_line), 'r-', linewidth=2, alpha=0.7)

    ax.set_xlabel(f'{protein} (binned)')
    ax.set_ylabel(f'{gene} (binned)')
    ax.set_title('Spatial Correlation')

    # Add metrics annotation
    metrics_text = (
        f"Cell Pearson r: {result['cell_pearson_r']:.3f}\n"
        f"Cell Spearman r: {result['cell_spearman_r']:.3f}\n"
        f"Grid Pearson r: {result['grid_pearson_r']:.3f}\n"
        f"SSIM: {result['ssim']:.3f}\n"
        f"n={result['n_cells']:,} cells"
    )
    ax.text(0.02, 0.98, metrics_text, transform=ax.transAxes,
           fontsize=9, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    # Stage color indicator
    stage_color = STAGE_COLORS.get(stage, 'gray')

    plt.suptitle(f'{sample_id} | {stage.title()} | {protein} vs {gene}',
                fontsize=12, fontweight='bold', color=stage_color)
    plt.tight_layout()

    # Save
    panel_dir = output_dir / 'sample_panels'
    panel_dir.mkdir(exist_ok=True)

    filename = f'{sample_id}_{protein}_{gene}.png'
    fig.savefig(panel_dir / filename, dpi=150, bbox_inches='tight')
    plt.close(fig)


def generate_protein_summary(all_results: List[Dict], protein: str, gene: str, output_dir: Path):
    """Generate summary figure for one protein across all samples."""
    logger.info(f"  Generating summary for {protein}...")

    # Filter to this protein
    protein_results = [r for r in all_results if r['protein'] == protein and 'error' not in r]

    if not protein_results:
        logger.warning(f"    No valid results for {protein}")
        return

    # Sort by stage
    stage_order = {'normal': 0, 'control': 1, 'metaplasia': 2, 'cancer': 3}

    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, 3, figure=fig, height_ratios=[1, 1])

    # Panel 1: Cell-level Pearson by stage
    ax = fig.add_subplot(gs[0, 0])

    # Get unique stages present
    stages = sorted(set(r.get('stage', 'unknown') for r in protein_results),
                   key=lambda x: stage_order.get(x, 99))

    stage_correlations = {s: [] for s in stages}
    for r in protein_results:
        stage = r.get('stage', 'unknown')
        if not np.isnan(r['cell_pearson_r']):
            stage_correlations[stage].append(r['cell_pearson_r'])

    positions = list(range(len(stages)))
    bp_data = [stage_correlations.get(s, []) for s in stages]

    if any(bp_data):
        bp = ax.boxplot(bp_data, positions=positions, patch_artist=True)
        for patch, stage in zip(bp['boxes'], stages):
            patch.set_facecolor(STAGE_COLORS.get(stage, 'gray'))
            patch.set_alpha(0.6)

    ax.set_xticks(positions)
    ax.set_xticklabels([s.title() for s in stages])
    ax.set_ylabel('Pearson r')
    ax.set_title('Cell-level Correlation by Stage')
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)

    # Panel 2: SSIM by stage
    ax = fig.add_subplot(gs[0, 1])

    ssim_by_stage = {s: [] for s in stages}
    for r in protein_results:
        stage = r.get('stage', 'unknown')
        if not np.isnan(r['ssim']):
            ssim_by_stage[stage].append(r['ssim'])

    bp_data = [ssim_by_stage.get(s, []) for s in stages]
    if any(bp_data):
        bp = ax.boxplot(bp_data, positions=positions, patch_artist=True)
        for patch, stage in zip(bp['boxes'], stages):
            patch.set_facecolor(STAGE_COLORS.get(stage, 'gray'))
            patch.set_alpha(0.6)

    ax.set_xticks(positions)
    ax.set_xticklabels([s.title() for s in stages])
    ax.set_ylabel('SSIM')
    ax.set_title('Spatial Similarity by Stage')
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)

    # Panel 3: Sample scatter (Pearson vs SSIM)
    ax = fig.add_subplot(gs[0, 2])

    for r in protein_results:
        if not np.isnan(r['cell_pearson_r']) and not np.isnan(r['ssim']):
            stage = r.get('stage', 'unknown')
            ax.scatter(r['cell_pearson_r'], r['ssim'],
                      c=STAGE_COLORS.get(stage, 'gray'),
                      s=60, alpha=0.7, edgecolors='black', linewidth=0.5)

    ax.set_xlabel('Pearson r')
    ax.set_ylabel('SSIM')
    ax.set_title('Pearson vs SSIM by Sample')

    # Add legend
    for stage in stages:
        ax.scatter([], [], c=STAGE_COLORS.get(stage, 'gray'), label=stage.title())
    ax.legend(loc='best', fontsize=8)

    # Panel 4: Mean expression by stage (protein)
    ax = fig.add_subplot(gs[1, 0])

    protein_means = {s: [] for s in stages}
    for r in protein_results:
        stage = r.get('stage', 'unknown')
        protein_means[stage].append(r['protein_mean'])

    bp_data = [protein_means.get(s, []) for s in stages]
    if any(bp_data):
        bp = ax.boxplot(bp_data, positions=positions, patch_artist=True)
        for patch, stage in zip(bp['boxes'], stages):
            patch.set_facecolor(STAGE_COLORS.get(stage, 'gray'))
            patch.set_alpha(0.6)

    ax.set_xticks(positions)
    ax.set_xticklabels([s.title() for s in stages])
    ax.set_ylabel('Mean Intensity')
    ax.set_title(f'{protein} Protein Expression')

    # Panel 5: Mean expression by stage (RNA)
    ax = fig.add_subplot(gs[1, 1])

    rna_means = {s: [] for s in stages}
    for r in protein_results:
        stage = r.get('stage', 'unknown')
        rna_means[stage].append(r['rna_mean'])

    bp_data = [rna_means.get(s, []) for s in stages]
    if any(bp_data):
        bp = ax.boxplot(bp_data, positions=positions, patch_artist=True)
        for patch, stage in zip(bp['boxes'], stages):
            patch.set_facecolor(STAGE_COLORS.get(stage, 'gray'))
            patch.set_alpha(0.6)

    ax.set_xticks(positions)
    ax.set_xticklabels([s.title() for s in stages])
    ax.set_ylabel('Mean Expression')

    gene_label = gene
    if protein == 'PanCK':
        gene_label = f'{gene} (proxy)'
    ax.set_title(f'{gene_label} RNA Expression')

    # Panel 6: Summary statistics
    ax = fig.add_subplot(gs[1, 2])
    ax.axis('off')

    # Compute overall statistics
    pearson_vals = [r['cell_pearson_r'] for r in protein_results if not np.isnan(r['cell_pearson_r'])]
    ssim_vals = [r['ssim'] for r in protein_results if not np.isnan(r['ssim'])]

    summary_text = f"""
{protein} → {gene} Summary
{'='*35}

Samples analyzed: {len(protein_results)}
Zero-variance cases: {len([r for r in all_results if r['protein'] == protein and r.get('note') == 'zero_variance_protein' or r.get('note') == 'zero_variance_rna'])}

Cell Pearson r:
  Mean: {np.mean(pearson_vals):.3f} (±{np.std(pearson_vals):.3f})
  Range: [{np.min(pearson_vals):.3f}, {np.max(pearson_vals):.3f}]

SSIM:
  Mean: {np.mean(ssim_vals):.3f} (±{np.std(ssim_vals):.3f})
  Range: [{np.min(ssim_vals):.3f}, {np.max(ssim_vals):.3f}]
"""

    if protein == 'PanCK':
        summary_text += """
⚠️ NOTE: PanCK → EPCAM is a PROXY mapping.
Cytokeratins (KRT8/18/19) are not in the RNA panel.
EPCAM is used as an epithelial marker alternative.
"""

    ax.text(0.1, 0.95, summary_text, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    plt.suptitle(f'{protein} vs {gene} - Summary Across All Samples',
                fontsize=14, fontweight='bold')
    plt.tight_layout()

    # Save
    summary_dir = output_dir / 'protein_summary'
    summary_dir.mkdir(exist_ok=True)

    filename = f'{protein}_{gene}_summary.png'
    fig.savefig(summary_dir / filename, dpi=200, bbox_inches='tight')
    fig.savefig(summary_dir / filename.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close(fig)


def generate_summary_heatmap(all_results: List[Dict], output_dir: Path):
    """Generate summary heatmap of correlations across all samples and proteins."""
    logger.info("Generating summary heatmap...")

    # Create DataFrames
    pearson_data = {}
    ssim_data = {}

    for r in all_results:
        if 'error' in r:
            continue

        sample_id = r['sample_id']
        protein = r['protein']

        if sample_id not in pearson_data:
            pearson_data[sample_id] = {}
            ssim_data[sample_id] = {}

        pearson_data[sample_id][protein] = r['cell_pearson_r']
        ssim_data[sample_id][protein] = r['ssim']

    pearson_df = pd.DataFrame(pearson_data).T
    ssim_df = pd.DataFrame(ssim_data).T

    # Sort columns by protein order
    protein_order = [p for p in PROTEIN_TO_GENE.keys() if p in pearson_df.columns]
    pearson_df = pearson_df[protein_order]
    ssim_df = ssim_df[protein_order]

    fig, axes = plt.subplots(1, 2, figsize=(16, max(8, len(pearson_df) * 0.4)))

    # Panel 1: Pearson heatmap
    ax = axes[0]
    sns.heatmap(pearson_df, cmap='RdBu_r', center=0, ax=ax,
                vmin=-1, vmax=1, annot=True, fmt='.2f', annot_kws={'fontsize': 7},
                cbar_kws={'label': 'Pearson r'})
    ax.set_title('Cell-level Pearson Correlation', fontsize=12)
    ax.set_xlabel('Protein')
    ax.set_ylabel('Sample')

    # Panel 2: SSIM heatmap
    ax = axes[1]
    sns.heatmap(ssim_df, cmap='YlGnBu', ax=ax,
                vmin=0, vmax=1, annot=True, fmt='.2f', annot_kws={'fontsize': 7},
                cbar_kws={'label': 'SSIM'})
    ax.set_title('Spatial Similarity (SSIM)', fontsize=12)
    ax.set_xlabel('Protein')
    ax.set_ylabel('Sample')

    plt.suptitle('Protein-RNA Correlation Summary', fontsize=14, fontweight='bold')
    plt.tight_layout()

    fig.savefig(output_dir / 'summary_heatmap.png', dpi=200, bbox_inches='tight')
    fig.savefig(output_dir / 'summary_heatmap.pdf', bbox_inches='tight')
    plt.close(fig)

    # Save matrices as CSV
    pearson_df.to_csv(output_dir / 'correlation_matrix.csv')
    ssim_df.to_csv(output_dir / 'ssim_matrix.csv')


# =============================================================================
# Report Generation
# =============================================================================

def generate_report(output_dir: Path, all_results: List[Dict], elapsed_seconds: float):
    """Generate markdown report."""
    logger.info("Generating report...")

    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Count statistics
    total_pairs = len(all_results)
    valid_pairs = len([r for r in all_results if 'error' not in r])
    zero_var_cases = len([r for r in all_results if r.get('note') in ['zero_variance_protein', 'zero_variance_rna']])

    # Get unique samples and proteins
    samples = set(r['sample_id'] for r in all_results if 'error' not in r)
    proteins = set(r['protein'] for r in all_results if 'error' not in r)

    # Compute overall correlations
    pearson_vals = [r['cell_pearson_r'] for r in all_results if 'error' not in r and not np.isnan(r['cell_pearson_r'])]
    ssim_vals = [r['ssim'] for r in all_results if 'error' not in r and not np.isnan(r['ssim'])]

    report = f"""# Protein-RNA Spatial Correlation Report

Generated: {timestamp}
Runtime: {elapsed_seconds/60:.1f} minutes

## Overview

This analysis compares protein immunofluorescence intensity with RNA expression
for 15 matched protein-gene pairs across all samples.

| Metric | Value |
|--------|-------|
| Total protein-sample pairs | {total_pairs} |
| Valid pairs analyzed | {valid_pairs} |
| Samples | {len(samples)} |
| Proteins | {len(proteins)} |
| Zero-variance cases | {zero_var_cases} |

## Protein-Gene Mapping

| Protein | RNA Gene | Status |
|---------|----------|--------|
| CD4 | CD4 | Exact match |
| CD68 | CD68 | Exact match |
| FOXP3 | FOXP3 | Exact match |
| CD3 | CD3D | Primary chain |
| CD8 | CD8A | Primary chain |
| HLA-DR | HLA-DRA | MHC-II alpha |
| KI67 | MKI67 | Standard name |
| PD1 | PDCD1 | Standard name |
| PDL1 | CD274 | Standard name |
| aSMA | ACTA2 | Standard name |
| CD31 | PECAM1 | Standard name |
| **PanCK** | **EPCAM** | **⚠️ PROXY** |
| CD20 | MS4A1 | B cell marker |
| CD45 | PTPRC | Pan-immune |
| CD11c | ITGAX | DC marker |

**Note on PanCK→EPCAM**: Cytokeratins (KRT8, KRT18, KRT19) are not in the RNA panel.
EPCAM is used as an epithelial marker proxy. Correlations should be interpreted with caution.

## Overall Correlation Statistics

### Cell-level Pearson Correlation

| Statistic | Value |
|-----------|-------|
| Mean | {np.mean(pearson_vals):.3f} |
| Std Dev | {np.std(pearson_vals):.3f} |
| Min | {np.min(pearson_vals):.3f} |
| Max | {np.max(pearson_vals):.3f} |
| Median | {np.median(pearson_vals):.3f} |

### Spatial Similarity (SSIM)

| Statistic | Value |
|-----------|-------|
| Mean | {np.mean(ssim_vals):.3f} |
| Std Dev | {np.std(ssim_vals):.3f} |
| Min | {np.min(ssim_vals):.3f} |
| Max | {np.max(ssim_vals):.3f} |
| Median | {np.median(ssim_vals):.3f} |

## Zero-Variance Cases

"""

    # List zero-variance cases
    zero_var_results = [r for r in all_results if r.get('note') in ['zero_variance_protein', 'zero_variance_rna']]
    if zero_var_results:
        report += "| Sample | Protein | Issue |\n|--------|---------|-------|\n"
        for r in zero_var_results:
            report += f"| {r['sample_id']} | {r['protein']} | {r.get('note', 'unknown')} |\n"
    else:
        report += "No zero-variance cases detected.\n"

    report += """

## Output Files

```
results/protein_rna_correlation/
├── sample_panels/           # Individual sample-protein panels
├── protein_summary/         # Per-protein summary figures
├── correlation_details.csv  # All metrics for all pairs
├── correlation_matrix.csv   # Sample × Protein Pearson r
├── ssim_matrix.csv          # Sample × Protein SSIM
├── zero_variance_cases.csv  # Cases with no correlation
├── summary_heatmap.png      # Visual correlation matrix
└── PROTEIN_RNA_REPORT.md
```

## Methods

### Correlation Metrics

1. **Cell-level Pearson r**: Direct correlation between protein intensity and RNA expression for all cells
2. **Cell-level Spearman r**: Rank-based correlation (robust to outliers)
3. **Grid Pearson r**: Correlation on spatially-binned averages (50×50 grid)
4. **SSIM**: Structural Similarity Index measuring spatial pattern similarity

### Normalization

- Protein: Raw intensity values from immunofluorescence
- RNA: Normalized counts from merged_corrected.h5ad
- For SSIM: Both modalities z-score normalized before comparison

### Grid Binning

- Grid size: 50×50 bins
- Binning statistic: Mean value per bin
- Empty bins: Set to 0
"""

    with open(output_dir / 'PROTEIN_RNA_REPORT.md', 'w') as f:
        f.write(report)

    logger.info("  Report saved")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Protein-RNA spatial correlation analysis')
    parser.add_argument('--skip-env-check', action='store_true',
                       help='Skip conda environment check')
    parser.add_argument('--skip-panels', action='store_true',
                       help='Skip individual panel generation (faster)')
    parser.add_argument('--proteins', type=str, nargs='+',
                       help='Analyze specific proteins only')
    parser.add_argument('--samples', type=str, nargs='+',
                       help='Analyze specific samples only')
    args = parser.parse_args()

    # Check environment
    check_environment(args.skip_env_check)

    # Create output directories
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / 'sample_panels').mkdir(exist_ok=True)
    (OUTPUT_DIR / 'protein_summary').mkdir(exist_ok=True)

    # Validate input
    if not INPUT_PATH.exists():
        logger.error(f"Input file not found: {INPUT_PATH}")
        sys.exit(1)

    import time
    start_time = time.time()

    # Load data
    logger.info(f"Loading data from {INPUT_PATH}...")
    adata = sc.read_h5ad(INPUT_PATH)
    logger.info(f"  Loaded {adata.n_obs:,} cells, {adata.n_vars} genes")

    # Get samples to analyze
    all_samples = sorted(adata.obs['sample_id'].unique())
    if args.samples:
        samples = [s for s in args.samples if s in all_samples]
    else:
        samples = all_samples

    logger.info(f"  Analyzing {len(samples)} samples")

    # Get proteins to analyze
    if args.proteins:
        proteins = {p: PROTEIN_TO_GENE[p] for p in args.proteins if p in PROTEIN_TO_GENE}
    else:
        proteins = PROTEIN_TO_GENE

    logger.info(f"  Analyzing {len(proteins)} proteins")

    # Compute global bounds (for consistent grids - not used in current implementation
    # but available for future use)
    global_bounds = compute_global_bounds(adata)
    logger.info(f"  Global bounds: x=[{global_bounds[0]:.0f}, {global_bounds[1]:.0f}], "
               f"y=[{global_bounds[2]:.0f}, {global_bounds[3]:.0f}]")

    # Get sample stages
    sample_stages = adata.obs.groupby('sample_id')['stage'].first().to_dict()

    # Run analysis
    all_results = []

    total_pairs = len(samples) * len(proteins)
    pair_count = 0

    for sample_id in samples:
        logger.info(f"Processing {sample_id}...")

        for protein, gene in proteins.items():
            pair_count += 1
            if pair_count % 50 == 0:
                logger.info(f"  Progress: {pair_count}/{total_pairs} pairs")

            result = analyze_sample_protein(adata, sample_id, protein, gene,
                                           global_bounds, GRID_SIZE)

            # Add stage to result
            result['stage'] = sample_stages.get(sample_id, 'unknown')

            all_results.append(result)

            # Generate panel if not skipped
            if not args.skip_panels and 'error' not in result:
                generate_sample_panel(result, OUTPUT_DIR, result['stage'])

            # Clean up grid data to save memory
            for key in ['_protein_grid', '_rna_grid', '_coords', '_bounds']:
                if key in result:
                    del result[key]

    # Save detailed results
    logger.info("Saving detailed results...")

    # Filter out internal keys
    results_clean = []
    for r in all_results:
        r_clean = {k: v for k, v in r.items() if not k.startswith('_')}
        results_clean.append(r_clean)

    results_df = pd.DataFrame(results_clean)
    results_df.to_csv(OUTPUT_DIR / 'correlation_details.csv', index=False)

    # Save zero-variance cases
    zero_var = results_df[results_df['note'].isin(['zero_variance_protein', 'zero_variance_rna'])]
    zero_var.to_csv(OUTPUT_DIR / 'zero_variance_cases.csv', index=False)

    # Re-run analysis for summary figures (need grids again)
    logger.info("Generating summary figures...")

    # Reload results with grids for protein summaries
    for protein, gene in proteins.items():
        protein_results = []
        for sample_id in samples:
            result = analyze_sample_protein(adata, sample_id, protein, gene,
                                           global_bounds, GRID_SIZE)
            result['stage'] = sample_stages.get(sample_id, 'unknown')
            protein_results.append(result)

        generate_protein_summary(protein_results, protein, gene, OUTPUT_DIR)

        # Clear memory
        del protein_results
        gc.collect()

    # Generate summary heatmap (using saved results)
    all_results_for_heatmap = []
    for sample_id in samples:
        for protein, gene in proteins.items():
            result = analyze_sample_protein(adata, sample_id, protein, gene,
                                           global_bounds, GRID_SIZE)
            result['stage'] = sample_stages.get(sample_id, 'unknown')
            all_results_for_heatmap.append(result)

    generate_summary_heatmap(all_results_for_heatmap, OUTPUT_DIR)

    # Generate report
    elapsed = time.time() - start_time
    generate_report(OUTPUT_DIR, results_clean, elapsed)

    logger.info("="*60)
    logger.info(f"Protein-RNA correlation analysis complete in {elapsed/60:.1f} minutes")
    logger.info(f"Output directory: {OUTPUT_DIR}")
    logger.info(f"Total panels: {len([r for r in results_clean if 'error' not in r])}")
    logger.info("="*60)


if __name__ == '__main__':
    main()
