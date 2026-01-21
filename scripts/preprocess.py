#!/usr/bin/env python3
"""
Spatial Biology Hackathon 2026 - Preprocessing Pipeline
========================================================

QC filtering and preprocessing for spatial transcriptomics data.
Based on scanpy workflow patterns from omicverse/scSLAT.

Author: Max Van Belkum
Date: 2026-01-20
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import yaml
from pathlib import Path
from typing import Optional, Dict, Tuple
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Configuration
# =============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
CONFIG_PATH = PROJECT_ROOT / "config"

# Default thresholds (platform-specific in config/qc_thresholds.yaml)
DEFAULT_THRESHOLDS = {
    'Visium': {
        'min_genes': 200,
        'max_genes': 8000,
        'min_counts': 500,
        'max_counts': 50000,
        'max_pct_mt': 25.0,
        'min_cells_per_gene': 3,
    },
    'G4X': {
        'min_genes': 50,
        'max_genes': 500,
        'min_counts': 100,
        'max_counts': 10000,
        'max_pct_mt': 25.0,
        'min_cells_per_gene': 3,
    }
}


def load_thresholds(platform: str, config_path: Optional[Path] = None) -> Dict:
    """
    Load QC thresholds from config or use defaults.

    Args:
        platform: Platform name (Visium, G4X, etc.)
        config_path: Path to qc_thresholds.yaml

    Returns:
        Dictionary with threshold values
    """
    if config_path is None:
        config_path = CONFIG_PATH / "qc_thresholds.yaml"

    if config_path.exists():
        with open(config_path) as f:
            config = yaml.safe_load(f)
        if platform in config:
            return config[platform]

    # Fall back to defaults
    return DEFAULT_THRESHOLDS.get(platform, DEFAULT_THRESHOLDS['Visium'])


# =============================================================================
# Preprocessing Functions
# =============================================================================

def calculate_qc(adata: ad.AnnData) -> ad.AnnData:
    """
    Calculate QC metrics without filtering.

    Adds:
        - n_genes_by_counts
        - total_counts
        - pct_counts_mt (if MT genes present)

    Args:
        adata: AnnData object

    Returns:
        AnnData with QC metrics added
    """
    # Adjust percent_top for samples with few features (G4X has ~387)
    n_features = adata.n_vars
    if n_features < 100:
        percent_top = [10, 20, 50]
    elif n_features < 500:
        percent_top = [50, 100, 200]
    else:
        percent_top = [50, 100, 200, 500]

    # Basic QC metrics
    sc.pp.calculate_qc_metrics(adata, percent_top=percent_top, inplace=True)

    # Mitochondrial percentage
    mito_genes = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
    n_mito = mito_genes.sum()

    if n_mito > 0:
        adata.var['mt'] = mito_genes
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
        print(f"  Found {n_mito} mitochondrial genes")
    else:
        adata.obs['pct_counts_mt'] = 0.0
        print(f"  No mitochondrial genes found")

    return adata


def filter_cells(
    adata: ad.AnnData,
    min_genes: int = 200,
    max_genes: int = 8000,
    min_counts: int = 500,
    max_counts: int = 50000,
    max_pct_mt: float = 25.0,
) -> ad.AnnData:
    """
    Filter cells based on QC metrics.

    Args:
        adata: AnnData with QC metrics
        min_genes: Minimum genes per cell
        max_genes: Maximum genes per cell
        min_counts: Minimum UMIs per cell
        max_counts: Maximum UMIs per cell
        max_pct_mt: Maximum mitochondrial percentage

    Returns:
        Filtered AnnData (copy)
    """
    n_before = adata.n_obs

    # Ensure QC metrics exist
    if 'n_genes_by_counts' not in adata.obs.columns:
        adata = calculate_qc(adata)

    # Apply filters
    mask = (
        (adata.obs['n_genes_by_counts'] >= min_genes) &
        (adata.obs['n_genes_by_counts'] <= max_genes) &
        (adata.obs['total_counts'] >= min_counts) &
        (adata.obs['total_counts'] <= max_counts) &
        (adata.obs['pct_counts_mt'] <= max_pct_mt)
    )

    adata_filtered = adata[mask].copy()

    n_after = adata_filtered.n_obs
    pct_retained = (n_after / n_before) * 100 if n_before > 0 else 0

    print(f"  Cell filtering: {n_before:,} -> {n_after:,} ({pct_retained:.1f}% retained)")

    return adata_filtered


def filter_genes(
    adata: ad.AnnData,
    min_cells: int = 3,
) -> ad.AnnData:
    """
    Filter genes expressed in too few cells.

    Args:
        adata: AnnData object
        min_cells: Minimum cells expressing gene

    Returns:
        Filtered AnnData (copy)
    """
    n_before = adata.n_vars

    sc.pp.filter_genes(adata, min_cells=min_cells)

    n_after = adata.n_vars
    pct_retained = (n_after / n_before) * 100 if n_before > 0 else 0

    print(f"  Gene filtering: {n_before:,} -> {n_after:,} ({pct_retained:.1f}% retained)")

    return adata


def preprocess_spatial(
    adata: ad.AnnData,
    platform: str = "Visium",
    thresholds: Optional[Dict] = None,
    normalize: bool = False,
    log1p: bool = False,
) -> Tuple[ad.AnnData, Dict]:
    """
    Full preprocessing pipeline for spatial data.

    Args:
        adata: Raw AnnData object
        platform: Platform name for threshold lookup
        thresholds: Override thresholds (or loads from config)
        normalize: Whether to normalize counts
        log1p: Whether to log-transform

    Returns:
        Tuple of (processed AnnData, filtering stats dict)
    """
    print(f"Preprocessing ({platform})...")

    # Get thresholds
    if thresholds is None:
        thresholds = load_thresholds(platform)

    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Calculate QC
    adata = calculate_qc(adata)

    # Filter cells
    adata = filter_cells(
        adata,
        min_genes=thresholds['min_genes'],
        max_genes=thresholds['max_genes'],
        min_counts=thresholds['min_counts'],
        max_counts=thresholds.get('max_counts', 50000),
        max_pct_mt=thresholds['max_pct_mt'],
    )

    # Filter genes
    adata = filter_genes(adata, min_cells=thresholds['min_cells_per_gene'])

    # Optional normalization
    if normalize:
        sc.pp.normalize_total(adata, target_sum=1e4)
        print("  Normalized to 10,000 counts per cell")

    if log1p:
        sc.pp.log1p(adata)
        print("  Log-transformed")

    # Stats
    stats = {
        'cells_before': n_cells_before,
        'cells_after': adata.n_obs,
        'cells_retained_pct': (adata.n_obs / n_cells_before) * 100,
        'genes_before': n_genes_before,
        'genes_after': adata.n_vars,
        'genes_retained_pct': (adata.n_vars / n_genes_before) * 100,
    }

    print(f"  Final: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    return adata, stats


# =============================================================================
# Batch Processing
# =============================================================================

def preprocess_batch(
    adatas: Dict[str, ad.AnnData],
    platform_map: Optional[Dict[str, str]] = None,
) -> Tuple[Dict[str, ad.AnnData], pd.DataFrame]:
    """
    Preprocess multiple samples.

    Args:
        adatas: Dictionary of sample_name -> AnnData
        platform_map: Optional mapping of sample_name -> platform

    Returns:
        Tuple of (processed adatas dict, stats DataFrame)
    """
    processed = {}
    stats_list = []

    for sample_name, adata in adatas.items():
        print(f"\n{'='*50}")
        print(f"Processing: {sample_name}")
        print(f"{'='*50}")

        # Determine platform
        platform = "Visium"  # default
        if platform_map and sample_name in platform_map:
            platform = platform_map[sample_name]
        elif 'platform' in adata.obs.columns:
            platform = adata.obs['platform'].iloc[0]

        try:
            adata_proc, stats = preprocess_spatial(adata, platform=platform)
            stats['sample'] = sample_name
            stats['platform'] = platform
            processed[sample_name] = adata_proc
            stats_list.append(stats)
        except Exception as e:
            print(f"  ERROR: {e}")

    stats_df = pd.DataFrame(stats_list)

    return processed, stats_df


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

    from pathlib import Path
    from scripts import load_sample

    print("=" * 60)
    print("PREPROCESSING TEST")
    print("=" * 60)

    # Test on one sample
    try:
        adata = load_sample("YP03A")
        adata_proc, stats = preprocess_spatial(adata, platform="Visium")
        print(f"\nStats: {stats}")
    except Exception as e:
        print(f"Error: {e}")
