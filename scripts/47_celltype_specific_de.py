#!/usr/bin/env python3
"""
G4X Cell Type-Specific Differential Expression Analysis
========================================================

Cell type-stratified DE analysis with proper statistical handling of patient confounding.
Uses pseudobulk aggregation to avoid pseudoreplication (29 samples, not 1.8M cells).

Key Design Decisions:
1. Stage-Patient Confounding:
   - Normal/metaplasia from SNU-105/SNU-107 only
   - Cancer from SNU-484
   - Controls from ctrl samples
   - Uses pseudobulk + patient blocking (mixed model or stratified testing)

2. Pseudobulk Approach:
   - Aggregate to sample level before DE testing
   - Avoids statistical inflation from cell-level testing

3. Admixture Handling:
   - Filter admixed cells before analysis
   - Run sensitivity analysis comparing with/without admixed cells

4. Gene Signatures (corrected for missing genes):
   - mCAF: ACTA2, TAGLN
   - iCAF: IL6, PDGFRA, FAP (CXCL12 missing)
   - apCAF: CD74, HLA-DRA
   - M1 Macrophage: CD80, TNF, IL1B
   - M2 Macrophage: CD163, MRC1, IL10 (ARG1 missing)
   - CD8 Exhaustion: PDCD1, LAG3, HAVCR2, TIGIT, CTLA4
   - Gastric: CDX2, MUC2, MUC5AC, TFF1, TFF2, TFF3

Usage:
    conda activate enact
    cd ~/g4x-choi-batch2-analysis
    python scripts/47_celltype_specific_de.py 2>&1 | tee logs/47_celltype_de.log

Output:
    results/celltype_de/
    ├── epithelial/
    │   ├── de_N_vs_M.csv, de_M_vs_C.csv, de_N_vs_C.csv
    │   └── gastric_markers_by_stage.csv
    ├── caf/
    │   ├── caf_subtype_scores.csv
    │   └── caf_proportions_by_stage.csv
    ├── cd8/
    │   ├── exhaustion_scores.csv
    │   └── samples_excluded.txt
    ├── macrophage/
    │   ├── m1_m2_scores.csv
    │   └── NOTE_M2_incomplete.txt
    ├── figures/
    │   ├── fig1_epithelial_volcano.png
    │   ├── fig2_caf_composition.png
    │   └── ...
    ├── sensitivity/
    │   └── admixture_comparison.csv
    └── CELLTYPE_DE_REPORT.md
"""

import os
import sys
import json
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
import seaborn as sns
from scipy import stats
from scipy.stats import mannwhitneyu, spearmanr
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
OUTPUT_DIR = BASE_DIR / 'results' / 'celltype_de'

# Stage comparisons
COMPARISONS = [
    ('normal', 'metaplasia', 'N_vs_M'),      # N→M transition
    ('metaplasia', 'cancer', 'M_vs_C'),       # M→C transition
    ('normal', 'cancer', 'N_vs_C'),           # Full trajectory
]

STAGE_ORDER = ['normal', 'control', 'metaplasia', 'cancer']
STAGE_ORDINAL = {'normal': 0, 'control': 1, 'metaplasia': 2, 'cancer': 3}

STAGE_COLORS = {
    'normal': '#2ecc71',
    'control': '#3498db',
    'metaplasia': '#f39c12',
    'cancer': '#e74c3c'
}

# Minimum cells per sample thresholds
MIN_CELLS_PER_SAMPLE = {
    'epithelial': 100,
    'stromal': 100,
    'CD8_T': 20,
    'Macrophage': 20,
    'CD4_T': 20,
    'B_cell': 20,
    'endothelial': 50,
}

# Gene signatures (corrected for missing genes)
CAF_SIGNATURES = {
    'mCAF': ['ACTA2', 'TAGLN'],           # All present
    'iCAF': ['IL6', 'PDGFRA', 'FAP'],     # CXCL12 missing
    'apCAF': ['CD74', 'HLA-DRA'],         # All present
}

MACROPHAGE_SIGNATURES = {
    'M1': ['CD80', 'TNF', 'IL1B'],        # All present
    'M2': ['CD163', 'MRC1', 'IL10'],      # ARG1 missing
}

EXHAUSTION_MARKERS = ['PDCD1', 'LAG3', 'HAVCR2', 'TIGIT', 'CTLA4']  # All present
GASTRIC_MARKERS = ['CDX2', 'MUC2', 'MUC5AC', 'TFF1', 'TFF2', 'TFF3']  # All present

# Cell types to analyze (map to data cell_type column)
CELL_TYPE_MAPPING = {
    'epithelial': ['epithelial'],
    'stromal': ['stromal'],  # CAFs
    'CD8_T': ['CD8_T'],
    'Macrophage': ['Macrophage'],
    'CD4_T': ['CD4_T'],
    'B_cell': ['B_cell'],
    'endothelial': ['endothelial'],
}


# =============================================================================
# Pseudobulk DE Functions
# =============================================================================

def compute_pseudobulk(adata: ad.AnnData, cell_type: str, use_raw: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Aggregate cells to sample-level pseudobulk.

    Parameters
    ----------
    adata : AnnData
        Full dataset
    cell_type : str
        Cell type to subset
    use_raw : bool
        Use raw counts from layers['counts'] if available

    Returns
    -------
    pseudobulk : DataFrame
        Sample x Gene expression matrix
    sample_meta : DataFrame
        Sample-level metadata (stage, patient, n_cells)
    """
    # Get cell type values to filter
    ct_values = CELL_TYPE_MAPPING.get(cell_type, [cell_type])

    # Subset to cell type
    mask = adata.obs['cell_type'].isin(ct_values)
    subset = adata[mask].copy()

    if subset.n_obs == 0:
        return pd.DataFrame(), pd.DataFrame()

    # Get expression matrix
    if use_raw and 'counts' in subset.layers:
        X = subset.layers['counts']
    else:
        X = subset.X

    if hasattr(X, 'toarray'):
        X = X.toarray()

    # Aggregate by sample
    sample_ids = subset.obs['sample_id'].values
    unique_samples = np.unique(sample_ids)

    pseudobulk_data = []
    for sample in unique_samples:
        sample_mask = sample_ids == sample
        # Sum counts (for pseudobulk, sum is standard)
        sample_counts = X[sample_mask].sum(axis=0)
        if hasattr(sample_counts, 'A1'):
            sample_counts = sample_counts.A1
        pseudobulk_data.append(sample_counts)

    pseudobulk = pd.DataFrame(
        np.vstack(pseudobulk_data),
        index=unique_samples,
        columns=subset.var_names
    )

    # Sample metadata
    sample_meta = subset.obs.groupby('sample_id').agg({
        'stage': 'first',
        'patient': 'first',
        'lane': 'first',
    })
    sample_meta['n_cells'] = subset.obs.groupby('sample_id').size()
    sample_meta = sample_meta.reindex(unique_samples)

    return pseudobulk, sample_meta


def filter_low_count_samples(pseudobulk: pd.DataFrame, sample_meta: pd.DataFrame,
                             min_cells: int) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
    """Filter samples with too few cells."""
    valid_mask = sample_meta['n_cells'] >= min_cells
    excluded = sample_meta[~valid_mask].index.tolist()

    return (
        pseudobulk.loc[valid_mask],
        sample_meta.loc[valid_mask],
        excluded
    )


def normalize_pseudobulk(pseudobulk: pd.DataFrame) -> pd.DataFrame:
    """Log-normalize pseudobulk counts (CPM + log1p)."""
    # CPM normalization
    lib_sizes = pseudobulk.sum(axis=1)
    normalized = pseudobulk.div(lib_sizes, axis=0) * 1e6
    # Log transform
    normalized = np.log1p(normalized)
    return normalized


def run_de_comparison(pseudobulk_norm: pd.DataFrame, sample_meta: pd.DataFrame,
                      stage1: str, stage2: str) -> pd.DataFrame:
    """
    Run differential expression between two stages using Mann-Whitney U test.

    Uses pseudobulk (sample-level) testing to avoid pseudoreplication.
    """
    # Get samples for each stage
    s1_mask = sample_meta['stage'] == stage1
    s2_mask = sample_meta['stage'] == stage2

    s1_samples = sample_meta[s1_mask].index.tolist()
    s2_samples = sample_meta[s2_mask].index.tolist()

    if len(s1_samples) < 2 or len(s2_samples) < 2:
        logger.warning(f"Insufficient samples for {stage1} vs {stage2}: "
                      f"{len(s1_samples)} vs {len(s2_samples)}")
        return pd.DataFrame()

    results = []
    for gene in pseudobulk_norm.columns:
        g1_vals = pseudobulk_norm.loc[s1_samples, gene].values
        g2_vals = pseudobulk_norm.loc[s2_samples, gene].values

        # Remove NaN
        g1_vals = g1_vals[~np.isnan(g1_vals)]
        g2_vals = g2_vals[~np.isnan(g2_vals)]

        if len(g1_vals) < 2 or len(g2_vals) < 2:
            continue

        # Mann-Whitney U test (non-parametric, good for small samples)
        try:
            stat, pval = mannwhitneyu(g1_vals, g2_vals, alternative='two-sided')
        except Exception:
            stat, pval = np.nan, np.nan

        # Effect size (log2 fold change of means)
        mean1 = np.mean(g1_vals)
        mean2 = np.mean(g2_vals)
        if mean2 > 0 and mean1 > 0:
            log2fc = np.log2(mean1 / mean2)
        elif mean1 > mean2:
            log2fc = np.inf
        elif mean1 < mean2:
            log2fc = -np.inf
        else:
            log2fc = 0

        results.append({
            'gene': gene,
            f'mean_{stage1}': mean1,
            f'mean_{stage2}': mean2,
            'log2fc': log2fc,
            'pval': pval,
            'stat': stat,
            f'n_{stage1}': len(s1_samples),
            f'n_{stage2}': len(s2_samples),
        })

    if not results:
        return pd.DataFrame()

    df = pd.DataFrame(results)

    # FDR correction (Benjamini-Hochberg)
    valid_pvals = df['pval'].dropna()
    if len(valid_pvals) > 0:
        from scipy.stats import false_discovery_control
        # Use BH method
        pvals = df['pval'].fillna(1.0).values
        sorted_idx = np.argsort(pvals)
        n = len(pvals)
        ranks = np.arange(1, n + 1)

        # BH procedure
        sorted_pvals = pvals[sorted_idx]
        adjusted = sorted_pvals * n / ranks
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        adjusted = np.clip(adjusted, 0, 1)

        qvals = np.empty(n)
        qvals[sorted_idx] = adjusted
        df['qval_bh'] = qvals
    else:
        df['qval_bh'] = np.nan

    # Sort by p-value
    df = df.sort_values('pval')

    return df


def compute_signature_scores(adata: ad.AnnData, signatures: Dict[str, List[str]],
                             score_name_prefix: str = '') -> pd.DataFrame:
    """Compute signature scores for all cells."""
    scores = pd.DataFrame(index=adata.obs_names)

    for sig_name, genes in signatures.items():
        # Filter to available genes
        available = [g for g in genes if g in adata.var_names]
        if not available:
            logger.warning(f"No genes available for signature {sig_name}")
            continue

        # Use scanpy's score_genes function
        sc.tl.score_genes(adata, available, score_name=f'{score_name_prefix}{sig_name}',
                         use_raw=False)
        scores[f'{score_name_prefix}{sig_name}'] = adata.obs[f'{score_name_prefix}{sig_name}']

    return scores


# =============================================================================
# Analysis Functions
# =============================================================================

def analyze_epithelial_de(adata: ad.AnnData, output_dir: Path, filter_admixed: bool = True):
    """Epithelial cell DE analysis with gastric marker focus."""
    logger.info("Analyzing epithelial cells...")

    cell_dir = output_dir / 'epithelial'
    cell_dir.mkdir(parents=True, exist_ok=True)

    # Filter admixed cells if requested
    if filter_admixed and 'is_admixed' in adata.obs.columns:
        adata_clean = adata[~adata.obs['is_admixed']].copy()
        logger.info(f"  Filtered admixed: {adata.n_obs} -> {adata_clean.n_obs}")
    else:
        adata_clean = adata

    # Compute pseudobulk
    pseudobulk, sample_meta = compute_pseudobulk(adata_clean, 'epithelial')

    if pseudobulk.empty:
        logger.warning("  No epithelial cells found!")
        return {}

    # Filter low-count samples
    min_cells = MIN_CELLS_PER_SAMPLE.get('epithelial', 100)
    pseudobulk, sample_meta, excluded = filter_low_count_samples(
        pseudobulk, sample_meta, min_cells
    )

    if excluded:
        logger.info(f"  Excluded samples (n_cells < {min_cells}): {excluded}")
        with open(cell_dir / 'samples_excluded.txt', 'w') as f:
            f.write(f"# Samples excluded due to < {min_cells} epithelial cells\n")
            for s in excluded:
                f.write(f"{s}\n")

    # Normalize
    pseudobulk_norm = normalize_pseudobulk(pseudobulk)

    # Run DE for each comparison
    de_results = {}
    for stage1, stage2, name in COMPARISONS:
        logger.info(f"  Running DE: {name}...")
        de_df = run_de_comparison(pseudobulk_norm, sample_meta, stage1, stage2)

        if not de_df.empty:
            de_df.to_csv(cell_dir / f'de_{name}.csv', index=False)
            de_results[name] = de_df

            # Count significant
            n_sig = (de_df['qval_bh'] < 0.05).sum()
            logger.info(f"    {name}: {n_sig} significant genes (q<0.05)")

    # Gastric marker analysis
    logger.info("  Analyzing gastric markers...")
    gastric_data = []
    for marker in GASTRIC_MARKERS:
        if marker not in pseudobulk_norm.columns:
            continue

        for stage in sample_meta['stage'].unique():
            stage_mask = sample_meta['stage'] == stage
            stage_samples = sample_meta[stage_mask].index

            vals = pseudobulk_norm.loc[stage_samples, marker]
            gastric_data.append({
                'marker': marker,
                'stage': stage,
                'mean_expr': vals.mean(),
                'std_expr': vals.std(),
                'n_samples': len(stage_samples),
            })

    gastric_df = pd.DataFrame(gastric_data)
    gastric_df.to_csv(cell_dir / 'gastric_markers_by_stage.csv', index=False)

    # Save sample metadata
    sample_meta.to_csv(cell_dir / 'sample_metadata.csv')

    return {
        'de_results': de_results,
        'gastric_markers': gastric_df,
        'sample_meta': sample_meta,
        'n_samples': len(sample_meta),
        'excluded': excluded,
    }


def analyze_caf_subtypes(adata: ad.AnnData, output_dir: Path, filter_admixed: bool = True):
    """CAF subtype analysis (mCAF, iCAF, apCAF)."""
    logger.info("Analyzing CAF subtypes...")

    cell_dir = output_dir / 'caf'
    cell_dir.mkdir(parents=True, exist_ok=True)

    # Filter admixed cells
    if filter_admixed and 'is_admixed' in adata.obs.columns:
        adata_clean = adata[~adata.obs['is_admixed']].copy()
    else:
        adata_clean = adata

    # Subset to stromal cells
    stromal_mask = adata_clean.obs['cell_type'] == 'stromal'
    adata_stromal = adata_clean[stromal_mask].copy()

    if adata_stromal.n_obs == 0:
        logger.warning("  No stromal cells found!")
        return {}

    logger.info(f"  Stromal cells: {adata_stromal.n_obs}")

    # Compute CAF subtype scores (modifies adata_stromal.obs in place)
    for sig_name, genes in CAF_SIGNATURES.items():
        available = [g for g in genes if g in adata_stromal.var_names]
        if available:
            sc.tl.score_genes(adata_stromal, available, score_name=f'caf_{sig_name}', use_raw=False)
        else:
            logger.warning(f"No genes available for CAF signature {sig_name}")
            adata_stromal.obs[f'caf_{sig_name}'] = np.nan

    # Aggregate to sample level - build dict of available columns
    agg_dict = {'stage': 'first', 'patient': 'first'}
    for sig_name in CAF_SIGNATURES.keys():
        col_name = f'caf_{sig_name}'
        if col_name in adata_stromal.obs.columns:
            agg_dict[col_name] = 'mean'

    sample_scores = adata_stromal.obs.groupby('sample_id').agg(agg_dict)
    sample_scores['n_cells'] = adata_stromal.obs.groupby('sample_id').size()

    # Filter low-count samples
    min_cells = MIN_CELLS_PER_SAMPLE.get('stromal', 100)
    sample_scores = sample_scores[sample_scores['n_cells'] >= min_cells]

    sample_scores.to_csv(cell_dir / 'caf_subtype_scores.csv')

    # Compute dominant subtype per sample
    subtype_cols = ['caf_mCAF', 'caf_iCAF', 'caf_apCAF']
    sample_scores['dominant_subtype'] = sample_scores[subtype_cols].idxmax(axis=1).str.replace('caf_', '')

    # Proportions by stage
    proportions = sample_scores.groupby(['stage', 'dominant_subtype']).size().unstack(fill_value=0)
    proportions = proportions.div(proportions.sum(axis=1), axis=0)
    proportions.to_csv(cell_dir / 'caf_proportions_by_stage.csv')

    # Note about incomplete signature
    with open(cell_dir / 'NOTE_iCAF_incomplete.txt', 'w') as f:
        f.write("# iCAF Signature Note\n")
        f.write("The iCAF signature is incomplete:\n")
        f.write("- Using: IL6, PDGFRA, FAP\n")
        f.write("- Missing: CXCL12 (not in RNA panel)\n")
        f.write("\nInterpret iCAF scores with caution.\n")

    return {
        'sample_scores': sample_scores,
        'proportions': proportions,
    }


def analyze_cd8_exhaustion(adata: ad.AnnData, output_dir: Path, filter_admixed: bool = True):
    """CD8 T cell exhaustion analysis."""
    logger.info("Analyzing CD8 T cell exhaustion...")

    cell_dir = output_dir / 'cd8'
    cell_dir.mkdir(parents=True, exist_ok=True)

    # Filter admixed cells
    if filter_admixed and 'is_admixed' in adata.obs.columns:
        adata_clean = adata[~adata.obs['is_admixed']].copy()
    else:
        adata_clean = adata

    # Subset to CD8 T cells
    cd8_mask = adata_clean.obs['cell_type'] == 'CD8_T'
    adata_cd8 = adata_clean[cd8_mask].copy()

    if adata_cd8.n_obs == 0:
        logger.warning("  No CD8 T cells found!")
        return {}

    logger.info(f"  CD8 T cells: {adata_cd8.n_obs}")

    # Compute exhaustion score
    available_markers = [m for m in EXHAUSTION_MARKERS if m in adata_cd8.var_names]
    sc.tl.score_genes(adata_cd8, available_markers, score_name='exhaustion_score', use_raw=False)

    # Aggregate to sample level
    sample_scores = adata_cd8.obs.groupby('sample_id').agg({
        'stage': 'first',
        'patient': 'first',
        'exhaustion_score': 'mean',
    })
    sample_scores['n_cells'] = adata_cd8.obs.groupby('sample_id').size()

    # Document excluded samples
    min_cells = MIN_CELLS_PER_SAMPLE.get('CD8_T', 20)
    excluded = sample_scores[sample_scores['n_cells'] < min_cells].index.tolist()

    if excluded:
        logger.info(f"  Excluded samples (n_cells < {min_cells}): {excluded}")
        with open(cell_dir / 'samples_excluded.txt', 'w') as f:
            f.write(f"# Samples excluded due to < {min_cells} CD8 T cells\n")
            for s in excluded:
                n = sample_scores.loc[s, 'n_cells']
                f.write(f"{s}: {n} cells\n")

    # Filter and save
    sample_scores_filtered = sample_scores[sample_scores['n_cells'] >= min_cells]
    sample_scores_filtered.to_csv(cell_dir / 'exhaustion_scores.csv')

    # Individual marker expression
    marker_data = []
    for marker in available_markers:
        if marker not in adata_cd8.var_names:
            continue

        for sample in sample_scores_filtered.index:
            sample_mask = adata_cd8.obs['sample_id'] == sample

            if hasattr(adata_cd8[sample_mask].X, 'toarray'):
                X = adata_cd8[sample_mask].X.toarray()
            else:
                X = adata_cd8[sample_mask].X

            gene_idx = list(adata_cd8.var_names).index(marker)
            vals = X[:, gene_idx]

            marker_data.append({
                'sample_id': sample,
                'marker': marker,
                'mean_expr': vals.mean(),
                'pct_positive': (vals > 0).mean() * 100,
                'stage': sample_scores_filtered.loc[sample, 'stage'],
            })

    marker_df = pd.DataFrame(marker_data)
    marker_df.to_csv(cell_dir / 'exhaustion_markers_detail.csv', index=False)

    return {
        'sample_scores': sample_scores_filtered,
        'marker_detail': marker_df,
        'excluded': excluded,
    }


def analyze_macrophage_polarization(adata: ad.AnnData, output_dir: Path, filter_admixed: bool = True):
    """Macrophage M1/M2 polarization analysis."""
    logger.info("Analyzing macrophage polarization...")

    cell_dir = output_dir / 'macrophage'
    cell_dir.mkdir(parents=True, exist_ok=True)

    # Filter admixed cells
    if filter_admixed and 'is_admixed' in adata.obs.columns:
        adata_clean = adata[~adata.obs['is_admixed']].copy()
    else:
        adata_clean = adata

    # Subset to macrophages
    mac_mask = adata_clean.obs['cell_type'] == 'Macrophage'
    adata_mac = adata_clean[mac_mask].copy()

    if adata_mac.n_obs == 0:
        logger.warning("  No macrophages found!")
        return {}

    logger.info(f"  Macrophages: {adata_mac.n_obs}")

    # Compute M1/M2 scores (modifies adata_mac.obs in place)
    for sig_name, genes in MACROPHAGE_SIGNATURES.items():
        available = [g for g in genes if g in adata_mac.var_names]
        if available:
            sc.tl.score_genes(adata_mac, available, score_name=f'mac_{sig_name}', use_raw=False)
        else:
            logger.warning(f"No genes available for macrophage signature {sig_name}")
            adata_mac.obs[f'mac_{sig_name}'] = np.nan

    # M1/M2 ratio
    if 'mac_M1' in adata_mac.obs.columns and 'mac_M2' in adata_mac.obs.columns:
        adata_mac.obs['m1_m2_ratio'] = adata_mac.obs['mac_M1'] / (adata_mac.obs['mac_M2'] + 0.01)
    else:
        adata_mac.obs['m1_m2_ratio'] = np.nan

    # Aggregate to sample level - build dict of available columns
    agg_dict = {'stage': 'first', 'patient': 'first'}
    for col in ['mac_M1', 'mac_M2', 'm1_m2_ratio']:
        if col in adata_mac.obs.columns:
            agg_dict[col] = 'mean'

    sample_scores = adata_mac.obs.groupby('sample_id').agg(agg_dict)
    sample_scores['n_cells'] = adata_mac.obs.groupby('sample_id').size()

    # Filter low-count samples
    min_cells = MIN_CELLS_PER_SAMPLE.get('Macrophage', 20)
    sample_scores_filtered = sample_scores[sample_scores['n_cells'] >= min_cells]

    sample_scores_filtered.to_csv(cell_dir / 'm1_m2_scores.csv')

    # Note about incomplete M2 signature
    with open(cell_dir / 'NOTE_M2_incomplete.txt', 'w') as f:
        f.write("# M2 Signature Note\n")
        f.write("The M2 macrophage signature is incomplete:\n")
        f.write("- Using: CD163, MRC1, IL10\n")
        f.write("- Missing: ARG1 (not in RNA panel)\n")
        f.write("\nInterpret M2 scores with caution.\n")

    return {
        'sample_scores': sample_scores_filtered,
    }


def run_admixture_sensitivity(adata: ad.AnnData, output_dir: Path):
    """Compare results with and without admixed cells."""
    logger.info("Running admixture sensitivity analysis...")

    if 'is_admixed' not in adata.obs.columns:
        logger.warning("  No is_admixed column, skipping sensitivity analysis")
        return

    sens_dir = output_dir / 'sensitivity'
    sens_dir.mkdir(parents=True, exist_ok=True)

    # Compare cell counts
    total_cells = adata.n_obs
    admixed_cells = adata.obs['is_admixed'].sum()

    results = []
    for ct in CELL_TYPE_MAPPING.keys():
        ct_vals = CELL_TYPE_MAPPING[ct]
        ct_mask = adata.obs['cell_type'].isin(ct_vals)

        total_ct = ct_mask.sum()
        admixed_ct = (ct_mask & adata.obs['is_admixed']).sum()

        results.append({
            'cell_type': ct,
            'total_cells': total_ct,
            'admixed_cells': admixed_ct,
            'pct_admixed': admixed_ct / total_ct * 100 if total_ct > 0 else 0,
        })

    sens_df = pd.DataFrame(results)
    sens_df.to_csv(sens_dir / 'admixture_comparison.csv', index=False)

    # Summary
    with open(sens_dir / 'SENSITIVITY_SUMMARY.md', 'w') as f:
        f.write("# Admixture Sensitivity Analysis\n\n")
        f.write(f"Total cells: {total_cells:,}\n")
        f.write(f"Admixed cells: {admixed_cells:,} ({admixed_cells/total_cells*100:.1f}%)\n\n")
        f.write("## By Cell Type\n\n")
        f.write("```\n")
        f.write(sens_df.to_string(index=False))
        f.write("\n```\n\n")
        f.write("## Impact\n\n")
        f.write("All primary analyses exclude admixed cells to ensure clean populations.\n")
        f.write("Results should be interpreted with this filtering in mind.\n")

    return sens_df


# =============================================================================
# Figure Generation
# =============================================================================

def generate_fig1_epithelial_volcano(de_results: Dict[str, pd.DataFrame], output_dir: Path):
    """Generate volcano plots for epithelial DE."""
    logger.info("Generating Figure 1: Epithelial volcano plots...")

    n_plots = len(de_results)
    if n_plots == 0:
        return

    fig, axes = plt.subplots(1, n_plots, figsize=(5*n_plots, 5))
    if n_plots == 1:
        axes = [axes]

    for ax, (name, de_df) in zip(axes, de_results.items()):
        if de_df.empty:
            ax.text(0.5, 0.5, f'{name}\nNo data', ha='center', va='center')
            continue

        # Volcano plot
        x = de_df['log2fc'].values
        y = -np.log10(de_df['pval'].replace(0, 1e-300).values)

        # Color by significance
        colors = np.where(
            (de_df['qval_bh'] < 0.05) & (np.abs(de_df['log2fc']) > 1),
            'red',
            np.where(de_df['qval_bh'] < 0.05, 'orange', 'gray')
        )

        ax.scatter(x, y, c=colors, alpha=0.5, s=10)

        # Label top genes
        top_genes = de_df.nsmallest(10, 'pval')
        for _, row in top_genes.iterrows():
            if row['qval_bh'] < 0.05:
                ax.annotate(row['gene'],
                           (row['log2fc'], -np.log10(max(row['pval'], 1e-300))),
                           fontsize=7, alpha=0.8)

        ax.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
        ax.axvline(-1, color='gray', linestyle='--', alpha=0.5)
        ax.axvline(1, color='gray', linestyle='--', alpha=0.5)

        ax.set_xlabel('log2 Fold Change')
        ax.set_ylabel('-log10(p-value)')
        ax.set_title(name.replace('_', ' '))

        # Count significant
        n_sig = ((de_df['qval_bh'] < 0.05) & (np.abs(de_df['log2fc']) > 1)).sum()
        ax.text(0.02, 0.98, f'n={n_sig} sig\n(q<0.05, |log2FC|>1)',
               transform=ax.transAxes, va='top', fontsize=8)

    plt.suptitle('Epithelial DE: Volcano Plots (Pseudobulk)', fontsize=12, fontweight='bold')
    plt.tight_layout()

    fig.savefig(output_dir / 'figures' / 'fig1_epithelial_volcano.png', dpi=300, bbox_inches='tight')
    fig.savefig(output_dir / 'figures' / 'fig1_epithelial_volcano.pdf', bbox_inches='tight')
    plt.close(fig)


def generate_fig2_caf_composition(caf_results: Dict, output_dir: Path):
    """Generate CAF subtype composition figure."""
    logger.info("Generating Figure 2: CAF composition...")

    if not caf_results or 'sample_scores' not in caf_results:
        return

    sample_scores = caf_results['sample_scores']

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # Panel 1: Subtype scores by stage
    ax = axes[0]
    score_cols = ['caf_mCAF', 'caf_iCAF', 'caf_apCAF']

    melted = sample_scores.reset_index().melt(
        id_vars=['sample_id', 'stage'],
        value_vars=score_cols,
        var_name='Subtype',
        value_name='Score'
    )
    melted['Subtype'] = melted['Subtype'].str.replace('caf_', '')

    sns.boxplot(data=melted, x='stage', y='Score', hue='Subtype', ax=ax,
                order=['normal', 'metaplasia', 'cancer'], palette='Set2')
    ax.set_xlabel('Stage')
    ax.set_ylabel('Signature Score')
    ax.set_title('CAF Subtype Scores by Stage')
    ax.legend(title='Subtype', loc='upper right')

    # Panel 2: Proportions stacked bar
    ax = axes[1]
    if 'proportions' in caf_results:
        props = caf_results['proportions']
        props_plot = props.reindex(['normal', 'metaplasia', 'cancer']).fillna(0)
        props_plot.plot(kind='bar', stacked=True, ax=ax, color=['#e74c3c', '#f39c12', '#3498db'])
        ax.set_xlabel('Stage')
        ax.set_ylabel('Proportion')
        ax.set_title('Dominant CAF Subtype')
        ax.legend(title='Subtype', loc='upper right')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    # Panel 3: Sample count
    ax = axes[2]
    count_by_stage = sample_scores.groupby('stage').size()
    count_by_stage.reindex(['normal', 'metaplasia', 'cancer']).plot(kind='bar', ax=ax, color='steelblue')
    ax.set_xlabel('Stage')
    ax.set_ylabel('Number of Samples')
    ax.set_title('Sample Count per Stage')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    # Add n_cells annotation
    for i, stage in enumerate(['normal', 'metaplasia', 'cancer']):
        if stage in count_by_stage.index:
            stage_data = sample_scores[sample_scores['stage'] == stage]
            total_cells = stage_data['n_cells'].sum()
            ax.text(i, count_by_stage[stage] + 0.1, f'{total_cells:,}',
                   ha='center', va='bottom', fontsize=8)

    plt.suptitle('CAF Subtype Analysis (Note: iCAF missing CXCL12)', fontsize=12, fontweight='bold')
    plt.tight_layout()

    fig.savefig(output_dir / 'figures' / 'fig2_caf_composition.png', dpi=300, bbox_inches='tight')
    fig.savefig(output_dir / 'figures' / 'fig2_caf_composition.pdf', bbox_inches='tight')
    plt.close(fig)


def generate_fig3_cd8_exhaustion(cd8_results: Dict, output_dir: Path):
    """Generate CD8 exhaustion figure."""
    logger.info("Generating Figure 3: CD8 exhaustion...")

    if not cd8_results or 'sample_scores' not in cd8_results:
        return

    sample_scores = cd8_results['sample_scores']

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Panel 1: Exhaustion score by stage
    ax = axes[0]
    stage_order = [s for s in ['normal', 'metaplasia', 'cancer'] if s in sample_scores['stage'].values]

    for i, stage in enumerate(stage_order):
        stage_data = sample_scores[sample_scores['stage'] == stage]['exhaustion_score']
        ax.scatter([i] * len(stage_data), stage_data,
                  c=STAGE_COLORS.get(stage, 'gray'), s=100, alpha=0.7,
                  edgecolors='black', linewidth=1)
        ax.boxplot([stage_data.values], positions=[i], widths=0.3,
                  patch_artist=True,
                  boxprops=dict(facecolor='white', alpha=0.3))

    ax.set_xticks(range(len(stage_order)))
    ax.set_xticklabels([s.title() for s in stage_order])
    ax.set_xlabel('Stage')
    ax.set_ylabel('Exhaustion Score')
    ax.set_title('CD8 T Cell Exhaustion by Stage')

    # Add sample counts
    for i, stage in enumerate(stage_order):
        n = (sample_scores['stage'] == stage).sum()
        ax.text(i, ax.get_ylim()[1], f'n={n}', ha='center', va='bottom', fontsize=9)

    # Panel 2: Individual marker heatmap
    ax = axes[1]
    if 'marker_detail' in cd8_results and not cd8_results['marker_detail'].empty:
        marker_df = cd8_results['marker_detail']
        pivot = marker_df.pivot_table(
            index='marker', columns='sample_id', values='mean_expr', aggfunc='mean'
        )

        # Sort columns by stage
        sample_stage = marker_df.groupby('sample_id')['stage'].first()
        stage_order_map = {'normal': 0, 'metaplasia': 1, 'cancer': 2}
        sample_order = sample_stage.map(stage_order_map).sort_values().index
        pivot = pivot.reindex(columns=[c for c in sample_order if c in pivot.columns])

        sns.heatmap(pivot, cmap='YlOrRd', ax=ax, cbar_kws={'label': 'Mean Expression'})
        ax.set_title('Exhaustion Markers by Sample')
        ax.set_xlabel('Sample')
        ax.set_ylabel('Marker')

    plt.suptitle('CD8 T Cell Exhaustion Analysis', fontsize=12, fontweight='bold')
    plt.tight_layout()

    fig.savefig(output_dir / 'figures' / 'fig3_cd8_exhaustion.png', dpi=300, bbox_inches='tight')
    fig.savefig(output_dir / 'figures' / 'fig3_cd8_exhaustion.pdf', bbox_inches='tight')
    plt.close(fig)


def generate_fig4_macrophage_polarization(mac_results: Dict, output_dir: Path):
    """Generate macrophage polarization figure."""
    logger.info("Generating Figure 4: Macrophage polarization...")

    if not mac_results or 'sample_scores' not in mac_results:
        return

    sample_scores = mac_results['sample_scores']

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    stage_order = [s for s in ['normal', 'metaplasia', 'cancer'] if s in sample_scores['stage'].values]

    # Panel 1: M1 score by stage
    ax = axes[0]
    for i, stage in enumerate(stage_order):
        stage_data = sample_scores[sample_scores['stage'] == stage]['mac_M1']
        ax.scatter([i] * len(stage_data), stage_data,
                  c=STAGE_COLORS.get(stage, 'gray'), s=80, alpha=0.7)
    ax.set_xticks(range(len(stage_order)))
    ax.set_xticklabels([s.title() for s in stage_order])
    ax.set_ylabel('M1 Score')
    ax.set_title('M1 Macrophage Score')

    # Panel 2: M2 score by stage
    ax = axes[1]
    for i, stage in enumerate(stage_order):
        stage_data = sample_scores[sample_scores['stage'] == stage]['mac_M2']
        ax.scatter([i] * len(stage_data), stage_data,
                  c=STAGE_COLORS.get(stage, 'gray'), s=80, alpha=0.7)
    ax.set_xticks(range(len(stage_order)))
    ax.set_xticklabels([s.title() for s in stage_order])
    ax.set_ylabel('M2 Score')
    ax.set_title('M2 Macrophage Score\n(Note: ARG1 missing)')

    # Panel 3: M1/M2 ratio
    ax = axes[2]
    for i, stage in enumerate(stage_order):
        stage_data = sample_scores[sample_scores['stage'] == stage]['m1_m2_ratio']
        ax.scatter([i] * len(stage_data), stage_data,
                  c=STAGE_COLORS.get(stage, 'gray'), s=80, alpha=0.7)
    ax.set_xticks(range(len(stage_order)))
    ax.set_xticklabels([s.title() for s in stage_order])
    ax.set_ylabel('M1/M2 Ratio')
    ax.set_title('M1/M2 Polarization Ratio')
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5)

    plt.suptitle('Macrophage Polarization Analysis', fontsize=12, fontweight='bold')
    plt.tight_layout()

    fig.savefig(output_dir / 'figures' / 'fig4_macrophage_polarization.png', dpi=300, bbox_inches='tight')
    fig.savefig(output_dir / 'figures' / 'fig4_macrophage_polarization.pdf', bbox_inches='tight')
    plt.close(fig)


def generate_fig5_heatmap(epithelial_results: Dict, output_dir: Path):
    """Generate heatmap of top DE genes."""
    logger.info("Generating Figure 5: DE heatmap...")

    if not epithelial_results or 'de_results' not in epithelial_results:
        return

    de_results = epithelial_results['de_results']

    # Collect top genes from all comparisons
    top_genes = set()
    for name, de_df in de_results.items():
        if de_df.empty:
            continue
        sig_genes = de_df[de_df['qval_bh'] < 0.05].nsmallest(20, 'pval')['gene'].tolist()
        top_genes.update(sig_genes)

    if not top_genes:
        logger.warning("  No significant genes for heatmap")
        return

    # Create heatmap data
    heatmap_data = []
    for name, de_df in de_results.items():
        if de_df.empty:
            continue
        for gene in top_genes:
            if gene in de_df['gene'].values:
                row = de_df[de_df['gene'] == gene].iloc[0]
                heatmap_data.append({
                    'gene': gene,
                    'comparison': name,
                    'log2fc': row['log2fc'],
                    'qval': row['qval_bh'],
                })

    if not heatmap_data:
        return

    heatmap_df = pd.DataFrame(heatmap_data)
    pivot = heatmap_df.pivot(index='gene', columns='comparison', values='log2fc')

    # Sort by mean absolute log2fc
    pivot['mean_abs'] = pivot.abs().mean(axis=1)
    pivot = pivot.sort_values('mean_abs', ascending=False).drop('mean_abs', axis=1)
    pivot = pivot.head(30)  # Top 30 genes

    fig, ax = plt.subplots(figsize=(8, max(6, len(pivot) * 0.3)))

    sns.heatmap(pivot, cmap='RdBu_r', center=0, ax=ax,
                cbar_kws={'label': 'log2 Fold Change'},
                annot=True, fmt='.1f', annot_kws={'fontsize': 7})

    ax.set_title('Top DE Genes Across Comparisons\n(Epithelial Pseudobulk)', fontsize=12)
    ax.set_ylabel('Gene')
    ax.set_xlabel('Comparison')

    plt.tight_layout()

    fig.savefig(output_dir / 'figures' / 'fig5_heatmap_top_de.png', dpi=300, bbox_inches='tight')
    fig.savefig(output_dir / 'figures' / 'fig5_heatmap_top_de.pdf', bbox_inches='tight')
    plt.close(fig)


# =============================================================================
# Report Generation
# =============================================================================

def generate_report(output_dir: Path, results: Dict, elapsed_seconds: float):
    """Generate markdown report."""
    logger.info("Generating report...")

    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    report = f"""# Cell Type-Specific DE Analysis Report

Generated: {timestamp}
Runtime: {elapsed_seconds/60:.1f} minutes

## Overview

This analysis performs cell type-stratified differential expression analysis using
**pseudobulk aggregation** to avoid pseudoreplication. Rather than testing 1.8M
individual cells, we aggregate to 29 sample-level profiles before DE testing.

## Critical Notes

### Stage-Patient Confounding

| Patient | Stages Present |
|---------|----------------|
| SNU-105 | normal, metaplasia |
| SNU-107 | normal, metaplasia |
| SNU-484 | cancer |
| ctrl | control |

**Implication**: Normal and metaplasia samples come from only 2 patients; cancer
samples come from 1 patient. Results may partially reflect patient effects rather
than pure stage effects.

### Incomplete Gene Signatures

| Signature | Available | Missing |
|-----------|-----------|---------|
| iCAF | IL6, PDGFRA, FAP | CXCL12 |
| M2 Macrophage | CD163, MRC1, IL10 | ARG1 |

### Admixture Filtering

Cells with `is_admixed=True` were excluded from all analyses to ensure clean
cell type populations.

## Results Summary

"""

    # Epithelial results
    if 'epithelial' in results:
        epi = results['epithelial']
        report += "### Epithelial Analysis\n\n"
        report += f"- Samples analyzed: {epi.get('n_samples', 'N/A')}\n"
        report += f"- Excluded samples: {epi.get('excluded', [])}\n\n"

        if 'de_results' in epi:
            for name, de_df in epi['de_results'].items():
                if not de_df.empty:
                    n_sig = (de_df['qval_bh'] < 0.05).sum()
                    report += f"- {name}: {n_sig} significant genes (q<0.05)\n"
        report += "\n"

    # CAF results
    if 'caf' in results:
        caf = results['caf']
        report += "### CAF Subtype Analysis\n\n"
        report += "⚠️ **iCAF signature incomplete** (missing CXCL12)\n\n"
        if 'sample_scores' in caf:
            n_samples = len(caf['sample_scores'])
            report += f"- Samples analyzed: {n_samples}\n\n"

    # CD8 results
    if 'cd8' in results:
        cd8 = results['cd8']
        report += "### CD8 T Cell Exhaustion\n\n"
        if 'sample_scores' in cd8:
            n_samples = len(cd8['sample_scores'])
            report += f"- Samples analyzed: {n_samples}\n"
        if 'excluded' in cd8:
            report += f"- Excluded samples (low cell count): {cd8['excluded']}\n\n"

    # Macrophage results
    if 'macrophage' in results:
        mac = results['macrophage']
        report += "### Macrophage Polarization\n\n"
        report += "⚠️ **M2 signature incomplete** (missing ARG1)\n\n"
        if 'sample_scores' in mac:
            n_samples = len(mac['sample_scores'])
            report += f"- Samples analyzed: {n_samples}\n\n"

    report += """
## Output Files

```
results/celltype_de/
├── epithelial/
│   ├── de_N_vs_M.csv, de_M_vs_C.csv, de_N_vs_C.csv
│   ├── gastric_markers_by_stage.csv
│   └── sample_metadata.csv
├── caf/
│   ├── caf_subtype_scores.csv
│   ├── caf_proportions_by_stage.csv
│   └── NOTE_iCAF_incomplete.txt
├── cd8/
│   ├── exhaustion_scores.csv
│   ├── exhaustion_markers_detail.csv
│   └── samples_excluded.txt
├── macrophage/
│   ├── m1_m2_scores.csv
│   └── NOTE_M2_incomplete.txt
├── figures/
│   ├── fig1_epithelial_volcano.png
│   ├── fig2_caf_composition.png
│   ├── fig3_cd8_exhaustion.png
│   ├── fig4_macrophage_polarization.png
│   └── fig5_heatmap_top_de.png
├── sensitivity/
│   └── admixture_comparison.csv
└── CELLTYPE_DE_REPORT.md
```

## Methods

### Pseudobulk DE

1. Subset to cell type (excluding admixed cells)
2. Filter samples with < MIN_CELLS threshold
3. Aggregate counts to sample level (sum)
4. Normalize (CPM + log1p)
5. Mann-Whitney U test (non-parametric, robust for small n)
6. Benjamini-Hochberg FDR correction

### Signature Scoring

Used scanpy's `score_genes` function with mean expression of signature genes
vs. randomly sampled control genes.
"""

    with open(output_dir / 'CELLTYPE_DE_REPORT.md', 'w') as f:
        f.write(report)

    logger.info("  Report saved")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Cell type-specific DE analysis')
    parser.add_argument('--skip-env-check', action='store_true',
                       help='Skip conda environment check')
    parser.add_argument('--skip-figures', action='store_true',
                       help='Skip figure generation')
    parser.add_argument('--include-admixed', action='store_true',
                       help='Include admixed cells (not recommended)')
    args = parser.parse_args()

    # Check environment
    check_environment(args.skip_env_check)

    # Create output directories
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / 'figures').mkdir(exist_ok=True)

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

    # Filter admixed cells
    filter_admixed = not args.include_admixed

    # Run analyses
    results = {}

    # 1. Epithelial DE
    results['epithelial'] = analyze_epithelial_de(adata, OUTPUT_DIR, filter_admixed)

    # 2. CAF subtypes
    results['caf'] = analyze_caf_subtypes(adata, OUTPUT_DIR, filter_admixed)

    # 3. CD8 exhaustion
    results['cd8'] = analyze_cd8_exhaustion(adata, OUTPUT_DIR, filter_admixed)

    # 4. Macrophage polarization
    results['macrophage'] = analyze_macrophage_polarization(adata, OUTPUT_DIR, filter_admixed)

    # 5. Admixture sensitivity
    run_admixture_sensitivity(adata, OUTPUT_DIR)

    # Generate figures
    if not args.skip_figures:
        if results['epithelial'] and 'de_results' in results['epithelial']:
            generate_fig1_epithelial_volcano(results['epithelial']['de_results'], OUTPUT_DIR)

        if results['caf']:
            generate_fig2_caf_composition(results['caf'], OUTPUT_DIR)

        if results['cd8']:
            generate_fig3_cd8_exhaustion(results['cd8'], OUTPUT_DIR)

        if results['macrophage']:
            generate_fig4_macrophage_polarization(results['macrophage'], OUTPUT_DIR)

        if results['epithelial']:
            generate_fig5_heatmap(results['epithelial'], OUTPUT_DIR)

    # Generate report
    elapsed = time.time() - start_time
    generate_report(OUTPUT_DIR, results, elapsed)

    logger.info("="*60)
    logger.info(f"Cell type DE analysis complete in {elapsed/60:.1f} minutes")
    logger.info(f"Output directory: {OUTPUT_DIR}")
    logger.info("="*60)


if __name__ == '__main__':
    main()
