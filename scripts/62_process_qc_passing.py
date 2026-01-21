#!/usr/bin/env python3
"""
G4X Process QC-Passing Samples
==============================
Run full WNN + annotation pipeline on samples that passed QC.

Usage:
    conda activate enact
    python scripts/62_process_qc_passing.py [--parallel N] [--keep-admixed]

Arguments:
    --parallel      Number of parallel workers (default: 4)
    --keep-admixed  Flag admixed cells but don't remove them (for sensitivity analysis)

Output:
    results/qc_all_samples/final_processed/{sample}_final.h5ad
"""

import os
import sys
import shutil


def check_environment(skip_check: bool = False):
    """
    Verify running in correct conda environment.

    Parameters
    ----------
    skip_check : bool
        If True, skip the environment check entirely (use --skip-env-check flag)
    """
    if skip_check or os.environ.get("G4X_SKIP_ENV_CHECK"):
        return  # Allow override via flag or env var

    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    env_name = os.path.basename(conda_prefix) if conda_prefix else ""

    # Accept envs containing these patterns (e.g., 'enact-gpu', 'spatial-dev', 'scanpy3.10')
    valid_patterns = ['enact', 'spatial', 'scanpy', 'single-cell', 'scverse']
    is_valid = any(pattern in env_name.lower() for pattern in valid_patterns)

    if not is_valid:
        print("ERROR: This script requires a spatial analysis conda environment.")
        print(f"Current environment: {env_name or 'none'}")
        print(f"\nAccepted environment patterns: {', '.join(valid_patterns)}")
        print("\nOptions:")
        print("  1. Activate a valid env: conda activate enact")
        print("  2. Skip check: --skip-env-check")
        print("  3. Set env var: export G4X_SKIP_ENV_CHECK=1")
        sys.exit(1)


def check_disk_space(path, required_gb=5.0):
    """
    Check available disk space and warn/abort if insufficient.

    This implements a 3-tier disk space monitoring system:
    - Tier 1 (Initial): 5 GB required at pipeline start (comfortable buffer)
    - Tier 2 (Periodic): 2 GB required during processing (checked every 5 samples)
    - Tier 3 (Critical): 1 GB minimum (hard abort to prevent data corruption)

    The rationale:
    - Each processed h5ad is ~10-50 MB depending on cell count
    - Processing 32 samples could generate ~1-2 GB of output
    - 5 GB initial ensures room for intermediate files and logs
    - 2 GB periodic catches gradual exhaustion
    - 1 GB critical prevents mid-write failures

    Parameters
    ----------
    path : str or Path
        Directory to check
    required_gb : float
        Minimum required space in GB (default: 5.0 for initial check)

    Returns
    -------
    float
        Available space in GB
    """
    total, used, free = shutil.disk_usage(path)
    free_gb = free / (1024 ** 3)
    if free_gb < required_gb:
        print(f"WARNING: Low disk space! {free_gb:.1f} GB available, {required_gb:.1f} GB recommended")
        # Tier 3: Critical threshold - abort to prevent data loss
        if free_gb < 1.0:
            print("ERROR: Critically low disk space (<1 GB). Aborting to prevent data loss.")
            sys.exit(1)
    return free_gb


# Defer environment check until after argparse (see main())

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from scipy import sparse
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import argparse
import logging
import warnings
import gc

warnings.filterwarnings('ignore')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# =============================================================================
# Configuration
# =============================================================================

INPUT_DIR = Path("/home/user/g4x-choi-batch2-analysis/results/qc_all_samples/raw")
QC_SUMMARY = Path("/home/user/g4x-choi-batch2-analysis/results/qc_all_samples/sample_qc_summary.csv")
OUTPUT_DIR = Path("/home/user/g4x-choi-batch2-analysis/results/qc_all_samples/final_processed")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Marker definitions for annotation
LINEAGE_MARKERS = {
    'epithelial': ['PanCK'],
    'immune': ['CD45'],
    'stromal': ['aSMA'],
    'endothelial': ['CD31'],
}

CELL_TYPE_MARKERS = {
    # Immune subtypes
    'T_cell': ['CD3'],
    'CD4_T': ['CD3', 'CD4'],
    'CD8_T': ['CD3', 'CD8'],
    'Treg': ['CD3', 'CD4', 'FOXP3'],
    'B_cell': ['CD20'],
    'Macrophage': ['CD68'],
    'DC': ['CD11c', 'HLA-DR'],
    # Functional markers
    'Proliferating': ['KI67'],
    'PD1_positive': ['PD1'],
    'PDL1_positive': ['PDL1'],
}

# Conflict pairs for admixture detection
CONFLICT_PAIRS = [
    ('PanCK', 'CD45'),      # Epithelial + Immune
    ('PanCK', 'aSMA'),      # Epithelial + Stromal
    ('CD45', 'aSMA'),       # Immune + Stromal
    ('CD31', 'CD45'),       # Endothelial + Immune
    ('CD3', 'CD68'),        # T cell + Macrophage
    ('CD4', 'CD8'),         # CD4 + CD8
]

ADMIX_THRESHOLD = 0.7  # Percentile for "high" expression
ADMIX_CUTOFF = 0.3     # Flag cells above this score

# Cell-level QC thresholds (applied before integration)
# Note: filter_admixed can be overridden by --keep-admixed CLI flag
CELL_QC_THRESHOLDS = {
    'min_counts': 10,       # Minimum transcripts per cell
    'min_genes': 5,         # Minimum genes detected per cell
    'filter_admixed': True, # Remove admixed cells (set False via --keep-admixed for sensitivity)
}


# =============================================================================
# Processing Functions
# =============================================================================

def apply_cell_level_qc(adata, thresholds: dict = CELL_QC_THRESHOLDS) -> ad.AnnData:
    """
    Apply cell-level QC filtering BEFORE integration.

    Removes:
    - Low count cells (< min_counts)
    - Low gene cells (< min_genes)
    - Admixed cells (if filter_admixed=True)

    Returns filtered AnnData with QC metrics in uns.
    """
    n_start = adata.n_obs

    # Compute QC metrics if not present
    if 'n_counts' not in adata.obs:
        adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
    if 'n_genes' not in adata.obs:
        adata.obs['n_genes'] = np.array((adata.X > 0).sum(axis=1)).flatten()

    # Build filter mask
    keep_mask = np.ones(adata.n_obs, dtype=bool)

    # Filter low counts
    low_counts = adata.obs['n_counts'] < thresholds['min_counts']
    n_low_counts = low_counts.sum()
    keep_mask &= ~low_counts

    # Filter low genes
    low_genes = adata.obs['n_genes'] < thresholds['min_genes']
    n_low_genes = low_genes.sum()
    keep_mask &= ~low_genes

    # Filter admixed (if computed and enabled)
    n_admixed = 0
    if thresholds.get('filter_admixed', False) and 'is_admixed' in adata.obs:
        admixed = adata.obs['is_admixed'].astype(bool)
        n_admixed = admixed.sum()
        keep_mask &= ~admixed

    # Apply filter
    adata_filtered = adata[keep_mask].copy()

    # Store QC info
    n_removed = n_start - adata_filtered.n_obs
    adata_filtered.uns['cell_qc'] = {
        'n_start': n_start,
        'n_filtered': adata_filtered.n_obs,
        'n_removed_total': n_removed,
        'n_low_counts': int(n_low_counts),
        'n_low_genes': int(n_low_genes),
        'n_admixed': int(n_admixed),
        'pct_removed': (n_removed / n_start) * 100 if n_start > 0 else 0,
    }

    logger.debug(f"  Cell QC: {n_start:,} -> {adata_filtered.n_obs:,} ({n_removed:,} removed)")

    return adata_filtered


def compute_wnn_integration(adata, n_pcs_rna=30, n_pcs_prot=15, n_neighbors=20):
    """
    Compute simplified Weighted Nearest Neighbors integration of RNA and Protein.

    IMPORTANT: This is a simplified WNN implementation using variance-weighted
    concatenation of PCA embeddings. It differs from true WNN (Seurat v4/muon)
    which computes:
    1. Modality-specific k-NN graphs
    2. Per-cell weights based on local neighborhood consistency
    3. Cell-specific modality importance

    Our simplified approach:
    - Global variance-derived weights (not per-cell)
    - Direct weighted concatenation (not neighbor graph fusion)
    - RNA typically dominates due to higher dimensionality

    For production use, consider:
    - muon.pp.neighbors(mdata, method='wnn') for true WNN
    - MOFA+ for multi-modal factor analysis
    - totalVI for joint probabilistic modeling

    Returns AnnData with X_wnn embedding in obsm.
    """
    logger.debug("Computing WNN integration...")

    # RNA PCA
    if 'X_pca' not in adata.obsm:
        sc.pp.pca(adata, n_comps=min(n_pcs_rna, adata.n_vars - 1))

    rna_pca = adata.obsm['X_pca'][:, :min(n_pcs_rna, adata.obsm['X_pca'].shape[1])].copy()

    # Protein PCA
    protein_data = adata.obsm['protein'].copy()
    if sparse.issparse(protein_data):
        protein_data = protein_data.toarray()

    # Scale protein data
    scaler = StandardScaler()
    protein_scaled = scaler.fit_transform(protein_data)

    # Protein PCA
    n_prot_comps = min(n_pcs_prot, protein_data.shape[1] - 1)
    pca_prot = PCA(n_components=n_prot_comps)
    protein_pca = pca_prot.fit_transform(protein_scaled)

    adata.obsm['X_pca_protein'] = protein_pca

    # Compute modality weights based on within-modality variance
    rna_var = np.var(rna_pca, axis=0).sum()
    prot_var = np.var(protein_pca, axis=0).sum()
    total_var = rna_var + prot_var

    rna_weight = rna_var / total_var
    prot_weight = prot_var / total_var

    logger.debug(f"  Modality weights: RNA={rna_weight:.3f}, Protein={prot_weight:.3f}")

    # Scale PCAs to same range
    rna_scaled = StandardScaler().fit_transform(rna_pca)
    prot_scaled = StandardScaler().fit_transform(protein_pca)

    # Weighted combination
    X_wnn = np.hstack([
        rna_scaled * np.sqrt(rna_weight),
        prot_scaled * np.sqrt(prot_weight)
    ])

    adata.obsm['X_wnn'] = X_wnn
    adata.uns['wnn_weights'] = {'rna': rna_weight, 'protein': prot_weight}

    return adata


def compute_admixture_score(adata):
    """
    Compute admixture score based on conflicting marker co-expression.

    High admixture = potential segmentation artifact (doublet/multiplet).
    """
    logger.debug("Computing admixture scores...")

    protein_names = adata.uns.get('protein_names', [])
    protein_data = adata.obsm['protein']

    if sparse.issparse(protein_data):
        protein_data = protein_data.toarray()

    # Create protein DataFrame
    protein_df = pd.DataFrame(
        protein_data,
        columns=protein_names,
        index=adata.obs_names
    )

    # Compute percentile thresholds
    thresholds = protein_df.quantile(ADMIX_THRESHOLD)

    # Binary high expression
    high_expr = protein_df > thresholds

    # Count conflicts
    conflict_count = np.zeros(len(adata))

    for m1, m2 in CONFLICT_PAIRS:
        if m1 in high_expr.columns and m2 in high_expr.columns:
            conflict = high_expr[m1] & high_expr[m2]
            conflict_count += conflict.values.astype(int)

    # Normalize to [0, 1]
    max_conflicts = len(CONFLICT_PAIRS)
    admix_score = conflict_count / max_conflicts

    adata.obs['admixture_score'] = admix_score
    adata.obs['is_admixed'] = admix_score > ADMIX_CUTOFF

    n_admixed = adata.obs['is_admixed'].sum()
    pct_admixed = n_admixed / len(adata) * 100
    logger.debug(f"  Admixed cells: {n_admixed:,} ({pct_admixed:.1f}%)")

    return adata


def annotate_lineage(adata):
    """Assign major lineage based on protein markers."""
    logger.debug("Annotating lineages...")

    protein_names = adata.uns.get('protein_names', [])
    protein_data = adata.obsm['protein']

    if sparse.issparse(protein_data):
        protein_data = protein_data.toarray()

    protein_df = pd.DataFrame(
        protein_data,
        columns=protein_names,
        index=adata.obs_names
    )

    # Score each lineage
    lineage_scores = {}
    for lineage, markers in LINEAGE_MARKERS.items():
        available = [m for m in markers if m in protein_df.columns]
        if available:
            # Z-score then mean
            scores = protein_df[available].apply(lambda x: (x - x.mean()) / (x.std() + 1e-8))
            lineage_scores[lineage] = scores.mean(axis=1)

    # Assign lineage
    if lineage_scores:
        score_df = pd.DataFrame(lineage_scores)
        adata.obs['lineage'] = score_df.idxmax(axis=1)
        adata.obs['lineage_score'] = score_df.max(axis=1)
    else:
        adata.obs['lineage'] = 'unknown'
        adata.obs['lineage_score'] = 0.0

    return adata


def annotate_cell_types(adata):
    """Assign cell types within lineages."""
    logger.debug("Annotating cell types...")

    protein_names = adata.uns.get('protein_names', [])
    protein_data = adata.obsm['protein']

    if sparse.issparse(protein_data):
        protein_data = protein_data.toarray()

    protein_df = pd.DataFrame(
        protein_data,
        columns=protein_names,
        index=adata.obs_names
    )

    # Initialize cell type column
    adata.obs['cell_type'] = adata.obs['lineage'].copy()

    # Refine immune cells
    immune_mask = adata.obs['lineage'] == 'immune'

    if immune_mask.any():
        immune_idx = adata.obs_names[immune_mask]

        for cell_type, markers in CELL_TYPE_MARKERS.items():
            available = [m for m in markers if m in protein_df.columns]
            if available:
                # All markers must be above median for assignment
                all_positive = protein_df.loc[immune_idx, available].apply(
                    lambda x: x > x.median(), axis=0
                ).all(axis=1)

                positive_idx = immune_idx[all_positive]
                adata.obs.loc[positive_idx, 'cell_type'] = cell_type

    return adata


def process_sample(sample_info: dict, keep_admixed: bool = False) -> dict:
    """
    Process a single sample through the full pipeline.

    Parameters
    ----------
    sample_info : dict
        Contains 'sample_id' and 'path' keys
    keep_admixed : bool
        If True, flag admixed cells but don't remove them (for sensitivity analysis)

    Returns
    -------
    dict
        Processing result with status and metrics
    """
    sample_id = sample_info['sample_id']
    input_path = Path(sample_info['path'])
    output_path = OUTPUT_DIR / f"{sample_id}_final.h5ad"

    result = {'sample_id': sample_id, 'status': 'SUCCESS'}

    try:
        # Load
        adata = sc.read_h5ad(input_path)
        n_cells_raw = adata.n_obs

        # =====================================================================
        # CELL-LEVEL QC (BEFORE integration)
        # =====================================================================

        # Step 1: Compute admixture scores on raw data (needed for QC filtering)
        adata = compute_admixture_score(adata)
        n_admixed_raw = int(adata.obs['is_admixed'].sum())

        # Step 2: Apply cell-level QC filtering (removes low-count, low-gene, admixed)
        # If keep_admixed=True, cells are flagged but not removed (for sensitivity analysis)
        qc_thresholds = CELL_QC_THRESHOLDS.copy()
        if keep_admixed:
            qc_thresholds['filter_admixed'] = False
        adata = apply_cell_level_qc(adata, qc_thresholds)
        n_cells_post_qc = adata.n_obs

        logger.info(f"  {sample_id}: {n_cells_raw:,} -> {n_cells_post_qc:,} cells after QC")

        # =====================================================================
        # PREPROCESSING
        # =====================================================================

        # Normalize RNA
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Variable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=200, subset=False)

        # PCA
        sc.pp.pca(adata, n_comps=min(30, adata.n_vars - 1))

        # =====================================================================
        # INTEGRATION & CLUSTERING
        # =====================================================================

        # WNN integration (simplified weighted concatenation - see documentation)
        adata = compute_wnn_integration(adata)

        # Annotation
        adata = annotate_lineage(adata)
        adata = annotate_cell_types(adata)

        # Clustering on WNN
        sc.pp.neighbors(adata, use_rep='X_wnn', n_neighbors=15)
        sc.tl.leiden(adata, resolution=0.5)
        sc.tl.umap(adata)

        # Save
        adata.write(output_path)

        # Metrics
        result['n_cells_raw'] = n_cells_raw
        result['n_cells'] = adata.n_obs
        result['n_cells_removed'] = n_cells_raw - adata.n_obs
        result['pct_removed'] = (n_cells_raw - adata.n_obs) / n_cells_raw * 100
        result['n_admixed_raw'] = n_admixed_raw
        result['n_clusters'] = adata.obs['leiden'].nunique()

        # Lineage distribution
        lineage_counts = adata.obs['lineage'].value_counts()
        for lineage, count in lineage_counts.items():
            result[f'n_{lineage}'] = int(count)

        del adata
        gc.collect()

    except Exception as e:
        result['status'] = 'FAILED'
        result['error'] = str(e)
        logger.error(f"Failed to process {sample_id}: {e}")

    return result


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Process QC-passing G4X samples')
    parser.add_argument('--parallel', type=int, default=4, help='Number of parallel workers')
    parser.add_argument('--keep-admixed', action='store_true',
                        help='Flag admixed cells but keep them (for sensitivity analysis)')
    parser.add_argument('--skip-env-check', action='store_true',
                        help='Skip conda environment validation (use if running from equivalent env)')
    args = parser.parse_args()

    # Environment check (deferred to allow --skip-env-check flag)
    check_environment(skip_check=args.skip_env_check)

    logger.info("=" * 60)
    logger.info("G4X Process QC-Passing Samples")
    logger.info("=" * 60)

    if args.keep_admixed:
        logger.info("NOTE: --keep-admixed enabled - admixed cells flagged but NOT removed")
        logger.info("Use this for sensitivity analysis; default behavior filters admixed cells")

    # Load QC summary
    if not QC_SUMMARY.exists():
        logger.error(f"QC summary not found: {QC_SUMMARY}")
        logger.error("Run 61_comprehensive_qc.py first")
        return

    qc_df = pd.read_csv(QC_SUMMARY)
    passing = qc_df[qc_df['qc_status'] != 'FAIL']

    logger.info(f"QC-passing samples: {len(passing)}/{len(qc_df)}")

    # Load manifest for paths
    manifest = pd.read_csv(INPUT_DIR / "loading_manifest.csv")

    # Prepare processing list
    samples_to_process = []
    for _, row in passing.iterrows():
        sample_id = row['sample_id']
        manifest_row = manifest[manifest['sample_id'] == sample_id]
        if len(manifest_row) > 0:
            samples_to_process.append({
                'sample_id': sample_id,
                'path': manifest_row['path'].values[0]
            })

    logger.info(f"Processing {len(samples_to_process)} samples...")

    # Initial disk space check
    free_gb = check_disk_space(OUTPUT_DIR, required_gb=5.0)
    logger.info(f"Disk space available: {free_gb:.1f} GB")

    # Process samples
    results = []
    disk_check_interval = 5  # Check disk every N samples

    if args.parallel > 1:
        # Parallel processing
        with ProcessPoolExecutor(max_workers=args.parallel) as executor:
            futures = {
                executor.submit(process_sample, s, args.keep_admixed): s['sample_id']
                for s in samples_to_process
            }

            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing"):
                sample_id = futures[future]
                try:
                    result = future.result()
                    results.append(result)
                    logger.info(f"  {sample_id}: {result['status']}")
                except Exception as e:
                    results.append({'sample_id': sample_id, 'status': 'FAILED', 'error': str(e)})
                    logger.error(f"  {sample_id}: FAILED - {e}")
    else:
        # Sequential processing
        for i, sample_info in enumerate(tqdm(samples_to_process, desc="Processing")):
            # Periodic disk space check
            if i > 0 and i % disk_check_interval == 0:
                free_gb = check_disk_space(OUTPUT_DIR, required_gb=2.0)
                logger.debug(f"Disk check: {free_gb:.1f} GB free")

            result = process_sample(sample_info, args.keep_admixed)
            results.append(result)
            logger.info(f"  {result['sample_id']}: {result['status']}")

    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(OUTPUT_DIR / "processing_results.csv", index=False)

    # Summary
    n_success = (results_df['status'] == 'SUCCESS').sum()
    n_failed = (results_df['status'] == 'FAILED').sum()

    logger.info("=" * 60)
    logger.info("Processing Complete")
    logger.info(f"  Success: {n_success}")
    logger.info(f"  Failed: {n_failed}")

    if n_success > 0:
        total_cells = results_df[results_df['status'] == 'SUCCESS']['n_cells'].sum()
        logger.info(f"  Total cells processed: {total_cells:,}")

    logger.info(f"  Output: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
