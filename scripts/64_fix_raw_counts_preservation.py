#!/usr/bin/env python3
"""
64_fix_raw_counts_preservation.py

CRITICAL FIX: Re-process samples from raw files, preserving raw counts in layers.

This script fixes the bug in script 62 where raw counts were lost during normalization.
It reads from results/qc_all_samples/raw/ and outputs to results/qc_all_samples/final_processed_fixed/

The key fix:
    BEFORE normalization: adata.layers['counts'] = adata.X.copy()
    THEN normalize: sc.pp.normalize_total(); sc.pp.log1p()
    SAVE: Both X (normalized) and layers['counts'] (raw) are preserved

Usage:
    python scripts/64_fix_raw_counts_preservation.py [--parallel N] [--sample SAMPLE_ID]

Author: G4X Analysis Pipeline
Date: 2026-01-22
"""

import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from concurrent.futures import ProcessPoolExecutor, as_completed

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(f'logs/64_fix_raw_counts_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    ]
)
logger = logging.getLogger(__name__)

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
RAW_DIR = PROJECT_ROOT / "results" / "qc_all_samples" / "raw"
OUTPUT_DIR = PROJECT_ROOT / "results" / "qc_all_samples" / "final_processed_fixed"
QC_SUMMARY = PROJECT_ROOT / "results" / "qc_all_samples" / "sample_qc_summary.csv"


def validate_raw_counts(adata: ad.AnnData, sample_id: str) -> bool:
    """
    Validate that adata.X contains raw counts (integers).

    Args:
        adata: AnnData object
        sample_id: Sample identifier for logging

    Returns:
        True if valid raw counts, False otherwise
    """
    X = adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()

    # Check 1: Should be non-negative
    if (X < 0).any():
        logger.error(f"{sample_id}: Negative values found - not raw counts")
        return False

    # Check 2: Should be integer-like
    is_integer_like = np.allclose(X, np.round(X))
    if not is_integer_like:
        logger.error(f"{sample_id}: Non-integer values found - not raw counts")
        return False

    # Check 3: Max should be reasonable for counts (not log-transformed)
    x_max = X.max()
    if x_max < 10:
        logger.warning(f"{sample_id}: Max value {x_max:.1f} is very low - verify this is correct")

    logger.info(f"{sample_id}: Raw counts validated (max={x_max:.0f}, integer-like={is_integer_like})")
    return True


def process_sample_with_counts_preservation(raw_file: Path, output_dir: Path) -> dict:
    """
    Process a single sample, preserving raw counts in layers.

    This is the FIXED version of the processing in script 62.

    Args:
        raw_file: Path to raw h5ad file
        output_dir: Directory for output

    Returns:
        Dict with processing results
    """
    sample_id = raw_file.stem.replace('_raw', '')
    result = {'sample_id': sample_id, 'status': 'unknown', 'error': None}

    try:
        logger.info(f"Processing {sample_id}...")

        # Load raw data
        adata = sc.read_h5ad(raw_file)
        n_cells_raw = adata.n_obs

        # CRITICAL: Validate raw counts BEFORE any processing
        if not validate_raw_counts(adata, sample_id):
            result['status'] = 'failed'
            result['error'] = 'Input is not raw counts'
            return result

        # =====================================================================
        # CRITICAL FIX: Store raw counts BEFORE any normalization
        # =====================================================================
        logger.info(f"{sample_id}: Storing raw counts in .layers['counts']")
        adata.layers['counts'] = adata.X.copy()

        # Also store as integer type for explicit clarity
        if hasattr(adata.layers['counts'], 'toarray'):
            # Sparse matrix - keep sparse but ensure we have a copy
            pass
        else:
            adata.layers['counts'] = adata.layers['counts'].astype(np.int32)

        # =====================================================================
        # Standard QC filtering (same as script 62)
        # =====================================================================

        # Basic QC metrics (should already exist but recalculate for safety)
        sc.pp.calculate_qc_metrics(adata, inplace=True)

        # Filter cells with too few counts or genes
        min_counts = 10
        min_genes = 5
        sc.pp.filter_cells(adata, min_counts=min_counts)
        sc.pp.filter_cells(adata, min_genes=min_genes)

        n_cells_filtered = adata.n_obs
        logger.info(f"{sample_id}: {n_cells_raw:,} -> {n_cells_filtered:,} cells after filtering")

        # =====================================================================
        # Normalization (AFTER storing raw counts)
        # =====================================================================
        logger.info(f"{sample_id}: Normalizing...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # =====================================================================
        # Verify raw counts are preserved
        # =====================================================================
        counts_layer = adata.layers['counts']
        if hasattr(counts_layer, 'toarray'):
            counts_sample = counts_layer[:100, :50].toarray()
        else:
            counts_sample = counts_layer[:100, :50]

        counts_max = counts_sample.max()
        counts_integer = np.allclose(counts_sample, np.round(counts_sample))

        X_sample = adata.X[:100, :50]
        if hasattr(X_sample, 'toarray'):
            X_sample = X_sample.toarray()
        x_max = X_sample.max()

        logger.info(f"{sample_id}: X (normalized) max={x_max:.2f}, layers['counts'] max={counts_max:.0f}")

        if not counts_integer:
            result['status'] = 'failed'
            result['error'] = 'Raw counts were corrupted during processing'
            return result

        if counts_max < 10:
            logger.warning(f"{sample_id}: Raw counts max is low ({counts_max:.0f}) - verify input")

        # =====================================================================
        # Save with preserved counts
        # =====================================================================
        output_file = output_dir / f"{sample_id}_final_fixed.h5ad"
        adata.write(output_file)

        result['status'] = 'success'
        result['n_cells'] = n_cells_filtered
        result['counts_max'] = float(counts_max)
        result['normalized_max'] = float(x_max)
        result['output_file'] = str(output_file)

        logger.info(f"{sample_id}: Saved to {output_file}")

    except Exception as e:
        logger.error(f"{sample_id}: Error - {e}")
        result['status'] = 'failed'
        result['error'] = str(e)

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Fix raw counts preservation in processed samples"
    )
    parser.add_argument('--parallel', type=int, default=4,
                        help='Number of parallel workers (default: 4)')
    parser.add_argument('--sample', type=str, default=None,
                        help='Process only this sample (for testing)')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite existing output files')
    args = parser.parse_args()

    logger.info("=" * 60)
    logger.info("SCRIPT 64: Fix Raw Counts Preservation")
    logger.info("=" * 60)

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Get list of raw files
    raw_files = sorted(RAW_DIR.glob("*_raw.h5ad"))
    logger.info(f"Found {len(raw_files)} raw files in {RAW_DIR}")

    # Load QC summary to identify passing samples
    if QC_SUMMARY.exists():
        qc_df = pd.read_csv(QC_SUMMARY)
        passing_samples = set(qc_df[qc_df['qc_status'].isin(['PASS', 'WARN'])]['sample_id'])
        logger.info(f"Passing samples from QC: {len(passing_samples)}")
    else:
        logger.warning("QC summary not found - processing all samples")
        passing_samples = None

    # Filter to passing samples only
    if passing_samples:
        raw_files = [f for f in raw_files if f.stem.replace('_raw', '') in passing_samples]
        logger.info(f"Processing {len(raw_files)} passing samples")

    # Filter to single sample if specified
    if args.sample:
        raw_files = [f for f in raw_files if args.sample in f.stem]
        if not raw_files:
            logger.error(f"Sample {args.sample} not found")
            sys.exit(1)
        logger.info(f"Processing single sample: {args.sample}")

    # Skip already processed files unless --force
    if not args.force:
        existing = set(f.stem.replace('_final_fixed', '') for f in OUTPUT_DIR.glob("*_final_fixed.h5ad"))
        raw_files = [f for f in raw_files if f.stem.replace('_raw', '') not in existing]
        if existing:
            logger.info(f"Skipping {len(existing)} already processed samples (use --force to reprocess)")

    if not raw_files:
        logger.info("No samples to process")
        return

    logger.info(f"Will process {len(raw_files)} samples with {args.parallel} workers")

    # Process samples
    results = []

    if args.parallel == 1 or len(raw_files) == 1:
        # Sequential processing
        for raw_file in raw_files:
            result = process_sample_with_counts_preservation(raw_file, OUTPUT_DIR)
            results.append(result)
    else:
        # Parallel processing
        with ProcessPoolExecutor(max_workers=args.parallel) as executor:
            futures = {
                executor.submit(process_sample_with_counts_preservation, f, OUTPUT_DIR): f
                for f in raw_files
            }

            for future in as_completed(futures):
                result = future.result()
                results.append(result)

    # Summary
    logger.info("=" * 60)
    logger.info("PROCESSING SUMMARY")
    logger.info("=" * 60)

    success = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'failed']

    logger.info(f"Success: {len(success)}")
    logger.info(f"Failed: {len(failed)}")

    if failed:
        logger.error("Failed samples:")
        for r in failed:
            logger.error(f"  {r['sample_id']}: {r['error']}")

    # Save results
    results_df = pd.DataFrame(results)
    results_file = OUTPUT_DIR / "processing_results.csv"
    results_df.to_csv(results_file, index=False)
    logger.info(f"Results saved to {results_file}")

    logger.info("=" * 60)
    logger.info("NEXT STEP: Run 65_rebuild_merged_from_raw.py")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
