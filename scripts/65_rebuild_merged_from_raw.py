#!/usr/bin/env python3
"""
65_rebuild_merged_from_raw.py

Rebuild merged files with PROPER raw counts preservation.

This script:
1. Loads fixed processed files (from script 64) OR raw files directly
2. Merges all passing samples
3. Saves TRUE raw counts in merged_counts.h5ad
4. Normalizes and saves merged_normalized.h5ad
5. Applies batch correction and saves merged_corrected.h5ad
6. Validates raw counts are properly preserved at each step

The KEY DIFFERENCE from the original script 63:
- Raw counts are preserved in .layers['counts'] throughout
- merged_counts.h5ad actually contains raw counts (not log-normalized)
- Explicit validation at each step

Usage:
    python scripts/65_rebuild_merged_from_raw.py [--method harmony|scvi|none]

Prerequisites:
    - Either run script 64 first (recommended), OR
    - Use --from-raw flag to process directly from raw files

Author: G4X Analysis Pipeline
Date: 2026-01-22
"""

import argparse
import logging
import sys
import os
from pathlib import Path
from datetime import datetime
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Setup logging
LOG_FILE = f'logs/65_rebuild_merged_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
os.makedirs('logs', exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(LOG_FILE)
    ]
)
logger = logging.getLogger(__name__)

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
RAW_DIR = PROJECT_ROOT / "results" / "qc_all_samples" / "raw"
FIXED_DIR = PROJECT_ROOT / "results" / "qc_all_samples" / "final_processed_fixed"
OUTPUT_DIR = PROJECT_ROOT / "results" / "qc_all_samples" / "merged_fixed"
QC_SUMMARY = PROJECT_ROOT / "results" / "qc_all_samples" / "sample_qc_summary.csv"


def validate_counts(adata: ad.AnnData, location: str, expected_type: str) -> bool:
    """
    Validate counts at a specific location in the AnnData object.

    Args:
        adata: AnnData object
        location: 'X' or 'layers[counts]'
        expected_type: 'raw' (integers) or 'normalized' (log-transformed)

    Returns:
        True if validation passes
    """
    if location == 'X':
        data = adata.X
    elif location == "layers['counts']":
        if 'counts' not in adata.layers:
            logger.error(f"Validation failed: layers['counts'] not found")
            return False
        data = adata.layers['counts']
    else:
        raise ValueError(f"Unknown location: {location}")

    # Sample data for efficiency
    if hasattr(data, 'toarray'):
        sample = data[:min(1000, data.shape[0]), :min(100, data.shape[1])].toarray()
    else:
        sample = data[:min(1000, data.shape[0]), :min(100, data.shape[1])]

    x_max = sample.max()
    x_min = sample.min()
    x_mean = sample.mean()
    is_integer = np.allclose(sample, np.round(sample))

    if expected_type == 'raw':
        # Raw counts should be non-negative integers with max > 10
        if x_min < 0:
            logger.error(f"Validation failed at {location}: negative values found (min={x_min:.2f})")
            return False
        if not is_integer:
            logger.error(f"Validation failed at {location}: non-integer values found")
            return False
        if x_max < 10:
            logger.warning(f"Validation warning at {location}: max value is low ({x_max:.1f})")
        logger.info(f"Validation PASSED at {location}: raw counts (max={x_max:.0f}, mean={x_mean:.2f})")
        return True

    elif expected_type == 'normalized':
        # Normalized data should have max < 15 (log1p range) and not be integers
        if x_max > 20:
            logger.warning(f"Validation warning at {location}: max={x_max:.2f} seems high for normalized")
        if is_integer and x_max > 10:
            logger.error(f"Validation failed at {location}: appears to be raw counts, not normalized")
            return False
        logger.info(f"Validation PASSED at {location}: normalized (max={x_max:.2f}, mean={x_mean:.4f})")
        return True

    return True


def load_and_merge_samples(input_dir: Path, from_raw: bool = False) -> ad.AnnData:
    """
    Load and merge samples, ensuring raw counts are preserved.

    Args:
        input_dir: Directory containing h5ad files
        from_raw: If True, apply QC filtering to raw files

    Returns:
        Merged AnnData with raw counts in .layers['counts']
    """
    # Get file list
    if from_raw:
        files = sorted(input_dir.glob("*_raw.h5ad"))
    else:
        files = sorted(input_dir.glob("*_final_fixed.h5ad"))

    if not files:
        raise ValueError(f"No h5ad files found in {input_dir}")

    logger.info(f"Found {len(files)} files to merge")

    # Load QC summary for filtering
    if QC_SUMMARY.exists():
        qc_df = pd.read_csv(QC_SUMMARY)
        passing = set(qc_df[qc_df['qc_status'].isin(['PASS', 'WARN'])]['sample_id'])
        logger.info(f"Filtering to {len(passing)} passing samples from QC")
    else:
        passing = None
        logger.warning("QC summary not found - loading all samples")

    # Load samples
    adatas = []
    sample_ids = []

    for f in files:
        # Extract sample ID
        sample_id = f.stem.replace('_raw', '').replace('_final_fixed', '')

        # Skip failed samples
        if passing and sample_id not in passing:
            logger.info(f"Skipping {sample_id} (not in passing list)")
            continue

        logger.info(f"Loading {sample_id}...")
        adata = sc.read_h5ad(f)

        # If loading from raw, need to process
        if from_raw:
            # Validate input
            if not validate_counts(adata, 'X', 'raw'):
                logger.error(f"Skipping {sample_id} - invalid raw counts")
                continue

            # Store raw counts BEFORE any processing
            adata.layers['counts'] = adata.X.copy()

            # Basic QC filtering
            sc.pp.filter_cells(adata, min_counts=10)
            sc.pp.filter_cells(adata, min_genes=5)

            # Normalize
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)

        # Verify counts layer exists
        if 'counts' not in adata.layers:
            logger.error(f"Skipping {sample_id} - no counts layer found")
            continue

        # Add sample ID to obs
        adata.obs['sample_id'] = sample_id

        adatas.append(adata)
        sample_ids.append(sample_id)
        logger.info(f"  {sample_id}: {adata.n_obs:,} cells")

    if not adatas:
        raise ValueError("No samples loaded successfully")

    logger.info(f"Merging {len(adatas)} samples...")

    # Concatenate - use 'outer' join to keep all genes
    adata_merged = ad.concat(adatas, join='outer', label='sample_id', keys=sample_ids)
    adata_merged.obs_names_make_unique()

    # CRITICAL: Verify counts layer survived concatenation
    if 'counts' not in adata_merged.layers:
        logger.error("CRITICAL: counts layer lost during concatenation!")
        # Attempt to reconstruct
        logger.info("Attempting to reconstruct counts layer from individual samples...")
        # This shouldn't happen with anndata concat, but defensive coding
        raise RuntimeError("counts layer lost - cannot proceed")

    logger.info(f"Merged: {adata_merged.n_obs:,} cells, {adata_merged.n_vars} genes")

    return adata_merged


def apply_harmony_correction(adata: ad.AnnData, batch_key: str = 'lane') -> ad.AnnData:
    """Apply Harmony batch correction."""
    try:
        import scanpy.external as sce

        logger.info(f"Applying Harmony batch correction on {batch_key}...")

        # Run Harmony
        sce.pp.harmony_integrate(adata, key=batch_key)

        # Recompute neighbors and UMAP from corrected embedding
        sc.pp.neighbors(adata, use_rep='X_pca_harmony')
        sc.tl.umap(adata)

        logger.info("Harmony correction complete")
        return adata

    except Exception as e:
        logger.error(f"Harmony failed: {e}")
        return adata


def main():
    parser = argparse.ArgumentParser(
        description="Rebuild merged files with proper raw counts preservation"
    )
    parser.add_argument('--method', choices=['harmony', 'scvi', 'none'], default='harmony',
                        help='Batch correction method (default: harmony)')
    parser.add_argument('--batch-key', default='lane',
                        help='Batch key for correction (default: lane)')
    parser.add_argument('--from-raw', action='store_true',
                        help='Load directly from raw files (applies QC filtering)')
    parser.add_argument('--skip-correction', action='store_true',
                        help='Skip batch correction step')
    args = parser.parse_args()

    logger.info("=" * 70)
    logger.info("SCRIPT 65: Rebuild Merged Files with Raw Counts Preservation")
    logger.info("=" * 70)
    logger.info(f"Log file: {LOG_FILE}")

    # Determine input directory
    if args.from_raw:
        input_dir = RAW_DIR
        logger.info(f"Loading from RAW files: {input_dir}")
    else:
        input_dir = FIXED_DIR
        if not FIXED_DIR.exists() or not any(FIXED_DIR.glob("*.h5ad")):
            logger.error(f"Fixed directory not found or empty: {FIXED_DIR}")
            logger.error("Run script 64 first, or use --from-raw flag")
            sys.exit(1)
        logger.info(f"Loading from FIXED files: {input_dir}")

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # Step 1: Load and merge samples
    # =========================================================================
    logger.info("-" * 70)
    logger.info("STEP 1: Load and merge samples")
    logger.info("-" * 70)

    adata = load_and_merge_samples(input_dir, from_raw=args.from_raw)

    # Validate raw counts in layers
    if not validate_counts(adata, "layers['counts']", 'raw'):
        logger.error("CRITICAL: Raw counts validation failed after merge")
        sys.exit(1)

    # =========================================================================
    # Step 2: Save TRUE raw counts file
    # =========================================================================
    logger.info("-" * 70)
    logger.info("STEP 2: Save merged_counts.h5ad (TRUE raw counts)")
    logger.info("-" * 70)

    # Create a copy with raw counts in X
    adata_counts = adata.copy()
    adata_counts.X = adata_counts.layers['counts'].copy()

    # Remove normalized data from this file to avoid confusion
    del adata_counts.layers['counts']

    # Validate before saving
    if not validate_counts(adata_counts, 'X', 'raw'):
        logger.error("CRITICAL: Raw counts validation failed before saving merged_counts.h5ad")
        sys.exit(1)

    counts_file = OUTPUT_DIR / "merged_counts.h5ad"
    logger.info(f"Saving to {counts_file}...")
    adata_counts.write(counts_file)
    logger.info(f"Saved: {counts_file} ({counts_file.stat().st_size / 1e9:.2f} GB)")

    del adata_counts

    # =========================================================================
    # Step 3: Normalize if needed and save normalized file
    # =========================================================================
    logger.info("-" * 70)
    logger.info("STEP 3: Normalize and save merged_normalized.h5ad")
    logger.info("-" * 70)

    # Check if X is already normalized (from fixed files) or raw (from raw files)
    X_sample = adata.X[:1000, :100]
    if hasattr(X_sample, 'toarray'):
        X_sample = X_sample.toarray()
    is_normalized = X_sample.max() < 15 and not np.allclose(X_sample, np.round(X_sample))

    if is_normalized:
        logger.info("Data already normalized (from fixed files)")
    else:
        logger.info("Normalizing merged data...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # Validate normalized data
    if not validate_counts(adata, 'X', 'normalized'):
        logger.warning("Normalized data validation had issues - check logs")

    # Validate raw counts still in layers
    if not validate_counts(adata, "layers['counts']", 'raw'):
        logger.error("CRITICAL: Raw counts corrupted during normalization")
        sys.exit(1)

    # Compute HVGs and PCA for downstream analysis
    n_genes = adata.n_vars
    if n_genes <= 2000:
        logger.info(f"Targeted panel ({n_genes} genes) - using all genes")
        adata.var['highly_variable'] = True
    else:
        logger.info(f"Computing HVGs from {n_genes} genes...")
        sc.pp.highly_variable_genes(adata, n_top_genes=min(2000, n_genes - 1))

    n_pcs = min(30, n_genes - 1)
    logger.info(f"Computing PCA ({n_pcs} components)...")
    sc.pp.pca(adata, n_comps=n_pcs)

    logger.info("Computing neighbors and UMAP...")
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Save normalized file
    normalized_file = OUTPUT_DIR / "merged_normalized.h5ad"
    logger.info(f"Saving to {normalized_file}...")
    adata.write(normalized_file)
    logger.info(f"Saved: {normalized_file} ({normalized_file.stat().st_size / 1e9:.2f} GB)")

    # =========================================================================
    # Step 4: Batch correction
    # =========================================================================
    logger.info("-" * 70)
    logger.info("STEP 4: Batch correction and save merged_corrected.h5ad")
    logger.info("-" * 70)

    if args.skip_correction or args.method == 'none':
        logger.info("Skipping batch correction (--skip-correction or --method none)")
        adata_corrected = adata.copy()
    else:
        if args.method == 'harmony':
            adata_corrected = apply_harmony_correction(adata.copy(), args.batch_key)
        elif args.method == 'scvi':
            logger.info("scVI correction requested - checking for raw counts...")
            # scVI requires raw counts
            if not validate_counts(adata, "layers['counts']", 'raw'):
                logger.error("Cannot run scVI - raw counts not available")
                logger.info("Falling back to Harmony")
                adata_corrected = apply_harmony_correction(adata.copy(), args.batch_key)
            else:
                # scVI implementation would go here
                logger.warning("scVI not implemented in this script - using Harmony")
                adata_corrected = apply_harmony_correction(adata.copy(), args.batch_key)

    # Final validation of raw counts
    if not validate_counts(adata_corrected, "layers['counts']", 'raw'):
        logger.error("CRITICAL: Raw counts lost during batch correction")
        sys.exit(1)

    # Save corrected file
    corrected_file = OUTPUT_DIR / "merged_corrected.h5ad"
    logger.info(f"Saving to {corrected_file}...")
    adata_corrected.write(corrected_file)
    logger.info(f"Saved: {corrected_file} ({corrected_file.stat().st_size / 1e9:.2f} GB)")

    # =========================================================================
    # Step 5: Generate validation report
    # =========================================================================
    logger.info("-" * 70)
    logger.info("STEP 5: Generate validation report")
    logger.info("-" * 70)

    report = []
    report.append("# Merged Files Validation Report")
    report.append(f"\nGenerated: {datetime.now().isoformat()}")
    report.append(f"\n## Files Created\n")

    for name, fpath in [
        ("merged_counts.h5ad", counts_file),
        ("merged_normalized.h5ad", normalized_file),
        ("merged_corrected.h5ad", corrected_file)
    ]:
        adata_check = ad.read_h5ad(fpath, backed='r')
        X_sample = adata_check.X[:1000, :100]
        if hasattr(X_sample, 'toarray'):
            X_sample = X_sample.toarray()

        report.append(f"### {name}")
        report.append(f"- Shape: {adata_check.shape}")
        report.append(f"- X max: {X_sample.max():.2f}")
        report.append(f"- X mean: {X_sample.mean():.4f}")
        report.append(f"- X integer-like: {np.allclose(X_sample, np.round(X_sample))}")
        report.append(f"- Layers: {list(adata_check.layers.keys()) if adata_check.layers else 'None'}")

        if adata_check.layers and 'counts' in adata_check.layers:
            L_sample = adata_check.layers['counts'][:1000, :100]
            if hasattr(L_sample, 'toarray'):
                L_sample = L_sample.toarray()
            report.append(f"- layers['counts'] max: {L_sample.max():.0f}")
            report.append(f"- layers['counts'] integer-like: {np.allclose(L_sample, np.round(L_sample))}")

        report.append("")
        del adata_check

    report_file = OUTPUT_DIR / "VALIDATION_REPORT.md"
    with open(report_file, 'w') as f:
        f.write('\n'.join(report))
    logger.info(f"Validation report: {report_file}")

    # =========================================================================
    # Summary
    # =========================================================================
    logger.info("=" * 70)
    logger.info("REBUILD COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Output directory: {OUTPUT_DIR}")
    logger.info(f"Files created:")
    logger.info(f"  - merged_counts.h5ad (TRUE raw counts in X)")
    logger.info(f"  - merged_normalized.h5ad (normalized X, raw in layers['counts'])")
    logger.info(f"  - merged_corrected.h5ad (batch-corrected, raw in layers['counts'])")
    logger.info("")
    logger.info("NEXT STEP: Run 66_validate_raw_counts.py to verify")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
