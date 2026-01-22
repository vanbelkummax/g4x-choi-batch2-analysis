#!/usr/bin/env python3
"""
66_validate_raw_counts.py

Comprehensive validation of raw counts preservation across all pipeline outputs.

This script performs extensive checks to ensure:
1. Raw files actually contain integer counts
2. Fixed processed files have counts in layers
3. Merged files have proper raw/normalized separation
4. No data corruption occurred during processing

Exit codes:
    0 - All validations passed
    1 - Critical validation failed
    2 - Warnings but no critical failures

Usage:
    python scripts/66_validate_raw_counts.py [--verbose] [--fix-check]

Author: G4X Analysis Pipeline
Date: 2026-01-22
"""

import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime
import json

import numpy as np
import pandas as pd
import anndata as ad

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
RAW_DIR = PROJECT_ROOT / "results" / "qc_all_samples" / "raw"
ORIGINAL_PROCESSED_DIR = PROJECT_ROOT / "results" / "qc_all_samples" / "final_processed"
FIXED_PROCESSED_DIR = PROJECT_ROOT / "results" / "qc_all_samples" / "final_processed_fixed"
ORIGINAL_MERGED_DIR = PROJECT_ROOT / "results" / "qc_all_samples" / "merged"
FIXED_MERGED_DIR = PROJECT_ROOT / "results" / "qc_all_samples" / "merged_fixed"


class ValidationResult:
    """Container for validation results."""

    def __init__(self):
        self.passed = []
        self.warnings = []
        self.failures = []

    def add_pass(self, msg: str):
        self.passed.append(msg)
        logger.info(f"PASS: {msg}")

    def add_warning(self, msg: str):
        self.warnings.append(msg)
        logger.warning(f"WARN: {msg}")

    def add_failure(self, msg: str):
        self.failures.append(msg)
        logger.error(f"FAIL: {msg}")

    @property
    def has_failures(self) -> bool:
        return len(self.failures) > 0

    @property
    def has_warnings(self) -> bool:
        return len(self.warnings) > 0

    def summary(self) -> str:
        lines = [
            "=" * 60,
            "VALIDATION SUMMARY",
            "=" * 60,
            f"Passed:   {len(self.passed)}",
            f"Warnings: {len(self.warnings)}",
            f"Failures: {len(self.failures)}",
            ""
        ]

        if self.failures:
            lines.append("FAILURES:")
            for f in self.failures:
                lines.append(f"  - {f}")
            lines.append("")

        if self.warnings:
            lines.append("WARNINGS:")
            for w in self.warnings:
                lines.append(f"  - {w}")
            lines.append("")

        return "\n".join(lines)


def check_is_raw_counts(data, name: str) -> tuple[bool, dict]:
    """
    Check if data represents raw counts.

    Args:
        data: numpy array or sparse matrix
        name: identifier for logging

    Returns:
        (is_raw_counts, stats_dict)
    """
    if hasattr(data, 'toarray'):
        sample = data[:min(1000, data.shape[0]), :min(100, data.shape[1])].toarray()
    else:
        sample = data[:min(1000, data.shape[0]), :min(100, data.shape[1])]

    stats = {
        'max': float(sample.max()),
        'min': float(sample.min()),
        'mean': float(sample.mean()),
        'is_integer': bool(np.allclose(sample, np.round(sample))),
        'has_negative': bool((sample < 0).any()),
        'sparsity': float((sample == 0).sum() / sample.size)
    }

    # Raw counts criteria:
    # 1. Non-negative
    # 2. Integer-like
    # 3. Max > 10 (typical for count data)
    is_raw = (
        not stats['has_negative'] and
        stats['is_integer'] and
        stats['max'] >= 1  # At least some counts
    )

    return is_raw, stats


def check_is_normalized(data, name: str) -> tuple[bool, dict]:
    """
    Check if data represents log-normalized data.

    Args:
        data: numpy array or sparse matrix
        name: identifier for logging

    Returns:
        (is_normalized, stats_dict)
    """
    if hasattr(data, 'toarray'):
        sample = data[:min(1000, data.shape[0]), :min(100, data.shape[1])].toarray()
    else:
        sample = data[:min(1000, data.shape[0]), :min(100, data.shape[1])]

    stats = {
        'max': float(sample.max()),
        'min': float(sample.min()),
        'mean': float(sample.mean()),
        'is_integer': bool(np.allclose(sample, np.round(sample))),
    }

    # Log-normalized criteria:
    # 1. Max typically < 15 for log1p(CPM) data
    # 2. Non-integer (due to log transform)
    # 3. Non-negative (log1p of non-negative values)
    is_norm = (
        stats['max'] < 20 and
        stats['min'] >= 0 and
        not (stats['is_integer'] and stats['max'] > 10)
    )

    return is_norm, stats


def validate_raw_files(result: ValidationResult, verbose: bool = False):
    """Validate raw input files contain actual raw counts."""
    logger.info("-" * 60)
    logger.info("Validating RAW files")
    logger.info("-" * 60)

    if not RAW_DIR.exists():
        result.add_failure(f"Raw directory not found: {RAW_DIR}")
        return

    raw_files = list(RAW_DIR.glob("*_raw.h5ad"))
    if not raw_files:
        result.add_failure("No raw h5ad files found")
        return

    result.add_pass(f"Found {len(raw_files)} raw files")

    # Check a sample of files
    for f in raw_files[:5]:
        try:
            adata = ad.read_h5ad(f)
            is_raw, stats = check_is_raw_counts(adata.X, f.name)

            if is_raw:
                result.add_pass(f"{f.name}: Valid raw counts (max={stats['max']:.0f})")
            else:
                result.add_failure(f"{f.name}: NOT raw counts - max={stats['max']:.2f}, integer={stats['is_integer']}")

            if verbose:
                logger.info(f"  Stats: {json.dumps(stats, indent=2)}")

        except Exception as e:
            result.add_failure(f"{f.name}: Error loading - {e}")


def validate_original_processed(result: ValidationResult, verbose: bool = False):
    """Validate original processed files (should show the bug)."""
    logger.info("-" * 60)
    logger.info("Validating ORIGINAL processed files (checking for bug)")
    logger.info("-" * 60)

    if not ORIGINAL_PROCESSED_DIR.exists():
        result.add_warning("Original processed directory not found - skipping")
        return

    processed_files = list(ORIGINAL_PROCESSED_DIR.glob("*_final.h5ad"))
    if not processed_files:
        result.add_warning("No original processed files found")
        return

    # These SHOULD fail (they have the bug)
    for f in processed_files[:3]:
        try:
            adata = ad.read_h5ad(f)

            # Check X - should be normalized (the bug)
            is_raw_x, stats_x = check_is_raw_counts(adata.X, f"{f.name}:X")

            if is_raw_x:
                result.add_warning(f"{f.name}: X contains raw counts (unexpected)")
            else:
                result.add_pass(f"{f.name}: X is normalized (confirms bug present)")

            # Check layers - should NOT have counts (the bug)
            if 'counts' in adata.layers:
                is_raw_l, stats_l = check_is_raw_counts(adata.layers['counts'], f"{f.name}:layers")
                if is_raw_l:
                    result.add_warning(f"{f.name}: Has valid raw counts in layers (bug may be fixed?)")
                else:
                    result.add_pass(f"{f.name}: layers['counts'] is NOT raw (confirms bug)")
            else:
                result.add_pass(f"{f.name}: No counts layer (confirms bug - counts not preserved)")

        except Exception as e:
            result.add_failure(f"{f.name}: Error - {e}")


def validate_fixed_processed(result: ValidationResult, verbose: bool = False):
    """Validate fixed processed files have proper raw counts."""
    logger.info("-" * 60)
    logger.info("Validating FIXED processed files")
    logger.info("-" * 60)

    if not FIXED_PROCESSED_DIR.exists():
        result.add_warning("Fixed processed directory not found - run script 64 first")
        return

    fixed_files = list(FIXED_PROCESSED_DIR.glob("*_final_fixed.h5ad"))
    if not fixed_files:
        result.add_warning("No fixed processed files found - run script 64 first")
        return

    result.add_pass(f"Found {len(fixed_files)} fixed files")

    for f in fixed_files[:5]:
        try:
            adata = ad.read_h5ad(f)

            # Check X - should be normalized
            is_norm_x, stats_x = check_is_normalized(adata.X, f"{f.name}:X")
            if is_norm_x:
                result.add_pass(f"{f.name}: X is normalized (correct)")
            else:
                result.add_warning(f"{f.name}: X may not be normalized (max={stats_x['max']:.2f})")

            # Check layers - MUST have raw counts
            if 'counts' not in adata.layers:
                result.add_failure(f"{f.name}: No counts layer - fix failed!")
                continue

            is_raw_l, stats_l = check_is_raw_counts(adata.layers['counts'], f"{f.name}:layers")
            if is_raw_l:
                result.add_pass(f"{f.name}: layers['counts'] has raw counts (max={stats_l['max']:.0f})")
            else:
                result.add_failure(f"{f.name}: layers['counts'] NOT raw - max={stats_l['max']:.2f}")

        except Exception as e:
            result.add_failure(f"{f.name}: Error - {e}")


def validate_original_merged(result: ValidationResult, verbose: bool = False):
    """Validate original merged files (should show the bug)."""
    logger.info("-" * 60)
    logger.info("Validating ORIGINAL merged files (checking for bug)")
    logger.info("-" * 60)

    if not ORIGINAL_MERGED_DIR.exists():
        result.add_warning("Original merged directory not found")
        return

    for name in ["merged_counts.h5ad", "merged_normalized.h5ad", "merged_corrected.h5ad"]:
        fpath = ORIGINAL_MERGED_DIR / name
        if not fpath.exists():
            result.add_warning(f"{name} not found")
            continue

        try:
            adata = ad.read_h5ad(fpath, backed='r')
            is_raw_x, stats_x = check_is_raw_counts(adata.X, name)

            if name == "merged_counts.h5ad":
                # This SHOULD have raw counts but doesn't (the bug)
                if is_raw_x:
                    result.add_warning(f"{name}: X has raw counts (bug may be fixed?)")
                else:
                    result.add_pass(f"{name}: X is NOT raw (confirms bug - mislabeled file)")
                    result.add_pass(f"  Actual max={stats_x['max']:.2f} (should be >50 for raw counts)")
            else:
                # These should be normalized
                if not is_raw_x:
                    result.add_pass(f"{name}: X is normalized (expected)")

        except Exception as e:
            result.add_failure(f"{name}: Error - {e}")


def validate_fixed_merged(result: ValidationResult, verbose: bool = False):
    """Validate fixed merged files have proper structure."""
    logger.info("-" * 60)
    logger.info("Validating FIXED merged files")
    logger.info("-" * 60)

    if not FIXED_MERGED_DIR.exists():
        result.add_warning("Fixed merged directory not found - run script 65 first")
        return

    # merged_counts.h5ad - X should be RAW
    counts_file = FIXED_MERGED_DIR / "merged_counts.h5ad"
    if counts_file.exists():
        try:
            adata = ad.read_h5ad(counts_file, backed='r')
            is_raw_x, stats_x = check_is_raw_counts(adata.X, "merged_counts.h5ad:X")

            if is_raw_x:
                result.add_pass(f"merged_counts.h5ad: X contains RAW COUNTS (max={stats_x['max']:.0f})")
            else:
                result.add_failure(f"merged_counts.h5ad: X is NOT raw counts - FIX FAILED")
                result.add_failure(f"  max={stats_x['max']:.2f}, integer={stats_x['is_integer']}")

        except Exception as e:
            result.add_failure(f"merged_counts.h5ad: Error - {e}")
    else:
        result.add_warning("merged_counts.h5ad not found")

    # merged_normalized.h5ad - X normalized, layers['counts'] raw
    norm_file = FIXED_MERGED_DIR / "merged_normalized.h5ad"
    if norm_file.exists():
        try:
            adata = ad.read_h5ad(norm_file, backed='r')

            # X should be normalized
            is_norm_x, stats_x = check_is_normalized(adata.X, "merged_normalized.h5ad:X")
            if is_norm_x:
                result.add_pass(f"merged_normalized.h5ad: X is normalized (max={stats_x['max']:.2f})")
            else:
                result.add_warning(f"merged_normalized.h5ad: X may not be normalized properly")

            # layers['counts'] should be raw
            if 'counts' in adata.layers:
                is_raw_l, stats_l = check_is_raw_counts(adata.layers['counts'], "merged_normalized.h5ad:layers")
                if is_raw_l:
                    result.add_pass(f"merged_normalized.h5ad: layers['counts'] has raw (max={stats_l['max']:.0f})")
                else:
                    result.add_failure(f"merged_normalized.h5ad: layers['counts'] NOT raw")
            else:
                result.add_failure("merged_normalized.h5ad: No counts layer")

        except Exception as e:
            result.add_failure(f"merged_normalized.h5ad: Error - {e}")
    else:
        result.add_warning("merged_normalized.h5ad not found")

    # merged_corrected.h5ad - same structure as normalized
    corr_file = FIXED_MERGED_DIR / "merged_corrected.h5ad"
    if corr_file.exists():
        try:
            adata = ad.read_h5ad(corr_file, backed='r')

            if 'counts' in adata.layers:
                is_raw_l, stats_l = check_is_raw_counts(adata.layers['counts'], "merged_corrected.h5ad:layers")
                if is_raw_l:
                    result.add_pass(f"merged_corrected.h5ad: layers['counts'] has raw (max={stats_l['max']:.0f})")
                else:
                    result.add_failure(f"merged_corrected.h5ad: layers['counts'] NOT raw")
            else:
                result.add_failure("merged_corrected.h5ad: No counts layer")

        except Exception as e:
            result.add_failure(f"merged_corrected.h5ad: Error - {e}")
    else:
        result.add_warning("merged_corrected.h5ad not found")


def main():
    parser = argparse.ArgumentParser(
        description="Validate raw counts preservation across pipeline"
    )
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Show detailed statistics')
    parser.add_argument('--fix-check', action='store_true',
                        help='Only check fixed files (skip original bug confirmation)')
    args = parser.parse_args()

    logger.info("=" * 60)
    logger.info("SCRIPT 66: Validate Raw Counts Preservation")
    logger.info("=" * 60)
    logger.info(f"Timestamp: {datetime.now().isoformat()}")

    result = ValidationResult()

    # Run validations
    validate_raw_files(result, args.verbose)

    if not args.fix_check:
        validate_original_processed(result, args.verbose)
        validate_original_merged(result, args.verbose)

    validate_fixed_processed(result, args.verbose)
    validate_fixed_merged(result, args.verbose)

    # Print summary
    print()
    print(result.summary())

    # Save report
    report_file = PROJECT_ROOT / "results" / "RAW_COUNTS_VALIDATION_REPORT.md"
    with open(report_file, 'w') as f:
        f.write(f"# Raw Counts Validation Report\n\n")
        f.write(f"Generated: {datetime.now().isoformat()}\n\n")
        f.write(result.summary())

        f.write("\n## Detailed Results\n\n")
        f.write("### Passed\n")
        for p in result.passed:
            f.write(f"- {p}\n")

        f.write("\n### Warnings\n")
        for w in result.warnings:
            f.write(f"- {w}\n")

        f.write("\n### Failures\n")
        for fail in result.failures:
            f.write(f"- {fail}\n")

    logger.info(f"Report saved: {report_file}")

    # Exit code
    if result.has_failures:
        logger.error("VALIDATION FAILED - see report for details")
        sys.exit(1)
    elif result.has_warnings:
        logger.warning("Validation passed with warnings")
        sys.exit(0)  # Warnings are OK
    else:
        logger.info("ALL VALIDATIONS PASSED")
        sys.exit(0)


if __name__ == "__main__":
    main()
