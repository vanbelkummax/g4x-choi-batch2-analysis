#!/usr/bin/env python3
"""
G4X Full Dataset QC Pre-flight Check
=====================================
Verify all prerequisites are met before running the QC pipeline.

Usage:
    python scripts/preflight_check.py

Returns exit code 0 if ready, 1 if issues found.
"""

import sys
from pathlib import Path


def check_data_directories():
    """Verify data directories exist with expected samples."""
    print("\n1. Checking data directories...")

    DATA_ROOT = Path("/mnt/x/Choi_Batch_2_Tuesday")
    LANES = {
        'L001': 'g4-028-083-FC1-L001_5XGaAe5DB2dm7sRK',
        'L002': 'g4-028-083-FC1-L002_5XGaAe5DB2dm7sRK',
        'L003': 'g4-028-083-FC1-L003_5XGaAe5DB2dm7sRK',
        'L004': 'g4-028-083-FC1-L004_5XGaAe5DB2dm7sRK',
    }
    SAMPLES = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

    issues = []
    samples_found = 0

    for lane_id, lane_dir_name in LANES.items():
        lane_dir = DATA_ROOT / lane_dir_name
        if not lane_dir.exists():
            issues.append(f"Lane directory missing: {lane_dir}")
            continue

        lane_num = lane_id[-1]
        for prefix in SAMPLES:
            sample_id = f"{prefix}0{lane_num}"
            sample_dir = lane_dir / sample_id
            h5_file = sample_dir / "single_cell_data" / "feature_matrix.h5"

            if not sample_dir.exists():
                issues.append(f"Sample directory missing: {sample_id}")
            elif not h5_file.exists():
                issues.append(f"H5 file missing for {sample_id}")
            else:
                samples_found += 1

    print(f"   Found {samples_found}/32 samples")

    # Check Resolve baseline
    baseline_path = DATA_ROOT / "choi_preGC_b2_core_metrics.csv"
    if baseline_path.exists():
        print(f"   Resolve baseline: OK")
    else:
        issues.append(f"Resolve baseline missing: {baseline_path}")

    return issues


def check_environment():
    """Verify required packages are available."""
    print("\n2. Checking Python environment...")

    issues = []
    packages = [
        ('scanpy', 'sc'),
        ('anndata', 'ad'),
        ('squidpy', 'sq'),
        ('numpy', 'np'),
        ('pandas', 'pd'),
        ('h5py', 'h5py'),
        ('sklearn', 'sklearn'),
        ('scipy', 'scipy'),
        ('matplotlib', 'plt'),
        ('seaborn', 'sns'),
        ('tqdm', 'tqdm'),
    ]

    for pkg_name, _ in packages:
        try:
            __import__(pkg_name)
            print(f"   {pkg_name}: OK")
        except ImportError:
            issues.append(f"Missing package: {pkg_name}")

    return issues


def check_output_directories():
    """Verify output directories exist or can be created."""
    print("\n3. Checking output directories...")

    REPO_ROOT = Path("/home/user/g4x-choi-batch2-analysis")
    dirs_to_check = [
        REPO_ROOT / "results" / "qc_all_samples" / "raw",
        REPO_ROOT / "results" / "qc_all_samples" / "figures" / "per_sample",
        REPO_ROOT / "results" / "qc_all_samples" / "figures" / "cross_sample",
        REPO_ROOT / "results" / "qc_all_samples" / "figures" / "batch_effects",
        REPO_ROOT / "results" / "qc_all_samples" / "final_processed",
        REPO_ROOT / "logs",
    ]

    issues = []
    for d in dirs_to_check:
        try:
            d.mkdir(parents=True, exist_ok=True)
            print(f"   {d.relative_to(REPO_ROOT)}: OK")
        except Exception as e:
            issues.append(f"Cannot create directory {d}: {e}")

    return issues


def check_scripts():
    """Verify pipeline scripts exist and are valid Python."""
    print("\n4. Checking pipeline scripts...")

    REPO_ROOT = Path("/home/user/g4x-choi-batch2-analysis")
    scripts = [
        REPO_ROOT / "scripts" / "60_load_all_samples.py",
        REPO_ROOT / "scripts" / "61_comprehensive_qc.py",
        REPO_ROOT / "scripts" / "62_process_qc_passing.py",
    ]

    issues = []
    for script in scripts:
        if not script.exists():
            issues.append(f"Script missing: {script.name}")
            continue

        # Check syntax
        try:
            with open(script) as f:
                code = f.read()
            compile(code, script, 'exec')
            print(f"   {script.name}: OK")
        except SyntaxError as e:
            issues.append(f"Syntax error in {script.name}: {e}")

    return issues


def check_disk_space():
    """Verify sufficient disk space available."""
    print("\n5. Checking disk space...")

    import shutil

    output_path = Path("/home/user/g4x-choi-batch2-analysis/results")
    total, used, free = shutil.disk_usage(output_path)

    free_gb = free // (1024**3)
    required_gb = 15  # Estimate for full dataset

    print(f"   Available: {free_gb} GB")
    print(f"   Required: ~{required_gb} GB")

    if free_gb < required_gb:
        return [f"Low disk space: {free_gb} GB < {required_gb} GB required"]

    return []


def main():
    """Run all preflight checks."""
    print("=" * 60)
    print("G4X Full Dataset QC - Pre-flight Check")
    print("=" * 60)

    all_issues = []

    all_issues.extend(check_data_directories())
    all_issues.extend(check_environment())
    all_issues.extend(check_output_directories())
    all_issues.extend(check_scripts())
    all_issues.extend(check_disk_space())

    print("\n" + "=" * 60)

    if all_issues:
        print("ISSUES FOUND:")
        for issue in all_issues:
            print(f"  - {issue}")
        print(f"\nTotal issues: {len(all_issues)}")
        print("Please fix issues before running the pipeline.")
        return 1
    else:
        print("ALL CHECKS PASSED")
        print("\nReady to run:")
        print("  python scripts/60_load_all_samples.py 2>&1 | tee logs/60_loading.log")
        print("  python scripts/61_comprehensive_qc.py 2>&1 | tee logs/61_qc.log")
        print("  python scripts/62_process_qc_passing.py --parallel 8 2>&1 | tee logs/62_processing.log")
        return 0


if __name__ == "__main__":
    sys.exit(main())
