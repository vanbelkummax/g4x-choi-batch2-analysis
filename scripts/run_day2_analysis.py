#!/usr/bin/env python3
"""
Spatial Biology Hackathon 2026 - Day 2 Pipeline Runner
=======================================================

Runs the complete Day 2 analysis pipeline:
1. Clustering (04_clustering.py)
2. Spatial analysis (05_spatial_analysis.py)
3. Cell type annotation (06_annotation.py)
4. Comparative analysis (07_comparative.py)

Author: Max Van Belkum
Date: 2026-01-20
"""

import sys
import time
from pathlib import Path

# Add scripts to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from datetime import datetime


def run_full_pipeline():
    """Run complete Day 2 analysis pipeline."""

    start_time = time.time()

    print("="*70)
    print("SPATIAL BIOLOGY HACKATHON 2026 - DAY 2 ANALYSIS PIPELINE")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*70)

    # Step 1: Clustering (skip if already done)
    print("\n" + "="*70)
    print("STEP 1/4: CLUSTERING")
    print("="*70)

    clustered_dir = PROJECT_ROOT / "outputs" / "adata" / "clustered"
    if len(list(clustered_dir.glob("*_clustered.h5ad"))) >= 10:
        print("  Clustering already complete, skipping...")
    else:
        from scripts import run_clustering
        # Or import and run
        import subprocess
        subprocess.run([sys.executable, str(PROJECT_ROOT / "scripts" / "04_clustering.py")])

    # Step 2: Spatial Analysis
    print("\n" + "="*70)
    print("STEP 2/4: SPATIAL ANALYSIS")
    print("="*70)

    try:
        from scripts import run_spatial
    except ImportError:
        import subprocess
        subprocess.run([sys.executable, str(PROJECT_ROOT / "scripts" / "05_spatial_analysis.py")])

    # Step 3: Cell Type Annotation
    print("\n" + "="*70)
    print("STEP 3/4: CELL TYPE ANNOTATION")
    print("="*70)

    try:
        from scripts import run_annotation
    except ImportError:
        import subprocess
        subprocess.run([sys.executable, str(PROJECT_ROOT / "scripts" / "06_annotation.py")])

    # Step 4: Comparative Analysis
    print("\n" + "="*70)
    print("STEP 4/4: COMPARATIVE ANALYSIS")
    print("="*70)

    try:
        from scripts import run_comparative
    except ImportError:
        import subprocess
        subprocess.run([sys.executable, str(PROJECT_ROOT / "scripts" / "07_comparative.py")])

    # Summary
    elapsed = time.time() - start_time
    print("\n" + "="*70)
    print("PIPELINE COMPLETE")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total time: {elapsed/60:.1f} minutes")
    print("="*70)

    # List outputs
    print("\nGenerated outputs:")
    for subdir in ["clustered", "annotated", "spatial"]:
        dir_path = PROJECT_ROOT / "outputs" / "adata" / subdir
        if dir_path.exists():
            files = list(dir_path.glob("*.h5ad"))
            print(f"  {subdir}/: {len(files)} h5ad files")

    for subdir in ["clustering", "annotation", "spatial", "comparative"]:
        dir_path = PROJECT_ROOT / "outputs" / "figures" / subdir
        if dir_path.exists():
            files = list(dir_path.glob("*.png"))
            print(f"  figures/{subdir}/: {len(files)} figures")

    table_dir = PROJECT_ROOT / "outputs" / "tables"
    if table_dir.exists():
        files = list(table_dir.glob("*.csv"))
        print(f"  tables/: {len(files)} CSV files")


if __name__ == "__main__":
    run_full_pipeline()
