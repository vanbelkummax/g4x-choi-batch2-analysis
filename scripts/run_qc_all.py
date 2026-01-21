#!/usr/bin/env python3
"""
Spatial Biology Hackathon 2026 - Batch QC Runner
=================================================

Runs QC pipeline on all samples and generates reports.

Usage:
    python scripts/run_qc_all.py

Author: Max Van Belkum
Date: 2026-01-20
"""

import sys
from pathlib import Path

# Add project to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Import our modules
from scripts.data_loading import (
    load_all_samples,
    load_sample,
    PDAC_SAMPLES,
    G4X_SAMPLES,
)
from scripts.utils.qc_metrics import (
    calculate_qc_metrics,
    get_filtering_stats,
    metrics_to_dataframe,
    generate_qc_summary,
)
from scripts.utils.plotting import plot_qc_metrics, plot_spatial_qc
from scripts.preprocess import preprocess_spatial, load_thresholds

# Output paths
OUTPUT_DIR = PROJECT_ROOT / "outputs"
FIGURES_DIR = OUTPUT_DIR / "figures" / "qc"
ADATA_RAW_DIR = OUTPUT_DIR / "adata" / "raw"
ADATA_PROC_DIR = OUTPUT_DIR / "adata" / "preprocessed"
TABLES_DIR = OUTPUT_DIR / "tables"
REPORTS_DIR = OUTPUT_DIR / "reports"


def ensure_dirs():
    """Create output directories if they don't exist."""
    for d in [FIGURES_DIR, ADATA_RAW_DIR, ADATA_PROC_DIR, TABLES_DIR, REPORTS_DIR]:
        d.mkdir(parents=True, exist_ok=True)


def run_qc_pipeline():
    """
    Run full QC pipeline on all samples.

    Steps:
    1. Load all samples (8 PDAC + 2 G4X)
    2. Calculate QC metrics per sample
    3. Generate per-sample QC figures
    4. Save raw AnnData
    5. Preprocess (filter) samples
    6. Save preprocessed AnnData
    7. Generate combined metrics table
    8. Generate HTML summary report
    """
    print("=" * 70)
    print("SPATIAL BIOLOGY HACKATHON 2026 - QC PIPELINE")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    ensure_dirs()

    # Track results
    raw_metrics = []
    filtering_stats = []
    failed_samples = []

    # Define all samples
    all_samples = PDAC_SAMPLES + G4X_SAMPLES
    platform_map = {s: "Visium" for s in PDAC_SAMPLES}
    platform_map.update({s: "G4X" for s in G4X_SAMPLES})

    print(f"\nProcessing {len(all_samples)} samples: {all_samples}\n")

    for sample in all_samples:
        print(f"\n{'='*60}")
        print(f"SAMPLE: {sample}")
        print(f"{'='*60}")

        platform = platform_map[sample]

        try:
            # Load sample
            print(f"Loading {sample}...")
            adata = load_sample(sample, platform=platform)

            # Calculate QC metrics (before filtering)
            print("Calculating QC metrics...")
            adata_qc = adata.copy()
            # Adjust percent_top for samples with few features (G4X has ~387)
            n_features = adata_qc.n_vars
            if n_features < 100:
                percent_top = [10, 20, 50]
            elif n_features < 500:
                percent_top = [50, 100, 200]
            else:
                percent_top = [50, 100, 200, 500]
            sc.pp.calculate_qc_metrics(adata_qc, percent_top=percent_top, inplace=True)

            # Handle mitochondrial genes
            mito_genes = adata_qc.var_names.str.startswith('MT-') | adata_qc.var_names.str.startswith('mt-')
            if mito_genes.sum() > 0:
                adata_qc.var['mt'] = mito_genes
                sc.pp.calculate_qc_metrics(adata_qc, qc_vars=['mt'], inplace=True)
            else:
                adata_qc.obs['pct_counts_mt'] = 0.0

            # Get metrics
            metrics = calculate_qc_metrics(adata_qc, sample)
            raw_metrics.append(metrics)

            # Generate QC figure
            print("Generating QC plots...")
            thresholds = load_thresholds(platform)
            plot_qc_metrics(
                adata_qc,
                sample,
                output_path=FIGURES_DIR / f"{sample}_qc.png",
                thresholds=thresholds
            )

            # Generate spatial QC plot if coordinates available
            if 'spatial' in adata_qc.obsm:
                plot_spatial_qc(
                    adata_qc,
                    sample,
                    output_path=FIGURES_DIR / f"{sample}_spatial_qc.png"
                )

            # Save raw AnnData
            print("Saving raw AnnData...")
            adata_qc.write_h5ad(ADATA_RAW_DIR / f"{sample}.h5ad")

            # Preprocess (filter)
            print("Preprocessing...")
            adata_proc, stats = preprocess_spatial(adata, platform=platform)
            stats['sample'] = sample
            stats['platform'] = platform
            filtering_stats.append(stats)

            # Save preprocessed AnnData
            print("Saving preprocessed AnnData...")
            adata_proc.write_h5ad(ADATA_PROC_DIR / f"{sample}.h5ad")

            print(f"SUCCESS: {sample}")

        except Exception as e:
            print(f"FAILED: {sample} - {e}")
            failed_samples.append({'sample': sample, 'error': str(e)})
            import traceback
            traceback.print_exc()

    # Generate summary tables
    print(f"\n{'='*70}")
    print("GENERATING SUMMARY REPORTS")
    print(f"{'='*70}")

    # QC metrics table
    if raw_metrics:
        metrics_df = pd.DataFrame(raw_metrics)
        metrics_df.to_csv(TABLES_DIR / "qc_metrics.csv", index=False)
        print(f"Saved: {TABLES_DIR / 'qc_metrics.csv'}")

        # Print summary
        print("\n" + generate_qc_summary(metrics_df))

    # Filtering stats table
    if filtering_stats:
        filter_df = pd.DataFrame(filtering_stats)
        filter_df.to_csv(TABLES_DIR / "filtering_stats.csv", index=False)
        print(f"Saved: {TABLES_DIR / 'filtering_stats.csv'}")

    # Generate HTML report
    generate_html_report(raw_metrics, filtering_stats, failed_samples)

    # Final summary
    print(f"\n{'='*70}")
    print("PIPELINE COMPLETE")
    print(f"{'='*70}")
    print(f"Samples processed: {len(raw_metrics)}/{len(all_samples)}")
    print(f"Failed: {len(failed_samples)}")
    print(f"\nOutputs:")
    print(f"  Figures: {FIGURES_DIR}")
    print(f"  Raw data: {ADATA_RAW_DIR}")
    print(f"  Preprocessed: {ADATA_PROC_DIR}")
    print(f"  Tables: {TABLES_DIR}")
    print(f"  Report: {REPORTS_DIR / 'qc_summary.html'}")
    print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


def generate_html_report(raw_metrics, filtering_stats, failed_samples):
    """Generate HTML summary report."""
    metrics_df = pd.DataFrame(raw_metrics) if raw_metrics else pd.DataFrame()
    filter_df = pd.DataFrame(filtering_stats) if filtering_stats else pd.DataFrame()

    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>QC Summary - Spatial Hackathon 2026</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; background: #f5f5f5; }}
        h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; background: white; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
        tr:nth-child(even) {{ background-color: #f9f9f9; }}
        tr:hover {{ background-color: #f1f1f1; }}
        .summary-box {{ background: white; padding: 20px; border-radius: 8px; margin: 20px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .metric {{ font-size: 24px; font-weight: bold; color: #3498db; }}
        .label {{ color: #7f8c8d; font-size: 14px; }}
        .grid {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 20px; }}
        .figure-grid {{ display: grid; grid-template-columns: repeat(2, 1fr); gap: 20px; margin: 20px 0; }}
        .figure-card {{ background: white; padding: 10px; border-radius: 8px; text-align: center; }}
        .figure-card img {{ max-width: 100%; height: auto; }}
        .error {{ color: #e74c3c; }}
        .success {{ color: #27ae60; }}
    </style>
</head>
<body>
    <h1>Spatial Biology Hackathon 2026 - QC Summary</h1>
    <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

    <div class="summary-box">
        <div class="grid">
            <div>
                <div class="metric">{len(metrics_df)}</div>
                <div class="label">Samples Processed</div>
            </div>
            <div>
                <div class="metric">{metrics_df['n_cells'].sum():,}</div>
                <div class="label">Total Cells</div>
            </div>
            <div>
                <div class="metric">{metrics_df['median_genes'].median():,.0f}</div>
                <div class="label">Median Genes/Cell</div>
            </div>
            <div>
                <div class="metric">{len(failed_samples)}</div>
                <div class="label">Failed Samples</div>
            </div>
        </div>
    </div>

    <h2>QC Metrics (Before Filtering)</h2>
    {metrics_df.to_html(index=False, classes='table') if len(metrics_df) > 0 else '<p>No data</p>'}

    <h2>Filtering Statistics</h2>
    {filter_df.to_html(index=False, classes='table') if len(filter_df) > 0 else '<p>No data</p>'}

    <h2>QC Figures</h2>
    <div class="figure-grid">
"""

    # Add figure links
    for sample in metrics_df['sample'] if len(metrics_df) > 0 else []:
        html += f"""
        <div class="figure-card">
            <img src="../figures/qc/{sample}_qc.png" alt="{sample} QC">
            <p>{sample}</p>
        </div>
"""

    html += """
    </div>

    <h2>By Platform</h2>
"""

    if len(metrics_df) > 0:
        for platform in metrics_df['platform'].unique():
            subset = metrics_df[metrics_df['platform'] == platform]
            html += f"""
    <h3>{platform}</h3>
    <ul>
        <li>Samples: {len(subset)}</li>
        <li>Total cells: {subset['n_cells'].sum():,}</li>
        <li>Median genes/cell: {subset['median_genes'].median():,.0f}</li>
        <li>Median MT%: {subset['pct_mt_median'].median():.1f}%</li>
    </ul>
"""

    if failed_samples:
        html += """
    <h2 class="error">Failed Samples</h2>
    <ul>
"""
        for f in failed_samples:
            html += f"        <li>{f['sample']}: {f['error']}</li>\n"
        html += "    </ul>\n"

    html += """
</body>
</html>
"""

    report_path = REPORTS_DIR / "qc_summary.html"
    with open(report_path, 'w') as f:
        f.write(html)
    print(f"Saved: {report_path}")


if __name__ == "__main__":
    run_qc_pipeline()
