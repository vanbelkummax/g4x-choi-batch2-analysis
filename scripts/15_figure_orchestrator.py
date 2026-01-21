#!/usr/bin/env python3
"""
Figure Orchestrator for Spatial Hackathon 2026
==============================================

Runs all figure generation scripts and verifies outputs.
Generates a summary report of all figures.

Author: Max Van Belkum
Date: 2026-01-20
"""

import subprocess
import sys
from pathlib import Path
from datetime import datetime
import json

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
FIG_DIR = OUTPUT_DIR / "figures" / "showcase_v2"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Figure manifest - Complete list of all figures
FIGURE_MANIFEST = {
    # Methodology figures (NEW)
    "fig0_sample_exclusion": {
        "title": "Sample Exclusion Rationale",
        "script": "fig0_sample_exclusion.py",
        "category": "Methodology",
        "description": "Documents YP04A exclusion due to low cell count after QC"
    },

    # Traditional figures (from existing scripts - copied from showcase/)
    "fig1_sample_overview": {
        "title": "Sample Overview Grid",
        "script": None,  # Already generated
        "category": "Traditional",
        "description": "Spatial plots of all PDAC samples with cell types"
    },
    "fig2_cell_type_composition": {
        "title": "Cell Type Composition",
        "script": None,
        "category": "Traditional",
        "description": "Cell type proportions by response group"
    },

    # Significant Findings
    "fig3_centrality_analysis": {
        "title": "Centrality Analysis",
        "script": None,
        "category": "Significant Findings",
        "description": "Betweenness centrality comparison R vs NR"
    },
    "fig4_spatial_hub_cells": {
        "title": "Spatial Hub Cells",
        "script": None,
        "category": "Significant Findings",
        "description": "Spatial distribution of hub cells"
    },
    "fig5_persistent_homology": {
        "title": "Persistent Homology",
        "script": None,
        "category": "Significant Findings",
        "description": "Topological features comparison"
    },
    "fig6_betti_curves": {
        "title": "Betti Curves",
        "script": None,
        "category": "Significant Findings",
        "description": "Betti number curves by response"
    },

    # Other analyses
    "fig7_spatial_entropy": {
        "title": "Spatial Entropy",
        "script": None,
        "category": "Other Analyses",
        "description": "Spatial entropy distribution"
    },
    "fig8_mi_biomarkers": {
        "title": "MI Biomarkers",
        "script": None,
        "category": "Other Analyses",
        "description": "Mutual information biomarker discovery"
    },
    "fig9_treatment_response": {
        "title": "Treatment Response Biomarkers",
        "script": None,
        "category": "Other Analyses",
        "description": "Top biomarkers for treatment response"
    },
    "fig10_sample_gallery": {
        "title": "Sample Gallery",
        "script": None,
        "category": "Other Analyses",
        "description": "Gallery of all samples with QC metrics"
    },
    "fig11_acinar_comparison": {
        "title": "Acinar Cell Comparison",
        "script": None,
        "category": "Other Analyses",
        "description": "Acinar cell proportions by response"
    },
    "fig12_ratio_analysis": {
        "title": "Cell Type Ratio Analysis",
        "script": None,
        "category": "Other Analyses",
        "description": "Cell type ratio comparisons"
    },
    "fig13_ncrna_analysis": {
        "title": "ncRNA Analysis",
        "script": None,
        "category": "Other Analyses",
        "description": "Non-coding RNA analysis"
    },
    "fig14_deconvolution": {
        "title": "Deconvolution Analysis",
        "script": None,
        "category": "Other Analyses",
        "description": "Cell type deconvolution results"
    },
    "fig15_ligrec": {
        "title": "Ligand-Receptor Communication",
        "script": None,
        "category": "Other Analyses",
        "description": "Cell-cell communication analysis"
    },

    # NEW PCA figures
    "fig16_bulk_rna_pca": {
        "title": "Bulk RNA PCA",
        "script": "14_pca_analyses.py",
        "category": "PCA Analysis (NEW)",
        "description": "Pseudobulk PCA of top 500 HVGs"
    },
    "fig17_spatial_pca": {
        "title": "Spatial PCA",
        "script": "14_pca_analyses.py",
        "category": "PCA Analysis (NEW)",
        "description": "PCA of spatial and topology features"
    },
    "fig18_ncrna_pca": {
        "title": "ncRNA PCA",
        "script": "14_pca_analyses.py",
        "category": "PCA Analysis (NEW)",
        "description": "PCA of non-coding RNA expression"
    },
    "fig19_ncrna_integrated_pca": {
        "title": "ncRNA Integrated PCA",
        "script": "14_pca_analyses.py",
        "category": "PCA Analysis (NEW)",
        "description": "Integrated ncRNA and cell type PCA"
    },
    "fig20_multimodal_pca": {
        "title": "Multi-modal Integrated PCA",
        "script": "14_pca_analyses.py",
        "category": "PCA Analysis (NEW)",
        "description": "All features combined: expression + cell types + spatial"
    },
    "fig21_separation_methods": {
        "title": "Separation Methods Comparison",
        "script": "14_pca_analyses.py",
        "category": "PCA Analysis (NEW)",
        "description": "Comparison of PCA, LDA, RF, SVM for R vs NR"
    },
}


def run_script(script_name: str) -> bool:
    """Run a Python script and return success status."""
    script_path = PROJECT_ROOT / "scripts" / script_name
    if not script_path.exists():
        print(f"  ERROR: Script not found: {script_path}")
        return False

    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        if result.returncode != 0:
            print(f"  ERROR: Script failed with return code {result.returncode}")
            print(f"  STDERR: {result.stderr[:500]}")
            return False
        return True
    except subprocess.TimeoutExpired:
        print(f"  ERROR: Script timed out")
        return False
    except Exception as e:
        print(f"  ERROR: {e}")
        return False


def check_figure_exists(fig_name: str) -> dict:
    """Check if figure PNG and PDF exist."""
    png_path = FIG_DIR / f"{fig_name}.png"
    pdf_path = FIG_DIR / f"{fig_name}.pdf"

    result = {
        "name": fig_name,
        "png_exists": png_path.exists(),
        "pdf_exists": pdf_path.exists(),
        "png_size": png_path.stat().st_size if png_path.exists() else 0,
        "pdf_size": pdf_path.stat().st_size if pdf_path.exists() else 0,
    }
    return result


def generate_summary_report(results: list) -> str:
    """Generate a summary report of all figures."""
    report = []
    report.append("=" * 70)
    report.append("FIGURE GENERATION SUMMARY REPORT")
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("=" * 70)
    report.append("")

    # Count statistics
    total = len(results)
    png_ok = sum(1 for r in results if r['png_exists'])
    pdf_ok = sum(1 for r in results if r['pdf_exists'])
    total_size = sum(r['png_size'] + r['pdf_size'] for r in results)

    report.append(f"SUMMARY:")
    report.append(f"  Total figures: {total}")
    report.append(f"  PNG files: {png_ok}/{total}")
    report.append(f"  PDF files: {pdf_ok}/{total}")
    report.append(f"  Total size: {total_size / 1024 / 1024:.2f} MB")
    report.append("")

    # Group by category
    categories = {}
    for fig_name, meta in FIGURE_MANIFEST.items():
        cat = meta.get('category', 'Other')
        if cat not in categories:
            categories[cat] = []
        result = next((r for r in results if r['name'] == fig_name), None)
        categories[cat].append((fig_name, meta, result))

    for cat, figures in categories.items():
        report.append(f"\n{cat.upper()}")
        report.append("-" * 40)
        for fig_name, meta, result in figures:
            status = "OK" if result and result['png_exists'] else "MISSING"
            size = f"{result['png_size']/1024:.1f}KB" if result and result['png_exists'] else "N/A"
            report.append(f"  [{status:7}] {fig_name}")
            report.append(f"           {meta['title']}: {meta['description']}")
            if result and result['png_exists']:
                report.append(f"           Size: {size}")

    report.append("")
    report.append("=" * 70)
    report.append(f"Output directory: {FIG_DIR}")
    report.append("=" * 70)

    return "\n".join(report)


def main():
    """Main orchestration function."""
    print("=" * 70)
    print("SPATIAL HACKATHON 2026 - FIGURE ORCHESTRATOR")
    print("=" * 70)
    print("")

    # Track which scripts to run
    scripts_to_run = set()
    for meta in FIGURE_MANIFEST.values():
        if meta['script'] is not None:
            scripts_to_run.add(meta['script'])

    # Run each unique script
    print("STEP 1: Running figure generation scripts...")
    print("-" * 40)
    script_results = {}
    for script in sorted(scripts_to_run):
        print(f"\nRunning {script}...")
        success = run_script(script)
        script_results[script] = success
        print(f"  {'SUCCESS' if success else 'FAILED'}")

    # Check all figure outputs
    print("\n\nSTEP 2: Verifying figure outputs...")
    print("-" * 40)
    results = []
    for fig_name in FIGURE_MANIFEST.keys():
        result = check_figure_exists(fig_name)
        results.append(result)
        status = "OK" if result['png_exists'] else "MISSING"
        print(f"  [{status}] {fig_name}")

    # Generate and save summary report
    print("\n\nSTEP 3: Generating summary report...")
    print("-" * 40)
    report = generate_summary_report(results)
    print(report)

    # Save report
    report_path = FIG_DIR / "FIGURE_MANIFEST.txt"
    with open(report_path, 'w') as f:
        f.write(report)
    print(f"\nReport saved to: {report_path}")

    # Save JSON manifest
    json_manifest = {
        "generated": datetime.now().isoformat(),
        "output_dir": str(FIG_DIR),
        "figures": results,
        "manifest": FIGURE_MANIFEST
    }
    json_path = FIG_DIR / "figure_manifest.json"
    with open(json_path, 'w') as f:
        json.dump(json_manifest, f, indent=2)
    print(f"JSON manifest saved to: {json_path}")

    # Return success if all figures exist
    all_ok = all(r['png_exists'] for r in results)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
