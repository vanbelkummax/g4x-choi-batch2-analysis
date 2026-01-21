#!/usr/bin/env python3
"""
Script 20: cell2location Deconvolution
=======================================
Bayesian deconvolution with negative binomial model.
- Accounts for technical variation and platform effects
- Provides uncertainty estimates
- GPU-accelerated via PyTorch/scvi-tools

Reference: Kleshchevnikov et al. 2022, Nature Biotechnology
"""

import sys
sys.path.insert(0, '/home/user/spatial-hackathon-2026')

import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json

# cell2location imports
try:
    import cell2location
    from cell2location.models import RegressionModel
    from cell2location.utils import select_slide
    CELL2LOC_AVAILABLE = True
except ImportError:
    CELL2LOC_AVAILABLE = False
    print("Warning: cell2location not installed. Run: pip install cell2location")

from scripts.config import POLYMATHIC_DIR, PDAC_METADATA

# Derive SAMPLES and METADATA from config
SAMPLES = list(PDAC_METADATA.keys())
METADATA = PDAC_METADATA


# Configuration
OUTPUT_DIR = Path("/home/user/spatial-hackathon-2026/outputs/deconvolution")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

REFERENCE_DIR = Path("/home/user/data/hackathon/references")


def check_gpu():
    """Check GPU availability for cell2location."""
    import torch
    if torch.cuda.is_available():
        print(f"GPU available: {torch.cuda.get_device_name(0)}")
        print(f"CUDA version: {torch.version.cuda}")
        return "cuda:0"
    else:
        print("No GPU available, using CPU (will be slow)")
        return "cpu"


def load_reference_scrna(ref_path: Optional[Path] = None) -> sc.AnnData:
    """
    Load scRNA-seq reference for cell2location.

    Parameters
    ----------
    ref_path : Path, optional
        Path to h5ad file with scRNA-seq reference

    Returns
    -------
    AnnData with cell type annotations in obs['cell_type']
    """
    # Try default paths
    default_paths = [
        REFERENCE_DIR / "pdac_scrna_reference.h5ad",
        REFERENCE_DIR / "pdac_peng_2019.h5ad",
        Path("/home/user/work/polymax/references/pdac_scrna.h5ad"),
    ]

    if ref_path:
        default_paths.insert(0, ref_path)

    for path in default_paths:
        if path.exists():
            print(f"Loading reference: {path}")
            adata = sc.read_h5ad(path)

            # Ensure cell_type column exists
            if 'cell_type' not in adata.obs.columns:
                # Try common alternatives
                for col in ['celltype', 'CellType', 'cluster', 'annotation']:
                    if col in adata.obs.columns:
                        adata.obs['cell_type'] = adata.obs[col]
                        break

            if 'cell_type' in adata.obs.columns:
                print(f"  Cells: {adata.n_obs}")
                print(f"  Cell types: {adata.obs['cell_type'].nunique()}")
                return adata

    # No reference found - create synthetic
    print("No scRNA-seq reference found. Creating marker-based pseudo-reference...")
    return create_synthetic_reference()


def create_synthetic_reference() -> sc.AnnData:
    """
    Create synthetic reference based on PDAC marker genes.
    This is a fallback when no real scRNA-seq reference is available.
    """
    # PDAC cell type markers (from literature)
    markers = {
        "Ductal_cancer": ["KRT19", "KRT7", "EPCAM", "MUC1", "CEACAM5", "CEACAM6", "TFF1"],
        "Acinar": ["PRSS1", "CELA3A", "CPA1", "AMY2A", "PNLIP", "CTRB1", "PRSS2"],
        "Ductal_normal": ["KRT19", "KRT7", "CFTR", "AQP1", "SCTR"],
        "CAF_myCAF": ["ACTA2", "TAGLN", "MYH11", "ACTG2", "CNN1", "MYLK"],
        "CAF_iCAF": ["IL6", "CXCL1", "CXCL12", "LIF", "PDGFRA", "CFD", "DPT"],
        "Endothelial": ["PECAM1", "CDH5", "VWF", "CLDN5", "PLVAP", "ERG"],
        "Macrophage": ["CD68", "CD163", "CSF1R", "MRC1", "MARCO", "MSR1"],
        "T_cell": ["CD3D", "CD3E", "CD4", "CD8A", "TRAC", "TRBC2"],
        "B_cell": ["CD19", "MS4A1", "CD79A", "PAX5", "CD79B"],
        "Stellate": ["RGS5", "PDGFRB", "MCAM", "CSPG4", "ACTA2"],
    }

    # Get all marker genes
    all_markers = []
    for genes in markers.values():
        all_markers.extend(genes)
    all_markers = list(set(all_markers))

    # Create synthetic expression profiles
    n_cells_per_type = 100
    cell_types = list(markers.keys())
    n_cells = n_cells_per_type * len(cell_types)

    # Initialize expression matrix
    X = np.zeros((n_cells, len(all_markers)))
    obs_data = []

    for i, cell_type in enumerate(cell_types):
        start_idx = i * n_cells_per_type
        end_idx = start_idx + n_cells_per_type

        for j, gene in enumerate(all_markers):
            # High expression if marker for this type
            if gene in markers[cell_type]:
                X[start_idx:end_idx, j] = np.random.poisson(500, n_cells_per_type)
            else:
                # Low background expression
                X[start_idx:end_idx, j] = np.random.poisson(10, n_cells_per_type)

        obs_data.extend([cell_type] * n_cells_per_type)

    # Create AnnData
    adata = sc.AnnData(
        X=X.astype(np.float32),
        obs=pd.DataFrame({'cell_type': obs_data}, index=[f"cell_{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=all_markers)
    )

    print(f"Created synthetic reference with {adata.n_obs} cells, {len(cell_types)} cell types")
    print("WARNING: This is a simplified reference. Use real scRNA-seq for best results.")

    return adata


def train_reference_model(
    sc_adata: sc.AnnData,
    batch_key: Optional[str] = None,
    max_epochs: int = 250,
    device: str = "cuda:0"
) -> Tuple[sc.AnnData, RegressionModel]:
    """
    Train cell2location reference model on scRNA-seq data.

    Parameters
    ----------
    sc_adata : AnnData
        scRNA-seq reference data
    batch_key : str, optional
        Column for batch correction
    max_epochs : int
        Training epochs
    device : str
        CPU or CUDA device

    Returns
    -------
    Tuple of (annotated AnnData, trained model)
    """
    if not CELL2LOC_AVAILABLE:
        raise RuntimeError("cell2location not installed")

    print("\nTraining reference model...")

    # Prepare data
    sc_adata = sc_adata.copy()

    # Ensure raw counts
    if 'counts' in sc_adata.layers:
        sc_adata.X = sc_adata.layers['counts']

    # Filter genes
    sc.pp.filter_genes(sc_adata, min_cells=10)

    # Setup model
    cell2location.models.RegressionModel.setup_anndata(
        sc_adata,
        batch_key=batch_key,
        labels_key='cell_type'
    )

    # Create and train model
    ref_model = cell2location.models.RegressionModel(sc_adata)

    ref_model.train(
        max_epochs=max_epochs,
        use_gpu=(device != "cpu")
    )

    # Export posterior
    sc_adata = ref_model.export_posterior(
        sc_adata,
        sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
    )

    print("Reference model training complete")

    return sc_adata, ref_model


def run_cell2location(
    sp_adata: sc.AnnData,
    ref_adata: sc.AnnData,
    n_cells_per_location: int = 10,
    max_epochs: int = 30000,
    device: str = "cuda:0"
) -> sc.AnnData:
    """
    Run cell2location deconvolution on spatial data.

    Parameters
    ----------
    sp_adata : AnnData
        Spatial transcriptomics data
    ref_adata : AnnData
        Reference AnnData with trained signatures
    n_cells_per_location : int
        Expected cells per spot
    max_epochs : int
        Training epochs
    device : str
        CPU or CUDA device

    Returns
    -------
    AnnData with deconvolution results in obsm
    """
    if not CELL2LOC_AVAILABLE:
        raise RuntimeError("cell2location not installed")

    print("\nRunning cell2location deconvolution...")

    sp_adata = sp_adata.copy()

    # Get reference signatures
    if 'means_per_cluster_mu_fg' not in ref_adata.varm:
        raise ValueError("Reference model not trained. Run train_reference_model first.")

    inf_aver = ref_adata.varm['means_per_cluster_mu_fg']

    # Match genes
    common_genes = list(set(sp_adata.var_names) & set(ref_adata.var_names))
    print(f"Common genes: {len(common_genes)}")

    sp_adata = sp_adata[:, common_genes].copy()
    inf_aver = inf_aver.loc[common_genes, :]

    # Setup spatial model
    cell2location.models.Cell2location.setup_anndata(
        sp_adata,
        batch_key=None
    )

    # Create model
    sp_model = cell2location.models.Cell2location(
        sp_adata,
        cell_state_df=inf_aver,
        N_cells_per_location=n_cells_per_location,
        detection_alpha=200
    )

    # Train
    sp_model.train(
        max_epochs=max_epochs,
        use_gpu=(device != "cpu")
    )

    # Export results
    sp_adata = sp_model.export_posterior(
        sp_adata,
        sample_kwargs={'num_samples': 1000, 'batch_size': sp_adata.n_obs}
    )

    # Extract cell type abundances
    # Key results in: sp_adata.obsm['q05_cell_abundance_w_sf']

    print("Deconvolution complete")
    print(f"Results stored in: adata.obsm['q05_cell_abundance_w_sf']")

    return sp_adata


def visualize_deconvolution(
    sp_adata: sc.AnnData,
    sample_name: str,
    output_dir: Path,
    top_n: int = 6
):
    """
    Visualize cell2location deconvolution results.
    """
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)

    # Get abundances
    if 'q05_cell_abundance_w_sf' not in sp_adata.obsm:
        print("No deconvolution results found")
        return

    abundances = sp_adata.obsm['q05_cell_abundance_w_sf']
    cell_types = abundances.columns.tolist()

    # Convert to proportions
    proportions = abundances.div(abundances.sum(axis=1), axis=0)

    # 1. Overall composition pie chart
    fig, ax = plt.subplots(figsize=(10, 8))
    mean_props = proportions.mean()
    mean_props = mean_props[mean_props > 0.01]  # Filter < 1%
    mean_props.plot.pie(ax=ax, autopct='%1.1f%%')
    ax.set_title(f'{sample_name} - Cell Type Composition (cell2location)')
    ax.set_ylabel('')
    plt.tight_layout()
    plt.savefig(output_dir / f'{sample_name}_c2l_composition.png', dpi=150)
    plt.close()

    # 2. Spatial plots of top cell types
    top_types = mean_props.nlargest(top_n).index.tolist()

    # Check if spatial coordinates exist
    if 'spatial' in sp_adata.obsm:
        coords = sp_adata.obsm['spatial']
    else:
        print("No spatial coordinates found, skipping spatial visualization")
        return

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for i, ct in enumerate(top_types):
        if i >= len(axes):
            break

        ax = axes[i]
        scatter = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=proportions[ct],
            s=10,
            cmap='plasma'
        )
        ax.set_title(ct)
        ax.axis('equal')
        ax.axis('off')
        plt.colorbar(scatter, ax=ax, label='Proportion')

    # Hide unused axes
    for i in range(len(top_types), len(axes)):
        axes[i].axis('off')

    plt.suptitle(f'{sample_name} - Cell Type Spatial Distribution', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_dir / f'{sample_name}_c2l_spatial.png', dpi=150)
    plt.close()

    print(f"Saved visualizations for {sample_name}")


def run_simplified_deconvolution(sp_adata: sc.AnnData, markers: dict) -> pd.DataFrame:
    """
    Simplified marker-based deconvolution when cell2location is not available.
    Uses marker gene enrichment scores.
    """
    print("Running simplified marker-based deconvolution...")

    results = {}

    for cell_type, marker_genes in markers.items():
        # Filter to available genes
        available = [g for g in marker_genes if g in sp_adata.var_names]
        if not available:
            results[cell_type] = np.zeros(sp_adata.n_obs)
            continue

        # Calculate mean expression of markers
        expr = sp_adata[:, available].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray()

        # Z-score normalization
        mean_expr = expr.mean(axis=1)
        results[cell_type] = mean_expr

    # Convert to DataFrame
    df = pd.DataFrame(results, index=sp_adata.obs_names)

    # Normalize to proportions
    df = df.div(df.sum(axis=1), axis=0)
    df = df.fillna(0)

    return df


def main():
    """Run cell2location deconvolution on all samples."""
    print("=" * 60)
    print("cell2location Deconvolution Pipeline")
    print("=" * 60)

    # Check GPU
    device = check_gpu()

    # Check cell2location
    if not CELL2LOC_AVAILABLE:
        print("\ncell2location not installed. Running simplified deconvolution.")

    # Load reference
    ref_adata = load_reference_scrna()

    # Check if we have a proper reference or synthetic
    is_synthetic = ref_adata.n_obs < 500 or ref_adata.n_vars < 100

    if is_synthetic:
        print("\n" + "=" * 60)
        print("WARNING: Using synthetic/marker-based reference")
        print("For best results, provide a real PDAC scRNA-seq reference.")
        print("Suggested: Peng et al. 2019 (GSE111672)")
        print("=" * 60 + "\n")

    # Define markers for simplified deconvolution
    markers = {
        "Ductal_cancer": ["KRT19", "KRT7", "EPCAM", "MUC1", "CEACAM5"],
        "Acinar": ["PRSS1", "CELA3A", "CPA1", "AMY2A", "PNLIP"],
        "CAF_myCAF": ["ACTA2", "TAGLN", "MYH11", "ACTG2"],
        "CAF_iCAF": ["IL6", "CXCL1", "CXCL12", "LIF", "PDGFRA"],
        "Endothelial": ["PECAM1", "CDH5", "VWF", "CLDN5"],
        "Macrophage": ["CD68", "CD163", "CSF1R", "MRC1"],
        "T_cell": ["CD3D", "CD3E", "CD4", "CD8A"],
        "B_cell": ["CD19", "MS4A1", "CD79A"],
        "Stellate": ["RGS5", "PDGFRB", "MCAM"],
    }

    # Process samples
    all_results = {}

    for sample_id in SAMPLES:
        response = METADATA.get(sample_id, {}).get('response', 'Unknown')
        h5ad_path = POLYMATHIC_DIR / f"{sample_id}_polymathic.h5ad"

        if not h5ad_path.exists():
            print(f"Skipping {sample_id} (file not found)")
            continue

        print(f"\n{'=' * 40}")
        print(f"Processing {sample_id} ({response})")
        print(f"{'=' * 40}")

        try:
            # Load spatial data
            sp_adata = sc.read_h5ad(h5ad_path)
            print(f"Loaded: {sp_adata.n_obs} spots, {sp_adata.n_vars} genes")

            # Run deconvolution
            if CELL2LOC_AVAILABLE and not is_synthetic:
                # Full cell2location pipeline
                # First train reference (only once)
                if 'means_per_cluster_mu_fg' not in ref_adata.varm:
                    ref_adata, ref_model = train_reference_model(
                        ref_adata,
                        max_epochs=250,
                        device=device
                    )
                    # Save trained reference
                    ref_adata.write_h5ad(OUTPUT_DIR / "c2l_trained_reference.h5ad")

                # Run deconvolution
                sp_adata = run_cell2location(
                    sp_adata,
                    ref_adata,
                    n_cells_per_location=10,
                    max_epochs=30000,
                    device=device
                )

                # Extract proportions
                proportions = sp_adata.obsm['q05_cell_abundance_w_sf']
                proportions = proportions.div(proportions.sum(axis=1), axis=0)

            else:
                # Simplified marker-based deconvolution
                proportions = run_simplified_deconvolution(sp_adata, markers)

            # Store in adata
            sp_adata.obsm['deconv_proportions'] = proportions

            # Visualize
            visualize_deconvolution(sp_adata, sample_id, OUTPUT_DIR)

            # Save results
            sp_adata.write_h5ad(OUTPUT_DIR / f"{sample_id}_deconv.h5ad")
            proportions.to_csv(OUTPUT_DIR / f"{sample_id}_proportions.csv")

            # Store summary
            all_results[sample_id] = {
                'response': response,
                'mean_proportions': proportions.mean().to_dict(),
                'status': 'success'
            }

            print(f"✓ {sample_id} completed")

        except Exception as e:
            print(f"✗ {sample_id} failed: {e}")
            import traceback
            traceback.print_exc()
            all_results[sample_id] = {
                'response': response,
                'status': 'failed',
                'error': str(e)
            }

    # Save all results
    with open(OUTPUT_DIR / "c2l_all_results.json", 'w') as f:
        json.dump(all_results, f, indent=2)

    # Summary comparison R vs NR
    print("\n" + "=" * 60)
    print("SUMMARY: R vs NR Comparison")
    print("=" * 60)

    R_results = {k: v for k, v in all_results.items()
                 if v.get('response') == 'R' and v.get('status') == 'success'}
    NR_results = {k: v for k, v in all_results.items()
                  if v.get('response') == 'NR' and v.get('status') == 'success'}

    if R_results and NR_results:
        # Aggregate proportions
        R_props = pd.DataFrame([v['mean_proportions'] for v in R_results.values()])
        NR_props = pd.DataFrame([v['mean_proportions'] for v in NR_results.values()])

        comparison = pd.DataFrame({
            'R_mean': R_props.mean(),
            'NR_mean': NR_props.mean(),
        })
        comparison['diff'] = comparison['R_mean'] - comparison['NR_mean']
        comparison['enriched_in'] = np.where(comparison['diff'] > 0, 'R', 'NR')
        comparison = comparison.sort_values('diff', key=abs, ascending=False)

        print("\nCell type proportions:")
        print(comparison.round(3).to_string())

        comparison.to_csv(OUTPUT_DIR / "c2l_r_vs_nr_comparison.csv")

    print(f"\nResults saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
