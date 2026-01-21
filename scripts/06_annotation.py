#!/usr/bin/env python3
"""
Spatial Biology Hackathon 2026 - Cell Type Annotation
======================================================

Marker-based cell type annotation for:
- PDAC (pancreatic cancer) - 8 Visium samples
- Gastric cancer (G4X) - 2 multimodal samples

Approach: Score cells using canonical markers, then assign types.

Author: Max Van Belkum
Date: 2026-01-20
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Optional, Dict, List, Tuple
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Configuration
# =============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
ADATA_DIR = OUTPUT_DIR / "adata"
FIG_DIR = OUTPUT_DIR / "figures" / "annotation"
TABLE_DIR = OUTPUT_DIR / "tables"

# Create directories
FIG_DIR.mkdir(parents=True, exist_ok=True)
(ADATA_DIR / "annotated").mkdir(parents=True, exist_ok=True)

# =============================================================================
# Cell Type Markers
# =============================================================================

# PDAC markers (from Lau lab papers and PDAC single-cell atlases)
PDAC_MARKERS = {
    'Ductal_Epithelial': ['KRT19', 'KRT7', 'EPCAM', 'MUC1', 'KRT8', 'KRT18'],
    'Acinar': ['PRSS1', 'CPA1', 'CELA3A', 'CELA2A', 'REG1A', 'PNLIP'],
    'Endocrine': ['INS', 'GCG', 'SST', 'PPY', 'CHGA', 'SYP'],
    'Stellate_PSC': ['ACTA2', 'COL1A1', 'PDGFRB', 'VIM', 'DES', 'RGS5'],
    'CAF_iCAF': ['IL6', 'CXCL12', 'PDPN', 'LIF', 'HAS1'],  # Inflammatory CAF
    'CAF_myCAF': ['ACTA2', 'TAGLN', 'MYL9', 'POSTN'],  # Myofibroblastic CAF
    'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'PLVAP', 'CLDN5'],
    'T_cells': ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'],
    'B_cells': ['CD19', 'MS4A1', 'CD79A', 'CD79B'],
    'Myeloid': ['CD68', 'CD14', 'LYZ', 'ITGAM', 'CSF1R'],
    'Macrophage': ['CD68', 'CD163', 'MRC1', 'MSR1'],
    'NK_cells': ['NKG7', 'GNLY', 'NCAM1', 'KLRD1'],
}

# Gastric cancer markers (for G4X targeted panel)
GASTRIC_MARKERS = {
    'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1'],
    'Stromal': ['VIM', 'COL1A1', 'COL1A2', 'DCN', 'FAP'],
    'Endothelial': ['PECAM1', 'CDH5', 'VWF', 'PLVAP'],
    'Immune': ['PTPRC', 'CD45', 'CD3D', 'CD4', 'CD8A'],
    'Myeloid': ['CD68', 'CD14', 'LYZ'],
    'Smooth_Muscle': ['ACTA2', 'MYH11', 'TAGLN'],
}

# =============================================================================
# Marker Scoring Functions
# =============================================================================

def score_cell_types(
    adata: ad.AnnData,
    marker_dict: Dict[str, List[str]],
    layer: Optional[str] = None
) -> ad.AnnData:
    """
    Score cells for each cell type using marker genes.

    Uses sc.tl.score_genes which computes the average expression of a gene set
    minus the average of reference genes.

    Args:
        adata: AnnData with normalized expression
        marker_dict: Dictionary of cell_type -> list of marker genes
        layer: Layer to use for scoring (None = .X)

    Returns:
        AnnData with cell type scores in .obs
    """
    print(f"  Scoring {len(marker_dict)} cell types...")

    for cell_type, markers in marker_dict.items():
        # Filter to available markers
        available = [m for m in markers if m in adata.var_names]

        if len(available) == 0:
            print(f"    {cell_type}: No markers found, skipping")
            continue

        print(f"    {cell_type}: {len(available)}/{len(markers)} markers")

        # Score
        score_name = f'score_{cell_type}'
        try:
            sc.tl.score_genes(
                adata,
                gene_list=available,
                score_name=score_name,
                ctrl_size=min(50, adata.n_vars - len(available)),
                use_raw=False
            )
        except Exception as e:
            print(f"    Warning: Scoring failed for {cell_type}: {e}")
            adata.obs[score_name] = 0.0

    return adata


def assign_cell_types(
    adata: ad.AnnData,
    marker_dict: Dict[str, List[str]],
    min_score_diff: float = 0.1
) -> ad.AnnData:
    """
    Assign cell types based on highest score.

    Args:
        adata: AnnData with cell type scores
        marker_dict: Dictionary of cell types to get score names
        min_score_diff: Minimum score difference to assign a type

    Returns:
        AnnData with 'cell_type' in .obs
    """
    score_cols = [f'score_{ct}' for ct in marker_dict.keys() if f'score_{ct}' in adata.obs]

    if len(score_cols) == 0:
        print("  Warning: No scores found, cannot assign cell types")
        adata.obs['cell_type'] = 'Unknown'
        return adata

    # Get scores matrix
    scores_df = adata.obs[score_cols].copy()

    # Assign type with highest score
    max_scores = scores_df.max(axis=1)
    max_types = scores_df.idxmax(axis=1)

    # Clean up names (remove 'score_' prefix)
    cell_types = max_types.str.replace('score_', '')

    # Mark low-confidence assignments as 'Unknown'
    second_max = scores_df.apply(lambda x: x.nlargest(2).iloc[-1] if len(x.nlargest(2)) > 1 else 0, axis=1)
    low_conf = (max_scores - second_max) < min_score_diff
    cell_types[low_conf] = 'Low_Confidence'

    # Mark negative max scores as Unknown
    cell_types[max_scores < 0] = 'Unknown'

    adata.obs['cell_type'] = cell_types.astype('category')

    # Summary
    print("\n  Cell type assignments:")
    for ct, count in adata.obs['cell_type'].value_counts().items():
        print(f"    {ct}: {count} ({100*count/adata.n_obs:.1f}%)")

    return adata


# =============================================================================
# Visualization Functions
# =============================================================================

def plot_marker_dotplot(
    adata: ad.AnnData,
    marker_dict: Dict[str, List[str]],
    sample_name: str,
    groupby: str = 'cell_type'
) -> None:
    """Plot dotplot of markers by cell type."""

    # Collect available markers
    markers = []
    for ct, genes in marker_dict.items():
        markers.extend([g for g in genes if g in adata.var_names])
    markers = list(dict.fromkeys(markers))  # Deduplicate preserving order

    if len(markers) == 0:
        print(f"  No markers available for dotplot")
        return

    try:
        fig, ax = plt.subplots(figsize=(15, 8))
        sc.pl.dotplot(
            adata,
            var_names=markers[:30],  # Limit to 30 markers
            groupby=groupby,
            ax=ax,
            show=False,
            title=f'{sample_name} - Marker Expression'
        )
        plt.tight_layout()
        plt.savefig(FIG_DIR / f"{sample_name}_marker_dotplot.png", dpi=150, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"  Warning: Dotplot failed: {e}")


def plot_cell_types(
    adata: ad.AnnData,
    sample_name: str
) -> None:
    """Plot cell type assignments on UMAP and spatial."""

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # UMAP colored by cell type
    sc.pl.umap(adata, color='cell_type', ax=axes[0], show=False,
               title=f'{sample_name} - Cell Types (UMAP)')

    # Spatial plot
    if 'spatial' in adata.obsm:
        sc.pl.embedding(adata, basis='spatial', color='cell_type', ax=axes[1],
                       show=False, title=f'{sample_name} - Cell Types (Spatial)')
    else:
        axes[1].text(0.5, 0.5, 'No spatial coordinates', ha='center', va='center')
        axes[1].set_title('Spatial (N/A)')

    # Cell type proportions
    ct_counts = adata.obs['cell_type'].value_counts()
    axes[2].barh(ct_counts.index.astype(str), ct_counts.values)
    axes[2].set_xlabel('Count')
    axes[2].set_title(f'{sample_name} - Cell Type Proportions')

    plt.tight_layout()
    plt.savefig(FIG_DIR / f"{sample_name}_cell_types.png", dpi=150, bbox_inches='tight')
    plt.close()


def plot_score_heatmap(
    adata: ad.AnnData,
    marker_dict: Dict[str, List[str]],
    sample_name: str,
    groupby: str = 'leiden'
) -> None:
    """Plot heatmap of cell type scores by cluster."""

    score_cols = [f'score_{ct}' for ct in marker_dict.keys() if f'score_{ct}' in adata.obs]

    if len(score_cols) == 0:
        return

    # Average scores per cluster
    mean_scores = adata.obs.groupby(groupby)[score_cols].mean()

    # Clean up column names
    mean_scores.columns = [c.replace('score_', '') for c in mean_scores.columns]

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(mean_scores.T, cmap='RdBu_r', center=0, ax=ax,
                xticklabels=True, yticklabels=True)
    ax.set_title(f'{sample_name} - Cell Type Scores by Cluster')
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Cell Type Score')

    plt.tight_layout()
    plt.savefig(FIG_DIR / f"{sample_name}_score_heatmap.png", dpi=150, bbox_inches='tight')
    plt.close()


# =============================================================================
# Main Pipeline
# =============================================================================

def annotate_sample(
    adata: ad.AnnData,
    sample_name: str,
    platform: Optional[str] = None,
    save: bool = True
) -> ad.AnnData:
    """
    Full annotation pipeline for one sample.

    Args:
        adata: Clustered AnnData
        sample_name: Sample identifier
        platform: Platform type (Visium or G4X)
        save: Whether to save outputs

    Returns:
        Annotated AnnData
    """
    print(f"\n{'='*60}")
    print(f"Annotating: {sample_name}")
    print(f"{'='*60}")

    # Determine platform and markers
    if platform is None:
        platform = adata.obs['platform'].iloc[0] if 'platform' in adata.obs else 'Visium'

    if platform == 'G4X':
        marker_dict = GASTRIC_MARKERS
    else:
        marker_dict = PDAC_MARKERS

    print(f"  Platform: {platform}, using {len(marker_dict)} cell type definitions")

    # 1. Score cell types
    print("\n1. Scoring cell types...")
    adata = score_cell_types(adata, marker_dict)

    # 2. Assign cell types
    print("\n2. Assigning cell types...")
    adata = assign_cell_types(adata, marker_dict)

    # 3. Generate plots
    if save:
        print("\n3. Generating plots...")
        plot_cell_types(adata, sample_name)
        plot_score_heatmap(adata, marker_dict, sample_name)
        plot_marker_dotplot(adata, marker_dict, sample_name)

    # 4. Save
    if save:
        out_path = ADATA_DIR / "annotated" / f"{sample_name}_annotated.h5ad"
        adata.write_h5ad(out_path)
        print(f"\n  Saved to: {out_path}")

    return adata


def annotate_all_samples(
    input_dir: Optional[Path] = None,
    samples: Optional[List[str]] = None
) -> Dict[str, ad.AnnData]:
    """
    Annotate all clustered samples.

    Args:
        input_dir: Directory with clustered h5ad files
        samples: Specific samples to process (None = all)

    Returns:
        Dictionary of annotated AnnData objects
    """
    if input_dir is None:
        input_dir = ADATA_DIR / "clustered"

    # Find all clustered h5ad files
    h5ad_files = list(input_dir.glob("*_clustered.h5ad"))

    if samples:
        h5ad_files = [f for f in h5ad_files if f.stem.replace('_clustered', '') in samples]

    print(f"Found {len(h5ad_files)} samples to annotate")

    annotated = {}
    proportions_list = []

    for h5ad_path in sorted(h5ad_files):
        sample_name = h5ad_path.stem.replace('_clustered', '')

        try:
            # Load
            adata = sc.read_h5ad(h5ad_path)

            # Annotate
            adata = annotate_sample(adata, sample_name)
            annotated[sample_name] = adata

            # Collect proportions
            ct_props = adata.obs['cell_type'].value_counts(normalize=True)
            prop_dict = {'sample': sample_name, 'platform': adata.obs['platform'].iloc[0]}
            prop_dict.update(ct_props.to_dict())
            proportions_list.append(prop_dict)

        except Exception as e:
            print(f"  ERROR annotating {sample_name}: {e}")
            import traceback
            traceback.print_exc()

    # Save proportions table
    if proportions_list:
        props_df = pd.DataFrame(proportions_list).fillna(0)
        props_df.to_csv(TABLE_DIR / "cell_type_proportions.csv", index=False)
        print(f"\nCell type proportions saved to: {TABLE_DIR / 'cell_type_proportions.csv'}")

    return annotated


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Annotate cell types in clustered samples")
    parser.add_argument("--sample", type=str, help="Single sample to process")
    parser.add_argument("--input-dir", type=str, help="Input directory with clustered h5ad files")
    args = parser.parse_args()

    print("="*60)
    print("SPATIAL BIOLOGY HACKATHON 2026 - CELL TYPE ANNOTATION")
    print("="*60)

    if args.sample:
        # Single sample mode
        input_dir = Path(args.input_dir) if args.input_dir else ADATA_DIR / "clustered"
        adata = sc.read_h5ad(input_dir / f"{args.sample}_clustered.h5ad")
        annotate_sample(adata, args.sample)
    else:
        # Batch mode
        input_dir = Path(args.input_dir) if args.input_dir else None
        annotate_all_samples(input_dir)

    print("\n" + "="*60)
    print("ANNOTATION COMPLETE")
    print("="*60)
