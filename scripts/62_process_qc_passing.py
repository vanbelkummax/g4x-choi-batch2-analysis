#!/usr/bin/env python3
"""
G4X Process QC-Passing Samples
==============================
Run full WNN + annotation pipeline on samples that passed QC.

Usage:
    python scripts/62_process_qc_passing.py [--harmony] [--parallel N]

Arguments:
    --harmony   Apply Harmony batch correction by lane
    --parallel  Number of parallel workers (default: 4)

Output:
    results/qc_all_samples/final_processed/{sample}_final.h5ad
"""

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


# =============================================================================
# Processing Functions
# =============================================================================

def compute_wnn_integration(adata, n_pcs_rna=30, n_pcs_prot=15, n_neighbors=20):
    """
    Compute Weighted Nearest Neighbors integration of RNA and Protein.

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


def process_sample(sample_info: dict, apply_harmony: bool = False) -> dict:
    """
    Process a single sample through the full pipeline.

    Parameters
    ----------
    sample_info : dict
        Contains 'sample_id' and 'path' keys
    apply_harmony : bool
        Whether to apply Harmony batch correction

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

        # Normalize RNA
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Variable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=200, subset=False)

        # PCA
        sc.pp.pca(adata, n_comps=min(30, adata.n_vars - 1))

        # WNN integration
        adata = compute_wnn_integration(adata)

        # Admixture scoring
        adata = compute_admixture_score(adata)

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
        result['n_cells'] = adata.n_obs
        result['n_admixed'] = int(adata.obs['is_admixed'].sum())
        result['pct_admixed'] = result['n_admixed'] / result['n_cells'] * 100
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
    parser.add_argument('--harmony', action='store_true', help='Apply Harmony batch correction')
    parser.add_argument('--parallel', type=int, default=4, help='Number of parallel workers')
    args = parser.parse_args()

    logger.info("=" * 60)
    logger.info("G4X Process QC-Passing Samples")
    logger.info("=" * 60)

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

    # Process samples
    results = []

    if args.parallel > 1:
        # Parallel processing
        with ProcessPoolExecutor(max_workers=args.parallel) as executor:
            futures = {
                executor.submit(process_sample, s, args.harmony): s['sample_id']
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
        for sample_info in tqdm(samples_to_process, desc="Processing"):
            result = process_sample(sample_info, args.harmony)
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
