#!/usr/bin/env python3
"""
G4X Full Dataset Loading
========================
Load all 32 samples across 4 lanes into individual AnnData objects.

Usage:
    conda activate enact
    python scripts/60_load_all_samples.py

Output:
    results/qc_all_samples/raw/{sample}_raw.h5ad (32 files)
"""

import os
import sys


def check_environment():
    """Verify running in correct conda environment."""
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if "enact" not in conda_prefix:
        print("ERROR: This script requires the 'enact' conda environment.")
        print(f"Current environment: {conda_prefix or 'none'}")
        print("\nRun: conda activate enact")
        sys.exit(1)


check_environment()

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import h5py
from pathlib import Path
from tqdm import tqdm
import logging
import warnings

warnings.filterwarnings('ignore')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# =============================================================================
# Configuration
# =============================================================================

DATA_ROOT = Path("/mnt/x/Choi_Batch_2_Tuesday")
OUTPUT_DIR = Path("/home/user/g4x-choi-batch2-analysis/results/qc_all_samples/raw")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Lane directory mapping
LANES = {
    'L001': 'g4-028-083-FC1-L001_5XGaAe5DB2dm7sRK',
    'L002': 'g4-028-083-FC1-L002_5XGaAe5DB2dm7sRK',
    'L003': 'g4-028-083-FC1-L003_5XGaAe5DB2dm7sRK',
    'L004': 'g4-028-083-FC1-L004_5XGaAe5DB2dm7sRK',
}

# Sample prefixes (A-H)
SAMPLE_PREFIXES = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Sample metadata
SAMPLE_METADATA = {
    'A': {'stage': 'control', 'tissue': 'CRC', 'patient': 'ctrl'},
    'B': {'stage': 'normal', 'tissue': 'adjacent_normal', 'patient': 'SNU-105'},
    'C': {'stage': 'metaplasia', 'tissue': 'IM', 'patient': 'SNU-105'},
    'D': {'stage': 'cancer', 'tissue': 'GC', 'patient': 'SNU-105'},
    'E': {'stage': 'normal', 'tissue': 'adjacent_normal', 'patient': 'SNU-107'},
    'F': {'stage': 'metaplasia', 'tissue': 'IM', 'patient': 'SNU-107'},
    'G': {'stage': 'cancer', 'tissue': 'GC', 'patient': 'SNU-107'},
    'H': {'stage': 'cancer', 'tissue': 'GC', 'patient': 'SNU-484'},
}

# Protein markers
PROTEIN_MARKERS = [
    'ATPase', 'CD11c', 'CD20', 'CD3', 'CD31', 'CD4', 'CD45',
    'CD68', 'CD8', 'FOXP3', 'HLA-DR', 'Isotype', 'KI67',
    'PD1', 'PDL1', 'PanCK', 'aSMA'
]


# =============================================================================
# Loading Functions
# =============================================================================

def load_sample(sample_dir: Path, sample_id: str, lane: str) -> ad.AnnData:
    """
    Load a single G4X sample into AnnData format.

    Parameters
    ----------
    sample_dir : Path
        Path to sample directory containing single_cell_data/
    sample_id : str
        Sample identifier (e.g., 'A01')
    lane : str
        Lane identifier (e.g., 'L001')

    Returns
    -------
    AnnData
        Loaded sample with RNA in X and protein in obsm['protein']
    """
    h5_path = sample_dir / "single_cell_data" / "feature_matrix.h5"

    if not h5_path.exists():
        raise FileNotFoundError(f"No feature matrix for {sample_id}: {h5_path}")

    with h5py.File(h5_path, 'r') as f:
        # Load expression matrix
        X = f['X'][:]  # (n_cells, n_genes)

        # Load gene names
        gene_ids = [x.decode() if isinstance(x, bytes) else x
                    for x in f['var/gene_id'][:]]

        # Load cell IDs
        cell_ids = [x.decode() if isinstance(x, bytes) else x
                    for x in f['obs/cell_id'][:]]

        # Create AnnData
        adata = ad.AnnData(X=X.astype(np.float32))
        adata.var_names = gene_ids
        adata.obs_names = [f"{sample_id}_{cid}" for cid in cell_ids]

        # Load all obs columns
        for key in f['obs'].keys():
            if key != 'cell_id':
                data = f[f'obs/{key}'][:]
                if data.dtype == object or data.dtype.kind == 'S':
                    data = [x.decode() if isinstance(x, bytes) else x for x in data]
                adata.obs[key] = data

        # Load var columns
        for key in f['var'].keys():
            if key != 'gene_id':
                data = f[f'var/{key}'][:]
                if data.dtype == object or data.dtype.kind == 'S':
                    data = [x.decode() if isinstance(x, bytes) else x for x in data]
                adata.var[key] = data

    # Separate RNA and protein
    is_protein = adata.var_names.isin(PROTEIN_MARKERS)
    protein_data = adata[:, is_protein].X.copy()
    rna_data = adata[:, ~is_protein].X.copy()

    # Create RNA-only AnnData
    adata_rna = ad.AnnData(X=rna_data)
    adata_rna.var_names = adata.var_names[~is_protein].tolist()
    adata_rna.obs = adata.obs.copy()
    adata_rna.obs_names = adata.obs_names.tolist()

    # Store protein in obsm
    adata_rna.obsm['protein'] = protein_data
    adata_rna.uns['protein_names'] = adata.var_names[is_protein].tolist()

    # Add sample metadata
    prefix = sample_id[0]
    meta = SAMPLE_METADATA.get(prefix, {})
    adata_rna.obs['sample_id'] = sample_id
    adata_rna.obs['lane'] = lane
    adata_rna.obs['stage'] = meta.get('stage', 'unknown')
    adata_rna.obs['tissue'] = meta.get('tissue', 'unknown')
    adata_rna.obs['patient'] = meta.get('patient', 'unknown')

    # Basic QC metrics
    adata_rna.obs['n_counts'] = np.array(adata_rna.X.sum(axis=1)).flatten()
    adata_rna.obs['n_genes'] = np.array((adata_rna.X > 0).sum(axis=1)).flatten()
    adata_rna.obs['n_proteins'] = np.array((protein_data > 0).sum(axis=1)).flatten()

    return adata_rna


def get_all_sample_paths():
    """Get paths to all 32 samples."""
    samples = []

    for lane_id, lane_dir_name in LANES.items():
        lane_dir = DATA_ROOT / lane_dir_name
        if not lane_dir.exists():
            logger.warning(f"Lane directory not found: {lane_dir}")
            continue

        lane_num = lane_id[-1]  # '1', '2', '3', '4'

        for prefix in SAMPLE_PREFIXES:
            sample_id = f"{prefix}0{lane_num}"
            sample_dir = lane_dir / sample_id

            if sample_dir.exists():
                samples.append({
                    'sample_id': sample_id,
                    'lane': lane_id,
                    'path': sample_dir
                })
            else:
                logger.warning(f"Sample directory not found: {sample_dir}")

    return samples


def main():
    """Load all samples and save as individual h5ad files."""
    logger.info("=" * 60)
    logger.info("G4X Full Dataset Loading")
    logger.info("=" * 60)

    # Get all sample paths
    samples = get_all_sample_paths()
    logger.info(f"Found {len(samples)} samples")

    # Load each sample
    loaded = []
    failed = []

    for sample_info in tqdm(samples, desc="Loading samples"):
        sample_id = sample_info['sample_id']
        lane = sample_info['lane']
        path = sample_info['path']

        try:
            adata = load_sample(path, sample_id, lane)

            # Save individual file
            out_path = OUTPUT_DIR / f"{sample_id}_raw.h5ad"
            adata.write(out_path)

            loaded.append({
                'sample_id': sample_id,
                'lane': lane,
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'path': str(out_path)
            })

            logger.info(f"  {sample_id}: {adata.n_obs:,} cells, {adata.n_vars} genes")

        except Exception as e:
            logger.error(f"  {sample_id}: FAILED - {e}")
            failed.append({'sample_id': sample_id, 'error': str(e)})

    # Save manifest
    manifest = pd.DataFrame(loaded)
    manifest.to_csv(OUTPUT_DIR / "loading_manifest.csv", index=False)

    if failed:
        pd.DataFrame(failed).to_csv(OUTPUT_DIR / "loading_failures.csv", index=False)

    # Summary
    logger.info("=" * 60)
    logger.info(f"Successfully loaded: {len(loaded)}/32 samples")
    logger.info(f"Total cells: {manifest['n_cells'].sum():,}")
    logger.info(f"Output: {OUTPUT_DIR}")

    if failed:
        logger.warning(f"Failed samples: {[f['sample_id'] for f in failed]}")


if __name__ == "__main__":
    main()
