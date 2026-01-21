#!/usr/bin/env python3
"""
Spatial Biology Hackathon 2026 - Data Loading
==============================================

Unified data loaders for all spatial platforms:
- Visium (PDAC treatment response)
- G4X Multimodal (gastric cancer)

Author: Max Van Belkum
Date: 2026-01-20
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import h5py
from pathlib import Path
from typing import Optional, Dict, List, Tuple
import yaml
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Configuration
# =============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
DATA_ROOT = Path.home() / "data" / "hackathon"
CONFIG_PATH = PROJECT_ROOT / "config"

# Platform-specific paths
VISIUM_PATH = DATA_ROOT / "PDAC.Visium" / "spaceranger.out"
G4X_PATH = DATA_ROOT / "multimodal"

# Clinical metadata for PDAC samples
PDAC_METADATA = {
    "YP03A": {"patient": "YP03", "timepoint": "Pre", "response": "NR", "outcome": "Stable Disease"},
    "YP03C": {"patient": "YP03", "timepoint": "Post", "response": "NR", "outcome": "Stable Disease"},
    "YP04A": {"patient": "YP04", "timepoint": "Pre", "response": "NR", "outcome": "Progressive Disease"},
    "YP04C": {"patient": "YP04", "timepoint": "Post", "response": "NR", "outcome": "Progressive Disease"},
    "YP12A": {"patient": "YP12", "timepoint": "Pre", "response": "R", "outcome": "Stable Disease"},
    "YP12C": {"patient": "YP12", "timepoint": "Post", "response": "R", "outcome": "Stable Disease"},
    "YP15A": {"patient": "YP15", "timepoint": "Pre", "response": "R", "outcome": "Stable Disease"},
    "YP15C": {"patient": "YP15", "timepoint": "Post", "response": "R", "outcome": "Stable Disease"},
}

PDAC_SAMPLES = list(PDAC_METADATA.keys())
G4X_SAMPLES = ["A01", "B01"]

# =============================================================================
# Visium Loaders
# =============================================================================

def load_visium(sample_name: str) -> ad.AnnData:
    """
    Load standard Visium data from SpaceRanger output.

    Args:
        sample_name: Sample directory name (e.g., "YP03A")

    Returns:
        AnnData object with spatial coordinates and clinical metadata
    """
    sample_path = VISIUM_PATH / sample_name

    if not sample_path.exists():
        raise FileNotFoundError(f"Sample path not found: {sample_path}")

    adata = sc.read_visium(
        sample_path,
        count_file="filtered_feature_bc_matrix.h5",
        load_images=True
    )
    adata.var_names_make_unique()

    # Add sample identifier
    adata.obs["sample"] = sample_name
    adata.obs["platform"] = "Visium"

    # Add clinical metadata
    if sample_name in PDAC_METADATA:
        for key, value in PDAC_METADATA[sample_name].items():
            adata.obs[key] = value

    print(f"Loaded {sample_name}: {adata.n_obs} spots, {adata.n_vars} genes")
    return adata


def load_all_visium(samples: Optional[List[str]] = None) -> Dict[str, ad.AnnData]:
    """
    Load multiple Visium samples.

    Args:
        samples: List of sample names. If None, loads all PDAC samples.

    Returns:
        Dictionary mapping sample name to AnnData
    """
    if samples is None:
        samples = PDAC_SAMPLES

    adatas = {}
    for sample in samples:
        try:
            adatas[sample] = load_visium(sample)
        except Exception as e:
            print(f"Failed to load {sample}: {e}")

    print(f"\nLoaded {len(adatas)}/{len(samples)} Visium samples")
    return adatas


def concat_visium(adatas: Dict[str, ad.AnnData]) -> ad.AnnData:
    """Concatenate multiple Visium samples into one AnnData."""
    if len(adatas) == 0:
        return None

    combined = ad.concat(
        list(adatas.values()),
        join="outer",
        label="sample",
        keys=list(adatas.keys())
    )
    print(f"Combined: {combined.n_obs} total spots")
    return combined


# =============================================================================
# G4X Loaders (H5 Feature Matrix - FAST)
# =============================================================================

def load_g4x_h5(sample: str) -> ad.AnnData:
    """
    Load G4X data from H5 feature matrix (faster than CSV).

    The H5 file contains pre-aggregated cell x gene/protein matrices.

    Args:
        sample: Sample ID (A01 or B01)

    Returns:
        AnnData with RNA counts and cell metadata
    """
    sample_path = G4X_PATH / sample / "single_cell_data"
    h5_path = sample_path / "feature_matrix.h5"

    if not h5_path.exists():
        raise FileNotFoundError(f"Feature matrix not found: {h5_path}")

    # Read H5 file structure - G4X uses AnnData-like format with X, obs, var
    with h5py.File(h5_path, 'r') as f:
        if 'X' in f and 'obs' in f and 'var' in f:
            # AnnData-like H5 format
            X = f['X'][:]  # Shape: (n_cells, n_features)

            # Read var (genes)
            var_data = {}
            for key in f['var'].keys():
                data = f['var'][key][:]
                if data.dtype.kind == 'S':  # byte strings
                    data = np.array([x.decode() for x in data])
                var_data[key] = data

            # Use gene_id as index if available
            if 'gene_id' in var_data:
                var_index = var_data.pop('gene_id')
            else:
                var_index = [f"gene_{i}" for i in range(X.shape[1])]

            var_df = pd.DataFrame(var_data, index=var_index)

            # Read obs (cells) - numerical data only to avoid decode issues
            obs_data = {}
            for key in f['obs'].keys():
                data = f['obs'][key][:]
                if data.dtype.kind != 'S':  # skip byte strings for now
                    obs_data[key] = data

            obs_df = pd.DataFrame(obs_data)
            obs_df.index = [f"cell_{i}" for i in range(X.shape[0])]

            adata = ad.AnnData(X=X, obs=obs_df, var=var_df)

        elif 'matrix' in f:
            # 10x sparse matrix format
            grp = f['matrix']
            data = grp['data'][:]
            indices = grp['indices'][:]
            indptr = grp['indptr'][:]
            barcodes = [b.decode() for b in grp['barcodes'][:]]
            features = grp['features']
            gene_names = [g.decode() for g in features['name'][:]]
            feature_types = [t.decode() for t in features['feature_type'][:]]

            from scipy.sparse import csc_matrix
            shape = (len(barcodes), len(gene_names))
            X = csc_matrix((data, indices, indptr), shape=shape).T.tocsr()

            adata = ad.AnnData(
                X=X,
                obs=pd.DataFrame(index=barcodes),
                var=pd.DataFrame({'feature_type': feature_types}, index=gene_names)
            )
        else:
            raise ValueError(f"Unrecognized H5 format: {list(f.keys())}")

    adata.var_names_make_unique()

    # Load cell metadata
    metadata_path = sample_path / "cell_metadata.csv.gz"
    if metadata_path.exists():
        meta = pd.read_csv(metadata_path, compression='gzip', index_col=0)
        # Align index
        common_cells = adata.obs_names.intersection(meta.index.astype(str))
        if len(common_cells) > 0:
            adata = adata[common_cells].copy()
            meta = meta.loc[common_cells.astype(int) if meta.index.dtype != object else common_cells]
            for col in meta.columns:
                adata.obs[col] = meta[col].values

    # Add spatial coordinates if available
    if 'cell_x_coordinate' in adata.obs.columns:
        adata.obsm['spatial'] = adata.obs[['cell_x_coordinate', 'cell_y_coordinate']].values

    # Add metadata
    adata.obs["sample"] = sample
    adata.obs["platform"] = "G4X"

    print(f"Loaded G4X {sample}: {adata.n_obs} cells, {adata.n_vars} features")
    return adata


def load_g4x_rna_csv(sample: str) -> ad.AnnData:
    """
    Load G4X RNA from cell_by_transcript.csv.gz (backup method).

    Use load_g4x_h5() for faster loading.
    """
    sample_path = G4X_PATH / sample / "single_cell_data"
    csv_path = sample_path / "cell_by_transcript.csv.gz"

    if not csv_path.exists():
        raise FileNotFoundError(f"Transcript file not found: {csv_path}")

    print(f"Loading transcript CSV (this may take longer than H5)...")
    df = pd.read_csv(csv_path, compression='gzip', index_col=0)

    adata = ad.AnnData(df)
    adata.obs["sample"] = sample
    adata.obs["platform"] = "G4X"

    print(f"Loaded G4X RNA {sample}: {adata.n_obs} cells, {adata.n_vars} genes")
    return adata


def load_g4x_protein(sample: str) -> ad.AnnData:
    """Load G4X protein data."""
    sample_path = G4X_PATH / sample / "single_cell_data"
    csv_path = sample_path / "cell_by_protein.csv.gz"

    if not csv_path.exists():
        raise FileNotFoundError(f"Protein file not found: {csv_path}")

    df = pd.read_csv(csv_path, compression='gzip', index_col=0)

    adata = ad.AnnData(df)
    adata.obs["sample"] = sample
    adata.obs["platform"] = "G4X"
    adata.var["feature_type"] = "Protein"

    print(f"Loaded G4X Protein {sample}: {adata.n_obs} cells, {adata.n_vars} proteins")
    return adata


def load_all_g4x(samples: Optional[List[str]] = None) -> Dict[str, ad.AnnData]:
    """
    Load multiple G4X samples.

    Args:
        samples: List of sample names. If None, loads all G4X samples.

    Returns:
        Dictionary mapping sample name to AnnData
    """
    if samples is None:
        samples = G4X_SAMPLES

    adatas = {}
    for sample in samples:
        try:
            adatas[sample] = load_g4x_h5(sample)
        except Exception as e:
            print(f"Failed to load {sample}: {e}")
            # Try fallback
            try:
                adatas[sample] = load_g4x_rna_csv(sample)
            except Exception as e2:
                print(f"  Fallback also failed: {e2}")

    print(f"\nLoaded {len(adatas)}/{len(samples)} G4X samples")
    return adatas


# =============================================================================
# Unified Loaders
# =============================================================================

def load_sample(sample: str, platform: Optional[str] = None) -> ad.AnnData:
    """
    Load any sample by name, auto-detecting platform.

    Args:
        sample: Sample name (e.g., "YP03A" or "A01")
        platform: Platform hint ("Visium" or "G4X"). Auto-detected if None.

    Returns:
        AnnData object
    """
    if platform is None:
        if sample in PDAC_SAMPLES:
            platform = "Visium"
        elif sample in G4X_SAMPLES:
            platform = "G4X"
        else:
            raise ValueError(f"Unknown sample: {sample}")

    if platform == "Visium":
        return load_visium(sample)
    elif platform == "G4X":
        return load_g4x_h5(sample)
    else:
        raise ValueError(f"Unknown platform: {platform}")


def load_all_samples() -> Dict[str, ad.AnnData]:
    """Load all samples from all platforms."""
    adatas = {}

    print("=" * 60)
    print("Loading all Visium samples...")
    print("=" * 60)
    adatas.update(load_all_visium())

    print("\n" + "=" * 60)
    print("Loading all G4X samples...")
    print("=" * 60)
    adatas.update(load_all_g4x())

    print(f"\n{'=' * 60}")
    print(f"TOTAL: Loaded {len(adatas)} samples")
    print(f"{'=' * 60}")

    return adatas


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("SPATIAL BIOLOGY HACKATHON 2026 - DATA LOADING TEST")
    print("=" * 60)

    # Test single sample loading
    print("\n--- Testing Visium ---")
    try:
        adata_vis = load_visium("YP03A")
        print(f"  Success! Shape: {adata_vis.shape}")
        print(f"  Clinical: response={adata_vis.obs['response'].iloc[0]}, timepoint={adata_vis.obs['timepoint'].iloc[0]}")
    except Exception as e:
        print(f"  Error: {e}")

    print("\n--- Testing G4X ---")
    try:
        adata_g4x = load_g4x_h5("A01")
        print(f"  Success! Shape: {adata_g4x.shape}")
    except Exception as e:
        print(f"  Error: {e}")

    print("\n" + "=" * 60)
    print("Data loading test complete!")
    print("=" * 60)
