#!/usr/bin/env python3
"""
G4X Choi Batch 2 Analysis - Spatial Hackathon 2026
==================================================
32 samples, 2.3M cells, 387 genes + 17 proteins

Data: /mnt/x/Choi_Batch_2_Tuesday/
"""

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
import h5py
import json
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Configuration
DATA_ROOT = Path("/mnt/x/Choi_Batch_2_Tuesday")
OUTPUT_DIR = Path("/home/user/spatial-hackathon-2026/results/g4x_choi_batch2")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# WSL fast copy location (if available)
FAST_DATA = Path("/home/user/data/hackathon/g4x_choi_batch2")

def get_sample_dirs():
    """Get all sample directories."""
    samples = []
    # Find lane directories (named like g4-028-083-FC1-L00X_*)
    for lane_dir in sorted(DATA_ROOT.iterdir()):
        if lane_dir.is_dir() and lane_dir.name.startswith('g4-'):
            # Get sample subdirs (A01, B01, etc.)
            for sample_dir in sorted(lane_dir.iterdir()):
                if sample_dir.is_dir() and len(sample_dir.name) == 3:  # A01, B01, etc.
                    samples.append(sample_dir)
    return samples

def load_sample(sample_dir: Path, load_images: bool = False) -> ad.AnnData:
    """Load a single G4X sample into AnnData."""
    sample_id = sample_dir.name

    # Load feature matrix (H5) - G4X format has X, obs, var directly
    h5_path = sample_dir / "single_cell_data" / "feature_matrix.h5"
    if not h5_path.exists():
        raise FileNotFoundError(f"No feature matrix for {sample_id}")

    with h5py.File(h5_path, 'r') as f:
        # Load expression matrix
        X = f['X'][:]  # (n_cells, n_genes)

        # Load gene names
        gene_ids = [x.decode() if isinstance(x, bytes) else x for x in f['var/gene_id'][:]]

        # Load cell IDs
        cell_ids = [x.decode() if isinstance(x, bytes) else x for x in f['obs/cell_id'][:]]

        # Create AnnData
        adata = ad.AnnData(X=X.astype(np.float32))
        adata.var_names = gene_ids
        adata.obs_names = cell_ids

        # Load all obs columns
        for key in f['obs'].keys():
            if key != 'cell_id':
                data = f[f'obs/{key}'][:]
                if data.dtype == object:
                    data = [x.decode() if isinstance(x, bytes) else x for x in data]
                adata.obs[key] = data

        # Load var columns
        for key in f['var'].keys():
            if key != 'gene_id':
                data = f[f'var/{key}'][:]
                if data.dtype == object:
                    data = [x.decode() if isinstance(x, bytes) else x for x in data]
                adata.var[key] = data

    # Add spatial coordinates
    if 'cell_x' in adata.obs.columns and 'cell_y' in adata.obs.columns:
        adata.obsm['spatial'] = np.column_stack([
            adata.obs['cell_x'].values.astype(float),
            adata.obs['cell_y'].values.astype(float)
        ])

    # Extract protein intensities into obsm
    protein_cols = [c for c in adata.obs.columns if c.endswith('_intensity_mean')
                    and c not in ['nuclearstain_intensity_mean', 'cytoplasmicstain_intensity_mean']]
    if protein_cols:
        adata.obsm['protein'] = adata.obs[protein_cols].values.astype(float)
        adata.uns['protein_names'] = [c.replace('_intensity_mean', '') for c in protein_cols]

    # Add sample metadata
    adata.obs['sample_id'] = sample_id
    adata.obs['lane'] = sample_dir.parent.name

    # Load run metadata
    meta_json = sample_dir / "run_meta.json"
    if meta_json.exists():
        with open(meta_json) as f:
            run_meta = json.load(f)
            adata.uns['run_meta'] = run_meta

    return adata

def qc_sample(adata: ad.AnnData) -> ad.AnnData:
    """Run QC on a sample."""
    # Basic metrics - handle case where MT genes may not exist
    # Use percent_top=None to avoid "Positions outside range of features" error
    # that occurs with sparse G4X data where many cells have zero counts
    mt_genes = adata.var_names.str.startswith('MT-')
    if mt_genes.sum() > 0:
        adata.var['mt'] = mt_genes
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)
    else:
        sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)

    # Filter cells - be lenient for 387-gene panel
    sc.pp.filter_cells(adata, min_genes=5)
    sc.pp.filter_cells(adata, min_counts=10)

    # Filter genes
    sc.pp.filter_genes(adata, min_cells=5)

    return adata

def analyze_sample(adata: ad.AnnData) -> ad.AnnData:
    """Run standard analysis on a sample."""
    # Store raw counts
    adata.layers['counts'] = adata.X.copy()

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG - use cell_ranger flavor which is more robust for targeted panels
    try:
        n_genes = min(100, adata.n_vars - 1)
        sc.pp.highly_variable_genes(adata, n_top_genes=n_genes, flavor='cell_ranger')
    except Exception as e:
        print(f"    HVG warning: {e}, marking all genes as HVG")
        adata.var['highly_variable'] = True

    # PCA
    n_comps = min(50, adata.n_vars - 1, adata.n_obs - 1)
    sc.tl.pca(adata, n_comps=n_comps)

    # Neighbors and UMAP
    n_pcs = min(30, adata.obsm['X_pca'].shape[1])
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
    sc.tl.umap(adata)

    # Clustering
    sc.tl.leiden(adata, resolution=0.5, key_added='leiden_0.5')
    sc.tl.leiden(adata, resolution=1.0, key_added='leiden_1.0')

    return adata

def spatial_analysis(adata: ad.AnnData) -> ad.AnnData:
    """Run spatial analysis with squidpy."""
    if 'spatial' not in adata.obsm:
        print("    No spatial coordinates, skipping spatial analysis")
        return adata

    try:
        # Build spatial graph with k-nearest neighbors (more robust than Delaunay for large data)
        sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=15)

        # Neighborhood enrichment (fast)
        if 'leiden_0.5' in adata.obs.columns:
            sq.gr.nhood_enrichment(adata, cluster_key='leiden_0.5')
            print(f"    Neighborhood enrichment computed")

        # Spatial autocorrelation on top variable genes only (faster)
        hvg = adata.var_names[adata.var['highly_variable']].tolist()[:20] if 'highly_variable' in adata.var.columns else []
        if hvg:
            sq.gr.spatial_autocorr(adata, mode='moran', genes=hvg, n_perms=50, n_jobs=4)
            print(f"    Moran's I computed for {len(hvg)} genes")

    except Exception as e:
        print(f"    Spatial analysis warning: {e}")

    return adata

def main():
    """Main analysis pipeline."""
    print("=" * 60)
    print("G4X Choi Batch 2 Analysis")
    print("=" * 60)

    # Get samples
    samples = get_sample_dirs()
    print(f"Found {len(samples)} samples")

    # Process first few samples for initial analysis
    results = []

    for sample_dir in tqdm(samples[:8], desc="Processing samples"):  # First 8 samples
        sample_id = sample_dir.name
        print(f"\n--- Processing {sample_id} ---")

        try:
            # Load
            adata = load_sample(sample_dir)
            print(f"  Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

            # QC
            adata = qc_sample(adata)
            print(f"  After QC: {adata.n_obs} cells")

            # Analyze
            adata = analyze_sample(adata)

            # Spatial
            adata = spatial_analysis(adata)

            # Save
            out_path = OUTPUT_DIR / f"{sample_id}_analyzed.h5ad"
            adata.write(out_path)
            print(f"  Saved: {out_path}")

            results.append({
                'sample_id': sample_id,
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'n_clusters': adata.obs['leiden_0.5'].nunique(),
                'has_protein': 'protein' in adata.obsm
            })

        except Exception as e:
            print(f"  Error: {e}")
            results.append({
                'sample_id': sample_id,
                'error': str(e)
            })

    # Save summary
    summary_df = pd.DataFrame(results)
    summary_df.to_csv(OUTPUT_DIR / "analysis_summary.csv", index=False)
    print(f"\nSummary saved to {OUTPUT_DIR / 'analysis_summary.csv'}")

    return results

if __name__ == "__main__":
    main()
