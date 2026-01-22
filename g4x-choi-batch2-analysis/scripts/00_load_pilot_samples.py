#!/usr/bin/env python3
"""
00_load_pilot_samples.py - Load G4X pilot samples (E02, F02, G02)

Loads raw count data, preserves counts in .layers['counts'], runs basic QC.
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# === CONFIG ===
DATA_ROOT = Path('/mnt/x/Choi_Batch_2_Tuesday/g4-028-083-FC1-L002_5XGaAe5DB2dm7sRK')
OUTPUT = Path('/home/user/g4x-choi-batch2-analysis/results/pilot')
OUTPUT.mkdir(parents=True, exist_ok=True)

SAMPLES = {
    'E02': {'stage': 'Normal', 'path': DATA_ROOT / 'E02'},
    'F02': {'stage': 'Metaplasia', 'path': DATA_ROOT / 'F02'},
    'G02': {'stage': 'Cancer', 'path': DATA_ROOT / 'G02'},
}

# QC thresholds - lower for targeted 337-gene panel
MIN_GENES = 10  # Lower threshold for targeted panel
MIN_CELLS = 3
MAX_PCTL_COUNTS = 99  # Remove top 1% (potential doublets/artifacts)


def load_sample(sample_id: str, sample_info: dict) -> sc.AnnData:
    """Load a single G4X sample from raw data."""
    print(f"\n{'='*50}")
    print(f"Loading {sample_id} ({sample_info['stage']})")
    print(f"{'='*50}")

    data_path = sample_info['path'] / 'single_cell_data'

    # Load transcript counts
    print("  Loading cell_by_transcript.csv.gz...")
    counts_df = pd.read_csv(data_path / 'cell_by_transcript.csv.gz', index_col=0)

    # Remove control probes (NCP*, NCS*)
    gene_cols = [c for c in counts_df.columns if not c.startswith('NC')]
    counts_df = counts_df[gene_cols]
    print(f"  Genes after removing controls: {len(gene_cols)}")

    # Load metadata
    print("  Loading cell_metadata.csv.gz...")
    meta_df = pd.read_csv(data_path / 'cell_metadata.csv.gz', index_col=0)

    # Align indices
    common_cells = counts_df.index.intersection(meta_df.index)
    counts_df = counts_df.loc[common_cells]
    meta_df = meta_df.loc[common_cells]
    print(f"  Cells with both counts and metadata: {len(common_cells):,}")

    # Create AnnData with raw counts in .X
    adata = sc.AnnData(
        X=counts_df.values.astype(np.float32),
        obs=meta_df.copy(),
        var=pd.DataFrame(index=gene_cols)
    )
    adata.var_names = gene_cols
    adata.obs_names = common_cells.astype(str)

    # Store raw counts in layers IMMEDIATELY
    adata.layers['counts'] = adata.X.copy()

    # Extract spatial coordinates
    if 'cell_x' in meta_df.columns and 'cell_y' in meta_df.columns:
        adata.obsm['spatial'] = meta_df[['cell_x', 'cell_y']].values
        print(f"  Spatial coords: {adata.obsm['spatial'].shape}")

    # Add sample metadata
    adata.obs['sample_id'] = sample_id
    adata.obs['stage'] = sample_info['stage']

    print(f"  Initial: {adata.n_obs:,} cells x {adata.n_vars} genes")
    return adata


def run_qc(adata: sc.AnnData, sample_id: str) -> sc.AnnData:
    """Run QC metrics and filtering."""
    print(f"\n  Running QC for {sample_id}...")

    # Calculate QC metrics on raw counts
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    # QC metrics summary
    print(f"  Genes/cell: median={adata.obs['n_genes_by_counts'].median():.0f}, "
          f"range=[{adata.obs['n_genes_by_counts'].min():.0f}-{adata.obs['n_genes_by_counts'].max():.0f}]")
    print(f"  Counts/cell: median={adata.obs['total_counts'].median():.0f}")

    # Filtering
    n_before = adata.n_obs

    # Filter cells with too few genes
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)

    # Filter genes with too few cells
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)

    # Remove extreme count outliers (top 1%)
    count_thresh = np.percentile(adata.obs['total_counts'], MAX_PCTL_COUNTS)
    adata = adata[adata.obs['total_counts'] < count_thresh].copy()

    n_after = adata.n_obs
    n_genes = adata.n_vars
    print(f"  After QC: {n_after:,} cells x {n_genes} genes ({n_before - n_after:,} cells removed)")

    return adata


def normalize(adata: sc.AnnData) -> sc.AnnData:
    """Normalize counts (preserve raw in .layers['counts'])."""
    print("  Normalizing...")

    # .layers['counts'] already has raw counts from filtering
    # Normalize .X
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    return adata


def main():
    print("="*60)
    print("G4X Pilot Sample Loading")
    print("="*60)

    all_adatas = []

    for sample_id, sample_info in SAMPLES.items():
        # Load
        adata = load_sample(sample_id, sample_info)

        # QC (preserves layers)
        adata = run_qc(adata, sample_id)

        # Normalize (.X becomes log-normalized, layers['counts'] stays raw)
        adata = normalize(adata)

        # Save individual sample
        out_path = OUTPUT / f'{sample_id}_raw.h5ad'
        adata.write_h5ad(out_path)
        print(f"  Saved: {out_path}")

        all_adatas.append(adata)

    # Merge all samples
    print(f"\n{'='*50}")
    print("Merging samples...")
    print(f"{'='*50}")

    # Find common genes across all samples
    common_genes = set(all_adatas[0].var_names)
    for ad in all_adatas[1:]:
        common_genes &= set(ad.var_names)
    common_genes = sorted(common_genes)
    print(f"  Common genes: {len(common_genes)}")

    # Subset each to common genes and collect raw counts
    aligned_adatas = []
    raw_counts_list = []

    for ad in all_adatas:
        ad_aligned = ad[:, common_genes].copy()
        # Extract raw counts before they get lost in concat
        raw_counts_list.append(ad_aligned.layers['counts'].copy())
        aligned_adatas.append(ad_aligned)

    # Merge
    merged = sc.concat(aligned_adatas, label='sample_id', keys=list(SAMPLES.keys()))
    merged.obs_names_make_unique()

    # Restore raw counts layer by concatenating
    merged.layers['counts'] = np.vstack(raw_counts_list)

    merged_path = OUTPUT / 'merged_pilot.h5ad'
    merged.write_h5ad(merged_path)
    print(f"Saved merged: {merged_path}")

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for sample_id in SAMPLES.keys():
        n = (merged.obs['sample_id'] == sample_id).sum()
        print(f"  {sample_id}: {n:,} cells")
    print(f"  TOTAL: {merged.n_obs:,} cells x {merged.n_vars} genes")
    print(f"\nOutputs: {OUTPUT}")


if __name__ == '__main__':
    main()
