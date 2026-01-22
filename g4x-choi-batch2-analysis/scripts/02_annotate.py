#!/usr/bin/env python3
"""Marker-based cell type annotation for Block5 pilot."""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

SAMPLES = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}
INPUT = Path('/home/user/g4x-choi-batch2-analysis/results/qc_all_samples/merged/merged_corrected.h5ad')
OUTPUT = Path('/home/user/g4x-choi-batch2-analysis/results/pilot')
OUTPUT.mkdir(parents=True, exist_ok=True)

# Cell type markers (subset present in G4X panel)
MARKERS = {
    'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'CDH1'],
    'Pit_Mucous': ['MUC5AC', 'TFF1', 'GKN1'],
    'Goblet': ['MUC2', 'TFF3'],
    'Intestinal_Meta': ['CDX2', 'VIL1', 'MUC2'],
    'T_Cell': ['CD3D', 'CD3E', 'CD2'],
    'Macrophage': ['CD68', 'CD163', 'CSF1R'],
    'B_Cell': ['CD19', 'MS4A1', 'CD79A'],
    'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'VIM'],
    'CAF': ['FAP', 'ACTA2', 'PDGFRA'],
    'Endothelial': ['PECAM1', 'VWF', 'CDH5'],
}

# GPU UMAP
try:
    from cuml.manifold import UMAP as cumlUMAP
    import cupy as cp
    def gpu_umap(X):
        return cp.asnumpy(cumlUMAP(n_neighbors=15, min_dist=0.3, verbose=False).fit_transform(cp.asarray(X.astype(np.float32))))
except ImportError:
    import umap
    def gpu_umap(X):
        return umap.UMAP(n_neighbors=15, min_dist=0.3).fit_transform(X)

def score_cells(adata):
    """Score each cell for each cell type."""
    genes = set(adata.var_names)
    for ct, markers in MARKERS.items():
        valid = [m for m in markers if m in genes]
        if valid:
            sc.tl.score_genes(adata, valid, score_name=f's_{ct}')
        else:
            adata.obs[f's_{ct}'] = 0
    return adata

def assign_types(adata):
    """Assign cell type by max score."""
    score_cols = [c for c in adata.obs.columns if c.startswith('s_')]
    scores = adata.obs[score_cols].values
    types = [c[2:] for c in score_cols]  # Remove 's_' prefix
    adata.obs['cell_type'] = [types[i] for i in np.argmax(scores, axis=1)]
    adata.obs['type_score'] = np.max(scores, axis=1)
    return adata

def main():
    print("Loading data...")
    adata = sc.read_h5ad(INPUT)

    for sid, stage in SAMPLES.items():
        print(f"\n{sid} ({stage})")
        sample = adata[adata.obs['sample_id'] == sid].copy()

        # Cluster with optimal HVG (use 50 as default)
        sc.pp.highly_variable_genes(sample, n_top_genes=50, flavor='seurat_v3', layer='counts', span=0.3)
        hvg = sample[:, sample.var.highly_variable].copy()
        sc.pp.pca(hvg, n_comps=30)
        sc.pp.neighbors(hvg, n_neighbors=15)
        hvg.obsm['X_umap'] = gpu_umap(hvg.obsm['X_pca'])
        sc.tl.leiden(hvg, resolution=0.5)

        # Transfer to full data
        sample.obs['leiden'] = hvg.obs['leiden']
        sample.obsm['X_umap'] = hvg.obsm['X_umap']

        # Annotate
        sample = score_cells(sample)
        sample = assign_types(sample)

        # Save
        sample.write_h5ad(OUTPUT / f'{sid}_annotated.h5ad')

        # Summary
        print(f"  Clusters: {sample.obs['leiden'].nunique()}")
        print("  Cell types:")
        for ct, n in sample.obs['cell_type'].value_counts().head(5).items():
            print(f"    {ct}: {n:,} ({n/len(sample)*100:.1f}%)")

    print(f"\nSaved to: {OUTPUT}")

if __name__ == '__main__':
    main()
