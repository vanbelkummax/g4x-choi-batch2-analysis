#!/usr/bin/env python3
"""HVG optimization for Block5 pilot samples."""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
from pathlib import Path

# Config
SAMPLES = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}
HVG_CONFIGS = [20, 30, 50, 'all']
INPUT = Path('/home/user/g4x-choi-batch2-analysis/results/qc_all_samples/merged/merged_corrected.h5ad')
OUTPUT = Path('/home/user/g4x-choi-batch2-analysis/results/pilot')
OUTPUT.mkdir(parents=True, exist_ok=True)

# GPU UMAP if available
try:
    from cuml.manifold import UMAP as cumlUMAP
    import cupy as cp
    def gpu_umap(X):
        return cp.asnumpy(cumlUMAP(n_neighbors=15, min_dist=0.3, verbose=False).fit_transform(cp.asarray(X.astype(np.float32))))
    print("✓ GPU enabled")
except ImportError:
    import umap
    def gpu_umap(X):
        return umap.UMAP(n_neighbors=15, min_dist=0.3).fit_transform(X)
    print("✗ CPU mode")

def process(adata, n_hvg):
    """Process one HVG config."""
    ad = adata.copy()
    if n_hvg != 'all':
        sc.pp.highly_variable_genes(ad, n_top_genes=n_hvg, flavor='seurat_v3', layer='counts', span=0.3)
        ad = ad[:, ad.var.highly_variable].copy()

    n_pcs = min(30, ad.n_vars - 1)
    sc.pp.pca(ad, n_comps=n_pcs)
    sc.pp.neighbors(ad, n_neighbors=15, n_pcs=n_pcs)
    ad.obsm['X_umap'] = gpu_umap(ad.obsm['X_pca'])
    sc.tl.leiden(ad, resolution=0.5)

    sil = silhouette_score(ad.obsm['X_pca'], ad.obs['leiden']) if ad.obs['leiden'].nunique() > 1 else 0
    return {'n_hvg': n_hvg, 'n_clusters': ad.obs['leiden'].nunique(), 'silhouette': sil, 'umap': ad.obsm['X_umap'], 'leiden': ad.obs['leiden'].values}

def main():
    print("Loading data...")
    adata = sc.read_h5ad(INPUT)
    results = []

    for sid, stage in SAMPLES.items():
        print(f"\n{sid} ({stage})")
        sample = adata[adata.obs['sample_id'] == sid].copy()

        sample_results = {}
        for n_hvg in HVG_CONFIGS:
            r = process(sample, n_hvg)
            sample_results[str(n_hvg)] = r
            results.append({'sample': sid, 'stage': stage, 'n_hvg': str(n_hvg), 'clusters': r['n_clusters'], 'silhouette': r['silhouette']})
            print(f"  HVG={n_hvg}: {r['n_clusters']} clusters, sil={r['silhouette']:.3f}")

        # Save best config
        best = max(sample_results, key=lambda x: sample_results[x]['silhouette'])
        print(f"  ★ Best: HVG={best}")

    # Save results
    df = pd.DataFrame(results)
    df.to_csv(OUTPUT / 'hvg_results.csv', index=False)
    print(f"\nSaved: {OUTPUT / 'hvg_results.csv'}")
    print(df.pivot(index='sample', columns='n_hvg', values='silhouette').round(3))

if __name__ == '__main__':
    main()
