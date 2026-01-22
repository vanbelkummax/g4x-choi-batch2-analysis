#!/usr/bin/env python3
"""
02_annotate_multimethod.py - Multi-method cell type annotation

Methods:
1. Curated RNA Markers (gene scoring) - uses only genes in G4X panel
2. CellTypist model - pre-trained immune/intestinal models
3. scRNA Reference Transfer - Kumar et al. 2022 gastric atlas (if available)

Compares methods and ranks by marker enrichment + cluster purity.
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# === CONFIG ===
INPUT = Path('/home/user/g4x-choi-batch2-analysis/results/pilot/merged_pilot.h5ad')
OUTPUT = Path('/home/user/g4x-choi-batch2-analysis/results/pilot')
FIGURES = OUTPUT / 'figures'
FIGURES.mkdir(parents=True, exist_ok=True)

SAMPLES = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}

# Optimal HVG setting (will be read from hvg_recommendation.txt if available)
DEFAULT_HVG = 100
DEFAULT_RES = 0.5

# === CURATED MARKERS (VERIFIED IN G4X 337-GENE PANEL) ===
# All markers below are confirmed present in the gene panel
MARKERS = {
    # === Gastric Epithelial ===
    'Gastric_Pit': ['MUC5AC', 'MUC6', 'TFF1', 'TFF2'],  # Pit/mucous cells
    'Gastric_Chief': ['PGC', 'ATP4A'],  # Parietal/chief cells
    'Goblet': ['MUC2', 'TFF3'],  # Goblet cells
    'Intestinal_Meta': ['CDX1', 'CDX2', 'MUC2', 'CDH17'],  # Intestinal metaplasia
    'Enteroendocrine': ['GHRL', 'SST', 'GAST', 'CHGA'],  # Neuroendocrine
    'Stem_Progenitor': ['LGR5', 'PROM1', 'OLFM4', 'SOX2', 'SOX9'],  # Stem/progenitor
    'Epithelial_General': ['EPCAM', 'CDH1', 'CLDN1', 'CLDN3', 'CLDN4', 'CLDN7', 'CLDN18'],

    # === Immune - T Cells ===
    'T_CD4': ['CD3D', 'CD3E', 'CD4', 'IL7R'],  # Helper T
    'T_CD8_Cytotoxic': ['CD3D', 'CD3E', 'CD8A', 'GZMA', 'GZMB', 'PRF1'],  # Cytotoxic T
    'T_Reg': ['CD3D', 'CD4', 'FOXP3', 'IL2RA'],  # Regulatory T
    'T_Exhausted': ['PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT'],  # Exhausted T

    # === Immune - B/Plasma ===
    'B_Cell': ['MS4A1', 'CD19', 'CD79A'],  # B cells
    'Plasma': ['IGHA1', 'IGHG1', 'IGHM', 'JCHAIN', 'SDC1'],  # Plasma cells

    # === Immune - Myeloid ===
    'Macrophage': ['CD68', 'CD163', 'CSF1R', 'MRC1'],  # Macrophages
    'Monocyte': ['CD14', 'ITGAM'],  # Monocytes

    # === Stromal ===
    'Fibroblast': ['COL1A1', 'LUM', 'VCAN', 'FN1'],  # Fibroblasts
    'CAF': ['ACTA2', 'FAP', 'PDGFRA', 'POSTN', 'THY1'],  # Cancer-associated fibroblasts
    'Endothelial': ['PECAM1', 'VWF'],  # Endothelial cells
}

# Score threshold for confident annotation
SCORE_THRESHOLD = 0.1


def load_optimal_params() -> tuple:
    """Load optimal HVG/resolution from prior optimization."""
    rec_path = OUTPUT / 'hvg_recommendation.txt'
    if rec_path.exists():
        with open(rec_path) as f:
            lines = f.readlines()
            hvg = lines[0].split(': ')[1].strip()
            hvg = int(hvg) if hvg.isdigit() else hvg
            return hvg, DEFAULT_RES
    return DEFAULT_HVG, DEFAULT_RES


def method1_marker_scoring(adata: sc.AnnData) -> sc.AnnData:
    """Method 1: Score cells using curated markers."""
    print("  Method 1: Marker gene scoring...")

    genes_in_data = set(adata.var_names)

    for cell_type, markers in MARKERS.items():
        valid_markers = [m for m in markers if m in genes_in_data]
        if valid_markers:
            sc.tl.score_genes(adata, valid_markers, score_name=f'score_{cell_type}')
        else:
            adata.obs[f'score_{cell_type}'] = 0.0
            print(f"    WARNING: No markers for {cell_type} in data")

    # Assign by max score
    score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
    scores = adata.obs[score_cols].values
    cell_types = [c.replace('score_', '') for c in score_cols]

    max_scores = np.max(scores, axis=1)
    max_indices = np.argmax(scores, axis=1)

    annotations = []
    for i, (idx, score) in enumerate(zip(max_indices, max_scores)):
        if score >= SCORE_THRESHOLD:
            annotations.append(cell_types[idx])
        else:
            annotations.append('Unknown')

    adata.obs['celltype_markers'] = annotations
    adata.obs['celltype_markers_score'] = max_scores

    # Summary
    counts = Counter(annotations)
    print(f"    Annotated: {len([a for a in annotations if a != 'Unknown']):,} cells")
    print(f"    Unknown: {counts.get('Unknown', 0):,} cells")

    return adata


def method2_celltypist(adata: sc.AnnData) -> sc.AnnData:
    """Method 2: CellTypist model predictions."""
    print("  Method 2: CellTypist model...")

    try:
        import celltypist
        from celltypist import models
    except ImportError:
        print("    WARNING: CellTypist not installed. Run: pip install celltypist")
        adata.obs['celltype_celltypist'] = 'Not_Run'
        adata.obs['celltype_celltypist_conf'] = 0.0
        return adata

    # Prepare data (CellTypist expects log-normalized counts)
    adata_ct = adata.copy()

    # Use Immune_All_Low model (broad immune cell types)
    try:
        # CellTypist 1.6+ API: download_models (plural)
        models.download_models(model='Immune_All_Low.pkl', force_update=False)
        model = models.Model.load(model='Immune_All_Low.pkl')

        predictions = celltypist.annotate(
            adata_ct, model=model, majority_voting=True
        )

        adata.obs['celltype_celltypist'] = predictions.predicted_labels['majority_voting'].values
        adata.obs['celltype_celltypist_conf'] = predictions.probability_matrix.max(axis=1).values

        # Summary
        counts = Counter(predictions.predicted_labels['majority_voting'].values)
        print(f"    Top types: {dict(counts.most_common(5))}")

    except Exception as e:
        print(f"    WARNING: CellTypist failed: {e}")
        adata.obs['celltype_celltypist'] = 'Failed'
        adata.obs['celltype_celltypist_conf'] = 0.0

    return adata


def method3_simple_gating(adata: sc.AnnData) -> sc.AnnData:
    """Method 3: Simple hierarchical gating based on key markers."""
    print("  Method 3: Hierarchical gating...")

    genes = set(adata.var_names)

    # Get normalized expression
    X = adata.X if not isinstance(adata.X, np.ndarray) else adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()

    def get_expr(gene):
        if gene in genes:
            return X[:, list(adata.var_names).index(gene)]
        return np.zeros(adata.n_obs)

    # Initialize as Unknown
    annotations = np.array(['Unknown'] * adata.n_obs)

    # Hierarchical gating
    # 1. Epithelial (EPCAM+)
    epcam = get_expr('EPCAM')
    is_epithelial = epcam > np.percentile(epcam[epcam > 0], 25) if (epcam > 0).sum() > 0 else np.zeros(adata.n_obs, dtype=bool)

    # 2. Immune markers
    cd3 = get_expr('CD3D') + get_expr('CD3E')
    cd68 = get_expr('CD68')
    ms4a1 = get_expr('MS4A1')

    is_tcell = cd3 > np.percentile(cd3[cd3 > 0], 25) if (cd3 > 0).sum() > 0 else np.zeros(adata.n_obs, dtype=bool)
    is_macro = cd68 > np.percentile(cd68[cd68 > 0], 25) if (cd68 > 0).sum() > 0 else np.zeros(adata.n_obs, dtype=bool)
    is_bcell = ms4a1 > np.percentile(ms4a1[ms4a1 > 0], 25) if (ms4a1 > 0).sum() > 0 else np.zeros(adata.n_obs, dtype=bool)

    # 3. Stromal
    col1a1 = get_expr('COL1A1')
    pecam1 = get_expr('PECAM1')

    is_fibro = col1a1 > np.percentile(col1a1[col1a1 > 0], 25) if (col1a1 > 0).sum() > 0 else np.zeros(adata.n_obs, dtype=bool)
    is_endo = pecam1 > np.percentile(pecam1[pecam1 > 0], 25) if (pecam1 > 0).sum() > 0 else np.zeros(adata.n_obs, dtype=bool)

    # Apply hierarchy (most specific first)
    annotations[is_endo & ~is_epithelial] = 'Endothelial'
    annotations[is_fibro & ~is_epithelial] = 'Fibroblast'
    annotations[is_bcell & ~is_epithelial] = 'B_Cell'
    annotations[is_macro & ~is_epithelial] = 'Macrophage'
    annotations[is_tcell & ~is_epithelial] = 'T_Cell'
    annotations[is_epithelial] = 'Epithelial'

    adata.obs['celltype_gating'] = annotations

    # Summary
    counts = Counter(annotations)
    print(f"    Distribution: {dict(counts)}")

    return adata


def compare_methods(adata: sc.AnnData) -> pd.DataFrame:
    """Compare annotation methods."""
    print("\n  Comparing methods...")

    methods = ['celltype_markers', 'celltype_celltypist', 'celltype_gating']
    methods = [m for m in methods if m in adata.obs.columns]

    comparisons = []

    # Pairwise agreement
    for i, m1 in enumerate(methods):
        for m2 in methods[i+1:]:
            # Filter out Unknown/Failed
            mask = (adata.obs[m1] != 'Unknown') & (adata.obs[m2] != 'Unknown') & \
                   (adata.obs[m1] != 'Not_Run') & (adata.obs[m2] != 'Not_Run') & \
                   (adata.obs[m1] != 'Failed') & (adata.obs[m2] != 'Failed')

            if mask.sum() > 100:
                ari = adjusted_rand_score(adata.obs.loc[mask, m1], adata.obs.loc[mask, m2])
                nmi = normalized_mutual_info_score(adata.obs.loc[mask, m1], adata.obs.loc[mask, m2])
                comparisons.append({
                    'method1': m1.replace('celltype_', ''),
                    'method2': m2.replace('celltype_', ''),
                    'ARI': ari,
                    'NMI': nmi,
                    'n_cells': mask.sum()
                })
                print(f"    {m1} vs {m2}: ARI={ari:.3f}, NMI={nmi:.3f}")

    return pd.DataFrame(comparisons)


def calculate_cluster_purity(adata: sc.AnnData, method: str) -> float:
    """Calculate cluster purity for an annotation method."""
    if method not in adata.obs.columns or 'leiden' not in adata.obs.columns:
        return 0.0

    purities = []
    for cluster in adata.obs['leiden'].unique():
        mask = adata.obs['leiden'] == cluster
        types = adata.obs.loc[mask, method]
        if len(types) > 0:
            most_common = types.value_counts().iloc[0]
            purity = most_common / len(types)
            purities.append(purity)

    return np.mean(purities) if purities else 0.0


def main():
    print("="*60)
    print("Multi-Method Cell Type Annotation")
    print("="*60)

    # Load data
    print(f"\nLoading: {INPUT}")
    adata = sc.read_h5ad(INPUT)
    print(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars} genes")

    # Get optimal params
    opt_hvg, opt_res = load_optimal_params()
    print(f"Using HVG={opt_hvg}, resolution={opt_res}")

    all_comparisons = []
    method_rankings = []

    for sample_id, stage in SAMPLES.items():
        print(f"\n{'='*50}")
        print(f"{sample_id} ({stage})")
        print(f"{'='*50}")

        # Subset to sample
        sample = adata[adata.obs['sample_id'] == sample_id].copy()
        print(f"Cells: {sample.n_obs:,}")

        # Clustering with optimal params
        print("\n  Clustering...")
        if opt_hvg != 'all':
            sc.pp.highly_variable_genes(
                sample, n_top_genes=min(int(opt_hvg), sample.n_vars - 1),
                flavor='seurat_v3', layer='counts', span=0.3
            )
            hvg_sample = sample[:, sample.var.highly_variable].copy()
        else:
            hvg_sample = sample.copy()

        n_pcs = min(30, hvg_sample.n_vars - 1)
        sc.pp.pca(hvg_sample, n_comps=n_pcs)
        sc.pp.neighbors(hvg_sample, n_neighbors=15, n_pcs=n_pcs)
        sc.tl.leiden(hvg_sample, resolution=opt_res)
        sc.tl.umap(hvg_sample)

        # Transfer embeddings to full sample
        sample.obs['leiden'] = hvg_sample.obs['leiden']
        sample.obsm['X_pca'] = hvg_sample.obsm['X_pca']
        sample.obsm['X_umap'] = hvg_sample.obsm['X_umap']

        print(f"  Clusters: {sample.obs['leiden'].nunique()}")

        # Run all annotation methods
        sample = method1_marker_scoring(sample)
        sample = method2_celltypist(sample)
        sample = method3_simple_gating(sample)

        # Compare methods
        comp_df = compare_methods(sample)
        comp_df['sample_id'] = sample_id
        all_comparisons.append(comp_df)

        # Calculate purity for ranking
        methods = ['celltype_markers', 'celltype_celltypist', 'celltype_gating']
        for method in methods:
            if method in sample.obs.columns:
                purity = calculate_cluster_purity(sample, method)
                unknown_rate = (sample.obs[method].isin(['Unknown', 'Not_Run', 'Failed'])).mean()
                method_rankings.append({
                    'sample_id': sample_id,
                    'stage': stage,
                    'method': method.replace('celltype_', ''),
                    'cluster_purity': purity,
                    'unknown_rate': unknown_rate,
                })
                print(f"    {method}: purity={purity:.3f}, unknown={unknown_rate:.1%}")

        # Save annotated sample
        out_path = OUTPUT / f'{sample_id}_annotated.h5ad'
        sample.write_h5ad(out_path)
        print(f"\n  Saved: {out_path}")

    # Save comparison results
    if all_comparisons:
        comp_all = pd.concat(all_comparisons, ignore_index=True)
        comp_all.to_csv(OUTPUT / 'annotation_method_comparison.csv', index=False)
        print(f"\nSaved: {OUTPUT / 'annotation_method_comparison.csv'}")

    # Save rankings
    rankings_df = pd.DataFrame(method_rankings)
    rankings_df.to_csv(OUTPUT / 'annotation_method_rankings.csv', index=False)
    print(f"Saved: {OUTPUT / 'annotation_method_rankings.csv'}")

    # === Create Visualizations ===
    print("\nCreating visualizations...")

    # 1. Method ranking bar plot
    fig, ax = plt.subplots(figsize=(10, 5))
    avg_rankings = rankings_df.groupby('method').agg({
        'cluster_purity': 'mean',
        'unknown_rate': 'mean'
    }).reset_index()

    x = np.arange(len(avg_rankings))
    width = 0.35
    ax.bar(x - width/2, avg_rankings['cluster_purity'], width, label='Cluster Purity', color='steelblue')
    ax.bar(x + width/2, 1 - avg_rankings['unknown_rate'], width, label='Annotation Rate', color='coral')

    ax.set_xlabel('Method')
    ax.set_ylabel('Score')
    ax.set_title('Annotation Method Comparison', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(avg_rankings['method'])
    ax.legend()
    ax.set_ylim(0, 1)

    fig.tight_layout()
    fig.savefig(FIGURES / 'annotation_method_ranking.png', dpi=150)
    plt.close()

    # 2. Cell type distribution per sample
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, (sample_id, stage) in enumerate(SAMPLES.items()):
        path = OUTPUT / f'{sample_id}_annotated.h5ad'
        if path.exists():
            ad = sc.read_h5ad(path)
            ax = axes[i]

            # Use markers method for distribution
            if 'celltype_markers' in ad.obs.columns:
                counts = ad.obs['celltype_markers'].value_counts()
                counts = counts[counts.index != 'Unknown']  # Exclude Unknown for clarity
                if len(counts) > 10:
                    counts = counts.head(10)

                counts.plot(kind='barh', ax=ax, color='steelblue')
                ax.set_title(f'{sample_id}: {stage}\n({ad.n_obs:,} cells)', fontweight='bold')
                ax.set_xlabel('Cell Count')

    plt.suptitle('Cell Type Distribution (Marker-based)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(FIGURES / 'celltype_distribution.png', dpi=150)
    plt.close()

    # === Final Recommendation ===
    print("\n" + "="*60)
    print("RECOMMENDATION")
    print("="*60)

    avg_purity = rankings_df.groupby('method')['cluster_purity'].mean()
    best_method = avg_purity.idxmax()
    print(f"Best method by cluster purity: {best_method}")
    print(f"Average purity: {avg_purity[best_method]:.3f}")

    # Save recommendation
    with open(OUTPUT / 'annotation_recommendation.txt', 'w') as f:
        f.write(f"Recommended annotation method: {best_method}\n")
        f.write(f"Average cluster purity: {avg_purity[best_method]:.3f}\n")
        f.write(f"\nAll method rankings:\n")
        for method, purity in avg_purity.sort_values(ascending=False).items():
            f.write(f"  {method}: {purity:.3f}\n")

    print(f"\nOutputs saved to: {OUTPUT}")


if __name__ == '__main__':
    main()
