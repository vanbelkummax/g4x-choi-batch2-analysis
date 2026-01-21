#!/usr/bin/env python
"""
35_validate_vs_resolve.py - Compare our pipeline outputs against Resolve Biosciences official outputs

Compares:
1. Cell counts - Resolve core_metrics.csv vs our h5ad
2. Clustering - Resolve clustering_umap.csv.gz vs our Leiden clusters
3. Cell metadata - Resolve cell_metadata.csv.gz vs our h5ad obs
4. QC statistics - Compare filtering rates

Author: Max Van Belkum
Date: 2026-01-21
"""

import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path
import gzip
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# Configuration
# ============================================================================
RESOLVE_BASE = Path("/mnt/x/Choi_Batch_2_Tuesday")
RESOLVE_LANE = "g4-028-083-FC1-L001_5XGaAe5DB2dm7sRK"
OUR_BASE = Path.home() / "spatial-hackathon-2026/results/g4x_choi_batch2"

SAMPLES = ["A01", "B01", "C01", "D01", "E01", "F01", "G01", "H01"]

OUTPUT_DIR = OUR_BASE / "validation"
OUTPUT_DIR.mkdir(exist_ok=True)


def load_resolve_core_metrics():
    """Load Resolve's aggregated core metrics."""
    metrics_file = RESOLVE_BASE / "choi_preGC_b2_core_metrics.csv"
    df = pd.read_csv(metrics_file)
    df = df.set_index('sample_id')
    return df


def load_resolve_sample_metrics(sample: str):
    """Load Resolve's per-sample metrics."""
    sample_dir = RESOLVE_BASE / RESOLVE_LANE / sample
    metrics_file = sample_dir / "metrics/core_metrics.csv"
    return pd.read_csv(metrics_file)


def load_resolve_clustering(sample: str):
    """Load Resolve's clustering results."""
    sample_dir = RESOLVE_BASE / RESOLVE_LANE / sample
    cluster_file = sample_dir / "single_cell_data/clustering_umap.csv.gz"
    return pd.read_csv(cluster_file)


def load_resolve_cell_metadata(sample: str):
    """Load Resolve's cell metadata."""
    sample_dir = RESOLVE_BASE / RESOLVE_LANE / sample
    meta_file = sample_dir / "single_cell_data/cell_metadata.csv.gz"
    return pd.read_csv(meta_file)


def load_our_data(sample: str):
    """Load our processed h5ad files."""
    data = {}

    # WNN integrated
    wnn_file = OUR_BASE / "wnn_integrated" / f"{sample}_wnn.h5ad"
    if wnn_file.exists():
        data['wnn'] = sc.read_h5ad(wnn_file)

    # Segmentation QC
    qc_file = OUR_BASE / "segmentation_qc" / f"{sample}_qc.h5ad"
    if qc_file.exists():
        data['qc'] = sc.read_h5ad(qc_file)

    # Annotated
    annot_file = OUR_BASE / "annotated_v2" / f"{sample}_annotated.h5ad"
    if annot_file.exists():
        data['annotated'] = sc.read_h5ad(annot_file)

    return data


def compare_cell_counts():
    """Compare cell counts between Resolve and our pipeline."""
    print("\n" + "="*80)
    print("1. CELL COUNT COMPARISON")
    print("="*80)

    resolve_metrics = load_resolve_core_metrics()

    results = []
    for sample in SAMPLES:
        resolve_n = resolve_metrics.loc[sample, 'number_cells']

        our_data = load_our_data(sample)
        our_n = our_data['wnn'].n_obs if 'wnn' in our_data else 0

        diff = resolve_n - our_n
        pct_diff = (diff / resolve_n) * 100

        results.append({
            'sample': sample,
            'resolve_cells': resolve_n,
            'our_cells': our_n,
            'difference': diff,
            'pct_filtered': pct_diff
        })

        print(f"{sample}: Resolve={resolve_n:,} | Ours={our_n:,} | Filtered={diff:,} ({pct_diff:.1f}%)")

    df = pd.DataFrame(results)
    print(f"\nTotal: Resolve={df['resolve_cells'].sum():,} | Ours={df['our_cells'].sum():,}")
    print(f"Overall filtering: {(1 - df['our_cells'].sum()/df['resolve_cells'].sum())*100:.1f}%")

    return df


def compare_clustering():
    """Compare clustering assignments between Resolve and our pipeline."""
    print("\n" + "="*80)
    print("2. CLUSTERING COMPARISON")
    print("="*80)

    results = []
    for sample in SAMPLES:
        # Load Resolve clustering
        resolve_clust = load_resolve_clustering(sample)
        resolve_clusters = resolve_clust['leiden_0.400'].value_counts()

        # Load our data
        our_data = load_our_data(sample)
        if 'wnn' not in our_data:
            continue
        adata = our_data['wnn']

        # Check for matching cells
        resolve_labels = set(resolve_clust['label'].values)
        our_labels = set(adata.obs_names)
        common = resolve_labels & our_labels

        # Get our clustering
        our_cluster_col = 'leiden_wnn_0_5' if 'leiden_wnn_0_5' in adata.obs.columns else 'leiden_0.5'
        our_clusters = adata.obs[our_cluster_col].value_counts()

        results.append({
            'sample': sample,
            'resolve_n_clusters': resolve_clusters.nunique(),
            'our_n_clusters': our_clusters.nunique(),
            'common_cells': len(common),
            'resolve_cells': len(resolve_labels),
            'our_cells': len(our_labels)
        })

        print(f"{sample}: Resolve clusters={len(resolve_clusters)} | Our clusters={len(our_clusters)} | Common cells={len(common):,}")

    return pd.DataFrame(results)


def compare_cluster_concordance(sample: str):
    """Compare cluster assignments for shared cells."""
    resolve_clust = load_resolve_clustering(sample)
    resolve_clust = resolve_clust.set_index('label')

    our_data = load_our_data(sample)
    if 'wnn' not in our_data:
        return None
    adata = our_data['wnn']

    # Find common cells
    common_cells = list(set(resolve_clust.index) & set(adata.obs_names))
    if len(common_cells) == 0:
        return None

    # Get cluster assignments
    resolve_assign = resolve_clust.loc[common_cells, 'leiden_0.400']
    our_cluster_col = 'leiden_wnn_0_5' if 'leiden_wnn_0_5' in adata.obs.columns else 'leiden_0.5'
    our_assign = adata.obs.loc[common_cells, our_cluster_col]

    # Create contingency table
    contingency = pd.crosstab(resolve_assign, our_assign)

    return {
        'n_common': len(common_cells),
        'contingency': contingency,
        'resolve_clusters': resolve_assign.nunique(),
        'our_clusters': our_assign.nunique()
    }


def compare_qc_stats():
    """Compare QC statistics."""
    print("\n" + "="*80)
    print("3. QC STATISTICS COMPARISON")
    print("="*80)

    resolve_metrics = load_resolve_core_metrics()

    results = []
    for sample in SAMPLES:
        resolve_row = resolve_metrics.loc[sample]

        # Resolve stats
        resolve_median_tx = resolve_row['median_transcripts_per_cell']
        resolve_median_genes = resolve_row['median_unique_genes_per_cell']
        resolve_pct_empty = resolve_row['pct_empty_cells']

        # Our stats
        our_data = load_our_data(sample)
        if 'wnn' not in our_data:
            continue
        adata = our_data['wnn']

        our_median_tx = adata.obs['total_counts'].median()
        our_median_genes = adata.obs['n_genes'].median()

        results.append({
            'sample': sample,
            'resolve_median_tx': resolve_median_tx,
            'our_median_tx': our_median_tx,
            'resolve_median_genes': resolve_median_genes,
            'our_median_genes': our_median_genes,
            'resolve_pct_empty': resolve_pct_empty
        })

        print(f"{sample}:")
        print(f"  Median transcripts: Resolve={resolve_median_tx:.0f} | Ours={our_median_tx:.0f}")
        print(f"  Median genes: Resolve={resolve_median_genes:.0f} | Ours={our_median_genes:.0f}")

    return pd.DataFrame(results)


def compare_cell_metadata(sample: str):
    """Compare cell-level metadata for a sample."""
    resolve_meta = load_resolve_cell_metadata(sample)
    resolve_meta = resolve_meta.set_index('label')

    our_data = load_our_data(sample)
    if 'wnn' not in our_data:
        return None
    adata = our_data['wnn']

    # Find common cells
    common_cells = list(set(resolve_meta.index) & set(adata.obs_names))

    # Compare spatial coordinates
    resolve_x = resolve_meta.loc[common_cells, 'cell_x']
    resolve_y = resolve_meta.loc[common_cells, 'cell_y']
    our_x = adata.obs.loc[common_cells, 'cell_x']
    our_y = adata.obs.loc[common_cells, 'cell_y']

    coord_corr_x = np.corrcoef(resolve_x, our_x)[0, 1]
    coord_corr_y = np.corrcoef(resolve_y, our_y)[0, 1]

    # Compare nuclei area
    resolve_area = resolve_meta.loc[common_cells, 'nuclei_area']
    our_area = adata.obs.loc[common_cells, 'nuclei_area']
    area_corr = np.corrcoef(resolve_area, our_area)[0, 1]

    return {
        'n_common': len(common_cells),
        'coord_corr_x': coord_corr_x,
        'coord_corr_y': coord_corr_y,
        'area_corr': area_corr
    }


def compare_annotations():
    """Compare cell type annotations (where available)."""
    print("\n" + "="*80)
    print("4. ANNOTATION COMPARISON")
    print("="*80)

    results = []
    for sample in SAMPLES:
        our_data = load_our_data(sample)
        if 'annotated' not in our_data:
            continue
        adata = our_data['annotated']

        # Our annotations
        if 'lineage' in adata.obs.columns:
            lineage_counts = adata.obs['lineage'].value_counts(normalize=True) * 100

            results.append({
                'sample': sample,
                'n_cells': adata.n_obs,
                'pct_epithelial': lineage_counts.get('Epithelial', 0),
                'pct_immune': lineage_counts.get('Immune', 0),
                'pct_stromal': lineage_counts.get('Stromal', 0),
                'pct_endothelial': lineage_counts.get('Endothelial', 0)
            })

            print(f"{sample}: Epi={lineage_counts.get('Epithelial', 0):.1f}% | Imm={lineage_counts.get('Immune', 0):.1f}% | Str={lineage_counts.get('Stromal', 0):.1f}%")

    return pd.DataFrame(results)


def generate_validation_report():
    """Generate comprehensive validation report."""
    print("\n" + "="*80)
    print("GENERATING VALIDATION REPORT")
    print("="*80)

    # Run all comparisons
    cell_counts = compare_cell_counts()
    clustering = compare_clustering()
    qc_stats = compare_qc_stats()
    annotations = compare_annotations()

    # Cell metadata correlations
    print("\n" + "="*80)
    print("5. CELL METADATA CORRELATIONS")
    print("="*80)

    meta_results = []
    for sample in SAMPLES:
        result = compare_cell_metadata(sample)
        if result:
            meta_results.append({
                'sample': sample,
                **result
            })
            print(f"{sample}: coord_corr={result['coord_corr_x']:.4f},{result['coord_corr_y']:.4f} | area_corr={result['area_corr']:.4f}")

    meta_df = pd.DataFrame(meta_results)

    # Save results
    cell_counts.to_csv(OUTPUT_DIR / "cell_count_comparison.csv", index=False)
    clustering.to_csv(OUTPUT_DIR / "clustering_comparison.csv", index=False)
    qc_stats.to_csv(OUTPUT_DIR / "qc_stats_comparison.csv", index=False)
    annotations.to_csv(OUTPUT_DIR / "annotation_summary.csv", index=False)
    meta_df.to_csv(OUTPUT_DIR / "metadata_correlation.csv", index=False)

    # Generate summary
    print("\n" + "="*80)
    print("VALIDATION SUMMARY")
    print("="*80)

    total_resolve = cell_counts['resolve_cells'].sum()
    total_ours = cell_counts['our_cells'].sum()
    overall_filter_pct = (1 - total_ours/total_resolve) * 100

    print(f"""
┌─────────────────────────────────────────────────────────────────────┐
│                    RESOLVE vs OUR PIPELINE                         │
├─────────────────────────────────────────────────────────────────────┤
│ Total Cells:                                                        │
│   Resolve:     {total_resolve:>10,}                                         │
│   Ours:        {total_ours:>10,}                                         │
│   Filtered:    {total_resolve - total_ours:>10,} ({overall_filter_pct:.1f}%)                                │
├─────────────────────────────────────────────────────────────────────┤
│ Median Transcripts/Cell:                                            │
│   Resolve:     {qc_stats['resolve_median_tx'].mean():.1f}                                              │
│   Ours:        {qc_stats['our_median_tx'].mean():.1f}                                              │
├─────────────────────────────────────────────────────────────────────┤
│ Metadata Correlation (spatial coords):                              │
│   X coord:     {meta_df['coord_corr_x'].mean():.4f}                                            │
│   Y coord:     {meta_df['coord_corr_y'].mean():.4f}                                            │
│   Area:        {meta_df['area_corr'].mean():.4f}                                            │
├─────────────────────────────────────────────────────────────────────┤
│ Our Lineage Annotations:                                            │
│   Epithelial:  {annotations['pct_epithelial'].mean():.1f}%                                            │
│   Immune:      {annotations['pct_immune'].mean():.1f}%                                            │
│   Stromal:     {annotations['pct_stromal'].mean():.1f}%                                            │
└─────────────────────────────────────────────────────────────────────┘
""")

    # Save summary
    summary = {
        'total_resolve_cells': total_resolve,
        'total_our_cells': total_ours,
        'overall_filter_pct': overall_filter_pct,
        'mean_coord_corr_x': meta_df['coord_corr_x'].mean(),
        'mean_coord_corr_y': meta_df['coord_corr_y'].mean(),
        'mean_area_corr': meta_df['area_corr'].mean(),
        'mean_pct_epithelial': annotations['pct_epithelial'].mean(),
        'mean_pct_immune': annotations['pct_immune'].mean(),
        'mean_pct_stromal': annotations['pct_stromal'].mean()
    }
    pd.Series(summary).to_csv(OUTPUT_DIR / "validation_summary.csv")

    print(f"\nResults saved to: {OUTPUT_DIR}")

    return {
        'cell_counts': cell_counts,
        'clustering': clustering,
        'qc_stats': qc_stats,
        'annotations': annotations,
        'metadata': meta_df
    }


if __name__ == "__main__":
    results = generate_validation_report()
