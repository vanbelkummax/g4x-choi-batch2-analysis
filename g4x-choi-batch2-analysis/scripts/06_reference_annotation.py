#!/usr/bin/env python3
"""
06_reference_annotation.py - Reference-based cell type annotation

Downloads Kumar et al. GSE183904 gastric cancer atlas and uses KNN-based
label transfer to annotate G4X cells. Works with transformed data.

Reference: "Single-Cell Atlas of Lineage States, Tumor Microenvironment,
and Subtype-Specific Expression Programs in Gastric Cancer"
Cancer Discovery 2022 (PMID: 34642171)
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.sparse import issparse
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import LabelEncoder
from collections import Counter
import warnings
import gzip
import tarfile
import os
import urllib.request

warnings.filterwarnings('ignore')

# === CONFIG ===
BASE = Path('/home/user/g4x-choi-batch2-analysis')
RESULTS = BASE / 'results/pilot'
FIGURES = RESULTS / 'figures'
REFERENCE_DIR = BASE / 'reference'

SAMPLES = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}

# GEO accession for Kumar et al. 2022 gastric cancer atlas
GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183904/suppl/"
GEO_FILES = [
    "GSE183904_matrix.mtx.gz",
    "GSE183904_barcodes.tsv.gz",
    "GSE183904_genes.tsv.gz",
    "GSE183904_metadata.csv.gz"
]

# KNN parameters
N_NEIGHBORS = 15
N_PCS = 30


def download_reference():
    """Download Kumar et al. 2022 reference dataset from GEO."""
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

    # Check if already downloaded
    ref_h5ad = REFERENCE_DIR / 'kumar_gastric_reference.h5ad'
    if ref_h5ad.exists():
        print(f"  Reference already exists: {ref_h5ad}")
        return ref_h5ad

    print("  Downloading Kumar et al. GSE183904 reference...")

    # Try to download supplementary files
    downloaded = []
    for fname in GEO_FILES:
        fpath = REFERENCE_DIR / fname
        if not fpath.exists():
            url = GEO_URL + fname
            print(f"    Downloading {fname}...")
            try:
                urllib.request.urlretrieve(url, fpath)
                downloaded.append(fpath)
            except Exception as e:
                print(f"    WARNING: Could not download {fname}: {e}")

    # Check if we got all files
    mtx_path = REFERENCE_DIR / 'GSE183904_matrix.mtx.gz'
    barcodes_path = REFERENCE_DIR / 'GSE183904_barcodes.tsv.gz'
    genes_path = REFERENCE_DIR / 'GSE183904_genes.tsv.gz'
    meta_path = REFERENCE_DIR / 'GSE183904_metadata.csv.gz'

    if not all(p.exists() for p in [mtx_path, barcodes_path, genes_path, meta_path]):
        print("  WARNING: Could not download all reference files")
        print("  Falling back to synthetic reference based on marker signatures")
        return create_marker_reference()

    # Load as AnnData
    print("  Loading reference data...")
    from scipy.io import mmread

    # Read MTX
    X = mmread(gzip.open(mtx_path, 'rb')).T.tocsr()

    # Read barcodes and genes
    with gzip.open(barcodes_path, 'rt') as f:
        barcodes = [line.strip() for line in f]
    with gzip.open(genes_path, 'rt') as f:
        genes = [line.strip().split('\t')[0] for line in f]

    # Read metadata
    meta = pd.read_csv(meta_path, compression='gzip', index_col=0)

    # Create AnnData
    ref = sc.AnnData(X=X)
    ref.obs_names = barcodes[:X.shape[0]]
    ref.var_names = genes[:X.shape[1]]

    # Add metadata
    ref.obs = ref.obs.join(meta, how='left')

    # Preprocess
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)

    # Save
    ref.write_h5ad(ref_h5ad)
    print(f"  Saved reference: {ref_h5ad}")

    return ref_h5ad


def create_marker_reference():
    """Create a marker-based pseudo-reference when GEO download fails."""
    ref_h5ad = REFERENCE_DIR / 'marker_reference.h5ad'

    if ref_h5ad.exists():
        return ref_h5ad

    print("  Creating marker-based pseudo-reference...")

    # Cell types and their canonical markers
    # Based on Kumar et al. 2022 and general gastric cancer literature
    CELL_TYPE_MARKERS = {
        # Epithelial
        'Gastric_Pit_Mucous': ['MUC5AC', 'MUC6', 'TFF1', 'TFF2', 'EPCAM'],
        'Gastric_Chief_Parietal': ['PGC', 'ATP4A', 'EPCAM'],
        'Goblet': ['MUC2', 'TFF3', 'EPCAM', 'SPINK4'],
        'Intestinal_Enterocyte': ['CDX1', 'CDX2', 'CDH17', 'EPCAM'],
        'Enteroendocrine': ['GHRL', 'SST', 'GAST', 'CHGA'],
        'Stem_Progenitor': ['LGR5', 'PROM1', 'OLFM4', 'SOX9'],
        'Proliferating_Epithelial': ['MKI67', 'TOP2A', 'EPCAM'],

        # T cells
        'T_CD4_Naive': ['CD3D', 'CD3E', 'CD4', 'CCR7', 'IL7R'],
        'T_CD4_Memory': ['CD3D', 'CD3E', 'CD4', 'IL7R'],
        'T_CD8_Cytotoxic': ['CD3D', 'CD3E', 'CD8A', 'GZMA', 'GZMB', 'PRF1'],
        'T_CD8_Exhausted': ['CD3D', 'CD8A', 'PDCD1', 'CTLA4', 'LAG3', 'HAVCR2'],
        'T_Reg': ['CD3D', 'CD4', 'FOXP3', 'IL2RA'],
        'T_NK': ['NKG7', 'GNLY', 'KLRD1', 'CD3D'],

        # B cells
        'B_Cell': ['MS4A1', 'CD19', 'CD79A', 'CD79B'],
        'Plasma_Cell': ['IGHA1', 'IGHG1', 'IGHM', 'JCHAIN', 'SDC1'],

        # Myeloid
        'Macrophage_M1': ['CD68', 'CD163', 'CSF1R', 'IL1B', 'TNF'],
        'Macrophage_M2': ['CD68', 'CD163', 'MRC1', 'CD206'],
        'Monocyte': ['CD14', 'ITGAM', 'S100A9'],
        'Dendritic_Cell': ['CD1C', 'CLEC10A', 'HLA-DRA'],
        'Mast_Cell': ['KIT', 'TPSAB1', 'CPA3'],

        # Stromal
        'Fibroblast': ['COL1A1', 'LUM', 'VCAN', 'DCN'],
        'CAF_Myofibroblast': ['ACTA2', 'TAGLN', 'MYH11'],
        'CAF_Inflammatory': ['FAP', 'PDPN', 'THY1', 'POSTN'],
        'Endothelial': ['PECAM1', 'VWF', 'CDH5'],
        'Pericyte': ['PDGFRB', 'RGS5', 'ACTA2'],
    }

    # Get all unique markers
    all_markers = set()
    for markers in CELL_TYPE_MARKERS.values():
        all_markers.update(markers)
    all_markers = sorted(all_markers)

    # Create pseudo-expression profiles (centroids)
    n_types = len(CELL_TYPE_MARKERS)
    n_genes = len(all_markers)

    # Create expression matrix
    X = np.zeros((n_types, n_genes))
    cell_types = list(CELL_TYPE_MARKERS.keys())

    for i, (ct, markers) in enumerate(CELL_TYPE_MARKERS.items()):
        for marker in markers:
            if marker in all_markers:
                j = all_markers.index(marker)
                X[i, j] = 2.0 + np.random.uniform(0.5, 1.5)  # High expression

    # Add some background noise
    X += np.random.uniform(0, 0.1, X.shape)

    # Create AnnData
    ref = sc.AnnData(X=X)
    ref.var_names = all_markers
    ref.obs_names = cell_types
    ref.obs['cell_type'] = cell_types
    ref.obs['cell_type_coarse'] = [ct.split('_')[0] for ct in cell_types]

    ref.write_h5ad(ref_h5ad)
    print(f"  Created marker reference with {n_types} cell types, {n_genes} markers")

    return ref_h5ad


def load_and_prepare_reference(ref_path: Path):
    """Load reference and prepare for label transfer."""
    print(f"  Loading reference: {ref_path}")
    ref = sc.read_h5ad(ref_path)

    # Check for cell type column
    type_col = None
    for col in ['cell_type', 'celltype', 'Cell_type', 'cluster', 'annotation']:
        if col in ref.obs.columns:
            type_col = col
            break

    if type_col is None:
        print(f"    WARNING: No cell type column found in reference")
        print(f"    Available columns: {list(ref.obs.columns)}")
        return None, None

    print(f"    Cell type column: {type_col}")
    print(f"    Reference shape: {ref.shape}")
    print(f"    Cell types: {ref.obs[type_col].nunique()}")

    return ref, type_col


def find_common_genes(query, reference):
    """Find genes common to both datasets."""
    query_genes = set(query.var_names)
    ref_genes = set(reference.var_names)
    common = sorted(query_genes & ref_genes)

    print(f"    Query genes: {len(query_genes)}")
    print(f"    Reference genes: {len(ref_genes)}")
    print(f"    Common genes: {len(common)}")

    return common


def knn_label_transfer(query, reference, type_col, n_neighbors=15, n_pcs=30):
    """Transfer labels from reference to query using KNN in PCA space."""
    print(f"  Performing KNN label transfer (k={n_neighbors})...")

    # Find common genes
    common_genes = find_common_genes(query, reference)

    if len(common_genes) < 50:
        print(f"    WARNING: Only {len(common_genes)} common genes - too few for reliable transfer")
        return None, None

    # Subset to common genes
    query_sub = query[:, common_genes].copy()
    ref_sub = reference[:, common_genes].copy()

    # Get expression matrices
    X_query = query_sub.X.toarray() if issparse(query_sub.X) else query_sub.X
    X_ref = ref_sub.X.toarray() if issparse(ref_sub.X) else ref_sub.X

    # Joint PCA
    print("    Computing joint PCA...")
    X_combined = np.vstack([X_ref, X_query])

    from sklearn.decomposition import PCA
    n_comps = min(n_pcs, min(X_combined.shape) - 1, len(common_genes) - 1)
    pca = PCA(n_components=n_comps)
    pca.fit(X_combined)

    X_ref_pca = pca.transform(X_ref)
    X_query_pca = pca.transform(X_query)

    # Train KNN on reference
    print(f"    Training KNN classifier...")
    labels = reference.obs[type_col].values
    le = LabelEncoder()
    labels_encoded = le.fit_transform(labels)

    knn = KNeighborsClassifier(n_neighbors=min(n_neighbors, len(labels)), weights='distance')
    knn.fit(X_ref_pca, labels_encoded)

    # Predict on query
    print(f"    Predicting labels for {query.n_obs:,} cells...")
    pred_encoded = knn.predict(X_query_pca)
    pred_labels = le.inverse_transform(pred_encoded)

    # Get prediction confidence (proportion of neighbors with same label)
    proba = knn.predict_proba(X_query_pca)
    confidence = np.max(proba, axis=1)

    print(f"    Mean confidence: {confidence.mean():.3f}")
    print(f"    Predictions: {Counter(pred_labels).most_common(5)}")

    return pred_labels, confidence


def annotate_sample(sample_id: str, ref, type_col):
    """Annotate a single sample with reference labels."""
    input_path = RESULTS / f'{sample_id}_annotated.h5ad'

    if not input_path.exists():
        print(f"  WARNING: {input_path} not found")
        return None

    print(f"\n  Annotating {sample_id}...")
    adata = sc.read_h5ad(input_path)
    print(f"    Cells: {adata.n_obs:,}")

    # Perform label transfer
    labels, confidence = knn_label_transfer(
        adata, ref, type_col,
        n_neighbors=N_NEIGHBORS,
        n_pcs=N_PCS
    )

    if labels is None:
        adata.obs['celltype_reference'] = 'Transfer_Failed'
        adata.obs['celltype_reference_conf'] = 0.0
    else:
        adata.obs['celltype_reference'] = labels
        adata.obs['celltype_reference_conf'] = confidence

        # Create coarse labels
        coarse = [ct.split('_')[0] if '_' in ct else ct for ct in labels]
        adata.obs['celltype_reference_coarse'] = coarse

    # Save
    output_path = RESULTS / f'{sample_id}_reference_annotated.h5ad'
    adata.write_h5ad(output_path)
    print(f"    Saved: {output_path}")

    return adata


def compare_annotations(adata, sample_id, stage):
    """Compare marker-based and reference-based annotations."""
    results = []

    if 'celltype_markers' not in adata.obs.columns:
        return pd.DataFrame()
    if 'celltype_reference' not in adata.obs.columns:
        return pd.DataFrame()

    # Get proportions for each method
    marker_counts = adata.obs['celltype_markers'].value_counts(normalize=True)
    ref_counts = adata.obs['celltype_reference'].value_counts(normalize=True)

    # Goblet comparison (key validation)
    goblet_marker = marker_counts.get('Goblet', 0)
    goblet_ref = ref_counts.get('Goblet', 0)

    results.append({
        'sample_id': sample_id,
        'stage': stage,
        'cell_type': 'Goblet',
        'marker_pct': goblet_marker * 100,
        'reference_pct': goblet_ref * 100,
        'difference': abs(goblet_marker - goblet_ref) * 100
    })

    # All cell types comparison
    all_types = set(marker_counts.index) | set(ref_counts.index)
    for ct in all_types:
        if ct in ['Unknown', 'Transfer_Failed']:
            continue
        results.append({
            'sample_id': sample_id,
            'stage': stage,
            'cell_type': ct,
            'marker_pct': marker_counts.get(ct, 0) * 100,
            'reference_pct': ref_counts.get(ct, 0) * 100,
            'difference': abs(marker_counts.get(ct, 0) - ref_counts.get(ct, 0)) * 100
        })

    return pd.DataFrame(results)


def create_visualizations(annotated_samples: dict):
    """Create comparison visualizations."""
    FIGURES.mkdir(parents=True, exist_ok=True)

    # 1. Reference annotation UMAP per sample
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for i, (sample_id, (adata, stage)) in enumerate(annotated_samples.items()):
        ax = axes[i]

        if 'celltype_reference' in adata.obs.columns and 'X_umap' in adata.obsm:
            # Get top 10 cell types for color
            top_types = adata.obs['celltype_reference'].value_counts().head(10).index.tolist()
            colors = adata.obs['celltype_reference'].apply(
                lambda x: x if x in top_types else 'Other'
            )

            sc.pl.umap(adata, color='celltype_reference', ax=ax, show=False,
                      legend_loc='on data', legend_fontsize=6,
                      title=f'{sample_id}: {stage}')

    plt.suptitle('Reference-Based Cell Type Annotation', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(FIGURES / 'reference_annotation_umap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {FIGURES / 'reference_annotation_umap.png'}")

    # 2. Marker vs Reference comparison bar chart
    comparison_data = []
    for sample_id, (adata, stage) in annotated_samples.items():
        comp_df = compare_annotations(adata, sample_id, stage)
        comparison_data.append(comp_df)

    if comparison_data:
        comp_all = pd.concat(comparison_data, ignore_index=True)
        comp_all.to_csv(RESULTS / 'reference_annotations.csv', index=False)

        # Goblet validation plot
        goblet_df = comp_all[comp_all['cell_type'] == 'Goblet']
        if len(goblet_df) > 0:
            fig, ax = plt.subplots(figsize=(10, 6))

            x = np.arange(len(goblet_df))
            width = 0.35

            ax.bar(x - width/2, goblet_df['marker_pct'], width,
                   label='Marker-based', color='steelblue')
            ax.bar(x + width/2, goblet_df['reference_pct'], width,
                   label='Reference-based', color='coral')

            ax.set_xlabel('Sample')
            ax.set_ylabel('Goblet Cell %')
            ax.set_title('Goblet Cell Proportion: Marker vs Reference Validation',
                        fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels([f"{row['sample_id']}\n({row['stage']})"
                               for _, row in goblet_df.iterrows()])
            ax.legend()

            # Add values on bars
            for i, (_, row) in enumerate(goblet_df.iterrows()):
                ax.annotate(f"{row['marker_pct']:.1f}%",
                           xy=(i - width/2, row['marker_pct']),
                           ha='center', va='bottom', fontsize=9)
                ax.annotate(f"{row['reference_pct']:.1f}%",
                           xy=(i + width/2, row['reference_pct']),
                           ha='center', va='bottom', fontsize=9)

            plt.tight_layout()
            fig.savefig(FIGURES / 'goblet_validation.png', dpi=150)
            plt.close()
            print(f"  Saved: {FIGURES / 'goblet_validation.png'}")

    # 3. Confidence distribution
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    for i, (sample_id, (adata, stage)) in enumerate(annotated_samples.items()):
        ax = axes[i]

        if 'celltype_reference_conf' in adata.obs.columns:
            conf = adata.obs['celltype_reference_conf']
            ax.hist(conf, bins=50, color='steelblue', edgecolor='white')
            ax.axvline(conf.mean(), color='red', linestyle='--',
                      label=f'Mean: {conf.mean():.2f}')
            ax.set_xlabel('Prediction Confidence')
            ax.set_ylabel('Cell Count')
            ax.set_title(f'{sample_id}: {stage}')
            ax.legend()

    plt.suptitle('Reference Label Transfer Confidence', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(FIGURES / 'reference_confidence.png', dpi=150)
    plt.close()
    print(f"  Saved: {FIGURES / 'reference_confidence.png'}")


def main():
    print("=" * 60)
    print("Reference-Based Cell Type Annotation")
    print("Kumar et al. 2022 Gastric Cancer Atlas")
    print("=" * 60)

    # Step 1: Download/create reference
    print("\n[1/3] Preparing reference dataset...")
    ref_path = download_reference()

    if ref_path is None:
        print("ERROR: Could not prepare reference dataset")
        return

    # Load reference
    ref, type_col = load_and_prepare_reference(ref_path)

    if ref is None:
        print("ERROR: Could not load reference")
        return

    # Step 2: Annotate each sample
    print("\n[2/3] Annotating samples...")
    annotated = {}

    for sample_id, stage in SAMPLES.items():
        adata = annotate_sample(sample_id, ref, type_col)
        if adata is not None:
            annotated[sample_id] = (adata, stage)

    # Step 3: Create visualizations
    print("\n[3/3] Creating visualizations...")
    create_visualizations(annotated)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    for sample_id, (adata, stage) in annotated.items():
        if 'celltype_reference' in adata.obs.columns:
            counts = adata.obs['celltype_reference'].value_counts()
            print(f"\n{sample_id} ({stage}):")
            print(f"  Top 5 cell types:")
            for ct, n in counts.head(5).items():
                pct = n / adata.n_obs * 100
                print(f"    {ct}: {n:,} ({pct:.1f}%)")

    print(f"\nOutputs saved to: {RESULTS}")


if __name__ == '__main__':
    main()
