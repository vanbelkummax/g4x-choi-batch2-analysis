#!/usr/bin/env python3
"""
06_reference_annotation.py - Reference-based cell type annotation

Builds reference from Kumar et al. GSE183904 Normal gastric samples,
annotates them with markers, then transfers labels to G4X data via KNN.

Reference: "Single-Cell Atlas of Lineage States, Tumor Microenvironment,
and Subtype-Specific Expression Programs in Gastric Cancer"
Cancer Discovery 2022 (PMID: 34642171)
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.sparse import csr_matrix, issparse
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA
from collections import Counter
import warnings
import gzip
from concurrent.futures import ProcessPoolExecutor, as_completed
import gc

warnings.filterwarnings('ignore')

# === CONFIG ===
BASE = Path('/home/user/g4x-choi-batch2-analysis')
RESULTS = BASE / 'results/pilot'
FIGURES = RESULTS / 'figures'
REFERENCE_DIR = BASE / 'reference/GSE183904'

SAMPLES = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}

# Normal samples from Kumar et al. (for building reference)
# GSM5573466_sample1 = Normal, GSM5573469_sample4 = Normal, etc.
NORMAL_SAMPLES = [
    ('GSM5573466_sample1.csv.gz', 'sample1_normal'),
    ('GSM5573469_sample4.csv.gz', 'sample4_normal'),
    ('GSM5573471_sample6.csv.gz', 'sample6_normal'),
    ('GSM5573474_sample9.csv.gz', 'sample9_normal'),
    ('GSM5573476_sample11.csv.gz', 'sample11_normal'),
    ('GSM5573486_sample21.csv.gz', 'sample21_normal'),
    ('GSM5573488_sample23.csv.gz', 'sample23_normal'),
    ('GSM5573490_sample25.csv.gz', 'sample25_normal'),
]

# Gastric-specific markers for reference annotation
GASTRIC_MARKERS = {
    # Epithelial
    'Gastric_Pit_Mucous': ['MUC5AC', 'MUC6', 'TFF1', 'TFF2', 'EPCAM', 'CDH1'],
    'Gastric_Chief_Parietal': ['PGC', 'ATP4A', 'ATP4B', 'GIF', 'EPCAM'],
    'Goblet': ['MUC2', 'TFF3', 'SPINK4', 'EPCAM'],
    'Enteroendocrine': ['CHGA', 'CHGB', 'GHRL', 'SST', 'GAST'],
    'Stem_Progenitor': ['LGR5', 'OLFM4', 'SOX9', 'ASCL2'],
    'Proliferating': ['MKI67', 'TOP2A', 'PCNA'],

    # T cells
    'T_CD4': ['CD3D', 'CD3E', 'CD4', 'IL7R', 'CCR7'],
    'T_CD8': ['CD3D', 'CD3E', 'CD8A', 'CD8B'],
    'T_Cytotoxic': ['GZMA', 'GZMB', 'GZMK', 'PRF1', 'NKG7'],
    'T_Reg': ['FOXP3', 'IL2RA', 'CTLA4'],
    'T_NK': ['KLRD1', 'NCAM1', 'GNLY', 'NKG7'],

    # B cells
    'B_Cell': ['MS4A1', 'CD19', 'CD79A', 'CD79B', 'PAX5'],
    'Plasma_Cell': ['JCHAIN', 'IGHA1', 'IGHG1', 'SDC1', 'MZB1'],

    # Myeloid
    'Macrophage': ['CD68', 'CD163', 'CSF1R', 'MRC1', 'MSR1'],
    'Monocyte': ['CD14', 'FCGR3A', 'S100A8', 'S100A9'],
    'Dendritic_Cell': ['CD1C', 'CLEC10A', 'FCER1A'],
    'Mast_Cell': ['KIT', 'TPSAB1', 'CPA3'],

    # Stromal
    'Fibroblast': ['COL1A1', 'COL1A2', 'LUM', 'DCN', 'VCAN'],
    'CAF': ['FAP', 'ACTA2', 'PDGFRA', 'POSTN', 'THY1'],
    'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'CLDN5'],
    'Pericyte': ['PDGFRB', 'RGS5', 'NOTCH3'],
}

# KNN parameters
N_NEIGHBORS = 15
N_PCS = 30


def load_sample(fpath: Path) -> sc.AnnData:
    """Load a single GSE183904 sample."""
    print(f"    Loading {fpath.name}...")

    with gzip.open(fpath, 'rt') as f:
        df = pd.read_csv(f, index_col=0)

    # Transpose: rows=cells, cols=genes
    df = df.T

    # Create AnnData
    adata = sc.AnnData(X=csr_matrix(df.values.astype(np.float32)))
    adata.obs_names = df.index.tolist()
    adata.var_names = df.columns.tolist()

    print(f"      Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    return adata


def build_reference():
    """Build reference AnnData from Normal gastric samples."""
    ref_h5ad = BASE / 'reference' / 'kumar_gastric_normal_reference.h5ad'

    if ref_h5ad.exists():
        print(f"  Loading existing reference: {ref_h5ad}")
        return sc.read_h5ad(ref_h5ad)

    print("  Building reference from Normal gastric samples...")

    # Load all normal samples
    adatas = []
    for fname, sample_id in NORMAL_SAMPLES:
        fpath = REFERENCE_DIR / fname
        if fpath.exists():
            ad = load_sample(fpath)
            ad.obs['sample_id'] = sample_id
            ad.obs['tissue'] = 'Normal'
            adatas.append(ad)
        else:
            print(f"    WARNING: {fpath} not found")

    if not adatas:
        print("  ERROR: No reference samples found")
        return None

    # Concatenate
    print(f"  Concatenating {len(adatas)} samples...")
    ref = sc.concat(adatas, join='outer')
    ref.obs_names_make_unique()
    print(f"  Reference: {ref.n_obs:,} cells x {ref.n_vars:,} genes")

    # Basic QC
    print("  QC filtering...")
    sc.pp.filter_cells(ref, min_genes=200)
    sc.pp.filter_genes(ref, min_cells=10)

    # Normalize
    print("  Normalizing...")
    ref.layers['counts'] = ref.X.copy()
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)

    # HVG
    print("  Finding HVGs...")
    sc.pp.highly_variable_genes(ref, n_top_genes=2000, flavor='seurat_v3', layer='counts')

    # PCA
    print("  Computing PCA...")
    sc.pp.pca(ref, n_comps=50)

    # Neighbors and clustering
    print("  Clustering...")
    sc.pp.neighbors(ref, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(ref, resolution=0.5)
    sc.tl.umap(ref)

    # Annotate with markers
    print("  Annotating with markers...")
    ref = annotate_reference_cells(ref)

    # Save
    ref.write_h5ad(ref_h5ad)
    print(f"  Saved reference: {ref_h5ad}")

    return ref


def annotate_reference_cells(ref: sc.AnnData) -> sc.AnnData:
    """Annotate reference cells using marker scoring."""
    genes_in_data = set(ref.var_names)

    # Score each cell type
    for cell_type, markers in GASTRIC_MARKERS.items():
        valid_markers = [m for m in markers if m in genes_in_data]
        if valid_markers:
            sc.tl.score_genes(ref, valid_markers, score_name=f'score_{cell_type}')
        else:
            ref.obs[f'score_{cell_type}'] = 0.0

    # Assign by max score
    score_cols = [c for c in ref.obs.columns if c.startswith('score_')]
    scores = ref.obs[score_cols].values
    cell_types = [c.replace('score_', '') for c in score_cols]

    max_indices = np.argmax(scores, axis=1)
    max_scores = np.max(scores, axis=1)

    annotations = []
    for idx, score in zip(max_indices, max_scores):
        if score > 0.1:
            annotations.append(cell_types[idx])
        else:
            annotations.append('Unknown')

    ref.obs['cell_type'] = annotations

    # Summary
    counts = Counter(annotations)
    print(f"    Reference cell type distribution:")
    for ct, n in sorted(counts.items(), key=lambda x: -x[1])[:10]:
        print(f"      {ct}: {n:,} ({n/ref.n_obs*100:.1f}%)")

    return ref


def knn_label_transfer(query, reference, n_neighbors=15, n_pcs=30):
    """Transfer labels from reference to query using KNN in PCA space."""
    print(f"  KNN label transfer (k={n_neighbors})...")

    # Find common genes
    query_genes = set(query.var_names)
    ref_genes = set(reference.var_names)
    common = sorted(query_genes & ref_genes)

    print(f"    Query genes: {len(query_genes)}")
    print(f"    Reference genes: {len(ref_genes)}")
    print(f"    Common genes: {len(common)}")

    if len(common) < 50:
        print(f"    WARNING: Only {len(common)} common genes - insufficient for transfer")
        return None, None

    # Subset to common genes
    query_sub = query[:, common].copy()
    ref_sub = reference[:, common].copy()

    # Get expression matrices
    X_query = query_sub.X.toarray() if issparse(query_sub.X) else query_sub.X
    X_ref = ref_sub.X.toarray() if issparse(ref_sub.X) else ref_sub.X

    # Joint PCA
    print("    Computing joint PCA...")
    X_combined = np.vstack([X_ref, X_query])
    n_comps = min(n_pcs, min(X_combined.shape) - 1, len(common) - 1)

    pca = PCA(n_components=n_comps)
    pca.fit(X_combined)

    X_ref_pca = pca.transform(X_ref)
    X_query_pca = pca.transform(X_query)

    # Train KNN
    print(f"    Training KNN on {reference.n_obs:,} reference cells...")
    labels = reference.obs['cell_type'].values
    le = LabelEncoder()
    labels_encoded = le.fit_transform(labels)

    knn = KNeighborsClassifier(n_neighbors=min(n_neighbors, len(labels)), weights='distance')
    knn.fit(X_ref_pca, labels_encoded)

    # Predict
    print(f"    Predicting labels for {query.n_obs:,} cells...")
    pred_encoded = knn.predict(X_query_pca)
    pred_labels = le.inverse_transform(pred_encoded)

    # Confidence
    proba = knn.predict_proba(X_query_pca)
    confidence = np.max(proba, axis=1)

    print(f"    Mean confidence: {confidence.mean():.3f}")
    print(f"    Top predictions: {Counter(pred_labels).most_common(5)}")

    return pred_labels, confidence


def annotate_sample(sample_id: str, reference):
    """Annotate a single G4X sample with reference labels."""
    input_path = RESULTS / f'{sample_id}_annotated.h5ad'

    if not input_path.exists():
        print(f"  WARNING: {input_path} not found")
        return None

    print(f"\n  Annotating {sample_id}...")
    adata = sc.read_h5ad(input_path)
    print(f"    Cells: {adata.n_obs:,}")

    # Transfer labels
    labels, confidence = knn_label_transfer(adata, reference, N_NEIGHBORS, N_PCS)

    if labels is None:
        adata.obs['celltype_reference'] = 'Transfer_Failed'
        adata.obs['celltype_reference_conf'] = 0.0
    else:
        adata.obs['celltype_reference'] = labels
        adata.obs['celltype_reference_conf'] = confidence

        # Coarse labels
        coarse = [ct.split('_')[0] if '_' in ct else ct for ct in labels]
        adata.obs['celltype_reference_coarse'] = coarse

    # Save
    output_path = RESULTS / f'{sample_id}_reference_annotated.h5ad'
    adata.write_h5ad(output_path)
    print(f"    Saved: {output_path}")

    return adata


def create_visualizations(annotated: dict, reference):
    """Create visualization figures."""
    FIGURES.mkdir(parents=True, exist_ok=True)

    # 1. Reference UMAP
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(reference, color='cell_type', ax=ax, show=False,
              legend_loc='on data', legend_fontsize=6,
              title='Kumar et al. Normal Gastric Reference')
    plt.tight_layout()
    fig.savefig(FIGURES / 'reference_umap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {FIGURES / 'reference_umap.png'}")

    # 2. G4X samples with reference annotations
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for i, (sample_id, (adata, stage)) in enumerate(annotated.items()):
        ax = axes[i]
        if 'celltype_reference' in adata.obs.columns and 'X_umap' in adata.obsm:
            sc.pl.umap(adata, color='celltype_reference', ax=ax, show=False,
                      legend_loc='on data', legend_fontsize=5,
                      title=f'{sample_id}: {stage}')

    plt.suptitle('G4X Samples - Reference-Based Annotation', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(FIGURES / 'reference_annotation_umap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {FIGURES / 'reference_annotation_umap.png'}")

    # 3. Comparison: Marker vs Reference
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, (sample_id, (adata, stage)) in enumerate(annotated.items()):
        ax = axes[i]

        if 'celltype_markers' in adata.obs.columns and 'celltype_reference' in adata.obs.columns:
            # Get proportions
            marker_props = adata.obs['celltype_markers'].value_counts(normalize=True)
            ref_props = adata.obs['celltype_reference'].value_counts(normalize=True)

            # Key types
            key_types = ['Goblet', 'Fibroblast', 'Plasma_Cell', 'Macrophage', 'T_CD8',
                        'Gastric_Pit_Mucous', 'Gastric_Chief_Parietal', 'Enteroendocrine']
            x = np.arange(len(key_types))
            width = 0.35

            # Map marker types to reference types
            marker_vals = []
            for kt in key_types:
                if kt == 'Gastric_Pit_Mucous':
                    val = marker_props.get('Gastric_Pit', 0)
                elif kt == 'Gastric_Chief_Parietal':
                    val = marker_props.get('Gastric_Chief', 0)
                else:
                    val = marker_props.get(kt, 0)
                marker_vals.append(val * 100)

            ref_vals = [ref_props.get(kt, 0) * 100 for kt in key_types]

            ax.bar(x - width/2, marker_vals, width, label='Marker', color='steelblue')
            ax.bar(x + width/2, ref_vals, width, label='Reference', color='coral')

            ax.set_xlabel('Cell Type')
            ax.set_ylabel('Proportion (%)')
            ax.set_title(f'{sample_id}: {stage}', fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels(key_types, rotation=45, ha='right', fontsize=7)
            ax.legend()

    plt.suptitle('Marker vs Reference Annotation', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(FIGURES / 'validation_comparison.png', dpi=150)
    plt.close()
    print(f"  Saved: {FIGURES / 'validation_comparison.png'}")

    # 4. Confidence histogram
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    for i, (sample_id, (adata, stage)) in enumerate(annotated.items()):
        ax = axes[i]
        if 'celltype_reference_conf' in adata.obs.columns:
            conf = adata.obs['celltype_reference_conf']
            ax.hist(conf, bins=50, color='steelblue', edgecolor='white')
            ax.axvline(conf.mean(), color='red', linestyle='--', label=f'Mean: {conf.mean():.2f}')
            ax.set_xlabel('Confidence')
            ax.set_ylabel('Cell Count')
            ax.set_title(f'{sample_id}: {stage}')
            ax.legend()

    plt.suptitle('Label Transfer Confidence', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(FIGURES / 'reference_confidence.png', dpi=150)
    plt.close()
    print(f"  Saved: {FIGURES / 'reference_confidence.png'}")


def save_comparison_csv(annotated: dict):
    """Save detailed comparison CSV."""
    rows = []

    for sample_id, (adata, stage) in annotated.items():
        if 'celltype_reference' not in adata.obs.columns:
            continue

        # Reference proportions
        ref_counts = adata.obs['celltype_reference'].value_counts()
        for ct, n in ref_counts.items():
            rows.append({
                'sample_id': sample_id,
                'stage': stage,
                'method': 'reference',
                'cell_type': ct,
                'count': n,
                'proportion': n / adata.n_obs
            })

        # Marker proportions
        if 'celltype_markers' in adata.obs.columns:
            marker_counts = adata.obs['celltype_markers'].value_counts()
            for ct, n in marker_counts.items():
                rows.append({
                    'sample_id': sample_id,
                    'stage': stage,
                    'method': 'markers',
                    'cell_type': ct,
                    'count': n,
                    'proportion': n / adata.n_obs
                })

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / 'reference_annotations.csv', index=False)
    print(f"\nSaved: {RESULTS / 'reference_annotations.csv'}")


def main():
    print("=" * 60)
    print("Reference-Based Cell Type Annotation")
    print("Kumar et al. 2022 Normal Gastric Reference")
    print("=" * 60)

    # Step 1: Build reference
    print("\n[1/3] Building reference from Normal gastric samples...")
    reference = build_reference()

    if reference is None:
        print("ERROR: Could not build reference")
        return

    print(f"\n  Reference ready: {reference.n_obs:,} cells, "
          f"{reference.obs['cell_type'].nunique()} cell types")

    # Step 2: Annotate G4X samples
    print("\n[2/3] Annotating G4X samples...")
    annotated = {}

    for sample_id, stage in SAMPLES.items():
        adata = annotate_sample(sample_id, reference)
        if adata is not None:
            annotated[sample_id] = (adata, stage)

    # Step 3: Visualize
    print("\n[3/3] Creating visualizations...")
    create_visualizations(annotated, reference)
    save_comparison_csv(annotated)

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

    # Goblet validation
    print("\n" + "-" * 40)
    print("GOBLET VALIDATION (Marker: E02=6%, F02=7%, G02=48%)")
    print("-" * 40)
    for sample_id, (adata, stage) in annotated.items():
        if 'celltype_reference' in adata.obs.columns:
            goblet_n = (adata.obs['celltype_reference'] == 'Goblet').sum()
            goblet_pct = goblet_n / adata.n_obs * 100
            print(f"  {sample_id} ({stage}): {goblet_pct:.1f}% Goblet")

    print(f"\nOutputs: {RESULTS}")


if __name__ == '__main__':
    main()
