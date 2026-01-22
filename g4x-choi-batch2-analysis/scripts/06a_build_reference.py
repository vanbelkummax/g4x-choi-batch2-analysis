#!/usr/bin/env python3
"""
Build annotated scRNA-seq reference from GSE134520
Zhang et al. 2019 Cell Reports - Gastric premalignant lesions
NAG → CAG → IM → EGC progression

This creates a reference AnnData with cell type labels for transfer to G4X data
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import gzip

sc.settings.verbosity = 2

# Paths
REF_DIR = Path("/home/user/g4x-choi-batch2-analysis/reference")
OUTPUT_DIR = Path("/home/user/g4x-choi-batch2-analysis/output")
OUTPUT_DIR.mkdir(exist_ok=True)

# Sample metadata
SAMPLES = {
    'GSM3954946_processed_NAG1.txt.gz': ('NAG', 'NAG1'),
    'GSM3954947_processed_NAG2.txt.gz': ('NAG', 'NAG2'),
    'GSM3954948_processed_NAG3.txt.gz': ('NAG', 'NAG3'),
    'GSM3954949_processed_CAG1.txt.gz': ('CAG', 'CAG1'),
    'GSM3954950_processed_CAG2.txt.gz': ('CAG', 'CAG2'),
    'GSM3954951_processed_CAG3.txt.gz': ('CAG', 'CAG3'),
    'GSM3954952_processed_IMW1.txt.gz': ('IM_Wild', 'IMW1'),
    'GSM3954953_processed_IMW2.txt.gz': ('IM_Wild', 'IMW2'),
    'GSM3954954_processed_IMS1.txt.gz': ('IM_Severe', 'IMS1'),
    'GSM3954955_processed_IMS2.txt.gz': ('IM_Severe', 'IMS2'),
    'GSM3954956_processed_IMS3.txt.gz': ('IM_Severe', 'IMS3'),
    'GSM3954957_processed_IMS4.txt.gz': ('IM_Severe', 'IMS4'),
    'GSM3954958_processed_EGC.txt.gz': ('EGC', 'EGC1'),
}

# Cell type markers from Zhang et al. 2019
MARKERS = {
    # Epithelial - Gastric
    'PMC': ['GKN1', 'GKN2', 'TFF1', 'MUC5AC', 'DPCR1', 'PSCA', 'MUC1', 'TFF2'],
    'GMC': ['PGC', 'MUC6', 'PRR4', 'C6orf58', 'LYZ'],
    'Chief': ['PGA3', 'PGA4', 'CHIA', 'LIPF', 'PDIA2'],
    'Neck_like': ['CXCL3', 'LYZ', 'PGC', 'CXCL2'],
    'PC': ['ATP4A', 'ATP4B'],  # Parietal cells - if present

    # Epithelial - Intestinal (metaplasia)
    'Enterocyte': ['FABP1', 'RBP2', 'FABP2', 'ANPEP', 'APOA4', 'APOA1', 'SI'],
    'Goblet': ['SPINK4', 'TFF3', 'ITLN1', 'ZG16', 'MUC2', 'WFDC2'],
    'Enteroendocrine': ['CHGA', 'PCSK1N', 'SCG5', 'CHGB', 'TPH1'],

    # Stem/Progenitor
    'MSC': ['OLFM4', 'REG1A', 'CLDN4', 'TSPAN8', 'LGR5', 'SOX9'],

    # Cancer
    'Cancer': ['REG4', 'CLDN7', 'TSPAN8', 'KRT18', 'LGALS3', 'CEACAM6'],

    # Immune
    'T_cell': ['CCL5', 'CD52', 'CCL4', 'KLRB1', 'GZMA', 'CD3D', 'CD3E'],
    'B_cell': ['IGJ', 'IGLL5', 'MZB1', 'CD79A', 'MS4A1'],
    'Macrophage': ['HLA-DRA', 'HLA-DPB1', 'C1QA', 'C1QB', 'CD68'],
    'Mast': ['TPSAB1', 'CPA3', 'HPGDS'],

    # Stromal
    'Fibroblast': ['CXCL14', 'LUM', 'DCN', 'COL1A1', 'COL1A2', 'POSTN'],
    'Endothelial': ['PLVAP', 'CD320', 'SPARCL1', 'CLDN5', 'PECAM1', 'VWF'],
    'Smooth_Muscle': ['RGS5', 'ACTA2', 'TAGLN', 'MYL9'],
}

def load_sample(filepath):
    """Load a single sample from GSE134520 format (genes x cells)"""
    print(f"  Loading {filepath.name}...")

    # Read transposed (genes as rows, cells as columns)
    df = pd.read_csv(filepath, sep='\t', index_col=0, compression='gzip')

    # Transpose to cells x genes
    df = df.T

    # Create AnnData
    adata = sc.AnnData(df.values.astype(np.float32))
    adata.obs_names = df.index
    adata.var_names = df.columns

    return adata

def main():
    print("="*60)
    print("Building GSE134520 Reference")
    print("="*60)

    # Load all samples
    print("\n1. Loading samples...")
    adatas = []
    for filename, (condition, sample_id) in SAMPLES.items():
        filepath = REF_DIR / filename
        if filepath.exists():
            adata = load_sample(filepath)
            adata.obs['condition'] = condition
            adata.obs['sample'] = sample_id
            adatas.append(adata)
            print(f"     {sample_id}: {adata.n_obs} cells")

    # Concatenate
    print("\n2. Concatenating samples...")
    adata = sc.concat(adatas, join='outer')
    adata.obs_names_make_unique()
    print(f"   Total: {adata.n_obs} cells, {adata.n_vars} genes")

    # Basic QC
    print("\n3. Quality control...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

    # Filter
    adata = adata[adata.obs['pct_counts_mt'] < 20, :].copy()
    print(f"   After QC: {adata.n_obs} cells, {adata.n_vars} genes")

    # Store raw
    adata.raw = adata.copy()

    # Normalize and log
    print("\n4. Normalizing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG
    print("\n5. Finding variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key='sample')

    # Scale and PCA
    print("\n6. PCA...")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50)

    # Batch correction with Harmony
    print("\n7. Batch correction (Harmony)...")
    try:
        import harmonypy
        sc.external.pp.harmony_integrate(adata, key='sample')
        use_rep = 'X_pca_harmony'
    except:
        print("   Harmony not available, using uncorrected PCA")
        use_rep = 'X_pca'

    # Neighbors and UMAP
    print("\n8. Computing neighbors and UMAP...")
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=15)
    sc.tl.umap(adata)

    # Clustering
    print("\n9. Clustering...")
    sc.tl.leiden(adata, resolution=0.8)

    # Cell type annotation using markers
    print("\n10. Annotating cell types...")

    # Score each cell type
    for cell_type, markers in MARKERS.items():
        valid_markers = [m for m in markers if m in adata.var_names]
        if valid_markers:
            sc.tl.score_genes(adata, valid_markers, score_name=f'score_{cell_type}')

    # Assign cell types based on max score
    score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
    if score_cols:
        scores = adata.obs[score_cols].values
        max_idx = np.argmax(scores, axis=1)
        max_scores = np.max(scores, axis=1)

        cell_types = [score_cols[i].replace('score_', '') for i in max_idx]
        # Low confidence threshold
        cell_types = [ct if max_scores[i] > 0.1 else 'Unknown'
                     for i, ct in enumerate(cell_types)]

        adata.obs['cell_type'] = pd.Categorical(cell_types)
        adata.obs['cell_type_score'] = max_scores

    # Summary
    print("\n" + "="*60)
    print("REFERENCE SUMMARY")
    print("="*60)
    print(f"Total cells: {adata.n_obs:,}")
    print(f"Total genes: {adata.n_vars:,}")
    print("\nCells per condition:")
    print(adata.obs['condition'].value_counts())
    print("\nCell types:")
    print(adata.obs['cell_type'].value_counts())

    # Save
    output_path = OUTPUT_DIR / "GSE134520_reference.h5ad"
    adata.write(output_path)
    print(f"\n✅ Saved reference to {output_path}")

    # Quick visualization
    print("\n11. Generating QC plots...")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    sc.pl.umap(adata, color='condition', ax=axes[0], show=False, title='Condition')
    sc.pl.umap(adata, color='cell_type', ax=axes[1], show=False, title='Cell Type',
               legend_loc='right margin')
    sc.pl.umap(adata, color='leiden', ax=axes[2], show=False, title='Leiden Clusters')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'reference_umap.png', dpi=150, bbox_inches='tight')
    print(f"✅ Saved UMAP to {OUTPUT_DIR}/reference_umap.png")

    return adata

if __name__ == "__main__":
    adata = main()
