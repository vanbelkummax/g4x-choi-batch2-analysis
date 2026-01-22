#!/usr/bin/env python3
"""
Label transfer from GSE134520 reference to G4X data
Uses scanpy.tl.ingest for label transfer
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 2

# Paths
DATA_DIR = Path("/home/user/g4x-choi-batch2-analysis")
OUTPUT_DIR = DATA_DIR / "output"

print("="*60)
print("Label Transfer: GSE134520 → G4X")
print("="*60)

# Load reference
print("\n1. Loading reference (GSE134520)...")
ref = sc.read_h5ad(OUTPUT_DIR / "GSE134520_reference.h5ad")
print(f"   Reference: {ref.n_obs:,} cells, {ref.n_vars:,} genes")
print(f"   Cell types: {ref.obs['cell_type'].nunique()}")

# Load G4X query data
print("\n2. Loading G4X query data...")
query = sc.read_h5ad(OUTPUT_DIR / "pilot_preprocessed.h5ad")
print(f"   Query: {query.n_obs:,} cells, {query.n_vars:,} genes")

# Find common genes
print("\n3. Finding common genes...")
common_genes = ref.var_names.intersection(query.var_names)
print(f"   Common genes: {len(common_genes)}")

if len(common_genes) < 100:
    print("   ⚠️ Warning: Very few common genes. Results may be unreliable.")

# Subset to common genes
ref_common = ref[:, common_genes].copy()
query_common = query[:, common_genes].copy()

# Normalize query if not already
print("\n4. Preparing query data...")
if query_common.X.max() > 50:  # Likely not log-normalized
    sc.pp.normalize_total(query_common, target_sum=1e4)
    sc.pp.log1p(query_common)

# Re-run PCA on reference with common genes
print("\n5. Running PCA on reference...")
sc.pp.scale(ref_common, max_value=10)
sc.tl.pca(ref_common, n_comps=30)

# Ingest: map query onto reference
print("\n6. Label transfer via ingest...")
sc.tl.ingest(query_common, ref_common, obs='cell_type')

# Copy results back to original query
query.obs['cell_type_transferred'] = query_common.obs['cell_type']
query.obsm['X_umap_transferred'] = query_common.obsm['X_umap']

# Summary
print("\n" + "="*60)
print("LABEL TRANSFER COMPLETE")
print("="*60)

for sample in ['E02_Normal', 'F02_Metaplasia', 'G02_Cancer']:
    print(f"\n{sample}:")
    mask = query.obs['sample'] == sample
    counts = query.obs.loc[mask, 'cell_type_transferred'].value_counts()
    total = counts.sum()
    for ct, n in counts.head(8).items():
        print(f"   {ct}: {n:,} ({100*n/total:.1f}%)")

# Save
output_path = OUTPUT_DIR / "pilot_annotated.h5ad"
query.write(output_path)
print(f"\n✅ Saved to {output_path}")

# Visualization
print("\n7. Generating figures...")

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Row 1: Cell types per sample
for i, sample in enumerate(['E02_Normal', 'F02_Metaplasia', 'G02_Cancer']):
    ax = axes[0, i]
    sample_data = query[query.obs['sample'] == sample].copy()
    sc.pl.umap(sample_data, color='cell_type_transferred', ax=ax, show=False,
               title=f'{sample}\n({sample_data.n_obs:,} cells)',
               legend_loc='right margin' if i == 2 else 'none')

# Row 2: Key cell type distributions
ax = axes[1, 0]
# Stacked bar chart of cell types
ct_counts = pd.crosstab(query.obs['sample'], query.obs['cell_type_transferred'])
ct_pct = ct_counts.div(ct_counts.sum(axis=1), axis=0) * 100
ct_pct.plot(kind='bar', stacked=True, ax=ax, legend=False)
ax.set_xlabel('')
ax.set_ylabel('Percentage')
ax.set_title('Cell Type Composition')
ax.tick_params(axis='x', rotation=45)

# Enterocyte and Goblet focus
ax = axes[1, 1]
key_types = ['Enterocyte', 'Goblet', 'PMC', 'GMC']
for ct in key_types:
    if ct in ct_pct.columns:
        ax.bar(ct_pct.index, ct_pct[ct], label=ct, alpha=0.7)
ax.set_ylabel('Percentage')
ax.set_title('Key Epithelial Types')
ax.legend()
ax.tick_params(axis='x', rotation=45)

# Cancer vs non-cancer
ax = axes[1, 2]
cancer_types = ['Cancer', 'MSC']  # Cancer and stem cells
for ct in cancer_types:
    if ct in ct_pct.columns:
        ax.bar(ct_pct.index, ct_pct[ct], label=ct, alpha=0.7)
ax.set_ylabel('Percentage')
ax.set_title('Cancer/Stem Cells')
ax.legend()
ax.tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'fig_label_transfer.png', dpi=150, bbox_inches='tight')
plt.savefig(OUTPUT_DIR / 'fig_label_transfer.pdf', bbox_inches='tight')
print(f"✅ Saved figures to {OUTPUT_DIR}/fig_label_transfer.png")

print("\n✅ Label transfer complete!")
