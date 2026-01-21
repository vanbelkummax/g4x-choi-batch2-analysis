#!/usr/bin/env python
"""
G4X Multimodal Statistical Analysis
====================================
Gastric cancer spatial transcriptomics with protein data

Analyses:
1. Cell type composition comparison (chi-square, Fisher's exact)
2. Protein-RNA correlation patterns
3. Spatially variable genes (Moran's I)
4. Neighborhood enrichment
5. Immune infiltration comparison
"""

import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import chi2_contingency, fisher_exact, spearmanr, pearsonr
import warnings
warnings.filterwarnings('ignore')

# Settings
sc.settings.verbosity = 1
np.random.seed(42)

print("=" * 70)
print("G4X MULTIMODAL STATISTICAL ANALYSIS - GASTRIC CANCER")
print("=" * 70)

# =============================================================================
# LOAD DATA
# =============================================================================
print("\n[1] LOADING DATA...")
data_dir = "/home/user/spatial-hackathon-2026/outputs/adata/g4x/"

adata_a01 = sc.read_h5ad(f"{data_dir}/A01_g4x.h5ad")
adata_b01 = sc.read_h5ad(f"{data_dir}/B01_g4x.h5ad")

print(f"A01: {adata_a01.n_obs} cells, {adata_a01.n_vars} genes")
print(f"B01: {adata_b01.n_obs} cells, {adata_b01.n_vars} genes")

# Check what's available
print("\nA01 obs columns:", list(adata_a01.obs.columns)[:15])
print("A01 var columns:", list(adata_a01.var.columns)[:10])
print("A01 layers:", list(adata_a01.layers.keys()) if adata_a01.layers else "None")
print("A01 obsm:", list(adata_a01.obsm.keys()) if adata_a01.obsm else "None")

# Check for protein data
if 'protein' in adata_a01.obsm:
    print(f"A01 protein features: {adata_a01.obsm['protein'].shape[1]}")
if hasattr(adata_a01, 'uns') and 'protein_names' in adata_a01.uns:
    print(f"Protein names: {adata_a01.uns['protein_names'][:5]}...")

# =============================================================================
# CELL TYPE COMPOSITION COMPARISON
# =============================================================================
print("\n" + "=" * 70)
print("[2] CELL TYPE COMPOSITION COMPARISON")
print("=" * 70)

# Find cell type column
cell_type_col = None
for col in ['cell_type', 'leiden', 'cluster', 'celltype', 'annotation']:
    if col in adata_a01.obs.columns:
        cell_type_col = col
        break

if cell_type_col is None:
    # Use leiden clustering if no annotation
    cell_type_col = 'leiden' if 'leiden' in adata_a01.obs.columns else None

if cell_type_col:
    print(f"Using cell type column: {cell_type_col}")

    # Get counts
    a01_counts = adata_a01.obs[cell_type_col].value_counts()
    b01_counts = adata_b01.obs[cell_type_col].value_counts()

    print(f"\nA01 cell types:\n{a01_counts}")
    print(f"\nB01 cell types:\n{b01_counts}")

    # Create contingency table
    all_types = sorted(set(a01_counts.index) | set(b01_counts.index))
    contingency = pd.DataFrame(index=all_types, columns=['A01', 'B01'])
    for ct in all_types:
        contingency.loc[ct, 'A01'] = a01_counts.get(ct, 0)
        contingency.loc[ct, 'B01'] = b01_counts.get(ct, 0)

    contingency = contingency.fillna(0).astype(int)
    print(f"\nContingency table:\n{contingency}")

    # Chi-square test
    chi2, p_chi2, dof, expected = chi2_contingency(contingency.values)
    print(f"\nCHI-SQUARE TEST:")
    print(f"  Chi2 statistic: {chi2:.2f}")
    print(f"  Degrees of freedom: {dof}")
    print(f"  P-value: {p_chi2:.2e}")

    # Effect size (Cramer's V)
    n = contingency.values.sum()
    min_dim = min(contingency.shape[0] - 1, contingency.shape[1] - 1)
    cramers_v = np.sqrt(chi2 / (n * min_dim)) if min_dim > 0 else 0
    print(f"  Cramer's V: {cramers_v:.3f}")

    # Per-cell-type Fisher's exact tests (for 2x2 comparisons)
    print(f"\nPER-CELL-TYPE ENRICHMENT (Fisher's exact):")
    fisher_results = []
    for ct in all_types:
        # Create 2x2 table: this type vs others, A01 vs B01
        a01_this = contingency.loc[ct, 'A01']
        a01_other = contingency['A01'].sum() - a01_this
        b01_this = contingency.loc[ct, 'B01']
        b01_other = contingency['B01'].sum() - b01_this

        table_2x2 = [[a01_this, b01_this], [a01_other, b01_other]]
        odds_ratio, p_fisher = fisher_exact(table_2x2)

        # Proportion in each sample
        prop_a01 = a01_this / (a01_this + a01_other) * 100
        prop_b01 = b01_this / (b01_this + b01_other) * 100

        fisher_results.append({
            'cell_type': ct,
            'A01_count': a01_this,
            'B01_count': b01_this,
            'A01_pct': prop_a01,
            'B01_pct': prop_b01,
            'odds_ratio': odds_ratio,
            'p_value': p_fisher
        })

    fisher_df = pd.DataFrame(fisher_results)
    # Multiple testing correction (Bonferroni)
    fisher_df['p_adj'] = np.minimum(fisher_df['p_value'] * len(fisher_df), 1.0)
    fisher_df = fisher_df.sort_values('p_value')

    print(fisher_df.to_string(index=False))

    # Significant differences
    sig_types = fisher_df[fisher_df['p_adj'] < 0.05]
    if len(sig_types) > 0:
        print(f"\nSIGNIFICANT DIFFERENCES (p_adj < 0.05): {len(sig_types)}")
        for _, row in sig_types.iterrows():
            enriched = "A01" if row['odds_ratio'] > 1 else "B01"
            print(f"  {row['cell_type']}: enriched in {enriched} (OR={row['odds_ratio']:.2f}, p_adj={row['p_adj']:.2e})")
else:
    print("No cell type column found - skipping composition analysis")

# =============================================================================
# PROTEIN-RNA CORRELATION
# =============================================================================
print("\n" + "=" * 70)
print("[3] PROTEIN-RNA CORRELATION ANALYSIS")
print("=" * 70)

def analyze_protein_rna(adata, sample_name):
    """Analyze protein-RNA correlations for a sample."""

    # Check for protein data in different locations
    protein_data = None
    protein_names = None

    if 'protein' in adata.obsm:
        protein_data = adata.obsm['protein']
        if 'protein_names' in adata.uns:
            protein_names = adata.uns['protein_names']
    elif 'X_protein' in adata.obsm:
        protein_data = adata.obsm['X_protein']

    # Check var for protein markers
    if protein_data is None and 'feature_types' in adata.var.columns:
        protein_mask = adata.var['feature_types'] == 'Antibody Capture'
        if protein_mask.sum() > 0:
            protein_names = adata.var_names[protein_mask].tolist()
            protein_data = adata[:, protein_mask].X
            if hasattr(protein_data, 'toarray'):
                protein_data = protein_data.toarray()

    if protein_data is None:
        print(f"{sample_name}: No protein data found")
        return None

    if protein_names is None:
        protein_names = [f"Protein_{i}" for i in range(protein_data.shape[1])]

    print(f"\n{sample_name}: {len(protein_names)} proteins detected")

    # Get RNA data (use raw if available)
    if adata.raw is not None:
        rna_data = adata.raw.X
        rna_names = adata.raw.var_names
    else:
        rna_data = adata.X
        rna_names = adata.var_names

    if hasattr(rna_data, 'toarray'):
        rna_data = rna_data.toarray()
    if hasattr(protein_data, 'toarray'):
        protein_data = protein_data.toarray()

    # Find matching protein-gene pairs
    correlations = []

    for i, prot_name in enumerate(protein_names):
        # Clean protein name for matching
        gene_name = prot_name.replace('_TotalSeqC', '').replace('_', '-').upper()

        # Try to find matching gene
        matches = [g for g in rna_names if gene_name in g.upper() or g.upper() in gene_name]

        if len(matches) > 0:
            gene_idx = list(rna_names).index(matches[0])
            prot_vals = protein_data[:, i]
            rna_vals = rna_data[:, gene_idx]

            # Calculate correlation (Spearman is more robust)
            if np.std(prot_vals) > 0 and np.std(rna_vals) > 0:
                r_spearman, p_spearman = spearmanr(prot_vals, rna_vals)
                r_pearson, p_pearson = pearsonr(prot_vals, rna_vals)

                correlations.append({
                    'protein': prot_name,
                    'gene': matches[0],
                    'spearman_r': r_spearman,
                    'spearman_p': p_spearman,
                    'pearson_r': r_pearson,
                    'pearson_p': p_pearson
                })

    if len(correlations) > 0:
        corr_df = pd.DataFrame(correlations)
        corr_df = corr_df.sort_values('spearman_r', ascending=False)

        print(f"\nProtein-RNA correlations (matched pairs):")
        print(corr_df.to_string(index=False))

        # Summary statistics
        mean_r = corr_df['spearman_r'].mean()
        sig_pos = (corr_df['spearman_r'] > 0.3).sum()
        sig_neg = (corr_df['spearman_r'] < -0.3).sum()

        print(f"\nSummary:")
        print(f"  Mean Spearman r: {mean_r:.3f}")
        print(f"  Strong positive (r > 0.3): {sig_pos}")
        print(f"  Strong negative (r < -0.3): {sig_neg}")

        return corr_df
    else:
        print(f"No matching protein-gene pairs found")
        return None

# Run for both samples
corr_a01 = analyze_protein_rna(adata_a01, "A01")
corr_b01 = analyze_protein_rna(adata_b01, "B01")

# =============================================================================
# SPATIALLY VARIABLE GENES (MORAN'S I)
# =============================================================================
print("\n" + "=" * 70)
print("[4] SPATIALLY VARIABLE GENES (MORAN'S I)")
print("=" * 70)

def compute_morans_i(adata, sample_name, n_genes=50):
    """Compute Moran's I for spatially variable genes."""

    print(f"\n{sample_name}: Computing spatial autocorrelation...")

    # Ensure spatial coordinates
    if 'spatial' not in adata.obsm:
        print(f"  No spatial coordinates found")
        return None

    # Build spatial graph if needed
    if 'spatial_connectivities' not in adata.obsp:
        print(f"  Building spatial neighbors graph...")
        sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

    # Select highly variable genes if available
    if 'highly_variable' in adata.var.columns:
        genes = adata.var_names[adata.var['highly_variable']][:n_genes]
    else:
        # Use top expressed genes
        mean_expr = np.array(adata.X.mean(axis=0)).flatten()
        top_idx = np.argsort(mean_expr)[-n_genes:]
        genes = adata.var_names[top_idx]

    print(f"  Analyzing {len(genes)} genes...")

    # Compute Moran's I
    try:
        sq.gr.spatial_autocorr(
            adata,
            genes=genes,
            mode='moran',
            n_perms=100,
            n_jobs=4
        )

        # Get results
        moran_df = adata.uns['moranI'].copy()
        moran_df = moran_df.sort_values('I', ascending=False)

        print(f"\nTop 15 spatially variable genes (Moran's I):")
        print(moran_df.head(15).to_string())

        # Significant genes
        sig_genes = moran_df[moran_df['pval_norm'] < 0.05]
        print(f"\nSignificant SVGs (p < 0.05): {len(sig_genes)}")

        return moran_df

    except Exception as e:
        print(f"  Error computing Moran's I: {e}")
        return None

svg_a01 = compute_morans_i(adata_a01, "A01")
svg_b01 = compute_morans_i(adata_b01, "B01")

# Compare SVGs between samples
if svg_a01 is not None and svg_b01 is not None:
    print("\n" + "-" * 50)
    print("SVG COMPARISON BETWEEN SAMPLES")
    print("-" * 50)

    # Top SVGs in each
    top_a01 = set(svg_a01[svg_a01['pval_norm'] < 0.05].head(50).index)
    top_b01 = set(svg_b01[svg_b01['pval_norm'] < 0.05].head(50).index)

    shared = top_a01 & top_b01
    a01_specific = top_a01 - top_b01
    b01_specific = top_b01 - top_a01

    print(f"Top 50 significant SVGs:")
    print(f"  Shared: {len(shared)}")
    print(f"  A01-specific: {len(a01_specific)}")
    print(f"  B01-specific: {len(b01_specific)}")

    if len(shared) > 0:
        print(f"\nShared top SVGs: {list(shared)[:10]}")

# =============================================================================
# NEIGHBORHOOD ENRICHMENT
# =============================================================================
print("\n" + "=" * 70)
print("[5] NEIGHBORHOOD ENRICHMENT ANALYSIS")
print("=" * 70)

def compute_neighborhood_enrichment(adata, sample_name, cell_type_col):
    """Compute neighborhood enrichment between cell types."""

    print(f"\n{sample_name}: Computing neighborhood enrichment...")

    if cell_type_col not in adata.obs.columns:
        print(f"  No cell type column '{cell_type_col}' found")
        return None

    # Ensure spatial graph
    if 'spatial_connectivities' not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

    try:
        sq.gr.nhood_enrichment(
            adata,
            cluster_key=cell_type_col,
            n_perms=1000,
            seed=42
        )

        # Extract results
        zscore = adata.uns[f'{cell_type_col}_nhood_enrichment']['zscore']

        print(f"\nNeighborhood enrichment Z-scores:")
        print(pd.DataFrame(zscore,
                          index=adata.obs[cell_type_col].cat.categories,
                          columns=adata.obs[cell_type_col].cat.categories).round(2))

        # Find significant interactions
        print(f"\nSignificant interactions (|Z| > 2):")
        categories = list(adata.obs[cell_type_col].cat.categories)
        sig_interactions = []

        for i, cat1 in enumerate(categories):
            for j, cat2 in enumerate(categories):
                z = zscore[i, j]
                if abs(z) > 2:
                    interaction_type = "enriched" if z > 0 else "depleted"
                    sig_interactions.append({
                        'cell_type_1': cat1,
                        'cell_type_2': cat2,
                        'zscore': z,
                        'interaction': interaction_type
                    })

        if len(sig_interactions) > 0:
            sig_df = pd.DataFrame(sig_interactions)
            sig_df = sig_df.sort_values('zscore', key=abs, ascending=False)
            print(sig_df.head(20).to_string(index=False))

        return zscore

    except Exception as e:
        print(f"  Error: {e}")
        return None

if cell_type_col:
    nhood_a01 = compute_neighborhood_enrichment(adata_a01, "A01", cell_type_col)
    nhood_b01 = compute_neighborhood_enrichment(adata_b01, "B01", cell_type_col)

# =============================================================================
# IMMUNE INFILTRATION COMPARISON
# =============================================================================
print("\n" + "=" * 70)
print("[6] IMMUNE INFILTRATION ANALYSIS")
print("=" * 70)

# Define immune cell markers
immune_markers = {
    'T_cells': ['CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B'],
    'B_cells': ['CD19', 'CD79A', 'CD79B', 'MS4A1', 'BANK1'],
    'NK_cells': ['NKG7', 'GNLY', 'KLRD1', 'KLRB1', 'NCR1'],
    'Macrophages': ['CD68', 'CD163', 'CSF1R', 'CD14', 'MARCO'],
    'Dendritic': ['ITGAX', 'CD1C', 'CLEC9A', 'FCER1A'],
    'Mast_cells': ['TPSAB1', 'TPSB2', 'CPA3', 'KIT'],
    'Neutrophils': ['S100A8', 'S100A9', 'FCGR3B', 'CSF3R'],
    'Tregs': ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2']
}

def score_immune_cells(adata, sample_name):
    """Score immune cell signatures."""

    print(f"\n{sample_name}: Scoring immune signatures...")

    immune_scores = {}

    for cell_type, markers in immune_markers.items():
        # Find available markers
        available = [m for m in markers if m in adata.var_names]

        if len(available) >= 2:
            # Use scanpy score_genes
            sc.tl.score_genes(adata, available, score_name=f'{cell_type}_score')
            score = adata.obs[f'{cell_type}_score'].mean()
            immune_scores[cell_type] = {
                'mean_score': score,
                'markers_used': len(available),
                'std': adata.obs[f'{cell_type}_score'].std()
            }

    return immune_scores

# Score both samples
immune_a01 = score_immune_cells(adata_a01, "A01")
immune_b01 = score_immune_cells(adata_b01, "B01")

# Compare immune infiltration
print("\nIMMUNE SIGNATURE COMPARISON:")
print("-" * 60)

immune_comparison = []
for cell_type in immune_markers.keys():
    if cell_type in immune_a01 and cell_type in immune_b01:
        score_a01_col = f'{cell_type}_score'
        score_b01_col = f'{cell_type}_score'

        # Mann-Whitney U test
        u_stat, p_val = stats.mannwhitneyu(
            adata_a01.obs[score_a01_col],
            adata_b01.obs[score_b01_col],
            alternative='two-sided'
        )

        # Effect size (rank-biserial correlation)
        n1, n2 = len(adata_a01), len(adata_b01)
        r = 1 - (2 * u_stat) / (n1 * n2)

        immune_comparison.append({
            'cell_type': cell_type,
            'A01_mean': immune_a01[cell_type]['mean_score'],
            'B01_mean': immune_b01[cell_type]['mean_score'],
            'fold_change': immune_b01[cell_type]['mean_score'] / immune_a01[cell_type]['mean_score']
                          if immune_a01[cell_type]['mean_score'] != 0 else np.inf,
            'U_statistic': u_stat,
            'p_value': p_val,
            'effect_size_r': r
        })

immune_df = pd.DataFrame(immune_comparison)
# Multiple testing correction
immune_df['p_adj'] = np.minimum(immune_df['p_value'] * len(immune_df), 1.0)
immune_df = immune_df.sort_values('p_value')

print(immune_df.to_string(index=False))

# Significant differences
sig_immune = immune_df[immune_df['p_adj'] < 0.05]
if len(sig_immune) > 0:
    print(f"\nSIGNIFICANT IMMUNE DIFFERENCES (p_adj < 0.05):")
    for _, row in sig_immune.iterrows():
        higher = "A01" if row['A01_mean'] > row['B01_mean'] else "B01"
        print(f"  {row['cell_type']}: Higher in {higher} (FC={row['fold_change']:.2f}, p={row['p_adj']:.2e})")

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("SUMMARY OF STATISTICAL FINDINGS")
print("=" * 70)

print("""
KEY STATISTICAL RESULTS:

1. CELL TYPE COMPOSITION:
   - Chi-square test assesses overall composition differences
   - Fisher's exact tests identify which cell types differ
   - Bonferroni correction applied for multiple testing

2. PROTEIN-RNA CORRELATIONS:
   - Spearman correlation used (more robust to outliers)
   - Expected positive correlation if protein reflects RNA expression
   - Low correlation may indicate post-transcriptional regulation

3. SPATIALLY VARIABLE GENES:
   - Moran's I > 0 indicates spatial clustering (autocorrelation)
   - P-values from permutation testing
   - Shared vs sample-specific SVGs indicate biological differences

4. NEIGHBORHOOD ENRICHMENT:
   - Z-scores indicate spatial co-localization patterns
   - Positive Z: cell types cluster together (attraction)
   - Negative Z: cell types avoid each other (repulsion)

5. IMMUNE INFILTRATION:
   - Mann-Whitney U test compares score distributions
   - Effect size (rank-biserial r) indicates magnitude
   - Multiple testing correction applied
""")

# Save results
output_dir = "/home/user/spatial-hackathon-2026/outputs/statistics/"
import os
os.makedirs(output_dir, exist_ok=True)

if cell_type_col and 'fisher_df' in dir():
    fisher_df.to_csv(f"{output_dir}/g4x_cell_type_composition.csv", index=False)
if 'immune_df' in dir():
    immune_df.to_csv(f"{output_dir}/g4x_immune_comparison.csv", index=False)
if svg_a01 is not None:
    svg_a01.to_csv(f"{output_dir}/g4x_A01_morans_i.csv")
if svg_b01 is not None:
    svg_b01.to_csv(f"{output_dir}/g4x_B01_morans_i.csv")

print(f"\nResults saved to: {output_dir}")
print("=" * 70)
