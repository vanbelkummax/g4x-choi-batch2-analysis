#!/usr/bin/env python3
"""
Script 24: Multimodal Spatial Entropy Analysis
===============================================
Calculates and compares entropy metrics across R vs NR samples:
1. Cell type composition entropy (from deconvolution)
2. Gene expression entropy (from HVGs)
3. Combined multimodal entropy score

NOVEL: Integrates multiple modalities to quantify TME complexity.
"""

import sys
sys.path.insert(0, '/home/user/spatial-hackathon-2026')

import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import entropy, ttest_ind, mannwhitneyu
from pathlib import Path
import json

from scripts.config import POLYMATHIC_DIR, PDAC_METADATA

# Configuration
OUTPUT_DIR = Path("/home/user/spatial-hackathon-2026/outputs/entropy")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DECONV_DIR = Path("/home/user/spatial-hackathon-2026/outputs/deconvolution")


def load_deconv_results():
    """Load deconvolution results from previous analysis."""
    results_file = DECONV_DIR / "c2l_all_results.json"
    if results_file.exists():
        with open(results_file) as f:
            return json.load(f)
    return {}


def calculate_cell_type_entropy(proportions: dict) -> float:
    """
    Shannon entropy of cell type proportions.
    Higher entropy = more diverse TME composition.
    """
    props = np.array(list(proportions.values()))
    # Shift to positive for entropy calculation
    props = props - props.min() + 1e-10
    props = props / props.sum()
    return entropy(props)


def calculate_expression_entropy(adata: sc.AnnData, n_top_genes: int = 50) -> float:
    """
    Calculate gene expression entropy from top HVGs.
    Uses variance-based entropy as a measure of expression diversity.
    """
    # Get HVGs if available
    if 'highly_variable' in adata.var.columns:
        hvg_mask = adata.var['highly_variable']
        hvg_expr = adata[:, hvg_mask].X
    else:
        # Fallback: use top variable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=min(n_top_genes, adata.n_vars),
                                     flavor='seurat', inplace=False)
        hvg_expr = adata.X

    # Convert to dense if sparse
    if hasattr(hvg_expr, 'toarray'):
        hvg_expr = hvg_expr.toarray()

    # Calculate per-gene variance as pseudo-probability
    gene_vars = np.var(hvg_expr, axis=0)
    gene_vars = gene_vars[gene_vars > 0]  # Remove zeros

    if len(gene_vars) == 0:
        return 0.0

    # Normalize to probability distribution
    gene_probs = gene_vars / gene_vars.sum()
    return entropy(gene_probs)


def calculate_spot_level_entropy(adata: sc.AnnData, ct_col: str = 'cell_type') -> float:
    """
    Calculate spatial heterogeneity from spot-level cell type annotations.
    """
    if ct_col not in adata.obs.columns:
        return np.nan

    ct_counts = adata.obs[ct_col].value_counts()
    ct_probs = ct_counts / ct_counts.sum()
    return entropy(ct_probs)


def calculate_multimodal_entropy(
    sample_id: str,
    deconv_props: dict,
    adata: sc.AnnData = None,
    weights: tuple = (0.5, 0.5)  # (cell_type_weight, expression_weight)
) -> dict:
    """
    Calculate multimodal entropy combining multiple sources.

    Returns dict with individual and combined entropy scores.
    """
    results = {'sample_id': sample_id}

    # 1. Cell type composition entropy
    ct_entropy = calculate_cell_type_entropy(deconv_props)
    results['cell_type_entropy'] = ct_entropy

    # 2. Gene expression entropy (if adata available)
    if adata is not None:
        expr_entropy = calculate_expression_entropy(adata)
        results['expression_entropy'] = expr_entropy

        # 3. Spot-level entropy
        spot_entropy = calculate_spot_level_entropy(adata)
        results['spot_entropy'] = spot_entropy
    else:
        expr_entropy = np.nan
        results['expression_entropy'] = np.nan
        results['spot_entropy'] = np.nan

    # 4. Combined multimodal entropy
    if not np.isnan(expr_entropy):
        combined = weights[0] * ct_entropy + weights[1] * expr_entropy
    else:
        combined = ct_entropy

    results['combined_entropy'] = combined

    return results


def compare_entropy_R_vs_NR(entropy_results: list) -> dict:
    """
    Statistical comparison of entropy metrics between R and NR groups.
    """
    df = pd.DataFrame(entropy_results)

    # Add response status
    df['response'] = df['sample_id'].map(lambda x: PDAC_METADATA.get(x, {}).get('response', 'unknown'))

    R_data = df[df['response'] == 'R']
    NR_data = df[df['response'] == 'NR']

    if len(R_data) < 2 or len(NR_data) < 2:
        return {'error': 'Insufficient samples for comparison'}

    comparisons = {}

    for metric in ['cell_type_entropy', 'expression_entropy', 'combined_entropy']:
        R_vals = R_data[metric].dropna().values
        NR_vals = NR_data[metric].dropna().values

        if len(R_vals) < 2 or len(NR_vals) < 2:
            continue

        # t-test
        t_stat, t_pval = ttest_ind(R_vals, NR_vals)

        # Mann-Whitney U test
        u_stat, mwu_pval = mannwhitneyu(R_vals, NR_vals, alternative='two-sided')

        # Effect size (Cohen's d)
        pooled_std = np.sqrt((np.var(R_vals) + np.var(NR_vals)) / 2)
        cohens_d = (np.mean(R_vals) - np.mean(NR_vals)) / pooled_std if pooled_std > 0 else 0

        comparisons[metric] = {
            'R_mean': float(np.mean(R_vals)),
            'R_std': float(np.std(R_vals)),
            'NR_mean': float(np.mean(NR_vals)),
            'NR_std': float(np.std(NR_vals)),
            't_statistic': float(t_stat),
            't_pvalue': float(t_pval),
            'mwu_statistic': float(u_stat),
            'mwu_pvalue': float(mwu_pval),
            'cohens_d': float(cohens_d),
            'direction': 'R > NR' if cohens_d > 0 else 'NR > R'
        }

    return comparisons


def main():
    print("=" * 60)
    print("Multimodal Spatial Entropy Analysis")
    print("=" * 60)

    # Load deconvolution results
    deconv_results = load_deconv_results()

    if not deconv_results:
        print("No deconvolution results found. Run script 20 first.")
        return

    print(f"Loaded deconvolution for {len(deconv_results)} samples")

    # Calculate entropy for each sample
    all_results = []

    for sample_id, sample_data in deconv_results.items():
        if sample_data.get('status') != 'success':
            print(f"Skipping {sample_id} (failed deconvolution)")
            continue

        print(f"\nProcessing {sample_id}...")

        # Load adata if available
        adata_path = POLYMATHIC_DIR / f"{sample_id}_polymathic.h5ad"
        adata = None

        if adata_path.exists():
            try:
                adata = sc.read_h5ad(adata_path)
                print(f"  Loaded adata: {adata.n_obs} spots, {adata.n_vars} genes")
            except Exception as e:
                print(f"  Could not load adata: {e}")

        # Calculate multimodal entropy
        props = sample_data.get('mean_proportions', {})
        result = calculate_multimodal_entropy(sample_id, props, adata)
        result['response'] = PDAC_METADATA.get(sample_id, {}).get('response', 'unknown')

        all_results.append(result)

        print(f"  Cell type entropy: {result['cell_type_entropy']:.4f}")
        print(f"  Expression entropy: {result.get('expression_entropy', np.nan):.4f}")
        print(f"  Combined entropy: {result['combined_entropy']:.4f}")

    # Save individual results
    df = pd.DataFrame(all_results)
    df.to_csv(OUTPUT_DIR / "multimodal_entropy_results.csv", index=False)
    print(f"\nSaved results to {OUTPUT_DIR / 'multimodal_entropy_results.csv'}")

    # Statistical comparison
    print("\n" + "=" * 60)
    print("R vs NR Entropy Comparison")
    print("=" * 60)

    comparisons = compare_entropy_R_vs_NR(all_results)

    if 'error' in comparisons:
        print(f"Error: {comparisons['error']}")
    else:
        for metric, stats in comparisons.items():
            print(f"\n{metric}:")
            print(f"  R: {stats['R_mean']:.4f} ± {stats['R_std']:.4f}")
            print(f"  NR: {stats['NR_mean']:.4f} ± {stats['NR_std']:.4f}")
            print(f"  Cohen's d: {stats['cohens_d']:.3f} ({stats['direction']})")
            print(f"  p-value (t-test): {stats['t_pvalue']:.4f}")
            print(f"  p-value (MWU): {stats['mwu_pvalue']:.4f}")

        # Save comparison results
        with open(OUTPUT_DIR / "entropy_comparison.json", 'w') as f:
            json.dump(comparisons, f, indent=2)

        print(f"\nSaved comparison to {OUTPUT_DIR / 'entropy_comparison.json'}")

    # Summary table
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    summary_df = df[['sample_id', 'response', 'cell_type_entropy',
                     'expression_entropy', 'combined_entropy']]
    print(summary_df.to_string(index=False))

    return all_results, comparisons


if __name__ == "__main__":
    main()
