#!/usr/bin/env python3
"""
Merge RCTD Deconvolution Results with IM Scoring

Combines:
- RCTD cell type annotations
- RCTD weights (cell type proportions)
- IM and Hwang scores

Performs validation and cross-tabulation.
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import json
import warnings

warnings.filterwarnings('ignore')

# ============================================================================
# Configuration
# ============================================================================

CONFIG = {
    'g4x_path': 'results/pilot/merged_pilot.h5ad',
    'rctd_annotations_path': 'output/g4x_rctd_annotations.csv',
    'rctd_weights_path': 'output/g4x_rctd_weights.csv',
    'im_scores_path': 'output/g4x_im_scores.csv',
    'output_dir': 'output',
    'output_h5ad': 'output/pilot_final_annotated.h5ad',
}

# ============================================================================
# Main Pipeline
# ============================================================================

def main():
    print("=" * 60)
    print("Merge RCTD + IM Results")
    print("=" * 60)

    output_dir = Path(CONFIG['output_dir'])

    # Load original data
    print("\n--- Loading G4X Data ---")
    adata = sc.read_h5ad(CONFIG['g4x_path'])
    print(f"Shape: {adata.shape}")
    print(f"Original obs columns: {list(adata.obs.columns)[:10]}...")

    # Load RCTD annotations
    print("\n--- Loading RCTD Annotations ---")
    rctd_path = Path(CONFIG['rctd_annotations_path'])

    if not rctd_path.exists():
        print(f"Warning: RCTD annotations not found at {rctd_path}")
        print("Run 07_rctd_deconvolution.R first!")
        rctd_df = None
    else:
        rctd_df = pd.read_csv(rctd_path)
        print(f"RCTD annotations: {len(rctd_df)} cells")
        print(f"Columns: {list(rctd_df.columns)}")

        # Merge RCTD annotations
        rctd_df = rctd_df.set_index('cell_id')

        # Add RCTD columns to adata with NA for missing cells
        for col in rctd_df.columns:
            if col in ['cell_type_rctd', 'second_type', 'spot_class']:
                # Reindex to match adata, filling missing with NA
                adata.obs[col] = rctd_df[col].reindex(adata.obs_names).values

        n_annotated = adata.obs['cell_type_rctd'].notna().sum()
        print(f"RCTD annotated: {n_annotated:,} / {adata.n_obs:,} cells ({100*n_annotated/adata.n_obs:.1f}%)")
        print(f"RCTD cell types: {adata.obs['cell_type_rctd'].value_counts().to_dict()}")

    # Load RCTD weights
    print("\n--- Loading RCTD Weights ---")
    weights_path = Path(CONFIG['rctd_weights_path'])

    if not weights_path.exists():
        print(f"Warning: RCTD weights not found at {weights_path}")
    else:
        weights_df = pd.read_csv(weights_path)
        print(f"RCTD weights: {weights_df.shape}")

        # Add weights to obsm
        weights_df = weights_df.set_index('cell_id')
        weight_cols = [c for c in weights_df.columns if c != 'cell_id']

        # Reindex to match adata, filling missing with 0
        weights_matrix = weights_df[weight_cols].reindex(adata.obs_names).fillna(0).values
        adata.obsm['rctd_weights'] = weights_matrix

        # Store weight column names
        adata.uns['rctd_cell_types'] = weight_cols

        print(f"Weight columns: {weight_cols[:5]}...")

    # Load IM scores
    print("\n--- Loading IM Scores ---")
    im_path = Path(CONFIG['im_scores_path'])

    if not im_path.exists():
        print(f"Warning: IM scores not found at {im_path}")
        print("Run 08_im_scoring.py first!")
    else:
        im_df = pd.read_csv(im_path)
        print(f"IM scores: {len(im_df)} cells")

        # Merge IM scores
        im_df = im_df.set_index('cell_id')

        score_cols = ['im_intestinal', 'im_gastric', 'im_stem_repair',
                     'hwang_cancer', 'im_ratio']

        for col in score_cols:
            if col in im_df.columns:
                adata.obs[col] = im_df.loc[adata.obs_names, col].values

        print(f"Added IM score columns: {score_cols}")

    # ============================================================================
    # Cross-Validation
    # ============================================================================

    if rctd_df is not None and 'im_intestinal' in adata.obs:
        print("\n" + "=" * 60)
        print("Cross-Validation: RCTD vs IM")
        print("=" * 60)

        # Dynamic intestinal types detection
        all_types = set(adata.obs["cell_type_rctd"].dropna().unique())
        intestinal_patterns = ['entero', 'goblet', 'intestin', 'absorptive']
        intestinal_types = [t for t in all_types if any(
            x in t.lower() for x in intestinal_patterns
        )]

        print(f"\nDetected intestinal types: {intestinal_types}")

        # IM score by RCTD cell type
        print("\n--- IM Scores by RCTD Cell Type ---")
        if 'cell_type_rctd' in adata.obs:
            im_by_type = adata.obs.groupby('cell_type_rctd')['im_intestinal'].agg(
                ['mean', 'std', 'count']
            ).sort_values('mean', ascending=False)
            print(im_by_type.round(3).head(10).to_string())

        # Hwang score by RCTD cell type
        print("\n--- Hwang Cancer Score by RCTD Cell Type ---")
        if 'cell_type_rctd' in adata.obs:
            hwang_by_type = adata.obs.groupby('cell_type_rctd')['hwang_cancer'].agg(
                ['mean', 'std', 'count']
            ).sort_values('mean', ascending=False)
            print(hwang_by_type.round(3).head(10).to_string())

        # Cross-tab: RCTD intestinal types vs high IM
        print("\n--- Cross-Tab: RCTD Intestinal Types vs High IM Score ---")
        im_threshold = adata.obs['im_intestinal'].quantile(0.75)
        adata.obs['high_im'] = adata.obs['im_intestinal'] > im_threshold

        if len(intestinal_types) > 0:
            adata.obs['is_intestinal_rctd'] = adata.obs['cell_type_rctd'].isin(intestinal_types)
            crosstab = pd.crosstab(
                adata.obs['is_intestinal_rctd'],
                adata.obs['high_im'],
                margins=True
            )
            print(crosstab.to_string())

            # Agreement rate
            agreement = ((adata.obs['is_intestinal_rctd'] == adata.obs['high_im'])).mean()
            print(f"\nRCTD-IM Agreement Rate: {agreement:.1%}")

    # ============================================================================
    # Per-Sample Summary
    # ============================================================================

    print("\n" + "=" * 60)
    print("Per-Sample Summary")
    print("=" * 60)

    # Cell type distribution per sample
    if 'cell_type_rctd' in adata.obs:
        print("\n--- Cell Type Distribution by Sample ---")
        crosstab = pd.crosstab(
            adata.obs['sample_id'],
            adata.obs['cell_type_rctd'],
            normalize='index'
        ).round(3) * 100

        # Show top 5 cell types per sample
        for sample in crosstab.index:
            top_types = crosstab.loc[sample].nlargest(5)
            print(f"\n{sample}:")
            for ct, pct in top_types.items():
                print(f"  {ct}: {pct:.1f}%")

    # ============================================================================
    # Save Final Annotated Data
    # ============================================================================

    print("\n--- Saving Final Annotated Data ---")

    # Save h5ad
    adata.write_h5ad(CONFIG['output_h5ad'])
    print(f"Saved: {CONFIG['output_h5ad']}")

    # Save summary CSV
    summary_cols = ['sample_id', 'stage']
    if 'cell_type_rctd' in adata.obs:
        summary_cols.extend(['cell_type_rctd', 'spot_class'])
    if 'im_intestinal' in adata.obs:
        summary_cols.extend(['im_intestinal', 'im_gastric', 'im_stem_repair',
                            'hwang_cancer', 'im_ratio'])

    summary_df = adata.obs[summary_cols].copy()
    summary_df['cell_id'] = adata.obs_names
    summary_df = summary_df[['cell_id'] + summary_cols]

    summary_path = output_dir / 'pilot_final_summary.csv'
    summary_df.to_csv(summary_path, index=False)
    print(f"Saved: {summary_path}")

    # Save integration stats
    stats = {
        'n_cells': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'samples': adata.obs['sample_id'].value_counts().to_dict(),
    }

    if 'cell_type_rctd' in adata.obs:
        stats['rctd_types'] = adata.obs['cell_type_rctd'].value_counts().to_dict()
        stats['rctd_rejection_rate'] = float(
            (adata.obs['spot_class'] == 'reject').sum() / adata.n_obs
        ) if 'spot_class' in adata.obs else None

    if 'im_intestinal' in adata.obs:
        stats['im_means_by_stage'] = {
            stage: {
                'im_intestinal': float(adata.obs[adata.obs['stage'] == stage]['im_intestinal'].mean()),
                'hwang_cancer': float(adata.obs[adata.obs['stage'] == stage]['hwang_cancer'].mean()),
            }
            for stage in adata.obs['stage'].unique()
        }

    stats_path = output_dir / 'merge_stats.json'
    with open(stats_path, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"Saved: {stats_path}")

    print("\n" + "=" * 60)
    print("Merge Complete")
    print("=" * 60)

    # Final obs columns
    print(f"\nFinal obs columns: {list(adata.obs.columns)}")
    print(f"Final obsm keys: {list(adata.obsm.keys())}")

    return adata


if __name__ == '__main__':
    main()
