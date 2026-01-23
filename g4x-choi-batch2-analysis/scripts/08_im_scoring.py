#!/usr/bin/env python3
"""
IM (Intestinal Metaplasia) and Hwang Cancer Scoring for G4X Data

Computes:
1. IM Intestinal Score - markers of intestinal-type metaplasia
2. IM Gastric Score - markers of normal gastric mucosa
3. IM Stem/Repair Score - stem cell and tissue repair markers
4. IM Ratio - intestinal vs gastric balance
5. Hwang Cancer Score - 32-gene cancer signature (29/32 in panel)

Critical: Uses layers['counts'] then normalizes for accurate scoring
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
    'output_dir': 'output',
}

# Gene signatures
GENE_SETS = {
    # IM (Intestinal Metaplasia) Markers
    'im_intestinal': ['CDX2', 'CDX1', 'MUC2', 'TFF3', 'CDH17', 'GUCY2C',
                      'GPA33', 'CFTR', 'HNF4A'],

    'im_gastric': ['TFF1', 'TFF2', 'MUC5AC', 'MUC6', 'PGC', 'ATP4A', 'PSCA'],

    'im_stem_repair': ['LGR5', 'OLFM4', 'SOX9'],

    # Hwang 32-gene Cancer Signature (from Hwang lab papers)
    # 29/32 present in G4X 337-gene panel
    'hwang_cancer': [
        'ACTA2', 'AREG', 'ASCC2', 'BEST1', 'BRCA1', 'CREBBP', 'DDX5', 'EP300',
        'ESR1', 'FHL2', 'GNL3', 'HIPK2', 'HSF1', 'IGSF9', 'JUN', 'MSH6',
        'PARP1', 'PAWR', 'PCNA', 'PML', 'PPP2R5A', 'RPA1', 'SMAD3', 'SMARCA4',
        'TP53', 'TP63', 'WRN', 'WT1', 'WTAP'
    ],
}

# ============================================================================
# Helper Functions
# ============================================================================

def check_gene_presence(adata, gene_set, set_name):
    """Check which genes from a gene set are present in the data."""
    present = [g for g in gene_set if g in adata.var_names]
    missing = [g for g in gene_set if g not in adata.var_names]

    print(f"  {set_name}: {len(present)}/{len(gene_set)} genes present")
    if missing:
        print(f"    Missing: {missing[:5]}{'...' if len(missing) > 5 else ''}")

    return present, missing


def score_genes_safe(adata, gene_list, score_name, use_raw=False):
    """
    Safely score genes, handling cases where genes are missing.
    Uses sc.tl.score_genes which expects log-normalized data.
    """
    present_genes = [g for g in gene_list if g in adata.var_names]

    if len(present_genes) == 0:
        print(f"  Warning: No genes found for {score_name}, setting to NaN")
        adata.obs[score_name] = np.nan
        return 0

    if len(present_genes) < len(gene_list):
        print(f"  {score_name}: Using {len(present_genes)}/{len(gene_list)} genes")

    try:
        sc.tl.score_genes(
            adata,
            gene_list=present_genes,
            score_name=score_name,
            use_raw=use_raw
        )
        return len(present_genes)
    except Exception as e:
        print(f"  Error scoring {score_name}: {e}")
        adata.obs[score_name] = np.nan
        return 0


# ============================================================================
# Main Pipeline
# ============================================================================

def main():
    print("=" * 60)
    print("IM + Hwang Scoring Pipeline")
    print("=" * 60)

    # Create output directory
    output_dir = Path(CONFIG['output_dir'])
    output_dir.mkdir(exist_ok=True)

    # Load data
    print("\n--- Loading G4X Data ---")
    adata = sc.read_h5ad(CONFIG['g4x_path'])
    print(f"Shape: {adata.shape}")
    print(f"Samples: {adata.obs['sample_id'].value_counts().to_dict()}")
    print(f"Layers: {list(adata.layers.keys())}")

    # CRITICAL: Work on a copy with counts layer, then normalize
    print("\n--- Preparing for Scoring ---")
    print("Using layers['counts'] and normalizing...")

    adata_work = adata.copy()

    # Use raw counts from layers
    if 'counts' in adata_work.layers:
        adata_work.X = adata_work.layers['counts'].copy()
    else:
        print("Warning: No 'counts' layer found, using X as-is")

    # Normalize (sc.tl.score_genes expects log-normalized data)
    sc.pp.normalize_total(adata_work, target_sum=1e4)
    sc.pp.log1p(adata_work)

    # Check gene presence
    print("\n--- Checking Gene Presence ---")
    gene_stats = {}

    for set_name, genes in GENE_SETS.items():
        present, missing = check_gene_presence(adata_work, genes, set_name)
        gene_stats[set_name] = {
            'total': len(genes),
            'present': len(present),
            'present_genes': present,
            'missing_genes': missing
        }

    # Compute scores
    print("\n--- Computing Scores ---")

    # IM scores
    score_genes_safe(adata_work, GENE_SETS['im_intestinal'], 'im_intestinal')
    score_genes_safe(adata_work, GENE_SETS['im_gastric'], 'im_gastric')
    score_genes_safe(adata_work, GENE_SETS['im_stem_repair'], 'im_stem_repair')

    # Hwang cancer score
    score_genes_safe(adata_work, GENE_SETS['hwang_cancer'], 'hwang_cancer')

    # Compute IM ratio (intestinal / (intestinal + gastric))
    print("  Computing IM ratio...")
    im_int = adata_work.obs['im_intestinal'].values
    im_gas = adata_work.obs['im_gastric'].values

    # Shift to positive (scores can be negative)
    im_int_shifted = im_int - np.nanmin(im_int) + 0.01
    im_gas_shifted = im_gas - np.nanmin(im_gas) + 0.01

    adata_work.obs['im_ratio'] = im_int_shifted / (im_int_shifted + im_gas_shifted)

    # Copy scores back to original adata
    score_cols = ['im_intestinal', 'im_gastric', 'im_stem_repair',
                  'hwang_cancer', 'im_ratio']

    for col in score_cols:
        adata.obs[col] = adata_work.obs[col].values

    # ============================================================================
    # Summary Statistics
    # ============================================================================

    print("\n" + "=" * 60)
    print("Score Summary")
    print("=" * 60)

    # Overall stats
    print("\n--- Overall Score Statistics ---")
    for col in score_cols:
        vals = adata.obs[col].dropna()
        if len(vals) > 0:
            print(f"  {col}: mean={vals.mean():.3f}, std={vals.std():.3f}, "
                  f"range=[{vals.min():.3f}, {vals.max():.3f}]")

    # Per-sample stats
    print("\n--- Per-Sample Mean Scores ---")
    sample_stats = adata.obs.groupby('sample_id')[score_cols].mean()
    print(sample_stats.round(3).to_string())

    # Hwang score validation (should be higher in Cancer)
    print("\n--- Hwang Cancer Score by Stage ---")
    stage_hwang = adata.obs.groupby('stage')['hwang_cancer'].agg(['mean', 'std', 'count'])
    print(stage_hwang.round(3).to_string())

    # IM ratio by stage
    print("\n--- IM Ratio by Stage ---")
    stage_im = adata.obs.groupby('stage')['im_ratio'].agg(['mean', 'std', 'count'])
    print(stage_im.round(3).to_string())

    # ============================================================================
    # Save Results
    # ============================================================================

    print("\n--- Saving Results ---")

    # Save scores as CSV
    scores_df = adata.obs[['sample_id', 'stage'] + score_cols].copy()
    scores_df['cell_id'] = adata.obs_names
    scores_df = scores_df[['cell_id', 'sample_id', 'stage'] + score_cols]

    scores_path = output_dir / 'g4x_im_scores.csv'
    scores_df.to_csv(scores_path, index=False)
    print(f"Scores saved to: {scores_path}")

    # Save gene stats
    gene_stats_path = output_dir / 'im_gene_stats.json'
    with open(gene_stats_path, 'w') as f:
        json.dump(gene_stats, f, indent=2)
    print(f"Gene stats saved to: {gene_stats_path}")

    # Save summary stats
    summary = {
        'n_cells': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'samples': adata.obs['sample_id'].value_counts().to_dict(),
        'stages': adata.obs['stage'].value_counts().to_dict(),
        'score_means': {col: float(adata.obs[col].mean()) for col in score_cols},
        'score_by_stage': {
            stage: {col: float(adata.obs[adata.obs['stage'] == stage][col].mean())
                   for col in score_cols}
            for stage in adata.obs['stage'].unique()
        }
    }

    summary_path = output_dir / 'im_scoring_summary.json'
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Summary saved to: {summary_path}")

    print("\n" + "=" * 60)
    print("IM + Hwang Scoring Complete")
    print("=" * 60)

    return adata


if __name__ == '__main__':
    main()
