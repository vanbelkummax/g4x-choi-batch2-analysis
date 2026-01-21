#!/usr/bin/env python3
"""
G4X Merge and Batch Correction
==============================
Merge all processed samples and apply batch correction if needed.

This script:
1. Merges all QC-passing samples into a single AnnData
2. Computes batch metrics (LISI, silhouette) on merged data
3. Applies Harmony batch correction if metrics fail thresholds
4. Optionally uses scVI for more sophisticated correction
5. Optionally re-computes WNN using muon for true multi-modal integration

Usage:
    conda activate enact
    python scripts/63_merge_and_batch_correct.py [OPTIONS]

Arguments:
    --method        Batch correction method: 'harmony' (default), 'scvi', 'none'
    --force         Apply batch correction even if metrics pass
    --wnn-muon      Use muon for true WNN integration (requires muon installed)
    --batch-key     Column to use for batch correction (default: 'lane')

Output:
    results/qc_all_samples/merged/
    ├── merged_counts.h5ad        # True raw counts (before any normalization)
    ├── merged_normalized.h5ad    # Normalized/log1p data (pre-correction baseline)
    ├── merged_corrected.h5ad     # Batch-corrected data (if applied)
    ├── batch_assessment.csv      # Pre/post correction metrics
    └── MERGE_REPORT.md           # Summary report
"""

import os
import sys
import shutil


def check_environment(skip_check: bool = False):
    """
    Verify running in correct conda environment.

    Parameters
    ----------
    skip_check : bool
        If True, skip the environment check entirely (use --skip-env-check flag)
    """
    if skip_check or os.environ.get("G4X_SKIP_ENV_CHECK"):
        return  # Allow override via flag or env var

    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    env_name = os.path.basename(conda_prefix) if conda_prefix else ""

    # Accept envs containing these patterns (e.g., 'enact-gpu', 'spatial-dev', 'scanpy3.10')
    # Using substring matching to be more flexible
    valid_patterns = ['enact', 'spatial', 'scanpy', 'single-cell', 'scverse']
    is_valid = any(pattern in env_name.lower() for pattern in valid_patterns)

    if not is_valid:
        print("ERROR: This script requires a spatial analysis conda environment.")
        print(f"Current environment: {env_name or 'none'}")
        print(f"\nAccepted environment patterns: {', '.join(valid_patterns)}")
        print("\nOptions:")
        print("  1. Activate a valid env: conda activate enact")
        print("  2. Skip check: --skip-env-check")
        print("  3. Set env var: export G4X_SKIP_ENV_CHECK=1")
        sys.exit(1)


# Defer environment check until after argparse (see main())

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
import argparse
import logging
import warnings
import gc

warnings.filterwarnings('ignore')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# =============================================================================
# Configuration
# =============================================================================

INPUT_DIR = Path("/home/user/g4x-choi-batch2-analysis/results/qc_all_samples/final_processed")
OUTPUT_DIR = Path("/home/user/g4x-choi-batch2-analysis/results/qc_all_samples/merged")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Batch effect thresholds (same as 61)
BATCH_THRESHOLDS = {
    'max_silhouette_batch': 0.3,   # Low = good mixing
    'min_lisi': 3.0,               # High = good mixing
    'max_pc1_lane_variance': 0.20, # < 20% variance explained by lane
}


# =============================================================================
# Batch Assessment Functions
# =============================================================================

def compute_lisi(X: np.ndarray, batch_labels: np.ndarray, perplexity: int = 30) -> np.ndarray:
    """Compute Local Inverse Simpson Index (LISI) for batch mixing."""
    n_neighbors = min(perplexity * 3, len(X) - 1)
    nn = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto')
    nn.fit(X)
    _, indices = nn.kneighbors(X)

    lisi_scores = []
    for neighbors in indices:
        neighbor_batches = batch_labels[neighbors]
        counts = pd.Series(neighbor_batches).value_counts(normalize=True)
        simpson = (counts ** 2).sum()
        lisi = 1.0 / simpson
        lisi_scores.append(lisi)

    return np.array(lisi_scores)


def assess_batch_effects(adata: ad.AnnData, batch_key: str = 'lane') -> dict:
    """
    Compute batch effect metrics on merged data.

    Returns dict with:
        - mean_lisi: Average LISI score
        - silhouette_batch: Silhouette score for batch clustering
        - recommendation: 'PASS' or 'CORRECT'
    """
    logger.info(f"Assessing batch effects using '{batch_key}'...")

    # Ensure PCA exists
    if 'X_pca' not in adata.obsm:
        logger.info("  Computing PCA...")
        sc.pp.pca(adata, n_comps=30)

    batch_labels = adata.obs[batch_key].values
    pca_embedding = adata.obsm['X_pca']

    # Compute LISI
    logger.info("  Computing LISI...")
    lisi_scores = compute_lisi(pca_embedding, batch_labels)
    mean_lisi = np.mean(lisi_scores)

    # Compute silhouette
    logger.info("  Computing silhouette...")
    try:
        sample_size = min(10000, adata.n_obs)
        sil = silhouette_score(
            pca_embedding,
            batch_labels,
            sample_size=sample_size
        )
    except Exception as e:
        logger.warning(f"  Could not compute silhouette: {e}")
        sil = np.nan

    # Determine recommendation
    passes_lisi = mean_lisi >= BATCH_THRESHOLDS['min_lisi']
    passes_sil = sil < BATCH_THRESHOLDS['max_silhouette_batch']

    recommendation = 'PASS' if (passes_lisi and passes_sil) else 'CORRECT'

    metrics = {
        'mean_lisi': mean_lisi,
        'lisi_threshold': BATCH_THRESHOLDS['min_lisi'],
        'lisi_pass': passes_lisi,
        'silhouette_batch': sil,
        'silhouette_threshold': BATCH_THRESHOLDS['max_silhouette_batch'],
        'silhouette_pass': passes_sil,
        'recommendation': recommendation,
    }

    logger.info(f"  LISI: {mean_lisi:.2f} (threshold > {BATCH_THRESHOLDS['min_lisi']}) - {'PASS' if passes_lisi else 'FAIL'}")
    logger.info(f"  Silhouette: {sil:.3f} (threshold < {BATCH_THRESHOLDS['max_silhouette_batch']}) - {'PASS' if passes_sil else 'FAIL'}")
    logger.info(f"  Recommendation: {recommendation}")

    return metrics


# =============================================================================
# Batch Correction Methods
# =============================================================================

def apply_harmony(adata: ad.AnnData, batch_key: str = 'lane') -> ad.AnnData:
    """
    Apply Harmony batch correction.

    Harmony corrects the PCA embedding, preserving biological variation
    while removing batch effects.
    """
    logger.info("Applying Harmony batch correction...")

    try:
        import harmonypy
    except ImportError:
        logger.error("harmonypy not installed. Run: pip install harmonypy")
        return adata

    # Ensure PCA exists
    if 'X_pca' not in adata.obsm:
        sc.pp.pca(adata, n_comps=30)

    # Store original PCA
    adata.obsm['X_pca_uncorrected'] = adata.obsm['X_pca'].copy()

    # Run Harmony
    ho = harmonypy.run_harmony(
        adata.obsm['X_pca'],
        adata.obs,
        batch_key,
        max_iter_harmony=20
    )

    # Store corrected embedding
    adata.obsm['X_pca'] = ho.Z_corr.T
    adata.obsm['X_pca_harmony'] = ho.Z_corr.T

    # Recompute neighbors and UMAP on corrected embedding
    logger.info("  Recomputing neighbors and UMAP...")
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata)

    adata.uns['batch_correction'] = {
        'method': 'harmony',
        'batch_key': batch_key,
    }

    logger.info("  Harmony correction complete")
    return adata


def apply_scvi(adata: ad.AnnData, batch_key: str = 'lane') -> tuple[ad.AnnData, str]:
    """
    Apply scVI batch correction.

    scVI learns a latent representation that removes batch effects while
    preserving biological variation. More sophisticated than Harmony but
    requires more compute.

    IMPORTANT: scVI requires raw counts. This function expects raw counts to be
    stored in adata.layers['counts']. If not present, falls back to Harmony.

    Returns
    -------
    tuple[AnnData, str]
        (corrected adata, method_used) where method_used is 'scvi' or 'harmony_fallback'
    """
    logger.info("Applying scVI batch correction...")

    try:
        import scvi
    except ImportError:
        logger.error("scvi-tools not installed. Run: pip install scvi-tools")
        logger.warning("FALLBACK: Using Harmony instead of scVI")
        result = apply_harmony(adata, batch_key)
        return result, 'harmony_fallback_import'

    # Verify raw counts are available
    if 'counts' not in adata.layers:
        logger.error("Raw counts not found in adata.layers['counts']")
        logger.error("scVI requires raw counts, not normalized data!")
        logger.warning("FALLBACK: Using Harmony instead of scVI")
        result = apply_harmony(adata, batch_key)
        return result, 'harmony_fallback_no_counts'

    # Create a copy with raw counts in X for scVI
    logger.info("  Preparing data with raw counts for scVI...")
    adata_scvi = adata.copy()
    adata_scvi.X = adata_scvi.layers['counts'].copy()

    # Setup anndata for scVI (must be raw counts)
    scvi.model.SCVI.setup_anndata(
        adata_scvi,
        batch_key=batch_key,
        layer=None,  # Use X directly (which is now raw counts)
    )

    # Train model
    logger.info("  Training scVI model (this may take a while)...")
    model = scvi.model.SCVI(adata_scvi, n_latent=30, n_layers=2)
    model.train(max_epochs=100, early_stopping=True)

    # Get latent representation and transfer back to original adata
    adata.obsm['X_scVI'] = model.get_latent_representation()

    # Recompute neighbors and UMAP on scVI embedding
    logger.info("  Recomputing neighbors and UMAP...")
    sc.pp.neighbors(adata, use_rep='X_scVI')
    sc.tl.umap(adata)

    adata.uns['batch_correction'] = {
        'method': 'scvi',
        'batch_key': batch_key,
        'n_latent': 30,
        'n_layers': 2,
    }

    logger.info("  scVI correction complete")
    return adata, 'scvi'


def apply_muon_wnn(adata: ad.AnnData) -> ad.AnnData:
    """
    Re-compute WNN using muon for true multi-modal integration.

    This replaces the simplified WNN with Seurat v4-style WNN that:
    1. Computes modality-specific k-NN graphs
    2. Learns per-cell modality weights
    3. Fuses graphs into a unified representation
    """
    logger.info("Computing true WNN with muon...")

    try:
        import muon as mu
        from muon import MuData
    except ImportError:
        logger.error("muon not installed. Run: pip install muon")
        logger.info("Keeping simplified WNN.")
        return adata

    # Check if protein data exists
    if 'protein' not in adata.obsm:
        logger.warning("No protein data found in obsm['protein']. Cannot compute WNN.")
        return adata

    # Create MuData object
    logger.info("  Creating MuData object...")

    # RNA modality
    rna = adata.copy()

    # Protein modality
    protein_names = adata.uns.get('protein_names', [f'protein_{i}' for i in range(adata.obsm['protein'].shape[1])])
    protein = ad.AnnData(
        X=adata.obsm['protein'],
        obs=adata.obs.copy()
    )
    protein.var_names = protein_names

    # Ensure both have embeddings
    if 'X_pca' not in rna.obsm:
        sc.pp.pca(rna, n_comps=30)

    sc.pp.pca(protein, n_comps=min(15, protein.n_vars - 1))

    # Create MuData
    mdata = MuData({'rna': rna, 'protein': protein})

    # Compute WNN
    logger.info("  Computing WNN neighbors...")
    mu.pp.neighbors(mdata, method='wnn')

    # Extract WNN connectivities back to original AnnData
    adata.obsp['connectivities'] = mdata.obsp['connectivities']
    adata.obsp['distances'] = mdata.obsp['distances']
    adata.uns['neighbors'] = mdata.uns['neighbors']

    # Compute UMAP on WNN graph
    logger.info("  Computing UMAP...")
    sc.tl.umap(adata)

    adata.uns['wnn_method'] = 'muon'

    logger.info("  muon WNN complete")
    return adata


# =============================================================================
# Report Generation
# =============================================================================

def generate_merge_report(
    pre_metrics: dict,
    post_metrics: dict | None,
    n_samples: int,
    n_cells: int,
    method: str,
    save_path: Path
):
    """Generate merge and batch correction report."""

    report = f"""# G4X Merge and Batch Correction Report

**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}

## Summary

| Metric | Value |
|--------|-------|
| Samples merged | {n_samples} |
| Total cells | {n_cells:,} |
| Batch correction method | {method} |

## Pre-Correction Batch Assessment

| Metric | Value | Threshold | Status |
|--------|-------|-----------|--------|
| Mean LISI | {pre_metrics['mean_lisi']:.2f} | > {pre_metrics['lisi_threshold']} | {'PASS' if pre_metrics['lisi_pass'] else 'FAIL'} |
| Silhouette | {pre_metrics['silhouette_batch']:.3f} | < {pre_metrics['silhouette_threshold']} | {'PASS' if pre_metrics['silhouette_pass'] else 'FAIL'} |
| **Recommendation** | {pre_metrics['recommendation']} | | |

"""

    if post_metrics:
        report += f"""## Post-Correction Batch Assessment

| Metric | Value | Threshold | Status |
|--------|-------|-----------|--------|
| Mean LISI | {post_metrics['mean_lisi']:.2f} | > {post_metrics['lisi_threshold']} | {'PASS' if post_metrics['lisi_pass'] else 'FAIL'} |
| Silhouette | {post_metrics['silhouette_batch']:.3f} | < {post_metrics['silhouette_threshold']} | {'PASS' if post_metrics['silhouette_pass'] else 'FAIL'} |

## Improvement

| Metric | Before | After | Change |
|--------|--------|-------|--------|
| LISI | {pre_metrics['mean_lisi']:.2f} | {post_metrics['mean_lisi']:.2f} | {post_metrics['mean_lisi'] - pre_metrics['mean_lisi']:+.2f} |
| Silhouette | {pre_metrics['silhouette_batch']:.3f} | {post_metrics['silhouette_batch']:.3f} | {post_metrics['silhouette_batch'] - pre_metrics['silhouette_batch']:+.3f} |

"""
    else:
        report += """## Batch Correction

No batch correction applied (metrics passed thresholds or --method none specified).

"""

    report += """## Output Files

- `merged_counts.h5ad`: **True raw counts** (before ANY normalization) - use for scVI
- `merged_normalized.h5ad`: Normalized/log1p data (pre-correction baseline)
- `merged_corrected.h5ad`: Batch-corrected data (if applied)
- `batch_assessment.csv`: Pre/post correction metrics

## Next Steps

1. Review batch effect visualizations in `figures/batch_effects/`
2. Proceed to downstream analyses with `merged_corrected.h5ad`
3. For scVI: raw counts are available in `merged_counts.h5ad` or `adata.layers['counts']`
"""

    with open(save_path, 'w') as f:
        f.write(report)


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Merge and batch correct G4X samples')
    parser.add_argument('--method', type=str, default='harmony',
                        choices=['harmony', 'scvi', 'none'],
                        help='Batch correction method (default: harmony)')
    parser.add_argument('--force', action='store_true',
                        help='Apply batch correction even if metrics pass')
    parser.add_argument('--wnn-muon', action='store_true',
                        help='Use muon for true WNN integration')
    parser.add_argument('--batch-key', type=str, default='lane',
                        help='Column to use for batch correction (default: lane)')
    parser.add_argument('--skip-env-check', action='store_true',
                        help='Skip conda environment validation (use if running from equivalent env)')
    args = parser.parse_args()

    # Environment check (deferred to allow --skip-env-check flag)
    check_environment(skip_check=args.skip_env_check)

    logger.info("=" * 60)
    logger.info("G4X Merge and Batch Correction")
    logger.info("=" * 60)

    # Find processed samples
    processed_files = list(INPUT_DIR.glob("*_final.h5ad"))
    if not processed_files:
        logger.error(f"No processed files found in {INPUT_DIR}")
        logger.error("Run 62_process_qc_passing.py first")
        return

    logger.info(f"Found {len(processed_files)} processed samples")

    # Load and merge
    logger.info("Loading and merging samples...")
    adatas = []
    for f in processed_files:
        adata = sc.read_h5ad(f)
        adatas.append(adata)
        logger.info(f"  {f.stem}: {adata.n_obs:,} cells")

    # Concatenate
    adata_merged = ad.concat(adatas, join='outer')
    adata_merged.obs_names_make_unique()

    n_samples = len(processed_files)
    n_cells = adata_merged.n_obs

    logger.info(f"Merged: {n_cells:,} cells from {n_samples} samples")

    # ==========================================================================
    # CRITICAL: ALWAYS save raw counts BEFORE any normalization
    # This ensures:
    #   1. merged_counts.h5ad contains exactly what was merged (no transforms)
    #   2. scVI can access raw counts via adata.layers['counts']
    #
    # NOTE: We unconditionally store counts regardless of magnitude. The old
    # heuristic (X.max() > 100) could misclassify low-coverage samples as
    # "already normalized" and break scVI. Now we always preserve the original
    # data and let downstream tools decide what to do with it.
    # ==========================================================================

    # ALWAYS store original data in layers (required for scVI)
    logger.info("Storing original data in .layers['counts'] for scVI compatibility...")
    adata_merged.layers['counts'] = adata_merged.X.copy()

    # ALWAYS save raw/original counts file (before any transformation)
    logger.info("Saving original counts (merged_counts.h5ad)...")
    adata_merged.write(OUTPUT_DIR / "merged_counts.h5ad")

    # Check if normalization is needed
    # Use multiple heuristics: integer-like values, high max, no negative values
    x_max = adata_merged.X.max()
    x_min = adata_merged.X.min()
    is_integer_like = np.allclose(adata_merged.X.data, np.round(adata_merged.X.data)) if hasattr(adata_merged.X, 'data') else np.allclose(adata_merged.X, np.round(adata_merged.X))
    likely_needs_normalization = (x_max > 50) or (is_integer_like and x_min >= 0)

    if likely_needs_normalization:
        logger.info(f"Data appears to be counts (max={x_max:.1f}, integer-like={is_integer_like})")
        logger.info("Normalizing merged data...")
        sc.pp.normalize_total(adata_merged, target_sum=1e4)
        sc.pp.log1p(adata_merged)
    else:
        logger.info(f"Data may already be normalized (max={x_max:.2f}, integer-like={is_integer_like})")
        logger.info("Skipping normalization - data will be used as-is")
        logger.info("NOTE: .layers['counts'] still contains original data for scVI if needed")

    # Variable genes and PCA
    # For targeted panels (<2000 genes), skip HVG selection - use all genes
    n_genes = adata_merged.n_vars
    if n_genes <= 2000:
        logger.info(f"Targeted panel ({n_genes} genes) - using all genes for PCA")
        adata_merged.var['highly_variable'] = True
    else:
        logger.info(f"Computing HVGs from {n_genes} genes...")
        sc.pp.highly_variable_genes(adata_merged, n_top_genes=2000, subset=True)

    # PCA - cap components at n_genes - 1
    n_pcs = min(30, n_genes - 1)
    logger.info(f"Computing PCA with {n_pcs} components...")
    sc.pp.pca(adata_merged, n_comps=n_pcs)
    sc.pp.neighbors(adata_merged)
    sc.tl.umap(adata_merged)

    # Save normalized pre-correction baseline
    logger.info("Saving normalized data (merged_normalized.h5ad)...")
    adata_merged.write(OUTPUT_DIR / "merged_normalized.h5ad")

    # Assess batch effects
    pre_metrics = assess_batch_effects(adata_merged, args.batch_key)

    # Apply batch correction if needed
    post_metrics = None
    method_applied = 'none'
    fallback_reason = None  # Track if scVI fell back to Harmony

    if args.method != 'none' and (pre_metrics['recommendation'] == 'CORRECT' or args.force):
        if args.method == 'harmony':
            adata_merged = apply_harmony(adata_merged, args.batch_key)
            method_applied = 'harmony'
        elif args.method == 'scvi':
            # apply_scvi returns (adata, method_used) to track fallbacks
            adata_merged, actual_method = apply_scvi(adata_merged, args.batch_key)
            method_applied = actual_method

            # Log fallback explicitly
            if actual_method.startswith('harmony_fallback'):
                fallback_reason = actual_method.replace('harmony_fallback_', '')
                logger.warning(f"=" * 40)
                logger.warning(f"scVI FALLBACK OCCURRED")
                logger.warning(f"Reason: {fallback_reason}")
                logger.warning(f"Method used: Harmony (instead of scVI)")
                logger.warning(f"=" * 40)
                method_applied = f'harmony (scVI fallback: {fallback_reason})'

        # Re-assess after correction
        post_metrics = assess_batch_effects(adata_merged, args.batch_key)
    elif args.method == 'none':
        logger.info("Batch correction skipped (--method none)")
        method_applied = 'none'
    else:
        logger.info("Batch metrics pass thresholds - no correction needed")
        logger.info("Use --force to apply correction anyway")
        method_applied = 'none (passed)'

    # Optional: true WNN with muon
    if args.wnn_muon:
        adata_merged = apply_muon_wnn(adata_merged)

    # Save corrected (or just processed) data
    logger.info("Saving final merged data...")
    adata_merged.write(OUTPUT_DIR / "merged_corrected.h5ad")

    # Save metrics
    metrics_data = {'stage': ['pre-correction'], **{k: [v] for k, v in pre_metrics.items()}}
    if post_metrics:
        for k, v in post_metrics.items():
            metrics_data[k].append(v)
        metrics_data['stage'].append('post-correction')

    pd.DataFrame(metrics_data).to_csv(OUTPUT_DIR / "batch_assessment.csv", index=False)

    # Generate report
    generate_merge_report(
        pre_metrics=pre_metrics,
        post_metrics=post_metrics,
        n_samples=n_samples,
        n_cells=n_cells,
        method=method_applied,
        save_path=OUTPUT_DIR / "MERGE_REPORT.md"
    )

    # Summary
    logger.info("=" * 60)
    logger.info("Merge Complete")
    logger.info(f"  Samples: {n_samples}")
    logger.info(f"  Cells: {n_cells:,}")
    logger.info(f"  Method: {method_applied}")
    if post_metrics:
        logger.info(f"  LISI improvement: {pre_metrics['mean_lisi']:.2f} -> {post_metrics['mean_lisi']:.2f}")
    logger.info(f"  Output: {OUTPUT_DIR}")

    # Cleanup
    del adata_merged, adatas
    gc.collect()


if __name__ == "__main__":
    main()
