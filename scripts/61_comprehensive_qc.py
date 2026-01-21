#!/usr/bin/env python3
"""
G4X Comprehensive Quality Control
==================================
Multi-level QC analysis with batch effect assessment.

Usage:
    python scripts/61_comprehensive_qc.py

Output:
    results/qc_all_samples/
    ├── figures/per_sample/      # 32 QC panels
    ├── figures/cross_sample/    # Comparison plots
    ├── figures/batch_effects/   # Lane analysis
    ├── sample_qc_summary.csv    # Pass/fail + reasons
    └── QC_REPORT.md             # Comprehensive report
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from tqdm import tqdm
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

INPUT_DIR = Path("/home/user/g4x-choi-batch2-analysis/results/qc_all_samples/raw")
OUTPUT_DIR = Path("/home/user/g4x-choi-batch2-analysis/results/qc_all_samples")
FIG_DIR = OUTPUT_DIR / "figures"
FIG_DIR_SAMPLE = FIG_DIR / "per_sample"
FIG_DIR_CROSS = FIG_DIR / "cross_sample"
FIG_DIR_BATCH = FIG_DIR / "batch_effects"

for d in [FIG_DIR_SAMPLE, FIG_DIR_CROSS, FIG_DIR_BATCH]:
    d.mkdir(parents=True, exist_ok=True)

# Resolve baseline metrics for validation
RESOLVE_METRICS_PATH = Path("/mnt/x/Choi_Batch_2_Tuesday/choi_preGC_b2_core_metrics.csv")

# QC thresholds (evidence-based)
QC_THRESHOLDS = {
    'min_cells': 20_000,
    'min_median_transcripts_per_cell': 30,
    'min_median_genes_per_cell': 20,
    'max_pct_empty': 5.0,
    'min_pct_in_cells': 80.0,
}

# Batch effect thresholds
BATCH_THRESHOLDS = {
    'max_silhouette_batch': 0.3,  # Low = good mixing
    'min_lisi': 2.0,              # High = good mixing
    'max_pc1_lane_variance': 0.20,  # <20% variance explained by lane
}


# =============================================================================
# QC Functions
# =============================================================================

def compute_sample_qc_metrics(adata: ad.AnnData) -> dict:
    """Compute comprehensive QC metrics for a single sample."""
    metrics = {
        'sample_id': adata.obs['sample_id'].iloc[0],
        'lane': adata.obs['lane'].iloc[0],
        'stage': adata.obs['stage'].iloc[0],
        'patient': adata.obs['patient'].iloc[0],
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'median_counts': np.median(adata.obs['n_counts']),
        'median_genes': np.median(adata.obs['n_genes']),
        'mean_counts': np.mean(adata.obs['n_counts']),
        'std_counts': np.std(adata.obs['n_counts']),
        'pct_zero_count_cells': (adata.obs['n_counts'] == 0).mean() * 100,
        'q25_counts': np.percentile(adata.obs['n_counts'], 25),
        'q75_counts': np.percentile(adata.obs['n_counts'], 75),
    }

    # Cell area if available
    if 'cell_area' in adata.obs.columns:
        metrics['median_cell_area'] = np.median(adata.obs['cell_area'])
        metrics['std_cell_area'] = np.std(adata.obs['cell_area'])

    # Protein metrics
    if 'protein' in adata.obsm:
        protein_sum = adata.obsm['protein'].sum(axis=1)
        metrics['median_protein_counts'] = np.median(protein_sum)
        metrics['pct_protein_positive'] = (protein_sum > 0).mean() * 100

    return metrics


def apply_qc_thresholds(metrics: dict) -> tuple:
    """
    Apply QC thresholds and return pass/fail status with reasons.

    Returns
    -------
    tuple
        (status: 'PASS'|'WARN'|'FAIL', reasons: list of strings)
    """
    reasons = []
    status = 'PASS'

    # Check each threshold
    if metrics['n_cells'] < QC_THRESHOLDS['min_cells']:
        reasons.append(f"Low cell count: {metrics['n_cells']:,} < {QC_THRESHOLDS['min_cells']:,}")
        status = 'WARN'

    if metrics['median_counts'] < QC_THRESHOLDS['min_median_transcripts_per_cell']:
        reasons.append(f"Low median transcripts: {metrics['median_counts']:.1f} < {QC_THRESHOLDS['min_median_transcripts_per_cell']}")
        status = 'FAIL'

    if metrics['median_genes'] < QC_THRESHOLDS['min_median_genes_per_cell']:
        reasons.append(f"Low median genes: {metrics['median_genes']:.1f} < {QC_THRESHOLDS['min_median_genes_per_cell']}")
        status = 'FAIL'

    if metrics['pct_zero_count_cells'] > QC_THRESHOLDS['max_pct_empty']:
        reasons.append(f"High empty cells: {metrics['pct_zero_count_cells']:.1f}% > {QC_THRESHOLDS['max_pct_empty']}%")
        if status != 'FAIL':
            status = 'WARN'

    return status, reasons


def compute_lisi(X: np.ndarray, batch_labels: np.ndarray, perplexity: int = 30) -> np.ndarray:
    """
    Compute Local Inverse Simpson Index (LISI) for batch mixing assessment.

    Higher LISI = better batch mixing.
    LISI of 1 = all neighbors from same batch (bad)
    LISI of n_batches = perfect mixing (good)

    Parameters
    ----------
    X : np.ndarray
        Embedding matrix (n_cells, n_dims)
    batch_labels : np.ndarray
        Batch labels for each cell
    perplexity : int
        Number of neighbors to consider

    Returns
    -------
    np.ndarray
        LISI score for each cell
    """
    n_neighbors = min(perplexity * 3, len(X) - 1)
    nn = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto')
    nn.fit(X)
    _, indices = nn.kneighbors(X)

    lisi_scores = []
    for i, neighbors in enumerate(indices):
        neighbor_batches = batch_labels[neighbors]
        # Simpson index = sum(p^2), LISI = 1/Simpson
        counts = pd.Series(neighbor_batches).value_counts(normalize=True)
        simpson = (counts ** 2).sum()
        lisi = 1.0 / simpson
        lisi_scores.append(lisi)

    return np.array(lisi_scores)


def plot_sample_qc_panel(adata: ad.AnnData, sample_id: str, save_path: Path):
    """Generate 6-panel QC figure for a single sample."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'QC Panel: {sample_id}', fontsize=14, fontweight='bold')

    # 1. Spatial scatter: n_genes
    ax = axes[0, 0]
    if 'x_centroid' in adata.obs.columns and 'y_centroid' in adata.obs.columns:
        scatter = ax.scatter(
            adata.obs['x_centroid'],
            adata.obs['y_centroid'],
            c=adata.obs['n_genes'],
            s=1, cmap='viridis', alpha=0.5
        )
        plt.colorbar(scatter, ax=ax, label='n_genes')
        ax.set_title('Spatial: Genes per cell')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
    else:
        ax.text(0.5, 0.5, 'No spatial coords', ha='center', va='center')
        ax.set_title('Spatial: Genes per cell')

    # 2. Spatial scatter: n_counts
    ax = axes[0, 1]
    if 'x_centroid' in adata.obs.columns and 'y_centroid' in adata.obs.columns:
        scatter = ax.scatter(
            adata.obs['x_centroid'],
            adata.obs['y_centroid'],
            c=adata.obs['n_counts'],
            s=1, cmap='plasma', alpha=0.5
        )
        plt.colorbar(scatter, ax=ax, label='n_counts')
        ax.set_title('Spatial: Counts per cell')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
    else:
        ax.text(0.5, 0.5, 'No spatial coords', ha='center', va='center')
        ax.set_title('Spatial: Counts per cell')

    # 3. Violin: n_genes distribution
    ax = axes[0, 2]
    ax.violinplot(adata.obs['n_genes'], positions=[0], showmedians=True)
    ax.axhline(QC_THRESHOLDS['min_median_genes_per_cell'], color='r',
               linestyle='--', label=f"Threshold: {QC_THRESHOLDS['min_median_genes_per_cell']}")
    ax.set_title('Genes per cell distribution')
    ax.set_ylabel('n_genes')
    ax.legend()

    # 4. Scatter: n_counts vs n_genes
    ax = axes[1, 0]
    ax.scatter(adata.obs['n_counts'], adata.obs['n_genes'], s=1, alpha=0.3)
    ax.set_xlabel('n_counts')
    ax.set_ylabel('n_genes')
    ax.set_title('Counts vs Genes')

    # 5. Histogram: cell area (if available)
    ax = axes[1, 1]
    if 'cell_area' in adata.obs.columns:
        ax.hist(adata.obs['cell_area'], bins=50, edgecolor='black', alpha=0.7)
        ax.axvline(np.median(adata.obs['cell_area']), color='r',
                   linestyle='--', label=f"Median: {np.median(adata.obs['cell_area']):.1f}")
        ax.set_xlabel('Cell area (um2)')
        ax.set_ylabel('Count')
        ax.set_title('Cell area distribution')
        ax.legend()
    else:
        ax.text(0.5, 0.5, 'No cell area data', ha='center', va='center')
        ax.set_title('Cell area distribution')

    # 6. QC summary text
    ax = axes[1, 2]
    ax.axis('off')
    metrics = compute_sample_qc_metrics(adata)
    status, reasons = apply_qc_thresholds(metrics)

    summary_text = f"""
    Sample: {sample_id}
    Lane: {metrics['lane']}
    Stage: {metrics['stage']}

    Cells: {metrics['n_cells']:,}
    Median counts: {metrics['median_counts']:.1f}
    Median genes: {metrics['median_genes']:.1f}
    Empty cells: {metrics['pct_zero_count_cells']:.2f}%

    STATUS: {status}
    """

    if reasons:
        summary_text += "\nIssues:\n" + "\n".join(f"  - {r}" for r in reasons)

    color = {'PASS': 'green', 'WARN': 'orange', 'FAIL': 'red'}[status]
    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor=color, alpha=0.2))

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_cross_sample_comparison(all_metrics: pd.DataFrame, save_dir: Path):
    """Generate cross-sample comparison plots."""

    # 1. Bar plot: median counts by sample
    fig, ax = plt.subplots(figsize=(14, 8))
    stage_colors = {'control': 'gray', 'normal': 'green', 'metaplasia': 'orange', 'cancer': 'red', 'unknown': 'blue'}

    samples = all_metrics.sort_values('median_counts', ascending=True)
    colors = [stage_colors.get(s, 'blue') for s in samples['stage']]

    bars = ax.barh(range(len(samples)), samples['median_counts'], color=colors, alpha=0.7)

    # Mark failures
    for i, (_, row) in enumerate(samples.iterrows()):
        if row['qc_status'] == 'FAIL':
            ax.scatter(row['median_counts'] + 5, i, marker='x', color='red', s=100, zorder=5)

    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(samples['sample_id'])
    ax.set_xlabel('Median counts per cell')
    ax.set_title('Median Counts by Sample (X = QC Fail)')
    ax.axvline(QC_THRESHOLDS['min_median_transcripts_per_cell'], color='red',
               linestyle='--', label=f"Threshold: {QC_THRESHOLDS['min_median_transcripts_per_cell']}")

    # Custom legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=s) for s, c in stage_colors.items() if s != 'unknown']
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()
    plt.savefig(save_dir / "median_counts_by_sample.png", dpi=150)
    plt.close()

    # 2. Box plot: metrics by lane
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for ax, metric, title in zip(
        axes,
        ['median_counts', 'median_genes', 'n_cells'],
        ['Median Counts', 'Median Genes', 'Cell Count']
    ):
        all_metrics.boxplot(column=metric, by='lane', ax=ax)
        ax.set_title(f'{title} by Lane')
        ax.set_xlabel('Lane')

    plt.suptitle('')
    plt.tight_layout()
    plt.savefig(save_dir / "metrics_by_lane.png", dpi=150)
    plt.close()

    # 3. Heatmap: sample metrics
    fig, ax = plt.subplots(figsize=(12, 10))
    pivot = all_metrics.set_index('sample_id')[['median_counts', 'median_genes', 'n_cells', 'pct_zero_count_cells']]
    pivot_norm = (pivot - pivot.mean()) / pivot.std()
    sns.heatmap(pivot_norm, ax=ax, cmap='RdBu_r', center=0, annot=False)
    ax.set_title('Sample QC Metrics (Z-scored)')
    plt.tight_layout()
    plt.savefig(save_dir / "sample_metrics_heatmap.png", dpi=150)
    plt.close()


def plot_batch_effects(adata_combined: ad.AnnData, save_dir: Path):
    """Generate batch effect assessment plots."""

    # 1. UMAP colored by lane vs stage
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    if 'X_umap' in adata_combined.obsm:
        for lane in sorted(adata_combined.obs['lane'].unique()):
            mask = adata_combined.obs['lane'] == lane
            ax.scatter(
                adata_combined.obsm['X_umap'][mask, 0],
                adata_combined.obsm['X_umap'][mask, 1],
                s=1, alpha=0.3, label=lane
            )
        ax.legend(markerscale=5)
        ax.set_title('UMAP colored by Lane')
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')

    ax = axes[1]
    if 'X_umap' in adata_combined.obsm:
        stage_colors = {'control': 'gray', 'normal': 'green', 'metaplasia': 'orange', 'cancer': 'red'}
        for stage in sorted(adata_combined.obs['stage'].unique()):
            mask = adata_combined.obs['stage'] == stage
            ax.scatter(
                adata_combined.obsm['X_umap'][mask, 0],
                adata_combined.obsm['X_umap'][mask, 1],
                s=1, alpha=0.3, label=stage,
                c=stage_colors.get(stage, 'blue')
            )
        ax.legend(markerscale=5)
        ax.set_title('UMAP colored by Stage')
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')

    plt.tight_layout()
    plt.savefig(save_dir / "umap_batch_vs_biology.png", dpi=150)
    plt.close()

    # 2. LISI distribution by lane
    if 'lisi_lane' in adata_combined.obs.columns:
        fig, ax = plt.subplots(figsize=(8, 6))
        adata_combined.obs.boxplot(column='lisi_lane', by='lane', ax=ax)
        ax.axhline(BATCH_THRESHOLDS['min_lisi'], color='red', linestyle='--',
                   label=f"Good mixing threshold: {BATCH_THRESHOLDS['min_lisi']}")
        ax.set_title('LISI Score by Lane (Higher = Better Mixing)')
        ax.set_xlabel('Lane')
        ax.set_ylabel('LISI')
        ax.legend()
        plt.suptitle('')
        plt.tight_layout()
        plt.savefig(save_dir / "lisi_by_lane.png", dpi=150)
        plt.close()

    # 3. PCA variance by lane
    if 'X_pca' in adata_combined.obsm:
        fig, ax = plt.subplots(figsize=(8, 6))

        # Compute variance explained by lane for each PC
        pc_data = pd.DataFrame(
            adata_combined.obsm['X_pca'][:, :10],
            columns=[f'PC{i+1}' for i in range(10)]
        )
        pc_data['lane'] = adata_combined.obs['lane'].values

        # ANOVA F-statistic for each PC
        f_stats = []
        for pc in [f'PC{i+1}' for i in range(10)]:
            groups = [pc_data[pc_data['lane'] == l][pc].values
                      for l in pc_data['lane'].unique()]
            f, p = stats.f_oneway(*groups)
            f_stats.append({'PC': pc, 'F_statistic': f, 'p_value': p})

        f_df = pd.DataFrame(f_stats)
        ax.bar(f_df['PC'], f_df['F_statistic'])
        ax.axhline(10, color='red', linestyle='--', label='Concern threshold')
        ax.set_xlabel('Principal Component')
        ax.set_ylabel('F-statistic (Lane)')
        ax.set_title('Lane Effect on Principal Components')
        ax.legend()

        plt.tight_layout()
        plt.savefig(save_dir / "pca_lane_variance.png", dpi=150)
        plt.close()


def generate_qc_report(all_metrics: pd.DataFrame, batch_metrics: dict, save_path: Path):
    """Generate comprehensive QC report in Markdown format."""

    n_pass = (all_metrics['qc_status'] == 'PASS').sum()
    n_warn = (all_metrics['qc_status'] == 'WARN').sum()
    n_fail = (all_metrics['qc_status'] == 'FAIL').sum()

    report = f"""# G4X Full Dataset QC Report

**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}

## Summary

| Metric | Value |
|--------|-------|
| Total samples | {len(all_metrics)} |
| PASS | {n_pass} |
| WARN | {n_warn} |
| FAIL | {n_fail} |
| Total cells (all) | {all_metrics['n_cells'].sum():,} |
| Total cells (passing) | {all_metrics[all_metrics['qc_status'] != 'FAIL']['n_cells'].sum():,} |

## QC Thresholds

| Threshold | Value |
|-----------|-------|
| min_cells | {QC_THRESHOLDS['min_cells']:,} |
| min_median_transcripts | {QC_THRESHOLDS['min_median_transcripts_per_cell']} |
| min_median_genes | {QC_THRESHOLDS['min_median_genes_per_cell']} |
| max_pct_empty | {QC_THRESHOLDS['max_pct_empty']}% |

## Failed Samples

"""

    failed = all_metrics[all_metrics['qc_status'] == 'FAIL']
    if len(failed) > 0:
        for _, row in failed.iterrows():
            report += f"### {row['sample_id']}\n"
            report += f"- **Lane:** {row['lane']}\n"
            report += f"- **Stage:** {row['stage']}\n"
            report += f"- **Cells:** {row['n_cells']:,}\n"
            report += f"- **Median counts:** {row['median_counts']:.1f}\n"
            report += f"- **Median genes:** {row['median_genes']:.1f}\n"
            report += f"- **Reasons:** {row['qc_reasons']}\n\n"
    else:
        report += "No samples failed QC.\n\n"

    report += """## Warned Samples

"""
    warned = all_metrics[all_metrics['qc_status'] == 'WARN']
    if len(warned) > 0:
        for _, row in warned.iterrows():
            report += f"- **{row['sample_id']}**: {row['qc_reasons']}\n"
    else:
        report += "No warnings.\n"

    report += f"""

## Batch Effect Assessment

| Metric | Value | Threshold | Status |
|--------|-------|-----------|--------|
| Mean LISI (lane) | {batch_metrics.get('mean_lisi', 0):.2f} | > {BATCH_THRESHOLDS['min_lisi']} | {'PASS' if batch_metrics.get('mean_lisi', 0) > BATCH_THRESHOLDS['min_lisi'] else 'CONCERN'} |
| Silhouette (lane) | {batch_metrics.get('silhouette_lane', 0):.3f} | < {BATCH_THRESHOLDS['max_silhouette_batch']} | {'PASS' if batch_metrics.get('silhouette_lane', 1) < BATCH_THRESHOLDS['max_silhouette_batch'] else 'CONCERN'} |

## Lane Statistics

"""

    for lane in sorted(all_metrics['lane'].unique()):
        lane_data = all_metrics[all_metrics['lane'] == lane]
        report += f"### {lane}\n"
        report += f"- Samples: {len(lane_data)}\n"
        report += f"- Total cells: {lane_data['n_cells'].sum():,}\n"
        report += f"- Mean median counts: {lane_data['median_counts'].mean():.1f}\n"
        report += f"- Pass/Warn/Fail: {(lane_data['qc_status']=='PASS').sum()}/{(lane_data['qc_status']=='WARN').sum()}/{(lane_data['qc_status']=='FAIL').sum()}\n\n"

    report += """## Recommendations

1. **Exclude** all FAIL samples from downstream analysis
2. **Monitor** WARN samples for unusual patterns
3. **Consider Harmony** batch correction if LISI < 2.0
4. **Verify** Lane 4 cell sizes don't confound results

## Files Generated

- `sample_qc_summary.csv`: All metrics and status
- `figures/per_sample/`: Individual QC panels
- `figures/cross_sample/`: Comparison plots
- `figures/batch_effects/`: Batch assessment
"""

    with open(save_path, 'w') as f:
        f.write(report)


# =============================================================================
# Main
# =============================================================================

def main():
    """Run comprehensive QC analysis."""
    logger.info("=" * 60)
    logger.info("G4X Comprehensive QC Analysis")
    logger.info("=" * 60)

    # Load manifest
    manifest_path = INPUT_DIR / "loading_manifest.csv"
    if not manifest_path.exists():
        logger.error(f"Manifest not found: {manifest_path}")
        logger.error("Run 60_load_all_samples.py first")
        return

    manifest = pd.read_csv(manifest_path)
    logger.info(f"Found {len(manifest)} samples in manifest")

    # Load Resolve baseline for validation
    resolve_baseline = None
    if RESOLVE_METRICS_PATH.exists():
        resolve_baseline = pd.read_csv(RESOLVE_METRICS_PATH)
        logger.info(f"Loaded Resolve baseline: {len(resolve_baseline)} samples")

    # Process each sample
    all_metrics = []
    combined_data = []

    for _, row in tqdm(manifest.iterrows(), total=len(manifest), desc="Processing samples"):
        sample_id = row['sample_id']
        h5ad_path = Path(row['path'])

        try:
            # Load sample
            adata = sc.read_h5ad(h5ad_path)

            # Compute metrics
            metrics = compute_sample_qc_metrics(adata)
            status, reasons = apply_qc_thresholds(metrics)
            metrics['qc_status'] = status
            metrics['qc_reasons'] = "; ".join(reasons) if reasons else "OK"

            # Compare to Resolve baseline
            if resolve_baseline is not None:
                resolve_row = resolve_baseline[resolve_baseline['sample_id'] == sample_id]
                if len(resolve_row) > 0:
                    metrics['resolve_median_trans'] = resolve_row['median_transcripts_per_cell'].values[0]
                    metrics['resolve_pct_empty'] = resolve_row['pct_empty_cells'].values[0]

            all_metrics.append(metrics)

            # Generate per-sample QC panel
            plot_sample_qc_panel(adata, sample_id, FIG_DIR_SAMPLE / f"{sample_id}_qc.png")

            # Store for combined analysis (subsample for memory)
            if adata.n_obs > 10000:
                sc.pp.subsample(adata, n_obs=10000)
            combined_data.append(adata)

            del adata
            gc.collect()

        except Exception as e:
            logger.error(f"Error processing {sample_id}: {e}")
            all_metrics.append({
                'sample_id': sample_id,
                'qc_status': 'ERROR',
                'qc_reasons': str(e)
            })

    # Create metrics DataFrame
    metrics_df = pd.DataFrame(all_metrics)
    metrics_df.to_csv(OUTPUT_DIR / "sample_qc_summary.csv", index=False)

    # Cross-sample plots
    logger.info("Generating cross-sample plots...")
    plot_cross_sample_comparison(metrics_df, FIG_DIR_CROSS)

    # Batch effect analysis
    logger.info("Performing batch effect analysis...")
    batch_metrics = {}

    if combined_data:
        # Concatenate (subsampled) data
        adata_combined = ad.concat(combined_data, join='outer')
        adata_combined.obs_names_make_unique()

        # Normalize and compute embeddings
        sc.pp.normalize_total(adata_combined, target_sum=1e4)
        sc.pp.log1p(adata_combined)
        sc.pp.highly_variable_genes(adata_combined, n_top_genes=200, subset=True)
        sc.pp.pca(adata_combined, n_comps=30)
        sc.pp.neighbors(adata_combined)
        sc.tl.umap(adata_combined)

        # Compute LISI
        lane_labels = adata_combined.obs['lane'].values
        lisi_scores = compute_lisi(adata_combined.obsm['X_pca'], lane_labels)
        adata_combined.obs['lisi_lane'] = lisi_scores
        batch_metrics['mean_lisi'] = np.mean(lisi_scores)

        # Compute silhouette
        try:
            sil = silhouette_score(
                adata_combined.obsm['X_pca'],
                adata_combined.obs['lane'],
                sample_size=min(10000, adata_combined.n_obs)
            )
            batch_metrics['silhouette_lane'] = sil
        except Exception as e:
            logger.warning(f"Could not compute silhouette: {e}")
            batch_metrics['silhouette_lane'] = np.nan

        # Batch effect plots
        plot_batch_effects(adata_combined, FIG_DIR_BATCH)

        del adata_combined
        gc.collect()

    # Generate report
    logger.info("Generating QC report...")
    generate_qc_report(metrics_df, batch_metrics, OUTPUT_DIR / "QC_REPORT.md")

    # Summary
    logger.info("=" * 60)
    logger.info("QC Analysis Complete")
    logger.info(f"  PASS: {(metrics_df['qc_status'] == 'PASS').sum()}")
    logger.info(f"  WARN: {(metrics_df['qc_status'] == 'WARN').sum()}")
    logger.info(f"  FAIL: {(metrics_df['qc_status'] == 'FAIL').sum()}")
    logger.info(f"  Output: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
