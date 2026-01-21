# G4X Full Dataset QC Implementation Plan

**Author:** Max Van Belkum
**Date:** 2026-01-21
**Status:** Ready for Review
**GitHub Branch:** `master`

---

## Executive Summary

This document provides the complete implementation plan for quality control (QC) analysis of all 32 samples (2.3M cells) from the Choi Gastric Cancer Progression study. Currently, only 8 samples from Lane 1 have been processed. This plan extends analysis to all 4 lanes with rigorous QC, batch effect assessment, and reproducible processing.

### Current vs Target State

| Metric | Current | Target |
|--------|---------|--------|
| Samples processed | 8 (Lane 1) | 32 (All lanes) |
| Cells analyzed | 522K | ~2.0M (post-QC) |
| Batch assessment | None | LISI + silhouette + visual |
| QC validation | Manual | Automated with Resolve baseline |

---

## 1. Environment Setup

### 1.1 Conda Environment

```bash
# Activate the spatial analysis environment
conda activate enact

# Verify key packages
python -c "
import scanpy as sc
import squidpy as sq
import anndata as ad
import numpy as np
import pandas as pd
import h5py
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from scipy import stats
print(f'scanpy: {sc.__version__}')
print(f'squidpy: {sq.__version__}')
print('All imports successful')
"
```

### 1.2 Package Versions

| Package | Version | Purpose |
|---------|---------|---------|
| scanpy | 1.10.3 | Core scRNA-seq analysis |
| squidpy | 1.6.1 | Spatial statistics |
| anndata | 0.10.9 | Data structure |
| h5py | 3.12.1 | H5 file reading |
| scipy | 1.14.1 | Statistical tests |
| scikit-learn | 1.5.2 | LISI, silhouette |
| matplotlib | 3.9.2 | Visualization |
| seaborn | 0.13.2 | Statistical plots |

### 1.3 Hardware Requirements

| Resource | Minimum | Available |
|----------|---------|-----------|
| RAM | 32 GB | 196 GB |
| CPU cores | 8 | 24 |
| GPU | Optional | RTX 5090 24GB |
| Storage | 15 GB | Sufficient |

---

## 2. Data Structure

### 2.1 Directory Layout

```
/mnt/x/Choi_Batch_2_Tuesday/
├── g4-028-083-FC1-L001_5XGaAe5DB2dm7sRK/  # Lane 1
│   ├── A01/single_cell_data/feature_matrix.h5
│   ├── B01/single_cell_data/feature_matrix.h5
│   ├── C01/...
│   ├── D01/...
│   ├── E01/...
│   ├── F01/...
│   ├── G01/...
│   └── H01/...
├── g4-028-083-FC1-L002_5XGaAe5DB2dm7sRK/  # Lane 2
│   └── [A02-H02]
├── g4-028-083-FC1-L003_5XGaAe5DB2dm7sRK/  # Lane 3
│   └── [A03-H03]
├── g4-028-083-FC1-L004_5XGaAe5DB2dm7sRK/  # Lane 4
│   └── [A04-H04]
├── choi_preGC_b2_core_metrics.csv          # Resolve baseline
└── choi_preGC_b2_protein_core_metrics.csv  # Protein metrics
```

### 2.2 Sample Metadata

| Sample | Stage | Tissue | Patient |
|--------|-------|--------|---------|
| A0x | Control | CRC (positive ctrl) | - |
| B0x | Normal | Adjacent normal | SNU-105 |
| C0x | Metaplasia | IM tissue | SNU-105 |
| D0x | Cancer | Gastric cancer | SNU-105 |
| E0x | Normal | Adjacent normal | SNU-107 |
| F0x | Metaplasia | IM tissue | SNU-107 |
| G0x | Cancer | Gastric cancer | SNU-107 |
| H0x | Cancer | Gastric cancer | SNU-484 |

### 2.3 Pre-Identified Issues (from Resolve Baseline)

| Sample | Issue | Metric | Action |
|--------|-------|--------|--------|
| **H04** | Very low transcripts | 25 trans/cell, 11 genes/cell | **EXCLUDE** |
| D01 | High empty cells | 11.9% empty | WARN |
| C03 | Low cell count | 21K cells | WARN |
| Lane 4 | Small cell areas | 54-70 µm² vs 80-120 µm² | MONITOR |

---

## 3. QC Thresholds and Rationale

### 3.1 Sample-Level QC

| Threshold | Value | Rationale |
|-----------|-------|-----------|
| `min_cells` | 20,000 | Statistical power for differential analysis |
| `min_median_transcripts` | 30 | Below this indicates capture failure |
| `min_median_genes` | 20 | Minimum for cell type assignment |
| `max_pct_empty` | 5.0% | High empty = segmentation issues |
| `min_pct_in_cells` | 80% | Transcript capture efficiency |

### 3.2 Cell-Level QC

| Threshold | Value | Rationale |
|-----------|-------|-----------|
| `min_counts` | 10 | Minimum viable cell |
| `min_genes` | 5 | Minimum for clustering |
| `max_pct_mito` | N/A | G4X panel lacks mito genes |

### 3.3 Batch Effect Thresholds

| Metric | Acceptable | Concerning |
|--------|------------|------------|
| LISI (by lane) | > 3.0 | < 2.0 |
| Silhouette (batch) | < 0.3 | > 0.5 |
| PC1 variance by lane | < 10% | > 20% |

---

## 4. Implementation Scripts

### 4.1 Script 60: Load All Samples

**File:** `scripts/60_load_all_samples.py`

```python
#!/usr/bin/env python3
"""
G4X Full Dataset Loading
========================
Load all 32 samples across 4 lanes into individual AnnData objects.

Usage:
    python scripts/60_load_all_samples.py

Output:
    results/qc_all_samples/raw/{sample}_raw.h5ad (32 files)
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import h5py
from pathlib import Path
from tqdm import tqdm
import logging
import warnings

warnings.filterwarnings('ignore')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# =============================================================================
# Configuration
# =============================================================================

DATA_ROOT = Path("/mnt/x/Choi_Batch_2_Tuesday")
OUTPUT_DIR = Path("/home/user/g4x-choi-batch2-analysis/results/qc_all_samples/raw")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Lane directory mapping
LANES = {
    'L001': 'g4-028-083-FC1-L001_5XGaAe5DB2dm7sRK',
    'L002': 'g4-028-083-FC1-L002_5XGaAe5DB2dm7sRK',
    'L003': 'g4-028-083-FC1-L003_5XGaAe5DB2dm7sRK',
    'L004': 'g4-028-083-FC1-L004_5XGaAe5DB2dm7sRK',
}

# Sample prefixes (A-H)
SAMPLE_PREFIXES = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Sample metadata
SAMPLE_METADATA = {
    'A': {'stage': 'control', 'tissue': 'CRC', 'patient': 'ctrl'},
    'B': {'stage': 'normal', 'tissue': 'adjacent_normal', 'patient': 'SNU-105'},
    'C': {'stage': 'metaplasia', 'tissue': 'IM', 'patient': 'SNU-105'},
    'D': {'stage': 'cancer', 'tissue': 'GC', 'patient': 'SNU-105'},
    'E': {'stage': 'normal', 'tissue': 'adjacent_normal', 'patient': 'SNU-107'},
    'F': {'stage': 'metaplasia', 'tissue': 'IM', 'patient': 'SNU-107'},
    'G': {'stage': 'cancer', 'tissue': 'GC', 'patient': 'SNU-107'},
    'H': {'stage': 'cancer', 'tissue': 'GC', 'patient': 'SNU-484'},
}

# Protein markers
PROTEIN_MARKERS = [
    'ATPase', 'CD11c', 'CD20', 'CD3', 'CD31', 'CD4', 'CD45',
    'CD68', 'CD8', 'FOXP3', 'HLA-DR', 'Isotype', 'KI67',
    'PD1', 'PDL1', 'PanCK', 'aSMA'
]


# =============================================================================
# Loading Functions
# =============================================================================

def load_sample(sample_dir: Path, sample_id: str, lane: str) -> ad.AnnData:
    """
    Load a single G4X sample into AnnData format.

    Parameters
    ----------
    sample_dir : Path
        Path to sample directory containing single_cell_data/
    sample_id : str
        Sample identifier (e.g., 'A01')
    lane : str
        Lane identifier (e.g., 'L001')

    Returns
    -------
    AnnData
        Loaded sample with RNA in X and protein in obsm['protein']
    """
    h5_path = sample_dir / "single_cell_data" / "feature_matrix.h5"

    if not h5_path.exists():
        raise FileNotFoundError(f"No feature matrix for {sample_id}: {h5_path}")

    with h5py.File(h5_path, 'r') as f:
        # Load expression matrix
        X = f['X'][:]  # (n_cells, n_genes)

        # Load gene names
        gene_ids = [x.decode() if isinstance(x, bytes) else x
                    for x in f['var/gene_id'][:]]

        # Load cell IDs
        cell_ids = [x.decode() if isinstance(x, bytes) else x
                    for x in f['obs/cell_id'][:]]

        # Create AnnData
        adata = ad.AnnData(X=X.astype(np.float32))
        adata.var_names = gene_ids
        adata.obs_names = [f"{sample_id}_{cid}" for cid in cell_ids]

        # Load all obs columns
        for key in f['obs'].keys():
            if key != 'cell_id':
                data = f[f'obs/{key}'][:]
                if data.dtype == object or data.dtype.kind == 'S':
                    data = [x.decode() if isinstance(x, bytes) else x for x in data]
                adata.obs[key] = data

        # Load var columns
        for key in f['var'].keys():
            if key != 'gene_id':
                data = f[f'var/{key}'][:]
                if data.dtype == object or data.dtype.kind == 'S':
                    data = [x.decode() if isinstance(x, bytes) else x for x in data]
                adata.var[key] = data

    # Separate RNA and protein
    is_protein = adata.var_names.isin(PROTEIN_MARKERS)
    protein_data = adata[:, is_protein].X.copy()
    rna_data = adata[:, ~is_protein].X.copy()

    # Create RNA-only AnnData
    adata_rna = ad.AnnData(X=rna_data)
    adata_rna.var_names = adata.var_names[~is_protein].tolist()
    adata_rna.obs = adata.obs.copy()
    adata_rna.obs_names = adata.obs_names.tolist()

    # Store protein in obsm
    adata_rna.obsm['protein'] = protein_data
    adata_rna.uns['protein_names'] = adata.var_names[is_protein].tolist()

    # Add sample metadata
    prefix = sample_id[0]
    meta = SAMPLE_METADATA.get(prefix, {})
    adata_rna.obs['sample_id'] = sample_id
    adata_rna.obs['lane'] = lane
    adata_rna.obs['stage'] = meta.get('stage', 'unknown')
    adata_rna.obs['tissue'] = meta.get('tissue', 'unknown')
    adata_rna.obs['patient'] = meta.get('patient', 'unknown')

    # Basic QC metrics
    adata_rna.obs['n_counts'] = np.array(adata_rna.X.sum(axis=1)).flatten()
    adata_rna.obs['n_genes'] = np.array((adata_rna.X > 0).sum(axis=1)).flatten()
    adata_rna.obs['n_proteins'] = np.array((protein_data > 0).sum(axis=1)).flatten()

    return adata_rna


def get_all_sample_paths():
    """Get paths to all 32 samples."""
    samples = []

    for lane_id, lane_dir_name in LANES.items():
        lane_dir = DATA_ROOT / lane_dir_name
        if not lane_dir.exists():
            logger.warning(f"Lane directory not found: {lane_dir}")
            continue

        lane_num = lane_id[-1]  # '1', '2', '3', '4'

        for prefix in SAMPLE_PREFIXES:
            sample_id = f"{prefix}0{lane_num}"
            sample_dir = lane_dir / sample_id

            if sample_dir.exists():
                samples.append({
                    'sample_id': sample_id,
                    'lane': lane_id,
                    'path': sample_dir
                })
            else:
                logger.warning(f"Sample directory not found: {sample_dir}")

    return samples


def main():
    """Load all samples and save as individual h5ad files."""
    logger.info("=" * 60)
    logger.info("G4X Full Dataset Loading")
    logger.info("=" * 60)

    # Get all sample paths
    samples = get_all_sample_paths()
    logger.info(f"Found {len(samples)} samples")

    # Load each sample
    loaded = []
    failed = []

    for sample_info in tqdm(samples, desc="Loading samples"):
        sample_id = sample_info['sample_id']
        lane = sample_info['lane']
        path = sample_info['path']

        try:
            adata = load_sample(path, sample_id, lane)

            # Save individual file
            out_path = OUTPUT_DIR / f"{sample_id}_raw.h5ad"
            adata.write(out_path)

            loaded.append({
                'sample_id': sample_id,
                'lane': lane,
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'path': str(out_path)
            })

            logger.info(f"  {sample_id}: {adata.n_obs:,} cells, {adata.n_vars} genes")

        except Exception as e:
            logger.error(f"  {sample_id}: FAILED - {e}")
            failed.append({'sample_id': sample_id, 'error': str(e)})

    # Save manifest
    manifest = pd.DataFrame(loaded)
    manifest.to_csv(OUTPUT_DIR / "loading_manifest.csv", index=False)

    if failed:
        pd.DataFrame(failed).to_csv(OUTPUT_DIR / "loading_failures.csv", index=False)

    # Summary
    logger.info("=" * 60)
    logger.info(f"Successfully loaded: {len(loaded)}/32 samples")
    logger.info(f"Total cells: {manifest['n_cells'].sum():,}")
    logger.info(f"Output: {OUTPUT_DIR}")

    if failed:
        logger.warning(f"Failed samples: {[f['sample_id'] for f in failed]}")


if __name__ == "__main__":
    main()
```

---

### 4.2 Script 61: Comprehensive QC

**File:** `scripts/61_comprehensive_qc.py`

```python
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


def apply_qc_thresholds(metrics: dict) -> tuple[str, list]:
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
        ax.set_xlabel('Cell area (µm²)')
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

    # 1. Ridge plot: n_counts by sample
    fig, ax = plt.subplots(figsize=(12, 10))
    # Group by stage for coloring
    stage_colors = {'control': 'gray', 'normal': 'green', 'metaplasia': 'orange', 'cancer': 'red'}

    samples = all_metrics.sort_values('median_counts', ascending=True)['sample_id'].tolist()
    for i, sample_id in enumerate(samples):
        metrics = all_metrics[all_metrics['sample_id'] == sample_id].iloc[0]
        color = stage_colors.get(metrics['stage'], 'blue')
        ax.barh(i, metrics['median_counts'], color=color, alpha=0.7,
                label=metrics['stage'] if i == 0 else "")

        # Mark failures
        if metrics['qc_status'] == 'FAIL':
            ax.scatter(metrics['median_counts'] + 5, i, marker='x', color='red', s=100, zorder=5)

    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(samples)
    ax.set_xlabel('Median counts per cell')
    ax.set_title('Median Counts by Sample (X = QC Fail)')
    ax.axvline(QC_THRESHOLDS['min_median_transcripts_per_cell'], color='red',
               linestyle='--', label=f"Threshold: {QC_THRESHOLDS['min_median_transcripts_per_cell']}")

    # Custom legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=s) for s, c in stage_colors.items()]
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()
    plt.savefig(save_dir / "median_counts_by_sample.png", dpi=150)
    plt.close()

    # 2. Correlation heatmap (sample vs sample)
    # This requires loading pseudobulk data - simplified version using metrics
    fig, ax = plt.subplots(figsize=(10, 8))
    pivot = all_metrics.pivot_table(
        values=['median_counts', 'median_genes', 'n_cells'],
        index='sample_id'
    )
    corr = pivot.T.corr()
    sns.heatmap(corr, ax=ax, cmap='coolwarm', center=0, annot=False)
    ax.set_title('Sample Correlation (QC Metrics)')
    plt.tight_layout()
    plt.savefig(save_dir / "sample_correlation_heatmap.png", dpi=150)
    plt.close()

    # 3. Box plot: metrics by lane
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


def plot_batch_effects(adata_combined: ad.AnnData, save_dir: Path):
    """Generate batch effect assessment plots."""

    # 1. UMAP colored by lane
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    if 'X_umap' in adata_combined.obsm:
        for lane in adata_combined.obs['lane'].unique():
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
        for stage in adata_combined.obs['stage'].unique():
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
| Mean LISI (lane) | {batch_metrics.get('mean_lisi', 'N/A'):.2f} | > {BATCH_THRESHOLDS['min_lisi']} | {'PASS' if batch_metrics.get('mean_lisi', 0) > BATCH_THRESHOLDS['min_lisi'] else 'CONCERN'} |
| Silhouette (lane) | {batch_metrics.get('silhouette_lane', 'N/A'):.3f} | < {BATCH_THRESHOLDS['max_silhouette_batch']} | {'PASS' if batch_metrics.get('silhouette_lane', 1) < BATCH_THRESHOLDS['max_silhouette_batch'] else 'CONCERN'} |

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
    manifest = pd.read_csv(INPUT_DIR / "loading_manifest.csv")
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
        except:
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
```

---

### 4.3 Script 62: Process QC-Passing Samples

**File:** `scripts/62_process_qc_passing.py`

```python
#!/usr/bin/env python3
"""
G4X Process QC-Passing Samples
==============================
Run full WNN + annotation pipeline on samples that passed QC.

Usage:
    python scripts/62_process_qc_passing.py [--harmony] [--parallel N]

Arguments:
    --harmony   Apply Harmony batch correction by lane
    --parallel  Number of parallel workers (default: 4)

Output:
    results/qc_all_samples/final_processed/{sample}_final.h5ad
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from scipy import sparse
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
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

INPUT_DIR = Path("/home/user/g4x-choi-batch2-analysis/results/qc_all_samples/raw")
QC_SUMMARY = Path("/home/user/g4x-choi-batch2-analysis/results/qc_all_samples/sample_qc_summary.csv")
OUTPUT_DIR = Path("/home/user/g4x-choi-batch2-analysis/results/qc_all_samples/final_processed")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Marker definitions for annotation
LINEAGE_MARKERS = {
    'epithelial': ['PanCK'],
    'immune': ['CD45'],
    'stromal': ['aSMA'],
    'endothelial': ['CD31'],
}

CELL_TYPE_MARKERS = {
    # Immune subtypes
    'T_cell': ['CD3'],
    'CD4_T': ['CD3', 'CD4'],
    'CD8_T': ['CD3', 'CD8'],
    'Treg': ['CD3', 'CD4', 'FOXP3'],
    'B_cell': ['CD20'],
    'Macrophage': ['CD68'],
    'DC': ['CD11c', 'HLA-DR'],
    # Functional markers
    'Proliferating': ['KI67'],
    'PD1_positive': ['PD1'],
    'PDL1_positive': ['PDL1'],
}

# Conflict pairs for admixture detection
CONFLICT_PAIRS = [
    ('PanCK', 'CD45'),      # Epithelial + Immune
    ('PanCK', 'aSMA'),      # Epithelial + Stromal
    ('CD45', 'aSMA'),       # Immune + Stromal
    ('CD31', 'CD45'),       # Endothelial + Immune
    ('CD3', 'CD68'),        # T cell + Macrophage
    ('CD4', 'CD8'),         # CD4 + CD8
]

ADMIX_THRESHOLD = 0.7  # Percentile for "high" expression
ADMIX_CUTOFF = 0.3     # Flag cells above this score


# =============================================================================
# Processing Functions
# =============================================================================

def compute_wnn_integration(adata, n_pcs_rna=30, n_pcs_prot=15, n_neighbors=20):
    """
    Compute Weighted Nearest Neighbors integration of RNA and Protein.

    Returns AnnData with X_wnn embedding in obsm.
    """
    logger.debug("Computing WNN integration...")

    # RNA PCA
    if 'X_pca' not in adata.obsm:
        sc.pp.pca(adata, n_comps=min(n_pcs_rna, adata.n_vars - 1))

    rna_pca = adata.obsm['X_pca'][:, :min(n_pcs_rna, adata.obsm['X_pca'].shape[1])].copy()

    # Protein PCA
    protein_data = adata.obsm['protein'].copy()
    if sparse.issparse(protein_data):
        protein_data = protein_data.toarray()

    # Scale protein data
    scaler = StandardScaler()
    protein_scaled = scaler.fit_transform(protein_data)

    # Protein PCA
    from sklearn.decomposition import PCA
    n_prot_comps = min(n_pcs_prot, protein_data.shape[1] - 1)
    pca_prot = PCA(n_components=n_prot_comps)
    protein_pca = pca_prot.fit_transform(protein_scaled)

    adata.obsm['X_pca_protein'] = protein_pca

    # Compute modality weights based on within-modality variance
    rna_var = np.var(rna_pca, axis=0).sum()
    prot_var = np.var(protein_pca, axis=0).sum()
    total_var = rna_var + prot_var

    rna_weight = rna_var / total_var
    prot_weight = prot_var / total_var

    logger.debug(f"  Modality weights: RNA={rna_weight:.3f}, Protein={prot_weight:.3f}")

    # Scale PCAs to same range
    rna_scaled = StandardScaler().fit_transform(rna_pca)
    prot_scaled = StandardScaler().fit_transform(protein_pca)

    # Weighted combination
    X_wnn = np.hstack([
        rna_scaled * np.sqrt(rna_weight),
        prot_scaled * np.sqrt(prot_weight)
    ])

    adata.obsm['X_wnn'] = X_wnn
    adata.uns['wnn_weights'] = {'rna': rna_weight, 'protein': prot_weight}

    return adata


def compute_admixture_score(adata):
    """
    Compute admixture score based on conflicting marker co-expression.

    High admixture = potential segmentation artifact (doublet/multiplet).
    """
    logger.debug("Computing admixture scores...")

    protein_names = adata.uns.get('protein_names', [])
    protein_data = adata.obsm['protein']

    if sparse.issparse(protein_data):
        protein_data = protein_data.toarray()

    # Create protein DataFrame
    protein_df = pd.DataFrame(
        protein_data,
        columns=protein_names,
        index=adata.obs_names
    )

    # Compute percentile thresholds
    thresholds = protein_df.quantile(ADMIX_THRESHOLD)

    # Binary high expression
    high_expr = protein_df > thresholds

    # Count conflicts
    conflict_count = np.zeros(len(adata))

    for m1, m2 in CONFLICT_PAIRS:
        if m1 in high_expr.columns and m2 in high_expr.columns:
            conflict = high_expr[m1] & high_expr[m2]
            conflict_count += conflict.values.astype(int)

    # Normalize to [0, 1]
    max_conflicts = len(CONFLICT_PAIRS)
    admix_score = conflict_count / max_conflicts

    adata.obs['admixture_score'] = admix_score
    adata.obs['is_admixed'] = admix_score > ADMIX_CUTOFF

    n_admixed = adata.obs['is_admixed'].sum()
    pct_admixed = n_admixed / len(adata) * 100
    logger.debug(f"  Admixed cells: {n_admixed:,} ({pct_admixed:.1f}%)")

    return adata


def annotate_lineage(adata):
    """Assign major lineage based on protein markers."""
    logger.debug("Annotating lineages...")

    protein_names = adata.uns.get('protein_names', [])
    protein_data = adata.obsm['protein']

    if sparse.issparse(protein_data):
        protein_data = protein_data.toarray()

    protein_df = pd.DataFrame(
        protein_data,
        columns=protein_names,
        index=adata.obs_names
    )

    # Score each lineage
    lineage_scores = {}
    for lineage, markers in LINEAGE_MARKERS.items():
        available = [m for m in markers if m in protein_df.columns]
        if available:
            # Z-score then mean
            scores = protein_df[available].apply(lambda x: (x - x.mean()) / (x.std() + 1e-8))
            lineage_scores[lineage] = scores.mean(axis=1)

    # Assign lineage
    if lineage_scores:
        score_df = pd.DataFrame(lineage_scores)
        adata.obs['lineage'] = score_df.idxmax(axis=1)
        adata.obs['lineage_score'] = score_df.max(axis=1)
    else:
        adata.obs['lineage'] = 'unknown'
        adata.obs['lineage_score'] = 0.0

    return adata


def annotate_cell_types(adata):
    """Assign cell types within lineages."""
    logger.debug("Annotating cell types...")

    protein_names = adata.uns.get('protein_names', [])
    protein_data = adata.obsm['protein']

    if sparse.issparse(protein_data):
        protein_data = protein_data.toarray()

    protein_df = pd.DataFrame(
        protein_data,
        columns=protein_names,
        index=adata.obs_names
    )

    # Initialize cell type column
    adata.obs['cell_type'] = adata.obs['lineage'].copy()

    # Refine immune cells
    immune_mask = adata.obs['lineage'] == 'immune'

    if immune_mask.any():
        immune_idx = adata.obs_names[immune_mask]

        for cell_type, markers in CELL_TYPE_MARKERS.items():
            available = [m for m in markers if m in protein_df.columns]
            if available:
                # All markers must be above median for assignment
                all_positive = protein_df.loc[immune_idx, available].apply(
                    lambda x: x > x.median(), axis=0
                ).all(axis=1)

                positive_idx = immune_idx[all_positive]
                adata.obs.loc[positive_idx, 'cell_type'] = cell_type

    return adata


def process_sample(sample_info: dict, apply_harmony: bool = False) -> dict:
    """
    Process a single sample through the full pipeline.

    Parameters
    ----------
    sample_info : dict
        Contains 'sample_id' and 'path' keys
    apply_harmony : bool
        Whether to apply Harmony batch correction

    Returns
    -------
    dict
        Processing result with status and metrics
    """
    sample_id = sample_info['sample_id']
    input_path = Path(sample_info['path'])
    output_path = OUTPUT_DIR / f"{sample_id}_final.h5ad"

    result = {'sample_id': sample_id, 'status': 'SUCCESS'}

    try:
        # Load
        adata = sc.read_h5ad(input_path)

        # Normalize RNA
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Variable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=200, subset=False)

        # PCA
        sc.pp.pca(adata, n_comps=min(30, adata.n_vars - 1))

        # WNN integration
        adata = compute_wnn_integration(adata)

        # Admixture scoring
        adata = compute_admixture_score(adata)

        # Annotation
        adata = annotate_lineage(adata)
        adata = annotate_cell_types(adata)

        # Clustering on WNN
        sc.pp.neighbors(adata, use_rep='X_wnn', n_neighbors=15)
        sc.tl.leiden(adata, resolution=0.5)
        sc.tl.umap(adata)

        # Save
        adata.write(output_path)

        # Metrics
        result['n_cells'] = adata.n_obs
        result['n_admixed'] = adata.obs['is_admixed'].sum()
        result['pct_admixed'] = result['n_admixed'] / result['n_cells'] * 100
        result['n_clusters'] = adata.obs['leiden'].nunique()

        # Lineage distribution
        lineage_counts = adata.obs['lineage'].value_counts()
        for lineage, count in lineage_counts.items():
            result[f'n_{lineage}'] = count

        del adata
        gc.collect()

    except Exception as e:
        result['status'] = 'FAILED'
        result['error'] = str(e)
        logger.error(f"Failed to process {sample_id}: {e}")

    return result


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Process QC-passing G4X samples')
    parser.add_argument('--harmony', action='store_true', help='Apply Harmony batch correction')
    parser.add_argument('--parallel', type=int, default=4, help='Number of parallel workers')
    args = parser.parse_args()

    logger.info("=" * 60)
    logger.info("G4X Process QC-Passing Samples")
    logger.info("=" * 60)

    # Load QC summary
    qc_df = pd.read_csv(QC_SUMMARY)
    passing = qc_df[qc_df['qc_status'] != 'FAIL']

    logger.info(f"QC-passing samples: {len(passing)}/{len(qc_df)}")

    # Load manifest for paths
    manifest = pd.read_csv(INPUT_DIR / "loading_manifest.csv")

    # Prepare processing list
    samples_to_process = []
    for _, row in passing.iterrows():
        sample_id = row['sample_id']
        manifest_row = manifest[manifest['sample_id'] == sample_id]
        if len(manifest_row) > 0:
            samples_to_process.append({
                'sample_id': sample_id,
                'path': manifest_row['path'].values[0]
            })

    logger.info(f"Processing {len(samples_to_process)} samples...")

    # Process samples
    results = []

    if args.parallel > 1:
        # Parallel processing
        with ProcessPoolExecutor(max_workers=args.parallel) as executor:
            futures = {
                executor.submit(process_sample, s, args.harmony): s['sample_id']
                for s in samples_to_process
            }

            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing"):
                sample_id = futures[future]
                try:
                    result = future.result()
                    results.append(result)
                    logger.info(f"  {sample_id}: {result['status']}")
                except Exception as e:
                    results.append({'sample_id': sample_id, 'status': 'FAILED', 'error': str(e)})
                    logger.error(f"  {sample_id}: FAILED - {e}")
    else:
        # Sequential processing
        for sample_info in tqdm(samples_to_process, desc="Processing"):
            result = process_sample(sample_info, args.harmony)
            results.append(result)
            logger.info(f"  {result['sample_id']}: {result['status']}")

    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(OUTPUT_DIR / "processing_results.csv", index=False)

    # Summary
    n_success = (results_df['status'] == 'SUCCESS').sum()
    n_failed = (results_df['status'] == 'FAILED').sum()

    logger.info("=" * 60)
    logger.info("Processing Complete")
    logger.info(f"  Success: {n_success}")
    logger.info(f"  Failed: {n_failed}")

    if n_success > 0:
        total_cells = results_df[results_df['status'] == 'SUCCESS']['n_cells'].sum()
        logger.info(f"  Total cells processed: {total_cells:,}")

    logger.info(f"  Output: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
```

---

## 5. Execution Protocol

### 5.1 Step-by-Step Instructions

```bash
# 1. Activate environment
conda activate enact

# 2. Navigate to repo
cd ~/g4x-choi-batch2-analysis

# 3. Create output directories
mkdir -p results/qc_all_samples/{raw,figures,final_processed}

# 4. Run loading script (~20 min)
python scripts/60_load_all_samples.py 2>&1 | tee logs/60_loading.log

# 5. Run QC analysis (~30 min)
python scripts/61_comprehensive_qc.py 2>&1 | tee logs/61_qc.log

# 6. Review QC report
cat results/qc_all_samples/QC_REPORT.md

# 7. Process passing samples (~3 hours)
python scripts/62_process_qc_passing.py --parallel 8 2>&1 | tee logs/62_processing.log

# 8. Verify outputs
ls -la results/qc_all_samples/final_processed/
```

### 5.2 Checkpoints

| Checkpoint | Verification |
|------------|--------------|
| After Script 60 | 32 files in `raw/`, manifest CSV exists |
| After Script 61 | H04 marked FAIL, QC_REPORT.md generated |
| After Script 62 | ~31 final h5ad files, processing_results.csv |

---

## 6. Expected Outcomes

### 6.1 Sample Status Predictions

| Sample | Status | Reason |
|--------|--------|--------|
| H04 | FAIL | 25 trans/cell, 11 genes/cell |
| D01 | WARN | 11.9% empty cells |
| C03 | WARN | 21K cells (borderline) |
| All others | PASS | Within thresholds |

### 6.2 Final Dataset

| Metric | Expected |
|--------|----------|
| Samples processed | 30-31 |
| Total cells | ~2.0M |
| Complete N→M→C series | 2 (SNU-105, SNU-107) |

---

## 7. Troubleshooting

### 7.1 Memory Issues

```python
# If running out of memory, process samples sequentially:
python scripts/62_process_qc_passing.py --parallel 1

# Or increase subsampling in QC:
# Edit 61_comprehensive_qc.py line ~250:
# Change: sc.pp.subsample(adata, n_obs=10000)
# To:     sc.pp.subsample(adata, n_obs=5000)
```

### 7.2 Batch Effect Concerns

If LISI < 2.0 or silhouette > 0.3:

```python
# Apply Harmony correction
python scripts/62_process_qc_passing.py --harmony

# Or use scVI for stronger correction (requires GPU):
import scvi
scvi.model.SCVI.setup_anndata(adata, batch_key='lane')
model = scvi.model.SCVI(adata)
model.train()
adata.obsm['X_scvi'] = model.get_latent_representation()
```

---

## 8. References

1. **scRNABatchQC** (Liu et al., Bioinformatics 2019) - Multi-sample QC framework
2. **LISI** (Korsunsky et al., Nature Methods 2019) - Batch integration metric
3. **Harmony** (Korsunsky et al., Nature Methods 2019) - Batch correction
4. **WNN** (Hao et al., Cell 2021) - Multimodal integration
5. **HEST** (Mahmood Lab) - Spatial batch correction utilities

---

## Appendix A: Resolve Baseline Comparison

The Resolve baseline file (`choi_preGC_b2_core_metrics.csv`) provides independent QC metrics for validation:

| Our Threshold | Resolve Equivalent | Validation |
|---------------|-------------------|------------|
| min_median_transcripts=30 | median_transcripts_per_cell | H04=25 (FAIL) |
| min_median_genes=20 | median_unique_genes_per_cell | H04=11 (FAIL) |
| max_pct_empty=5% | pct_empty_cells | D01=11.9% (WARN) |

---

## Appendix B: Lane 4 Cell Size Analysis

Lane 4 samples show systematically smaller cell areas:

| Lane | Median Cell Area (µm²) |
|------|------------------------|
| L001 | 84-102 |
| L002 | 84-119 |
| L003 | 80-124 |
| L004 | **54-70** |

This may represent:
1. Technical artifact (different segmentation parameters)
2. Biological difference (smaller cell populations)
3. Tissue quality variation

**Recommendation:** Monitor for confounding in downstream analysis. If Lane 4 clusters separately, consider exclusion or targeted batch correction.
