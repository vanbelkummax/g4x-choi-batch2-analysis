# G4X Choi Batch 2 Analysis Guide

## Overview

This analysis examines gastric cancer progression using 3 samples from the G4X platform:
- **E02**: Normal gastric tissue (27,775 cells)
- **F02**: Metaplasia (30,229 cells)
- **G02**: Cancer (76,463 cells)

**Total: 134,467 cells across 337 targeted genes**

---

## Directory Structure

```
output/
├── ANALYSIS_GUIDE.md          # This file
├── figures/                   # All visualizations
├── data/                      # CSV results and recommendations
└── scripts/                   # Python scripts that generated everything
```

---

## Figure-to-Script Mapping

### 1. Data Loading & QC
**Script:** `scripts/00_load_pilot_samples.py`

| Output | Description |
|--------|-------------|
| (no figures) | Loads raw G4X data, applies minimal QC (min_genes=10, min_cells=3), saves h5ad files |

**What it does:**
- Reads cell-by-transcript matrices from G4X output
- Removes low-quality cells (few genes detected)
- Removes genes detected in very few cells
- Normalizes and log-transforms expression
- Saves individual sample files + merged dataset

---

### 2. HVG Optimization
**Script:** `scripts/01_hvg_optimize_parallel.py`

| Figure | Description |
|--------|-------------|
| `hvg_silhouette_heatmap.png` | Heatmap showing silhouette scores across HVG counts (30, 50, 100, 200, all) and clustering resolutions (0.3, 0.5, 0.8, 1.0) |

| Data File | Description |
|-----------|-------------|
| `hvg_optimization.csv` | All 60 configurations with silhouette and Calinski-Harabasz scores |
| `hvg_recommendation.txt` | Best settings per sample |

**What it does:**
- Tests 5 HVG settings × 4 resolutions = 20 configs per sample
- Uses GPU-accelerated PCA (PyTorch)
- Computes silhouette score (cluster separation quality)
- Runs 20 workers in parallel for speed

**How to interpret:**
- Higher silhouette = better cluster separation
- Green cells in heatmap = good configurations
- **Key finding:** Cancer (G02) shows 2x better separation (0.193) than normal/metaplasia (~0.11)
- This suggests cancer develops more distinct cellular subpopulations

---

### 3. Cell Type Annotation
**Script:** `scripts/02_annotate_multimethod.py`

| Figure | Description |
|--------|-------------|
| `annotation_method_ranking.png` | Comparison of annotation methods |
| `celltype_distribution.png` | Cell type counts per sample |
| `celltype_proportions.png` | Stacked bar chart of cell type proportions |

| Data File | Description |
|-----------|-------------|
| `annotation_method_comparison.csv` | Metrics for each annotation method |
| `annotation_method_rankings.csv` | Ranked methods by performance |
| `annotation_recommendation.txt` | Best method recommendation |

**What it does:**
- Implements 3 annotation methods:
  1. **Marker scoring** - Uses curated gastric markers (best performer)
  2. **CellTypist** - Pre-trained ML model
  3. **Hierarchical gating** - Threshold-based classification
- Compares methods by cluster purity and unknown rate

**How to interpret:**
- Marker-based annotation works best for this targeted panel
- Look for shifts in cell type proportions across disease stages
- High proportion of a cell type in cancer but not normal = disease-associated

---

### 4. Spatial Visualization
**Script:** `scripts/03_spatial_viz.py`

| Figure | Description |
|--------|-------------|
| `E02_umap.png` | UMAP embedding colored by cell type (Normal) |
| `F02_umap.png` | UMAP embedding colored by cell type (Metaplasia) |
| `G02_umap.png` | UMAP embedding colored by cell type (Cancer) |
| `E02_spatial.png` | Spatial plot with cell positions (Normal) |
| `F02_spatial.png` | Spatial plot with cell positions (Metaplasia) |
| `G02_spatial.png` | Spatial plot with cell positions (Cancer) |
| `progression_3panel.png` | Side-by-side comparison: Normal → Metaplasia → Cancer |

**What it does:**
- Creates UMAP dimensionality reduction for visualization
- Plots cells in their original spatial coordinates
- Colors by cell type annotation

**How to interpret:**
- **UMAP:** Cells close together = similar expression profiles
- **Spatial plots:** Shows tissue architecture and cell organization
- **Progression panel:** Look for changes in cell distribution across stages
- Clusters that appear in cancer but not normal = disease-specific populations

---

### 5. Differential Expression Analysis
**Script:** `scripts/04_exploratory_deg.py`

| Figure | Description |
|--------|-------------|
| `deg_volcano_global_Cancer_vs_Normal.png` | Volcano plot: Cancer vs Normal |
| `deg_volcano_global_Metaplasia_vs_Normal.png` | Volcano plot: Metaplasia vs Normal |
| `deg_volcano_global_Cancer_vs_Metaplasia.png` | Volcano plot: Cancer vs Metaplasia |
| `deg_heatmap_top.png` | Heatmap of top differentially expressed genes |

| Data File | Description |
|-----------|-------------|
| `deg_exploratory_full.csv` | All DEG results (global + cell-type specific) |
| `deg_exploratory_top.csv` | Top DEGs filtered by fold-change |

**What it does:**
- Compares gene expression between disease stages
- Uses Wilcoxon rank-sum test
- Calculates log2 fold-change
- Identifies both global and cell-type-specific changes

**How to interpret:**
- **Volcano plots:**
  - X-axis = log2 fold-change (positive = upregulated in condition 2)
  - Y-axis = -log10(p-value) (higher = more significant)
  - Top-right quadrant = significantly upregulated
  - Top-left quadrant = significantly downregulated
- **Key DEGs found:**
  - COL1A1 ↑3.5x in cancer = CAF/stromal infiltration
  - MUC2 ↑3.0x in cancer = goblet cell differentiation
  - CLDN3 ↑2.9x in cancer = tight junction remodeling
  - PGC, ATP4A ↓ in cancer = loss of gastric identity

**CAUTION:** N=1 patient - p-values are exploratory only. Rankings by fold-change are meaningful.

---

### 6. Spatial Statistics
**Script:** `scripts/05_spatial_statistics.py`

| Figure | Description |
|--------|-------------|
| `spatial_neighborhood_enrichment.png` | Heatmap showing which cell types co-localize |
| `spatial_co_occurrence.png` | Co-occurrence probability curves |

| Data File | Description |
|-----------|-------------|
| `spatial_colocalization_summary.csv` | Z-scores for all cell type pairs |

**What it does:**
- Builds spatial neighbor graph (15 nearest neighbors)
- Computes neighborhood enrichment (which cell types are neighbors more than expected by chance)
- Calculates co-occurrence statistics

**How to interpret:**
- **Neighborhood enrichment heatmap:**
  - Red = cell types found together more than expected (co-localization)
  - Blue = cell types found together less than expected (avoidance)
  - Diagonal = self-association
- **Z-scores:**
  - |Z| > 2 = significant enrichment/depletion
  - Higher Z = stronger association
- **Key findings:**
  - Cancer shows Z=145 for clusters 2↔4 = extremely strong co-localization
  - This suggests specific cell type niches form in tumor microenvironment

---

## Key Biological Findings

### 1. Cancer Has More Distinct Subpopulations
- Silhouette score: Cancer (0.193) >> Normal/Metaplasia (~0.11)
- Interpretation: Tumor evolution creates distinct cellular states

### 2. Stromal Infiltration in Cancer
- COL1A1 (collagen) upregulated 3.5x
- Cancer-associated fibroblasts (CAFs) expand in tumor

### 3. Epithelial Remodeling
- MUC2 (mucin) upregulated = goblet cell differentiation
- CLDN3 (claudin) upregulated = tight junction changes
- PGC, ATP4A downregulated = loss of normal gastric cell identity

### 4. Spatial Reorganization
- Normal tissue: moderate cell type mixing
- Cancer: strong cell type clustering (distinct niches)

---

## Reproducing the Analysis

```bash
cd /home/user/g4x-choi-batch2-analysis
conda activate enact

# Run full pipeline
python scripts/00_load_pilot_samples.py      # ~2 min
python scripts/01_hvg_optimize_parallel.py   # ~6 min (parallel)
python scripts/02_annotate_multimethod.py    # ~3 min
python scripts/03_spatial_viz.py             # ~2 min
python scripts/04_exploratory_deg.py         # ~3 min
python scripts/05_spatial_statistics.py      # ~2 min
```

---

## Limitations

1. **N=1 patient** - All statistics are exploratory; no biological replicates
2. **337-gene panel** - Limited to targeted genes, not whole transcriptome
3. **Annotation uncertainty** - Some cells may be misclassified
4. **No protein validation** - G4X has IF data we haven't integrated yet

---

## Next Steps

1. **Integrate protein data** - Validate top DEGs at protein level
2. **Trajectory analysis** - Model progression Normal → Metaplasia → Cancer
3. **Expand to full cohort** - Analyze all available samples
4. **Literature validation** - Cross-reference findings with published gastric cancer studies

---

## Contact

Analysis performed: 2026-01-22
Platform: G4X (targeted spatial transcriptomics)
Environment: conda activate enact (Python 3.9 + scanpy + squidpy)
