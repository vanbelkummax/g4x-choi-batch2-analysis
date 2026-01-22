# G4X Pilot Analysis - Session Summary

**Date:** 2026-01-22
**Repository:** `/home/user/g4x-choi-batch2-analysis/`
**Branch:** `pilot-clean`
**Environment:** `conda activate enact`

---

## Completed Work

### 1. Data Loading (00_load_pilot_samples.py) ✅
- Loaded 3 gastric cancer progression samples from G4X platform
- Removed scrublet (unreliable for targeted 337-gene panels)
- Applied minimal QC (min_genes=10, min_cells=3)

| Sample | Stage | Cells After QC | Genes |
|--------|-------|----------------|-------|
| E02 | Normal | 27,775 | 338 |
| F02 | Metaplasia | 30,229 | 337 |
| G02 | Cancer | 76,463 | 338 |
| **TOTAL** | - | **134,467** | 337 |

### 2. HVG Optimization (01_hvg_optimize.py) 🔄
- Script created with grid search: HVG=[30,50,100,200,'all'] x Resolution=[0.3,0.5,0.8,1.0]
- Uses silhouette score + Calinski-Harabasz for evaluation
- **Status:** Script works but takes time with 134K cells - needs to be run separately

### 3. Cell Type Annotation (02_annotate_multimethod.py) ✅
Three methods implemented:
1. **Marker Scoring** - Using curated gastric-specific markers (WORKING)
2. **CellTypist** - Fixed API but needs testing (was failing)
3. **Hierarchical Gating** - Simple threshold-based (WORKING)

**Marker-based results (cluster purity ~0.4-0.6):**
- Identifies major populations: Epithelial, T cells, B cells, Macrophages, Fibroblasts
- Low unknown rate (<1%)

### 4. Visualization (03_spatial_viz.py) ✅
Generated 16 figures in `results/pilot/figures/`:
- UMAP plots per sample (E02, F02, G02)
- Spatial plots per sample
- 3-panel progression figure
- Cell type proportions
- Annotation method comparison

### 5. DEG Analysis (04_exploratory_deg.py) ✅
**Key findings (exploratory, n=1 patient):**

| Comparison | Up | Down |
|------------|-----|------|
| Normal→Metaplasia | GAST, TFF1, TFF3 | CBLIF, ATP4A, PGC |
| Normal→Cancer | COL1A1, MUC2, CLDN3 | PGC, TFF2, IGHA1 |
| Metaplasia→Cancer | COL1A1, CLDN3, MUC2 | TFF2, TFF1, GAST |

**Biological interpretation:**
- Cancer has increased CAF infiltration (COL1A1)
- Loss of gastric identity (PGC, ATP4A, CBLIF)
- Increased mucin (MUC2) and tight junctions (CLDN3)

---

## Output Files

```
results/pilot/
├── E02_raw.h5ad          # QC'd, normalized
├── F02_raw.h5ad
├── G02_raw.h5ad
├── merged_pilot.h5ad     # Combined, ready for analysis
├── E02_annotated.h5ad    # With cell type annotations
├── F02_annotated.h5ad
├── G02_annotated.h5ad
├── deg_exploratory_full.csv   # All DEGs
├── deg_exploratory_top.csv    # Top DEGs
├── annotation_method_comparison.csv
├── annotation_method_rankings.csv
└── figures/
    ├── *_umap.png        # 3 UMAP plots
    ├── *_spatial.png     # 3 spatial plots
    ├── progression_3panel.png
    ├── celltype_proportions.png
    ├── deg_volcano_*.png # 3 volcano plots
    └── deg_heatmap_top.png
```

---

## Next Session Tasks

### 1. Run HVG Optimization (PRIORITY)
```bash
cd /home/user/g4x-choi-batch2-analysis
conda activate enact
python scripts/01_hvg_optimize.py 2>&1 | tee hvg_output.log
```
- Takes ~10-15 min with 134K cells
- Will produce `hvg_optimization.csv` and `hvg_recommendation.txt`

### 2. Fix & Re-run CellTypist
The API was updated. Current fix is in place but untested:
```bash
python scripts/02_annotate_multimethod.py
```

### 3. Consider Additional Analysis
- **Spatial statistics** (neighborhood enrichment, Ripley's K)
- **Trajectory analysis** (if Normal→Metaplasia→Cancer progression)
- **scRNA reference mapping** (Kumar et al. 2022 gastric atlas - GSE183904)

### 4. Biological Validation
Top DEG genes (COL1A1, MUC2, GAST, TFF1/2) should be:
- Validated in protein data (G4X has IF channel)
- Cross-referenced with literature

---

## Quick Commands

```bash
# Activate environment
cd /home/user/g4x-choi-batch2-analysis
conda activate enact

# Check outputs
ls results/pilot/
ls results/pilot/figures/

# View DEG results
head -20 results/pilot/deg_exploratory_top.csv

# Load data in Python
import scanpy as sc
adata = sc.read_h5ad('results/pilot/merged_pilot.h5ad')
print(adata)  # 134,467 x 337
```

---

## Important Notes

1. **Statistical caveat:** N=1 patient - all p-values are meaningless. Rankings by fold-change only.
2. **Gene panel:** 337 targeted genes - not whole transcriptome. Limits some analyses.
3. **Cell types:** Marker-based annotation works well; CellTypist needs re-testing.
4. **GPU available:** RTX 5090 24GB - underutilized in current analysis.
