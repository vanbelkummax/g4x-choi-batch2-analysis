# G4X Gastric Cancer: Detailed QC Checklist

**Date:** 2026-01-21 | **Dataset:** Choi Batch 2 (8 samples, 522K cells)

---

## Phase 1: Per-Sample QC (MUST DO FOR EACH SAMPLE)

### 1.1 Cell-Level Metrics
- [ ] **A01** - Violin plots: total_counts, n_genes_by_counts, pct_counts_mt
- [ ] **B01** - Violin plots: total_counts, n_genes_by_counts, pct_counts_mt
- [ ] **C01** - Violin plots: total_counts, n_genes_by_counts, pct_counts_mt
- [ ] **D01** - Violin plots: total_counts, n_genes_by_counts, pct_counts_mt
- [ ] **E01** - Violin plots: total_counts, n_genes_by_counts, pct_counts_mt
- [ ] **F01** - Violin plots: total_counts, n_genes_by_counts, pct_counts_mt
- [ ] **G01** - Violin plots: total_counts, n_genes_by_counts, pct_counts_mt
- [ ] **H01** - Violin plots: total_counts, n_genes_by_counts, pct_counts_mt

### 1.2 Spatial Distribution QC
- [ ] **A01** - Spatial scatter: total_counts heatmap
- [ ] **B01** - Spatial scatter: total_counts heatmap
- [ ] **C01** - Spatial scatter: total_counts heatmap
- [ ] **D01** - Spatial scatter: total_counts heatmap
- [ ] **E01** - Spatial scatter: total_counts heatmap
- [ ] **F01** - Spatial scatter: total_counts heatmap
- [ ] **G01** - Spatial scatter: total_counts heatmap
- [ ] **H01** - Spatial scatter: total_counts heatmap

### 1.3 Gene Detection QC
- [ ] **A01** - Histogram: genes detected per cell
- [ ] **B01** - Histogram: genes detected per cell
- [ ] **C01** - Histogram: genes detected per cell
- [ ] **D01** - Histogram: genes detected per cell
- [ ] **E01** - Histogram: genes detected per cell
- [ ] **F01** - Histogram: genes detected per cell
- [ ] **G01** - Histogram: genes detected per cell
- [ ] **H01** - Histogram: genes detected per cell

---

## Phase 2: Cross-Sample Comparisons

### 2.1 Distribution Comparisons
- [ ] Ridge plot: total_counts across all 8 samples
- [ ] Ridge plot: n_genes across all 8 samples
- [ ] Box plot: cells per sample with outlier detection
- [ ] Stacked bar: cell type proportions per sample

### 2.2 Batch Effect Assessment
- [ ] UMAP colored by sample_id (check mixing)
- [ ] UMAP colored by stage (N/M/C separation)
- [ ] PCA variance by sample vs biological factors
- [ ] Silhouette score: sample vs cell_type clustering

---

## Phase 3: RNA QC Metrics

### 3.1 Count Distribution
- [ ] Scatter: total_counts vs n_genes (per sample, colored)
- [ ] Histogram: log10(total_counts) distribution
- [ ] Density plot: counts per gene across samples

### 3.2 Gene-Level QC
- [ ] Top 20 highest expressed genes (bar plot)
- [ ] Cumulative counts: top N genes explain X% variance
- [ ] Zero-inflation: % cells expressing each gene (histogram)
- [ ] Spatial autocorrelation (Moran's I) for top variable genes

### 3.3 Filtering Thresholds
- [ ] Document: min_counts threshold applied
- [ ] Document: min_genes threshold applied
- [ ] Document: max_pct_mt threshold applied
- [ ] Before/after filtering cell counts table

---

## Phase 4: Protein QC Metrics

### 4.1 Protein Signal Quality
- [ ] Violin plots: all 17 protein markers per sample
- [ ] Isotype control assessment (background levels)
- [ ] Heatmap: protein-protein correlations

### 4.2 RNA-Protein Concordance
- [ ] Scatter: CD8A (RNA) vs CD8 (protein) with correlation
- [ ] Scatter: PDCD1 (RNA) vs PD1 (protein) with correlation
- [ ] Scatter: MS4A1 (RNA) vs CD20 (protein) with correlation
- [ ] Scatter: MKI67 (RNA) vs KI67 (protein) with correlation
- [ ] Summary table: all RNA-protein pair correlations

### 4.3 Multi-Modal Integration QC
- [ ] UMAP: RNA-only embedding
- [ ] UMAP: Protein-only embedding
- [ ] UMAP: WNN combined embedding
- [ ] Silhouette comparison across modalities

---

## Phase 5: Spatial QC

### 5.1 Tissue Coverage
- [ ] Spatial plot: cell density heatmap per sample
- [ ] Spatial plot: identify tissue holes/artifacts
- [ ] Convex hull area per sample (coverage metric)

### 5.2 Cell Type Spatial Distribution
- [ ] Spatial plot: cell types colored (each sample)
- [ ] Neighborhood enrichment heatmap per sample
- [ ] Ripley's L function: cell type clustering assessment

### 5.3 Segmentation Quality
- [ ] Histogram: cell area distribution
- [ ] Spatial plot: flag unusually large/small cells
- [ ] Doublet assessment via cell size + transcript count

---

## Phase 6: Annotation QC

### 6.1 Cell Type Validation
- [ ] Dot plot: canonical markers by cell type
- [ ] Violin plot: marker expression per cell type
- [ ] Confusion assessment: ambiguous cell assignments

### 6.2 Stage Annotation Validation
- [ ] Composition bar chart: cell types per stage (N/M/C)
- [ ] Chi-square test: cell type ~ stage association
- [ ] Marker expression: CDX2, MUC2 by stage (intestinalization)

---

## Phase 7: Summary Statistics

### 7.1 Final Dataset Metrics
- [ ] Table: cells per sample (before/after QC)
- [ ] Table: median genes per cell per sample
- [ ] Table: median counts per cell per sample
- [ ] Table: cell type counts per sample

### 7.2 QC Report Generation
- [ ] Compile all figures into HTML report
- [ ] Flag any samples/cells requiring attention
- [ ] Document all filtering decisions

---

## QC Thresholds Applied

| Metric | Threshold | Rationale |
|--------|-----------|-----------|
| min_counts | TBD | Minimum UMI per cell |
| min_genes | TBD | Minimum genes detected |
| max_pct_mt | TBD | Maximum mitochondrial % |
| min_cells | TBD | Minimum cells per gene |

---

## Output Files

All QC figures saved to: `results/qc_figures/`
- `{sample}_violin_qc.png` - Per-sample violin plots
- `{sample}_spatial_qc.png` - Per-sample spatial QC
- `cross_sample_comparison.png` - Batch comparison
- `rna_protein_concordance.png` - Multi-modal QC
- `qc_summary_report.html` - Compiled report

---

## Script

Run: `python scripts/50_comprehensive_qc.py`
