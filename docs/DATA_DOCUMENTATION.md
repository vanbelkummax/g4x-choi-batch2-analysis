# G4X Gastric Cancer Dataset - Comprehensive Documentation

## Executive Summary

| Metric | Value |
|--------|-------|
| **Study** | Choi Lab Gastric Pre-Cancer Batch 2 |
| **Platform** | G4X (Resolve Biosciences) |
| **Analysis Focus** | Lane 1: 8 samples, 522,188 cells |
| **Modalities** | RNA (387 genes) + Protein (17 markers) |
| **Processing** | Caretta v25.08.0 |

---

## Sample Inventory

### Lane 1 (Primary Analysis Set)

| Sample | Cells | Transcripts | Median Trans/Cell | Median Genes/Cell | Tissue Area (mm²) |
|--------|------:|------------:|------------------:|------------------:|------------------:|
| **A01** | 127,792 | 9,256,166 | 42 | 28 | 14.72 |
| **B01** | 30,194 | 2,799,376 | 69 | 36 | 5.01 |
| **C01** | 38,588 | 2,508,102 | 47 | 29 | 6.06 |
| **D01** | 31,711 | 2,627,133 | 64 | 35 | 4.61 |
| **E01** | 34,353 | 2,804,629 | 60 | 35 | 5.46 |
| **F01** | 62,875 | 2,940,220 | 34 | 26 | 7.94 |
| **G01** | 158,727 | 15,816,882 | 81 | 46 | 17.78 |
| **H01** | 85,842 | 9,199,621 | 84 | 42 | 9.84 |
| **TOTAL** | **570,082** | **47,952,129** | **60** (mean) | **35** (mean) | **71.42** |

### Quality Metrics

| Sample | % Empty Cells | % Transcripts in Cells | Q20 Score |
|--------|-------------:|-----------------------:|----------:|
| A01 | 0.77% | 94.12% | 91.16% |
| B01 | 2.22% | 89.50% | 89.77% |
| C01 | 1.74% | 89.71% | 89.25% |
| D01 | 11.91% | 88.17% | 89.74% |
| E01 | 0.63% | 89.06% | 89.03% |
| F01 | 1.52% | 93.75% | 90.79% |
| G01 | 0.08% | 95.22% | 88.13% |
| H01 | 0.25% | 93.94% | 89.41% |

---

## Molecular Panels

### RNA Panel (387 Genes)

**Gene Categories:**
- **Immune markers:** CD3D, CD4, CD8A, FOXP3, MS4A1 (CD20), PTPRC (CD45), CD68, ITGAX (CD11c)
- **Checkpoint:** PDCD1 (PD1), CD274 (PDL1), LAG3, HAVCR2 (TIM3), TIGIT
- **Epithelial:** EPCAM, KRT8, KRT18, KRT19, CDX2, MUC2, MUC5AC
- **Stromal:** ACTA2, COL1A1, FAP, PDGFRA, PDGFRB
- **Gastric markers:** TFF1, TFF2, TFF3, MUC6, PGA3, PGC
- **Cytokines/Chemokines:** IL6, CXCL9, CXCL10, CXCL13, CCL18, IFNG

### Protein Panel (17 Markers)

| Marker | Target | Cell Type Association | Mean SNR |
|--------|--------|----------------------|----------|
| **CD3** | T cells (pan) | All T lymphocytes | High |
| **CD4** | Helper T cells | CD4+ T cells, some macrophages | High |
| **CD8** | Cytotoxic T cells | CD8+ T cells | High |
| **FOXP3** | Tregs | Regulatory T cells | Moderate |
| **CD20** | B cells | B lymphocytes | High |
| **CD45** | Pan-leukocyte | All immune cells | High |
| **CD68** | Macrophages | Macrophages, monocytes | High |
| **CD11c** | Dendritic cells | DCs, some macrophages | High |
| **HLA-DR** | APCs | Antigen-presenting cells | High |
| **CD31** | Endothelial | Blood vessels | High |
| **PanCK** | Epithelial | Epithelial/tumor cells | **Highest (96.0)** |
| **aSMA** | Stromal | Myofibroblasts, CAFs | High |
| **KI67** | Proliferation | Dividing cells | Moderate |
| **PD1** | Checkpoint | Exhausted T cells | Moderate |
| **PDL1** | Checkpoint | Tumor/immune cells | Moderate |
| **ATPase** | Metabolic | Various | Variable |
| **Isotype** | Negative control | None (background) | Low |

---

## Multimodal Integration Results

### WNN (Weighted Nearest Neighbor) Integration

| Sample | RNA Weight | Protein Weight | WNN Clusters | % Flagged (Admixture) |
|--------|----------:|---------------:|-------------:|----------------------:|
| A01 | 10.3% | 89.7% | 11 | 24.2% |
| B01 | 11.6% | 88.4% | 11 | 48.5% |
| C01 | 11.5% | 88.5% | 12 | 45.5% |
| D01 | 11.0% | 89.0% | 11 | 48.8% |
| E01 | 11.3% | 88.7% | 10 | 40.8% |
| F01 | 12.7% | 87.3% | 11 | 52.9% |
| G01 | 10.7% | 89.3% | 11 | 33.9% |
| H01 | 11.5% | 88.5% | 11 | 52.1% |

**Key Insight:** Protein data dominates integration (~89% weight) due to higher signal clarity.

### RNA-Protein Correlation

| Marker | Spearman r | Interpretation |
|--------|----------:|----------------|
| HLA-DR | **0.31** | Best correlation |
| aSMA | 0.18 | Good correlation |
| CD68 | 0.17 | Good correlation |
| CD45 | 0.15 | Moderate |
| PanCK | 0.12 | Moderate |
| CD3 | 0.08 | Low |
| FOXP3 | 0.01 | Very low |
| KI67 | **0.003** | No correlation |

**Note:** Low correlations suggest post-transcriptional regulation and protein stability differences.

---

## Cell Type Annotation

### Hierarchical Annotation Strategy

```
Level 1: Lineage
├── Epithelial (PanCK+)
├── Immune (CD45+)
│   ├── T cells (CD3+)
│   │   ├── CD8+ T cells (CD8+)
│   │   │   ├── Exhausted (PD1+, LAG3+, TIGIT+)
│   │   │   └── Active
│   │   ├── CD4+ T cells (CD4+, FOXP3-)
│   │   └── Tregs (CD4+, FOXP3+)
│   ├── B cells (CD20+)
│   ├── Macrophages (CD68+)
│   └── Dendritic cells (CD11c+, HLA-DR+)
├── Stromal (aSMA+ or COL1A1+)
│   └── Fibroblasts / CAFs
└── Endothelial (CD31+)
```

### Cell Type Distribution (Lane 1)

| Cell Type | Total Cells | % of Total |
|-----------|------------:|-----------:|
| Epithelial | 106,726 | 18.7% |
| Fibroblast/Stromal | 69,626 | 12.2% |
| Endothelial | 70,653 | 12.4% |
| Macrophages | 14,073 | 2.5% |
| CD8+ T cells | 19,422 | 3.4% |
| CD4+ T cells | 11,405 | 2.0% |
| B cells | 16,470 | 2.9% |
| Tregs | 4,523 | 0.8% |
| Dendritic cells | 8,742 | 1.5% |
| Other Immune | 23,517 | 4.1% |

---

## Immune Microenvironment Analysis

### Per-Sample Immune Metrics

| Sample | Total Immune % | CD8 % | Treg % | CD8/Treg Ratio | PD1+ T % | PDL1+ Tumor % | CD8 Exhausted % |
|--------|---------------:|------:|-------:|---------------:|---------:|--------------:|----------------:|
| **A01** | 20.0% | 0.94% | 0.60% | 1.54 | 35.4% | 18.8% | 25.0% |
| **B01** | 20.0% | 3.13% | 0.07% | **39.3** | 8.8% | 19.7% | 25.0% |
| **C01** | 20.0% | 2.69% | 0.68% | 3.92 | 15.5% | 34.4% | 25.0% |
| **D01** | 20.0% | 3.28% | 0.50% | 6.49 | 28.5% | 22.6% | 25.0% |
| **E01** | 20.0% | 3.85% | 0.77% | 4.96 | 26.0% | 21.8% | 25.0% |
| **F01** | 20.0% | 3.78% | 1.43% | 2.63 | 35.0% | 18.9% | 25.0% |
| **G01** | 20.0% | 3.50% | 1.77% | 1.96 | 40.0% | 28.5% | 25.0% |
| **H01** | 20.0% | 4.39% | 0.82% | 5.31 | 20.7% | 26.8% | 25.0% |

### Sample Groupings

| Group | Criteria | Samples |
|-------|----------|---------|
| **High CD8/Treg** | Ratio > 5 | B01, D01, H01 |
| **Low CD8/Treg** | Ratio < 2 | A01, F01, G01 |
| **High PDL1** | >25% PDL1+ tumor | C01, G01, H01 |
| **Low PDL1** | <20% PDL1+ tumor | A01, B01, F01 |
| **High Exhaustion** | >25% exhausted CD8 | A01, F01, G01 |

---

## Spatial Analysis Results

### Neighborhood Enrichment Z-scores

**Strong Co-localization (z > 3):**
- CD8 T ↔ Other Immune: z = 14.61
- CD4 T ↔ Tregs: z = 5.49
- B cells ↔ CD8 T cells: z = 4.37

**Strong Avoidance (z < -3):**
- Epithelial ↔ Stromal: z = -17.0
- Epithelial ↔ CD8 T: z = -10.89
- Epithelial ↔ Macrophages: z = -5.14

**Self-Aggregation:**
- Epithelial: z = 44.25 (highest)
- Endothelial: z = 28.50
- Stromal: z = 21.84

---

## Data Files

### Directory Structure

```
spatial-hackathon-2026/results/g4x_choi_batch2/
├── A01_analyzed.h5ad          # Raw analyzed (357 genes)
├── B01_analyzed.h5ad
├── ...
├── annotated_v2/
│   ├── A01_annotated.h5ad     # With cell type annotations
│   ├── ...
│   └── annotation_summary.csv
├── wnn_integrated/
│   ├── A01_wnn.h5ad           # WNN integrated
│   ├── ...
│   └── wnn_summary.csv
└── figures/
    └── *.png
```

### AnnData Object Structure

```python
adata.obs columns:
- Sample metadata: sample_id, lane
- QC metrics: n_genes, n_counts, total_counts
- Protein intensities: CD3_intensity_mean, CD8_intensity_mean, ...
- Spatial: cell_x, cell_y
- Clustering: leiden_0.5, leiden_1.0, leiden_wnn_0_5, leiden_wnn_1_0
- WNN: rna_weight, protein_weight, admixture_score, admixture_flag
- Annotation: lineage, cell_type, cell_type_detailed, confidence scores

adata.obsm keys:
- X_pca: PCA embedding
- X_umap: UMAP embedding
- protein: Protein expression matrix (17 markers)
- spatial: Spatial coordinates (x, y)

adata.var: 357 genes with standard QC metrics
```

---

## Known Limitations

1. **Panel Size:** 387 genes (vs ~20,000 whole transcriptome)
   - Some cell subtypes may be under-resolved
   - Limited to targeted gene panel

2. **RNA-Protein Discordance:** Mean r = 0.088
   - Protein-based annotation more reliable
   - RNA may capture transcriptional states not yet translated

3. **Admixture Flagging:** 24-53% cells flagged per sample
   - Mixed signals at cell boundaries
   - Conservative QC approach

4. **Clinical Metadata:** Progression stage (Normal/Metaplasia/Cancer) annotations pending from Choi Lab

---

## Citation

If using this data, please cite:
- Choi Lab (data generation)
- Resolve Biosciences G4X platform
- Analysis pipeline: Van Belkum M, 2026

---

## Contact

Max Van Belkum
MD-PhD Student, Vanderbilt University
GitHub: [@vanbelkummax](https://github.com/vanbelkummax)
