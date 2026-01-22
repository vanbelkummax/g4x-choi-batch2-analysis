# Protein-RNA Spatial Correlation Report

Generated: 2026-01-22 01:59:20
Runtime: 4.3 minutes

## Overview

This analysis compares protein immunofluorescence intensity with RNA expression
for 15 matched protein-gene pairs across all samples.

| Metric | Value |
|--------|-------|
| Total protein-sample pairs | 435 |
| Valid pairs analyzed | 406 |
| Samples | 29 |
| Proteins | 14 |
| Zero-variance cases | 12 |

## Protein-Gene Mapping

| Protein | RNA Gene | Status |
|---------|----------|--------|
| CD4 | CD4 | Exact match |
| CD68 | CD68 | Exact match |
| FOXP3 | FOXP3 | Exact match |
| CD3 | CD3D | Primary chain |
| CD8 | CD8A | Primary chain |
| HLA-DR | HLA-DRA | MHC-II alpha |
| KI67 | MKI67 | Standard name |
| PD1 | PDCD1 | Standard name |
| PDL1 | CD274 | Standard name |
| aSMA | ACTA2 | Standard name |
| CD31 | PECAM1 | Standard name |
| **PanCK** | **EPCAM** | **⚠️ PROXY** |
| CD20 | MS4A1 | B cell marker |
| CD45 | PTPRC | Pan-immune |
| CD11c | ITGAX | DC marker |

**Note on PanCK→EPCAM**: Cytokeratins (KRT8, KRT18, KRT19) are not in the RNA panel.
EPCAM is used as an epithelial marker proxy. Correlations should be interpreted with caution.

## Overall Correlation Statistics

### Cell-level Pearson Correlation

| Statistic | Value |
|-----------|-------|
| Mean | 0.095 |
| Std Dev | 0.111 |
| Min | -0.056 |
| Max | 0.435 |
| Median | 0.043 |

### Spatial Similarity (SSIM)

| Statistic | Value |
|-----------|-------|
| Mean | 0.501 |
| Std Dev | 0.164 |
| Min | 0.020 |
| Max | 0.807 |
| Median | 0.536 |

## Zero-Variance Cases

| Sample | Protein | Issue |
|--------|---------|-------|
| A02 | PanCK | zero_variance_rna |
| B02 | PanCK | zero_variance_rna |
| B03 | PanCK | zero_variance_rna |
| C03 | PanCK | zero_variance_rna |
| D01 | PanCK | zero_variance_rna |
| D03 | KI67 | zero_variance_rna |
| D03 | PanCK | zero_variance_rna |
| E01 | PanCK | zero_variance_rna |
| F02 | PanCK | zero_variance_rna |
| G03 | PanCK | zero_variance_rna |
| G04 | PanCK | zero_variance_rna |
| H03 | PanCK | zero_variance_rna |


## Output Files

```
results/protein_rna_correlation/
├── sample_panels/           # Individual sample-protein panels
├── protein_summary/         # Per-protein summary figures
├── correlation_details.csv  # All metrics for all pairs
├── correlation_matrix.csv   # Sample × Protein Pearson r
├── ssim_matrix.csv          # Sample × Protein SSIM
├── zero_variance_cases.csv  # Cases with no correlation
├── summary_heatmap.png      # Visual correlation matrix
└── PROTEIN_RNA_REPORT.md
```

## Methods

### Correlation Metrics

1. **Cell-level Pearson r**: Direct correlation between protein intensity and RNA expression for all cells
2. **Cell-level Spearman r**: Rank-based correlation (robust to outliers)
3. **Grid Pearson r**: Correlation on spatially-binned averages (50×50 grid)
4. **SSIM**: Structural Similarity Index measuring spatial pattern similarity

### Normalization

- Protein: Raw intensity values from immunofluorescence
- RNA: Normalized counts from merged_corrected.h5ad
- For SSIM: Both modalities z-score normalized before comparison

### Grid Binning

- Grid size: 50×50 bins
- Binning statistic: Mean value per bin
- Empty bins: Set to 0
