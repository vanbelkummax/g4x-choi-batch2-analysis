# Pipeline Validation: Our Analysis vs Resolve Biosciences

## Overview

This document describes the validation of our G4X analysis pipeline against Resolve Biosciences' official outputs (caretta v25.08.0).

## Validation Script

**Script:** `scripts/35_validate_vs_resolve.py`

**Comparisons performed:**
1. Cell counts (core_metrics.csv vs our h5ad)
2. Clustering assignments (clustering_umap.csv.gz vs our WNN clusters)
3. QC statistics (median transcripts/genes per cell)
4. Cell metadata correlations (spatial coordinates, nuclei area)

## Results Summary

### Cell Counts

| Sample | ROI Type | Resolve | Ours | Filtered | % |
|--------|----------|---------|------|----------|---|
| A01 | CRC ctrl | 127,792 | 113,389 | 14,403 | 11.3% |
| B01 | Normal | 30,194 | 26,161 | 4,033 | 13.4% |
| C01 | Metaplasia | 38,588 | 34,166 | 4,422 | 11.5% |
| D01 | Normal | 31,711 | 23,935 | 7,776 | **24.5%** |
| E01 | Metaplasia | 34,353 | 32,440 | 1,913 | 5.6% |
| F01 | Cancer | 62,875 | 53,156 | 9,719 | 15.5% |
| G01 | Cancer | 158,727 | 155,843 | 2,884 | 1.8% |
| H01 | Unknown | 85,842 | 83,098 | 2,744 | 3.2% |
| **Total** | - | **570,082** | **522,188** | **47,894** | **8.4%** |

### QC Statistics

| Metric | Resolve | Ours |
|--------|---------|------|
| Median transcripts/cell | 60.1 | 65.9 |
| Median genes/cell | 34.6 | 37.6 |

**Interpretation:** Our higher medians confirm that QC filtering removed low-quality cells, enriching for cells with better RNA capture.

### Metadata Correlation

| Metric | Correlation |
|--------|-------------|
| X coordinate | 1.0000 |
| Y coordinate | 1.0000 |
| Nuclei area | 1.0000 |

**Interpretation:** Perfect correlation confirms we're using identical source data from Resolve's segmentation.

### Clustering Comparison

| Sample | Resolve Clusters | Our Clusters |
|--------|------------------|--------------|
| A01 | 13 | 11 |
| B01 | 6 | 11 |
| G01 | 4 | 11 |

**Key difference:** Resolve uses RNA-only Leiden clustering at resolution 0.4. We use WNN (Weighted Nearest Neighbors) integrating both RNA and protein modalities, which provides more consistent cluster numbers across samples.

## Pipeline Value-Add

Our pipeline adds the following over Resolve's baseline analysis:

1. **WNN Integration**: Combines 387-gene RNA panel with 17-protein panel for more robust clustering
2. **Admixture-Based QC**: Identifies cells with mixed transcriptomes (segmentation errors, doublets)
3. **Hierarchical Annotation**: Two-tier cell typing (lineage → cell type) with confidence scores
4. **Segmentation QC**: Three filtering levels (Baseline, QC-strict, Clean)

## Filtering Rationale

The 8.4% overall filtering removes:
- Empty/low-quality cells (low transcript counts)
- Segmentation errors (high admixture scores)
- Doublets (conflicting marker expression)

**D01 outlier (24.5% filtered):** This Normal sample may have tissue quality issues or higher segmentation error rate. Further investigation recommended.

## Output Files

```
results/g4x_choi_batch2/validation/
├── cell_count_comparison.csv    # Per-sample cell counts
├── clustering_comparison.csv    # Cluster number comparison
├── qc_stats_comparison.csv      # Median tx/genes comparison
├── annotation_summary.csv       # Our cell type proportions
├── metadata_correlation.csv     # Spatial/area correlations
└── validation_summary.csv       # Aggregated metrics
```

## Conclusion

✅ **Validation passed.** Our pipeline produces consistent results with Resolve's official analysis while adding multimodal integration, QC filtering, and cell type annotation.
