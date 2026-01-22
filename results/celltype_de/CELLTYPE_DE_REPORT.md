# Cell Type-Specific DE Analysis Report

Generated: 2026-01-22 01:54:34
Runtime: 0.8 minutes

## Overview

This analysis performs cell type-stratified differential expression analysis using
**pseudobulk aggregation** to avoid pseudoreplication. Rather than testing 1.8M
individual cells, we aggregate to 29 sample-level profiles before DE testing.

## Critical Notes

### Stage-Patient Confounding

| Patient | Stages Present |
|---------|----------------|
| SNU-105 | normal, metaplasia |
| SNU-107 | normal, metaplasia |
| SNU-484 | cancer |
| ctrl | control |

**Implication**: Normal and metaplasia samples come from only 2 patients; cancer
samples come from 1 patient. Results may partially reflect patient effects rather
than pure stage effects.

### Incomplete Gene Signatures

| Signature | Available | Missing |
|-----------|-----------|---------|
| iCAF | IL6, PDGFRA, FAP | CXCL12 |
| M2 Macrophage | CD163, MRC1, IL10 | ARG1 |

### Admixture Filtering

Cells with `is_admixed=True` were excluded from all analyses to ensure clean
cell type populations.

## Results Summary

### Epithelial Analysis

- Samples analyzed: 29
- Excluded samples: []

- N_vs_M: 0 significant genes (q<0.05)
- M_vs_C: 0 significant genes (q<0.05)
- N_vs_C: 0 significant genes (q<0.05)

### CAF Subtype Analysis

⚠️ **iCAF signature incomplete** (missing CXCL12)

- Samples analyzed: 29

### CD8 T Cell Exhaustion

- Samples analyzed: 29
- Excluded samples (low cell count): []

### Macrophage Polarization

⚠️ **M2 signature incomplete** (missing ARG1)

- Samples analyzed: 29


## Output Files

```
results/celltype_de/
├── epithelial/
│   ├── de_N_vs_M.csv, de_M_vs_C.csv, de_N_vs_C.csv
│   ├── gastric_markers_by_stage.csv
│   └── sample_metadata.csv
├── caf/
│   ├── caf_subtype_scores.csv
│   ├── caf_proportions_by_stage.csv
│   └── NOTE_iCAF_incomplete.txt
├── cd8/
│   ├── exhaustion_scores.csv
│   ├── exhaustion_markers_detail.csv
│   └── samples_excluded.txt
├── macrophage/
│   ├── m1_m2_scores.csv
│   └── NOTE_M2_incomplete.txt
├── figures/
│   ├── fig1_epithelial_volcano.png
│   ├── fig2_caf_composition.png
│   ├── fig3_cd8_exhaustion.png
│   ├── fig4_macrophage_polarization.png
│   └── fig5_heatmap_top_de.png
├── sensitivity/
│   └── admixture_comparison.csv
└── CELLTYPE_DE_REPORT.md
```

## Methods

### Pseudobulk DE

1. Subset to cell type (excluding admixed cells)
2. Filter samples with < MIN_CELLS threshold
3. Aggregate counts to sample level (sum)
4. Normalize (CPM + log1p)
5. Mann-Whitney U test (non-parametric, robust for small n)
6. Benjamini-Hochberg FDR correction

### Signature Scoring

Used scanpy's `score_genes` function with mean expression of signature genes
vs. randomly sampled control genes.
