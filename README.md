# G4X Choi Batch 2 Analysis

Comprehensive multimodal spatial transcriptomics analysis of G4X data from the Choi Lab (Batch 2).

## Dataset Overview

| Metric | Value |
|--------|-------|
| **Platform** | G4X (imaging-based spatial transcriptomics) |
| **Software** | Caretta v25.08.0 |
| **Total Samples** | 32 (4 lanes × 8 samples) |
| **Total Cells** | 2,308,968 |
| **Total Transcripts** | 192,411,241 |
| **Transcript Panel** | 387 genes |
| **Protein Panel** | 17 markers (immune-focused) |

### Protein Panel (17 markers)

| Marker | Cell Type Association |
|--------|----------------------|
| CD3 | T cells (pan) |
| CD4 | Helper T cells |
| CD8 | Cytotoxic T cells |
| FOXP3 | Regulatory T cells |
| CD20 | B cells |
| CD45 | Pan-leukocyte |
| CD68 | Macrophages |
| CD11c | Dendritic cells |
| HLA-DR | Antigen-presenting cells |
| CD31 | Endothelial cells |
| PanCK | Epithelial cells |
| aSMA | Stromal/myofibroblasts |
| KI67 | Proliferating cells |
| PD1 | Checkpoint (T cells) |
| PDL1 | Checkpoint (tumor/immune) |
| ATPase | Metabolic |
| Isotype | Negative control |

### Sample Summary

| Sample | Cells | Med. Trans/Cell | Med. Genes/Cell |
|--------|------:|----------------:|----------------:|
| A01 | 127,792 | 42 | 28 |
| B01 | 30,194 | 69 | 36 |
| C01 | 38,588 | 47 | 29 |
| D01 | 31,711 | 64 | 35 |
| E01 | 34,353 | 60 | 35 |
| F01 | 62,875 | 34 | 26 |
| G01 | 158,727 | 81 | 46 |
| H01 | 85,842 | 84 | 42 |
| A02 | 120,260 | 48 | 32 |
| B02 | 37,511 | 64 | 38 |
| C02 | 52,588 | 55 | 35 |
| D02 | 56,035 | 48 | 32 |
| E02 | 31,790 | 51 | 33 |
| F02 | 33,857 | 56 | 35 |
| G02 | 79,634 | 89 | 50 |
| H02 | 87,215 | 109 | 55 |
| A03 | 145,070 | 48 | 31 |
| B03 | 38,095 | 57 | 33 |
| C03 | 21,231 | 47 | 30 |
| D03 | 33,227 | 48 | 31 |
| E03 | 86,003 | 63 | 41 |
| F03 | 82,704 | 50 | 33 |
| G03 | 45,583 | 67 | 40 |
| H03 | 35,585 | 126 | 36 |
| A04 | 130,891 | 65 | 38 |
| B04 | 119,211 | 67 | 29 |
| C04 | 89,387 | 74 | 33 |
| D04 | 169,683 | 75 | 32 |
| E04 | 68,014 | 67 | 34 |
| F04 | 78,831 | 62 | 30 |
| G04 | 64,147 | 53 | 29 |
| H04 | 32,334 | 25 | 11 |

## Analysis Pipeline

### 1. Quality Control
- Cell/sample QC metrics
- Protein signal-to-noise ratio
- Transcript quality assessment

### 2. Cell Type Annotation
- Protein marker-based gating
- Unsupervised clustering
- Cell type identification

### 3. Multimodal Integration
- RNA-protein correlation
- Cross-modality validation
- Concordance analysis

### 4. Spatial Analysis
- Neighborhood enrichment
- Cell-cell interactions
- Spatial niche identification

### 5. Hypothesis-Driven Analysis
- Immune infiltration patterns
- TME composition
- Checkpoint expression analysis

## Documentation

| Document | Description |
|----------|-------------|
| [📊 DATA_DOCUMENTATION.md](docs/DATA_DOCUMENTATION.md) | Comprehensive dataset documentation with sample metrics, QC, and cell type distributions |
| [📋 PROJECT_IDEAS_RANKED.md](docs/PROJECT_IDEAS_RANKED.md) | Ranked project ideas with feasibility/significance scoring |
| [📚 G4X_Project_Ideas_Evidence_Based.md](docs/G4X_Project_Ideas_Evidence_Based.md) | Literature-grounded project prioritization with Polymath/Vanderbilt references |
| [🖥️ G4X_Project_Ideas.html](docs/G4X_Project_Ideas.html) | Interactive React dashboard for project exploration |

## Repository Structure

```
g4x-choi-batch2-analysis/
├── scripts/           # Analysis scripts
├── data/              # Processed data
├── docs/              # Project documentation
│   ├── DATA_DOCUMENTATION.md
│   ├── PROJECT_IDEAS_RANKED.md
│   ├── G4X_Project_Ideas_Evidence_Based.md
│   └── G4X_Project_Ideas.html
├── results/
│   ├── figures/       # Publication figures
│   └── tables/        # Summary tables
├── notebooks/         # Exploratory notebooks
└── README.md
```

## Data Location

Raw data: `/mnt/x/Choi_Batch_2_Tuesday/`

## Requirements

- Python 3.10+
- scanpy, squidpy, anndata
- pandas, numpy, scipy
- matplotlib, seaborn

## Usage

```bash
# Run full analysis pipeline
python scripts/01_qc_analysis.py
python scripts/02_cell_type_annotation.py
python scripts/03_multimodal_integration.py
python scripts/04_spatial_analysis.py
python scripts/05_hypothesis_analysis.py
```

## Key Findings

### Quality Control
- **32 samples** analyzed, 2.3M total cells
- Mean 62.3 transcripts/cell, 34.3 genes/cell
- All protein markers show high SNR (>10), with PanCK highest (96.0)
- 2 samples flagged for low transcript detection (F01: 34, H04: 25)

### Cell Type Composition
| Cell Type | Mean % | Range |
|-----------|-------:|------:|
| Epithelial | 18.7% | 17.3-19.9% |
| Stromal | 12.2% | 5.2-16.1% |
| Endothelial | 12.4% | 9.4-15.6% |
| CD8+ T cells | 3.4% | 0.9-9.1% |
| B cells | 2.9% | 1.3-4.9% |
| Macrophages | 2.5% | 0.9-4.8% |
| CD4+ T cells | 2.0% | 0.6-3.4% |
| Tregs | 0.8% | 0.1-2.0% |

### Multimodal Integration (RNA-Protein Correlation)
- Mean Spearman r = 0.088 (typical for imaging platforms)
- **Best correlations:** HLA-DR (0.31), aSMA (0.18), CD68 (0.17)
- **Weakest correlations:** KI67 (0.003), FOXP3 (0.01), PD1 (0.02)
- Low correlation may indicate post-transcriptional regulation

### Spatial Analysis - Neighborhood Enrichment

**Strong co-localization (positive z-score):**
- B cells ↔ CD8 T cells (z=4.37)
- CD8 T ↔ Immune_other (z=14.61)
- CD4 T ↔ Tregs (z=5.49)
- All immune cells cluster together

**Strong avoidance (negative z-score):**
- Epithelial ↔ Stromal (z=-17.0)
- Epithelial ↔ Immune (z=-5 to -11)
- Clear tissue compartmentalization

**Self-aggregation:**
- All cell types show homotypic clustering
- Epithelial cells: z=44.25 (highest)
- Endothelial cells: z=28.50

### Immune Microenvironment
- **Total immune:** 20.0% across samples
- **CD8/Treg ratio:** 7.3 ± 9.5 (favorable in most samples)
- **PD1+ T cells:** 29.0% ± 13.2%
- **PDL1+ tumor:** 25.3% ± 9.0%
- **CD8 exhaustion:** 25% of CD8+ T cells show exhaustion phenotype

### Notable Samples
| Category | Samples | Characteristic |
|----------|---------|----------------|
| High CD8/Treg ratio | B01 (39.3), H03 (38.4) | Favorable immune |
| Low CD8/Treg ratio | A03 (0.9), A04 (1.3) | Treg-dominated |
| High PDL1 tumor | 11 samples >30% | Checkpoint expression |

## Figures

| Figure | Description |
|--------|-------------|
| `qc_overview.png` | Sample QC metrics |
| `protein_snr.png` | Protein signal-to-noise |
| `cell_type_composition.png` | Cell type proportions |
| `marker_expression_by_celltype.png` | Marker validation |
| `protein_rna_correlation.png` | Multimodal concordance |
| `neighborhood_enrichment.png` | Cell-cell interactions |
| `spatial_maps.png` | Spatial distributions |
| `immune_infiltration.png` | TME composition |
| `exhaustion_analysis.png` | T cell exhaustion |

## Author

Max Van Belkum, MD-PhD Student, Vanderbilt University

## Date

January 2026
