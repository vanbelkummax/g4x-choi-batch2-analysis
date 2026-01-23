# G4X Gastric Cancer Progression Analysis

## Overview

This report summarizes a spatial transcriptomics analysis of gastric cancer progression using the G4X platform. We analyzed three tissue samples from the same patient representing the progression from normal gastric tissue to cancer.

---

## The Data

| Sample | Tissue Type | Cells | Description |
|--------|-------------|------:|-------------|
| E02 | Normal | 27,775 | Normal gastric mucosa |
| F02 | Metaplasia | 30,229 | Intestinal metaplasia (precancerous) |
| G02 | Cancer | 76,463 | Gastric adenocarcinoma |
| **Total** | | **134,467** | |

**Platform:** G4X targeted spatial transcriptomics (337 genes)

---

## What We Did

### Phase 1: Initial Analysis (Marker-Based)

We first performed standard single-cell analysis:

1. **Quality Control** - Removed low-quality cells
2. **Clustering Optimization** - Tested multiple parameter combinations to find optimal clustering
3. **Cell Type Annotation** - Used curated gastric-specific markers to identify cell types
4. **Differential Expression** - Compared gene expression between disease stages
5. **Spatial Statistics** - Analyzed how cell types are organized in tissue

**Key findings from Phase 1:**
- Cancer sample shows better cluster separation (silhouette 0.193) than normal (0.11)
- Top upregulated genes in cancer: **COL1A1** (3.5×), **MUC2** (3.0×), **CLDN3** (2.9×)
- Loss of normal gastric markers: PGC, ATP4A (parietal/chief cell markers)
- Increased stromal infiltration in cancer

### Phase 2: RCTD Deconvolution

We then applied RCTD (Robust Cell Type Decomposition), a reference-based method that transfers cell type labels from a published single-cell RNA-seq atlas to our spatial data.

**Reference used:** GSE134520 (Zhang et al. 2019)
- 44,012 cells from gastric tissue
- Conditions: Normal, Chronic Atrophic Gastritis, Intestinal Metaplasia, Early Gastric Cancer
- 18 cell types annotated

**RCTD Results:**
- 117,887 cells (87.7%) successfully annotated
- 18 cell types identified
- 14% rejection rate (cells that couldn't be confidently assigned)

### Phase 3: Validation and Scoring

We computed additional scores to validate findings:

1. **Intestinal Metaplasia (IM) Score** - Based on CDX2, MUC2, TFF3 expression
2. **Hwang Cancer Signature** - 32-gene signature from published gastric cancer studies

---

## Key Findings

### Finding 1: Massive Intestinal Metaplasia in Cancer (MUC2 Expression)

MUC2 is a goblet cell marker and indicator of intestinal metaplasia.

| Stage | MUC2+ Cells | Total Cells | Percentage |
|-------|------------:|------------:|-----------:|
| Normal | 2,916 | 27,775 | 10.5% |
| Metaplasia | 4,249 | 30,229 | 14.1% |
| Cancer | 52,205 | 76,463 | **68.3%** |

**Interpretation:** The cancer sample shows 68% MUC2+ cells compared to only 10% in normal tissue. This represents widespread intestinal metaplasia, where normal gastric epithelium is replaced by intestinal-type cells. This aligns with our Phase 1 DEG analysis showing MUC2 as a top upregulated gene.

---

### Finding 2: RCTD "Goblet" Annotation Fails

| Stage | RCTD "Goblet" | Total Cells | Percentage |
|-------|------------:|------------:|-----------:|
| Normal | 265 | 27,775 | 0.95% |
| Metaplasia | 276 | 30,229 | 0.91% |
| Cancer | 799 | 76,463 | 1.04% |

**Interpretation:** RCTD labels only ~1% of cells as "Goblet" across all stages, completely missing the 68% MUC2+ cells in cancer. This failure occurs because the reference dataset contained only 623 goblet cells (1.4% of 44,012), providing insufficient training data for accurate deconvolution.

**Recommendation:** Do not use RCTD "Goblet" labels. Use direct MUC2 expression instead.

---

### Finding 3: RCTD "Cancer" Is Not a Cancer Indicator

| Stage | RCTD "Cancer" | Total Cells | Percentage |
|-------|------------:|------------:|-----------:|
| Normal | 3,069 | 27,775 | 13.73% |
| Metaplasia | 3,257 | 30,229 | 13.24% |
| Cancer | 8,756 | 76,463 | 12.34% |

**Interpretation:** The RCTD "Cancer" cell type is stable at ~13% across ALL samples, including normal tissue. This label does not indicate malignancy. It represents cells matching the transcriptional profile of cells from the reference's Early Gastric Cancer (EGC) and Intestinal Metaplasia-Severe samples. The reference "Cancer" cells (4,377 total) primarily came from pre-malignant conditions, not true cancer tissue.

**Recommendation:** Do not use RCTD "Cancer" labels for malignancy detection.

---

### Finding 4: Hwang Cancer Signature Validates Progression

The Hwang 32-gene cancer signature (29/32 genes present in our panel) provides a validated malignancy indicator.

| Stage | Cells | Mean Score | Status |
|-------|------:|----------:|--------|
| Normal | 27,775 | -0.052 | Baseline |
| Metaplasia | 30,229 | -0.051 | No change |
| Cancer | 76,463 | **+0.101** | Elevated |

**Interpretation:** Normal and metaplasia samples show similar baseline scores (~-0.05), while the cancer sample shows elevated scores (+0.101). This signature correctly captures cancer progression, unlike the RCTD "Cancer" cell type.

**Recommendation:** Use the Hwang cancer signature for malignancy assessment.

---

## Integration with Earlier Findings

Our RCTD and scoring analysis aligns with and extends the initial marker-based analysis:

| Finding | Phase 1 (Markers) | Phase 2 (RCTD/Scores) |
|---------|-------------------|----------------------|
| **MUC2 upregulation** | 3.0× fold-change in DEG | 68% MUC2+ cells in cancer |
| **Stromal infiltration** | COL1A1 3.5× upregulated | Fibroblast 10.4% of cells |
| **Loss of gastric identity** | PGC, ATP4A downregulated | Gastric markers decrease |
| **Cancer progression** | Distinct clustering | Hwang score +0.101 |

---

## Limitations

1. **N=1 Patient** - All samples from one individual; no biological replicates
2. **337-Gene Panel** - Limited to targeted genes, not whole transcriptome
3. **Reference Limitations** - GSE134520 had sparse goblet cells (1.4%)
4. **No Protein Validation** - G4X has IF protein data we haven't integrated

---

## Files Included

### Figures
| File | Description |
|------|-------------|
| `fig1_muc2_analysis.png` | MUC2 expression analysis with interpretation |
| `fig2_rctd_goblet_failure.png` | RCTD Goblet annotation failure |
| `fig3_rctd_cancer_misleading.png` | RCTD "Cancer" label analysis |
| `fig4_hwang_validated.png` | Hwang cancer signature validation |
| `fig_combined_analysis.png/pdf` | Combined 4-panel summary |
| `fig_integrated_panel.png` | Spatial visualization panel |
| `fig_qc_panel.png` | QC metrics panel |
| `fig_celltype_focus.png` | Top cell types spatial distribution |

### Data Files
| File | Description |
|------|-------------|
| `g4x_cell_labels.csv` | Cell ID and RCTD label (2 columns) |
| `g4x_cell_labels_enhanced.csv` | Labels with MUC2 expression |
| `g4x_rctd_annotations.csv` | Full RCTD results |
| `g4x_rctd_weights.csv` | RCTD cell type weights |
| `g4x_im_scores.csv` | IM and Hwang scores per cell |
| `pilot_final_summary.csv` | Combined summary table |

### Scripts
| Script | Purpose |
|--------|---------|
| `00_load_pilot_samples.py` | Load and QC G4X data |
| `01_hvg_optimize_parallel.py` | Optimize clustering parameters |
| `02_annotate_multimethod.py` | Multi-method cell type annotation |
| `03_spatial_viz.py` | Generate spatial visualizations |
| `04_exploratory_deg.py` | Differential expression analysis |
| `05_spatial_statistics.py` | Neighborhood enrichment |
| `06a_build_reference.py` | Build GSE134520 reference |
| `07_rctd_deconvolution.R` | Run RCTD deconvolution |
| `08_im_scoring.py` | Compute IM and Hwang scores |
| `09_merge_results.py` | Merge all annotations |
| `10_publication_figures.py` | Generate final figures |

---

## Conclusions

1. **MUC2 expression** is the most reliable indicator of intestinal metaplasia, showing a dramatic increase from 10.5% (normal) to 68.3% (cancer).

2. **RCTD deconvolution has limitations** for this dataset due to reference imbalances:
   - "Goblet" annotation misses 98% of actual goblet/IM cells
   - "Cancer" annotation does not indicate malignancy

3. **The Hwang 32-gene cancer signature** correctly captures cancer progression and should be used for malignancy assessment.

4. **The cancer sample shows widespread intestinal metaplasia**, consistent with the well-established model of gastric carcinogenesis through the metaplasia-dysplasia-carcinoma sequence.

---

## Methods Summary

- **Platform:** G4X targeted spatial transcriptomics (337 genes)
- **Cells analyzed:** 134,467 across 3 samples
- **Reference for RCTD:** GSE134520 (44,012 cells, 18 cell types)
- **Software:** Python (scanpy, squidpy), R (spacexr/RCTD)
- **Date:** January 2026

---

*Analysis performed in conda environment `enact`*
*Repository: `/home/user/g4x-choi-batch2-analysis/`*
