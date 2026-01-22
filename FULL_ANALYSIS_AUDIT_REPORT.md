# G4X Gastric Cancer Progression Analysis
## Full Audit Report

> ## 🚨 CRITICAL DATA QUALITY WARNING
>
> **The results in this report are PRELIMINARY and NOT TRUSTWORTHY.**
>
> - **CRITICAL BUG:** `merged_counts.h5ad` contains log-normalized data, NOT raw counts
> - **Zero significant findings** across all analyses
> - **Insufficient statistical power** (N=7-10 per stage)
> - **Reanalysis from true raw counts required**
>
> See [`DATA_QUALITY_NOTICE.md`](DATA_QUALITY_NOTICE.md) and [`CRITICAL_AUDIT_RAW_COUNTS_BUG.md`](CRITICAL_AUDIT_RAW_COUNTS_BUG.md) for details.
>
> **Trustworthy data:** `results/qc_all_samples/raw/*.h5ad` (individual sample raw counts)

---

**Project:** Choi_GC_preGC_batch2
**Platform:** G4X (Resolve Biosciences)
**Analysis Period:** January 2026
**Repository:** https://github.com/vanbelkummax/g4x-choi-batch2-analysis
**Current Commit:** `70141c0`

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Dataset Description](#2-dataset-description)
3. [Quality Control Pipeline](#3-quality-control-pipeline)
4. [Analysis Scripts & Results](#4-analysis-scripts--results)
5. [Numerical Results Summary](#5-numerical-results-summary)
6. [Statistical Significance Assessment](#6-statistical-significance-assessment)
7. [Limitations & Caveats](#7-limitations--caveats)
8. [Recommended Next Steps](#8-recommended-next-steps)
9. [File Manifest](#9-file-manifest)
10. [Reproducibility](#10-reproducibility)

---

## 1. Executive Summary

### Key Findings

| Analysis | Result | Statistical Significance |
|----------|--------|-------------------------|
| Bulk stage separation (PCA) | Stages overlap | p = 0.58 (Kruskal-Wallis) |
| Epithelial DE (N→M→C) | 0/341 genes | ❌ None at q < 0.05 or q < 0.10 |
| CAF subtype shifts | No trend | ❌ Not significant |
| CD8 exhaustion by stage | No trend | ❌ Not significant |
| Macrophage M1/M2 ratio | No trend | ❌ Not significant |
| Protein-RNA correlation | r = 0.095 | ⚠️ Poor concordance |

### Bottom Line

**No statistically significant transcriptional differences were detected between disease stages (Normal → Metaplasia → Cancer) using pseudobulk approaches.** This is attributable to:

1. Stage-patient partial confounding (though 2 patients have all 3 stages)
2. Insufficient power with N=29 samples for pseudobulk DE
3. Bulk measurements masking cell-level heterogeneity

### Positive Findings

- QC pipeline successfully processed 1.83M cells across 29 samples
- Batch correction (Harmony) improved integration (LISI: 2.46 → 2.66)
- Cell type annotation identified 14 populations
- Protein-RNA spatial patterns show moderate overlap (SSIM ~0.50)
- **Two patients (SNU-105, SNU-107) have all three disease stages**, enabling within-patient comparisons

---

## 2. Dataset Description

### Platform Specifications

| Parameter | Value |
|-----------|-------|
| Platform | G4X (Resolve Biosciences) |
| RNA genes | 387 (341 after QC filtering) |
| Protein markers | 17 (excluding controls) |
| Resolution | Single-cell spatial |

### Sample Summary

| Metric | Count |
|--------|-------|
| Total samples loaded | 32 |
| Samples passing QC | 29 |
| Samples excluded | 3 (C02, H02, H04) |
| Total cells | 1,835,026 |

### Samples Per Stage

| Stage | N Samples | N Cells | % Total |
|-------|-----------|---------|---------|
| Cancer | 10 | 659,029 | 35.9% |
| Control | 4 | 439,670 | 24.0% |
| Metaplasia | 7 | 347,862 | 19.0% |
| Normal | 8 | 388,465 | 21.2% |

### Samples Per Patient

| Patient | N Cells | Stages Present |
|---------|---------|----------------|
| SNU-105 | 558,052 | Normal, Metaplasia, Cancer |
| SNU-107 | 727,282 | Normal, Metaplasia, Cancer |
| SNU-484 | 110,022 | Cancer only |
| ctrl | 439,670 | Control only |

### Stage-Patient Distribution (Cell Counts)

|         | Cancer | Control | Metaplasia | Normal |
|---------|--------|---------|------------|--------|
| SNU-105 | 241,149 | 0 | 127,313 | 189,590 |
| SNU-107 | 307,858 | 0 | 220,549 | 198,875 |
| SNU-484 | 110,022 | 0 | 0 | 0 |
| ctrl | 0 | 439,670 | 0 | 0 |

### Cell Type Distribution

| Cell Type | N Cells | % Total |
|-----------|---------|---------|
| Epithelial | 665,200 | 36.3% |
| Stromal | 376,855 | 20.5% |
| Endothelial | 280,361 | 15.3% |
| PDL1_positive | 256,298 | 14.0% |
| PD1_positive | 87,572 | 4.8% |
| Proliferating | 48,737 | 2.7% |
| Macrophage | 31,159 | 1.7% |
| Immune (other) | 27,518 | 1.5% |
| B_cell | 19,973 | 1.1% |
| DC | 19,115 | 1.0% |
| CD8_T | 11,131 | 0.6% |
| CD4_T | 4,655 | 0.3% |
| T_cell | 4,348 | 0.2% |
| Treg | 2,104 | 0.1% |

### Protein Panel

| Marker | Target | Available in RNA? |
|--------|--------|-------------------|
| CD4 | T helper cells | ✓ CD4 |
| CD8 | Cytotoxic T cells | ✓ CD8A |
| CD3 | Pan-T cells | ✓ CD3D |
| CD68 | Macrophages | ✓ CD68 |
| FOXP3 | Tregs | ✓ FOXP3 |
| CD20 | B cells | ✓ MS4A1 |
| CD45 | Pan-immune | ✓ PTPRC |
| CD11c | Dendritic cells | ✓ ITGAX |
| HLA-DR | MHC-II | ✓ HLA-DRA |
| PD1 | Checkpoint | ✓ PDCD1 |
| PDL1 | Checkpoint ligand | ✓ CD274 |
| KI67 | Proliferation | ✓ MKI67 |
| aSMA | CAFs/smooth muscle | ✓ ACTA2 |
| CD31 | Endothelial | ✓ PECAM1 |
| PanCK | Epithelial | ⚠️ EPCAM (proxy) |
| ATPase | Metabolic | ✗ Not mapped |
| Isotype | Control | ✗ Control |

---

## 3. Quality Control Pipeline

### Scripts Executed

| Script | Purpose | Runtime | Output |
|--------|---------|---------|--------|
| `60_load_all_samples.py` | Load 32 raw samples | ~20 min | `results/qc_all_samples/raw/` |
| `61_comprehensive_qc.py` | Sample-level QC | ~30 min | QC metrics, exclusion flags |
| `62_process_qc_passing.py` | Cell-level QC, WNN | ~3 hr | Processed h5ad files |
| `63_merge_and_batch_correct.py` | Merge + Harmony | ~45 min | `merged_corrected.h5ad` |

### QC Metrics Applied

**Sample-level exclusion criteria:**
- Samples with < 1,000 cells
- Samples with > 50% mitochondrial reads
- Technical failures (H04 excluded)

**Cell-level filtering:**
- Minimum genes per cell: 50
- Maximum genes per cell: 5,000
- Mitochondrial fraction < 20%

### Batch Correction Assessment

| Metric | Pre-Harmony | Post-Harmony | Interpretation |
|--------|-------------|--------------|----------------|
| LISI score | 2.462 | 2.662 | ↑ Improved mixing |
| Silhouette | -0.0013 | -0.0060 | ↓ Less batch separation |

### Excluded Samples

| Sample | Reason |
|--------|--------|
| C02 | Technical failure |
| H02 | Low cell count |
| H04 | QC failure |

---

## 4. Analysis Scripts & Results

### 4.1 PCA Deep Dive (Script 46)

**Script:** `scripts/46_pca_deep_dive.py`
**Output:** `results/pca_deep_dive/`

#### Variance Decomposition

| Factor | Variance Explained (PC1) | Type |
|--------|-------------------------|------|
| Cell type | 19.2% | Biological |
| Lane | 5.6% | Technical |
| Stage | 1.3% | Biological |
| Patient | 1.2% | Technical |

#### Progression Statistics (Pseudobulk PCA)

| Test | Statistic | P-value | Interpretation |
|------|-----------|---------|----------------|
| Kruskal-Wallis (PC1 by stage) | 1.97 | **0.578** | No stage separation |
| Spearman correlation (PC1 vs progression) | ρ = -0.082 | **0.673** | No monotonic trend |

**Conclusion:** Disease stages do not separate in principal component space at the sample level.

---

### 4.2 Cell Type-Specific DE (Script 47)

**Script:** `scripts/47_celltype_specific_de.py`
**Output:** `results/celltype_de/`

#### Approach
- Pseudobulk aggregation: 1.83M cells → 29 sample-level profiles
- Statistical test: Mann-Whitney U (two-group comparisons)
- Multiple testing correction: Benjamini-Hochberg FDR

#### Differential Expression Results

| Comparison | Genes Tested | Significant (q < 0.05) | Significant (q < 0.10) |
|------------|--------------|------------------------|------------------------|
| Normal → Metaplasia | 341 | **0** | **0** |
| Metaplasia → Cancer | 341 | **0** | **0** |
| Normal → Cancer | 341 | **0** | **0** |

#### Top Genes by Uncorrected P-value (N vs C)

| Gene | log2FC | P-value (raw) | Q-value (BH) |
|------|--------|---------------|--------------|
| (All genes q > 0.10) | - | - | - |

**Note:** Even the top-ranked genes do not approach significance after FDR correction.

#### CAF Subtype Analysis

| Subtype | Signature Genes | Mean Score (N) | Mean Score (M) | Mean Score (C) |
|---------|-----------------|----------------|----------------|----------------|
| mCAF | ACTA2, TAGLN | 0.42 | 0.45 | 0.44 |
| iCAF | IL6, PDGFRA, FAP | 0.31 | 0.33 | 0.32 |
| apCAF | CD74, HLA-DRA | 0.28 | 0.27 | 0.29 |

**Note:** iCAF signature incomplete (CXCL12 not in panel)

#### CD8 Exhaustion Analysis

- Cells analyzed: 11,131 CD8 T cells
- Exhaustion markers: PDCD1, LAG3, HAVCR2, TIGIT, CTLA4

| Stage | Mean Exhaustion Score | Std |
|-------|----------------------|-----|
| Normal | 0.18 | 0.12 |
| Metaplasia | 0.19 | 0.11 |
| Cancer | 0.20 | 0.13 |

**Statistical test:** Kruskal-Wallis p > 0.05 (not significant)

#### Macrophage Polarization

- Cells analyzed: 31,159 macrophages
- M1 markers: CD80, TNF, IL1B
- M2 markers: CD163, MRC1, IL10 (ARG1 missing)

| Stage | M1 Score | M2 Score | M1/M2 Ratio |
|-------|----------|----------|-------------|
| Normal | 0.24 | 0.31 | 0.77 |
| Metaplasia | 0.25 | 0.30 | 0.83 |
| Cancer | 0.26 | 0.29 | 0.90 |

**Statistical test:** Kruskal-Wallis p > 0.05 (not significant)

---

### 4.3 Protein-RNA Spatial Correlation (Script 48)

**Script:** `scripts/48_protein_rna_spatial_correlation.py`
**Output:** `results/protein_rna_correlation/`

#### Overall Statistics

| Metric | N | Mean | Std | Min | Max | Median |
|--------|---|------|-----|-----|-----|--------|
| Cell Pearson r | 394 | **0.095** | 0.111 | -0.056 | 0.435 | 0.043 |
| Cell Spearman ρ | 394 | **0.087** | 0.108 | -0.090 | 0.472 | 0.037 |
| Spatial SSIM | 406 | **0.501** | 0.164 | 0.020 | 0.807 | 0.536 |

#### Per-Protein Correlations

| Protein | Gene | Mean Pearson r | Mean SSIM | N Valid Samples |
|---------|------|----------------|-----------|-----------------|
| aSMA | ACTA2 | **0.264** | 0.644 | 29 |
| HLA-DR | HLA-DRA | **0.262** | 0.598 | 29 |
| CD31 | PECAM1 | **0.206** | 0.630 | 29 |
| CD68 | CD68 | 0.145 | 0.589 | 29 |
| CD8 | CD8A | 0.106 | 0.537 | 29 |
| CD4 | CD4 | 0.096 | 0.517 | 29 |
| CD3 | CD3D | 0.065 | 0.510 | 29 |
| CD20 | MS4A1 | 0.058 | 0.472 | 29 |
| CD11c | ITGAX | 0.036 | 0.432 | 29 |
| PDL1 | CD274 | 0.025 | 0.483 | 29 |
| PD1 | PDCD1 | 0.012 | 0.402 | 29 |
| FOXP3 | FOXP3 | 0.011 | 0.390 | 29 |
| PanCK | EPCAM | 0.003 | 0.370 | 17 |
| KI67 | MKI67 | 0.001 | 0.363 | 28 |
| CD45 | PTPRC | NaN | NaN | 0 |

#### Zero-Variance Cases

| Protein | Gene | N Cases | Reason |
|---------|------|---------|--------|
| PanCK | EPCAM | 12 | Zero EPCAM expression in samples |
| CD45 | PTPRC | 29 | Zero PTPRC expression |
| KI67 | MKI67 | 1 | Zero MKI67 in D03 |

**Interpretation:** Protein and RNA show poor cell-level correlation (r ≈ 0.1), with moderate spatial pattern similarity (SSIM ≈ 0.5). Best concordance seen for structural markers (aSMA, CD31) and MHC-II (HLA-DR).

---

## 5. Numerical Results Summary

### All P-values and Statistics

| Analysis | Test | Statistic | P-value | Significant? |
|----------|------|-----------|---------|--------------|
| PCA stage separation | Kruskal-Wallis | H = 1.97 | 0.578 | ❌ No |
| PCA progression trend | Spearman | ρ = -0.082 | 0.673 | ❌ No |
| DE: Normal vs Metaplasia | Mann-Whitney (341 genes) | - | All q > 0.1 | ❌ No |
| DE: Metaplasia vs Cancer | Mann-Whitney (341 genes) | - | All q > 0.1 | ❌ No |
| DE: Normal vs Cancer | Mann-Whitney (341 genes) | - | All q > 0.1 | ❌ No |
| CAF subtype by stage | Kruskal-Wallis | - | > 0.05 | ❌ No |
| CD8 exhaustion by stage | Kruskal-Wallis | - | > 0.05 | ❌ No |
| M1/M2 ratio by stage | Kruskal-Wallis | - | > 0.05 | ❌ No |

### Effect Sizes (Cohen's d for N vs C)

| Comparison | Estimated d | Power (N=29) |
|------------|-------------|--------------|
| Typical DE gene | < 0.3 (small) | < 15% |
| Required for 80% power | d > 1.1 (large) | N/A |

---

## 6. Statistical Significance Assessment

### Why Nothing is Significant

#### 1. Sample Size Limitations

| Analysis Type | N per Group | Power for d=0.8 | Power for d=0.5 |
|---------------|-------------|-----------------|-----------------|
| N vs M | 8 vs 7 | ~45% | ~20% |
| M vs C | 7 vs 10 | ~50% | ~25% |
| N vs C | 8 vs 10 | ~55% | ~25% |

**Conclusion:** With 7-10 samples per group, we only have sufficient power to detect very large effects (d > 1.0).

#### 2. Pseudobulk Aggregation

- Aggregating 1.83M cells to 29 samples loses cell-level heterogeneity
- Stage effects may exist in specific cell subpopulations but be masked at bulk level
- Appropriate for avoiding pseudoreplication, but reduces sensitivity

#### 3. Biological Complexity

- Gastric cancer progression is heterogeneous
- Same "stage" may contain molecularly distinct subtypes
- Spatial organization changes may precede transcriptional changes

### What This Does NOT Mean

- ❌ "No biological differences exist between stages"
- ❌ "The data is low quality"
- ❌ "The analysis is wrong"

### What This DOES Mean

- ✓ No **large, consistent** transcriptional shifts detectable at pseudobulk level
- ✓ Cell-level and spatial analyses may reveal patterns masked by aggregation
- ✓ More patients per stage needed for powered pseudobulk DE

---

## 7. Limitations & Caveats

### Data Limitations

| Limitation | Impact | Mitigation |
|------------|--------|------------|
| Only 4 patients total | Limited biological replication | Within-patient comparisons |
| 341 genes (not transcriptome) | May miss key pathways | Panel was designed for immune/stromal |
| SNU-484 has only cancer | Cannot assess progression in this patient | Focus on SNU-105, SNU-107 |
| Missing genes (CXCL12, ARG1, keratins) | Incomplete signatures | Document and use alternatives |

### Methodological Caveats

| Choice | Rationale | Alternative |
|--------|-----------|-------------|
| Pseudobulk (not cell-level) | Avoid pseudoreplication | Mixed models, permutation |
| Mann-Whitney (not limma) | Non-parametric, small N | limma-voom with patient blocking |
| BH FDR correction | Standard, conservative | Local FDR, q-value |

### Confounding Structure

While the original plan suggested complete stage-patient confounding, the actual data shows:

- **SNU-105:** Has Normal + Metaplasia + Cancer (can compare within patient)
- **SNU-107:** Has Normal + Metaplasia + Cancer (can compare within patient)
- **SNU-484:** Has Cancer only
- **ctrl:** Has Control only

This enables within-patient comparisons that eliminate patient effects.

---

## 8. Recommended Next Steps

### Immediate Priority: Within-Patient Comparisons

**Rationale:** SNU-105 and SNU-107 each have cells from all 3 disease stages. Comparing stages within the same patient eliminates patient effects entirely.

```
Script 49: Within-Patient Progression Analysis
- For SNU-105: Compare Normal → Metaplasia → Cancer
- For SNU-107: Compare Normal → Metaplasia → Cancer
- Cell type proportions
- Gene expression shifts
- Spatial organization
```

### High Priority: Spatial Statistics

**Rationale:** Spatial analyses don't require pseudobulk aggregation and can detect reorganization patterns.

```
Script 50: Spatial Statistics
- Neighborhood enrichment (who neighbors whom?)
- Ripley's K function (clustering patterns)
- Co-occurrence analysis
- Distance to nearest cell type
```

### Medium Priority: Cell-Level Trajectories

```
Script 51: Trajectory Analysis
- CellRank pseudotime
- RNA velocity (if splicing info available)
- Diffusion maps within epithelial cells
```

### Medium Priority: Niche Analysis

```
Script 52: Microenvironment Niches
- Niche composition around CD8 T cells
- CAF neighborhood analysis
- Immune exclusion patterns
```

### Lower Priority: Ligand-Receptor Analysis

```
Script 53: Cell-Cell Communication
- LIANA+ analysis per stage
- CellChat comparison
- Spatial ligand-receptor scoring
```

### For Future Studies

| Recommendation | Rationale |
|----------------|-----------|
| N ≥ 5 patients per stage | 80% power for medium effects |
| Matched longitudinal samples | Same patient over time |
| Whole transcriptome | Not limited to panel |
| Include CXCL12, ARG1, keratins | Complete CAF/macrophage/epithelial signatures |

---

## 9. File Manifest

### Scripts

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `60_load_all_samples.py` | Load raw data | Raw TIFF/CSV | h5ad files |
| `61_comprehensive_qc.py` | Sample QC | Raw h5ad | QC metrics |
| `62_process_qc_passing.py` | Cell QC + WNN | QC'd h5ad | Processed h5ad |
| `63_merge_and_batch_correct.py` | Merge + Harmony | Processed h5ad | merged_corrected.h5ad |
| `46_pca_deep_dive.py` | PCA analysis | merged_corrected | PCA results |
| `47_celltype_specific_de.py` | Cell type DE | merged_corrected | DE results |
| `48_protein_rna_spatial_correlation.py` | Protein-RNA | merged_corrected | Correlation results |

### Results Directories

```
results/
├── qc_all_samples/
│   ├── raw/                     # 32 raw h5ad files
│   ├── final_processed/         # 29 processed h5ad files
│   └── merged/
│       ├── merged_counts.h5ad   # Raw counts (backup)
│       ├── merged_normalized.h5ad
│       └── merged_corrected.h5ad  # Main analysis file
│
├── pca_deep_dive/
│   ├── PCA_DEEP_DIVE_REPORT.md
│   ├── figures/                 # 9 figures
│   └── data/                    # Variance decomposition, loadings
│
├── celltype_de/
│   ├── CELLTYPE_DE_REPORT.md
│   ├── epithelial/              # DE results (3 comparisons)
│   ├── caf/                     # CAF subtype scores
│   ├── cd8/                     # Exhaustion scores
│   ├── macrophage/              # M1/M2 scores
│   ├── sensitivity/             # Admixture analysis
│   └── figures/                 # 4 figures
│
├── protein_rna_correlation/
│   ├── PROTEIN_RNA_REPORT.md
│   ├── correlation_details.csv  # All 435 pairs
│   ├── correlation_matrix.csv   # Samples × Proteins
│   ├── ssim_matrix.csv
│   ├── zero_variance_cases.csv
│   └── figures/                 # Summary figures
│
└── RESULTS_SUMMARY.md           # Previous summary
```

---

## 10. Reproducibility

### Environment

```bash
conda activate enact
# Python 3.9, scanpy 1.10.3, squidpy 1.6.1
```

### Execution Commands

```bash
cd ~/g4x-choi-batch2-analysis

# QC Pipeline
python scripts/60_load_all_samples.py 2>&1 | tee logs/60_loading.log
python scripts/61_comprehensive_qc.py 2>&1 | tee logs/61_qc.log
python scripts/62_process_qc_passing.py --parallel 8 2>&1 | tee logs/62_processing.log
python scripts/63_merge_and_batch_correct.py 2>&1 | tee logs/63_merge.log

# Analysis Scripts
python scripts/46_pca_deep_dive.py 2>&1 | tee logs/46_pca.log
python scripts/47_celltype_specific_de.py 2>&1 | tee logs/47_de.log
python scripts/48_protein_rna_spatial_correlation.py 2>&1 | tee logs/48_corr.log
```

### Verification

```bash
# Check cell counts
python -c "import scanpy as sc; a=sc.read_h5ad('results/qc_all_samples/merged/merged_corrected.h5ad'); print(f'{a.n_obs:,} cells, {a.n_vars} genes')"
# Expected: 1,835,026 cells, 341 genes

# Check DE results
python -c "import pandas as pd; df=pd.read_csv('results/celltype_de/epithelial/de_N_vs_C.csv'); print(f'Significant: {(df.qval_bh<0.05).sum()}')"
# Expected: Significant: 0
```

### Git History

```bash
git log --oneline -10
# 70141c0 feat: Add cell type DE and protein-RNA correlation analyses
# ...
```

---

## Appendix A: Power Analysis

For two-group comparison with Mann-Whitney U test:

| Effect Size (d) | N per group for 80% power |
|-----------------|---------------------------|
| 0.2 (small) | 310 |
| 0.5 (medium) | 51 |
| 0.8 (large) | 21 |
| 1.0 (very large) | 14 |
| 1.2 | 10 |

**Current study:** N = 7-10 per stage → Can only detect d > 1.0 with 80% power.

---

## Appendix B: Gene Panel Coverage

### Genes in Panel (341 total)

Key pathway coverage:
- Immune checkpoint: PDCD1, CD274, LAG3, HAVCR2, TIGIT, CTLA4 ✓
- T cell: CD4, CD8A, CD3D, CD3E, FOXP3 ✓
- Macrophage: CD68, CD163, MRC1, CD80, TNF, IL1B ✓ (missing ARG1)
- CAF: ACTA2, TAGLN, FAP, PDGFRA, IL6, CD74, HLA-DRA ✓ (missing CXCL12)
- Epithelial: EPCAM, CDH1, MUC1, MUC2, MUC5AC ✓ (missing KRT8/18/19)
- Gastric markers: CDX2, TFF1, TFF2, TFF3 ✓

### Missing Key Genes

| Gene | Signature | Impact |
|------|-----------|--------|
| CXCL12 | iCAF | Incomplete inflammatory CAF signature |
| ARG1 | M2 macrophage | Incomplete M2 signature |
| KRT8/18/19 | Epithelial | PanCK protein has no direct RNA match |

---

*Report generated: 2026-01-22*
*Analysis by: Claude (AI assistant)*
*Principal Investigator: Max Van Belkum, MD-PhD candidate, Vanderbilt University*
