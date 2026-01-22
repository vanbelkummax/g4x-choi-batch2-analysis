# G4X Full Dataset QC Plan (v2 - Enhanced)

**Created:** 2026-01-21
**Status:** READY FOR APPROVAL

---

## Executive Summary

Processing all 32 samples with rigorous, standardized QC metrics based on Nature Biotechnology (2025) framework.

| Current State | Full Dataset | Gap |
|---------------|--------------|-----|
| 8 samples (Lane 1) | 32 samples | 24 samples unprocessed |
| 522K cells | 2.3M cells | 1.7M cells missing |

**Dataset:** 4.5mm G4X Multi-modal, Pre-GC custom panel, 16 protein panel, 28 ROIs, 13 patients

---

## Complete Sample Inventory

| Lane | Position | Patient | Stage | Cells | Notes |
|------|----------|---------|-------|-------|-------|
| **L001** | A01 | WD112575 | CRC Ctrl | 127,792 | ✅ Analyzed |
| | B01 | SNU-103 | N | 30,194 | ✅ Analyzed |
| | C01 | SNU-103 | M | 38,588 | ✅ Analyzed |
| | D01 | SNU-105 | N | 31,711 | ✅ 11.9% empty |
| | E01 | SNU-105 | M | 34,353 | ✅ Analyzed |
| | F01 | SNU-105 | C | 62,875 | ✅ Low depth (34 t/c) |
| | G01 | SNU-109 | C | 158,727 | ✅ Analyzed |
| | H01 | SNU-109 | X | 85,842 | ✅ Analyzed |
| **L002** | A02 | WD112575 | CRC Ctrl | 120,260 | |
| | B02 | SNU-106 | N | 37,511 | |
| | C02 | SNU-106 | C | 52,588 | |
| | D02 | SNU-106 | N | 56,035 | |
| | E02 | SNU-107 | N | 31,790 | |
| | F02 | SNU-107 | M | 33,857 | |
| | G02 | SNU-107 | C | 79,634 | |
| | H02 | SNU-107 | C | 87,215 | High depth (109 t/c) |
| **L003** | A03 | WD112575 | CRC Ctrl | 145,070 | |
| | B03 | SNU-108 | M | 38,095 | |
| | C03 | SNU-108 | C | 21,231 | Lowest cells |
| | D03 | SNU-108 | M | 33,227 | |
| | E03 | SNU-110 | C | 86,003 | |
| | F03 | SNU-110 | X | 82,704 | |
| | G03 | UTSW-4 | body | 45,583 | |
| | H03 | UTSW-4 | fundus | 35,585 | High depth (126 t/c) |
| **L004** | A04 | WD112575 | CRC Ctrl | 130,891 | Different patients |
| | B04 | Hu13 | C | 119,211 | Smaller cell areas |
| | C04 | Hu13 | M | 89,387 | (54-70 µm²) |
| | D04 | Hu13 | C | 169,683 | vs 85-125 µm² |
| | E04 | Hu15 | X | 68,014 | |
| | F04 | Hu15 | X | 78,831 | |
| | G04 | Hu15 | X | 64,147 | |
| | **H04** | **Hu15** | **X** | **32,334** | **25 t/c, 11 g/c** |

---

## Available Resources

| Resource | Available | Strategy |
|----------|-----------|----------|
| **RAM** | 196 GB | Load all samples simultaneously |
| **VRAM** | 23.5 GB | GPU-accelerated UMAP, neighbors |
| **CPUs** | 20 cores | Parallel QC computation |

---

## Standardized QC Framework

Based on Nature Biotechnology (2025) standardized metrics and protocols.io CosMx guidelines.

### Tier 1: Segmentation/Background Health

| Metric | Target | Warning | Fail |
|--------|--------|---------|------|
| **FTC** (Fraction Transcripts in Cells) | ≥0.80 | 0.60-0.80 | <0.60 |
| **Hard empty** (0 transcripts) | ~0% | >1% | >5% |
| **Effective empty** (<T_min) | <30% | 30-40% | >40% |

### Tier 2: Assay Sensitivity / Tissue Quality

| Metric | Target | Warning | Fail |
|--------|--------|---------|------|
| **Median transcripts/cell** | ≥50 | 30-50 | <30 |
| **Median genes/cell** | ≥20 | 15-20 | <15 |
| **10th percentile trans/cell** | ≥20 | 10-20 | <10 |

### Tier 3: Negative Probe / Noise Gate

| Metric | Target | Warning | Fail |
|--------|--------|---------|------|
| **Neg probe fraction per cell** | <0.05 | 0.05-0.10 | >0.10 |
| **Mean neg counts/cell** | Low | - | High |

### Tier 4: Segmentation Specificity (MECR)

| Metric | Target | Warning | Fail |
|--------|--------|---------|------|
| **MECR** (Mutually Exclusive Co-expression) | <0.05 | 0.05-0.10 | >0.10 |

**MECR pairs to test:**
- EPCAM vs PTPRC (epithelial vs immune)
- PanCK vs CD45 (epithelial vs immune - protein)
- CD3 vs CD68 (T-cell vs macrophage)
- COL1A1 vs EPCAM (stromal vs epithelial)

### Tier 5: FOV-Level QC

| Metric | Target | Warning | Fail |
|--------|--------|---------|------|
| **Mean counts/cell per FOV** | ≥100 | 50-100 | <50 |
| **Cell density consistency** | Stable | Variable | Patchy deserts |

---

## Multi-Threshold Sensitivity Testing

```python
T_MIN_VALUES = [20, 50, 100]  # Test all three

for t_min in T_MIN_VALUES:
    # Filter cells
    adata_filtered = adata[adata.obs['total_counts'] >= t_min]

    # Track: cell retention, cluster stability, marker expression
    results[t_min] = {
        'n_cells': len(adata_filtered),
        'pct_retained': len(adata_filtered) / len(adata) * 100,
        'n_clusters': len(adata_filtered.obs['leiden'].unique()),
        # ... marker retention
    }
```

**Goal:** Verify biology is stable across thresholds, not fragile to QC.

---

## Implementation: 3 Scripts

### Script 60: `60_load_all_samples.py`

**Purpose:** Load all 32 samples, compute standardized QC metrics

**Key outputs per sample:**
```python
qc_metrics = {
    'sample_id': str,
    'n_cells_total': int,
    'median_transcripts_per_cell': float,
    'p10_transcripts_per_cell': float,
    'p90_transcripts_per_cell': float,
    'median_genes_per_cell': float,
    'pct_hard_empty': float,           # 0 transcripts
    'pct_effective_empty_20': float,   # <20 transcripts
    'pct_effective_empty_50': float,   # <50 transcripts
    'pct_effective_empty_100': float,  # <100 transcripts
    'ftc': float,                      # Fraction transcripts in cells
    'median_cell_area_um2': float,
    'mecr_epcam_ptprc': float,         # If markers available
    # Negative probe metrics (if available)
    'mean_neg_counts_per_cell': float,
    'pct_cells_neg_fraction_gt_10': float,
}
```

**Output:**
- `results/qc_all_samples/raw/{sample}_raw.h5ad` (32 files)
- `results/qc_all_samples/qc_metrics_all_samples.csv`
- `results/qc_all_samples/loading_summary.csv`

---

### Script 61: `61_comprehensive_qc.py`

**Purpose:** Apply tiered QC, generate visualizations, assess batch effects

#### Phase 1: Build QC Scorecard

Apply tiered framework to classify each sample:
- **GREEN:** Passes all Tier 1-2 metrics
- **YELLOW:** Warning on any metric
- **RED:** Fails any metric

#### Phase 2: Sensitivity Analysis

Test T_min = 20, 50, 100 on all samples:
- Track cell retention curves
- Verify cluster stability
- Check marker retention

#### Phase 3: Visual QC (per sample, 8-panel figure)

1. Spatial scatter: total_counts heatmap
2. Spatial scatter: n_genes heatmap
3. Violin: transcript distribution (with T_min lines)
4. Scatter: counts vs genes (QC cloud)
5. Histogram: cell area distribution
6. Protein SNR bar chart
7. MECR heatmap (marker co-expression)
8. FOV-level counts/cell

#### Phase 4: Cross-Sample Comparisons

- Ridge plot: transcript distributions (all 32)
- Pseudobulk PCA: color by lane, shape by stage
- Sample correlation heatmap
- FTC vs median_transcripts scatter

#### Phase 5: Lane 4 Diagnosis

Three-step discipline for cell area difference:

1. **Diagnose:** Is it technical or compositional?
   - Compare cell type proportions Lane 4 vs others
   - Check segmentation parameters
   - Compare same patient CRC controls (A01-A04)

2. **Fix upstream if needed:**
   - If technical: resegmentation with harmonized parameters
   - If biological: document and proceed

3. **Model downstream:**
   - Include lane as batch covariate
   - Run sensitivity: with/without Lane 4

**Output:**
```
results/qc_all_samples/
├── figures/
│   ├── per_sample/           # 32 × 8-panel QC figures
│   ├── cross_sample/         # Comparison plots
│   ├── sensitivity/          # T_min sensitivity curves
│   └── lane4_diagnosis/      # Batch effect investigation
├── qc_scorecard.csv          # GREEN/YELLOW/RED per sample
├── sensitivity_results.csv   # T_min = 20/50/100 comparison
├── lane4_diagnosis.md        # Technical vs biological assessment
└── QC_REPORT.md              # Comprehensive report
```

---

### Script 62: `62_process_qc_passing.py`

**Purpose:** Run full pipeline on GREEN samples

**Pipeline per sample:**
1. Load raw h5ad
2. Apply chosen T_min filter
3. WNN integration
4. Admixture scoring
5. Segmentation QC
6. Hierarchical annotation
7. Save final h5ad

**Output:**
```
results/qc_all_samples/
├── final_processed/
│   ├── {sample}_final.h5ad   # GREEN samples only
│   └── ...
└── FINAL_REPORT.md
```

---

## Decision Tree for Sample Inclusion

```
┌─────────────────────────────────────────────────────┐
│                  For each sample                     │
└─────────────────────────────────────────────────────┘
                         │
                         ▼
            ┌─────────────────────────┐
            │   FTC ≥ 0.80?           │
            └─────────────────────────┘
                    │           │
                   YES         NO
                    │           │
                    ▼           ▼
    ┌────────────────────┐  ┌────────────────────┐
    │ Median trans ≥ 50? │  │ FTC 0.60-0.80?     │
    └────────────────────┘  └────────────────────┘
           │        │              │        │
          YES      NO             YES      NO
           │        │              │        │
           ▼        ▼              ▼        ▼
      ┌────────┐ ┌────────┐  ┌────────┐ ┌────────┐
      │ GREEN  │ │ YELLOW │  │ YELLOW │ │  RED   │
      │ (Pass) │ │ (Warn) │  │ (Warn) │ │ (Fail) │
      └────────┘ └────────┘  └────────┘ └────────┘
```

---

## Pre-Identified Issues (Updated)

| Sample | Metric | Value | Assessment |
|--------|--------|-------|------------|
| **H04** | trans/cell | 25 | **RED** - Likely fail |
| **H04** | genes/cell | 11 | **RED** - Likely fail |
| D01 | pct_empty | 11.9% | YELLOW - Investigate FTC |
| F01 | trans/cell | 34 | YELLOW - Below 50 |
| C03 | n_cells | 21K | YELLOW - Borderline |
| Lane 4 | cell_area | 54-70 µm² | **DIAGNOSE** - Tech vs Bio |

---

## Expected Outcomes

| Classification | Expected Count |
|----------------|----------------|
| GREEN (Pass) | ~26-28 |
| YELLOW (Warning) | ~3-4 |
| RED (Fail) | ~1-2 (H04 definite) |

### Predicted Final Cell Counts

| T_min | Cells Retained | % Retained |
|-------|----------------|------------|
| 20 | ~2.1M | ~95% |
| 50 | ~1.9M | ~85% |
| 100 | ~1.6M | ~70% |

---

## Execution Timeline

| Step | Script | Duration | Memory |
|------|--------|----------|--------|
| 1 | Create directories | 1 min | - |
| 2 | `60_load_all_samples.py` | ~15 min | ~25 GB |
| 3 | `61_comprehensive_qc.py` | ~30 min | ~50 GB |
| 4 | Review QC_REPORT.md | Manual | - |
| 5 | `62_process_qc_passing.py` | ~90 min | ~60 GB |
| **Total** | | **~2.5 hours** | **60 GB peak** |

---

## Verification Checklist

### Loading (Script 60)
- [ ] 32 raw h5ad files created
- [ ] qc_metrics_all_samples.csv has 32 rows
- [ ] All samples have spatial coordinates
- [ ] Protein data loaded (17 markers)

### QC (Script 61)
- [ ] H04 classified as RED
- [ ] Sensitivity curves generated for T_min = 20/50/100
- [ ] Lane 4 diagnosis completed
- [ ] MECR computed (if markers available)
- [ ] 32 QC panel figures generated
- [ ] QC_REPORT.md is comprehensive

### Pipeline (Script 62)
- [ ] Only GREEN samples processed
- [ ] WNN integration successful
- [ ] Cell type annotations present
- [ ] FINAL_REPORT.md generated

---

## Files to Create

| Script | Purpose |
|--------|---------|
| `scripts/60_load_all_samples.py` | Parallel loading + QC metrics |
| `scripts/61_comprehensive_qc.py` | Tiered QC + visualizations |
| `scripts/62_process_qc_passing.py` | Pipeline on passing samples |

## Files to Reuse

| Existing Script | What to Extract |
|-----------------|-----------------|
| `20_g4x_choi_batch2_analysis.py` | `load_sample()` |
| `30_wnn_integration.py` | `compute_weighted_integration()` |
| `31_segmentation_qc.py` | `compute_qc_scores()` |
| `32_hierarchical_annotation.py` | annotation logic |

---

## QC Scorecard Template

```
# G4X Choi Batch 2 - QC Scorecard

## Summary
- Total samples: 32
- GREEN (Pass): XX
- YELLOW (Warning): XX
- RED (Fail): XX

## Detailed Results

| Sample | FTC | Med Trans | Med Genes | Empty% | MECR | Area | Status |
|--------|-----|-----------|-----------|--------|------|------|--------|
| A01 | 0.94 | 42 | 28 | 0.8% | 0.03 | 85 | GREEN |
| ... | ... | ... | ... | ... | ... | ... | ... |
| H04 | 0.94 | 25 | 11 | 1.4% | ? | 69 | RED |

## Recommendations
1. Exclude H04 (fails trans + genes thresholds)
2. Include D01 with warning (high empty%, check FTC)
3. Use T_min = 50 as default (sensitivity stable)
4. Include lane as covariate (cell area difference in L4)
```

---

## Notes

- **Environment:** `conda activate enact`
- **FTC calculation:** May need to compute from raw data if not in core_metrics.csv
- **MECR:** Requires marker pair selection appropriate for this panel
- **Negative probes:** Check if available in G4X data format
