# CRITICAL AUDIT: Raw Counts Preservation Bug

**Date:** 2026-01-22
**Auditor:** Claude (Opus 4.5)
**Severity:** CRITICAL
**Status:** IDENTIFIED - FIX REQUIRED

---

## Executive Summary

A comprehensive technical audit of the G4X Choi Batch 2 analysis pipeline revealed a **critical data integrity bug**: raw counts were lost during per-sample processing (Script 62) and all downstream "counts" files actually contain log-normalized data.

| Component | Status | Notes |
|-----------|--------|-------|
| QC Thresholds | PASS | Appropriate and well-documented |
| Sample Filtering | PASS | H04 excluded, D01 flagged correctly |
| Cell Filtering | PASS | Standard QC metrics applied |
| **Raw Counts Preservation** | **CRITICAL FAIL** | Lost during processing |
| WNN Integration | WARN | Simplified, RNA-biased |
| Batch Correction | WARN | Marginal improvement |
| DE Analysis | WARN | Missing cross-comparison FDR |

---

## The Critical Bug

### Root Cause

**Script 62 (`62_process_qc_passing.py`)** normalizes data WITHOUT preserving raw counts:

```python
# Lines 489-490 in 62_process_qc_passing.py
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# Saves to final_processed/ - RAW COUNTS LOST HERE
```

**Script 63 (`63_merge_and_batch_correct.py`)** then loads already-normalized data:

```python
# Line 500-501: Loads from final_processed/ (already normalized)
adata = sc.read_h5ad(f)  # This is log-normalized, not raw counts

# Line 527: "Saves raw counts" but data is already log-transformed
adata_merged.layers['counts'] = adata_merged.X.copy()  # WRONG - this is normalized!
```

### Evidence

| File | X max | Mean | Integer-like | Actual Content |
|------|-------|------|--------------|----------------|
| `raw/F03_raw.h5ad` | 25 | 0.18 | YES | True raw counts |
| `final_processed/F03_final.h5ad` | 7.96 | 0.68 | NO | Log-normalized |
| `merged_counts.h5ad` | 8.11 | 0.57 | NO | **MISLABELED** - log-normalized |
| `merged_normalized.h5ad` | 8.11 | 0.57 | NO | Log-normalized |
| `merged_corrected.h5ad` | 8.11 | 0.57 | NO | Log-normalized |

All three merged files have **identical** data distributions because raw counts were already lost before merge.

### Impact

1. **scVI batch correction cannot work** - requires raw integer counts
2. **DESeq2/edgeR will fail** - require count matrices
3. **Any count-based method produces garbage** - input is log-normalized
4. **Results are reproducible but based on incorrect assumptions**

### Recovery Path

Raw counts ARE available in `results/qc_all_samples/raw/` (32 files, all integer counts verified).

---

## Additional Issues Found

### 1. WNN Integration is Simplified (MODERATE)

**Location:** Script 62, lines 241-311

**Issue:** Current implementation uses global variance weighting instead of true per-cell WNN:

```python
# Current (simplified):
rna_weight = rna_var / (rna_var + prot_var)  # Global weight
wnn_embedding = concat(rna_pca * rna_weight, prot_pca * prot_weight)

# Should be (true WNN):
# Per-cell modality weights based on neighborhood preservation
# Neighbor graph fusion from both modalities
```

**Impact:** RNA signal dominates (30 PCs vs 15 PCs = 2:1 bias). Multi-modal integration is suboptimal.

**Fix:** Use `muon.pp.wnn()` - function `apply_muon_wnn()` already exists in Script 63 but is not called.

### 2. Batch Correction Marginal (MODERATE)

| Metric | Pre-Correction | Post-Correction | Threshold | Status |
|--------|----------------|-----------------|-----------|--------|
| LISI | 2.46 | 2.66 | > 3.0 | FAIL |
| Silhouette | -0.001 | -0.006 | < 0.3 | PASS (degraded) |

Harmony improved LISI (+0.20) but **degraded** silhouette (-0.005). This suggests possible over-correction.

### 3. Missing Cross-Comparison FDR Correction (MODERATE)

**Location:** Script 47 (DE analysis)

**Issue:** FDR correction applied within each comparison (N vs M, M vs C, N vs C) but NOT across cell types.

**Impact:** If testing 5 cell types x 3 comparisons = 15 tests, family-wise error rate inflates.

**Fix:** Apply Benjamini-Hochberg across all comparisons combined.

### 4. Cell Type Annotation Overly Strict (MINOR)

**Location:** Script 62, lines 429-431

```python
# Current: ALL markers must be above median
all_positive = protein_df[markers].apply(lambda x: x > x.median()).all(axis=1)
```

**Impact:** Misses real cells where one marker is slightly below median.

**Fix:** Use majority voting or signature scoring instead of strict AND logic.

---

## What Was Done Correctly

1. **QC Thresholds:** Evidence-based (min 20K cells, 30 median counts, 20 median genes)
2. **Sample Exclusion:** H04 (failed) correctly excluded, D01 (warned) appropriately flagged
3. **Data Loading:** Proper RNA/protein separation from G4X format
4. **Pseudobulk Approach:** Correctly avoids pseudoreplication for DE
5. **Per-Comparison FDR:** Benjamini-Hochberg properly implemented within comparisons
6. **Disk Monitoring:** 3-tier system (5GB/2GB/1GB) prevents data loss
7. **Resume Capability:** Scripts can restart from checkpoints

---

## Variance Decomposition Results

The PCA variance decomposition reveals important structure:

| PC | Top Factor | R² | Interpretation |
|----|------------|-----|----------------|
| PC1 | cell_type | 19.2% | Biology dominates PC1 |
| PC1 | lane | 5.6% | Acceptable batch effect |
| PC1 | stage | 1.3% | Weak stage signal in PC1 |
| PC4 | patient | 29.5% | Strong patient effect |
| PC4 | stage | 25.5% | **Stage signal in PC4** |

**Key Insight:** Stage separation exists but is captured in **PC4**, not PC1. Previous pseudobulk analysis (p=0.578) only tested PC1.

---

## Recommended Fix Plan

### Phase 1: Critical Fix (MUST DO)

| Step | Script | Description | Time |
|------|--------|-------------|------|
| 1 | `64_fix_raw_counts_preservation.py` | Patch script 62 to preserve raw counts | 10 min |
| 2 | `65_rebuild_merged_from_raw.py` | Rebuild merged files from raw/ directory | 2-3 hrs |
| 3 | `66_validate_raw_counts.py` | Verify raw counts correctly preserved | 5 min |

### Phase 2: Statistical Rigor (SHOULD DO)

| Step | Script | Description | Time |
|------|--------|-------------|------|
| 4 | Patch script 47 | Add cross-comparison FDR correction | 15 min |
| 5 | Enable muon WNN | Call `apply_muon_wnn()` in script 63 | 30 min |

### Phase 3: Refinements (NICE TO HAVE)

| Step | Script | Description | Time |
|------|--------|-------------|------|
| 6 | Patch script 62 | Relax cell type annotation logic | 15 min |
| 7 | Documentation | Update all method descriptions | 20 min |

---

## Verification Checklist

After running fix scripts, verify:

- [ ] `merged_counts.h5ad` X.max() > 50 (raw count range)
- [ ] `merged_counts.h5ad` data is integer-like
- [ ] `merged_counts.h5ad` layers['counts'] matches X
- [ ] `merged_normalized.h5ad` X.max() < 15 (log1p range)
- [ ] `merged_corrected.h5ad` has layers['counts'] with integers
- [ ] scVI batch correction runs without error
- [ ] All 29 passing samples included

---

## Files Included in This Commit

```
CRITICAL_AUDIT_RAW_COUNTS_BUG.md    # This document
scripts/64_fix_raw_counts_preservation.py  # Patches script 62
scripts/65_rebuild_merged_from_raw.py      # Rebuilds merged files correctly
scripts/66_validate_raw_counts.py          # Validates the fix
```

---

## Lessons Learned

1. **Always preserve raw counts in `.layers['counts']` before ANY transformation**
2. **Validate data type (integer vs float) at each pipeline stage**
3. **Add explicit assertions for count-based methods**
4. **Name files accurately - "counts" should mean actual counts**
5. **Audit pipelines end-to-end, not just individual scripts**

---

## Contact

For questions about this audit, refer to the git history or the analysis logs in `logs/`.
