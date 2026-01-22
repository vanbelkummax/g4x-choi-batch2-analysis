# ⚠️ DATA QUALITY NOTICE

**Date:** 2026-01-22
**Status:** CRITICAL BUG IDENTIFIED - REANALYSIS REQUIRED

---

## 🚨 CRITICAL BUG: Raw Counts Lost

**The file `merged_counts.h5ad` does NOT contain raw counts.** It contains log-normalized data mislabeled as counts.

See [`CRITICAL_AUDIT_RAW_COUNTS_BUG.md`](CRITICAL_AUDIT_RAW_COUNTS_BUG.md) for full technical details.

| File | Claims to be | Actually contains |
|------|--------------|-------------------|
| `merged_counts.h5ad` | Raw counts | ❌ Log-normalized data |
| `merged_normalized.h5ad` | Normalized | Log-normalized data |
| `merged_corrected.h5ad` | Batch-corrected | Log-normalized data |

**True raw counts are ONLY in:** `results/qc_all_samples/raw/*.h5ad`

---

## Summary

The current analysis results should be considered **PRELIMINARY AND NOT TRUSTWORTHY** for the following reasons:

1. **CRITICAL: Raw counts lost during processing** (see above)
2. **No statistically significant findings** across all analyses
3. **Insufficient power** (N=7-10 samples per stage)
4. **Potential processing artifacts** from normalization/batch correction pipeline
5. **Stage-patient partial confounding** in study design

**Action Required:**
1. Run fix scripts (`64_`, `65_`, `66_`) to rebuild from true raw counts
2. Re-analyze with improved methodology

---

## Results Flagged as Unreliable

### ❌ NOT TRUSTWORTHY - Do Not Cite

| Directory | Analysis | Reason |
|-----------|----------|--------|
| `results/celltype_de/` | Pseudobulk DE | 0/341 significant genes; underpowered |
| `results/celltype_de/epithelial/de_*.csv` | Stage comparisons | All q > 0.1; no signal |
| `results/celltype_de/caf/` | CAF subtyping | No stage trends; incomplete signatures |
| `results/celltype_de/cd8/` | CD8 exhaustion | No significant differences |
| `results/celltype_de/macrophage/` | M1/M2 polarization | No significant differences |
| `results/protein_rna_correlation/` | Protein-RNA concordance | Mean r = 0.095 (very poor) |
| `results/pca_deep_dive/` | PCA stage separation | p = 0.58 (no separation) |

### ⚠️ USE WITH CAUTION

| Directory | Analysis | Caveat |
|-----------|----------|--------|
| `results/qc_all_samples/merged/merged_corrected.h5ad` | Batch-corrected data | Harmony correction applied; may mask biological signal |
| `results/qc_all_samples/merged/merged_normalized.h5ad` | Normalized data | Log-normalized; not suitable for DE |

### ✓ TRUSTWORTHY - Can Be Used

| Directory | Analysis | Notes |
|-----------|----------|-------|
| `results/qc_all_samples/raw/` | **Individual raw h5ad** | ✅ TRUE RAW COUNTS - Use these! |
| `results/qc_all_samples/sample_qc_summary.csv` | QC metrics | Sample-level statistics |

### ❌ MISLABELED - Do Not Trust

| File | Label | Actual Content |
|------|-------|----------------|
| `merged_counts.h5ad` | "Raw counts" | ❌ Log-normalized (MISLABELED) |

---

## What Went Wrong

### 1. Statistical Power
- Only 29 samples total across 4 stages
- N = 7-10 per stage insufficient for pseudobulk DE
- Would need N ≥ 50 per group for medium effect sizes

### 2. Pseudobulk Aggregation
- Collapsed 1.83M cells → 29 data points
- Lost all cell-level heterogeneity
- Stage effects may exist in subpopulations but be masked

### 3. Protein-RNA Discordance
- Mean correlation r = 0.095 indicates protein and RNA don't track
- Biological (translation lag) + technical factors
- Cannot use protein as proxy for RNA or vice versa

### 4. Missing Signature Genes
- iCAF: Missing CXCL12
- M2 macrophage: Missing ARG1
- Epithelial: Missing KRT8/18/19 (PanCK protein has no RNA match)

---

## Recommended Reanalysis Plan

### Start Fresh From Raw Counts

```bash
# Use this file for reanalysis:
results/qc_all_samples/merged/merged_counts.h5ad
```

### Priority Analyses (From Raw Counts)

1. **Within-Patient Comparisons**
   - SNU-105: Normal → Metaplasia → Cancer
   - SNU-107: Normal → Metaplasia → Cancer
   - Eliminates patient confounding entirely

2. **Cell-Level Analysis (Not Pseudobulk)**
   - Mixed models with patient as random effect
   - Permutation-based testing
   - Focus on specific cell populations

3. **Spatial Statistics**
   - Neighborhood enrichment
   - Ripley's K
   - Co-occurrence (don't require DE)

4. **Trajectory Analysis**
   - CellRank within epithelial cells
   - Pseudotime without stage labels

---

## Scripts to Rewrite

| Current Script | Issue | Recommended Fix |
|----------------|-------|-----------------|
| `47_celltype_specific_de.py` | Mann-Whitney underpowered | Use mixed models or within-patient design |
| `48_protein_rna_spatial_correlation.py` | Poor correlation expected | Focus on spatial patterns, not correlation |
| `46_pca_deep_dive.py` | Pseudobulk PCA insensitive | Cell-level PCA with subsampling |

---

## File Status Summary

```
results/
├── qc_all_samples/
│   ├── raw/                        ✅ TRUSTWORTHY (TRUE raw counts)
│   └── merged/
│       ├── merged_counts.h5ad      ❌ MISLABELED (actually log-normalized!)
│       ├── merged_normalized.h5ad  ⚠️ CAUTION (log-normalized)
│       └── merged_corrected.h5ad   ⚠️ CAUTION (batch corrected, log-normalized)
│
├── pca_deep_dive/                  ❌ NOT TRUSTWORTHY
├── celltype_de/                    ❌ NOT TRUSTWORTHY
└── protein_rna_correlation/        ❌ NOT TRUSTWORTHY
```

## Fix Scripts Available

Run these to rebuild from true raw counts:
```bash
python scripts/64_fix_raw_counts_preservation.py
python scripts/65_rebuild_merged_from_raw.py
python scripts/66_validate_raw_counts.py
```

---

## Contact

For questions about data quality or reanalysis plans, contact the analysis team.

**This notice should remain in the repository until reanalysis is complete.**
