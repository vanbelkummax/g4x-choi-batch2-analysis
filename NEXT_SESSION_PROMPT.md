# G4X Analysis - Next Session Prompt

**Date:** 2026-01-22
**Last Commit:** b2099cf
**Status:** QC Complete, Data Validated, Ready for Within-Patient Analysis

---

## Quick Context (Read This First)

You're working on a **G4X spatial transcriptomics** dataset analyzing **gastric cancer progression** (Normal → Metaplasia → Cancer).

### Current State
- **1.83M cells** across **29 samples** (3 patients + controls)
- **341 RNA genes** + **17 protein markers**
- QC, clustering (23 Leiden clusters), and UMAP complete
- **Data is VALID** for visualization and clustering

### The "Raw Counts Bug" - RESOLVED (Not a Problem)
We discovered that `merged_counts.h5ad` contains normalized data, not raw counts. However:
- **This does NOT affect** UMAP, PCA, clustering, cell type annotation, or Harmony batch correction
- **This would only matter** for scVI or DESeq2/edgeR, which we don't need
- **Raw counts exist** in `results/qc_all_samples/raw/` if ever needed
- **Conclusion:** Proceed with current data - it's valid

---

## Key Discovery: Within-Patient Progression

**Both SNU-105 and SNU-107 have complete Normal → Metaplasia → Cancer from the SAME patient.**

| Patient | Normal | Metaplasia | Cancer | Total |
|---------|--------|------------|--------|-------|
| SNU-105 | 189,590 | 127,313 | 241,149 | 558,052 |
| SNU-107 | 198,875 | 220,549 | 307,858 | 727,282 |

This is extremely valuable because:
- Eliminates patient-to-patient confounding
- Each patient serves as their own control
- Paired analysis has much higher statistical power

---

## What's Been Done

| Task | Status | Location |
|------|--------|----------|
| QC (32 samples) | ✅ Complete | `results/qc_all_samples/` |
| Sample filtering | ✅ 29 pass, 1 fail (H04), 2 excluded | `sample_qc_summary.csv` |
| Clustering | ✅ 23 Leiden clusters | `merged_corrected.h5ad` |
| UMAP | ✅ Complete | `merged_corrected.h5ad` |
| Batch correction | ✅ Harmony | `merged_corrected.h5ad` |
| Cell type annotation | ✅ Protein-based | In obs columns |
| Per-patient UMAPs | ✅ Complete | `results/per_patient_umap/` |
| PCA deep dive | ✅ Complete | `results/pca_deep_dive/` |
| Technical audit | ✅ Complete | `CRITICAL_AUDIT_RAW_COUNTS_BUG.md` |
| Data validity doc | ✅ Complete | `DATA_VALIDITY_STATUS.md` |

---

## What Hasn't Worked

| Analysis | Result | Why |
|----------|--------|-----|
| Pseudobulk DE (all cells) | 0 significant genes | Underpowered (N=7-10/group), cell type heterogeneity masks signal |
| Stage separation in PC1 | p=0.578 | Stage signal is in PC4, not PC1; cell type dominates (19.2% of variance) |
| Protein-RNA correlation | Mean r=0.095 | Expected for imaging platforms; post-transcriptional regulation |

**Key insight:** Global comparisons don't work. Need cell-type-specific and within-patient analyses.

---

## Recommended Next Steps (Priority Order)

### 1. Within-Patient Progression Analysis (HIGHEST PRIORITY)
```bash
# Script to create: scripts/49_within_patient_progression.py
# Compare: SNU-105 Normal → Metaplasia → Cancer
# Compare: SNU-107 Normal → Metaplasia → Cancer
```

What to analyze:
- Cell type proportion changes within each patient
- Cluster shifts during progression
- Spatial reorganization (immune infiltration/exclusion)
- Marker expression changes (CDX2, MUC2, MUC5AC for metaplasia)

### 2. Cell-Type-Specific DE
Instead of pseudobulk across all cells:
- Epithelial cells: Normal vs Cancer
- T cells: Normal vs Cancer (exhaustion markers)
- Fibroblasts: Normal vs Cancer (CAF emergence)

### 3. Spatial Analysis
- Immune exclusion patterns by stage
- Neighborhood composition changes
- TLS (tertiary lymphoid structure) detection

---

## Key File Locations

```
~/g4x-choi-batch2-analysis/
├── results/qc_all_samples/
│   ├── merged/
│   │   ├── merged_corrected.h5ad    # Main analysis file (1.83M cells)
│   │   ├── merged_normalized.h5ad   # Pre-correction
│   │   └── merged_counts.h5ad       # Actually normalized (mislabeled)
│   ├── raw/                         # TRUE raw counts (32 files)
│   └── sample_qc_summary.csv        # QC metrics
├── results/per_patient_umap/        # Patient-specific UMAPs
├── results/pca_deep_dive/           # PCA analysis
├── scripts/                         # Analysis scripts (60-66)
├── DATA_VALIDITY_STATUS.md          # What's valid
├── CRITICAL_AUDIT_RAW_COUNTS_BUG.md # Technical audit
└── TIERED_ANALYSIS_PLAN.md          # Full analysis roadmap
```

---

## Environment

```bash
conda activate enact
cd ~/g4x-choi-batch2-analysis
```

Key packages: scanpy 1.10.3, squidpy 1.6.1, muon 0.1.6

---

## Loading the Data

```python
import scanpy as sc

# Main analysis file
adata = sc.read_h5ad("results/qc_all_samples/merged/merged_corrected.h5ad")
print(f"Loaded: {adata.shape[0]:,} cells, {adata.shape[1]} genes")

# Key columns in adata.obs:
# - patient: SNU-105, SNU-107, SNU-484, ctrl
# - stage: normal, metaplasia, cancer, control
# - leiden: cluster labels (0-22)
# - cell_type: annotated cell types
# - [PROTEIN]_intensity_mean: 17 protein markers
# - cell_x, cell_y: spatial coordinates
```

---

## What NOT to Worry About

1. **Raw counts bug** - Doesn't affect your analyses
2. **Low protein-RNA correlation** - Normal for imaging platforms
3. **No DE genes in pseudobulk** - Expected; need stratified analysis
4. **LISI < 3.0 after Harmony** - Batch effects are modest anyway

---

## Summary for Next Session

**Start here:**
1. Load `merged_corrected.h5ad`
2. Subset to SNU-105 or SNU-107
3. Compare Normal → Metaplasia → Cancer within that patient
4. Look for:
   - Cell type proportion shifts
   - Cluster redistribution
   - Spatial pattern changes
   - Marker expression changes in specific cell types

**The data is valid. The opportunity is the within-patient paired design. Focus on that.**

---

## GitHub

Repository: https://github.com/vanbelkummax/g4x-choi-batch2-analysis
Latest commit: b2099cf (docs: Add data validity status and per-patient UMAP analysis)
