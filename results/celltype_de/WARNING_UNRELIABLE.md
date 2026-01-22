# ⚠️ WARNING: RESULTS NOT TRUSTWORTHY

**Status:** PRELIMINARY - DO NOT CITE

## Issues

1. **Zero significant genes** across all comparisons (q < 0.05)
2. **Underpowered**: N=7-10 samples per stage insufficient
3. **Pseudobulk approach** masks cell-level heterogeneity
4. **Incomplete signatures**: Missing CXCL12 (iCAF), ARG1 (M2)

## Recommendation

Re-analyze using:
- Within-patient comparisons (SNU-105, SNU-107 have all stages)
- Cell-level mixed models instead of pseudobulk
- Raw counts from `merged_counts.h5ad`

See `/DATA_QUALITY_NOTICE.md` for full details.
