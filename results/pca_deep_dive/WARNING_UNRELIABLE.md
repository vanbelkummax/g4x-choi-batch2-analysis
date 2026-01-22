# ⚠️ WARNING: RESULTS NOT TRUSTWORTHY

**Status:** PRELIMINARY - DO NOT CITE

## Issues

1. **No stage separation**: Kruskal-Wallis p = 0.578
2. **No progression trend**: Spearman ρ = -0.08, p = 0.67
3. **Cell type dominates variance** (19%) over stage (1.3%)
4. **Pseudobulk PCA** insensitive to cell-level changes

## Recommendation

Re-analyze using:
- Cell-level PCA with subsampling
- Within-patient trajectory analysis
- Variance decomposition per cell type
- Raw counts from `merged_counts.h5ad`

See `/DATA_QUALITY_NOTICE.md` for full details.
