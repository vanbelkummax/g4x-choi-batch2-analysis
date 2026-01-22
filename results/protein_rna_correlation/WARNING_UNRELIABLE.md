# ⚠️ WARNING: RESULTS NOT TRUSTWORTHY

**Status:** PRELIMINARY - DO NOT CITE

## Issues

1. **Very poor correlation**: Mean Pearson r = 0.095
2. **Protein-RNA discordance** is biologically expected but limits utility
3. **Zero-variance cases**: 12 pairs had no RNA signal (especially EPCAM/PTPRC)
4. **PanCK→EPCAM proxy** is not valid (keratins not in panel)

## What This Means

Protein and RNA spatial patterns show only **moderate overlap (SSIM ~0.5)** and **poor cell-level correlation (r ~0.1)**. This is expected due to:
- Translation lag (protein downstream of RNA)
- Post-transcriptional regulation
- Different detection sensitivities

## Recommendation

Do not use protein-RNA correlation for validation. Focus on:
- Spatial statistics (neighborhood enrichment)
- Within-modality analyses
- Raw counts from `merged_counts.h5ad`

See `/DATA_QUALITY_NOTICE.md` for full details.
