# Multimodal Spatial Integration: Method Comparison

## Overview

This document compares our G4X analysis pipeline against state-of-the-art multimodal integration methods for spatial transcriptomics.

## Method Landscape (2024-2025)

### Benchmark Summary

| Method | Type | ARI Score | Spatial-Aware | RNA+Protein | Notes |
|--------|------|-----------|---------------|-------------|-------|
| **WNN (Seurat)** | Metric learning | 0.71 | ❌ | ✅ | Our current approach |
| **totalVI** | VAE | 0.71 | ❌ | ✅ | Best for CITE-seq |
| **MOFA+** | Factor analysis | 0.63 | ❌ | ✅ | Merges some subtypes |
| **SpatialGlue** | GNN + attention | Higher | ✅ | ✅ | Spatial-specific, better anatomy |
| **SIMO** | Probabilistic | - | ✅ | ✅ | Multi-omics alignment |
| **Multigrate** | VAE | High | ❌ | ✅ | Best for complex combinations |

*ARI = Adjusted Rand Index on CBMC/HBMC benchmarks*

### Key References

- [Multitask benchmarking of single-cell multimodal omics integration methods](https://www.nature.com/articles/s41592-025-02856-3) - Nature Methods 2025
- [Benchmarking single-cell multi-modal data integrations](https://www.nature.com/articles/s41592-025-02737-9) - Nature Methods 2025
- [SpatialGlue: Deciphering spatial domains from spatial multi-omics](https://www.nature.com/articles/s41592-024-02316-4) - Nature Methods 2024
- [SIMO: Spatial integration of multi-omics single-cell data](https://www.nature.com/articles/s41467-025-56523-4) - Nature Communications 2025
- [Orthogonal multimodality integration (OMIC)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05773-y) - BMC Bioinformatics 2024

## Our Pipeline vs Competitors

### What We Use: WNN (Weighted Nearest Neighbors)

**Strengths:**
- Well-established (Seurat v4/v5)
- Variance-based modality weighting
- Good benchmark performance (ARI=0.71)
- Interpretable weights

**Weaknesses:**
- Not spatial-aware (ignores x,y coordinates)
- Computationally expensive at scale
- May merge some cell subtypes

### Potential Upgrades

#### 1. SpatialGlue (Recommended to Evaluate)
- **Why:** Designed specifically for spatial multi-omics
- **Advantage:** Uses spatial coordinates + dual attention mechanism
- **Result:** Better anatomical detail in benchmarks
- **Implementation:** `pip install SpatialGlue`

#### 2. totalVI
- **Why:** Strong benchmark performer for RNA+ADT
- **Advantage:** VAE provides uncertainty quantification
- **Consideration:** Not spatial-aware, similar performance to WNN

#### 3. SIMO
- **Why:** Probabilistic alignment across modalities
- **Advantage:** Handles missing data, extensible to more modalities
- **Status:** Newer, less benchmarked

## Our Unique Contributions

These are **not** in standard benchmarks and represent our value-add:

### 1. Admixture-Based QC (Novel)
```
Standard methods: Filter by counts/genes only
Our approach:    + Admixture score (mixed transcriptomes)
                 + Conflict score (marker disagreement)
                 + Neighborhood consistency
```
**Filters 8.4% additional cells** that pass standard QC but have segmentation errors.

### 2. Hierarchical Annotation
```
Standard: Cluster → Manual annotation
Ours:     Cluster → Lineage (4 types) → Cell type (detailed)
          + Confidence scores at each level
```

### 3. Progression-Matched Analysis
```
Standard: Cross-sectional comparison
Ours:     SNU-105 series (D01→E01→F01)
          Same patient: Normal → Metaplasia → Cancer
```

## Recommendations

### Short-term (This Project)
1. **Keep WNN** - Solid performance, well-validated
2. **Document admixture QC** - This is novel, needs benchmarking
3. **Validate annotations** - Compare to marker expression

### Future Improvements
1. **Try SpatialGlue** - May improve spatial domain detection
2. **Benchmark admixture QC** - Show it removes real artifacts
3. **Compare clustering** - ARI/NMI of WNN vs SpatialGlue vs totalVI

## Benchmark Plan

To properly validate our method:

```
┌─────────────────────────────────────────────────────────────────┐
│ BENCHMARK EXPERIMENT                                            │
├─────────────────────────────────────────────────────────────────┤
│ 1. Run on same data:                                            │
│    - WNN (our current)                                          │
│    - SpatialGlue                                                │
│    - totalVI                                                    │
│    - RNA-only Leiden (Resolve baseline)                         │
│                                                                 │
│ 2. Evaluate:                                                    │
│    - ARI vs known cell types (if available)                     │
│    - Marker expression coherence per cluster                    │
│    - Spatial domain continuity                                  │
│    - Biological interpretability                                │
│                                                                 │
│ 3. Admixture QC validation:                                     │
│    - Compare flagged cells to known doublets                    │
│    - Check if high-admixture cells are at boundaries            │
│    - Downstream impact on DE results                            │
└─────────────────────────────────────────────────────────────────┘
```

## Conclusion

**Current status:** WNN is competitive (tied for best ARI with totalVI), but SpatialGlue may be superior for spatial data.

**Our edge:** Admixture-based QC and hierarchical annotation are novel contributions not covered by existing benchmarks.

**Action:** Proceed with science (Option A) using current pipeline, then benchmark methods (Option B) to validate or upgrade.
