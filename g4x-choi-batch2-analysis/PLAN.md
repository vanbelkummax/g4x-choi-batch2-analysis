# Analysis Plan

## Goal
Optimize RNA clustering and annotate cell types for Normal → Metaplasia → Cancer progression.

## Phase 1: HVG Optimization
- Test 20, 30, 50, ALL genes
- Measure silhouette score
- Pick best per sample

## Phase 2: Cell Type Annotation
- Score cells using gastric marker signatures
- Assign by highest score
- Annotate clusters by majority vote

## Phase 3: Spatial Visualization
- Plot cells at (cell_x, cell_y) colored by type
- Compare progression: E02 → F02 → G02

## Key Markers

**Epithelial:** EPCAM, KRT8, CDH1
**Gastric:** MUC5AC, TFF1, CDX2 (metaplasia)
**Immune:** CD3D (T), CD68 (Mac), CD19 (B)
**Stromal:** COL1A1, FAP (CAF), PECAM1 (Endo)

## Success Criteria
- Silhouette > 0.1
- 5-20 clusters per sample
- Cell types form spatially coherent regions
