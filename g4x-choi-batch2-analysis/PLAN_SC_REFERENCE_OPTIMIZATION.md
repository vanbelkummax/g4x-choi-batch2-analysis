# Plan: Single-Cell Reference-Based Annotation & Resolution Optimization

**Status:** IN PROGRESS (2026-01-22)
**Branch:** `pilot-clean`

## Current Progress

### Completed
- [x] Downloaded GSE183904_RAW.tar (329 MB, 30 samples)
- [x] Created `06_reference_annotation.py` - KNN label transfer
- [x] Created `07_resolution_optimization.py` - Multi-res clustering
- [x] Pushed scripts to GitHub

### In Progress
- [ ] Building reference from 8 Normal gastric samples (~30-40 min)
- [ ] Annotating G4X samples via KNN

### Pending
- [ ] Resolution optimization
- [ ] Goblet validation

---

## Goal
Use established scRNA-seq reference datasets to:
1. Annotate our G4X cells with higher confidence
2. Validate marker-based findings (Goblet expansion 6%→48%)
3. Find optimal clustering resolution using biology-driven metrics
4. Enable cross-study comparison

---

## Reference Dataset: Kumar et al. GSE183904

**Paper:** "Single-Cell Atlas of Lineage States, Tumor Microenvironment, and Subtype-Specific Expression Programs in Gastric Cancer" (Cancer Discovery 2022)

| Attribute | Value |
|-----------|-------|
| GEO | GSE183904 |
| PMID | 34642171 |
| Samples | 40 total (Normal + Tumor + Peritoneal) |
| Format | Gene x Cell CSV matrices |
| Size | 329 MB |

### Normal Samples Used for Reference (8 samples)

| GSM ID | Sample | Tissue |
|--------|--------|--------|
| GSM5573466 | sample1 | Primary Gastric (Normal) |
| GSM5573469 | sample4 | Primary Gastric (Normal) |
| GSM5573471 | sample6 | Primary Gastric (Normal) |
| GSM5573474 | sample9 | Primary Gastric (Normal) |
| GSM5573476 | sample11 | Primary Gastric (Normal) |
| GSM5573486 | sample21 | Primary Gastric (Normal) |
| GSM5573488 | sample23 | Primary Gastric (Normal) |
| GSM5573490 | sample25 | Primary Gastric (Normal) |

---

## Implementation

### Phase 1: Build Reference (~10-15 min)
```
1. Load 8 Normal gastric samples
2. Concatenate into single reference
3. QC filtering (min_genes=200, min_cells=10)
4. Normalize and log-transform
5. HVG selection (2000 genes)
6. PCA (50 components)
7. Clustering + UMAP
8. Marker-based annotation
9. Save reference h5ad
```

### Phase 2: Label Transfer (~10 min)
```
For each G4X sample:
1. Find common genes (337 G4X ∩ ~20K reference)
2. Compute joint PCA
3. Train KNN classifier on reference (k=15)
4. Predict labels for G4X cells
5. Compute confidence scores
```

### Phase 3: Resolution Optimization (~15 min)
```
Resolutions: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0]

For each resolution:
    For each sample (E02, F02, G02):
        - Cluster at resolution
        - Compute NMI vs reference labels
        - Compute ARI vs reference labels
        - Compute cluster purity
        - Compute silhouette score

    Composite = 0.35*NMI + 0.25*ARI + 0.25*Purity + 0.15*Silhouette

Select resolution with highest mean composite
```

---

## Expected Outputs

| File | Description |
|------|-------------|
| `reference/kumar_gastric_normal_reference.h5ad` | Processed reference |
| `*_reference_annotated.h5ad` | G4X with reference labels |
| `reference_annotations.csv` | Comparison table |
| `resolution_optimization.csv` | Metrics at each resolution |
| `optimal_resolution.txt` | Recommended resolution |
| `figures/reference_umap.png` | Reference UMAP |
| `figures/validation_comparison.png` | Marker vs reference |
| `figures/resolution_optimization.png` | Resolution vs metrics |

---

## Validation Checklist

- [ ] Reference has >10K cells after QC
- [ ] Common genes >200
- [ ] Mean confidence >0.3
- [ ] Goblet trend consistent between methods
- [ ] Optimal resolution identified

---

## Time Estimate

| Step | Estimated |
|------|-----------|
| Load 8 Normal samples | ~5-10 min |
| Preprocess reference | ~5 min |
| Annotate G4X (3 samples) | ~10 min |
| Resolution optimization | ~15 min |
| **Total** | **~30-40 min** |

---

## Scripts

- `scripts/06_reference_annotation.py` - Reference building + KNN transfer
- `scripts/07_resolution_optimization.py` - Multi-res clustering + metrics

---

## Git Log

| Commit | Description |
|--------|-------------|
| `59f29d7a32` | feat: Add reference annotation and resolution optimization scripts |

---

*Last updated: 2026-01-22 15:50*
