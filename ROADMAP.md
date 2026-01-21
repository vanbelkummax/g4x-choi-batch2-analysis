# G4X Gastric Cancer Analysis: Project Roadmap

## Project Status

**Completed:**
- ✅ WNN Integration (8/8 samples)
- ✅ Segmentation QC with admixture scoring (8/8 samples)
- ✅ Hierarchical cell type annotation (8/8 samples)
- ✅ Validation vs Resolve (8.4% QC filtering, 1.0 correlation)

**Key Finding:** Epithelial cells decrease with progression (39% → 22%)

---

## Option A: Finish the Science

**Goal:** Characterize gastric cancer progression using spatial transcriptomics

**Timeline:** 1-2 days

**Output:** Poster/short paper on N→M→C progression

### A1. Progression Analysis Script
**Script:** `34_progression_analysis.py`

**Samples:** SNU-105 series (matched patient)
- D01: Normal gastric mucosa
- E01: Intestinal metaplasia
- F01: Gastric cancer

**Analyses:**
```
┌─────────────────────────────────────────────────────────────────┐
│ A1.1 Cell Type Proportion Changes                               │
│ ─────────────────────────────────────────────────────────────── │
│ - Bar plot: lineage % across N→M→C                              │
│ - Statistical test: chi-square or Fisher's exact                │
│ - Key question: Which cell types expand/contract?               │
├─────────────────────────────────────────────────────────────────┤
│ A1.2 Pseudobulk Differential Expression                         │
│ ─────────────────────────────────────────────────────────────── │
│ - Normal vs Metaplasia                                          │
│ - Metaplasia vs Cancer                                          │
│ - Normal vs Cancer                                              │
│ - Method: DESeq2 or edgeR on aggregated counts                  │
├─────────────────────────────────────────────────────────────────┤
│ A1.3 Spatial Neighborhood Analysis                              │
│ ─────────────────────────────────────────────────────────────── │
│ - Squidpy nhood_enrichment per stage                            │
│ - Question: Do immune cells get excluded from tumor?            │
│ - Co-occurrence analysis: which cells are neighbors?            │
├─────────────────────────────────────────────────────────────────┤
│ A1.4 Marker Gene Expression                                     │
│ ─────────────────────────────────────────────────────────────── │
│ - Known progression markers (CDX2, MUC2, MUC5AC)                │
│ - Spatial expression patterns                                   │
│ - Violin plots per stage                                        │
└─────────────────────────────────────────────────────────────────┘
```

### A2. Key Figures to Generate

| Figure | Content | Purpose |
|--------|---------|---------|
| Fig 1 | Spatial maps colored by cell type (N, M, C) | Overview |
| Fig 2 | Cell type proportion stacked bar | Composition changes |
| Fig 3 | Volcano plots (N vs M, M vs C) | DE genes |
| Fig 4 | Neighborhood enrichment heatmaps | Spatial organization |
| Fig 5 | Marker gene spatial expression | Validation |

### A3. Expected Findings

Based on preliminary data:
1. **Epithelial decrease:** 25% (Normal) → 27% (Metaplasia) → 22% (Cancer)
2. **Stromal increase:** Possible desmoplastic reaction (see G01: 37% fibroblasts)
3. **Immune changes:** TBD - may see exclusion in cancer

### A4. Deliverables

- [ ] `34_progression_analysis.py` - Main analysis script
- [ ] `results/g4x_choi_batch2/progression/` - Output directory
- [ ] Figures for poster/paper
- [ ] Summary statistics CSV

---

## Option B: Validate the Method

**Goal:** Benchmark our pipeline against alternatives and validate novel components

**Timeline:** 2-3 days

**Output:** Methods validation, potential methods paper

### B1. Clustering Comparison
**Script:** `36_clustering_benchmark.py`

**Methods to Compare:**
| Method | Implementation | Notes |
|--------|----------------|-------|
| WNN | Seurat/muon | Our current |
| RNA-only Leiden | Scanpy | Resolve baseline |
| SpatialGlue | spatialglue | Spatial-aware GNN |
| totalVI | scvi-tools | VAE-based |

**Metrics:**
- ARI (Adjusted Rand Index) vs reference
- NMI (Normalized Mutual Information)
- Silhouette score
- Marker coherence per cluster

### B2. Admixture QC Validation
**Script:** `37_admixture_validation.py`

**Questions to Answer:**
1. Are high-admixture cells at segmentation boundaries?
2. Do they have mixed marker expression?
3. Does removing them improve downstream DE?

**Validation Approach:**
```
┌─────────────────────────────────────────────────────────────────┐
│ B2.1 Spatial Distribution                                       │
│ ─────────────────────────────────────────────────────────────── │
│ - Plot admixture score spatially                                │
│ - Check: Are flagged cells at tissue edges/boundaries?          │
├─────────────────────────────────────────────────────────────────┤
│ B2.2 Marker Expression                                          │
│ ─────────────────────────────────────────────────────────────── │
│ - Compare marker profiles: clean vs flagged cells               │
│ - Flagged cells should have mixed lineage markers               │
├─────────────────────────────────────────────────────────────────┤
│ B2.3 Downstream Impact                                          │
│ ─────────────────────────────────────────────────────────────── │
│ - Run DE with/without QC filtering                              │
│ - Compare: Does QC improve signal-to-noise?                     │
└─────────────────────────────────────────────────────────────────┘
```

### B3. Annotation Validation
**Script:** `38_annotation_validation.py`

**Approach:**
1. Check marker gene expression per annotated cell type
2. Compare to published gastric cancer scRNA-seq atlases
3. Calculate annotation confidence metrics

### B4. Method Comparison Summary

**Generate comparison table:**
```
┌──────────────────┬─────────┬───────────┬──────────┬───────────┐
│ Method           │ ARI     │ Silhouette│ Runtime  │ Spatial   │
├──────────────────┼─────────┼───────────┼──────────┼───────────┤
│ RNA-only         │ ?       │ ?         │ Fast     │ No        │
│ WNN (ours)       │ ?       │ ?         │ Medium   │ No        │
│ SpatialGlue      │ ?       │ ?         │ Slow     │ Yes       │
│ totalVI          │ ?       │ ?         │ Medium   │ No        │
└──────────────────┴─────────┴───────────┴──────────┴───────────┘
```

### B5. Deliverables

- [ ] `36_clustering_benchmark.py`
- [ ] `37_admixture_validation.py`
- [ ] `38_annotation_validation.py`
- [ ] `results/g4x_choi_batch2/benchmarks/` - Comparison results
- [ ] METHOD_COMPARISON.md - Updated with results

---

## Recommended Execution Order

```
Week 1: Option A (Science)
├── Day 1: 34_progression_analysis.py
├── Day 2: Generate figures, interpret results
└── Day 3: Write up findings

Week 2: Option B (Methods)
├── Day 1: 36_clustering_benchmark.py
├── Day 2: 37_admixture_validation.py
├── Day 3: 38_annotation_validation.py
└── Day 4: Compile results, update documentation
```

**Rationale:** Science first validates method implicitly. If progression findings are biologically meaningful, the method is working.

---

## Dependencies

### For Option A
```bash
conda activate enact
# Already installed: scanpy, squidpy, muon
```

### For Option B (Additional)
```bash
# SpatialGlue
pip install SpatialGlue

# scvi-tools (for totalVI)
pip install scvi-tools
```

---

## Success Criteria

### Option A Success
- [ ] Clear N→M→C cell type changes documented
- [ ] At least 50 significant DE genes per comparison
- [ ] Spatial neighborhood differences between stages
- [ ] One compelling summary figure

### Option B Success
- [ ] WNN performs ≥ as well as alternatives (or identify better method)
- [ ] Admixture QC removes demonstrably bad cells
- [ ] Annotations match expected marker expression
- [ ] Clear recommendation for pipeline improvements

---

## Notes

- **Data location:** `/mnt/x/Choi_Batch_2_Tuesday/`
- **Results:** `~/spatial-hackathon-2026/results/g4x_choi_batch2/`
- **GitHub:** https://github.com/vanbelkummax/g4x-choi-batch2-analysis

See also:
- `VALIDATION.md` - Resolve comparison results
- `METHOD_COMPARISON.md` - Competitive landscape analysis
- `G4X_ANALYSIS_PLAN.md` - Original analysis plan
