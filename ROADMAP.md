# G4X Gastric Cancer Analysis: Project Roadmap

## Project Status

**Completed:**
- ✅ WNN Integration (8/8 samples)
- ✅ Segmentation QC with admixture scoring (8/8 samples)
- ✅ Hierarchical cell type annotation (8/8 samples)
- ✅ Validation vs Resolve (8.4% QC filtering, 1.0 correlation)

**Key Finding:** Epithelial cells decrease with progression (39% → 22%)

---

## Panel Availability (351 RNA genes, 17 proteins)

### Available Markers by Analysis Type

| Analysis | Available Markers | Feasibility |
|----------|-------------------|-------------|
| **CAF Subtyping** | mCAF: ACTA2, TAGLN, aSMA(prot) / iCAF: IL6, PDPN, PDGFRA / apCAF: CD74, HLA-DRA, HLA-DR(prot) | ✅ Full |
| **Cell-Cell Comm** | 23/26 key LR pairs (CCL2, CCL5, CXCL9, CXCL10, IL6, TGFB1, etc.) | ✅ Full |
| **Gastric Markers** | CDX2, MUC2, MUC5AC, MUC6, TFF1, TFF2, TFF3, PGC | ✅ Full |
| **IEX Signature** | Only TGFBI (1/4: DDR1, PAK4, DPEP1 missing) | ⚠️ Limited |
| **Proliferation** | PCNA, TOP2A, ERBB2, TP53, CDKN2A | ✅ Partial |

---

## COMPREHENSIVE ANALYSIS PLAN

### Phase 1: Core Progression Analysis
**Script:** `34_progression_analysis.py`
**Samples:** D01 (Normal) → E01 (Metaplasia) → F01 (Cancer)

```
┌─────────────────────────────────────────────────────────────────┐
│ 1.1 Cell Type Proportion Changes                                │
│ ─────────────────────────────────────────────────────────────── │
│ - Bar plot: lineage % across N→M→C                              │
│ - Statistical test: chi-square or Fisher's exact                │
│ - Key question: Which cell types expand/contract?               │
├─────────────────────────────────────────────────────────────────┤
│ 1.2 Pseudobulk Differential Expression                          │
│ ─────────────────────────────────────────────────────────────── │
│ - Normal vs Metaplasia                                          │
│ - Metaplasia vs Cancer                                          │
│ - Normal vs Cancer                                              │
│ - Method: scanpy rank_genes_groups (Wilcoxon) or DESeq2         │
├─────────────────────────────────────────────────────────────────┤
│ 1.3 Gastric Progression Markers                                 │
│ ─────────────────────────────────────────────────────────────── │
│ - CDX2: intestinal metaplasia marker                            │
│ - MUC2/MUC5AC: goblet cell / gastric foveolar                   │
│ - TFF1/TFF2/TFF3: trefoil factors                               │
│ - Spatial expression patterns + violin plots per stage          │
└─────────────────────────────────────────────────────────────────┘
```

**Deliverables:**
- [ ] Cell type proportion statistics
- [ ] DE gene lists (N vs M, M vs C, N vs C)
- [ ] Marker gene spatial plots

---

### Phase 2: Cell-Cell Communication Analysis (NEW)
**Script:** `35_cellcell_communication.py`
**Method:** LIANA+ (integrates CellPhoneDB, CellChat, NicheNet)

```
┌─────────────────────────────────────────────────────────────────┐
│ 2.1 LIANA+ Analysis                                             │
│ ─────────────────────────────────────────────────────────────── │
│ - Run LIANA rank_aggregate per sample/stage                     │
│ - Resource: 'consensus' database                                │
│ - Group by cell_type annotation                                 │
├─────────────────────────────────────────────────────────────────┤
│ 2.2 Stage Comparison                                            │
│ ─────────────────────────────────────────────────────────────── │
│ - Compare LR interactions: Normal vs Metaplasia vs Cancer       │
│ - Key question: Which interactions change with progression?     │
│ - Focus: Epithelial-Immune, CAF-Immune, CAF-Epithelial          │
├─────────────────────────────────────────────────────────────────┤
│ 2.3 Visualization                                               │
│ ─────────────────────────────────────────────────────────────── │
│ - Dotplot: CAF → Immune interactions                            │
│ - Chord diagram: top interactions per stage                     │
│ - Heatmap: LR pair magnitude changes                            │
└─────────────────────────────────────────────────────────────────┘
```

**Available LR Pairs (23/26):**
CCL2, CCL5, CXCL9, CXCL10, IL6, IL1B, TGFB1, CCR2, CCR5, CXCR3, CXCR4,
CD40-CD40LG, CD80/CD86-CD28/CTLA4, PD1-PDL1, VEGFA/B-KDR

**Deliverables:**
- [ ] LIANA results per stage
- [ ] Top changing interactions list
- [ ] Communication dotplots/chord diagrams

---

### Phase 3: CAF Subtyping (NEW)
**Script:** `36_caf_subtyping.py`
**Reference:** Song et al. 2025, Peng et al. 2022

```
┌─────────────────────────────────────────────────────────────────┐
│ 3.1 CAF Identification                                          │
│ ─────────────────────────────────────────────────────────────── │
│ - Filter fibroblast population from annotations                 │
│ - Verify with COL1A1, PDGFRA expression                         │
├─────────────────────────────────────────────────────────────────┤
│ 3.2 Subtype Scoring                                             │
│ ─────────────────────────────────────────────────────────────── │
│ - mCAF (myofibroblast): ACTA2, TAGLN + aSMA(protein)            │
│ - iCAF (inflammatory): IL6, PDPN, PDGFRA                        │
│ - apCAF (antigen-presenting): CD74, HLA-DRA + HLA-DR(protein)   │
│ - Method: sc.tl.score_genes() for each signature                │
├─────────────────────────────────────────────────────────────────┤
│ 3.3 Spatial Distribution                                        │
│ ─────────────────────────────────────────────────────────────── │
│ - Map CAF subtypes spatially                                    │
│ - Question: Are iCAFs near immune cells? mCAFs near epithelium? │
│ - Compare CAF subtype proportions N→M→C                         │
├─────────────────────────────────────────────────────────────────┤
│ 3.4 G01 Deep Dive                                               │
│ ─────────────────────────────────────────────────────────────── │
│ - G01 has 36.7% fibroblasts - subtype these                     │
│ - Correlate with immune infiltration                            │
└─────────────────────────────────────────────────────────────────┘
```

**Deliverables:**
- [ ] CAF subtype assignments
- [ ] Subtype proportion by stage
- [ ] Spatial CAF subtype maps

---

### Phase 4: Spatial Domain Detection (NEW)
**Script:** `37_spatial_domains.py`
**Method:** SpatialGlue (recommended) or SpatialPCA

```
┌─────────────────────────────────────────────────────────────────┐
│ 4.1 SpatialGlue Integration                                     │
│ ─────────────────────────────────────────────────────────────── │
│ - Input: RNA + Protein + Spatial coordinates                    │
│ - Uses dual attention mechanism for spatial-aware clustering    │
│ - Should give better anatomical detail than WNN                 │
├─────────────────────────────────────────────────────────────────┤
│ 4.2 Domain Characterization                                     │
│ ─────────────────────────────────────────────────────────────── │
│ - Identify spatially coherent regions                           │
│ - Compare to WNN clusters                                       │
│ - Check if spatial domains capture tissue architecture better   │
├─────────────────────────────────────────────────────────────────┤
│ 4.3 Spatially Variable Genes                                    │
│ ─────────────────────────────────────────────────────────────── │
│ - Moran's I via squidpy (sq.gr.spatial_autocorr)                │
│ - Identify genes with spatial patterning                        │
│ - Compare SVGs across stages                                    │
└─────────────────────────────────────────────────────────────────┘
```

**Dependencies:**
```bash
pip install SpatialGlue  # If not installed
```

**Deliverables:**
- [ ] Spatial domain assignments
- [ ] WNN vs SpatialGlue comparison
- [ ] Top SVG list per stage

---

### Phase 5: Trajectory/Pseudotime Analysis (NEW)
**Script:** `38_trajectory_analysis.py`
**Method:** SpatialPCA or CellRank

```
┌─────────────────────────────────────────────────────────────────┐
│ 5.1 SpatialPCA Pseudotime                                       │
│ ─────────────────────────────────────────────────────────────── │
│ - Spatially-aware dimension reduction                           │
│ - Infer pseudotime along spatial gradients                      │
│ - Can capture N→M→C transitions within tissues                  │
├─────────────────────────────────────────────────────────────────┤
│ 5.2 Epithelial Trajectory                                       │
│ ─────────────────────────────────────────────────────────────── │
│ - Focus on epithelial cells only                                │
│ - Normal gastric → Metaplasia → Cancer                          │
│ - Identify intermediate states                                  │
├─────────────────────────────────────────────────────────────────┤
│ 5.3 Trajectory-Associated Genes                                 │
│ ─────────────────────────────────────────────────────────────── │
│ - Genes changing along pseudotime                               │
│ - Validate with known progression markers                       │
└─────────────────────────────────────────────────────────────────┘
```

**Dependencies:**
```bash
pip install SpatialPCA  # R package, may need rpy2
# Alternative: Use scanpy PAGA + diffusion pseudotime
```

**Deliverables:**
- [ ] Pseudotime assignments
- [ ] Trajectory visualization
- [ ] Trajectory-associated gene list

---

### Phase 6: Spatial Statistics (Enhanced)
**Script:** `39_spatial_statistics.py`
**Method:** Squidpy

```
┌─────────────────────────────────────────────────────────────────┐
│ 6.1 Neighborhood Enrichment                                     │
│ ─────────────────────────────────────────────────────────────── │
│ - sq.gr.nhood_enrichment per stage                              │
│ - Compare neighborhood patterns N→M→C                           │
│ - Key: immune cell neighbors of epithelial cells                │
├─────────────────────────────────────────────────────────────────┤
│ 6.2 Co-occurrence Analysis                                      │
│ ─────────────────────────────────────────────────────────────── │
│ - sq.gr.co_occurrence                                           │
│ - Which cell types co-localize at each stage?                   │
│ - Changes in co-occurrence patterns                             │
├─────────────────────────────────────────────────────────────────┤
│ 6.3 Interaction Distances                                       │
│ ─────────────────────────────────────────────────────────────── │
│ - Compute pairwise distances between cell types                 │
│ - Compare immune-tumor distance across stages                   │
│ - Quantify immune exclusion                                     │
├─────────────────────────────────────────────────────────────────┤
│ 6.4 Ripley's Statistics                                         │
│ ─────────────────────────────────────────────────────────────── │
│ - sq.gr.ripley for spatial clustering patterns                  │
│ - Are certain cell types spatially clustered?                   │
└─────────────────────────────────────────────────────────────────┘
```

**Deliverables:**
- [ ] Neighborhood enrichment heatmaps per stage
- [ ] Co-occurrence plots
- [ ] Quantified immune exclusion metrics

---

### Phase 7: IEX Signature Check (Limited)
**Script:** Part of `34_progression_analysis.py`

**Status:** Only 1/4 IEX genes available (TGFBI)

```
┌─────────────────────────────────────────────────────────────────┐
│ 7.1 TGFBI Expression                                            │
│ ─────────────────────────────────────────────────────────────── │
│ - Plot TGFBI expression across stages                           │
│ - Spatial distribution of TGFBI                                 │
│ - Correlate with immune infiltration                            │
│ - Note: Full IEX signature not testable (DDR1, PAK4, DPEP1 NA)  │
└─────────────────────────────────────────────────────────────────┘
```

**Note:** Document this limitation in any publication. The IEX signature from Heiser et al. 2023 requires DDR1, TGFBI, PAK4, DPEP1 - only TGFBI is in our panel.

---

## Execution Order

```
Day 1: Phase 1 (Core Progression)
├── 34_progression_analysis.py
├── Cell type proportions, DE, gastric markers
└── Checkpoint: Clear N→M→C changes?

Day 2: Phase 2+3 (Communication + CAF)
├── 35_cellcell_communication.py (LIANA+)
├── 36_caf_subtyping.py
└── Checkpoint: Meaningful LR changes? CAF subtypes identified?

Day 3: Phase 4+5 (Spatial + Trajectory)
├── 37_spatial_domains.py (SpatialGlue)
├── 38_trajectory_analysis.py
└── Checkpoint: Better spatial structure? Trajectory makes sense?

Day 4: Phase 6 (Spatial Statistics)
├── 39_spatial_statistics.py
├── Compile all results
└── Generate summary figures
```

---

## Key Figures to Generate

| Figure | Content | Phase |
|--------|---------|-------|
| Fig 1 | Spatial maps colored by cell type (N, M, C) | 1 |
| Fig 2 | Cell type proportion stacked bar | 1 |
| Fig 3 | Volcano plots (N vs M, M vs C) | 1 |
| Fig 4 | Gastric marker spatial expression | 1 |
| Fig 5 | **LIANA dotplot: changing LR pairs** | 2 |
| Fig 6 | **CAF subtype spatial distribution** | 3 |
| Fig 7 | **Spatial domains vs WNN clusters** | 4 |
| Fig 8 | **Pseudotime trajectory** | 5 |
| Fig 9 | Neighborhood enrichment heatmaps | 6 |
| Fig 10 | Immune exclusion quantification | 6 |

---

## Success Criteria

- [ ] Clear N→M→C cell type changes documented
- [ ] At least 50 significant DE genes per comparison
- [ ] Cell-cell communication changes identified (LIANA)
- [ ] CAF subtypes characterized (mCAF/iCAF/apCAF proportions)
- [ ] Spatial domains improve on WNN clustering
- [ ] Epithelial trajectory captures progression
- [ ] Spatial statistics show immune exclusion in cancer

---

## Dependencies

```bash
conda activate enact

# Already installed
# scanpy, squidpy, muon, liana

# May need to install
pip install SpatialGlue  # Spatial domain detection
pip install scvi-tools   # totalVI (if benchmarking)
```

---

## References (from Polymath)

1. **Cell-cell communication:** Jin et al. 2025 "CellChat for systematic analysis" Nat Protoc
2. **CAF subtypes:** Song et al. 2025 "Antigen-presenting CAFs in gastric cancer"
3. **Spatial domains:** Shang & Zhou 2022 "SpatialPCA" Nat Commun
4. **IEX signature:** Heiser et al. 2023 "Molecular cartography" Cell (DDR1, TGFBI, PAK4, DPEP1)
5. **Spatial statistics:** Palla et al. 2022 "Squidpy" Nat Methods

---

## Notes

- **Data location:** `/mnt/x/Choi_Batch_2_Tuesday/`
- **Results:** `~/spatial-hackathon-2026/results/g4x_choi_batch2/`
- **GitHub:** https://github.com/vanbelkummax/g4x-choi-batch2-analysis

See also:
- `VALIDATION.md` - Resolve comparison results
- `METHOD_COMPARISON.md` - Competitive landscape analysis
