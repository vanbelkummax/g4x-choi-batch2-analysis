# G4X Gastric Cancer Progression: Tiered Analysis Plan

**Generated:** 2026-01-21
**Dataset:** Choi_GC_preGC_batch2 (G4X, Resolve Biosciences)
**Status:** QC Complete, Ready for Downstream Analysis

---

## Executive Summary

This plan prioritizes analyses based on:
1. **Scientific impact** - Novel insights into gastric cancer progression
2. **Data suitability** - Marker panel coverage (341 RNA + 17 protein)
3. **Methodological maturity** - Established vs. experimental approaches
4. **Resource efficiency** - Computational cost vs. insight yield

### Dataset Overview
| Metric | Value |
|--------|-------|
| Samples | 29 (of 32 passing QC) |
| Cells | 1,835,026 |
| Disease stages | Control, Normal, Metaplasia, Cancer |
| Progression series | SNU-105 (complete), SNU-107 (complete), SNU-484 (partial) |
| RNA genes | 341 (targeted panel) |
| Protein markers | 17 |

### Key Strengths
- **Matched progression series**: Same patient N→M→C enables true trajectory analysis
- **Multimodal data**: RNA + protein enables WNN integration and validation
- **High cell count**: 1.8M cells provides statistical power for rare populations
- **Spatial resolution**: Single-molecule G4X enables subcellular analysis

### Key Limitations
- **Targeted panel**: 341 genes limits transcriptome-wide discovery
- **IEX signature**: Only TGFBI available (1/4 markers)
- **No H&E**: Limits morphological validation

---

## 🔥 Tier 0: PCA Deep Dive (IMMEDIATE PRIORITY)

**User-requested priority analysis.** PCA-based exploration provides foundation for all downstream work.

### 0.1 Sample-Level Pseudobulk PCA
**Status:** Not started
**Script:** `46_pca_deep_dive.py` (to create)
**Priority:** 🔥 IMMEDIATE

**Rationale:** Understand global sample relationships before cell-level analyses. Pseudobulk aggregation reduces noise and reveals sample-level patterns.

**Methods:**
- Aggregate counts per sample (pseudobulk)
- PCA on log-normalized pseudobulk matrix
- Color by: stage, patient, lane (batch)
- Scree plot for variance explained

**Deliverables:**
- [ ] PC1 vs PC2 scatter (colored by stage)
- [ ] PC1 vs PC2 scatter (colored by patient)
- [ ] Scree plot (variance explained)
- [ ] Top 20 PC1/PC2 loading genes
- [ ] Figure: Sample trajectory in PCA space

---

### 0.2 Cell-Level PCA Analysis
**Status:** Not started
**Script:** `46_pca_deep_dive.py`
**Priority:** 🔥 IMMEDIATE

**Rationale:** Understand cell population structure and how it relates to disease stages.

**Methods:**
- Use existing `X_pca` from merged_corrected.h5ad
- Visualize PC1-PC4 colored by: cell type, stage, patient
- Density plots per stage in PC space
- PC loading analysis for biological interpretation

**Deliverables:**
- [ ] PC1 vs PC2 colored by cell type
- [ ] PC1 vs PC2 colored by stage (with density contours)
- [ ] PC loadings heatmap (top genes per PC)
- [ ] Stage separation quantification in PC space

---

### 0.3 Multimodal PCA Comparison
**Status:** Not started
**Script:** `46_pca_deep_dive.py`
**Priority:** 🔥 IMMEDIATE

**Rationale:** Compare information content of RNA vs Protein vs WNN integration.

**Methods:**
- PCA on RNA only
- PCA on Protein only
- PCA on WNN integrated space
- Compare: Which modality best separates stages?

**Deliverables:**
- [ ] Side-by-side PCA plots (RNA / Protein / WNN)
- [ ] Procrustes analysis (alignment between modalities)
- [ ] Variance explained comparison
- [ ] Modality contribution to stage separation

---

### 0.4 Batch Effect Visualization
**Status:** Not started
**Script:** `46_pca_deep_dive.py`
**Priority:** 🔥 IMMEDIATE

**Rationale:** Verify Harmony batch correction worked; understand residual batch effects.

**Methods:**
- PCA before Harmony (X_pca)
- PCA after Harmony (X_pca_harmony)
- Color by lane (batch variable)
- LISI scores visualization

**Deliverables:**
- [ ] Pre-correction PCA colored by lane
- [ ] Post-correction PCA colored by lane
- [ ] LISI improvement visualization
- [ ] Batch mixing quantification

---

### 0.5 MOFA+ Multi-Omics Factor Analysis
**Status:** Not started
**Script:** `47_mofa_analysis.py` (to create)
**Priority:** HIGH (after basic PCA)

**Rationale:** MOFA+ identifies shared and modality-specific factors, more sophisticated than simple PCA.

**Methods:**
- MOFA+ with RNA and Protein views
- Factor interpretation (which genes/proteins drive each factor?)
- Factor correlation with clinical variables (stage, patient)

**Deliverables:**
- [ ] MOFA factor loadings
- [ ] Factor-stage associations
- [ ] Shared vs modality-specific variance
- [ ] Top features per factor

**Polymath References:**
- MOFA+ (Argelaguet et al.) - Multi-omics factor analysis
- SMOPCA 2025 - Spatially aware multimodal PCA (doc_id: 47f141c8)

---

## Tier 1: High-Priority Core Analyses (Week 1-2)

These analyses directly address the primary hypothesis and have full marker support.

### 1.1 Cell Type Proportion Analysis Across Progression
**Status:** Partially complete (preliminary)
**Script:** `34_progression_analysis.py` (to create)
**Priority:** CRITICAL

**Rationale:** Establishes the cellular landscape changes N→M→C. Foundation for all downstream analyses.

**Methods:**
- Chi-squared tests for proportion differences
- Dirichlet regression for compositional analysis (recommended by Lau lab CRC atlas)
- Visualization: Stacked bar plots, alluvial diagrams

**Deliverables:**
- [ ] Cell type proportions per sample/stage
- [ ] Statistical significance of proportion changes
- [ ] Epithelial/immune/stromal ratio trends
- [ ] Figure: Alluvial diagram of cellular transitions

**Polymath References:**
- Lin et al. 2023 (PMID:36669472) - CRC multiplexed atlas methodology
- Liu et al. 2026 - Spatial omics comprehensive review

---

### 1.2 Gastric-to-Intestinal Metaplasia Marker Analysis
**Status:** Not started
**Script:** `34_progression_analysis.py`
**Priority:** CRITICAL

**Rationale:** Core biological question - spatial distribution of intestinalization markers.

**Available Markers:**
| Marker | Function | Panel Status |
|--------|----------|--------------|
| CDX2 | Intestinal TF | ✅ Available |
| MUC2 | Goblet cell mucin | ✅ Available |
| MUC5AC | Gastric mucin | ✅ Available |
| TFF1/TFF2/TFF3 | Trefoil factors | ✅ Available |
| SOX9 | Stem cell marker | ✅ Available |
| LGR5 | Stem cell marker | ✅ Available |

**Methods:**
- Spatial expression plots per marker
- Co-expression analysis (CDX2+/MUC2+ vs MUC5AC+)
- Gradient analysis from normal to metaplastic regions

**Deliverables:**
- [ ] Spatial marker expression maps per sample
- [ ] Intestinalization score per cell
- [ ] Transition zone characterization
- [ ] Figure: Marker gradient heatmaps

---

### 1.3 CAF Subtyping and Spatial Organization
**Status:** Planned in roadmap
**Script:** `36_caf_subtyping.py` (to create)
**Priority:** HIGH

**Rationale:** CAFs drive tumor progression. Recent Liu et al. 2025 paper identifies conserved spatial CAF subtypes across cancers.

**Available CAF Markers:**
| Subtype | Markers | Panel Status |
|---------|---------|--------------|
| myCAF | ACTA2, TAGLN, MYL9 | ✅ Full |
| iCAF | IL6, CXCL12, PDGFRA | ✅ Full |
| apCAF | CD74, HLA-DRA | ✅ Full |
| General | FAP, COL1A1, VIM | ✅ Full |

**Methods (from Liu et al. 2025):**
1. Score cells using CAF subtype signatures
2. Spatial domain detection for CAF neighborhoods
3. Neighborhood composition analysis (what cells surround each CAF subtype?)
4. CAF-immune interaction quantification

**Deliverables:**
- [ ] CAF subtype assignments per fibroblast
- [ ] Spatial distribution maps
- [ ] CAF neighborhood composition
- [ ] CAF proportion changes across progression
- [ ] Figure: CAF subtype spatial organization

**Polymath References:**
- Liu et al. 2025 - "Conserved spatial subtypes and cellular neighborhoods of CAFs" (doc_id: 5a68fd66)
- Yang et al. 2023 - "CAFs: from basic science to anticancer therapy" (doc_id: ce9f50d4)

---

### 1.4 Cell-Cell Communication Analysis
**Status:** Script exists (13_ligrec_communication.py)
**Script:** `35_cellcell_communication.py` (updated)
**Priority:** HIGH

**Rationale:** Understand how cell populations signal to each other across progression stages.

**Available LR Pairs (23/26):**
- Immune: CD274-PDCD1, CD80-CD28, CXCL9/10/11-CXCR3
- Stromal: COL1A1-ITGB1, FN1-ITGAV
- Growth: TGFB1-TGFBR1/2, VEGFA-KDR
- CAF-Immune: IL6-IL6R, CXCL12-CXCR4

**Methods:**
- LIANA+ with spatial constraints (distance-weighted)
- Stage-specific communication networks
- Sender-receiver analysis per cell type

**Deliverables:**
- [ ] LR interaction scores per stage
- [ ] Communication network visualizations
- [ ] Top differential interactions N→M→C
- [ ] CAF→Immune interaction dotplots
- [ ] Figure: Circos plots of communication patterns

**Polymath References:**
- CellChat and LIANA methodology papers (multiple in database)

---

## Tier 2: Advanced Spatial Analyses (Week 2-3)

These require more complex methods but offer unique spatial insights.

### 2.1 Spatial Statistics (Squidpy-based)
**Status:** Planned
**Script:** `39_spatial_statistics.py` (to create)
**Priority:** HIGH

**Methods:**
| Analysis | Function | Purpose |
|----------|----------|---------|
| Neighborhood enrichment | `sq.gr.nhood_enrichment()` | Cell type co-localization |
| Co-occurrence | `sq.gr.co_occurrence()` | Spatial correlation |
| Ripley's K/L | `sq.gr.ripley()` | Clustering patterns |
| Moran's I | `sq.gr.spatial_autocorr()` | Spatially variable genes |
| Centrality | `sq.gr.centrality_scores()` | Network topology |

**Deliverables:**
- [ ] Neighborhood enrichment heatmaps per stage
- [ ] Spatially variable gene list
- [ ] Cell type clustering statistics
- [ ] Figure: Co-occurrence matrices

**Cross-Domain Algorithm (Polymath):**
- **Community Detection** (graph theory) → Spatial cell communities
- **Louvain Algorithm** → Neighborhood clustering

---

### 2.2 Spatial Domain Detection
**Status:** Planned
**Script:** `37_spatial_domains.py` (to create)
**Priority:** MEDIUM-HIGH

**Rationale:** Identify tissue architecture domains that transcend individual cell types.

**Methods (from SMOPCA, SpatialGlue papers):**
1. **SpatialGlue** - Multimodal (RNA+protein) graph neural network
2. **GraphST** - Graph-based spatial domains
3. **MISO** - Multimodal integration (if needed)

**Deliverables:**
- [ ] Spatial domain assignments per cell
- [ ] Domain marker genes
- [ ] Domain transitions across progression
- [ ] Figure: Domain maps overlaid on tissue

**Polymath References:**
- Shang & Zhou 2022 - "Spatially aware dimension reduction" (doc_id: 49aa60f3)
- SMOPCA 2025 - Multimodal spatial domain detection (doc_id: 47f141c8)
- MISO 2025 - Multimodal spatial omics modeling (doc_id: 2766cb8f)

---

### 2.3 Immune Microenvironment Characterization
**Status:** Partially explored in pilot
**Script:** `38_immune_analysis.py` (to create)
**Priority:** HIGH

**Key Questions:**
- How does the CD8/Treg ratio change spatially across progression?
- Where are exhausted vs. functional T cells located?
- Is there immune exclusion in cancer samples?

**Available Immune Markers:**
| Population | Markers | Status |
|------------|---------|--------|
| CD8 T cells | CD8A, GZMB, PRF1 | ✅ |
| Tregs | FOXP3, IL2RA | ✅ |
| Macrophages | CD68, CD163, CD14 | ✅ |
| B cells | CD19, MS4A1 | ✅ |
| Exhaustion | PDCD1, CTLA4, LAG3, TIGIT | ✅ |

**Methods:**
- Exhaustion scoring per T cell
- Spatial mapping of immune states
- Immune exclusion analysis (tumor-immune distance)

**Deliverables:**
- [ ] Immune phenotype maps
- [ ] CD8/Treg spatial ratio
- [ ] Exhaustion gradient analysis
- [ ] Figure: Immune exclusion streamplots

---

## Tier 3: Trajectory and Dynamics (Week 3-4)

These analyses require careful interpretation given the cross-sectional design.

### 3.1 Pseudotime/CellRank Trajectory Analysis
**Status:** Planned
**Script:** `41_trajectory_analysis.py` (to create)
**Priority:** MEDIUM

**Rationale:** Model the N→M→C transition as a continuous trajectory.

**Caution:** Cross-sectional data; trajectories are inferred, not observed.

**Methods:**
- CellRank for fate probability estimation
- RNA velocity (if splice information available - CHECK)
- Pseudotime ordering within epithelial compartment

**Deliverables:**
- [ ] Pseudotime orderings
- [ ] Fate probability maps
- [ ] Transition driver genes
- [ ] Figure: UMAP with trajectory overlay

**Polymath Cross-Domain:**
- **Optimal Transport** (3830 mentions) → Cell fate trajectory via STORIES method
- Reference: "STORIES: learning cell fate landscapes" (doc_id: b91c3271)

---

### 3.2 TLS (Tertiary Lymphoid Structure) Detection
**Status:** Novel analysis
**Script:** `40_tls_detection.py` (to create)
**Priority:** MEDIUM

**Rationale:** TLS presence correlates with immunotherapy response (Cabrita et al., Nature 2020).

**Methods:**
1. Identify B+T cell aggregates
2. **Persistent Homology** for ring-shaped structure detection (H1 features)
3. Validate with CD21 (FDC marker) if available

**Deliverables:**
- [ ] TLS candidate regions
- [ ] TLS density per stage
- [ ] TLS maturity scoring
- [ ] Figure: TLS spatial maps

**Polymath Cross-Domain:**
- **Persistent Homology** (176 mentions, topology domain)
- Spatial uses: "tissue architecture analysis, tumor boundary detection"
- Reference: Mebane et al. 2025 - TLS spatial transcriptomics (from search results)

---

### 3.3 Multimodal PCA and Integration Validation
**Status:** Script exists (14_pca_analyses.py)
**Script:** `42_multimodal_validation.py` (to update)
**Priority:** MEDIUM

**Methods:**
- Pseudobulk PCA per sample
- MOFA+ for multi-omics factor analysis
- RNA-protein correlation validation

**Deliverables:**
- [ ] Sample-level PCA trajectories
- [ ] Key factors driving progression
- [ ] RNA-protein agreement metrics
- [ ] Figure: Sample trajectory plots

---

## Tier 4: Exploratory/Novel Methods (Week 4+)

High-risk, high-reward analyses requiring method development.

### 4.1 Spatial Entropy and Information Theory
**Status:** Novel
**Script:** `43_spatial_entropy.py` (to create)
**Priority:** LOW-MEDIUM

**Rationale:** Quantify spatial organization using information theory.

**Methods:**
- Spatial entropy of cell type distributions
- Mutual information between cell types
- Complexity measures across progression

**Polymath Cross-Domain:**
- **Entropy** (171 mentions, information_theory domain)
- Spatial uses: "gene selection, spatial information quantification"

---

### 4.2 Graph-Based Network Topology
**Status:** Planned in roadmap
**Script:** `44_network_topology.py` (to create)
**Priority:** LOW-MEDIUM

**Methods:**
- Delaunay triangulation of cell positions
- Betweenness centrality per cell type
- Network motif analysis

**Deliverables:**
- [ ] Cell connectivity graphs
- [ ] Hub cell identification
- [ ] Network topology changes N→M→C

---

### 4.3 Deep Dive: G01 (High Fibroblast Sample)
**Status:** Identified in QC
**Script:** `45_g01_deep_dive.py` (to create)
**Priority:** LOW

**Rationale:** G01 has 36.7% fibroblasts - potential desmoplastic phenotype.

**Methods:**
- Focused CAF analysis
- Comparison to other cancer samples
- Stromal-epithelial interface characterization

---

## Execution Schedule

| Week | Tier | Analyses | Scripts |
|------|------|----------|---------|
| **NOW** | **0** | **🔥 PCA Deep Dive** | **46_pca_deep_dive.py** |
| NOW | 0 | MOFA+ multi-omics factors | 47_mofa_analysis.py |
| 1 | 1 | Cell proportions, Gastric markers | 34_progression_analysis.py |
| 1-2 | 1 | CAF subtyping, Cell-cell communication | 35, 36 |
| 2 | 2 | Spatial statistics, Immune analysis | 38, 39 |
| 2-3 | 2 | Spatial domains | 37_spatial_domains.py |
| 3 | 3 | Trajectory, TLS detection | 40, 41 |
| 3-4 | 3 | Multimodal validation | 42 |
| 4+ | 4 | Exploratory analyses | 43, 44, 45 |

---

## Methods Reference Summary

### From Polymath Semantic Search
| Paper | Year | Key Method | Relevance |
|-------|------|------------|-----------|
| Liu et al. - Spatial CAF subtypes | 2025 | Spatial CAF subtyping | CAF analysis |
| Lin et al. - CRC multiplexed atlas | 2023 | 3D imaging, TLS | Methodology template |
| Shang & Zhou - SPACE | 2022 | Spatial domain detection | Domain analysis |
| SMOPCA | 2025 | Multimodal spatial PCA | Integration |
| MISO | 2025 | Multimodal integration | Domain detection |
| STORIES | 2025 | Optimal transport trajectory | Pseudotime |

### From Polymath Cross-Domain Algorithms
| Algorithm | Domain | Mentions | Spatial Application |
|-----------|--------|----------|---------------------|
| Optimal Transport | optimal_transport | 3830 | Cell fate trajectory |
| Persistent Homology | topology | 176 | TLS ring detection |
| Community Detection | graph_theory | 173 | Spatial communities |
| Louvain Algorithm | graph_theory | 128 | Neighborhood clustering |
| Entropy | information_theory | 171 | Spatial organization quantification |

### From Vanderbilt Professors
| Paper | Professor | Relevance |
|-------|-----------|-----------|
| CRC spatial atlas (PMID:35794563) | Ken Lau | CAF-TME crosstalk |
| 3D CRC atlas (PMID:36669472) | Ken Lau | Methodology, TLS analysis |
| DMBT1 in dysplasia (PMID:40026233) | Ken Lau | Progression markers |

---

## Output Files Summary

All outputs to: `results/downstream_analysis/`

```
results/downstream_analysis/
├── tier1_core/
│   ├── cell_proportions/
│   ├── gastric_markers/
│   ├── caf_subtyping/
│   └── cell_communication/
├── tier2_spatial/
│   ├── spatial_statistics/
│   ├── spatial_domains/
│   └── immune_analysis/
├── tier3_dynamics/
│   ├── trajectory/
│   ├── tls_detection/
│   └── multimodal_validation/
├── tier4_exploratory/
│   ├── spatial_entropy/
│   ├── network_topology/
│   └── deep_dives/
└── figures/
    ├── main/          # Publication-ready
    └── supplementary/  # Supporting figures
```

---

## Risk Assessment

| Analysis | Risk | Mitigation |
|----------|------|------------|
| Trajectory analysis | Cross-sectional data | Focus on matched patient series |
| TLS detection | May not have TLS in gastric | Report null result as valid finding |
| Spatial domains | Panel may limit resolution | Use multimodal (RNA+protein) |
| IEX signature | Only 1/4 markers | Report as limitation, use proxy |

---

## Next Immediate Steps

1. **Create output directory structure**
2. **Start with 34_progression_analysis.py** (cell proportions + gastric markers)
3. **Load merged_corrected.h5ad** and verify cell type annotations
4. **Generate first set of progression figures**

---

*Plan generated using Polymath knowledge base (2,308 papers, 50K algorithms, 1.2M Neo4j nodes)*
