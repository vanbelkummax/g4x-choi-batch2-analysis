# G4X Gastric Cancer: Evidence-Based Project Prioritization

**Generated:** 2026-01-21
**Data Sources:** Polymath (2,290 papers, 578K code chunks), Vanderbilt Professors (830 papers), Cross-domain algorithms (50K)

---

## Executive Summary

### Confirmed Literature Gaps (High Novelty Opportunities)

| Gap | Evidence | Polymath Hits |
|-----|----------|---------------|
| Gastric pre-cancer spatial transcriptomics | Only 1 review mentions it (Liu 2026) | 0 primary studies |
| TLS in gastric/stomach | 15 TLS spatial papers, 0 gastric-specific | Gap confirmed |
| Multimodal RNA+protein same-section gastric | No hits | Gap confirmed |
| CD8 exhaustion spatial in gastric | Methods exist (Hwang 2022), no gastric application | Gap confirmed |

### Key Vanderbilt Alignments

| Professor | Paper | PMID | Relevance |
|-----------|-------|------|-----------|
| **Hwang** | Conserved spatial CAF subtypes | 40154487 | Direct methodology for CAF subtyping |
| **Hwang** | ACTA2 predicts ICI response in gastric | 36508166 | We have aSMA protein! |
| **Hwang** | Spatial TIL analysis predicts survival | - | TIL methodology |
| **Lau** | CAF-TME crosstalk in CRC | 35794563 | Similar tissue, methods applicable |
| **Sarkar** | SIID joint imputation/deconvolution | 41249066 | Method for sparse spatial data |

### Cross-Domain Algorithm Opportunities

| Algorithm | Domain | Spatial Application | Polymath Mentions |
|-----------|--------|---------------------|-------------------|
| Persistent Homology | Topology | TLS/aggregate detection | 176 |
| Community Detection | Graph Theory | Spatial niches | 173 |
| Optimal Transport | Mathematics | Sample alignment | 3,830 |
| Entropy | Information Theory | Spatial gene selection | 171 |

---

## TIER 1: Maximum Impact (Do All)

### 1. First Spatial Atlas of Gastric Pre-Cancer Progression

**Score:** Significance ★★★★★ × Feasibility ★★★★★ × Novelty ★★★★★ = **125**

| Evidence | Source |
|----------|--------|
| Gap confirmed | Polymath: 0 primary gastric pre-cancer spatial studies |
| Only mention | Liu 2026 review: "spatial analyses identified a subset dominated by intestinal stem-cell lineages with transcriptomic convergence toward gastric cancer" |
| Methodology | Min 2019 (Lau lab): Kras-induced gastric dysplasia heterogeneity |

**Your Advantage:** 8 samples spanning Normal → Metaplasia → Cancer with 522K cells

**Implementation:**
```python
# Polymath-grounded analysis
mcp__polymath-v4__search_papers(query="gastric intestinal metaplasia single cell atlas", n=20)

# Key comparison
stages = adata.obs['stage'].unique()  # Normal, Metaplasia, Cancer
for stage in stages:
    subset = adata[adata.obs['stage'] == stage]
    # Cell type proportions, spatial statistics per stage
```

**Output:** Figure 2 of paper - cell composition + spatial organization per stage

---

### 2. CD8 Exhaustion Proximity Networks (Hwang Method Extension)

**Score:** Significance ★★★★★ × Feasibility ★★★★★ × Novelty ★★★★ = **100**

| Evidence | Source |
|----------|--------|
| Method validated | Hwang 2022 (CRC): "Spatial analysis of tumor-infiltrating lymphocytes predicts survival" |
| No gastric application | Polymath: 0 gastric CD8 exhaustion spatial studies |
| Your markers | PD1, PDL1, CD8 proteins + LAG3, HAVCR2, TIGIT RNA |

**Key Finding from Literature:**
> "PD1+ CD8 T cells to PDL1+ Tumor or Macrophage cells interaction score is often a better predictor of immunotherapy response than simple abundance counts" - Comprehensive Analysis Framework (Polymath)

**Implementation:**
```python
# Distance-based exhaustion scoring
exhausted_markers = ['PDCD1', 'LAG3', 'HAVCR2', 'TIGIT', 'CXCL13']
sc.tl.score_genes(adata, exhausted_markers, score_name='exhaustion_score')

# Proximity analysis
from scipy.spatial.distance import cdist
pdl1_positive = adata[adata.obs['PDL1_protein'] > threshold]
cd8_cells = adata[adata.obs['cell_type'] == 'CD8_T']
distances = cdist(cd8_cells.obsm['spatial'], pdl1_positive.obsm['spatial'])
```

**Output:** Figure 4 - CD8 exhaustion proximity score vs PDL1+ cell distance

---

### 3. TLS Detection in Gastric Pre-Cancer (Novel Territory)

**Score:** Significance ★★★★★ × Feasibility ★★★★ × Novelty ★★★★★ = **100**

| Evidence | Source |
|----------|--------|
| Gap confirmed | Polymath: 15 TLS spatial papers, 0 in gastric |
| Method exists | Tang 2025 (Cancer Cell): "Spatial transcriptomics reveals tryptophan metabolism restricting TLS maturation" |
| Key reference | Hwang 2025 CAF paper: "s4-CAFs interwoven with B/T cell aggregates and TLSs" |

**TLS Detection Criteria (from literature):**
- CD20+ B cell clusters
- Surrounded by CD4+/CD8+ T cells
- Often associated with s4-CAF subtype (Hwang)

**Cross-Domain Opportunity:**
> **Persistent Homology** (176 mentions): Use topological data analysis to detect TLS-like aggregates as "holes" in the cell-type spatial graph

**Implementation:**
```python
# B-T co-localization for TLS
b_cells = adata[adata.obs['cell_type'] == 'B_cell']
t_cells = adata[adata.obs['cell_type'].isin(['CD4_T', 'CD8_T'])]

# Squidpy co-occurrence
sq.gr.co_occurrence(adata, cluster_key="cell_type")
sq.pl.co_occurrence(adata, cluster_key="cell_type", clusters=['B_cell', 'CD4_T', 'CD8_T'])

# Persistent homology (novel approach)
# Detect 1-dimensional holes = B-T aggregates
from ripser import ripser
from persim import plot_diagrams
```

**Output:** Figure 5 - TLS detection algorithm + density per progression stage

---

### 4. Spatial CAF Subtyping (Hwang 2025 Replication + Gastric Extension)

**Score:** Significance ★★★★ × Feasibility ★★★★★ × Novelty ★★★★ = **80**

| Evidence | Source |
|----------|--------|
| Direct methodology | Hwang 2025 (PMID:40154487): 4 conserved spatial CAF subtypes |
| Gastric relevance | Hwang 2023 (PMID:36508166): "ACTA2 predicts ICI response in gastric cancer" |
| Your markers | ACTA2, TAGLN, IL6, PDPN, PDGFRA, CD74, HLA-DR ✓ |

**Hwang's 4 Spatial CAF Subtypes:**

| Subtype | Markers | Neighborhood | Your Panel |
|---------|---------|--------------|------------|
| s1-CAF | ACTA2high, dense | Tumor boundary | ✓ ACTA2, aSMA protein |
| s2-CAF | ECM genes | Stroma-rich | ✓ COL1A1, FN1 |
| s3-CAF | PDGFRBhigh | Myeloid niches | ✓ PDGFRA |
| s4-CAF | CD74, HLA-II | TLS-associated | ✓ CD74, HLA-DR protein |

**Key Quote:**
> "s4-CAFs were interwoven with clusters of immune cells, including B/T cell aggregates and TLSs, suggesting their potential role in orchestrating anti-tumor immunity" - Hwang 2025

**Implementation:**
```python
# Hwang CAF signature scoring
caf_signatures = {
    's1_mCAF': ['ACTA2', 'TAGLN', 'MYL9'],  # + aSMA protein
    's2_eCAF': ['COL1A1', 'COL3A1', 'FN1'],
    's3_iCAF': ['IL6', 'PDGFRA', 'CFD'],
    's4_apCAF': ['CD74', 'HLA-DRA', 'HLA-DRB1']  # + HLA-DR protein
}

for subtype, genes in caf_signatures.items():
    available = [g for g in genes if g in adata.var_names]
    sc.tl.score_genes(adata, available, score_name=subtype)
```

**Output:** Figure 3 - CAF subtype spatial maps + proportions by stage + survival

---

### 5. Immune Exclusion Gradient Quantification

**Score:** Significance ★★★★ × Feasibility ★★★★★ × Novelty ★★★★ = **80**

| Evidence | Source |
|----------|--------|
| Method validated | "Molecular cartography" paper: IEX (immune exclusion) 4-gene signature |
| Clinical relevance | "IEX had negative correlation with cytotoxic T cells and positive with Tregs" |
| Your data | Spatial coordinates + cell types = direct measurement |

**IEX Signature Genes (from Polymath):**
- TGFBI, DDR1, PAK4, DPEP1
- **Your panel:** Only TGFBI available (1/4)

**Alternative Approach:** Direct spatial distance measurement (more robust)

**Implementation:**
```python
# Immune-to-epithelial distance by stage
from scipy.spatial import cKDTree

epithelial = adata[adata.obs['lineage'] == 'Epithelial']
immune = adata[adata.obs['lineage'] == 'Immune']

tree = cKDTree(epithelial.obsm['spatial'])
distances, _ = tree.query(immune.obsm['spatial'], k=1)

# Compare per stage
for stage in ['Normal', 'Metaplasia', 'Cancer']:
    stage_distances = distances[immune.obs['stage'] == stage]
    print(f"{stage}: mean distance = {np.mean(stage_distances):.2f}")
```

**Output:** Figure 6 - Immune exclusion gradient visualization + statistics

---

## TIER 2: High Value If Time Permits

### 6. RNA-Protein Discordance as Biological Signal

**Score:** 75 | Gap: HIGH | Method: NOVEL

**Insight:** Your r=0.088 RNA-protein correlation is LOW. Instead of treating as QC failure, ask: **where do they disagree and why?**

| Evidence | Source |
|----------|--------|
| Biological meaning | Post-transcriptional regulation, protein stability |
| No spatial studies | Polymath: 0 papers on spatial RNA-protein discordance |

**Hypothesis:** High discordance regions = active post-transcriptional regulation = biologically important

---

### 7. Neighborhood Signatures for Stage Classification

**Score:** 64 | Method: Community Detection (173 mentions)

**Approach:** Use Louvain clustering on spatial neighborhood graphs to define "niche signatures"

```python
# Build spatial graph
sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=15)

# Neighborhood enrichment patterns as features
sq.gr.nhood_enrichment(adata, cluster_key="cell_type")

# Train classifier: neighborhood pattern → stage
from sklearn.ensemble import RandomForestClassifier
```

---

### 8. Treg/CD8 Ratio Spatial Gradients

**Score:** 48 | Direct measurement

**Your data:** CD8/Treg ratio = 7.3 average. Map spatially - does ratio change with distance from tumor?

---

### 9. PDL1+ Tumor Cell Clustering (Hot/Cold Regions)

**Score:** 48 | ICI prediction

**Question:** Are PDL1+ tumor cells clustered (immune hot) or dispersed (immune cold)?

```python
# Ripley's K for PDL1+ clustering
sq.gr.ripley(adata[adata.obs['PDL1_positive']], cluster_key='cell_type', mode='K')
```

---

### 10. Cross-Sample Spatial Domain Transfer

**Score:** 48 | Method: Optimal Transport (3,830 mentions)

**Question:** Do spatial niches transfer across patients? Use optimal transport for sample alignment.

---

## TIER 3: Quick Wins (Supplementary Material)

| # | Idea | Implementation | Output |
|---|------|----------------|--------|
| 11 | Protein QC benchmark | Document HLA-DR r=0.31 vs KI67 r=0.003 | Supp Table |
| 12 | Cell composition stats | Chi-square + effect sizes N/M/C | Table 1 |
| 13 | Neighborhood enrichment per stage | `sq.gr.nhood_enrichment()` × 3 | Supp Fig |
| 14 | Moran's I for key genes | `sq.gr.spatial_autocorr()` | Supp Fig |
| 15 | KI67 proliferation vs immune | Proliferation hotspots + immune exclusion | Supp Fig |

---

## TIER 4: Future Work (AI/ML)

| # | Idea | Feasibility Blocker | Cross-Domain Algorithm |
|---|------|---------------------|------------------------|
| 16 | GNN cell type prediction | Training infrastructure | Graph Convolution (109 mentions) |
| 17 | Contrastive RNA-protein learning | Novel method | Optimal Transport alignment |
| 18 | Stage classifier from spatial graphs | Sample size (n=8) | Community Detection |
| 19 | scGPT fine-tuning | 357 genes may be too few | - |
| 20 | Persistent homology for tissue architecture | Novel application | Persistent Homology (176 mentions) |

---

## Paper Strategy: "Gastric Immune Atlas"

**Recommended Title:**
> *Spatial multi-omic atlas of gastric pre-cancer progression reveals CAF niche reorganization and emergent tertiary lymphoid structures*

### Figure Plan

| Fig | Content | Ideas | Key Reference |
|-----|---------|-------|---------------|
| 1 | Study design + QC + validation | 12 | - |
| 2 | Cell composition N→M→C + spatial | 1 | Min 2019 |
| 3 | Spatial CAF subtypes (Hwang method) | 4 | Hwang 2025 |
| 4 | CD8 exhaustion proximity | 2 | Hwang 2022 |
| 5 | TLS detection in gastric | 3 | Tang 2025 |
| 6 | Immune exclusion gradients | 5 | Molecular cartography |

### Target Journals (ranked by fit)

1. **Cancer Cell** - Hwang publishes here, CAF focus
2. **Gut** - IF 24.5, gastric expertise
3. **Gastroenterology** - IF 29.4, clinical impact
4. **Nature Communications** - IF 16.6, broader reach

---

## Key Polymath Queries for Each Phase

```python
# Phase 1: Progression (run before implementing)
mcp__polymath-v4__search_papers(query="gastric metaplasia cancer progression single cell spatial", n=20)

# Phase 2: CAF subtyping
mcp__vanderbilt-professors__search_hwang_papers(query="CAF spatial subtype", limit=10)

# Phase 3: TLS detection
mcp__polymath-v4__search_papers(query="tertiary lymphoid structure detection algorithm spatial", n=15)

# Phase 4: Immune exclusion
mcp__polymath-v4__search_papers(query="immune exclusion score spatial tumor infiltrating", n=15)

# Cross-domain insights
mcp__polymath-v4__search_algorithms(bridges=True, query="cell aggregation detection topology")
```

---

## Validation Checklist Before Starting

```bash
# Verify CAF marker availability
python -c "
import scanpy as sc
adata = sc.read('results/g4x_choi_batch2/combined_wnn.h5ad')
caf_markers = ['ACTA2', 'TAGLN', 'IL6', 'PDGFRA', 'CD74', 'HLA-DRA']
for m in caf_markers:
    status = '✓' if m in adata.var_names else '✗'
    print(f'{status} {m}')
"

# Verify TLS markers
python -c "
import scanpy as sc
adata = sc.read('results/g4x_choi_batch2/combined_wnn.h5ad')
tls_markers = ['MS4A1', 'CD19', 'CD4', 'CD8A', 'CR2']  # CD20=MS4A1
for m in tls_markers:
    status = '✓' if m in adata.var_names else '✗'
    print(f'{status} {m}')
"
```

---

## Summary: What Changed from Original Brainstorm

| Original Assumption | Evidence Found | Action |
|---------------------|----------------|--------|
| "Zero gastric TLS papers" | Confirmed: 0 in Polymath | Priority UP |
| "CAF subtyping is novel" | Hwang 2025 defines 4 subtypes | Replication + extension |
| "CD8 exhaustion methods unclear" | Hwang 2022 has exact method | Direct application |
| "IEX signature available" | Only TGFBI (1/4 genes) | Use distance method instead |
| "TLS detection unclear" | s4-CAF associated with TLS (Hwang) | Link to CAF analysis |
| "Persistent homology applicable?" | 176 mentions, used for tissue architecture | Add as novel method |

**Key Insight:** Hwang lab papers provide direct methodological foundation. Frame as "Hwang CAF framework applied to gastric pre-cancer" for credibility.
