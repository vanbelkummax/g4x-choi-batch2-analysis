# G4X Gastric Cancer: Strategic Project Selection

## Executive Summary

**Dataset:** 522K cells, 8 samples, 387 RNA + 17 protein markers (G4X platform)
**Gap Found:** Zero spatial pre-cancer gastric studies in Polymath (2,300 papers)
**Recommended Focus:** Combine top 5 ideas into "Gastric Immune Atlas" paper

---

## Decision Matrix

```
                    FEASIBILITY
                 LOW    MED    HIGH
           ┌─────────────────────────┐
     HIGH  │  17     6,10   1,2,3,4,5│  ← PRIORITY
SIGNIFICANCE│  19    16,18   7,8,9   │
     LOW   │        20      11-15    │  ← Quick wins for supplement
           └─────────────────────────┘
```

---

## TIER 1: Core Paper (Do All 5)

### 1. Spatial Immune Architecture of Gastric Pre-Cancer Progression
**Score:** Sig⭐5 × Feas⭐5 × Nov⭐5 = **125**

| Aspect | Detail |
|--------|--------|
| Gap | Zero N→M→C spatial studies in stomach |
| Data | All 8 samples span progression stages |
| Output | Fig 2: Cell composition + spatial plots per stage |

**Implementation:**
```python
# Polymath grounding
mcp__polymath-v4__search_papers(query="gastric metaplasia immune infiltration spatial", n=10)

# Analysis (scripts/34_progression_analysis.py)
sq.gr.spatial_neighbors(adata, coord_type="generic")
sq.pl.spatial_scatter(adata, color="cell_type", shape="stage")
```

---

### 2. CD8 Exhaustion Proximity to PDL1+ Cells
**Score:** Sig⭐5 × Feas⭐5 × Nov⭐4 = **100**

| Aspect | Detail |
|--------|--------|
| Markers | PD1, PDL1, CD8 proteins ✓ |
| Finding | 25% exhaustion rate already detected |
| Output | Fig 4: Distance distributions, niche classification |

**Implementation:**
```python
# Polymath grounding
mcp__polymath-v4__search_papers(query="CD8 exhaustion PD1 proximity spatial tumor", n=10)

# Distance analysis
from scipy.spatial.distance import cdist
exhausted_cd8 = adata[adata.obs['phenotype'] == 'CD8_exhausted']
pdl1_cells = adata[adata.obs['PDL1_positive'] == True]
distances = cdist(exhausted_cd8.obsm['spatial'], pdl1_cells.obsm['spatial'])
```

---

### 3. TLS-like Structures in Gastric Metaplasia
**Score:** Sig⭐5 × Feas⭐4 × Nov⭐5 = **100**

| Aspect | Detail |
|--------|--------|
| Gap | Zero gastric TLS spatial papers in Polymath |
| Markers | CD20, CD4, CD8 for B-T aggregates ✓ |
| Output | Fig 5: TLS detection, density per stage |

**Implementation:**
```python
# Polymath grounding
mcp__polymath-v4__search_papers(query="tertiary lymphoid structures spatial transcriptomics", n=10)
mcp__polymath-v4__search_code(query="tertiary lymphoid structure detection", n=5)

# TLS detection (B-T co-localization)
b_cells = adata[adata.obs['cell_type'] == 'B_cell']
t_cells = adata[adata.obs['cell_type'].isin(['CD4_T', 'CD8_T'])]
sq.gr.co_occurrence(adata, cluster_key="cell_type")
```

---

### 4. Spatial CAF Subtyping (Hwang Framework)
**Score:** Sig⭐4 × Feas⭐5 × Nov⭐4 = **80**

| Aspect | Detail |
|--------|--------|
| Reference | Hwang 2025 (PMID:40154487) - 4 CAF subtypes |
| Markers | ACTA2, IL6, CD74, HLA-DR ✓ |
| Output | Fig 3: CAF subtype proportions by stage |

**Implementation:**
```python
# Polymath grounding - Hwang paper
mcp__vanderbilt-professors__search_hwang_papers(query="CAF spatial subtype", limit=5)

# CAF scoring (scripts/36_caf_subtyping.py)
caf_signatures = {
    'mCAF': ['ACTA2', 'TAGLN'],  # + aSMA protein
    'iCAF': ['IL6', 'PDPN', 'PDGFRA'],
    'apCAF': ['CD74', 'HLA-DRA']  # + HLA-DR protein
}
sc.tl.score_genes(adata, caf_signatures['mCAF'], score_name='mCAF_score')
```

---

### 5. Immune Exclusion Score
**Score:** Sig⭐4 × Feas⭐5 × Nov⭐4 = **80**

| Aspect | Detail |
|--------|--------|
| Question | Does cancer exclude immune cells? |
| Data | Spatial coords + cell types = direct measurement |
| Output | Fig 6: Distance gradients, exclusion heatmap |

**Implementation:**
```python
# Polymath grounding
mcp__polymath-v4__search_papers(query="immune exclusion spatial distance tumor", n=10)

# Distance to epithelium
epithelial = adata[adata.obs['lineage'] == 'Epithelial']
immune = adata[adata.obs['lineage'] == 'Immune']
# Compute per-stage distance distributions
```

---

## TIER 2: If Time Permits

| # | Idea | Score | Key Insight |
|---|------|-------|-------------|
| 6 | RNA-protein discordance as signal | 75 | r=0.088 is biology, not failure |
| 7 | Neighborhood signatures → stage classifier | 64 | ML on enrichment patterns |
| 8 | Treg/CD8 ratio spatial gradients | 48 | Already have 7.3 ratio |
| 9 | PDL1+ tumor clustering (hot/cold) | 48 | ICI response prediction |
| 10 | Cross-sample domain transfer | 48 | SpatialGlue validation |

---

## TIER 3: Quick Wins (Supplement Material)

| # | Idea | Implementation |
|---|------|----------------|
| 11 | Protein QC benchmark | Document HLA-DR r=0.31 vs KI67 r=0.003 |
| 12 | Cell composition stats | Chi-square + effect sizes for N/M/C |
| 13 | Neighborhood enrichment per stage | `sq.gr.nhood_enrichment()` × 3 stages |
| 14 | Moran's I for key genes | `sq.gr.spatial_autocorr()` |
| 15 | KI67 proliferation vs immune | Map hotspots, test exclusion |

---

## TIER 4: Future Work (AI/ML)

| # | Idea | Feasibility Blocker |
|---|------|---------------------|
| 16 | GNN cell type prediction | Training infrastructure |
| 17 | Contrastive RNA-protein learning | Novel method, high risk |
| 18 | Stage classifier from spatial graphs | Sample size (n=8) |
| 19 | scGPT fine-tuning | 357 genes may be too few |
| 20 | Pseudo-bulk → spatial deconvolution | Validation complexity |

---

## Paper Structure: "Gastric Immune Atlas"

**Title:** *Spatial immune microenvironment of gastric pre-cancer progression reveals exhausted T cell niches and emergent tertiary lymphoid structures*

### Figures

| Fig | Content | Ideas | Script |
|-----|---------|-------|--------|
| 1 | Study design + QC | 12 | `31_cell_type_annotation.py` |
| 2 | Cell composition N→M→C | 1 | `34_progression_analysis.py` |
| 3 | CAF spatial subtypes | 4 | `36_caf_subtyping.py` |
| 4 | CD8 exhaustion proximity | 2 | `35_cellcell_communication.py` |
| 5 | TLS detection | 3 | NEW: `40_tls_detection.py` |
| 6 | Immune exclusion gradients | 5 | `39_spatial_statistics.py` |

### Target Journals (ranked)

1. **Gut** - IF 24.5, gastric focus
2. **Gastroenterology** - IF 29.4, clinical impact
3. **Nature Communications** - IF 16.6, broader reach

---

## Immediate Next Steps

1. **Start with Idea #1** - progression analysis grounds everything else
2. **Run Polymath queries** before implementing each analysis
3. **Check Hwang papers** for CAF subtyping methodology
4. **Create `40_tls_detection.py`** - new script needed

### Quick Validation Commands

```bash
# Verify data availability
cd ~/spatial-hackathon-2026
python -c "import scanpy as sc; adata = sc.read('results/g4x_choi_batch2/combined_wnn.h5ad'); print(adata.obs['stage'].value_counts())"

# Verify marker availability for TLS
python -c "import scanpy as sc; adata = sc.read('results/g4x_choi_batch2/combined_wnn.h5ad'); print([g for g in ['CD20', 'CD4', 'CD8A'] if g in adata.var_names])"
```

---

## Key Polymath Queries (Pre-loaded)

```python
# Run these BEFORE implementing each analysis

# Gastric progression
mcp__polymath-v4__search_papers(query="gastric cancer metaplasia immune microenvironment progression", n=15)

# TLS detection methods
mcp__polymath-v4__search_code(query="tertiary lymphoid structure detection algorithm", n=10)

# CAF subtypes
mcp__polymath-v4__search_algorithms(query="fibroblast subtyping spatial", n=10)

# Cross-domain insights
mcp__polymath-v4__search_algorithms(bridges=True, query="immune spatial pattern")
```
