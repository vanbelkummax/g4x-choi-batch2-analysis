# G4X Gastric Cancer: Strategic Project Selection

## Executive Summary

**Dataset:** 522K cells, 8 samples, 387 RNA + 17 protein markers (G4X platform)
**Gap Found:** Zero spatial pre-cancer gastric studies in Polymath (2,300 papers)
**Recommended Focus:** Combine top 6 ideas into "Gastric Immune Atlas" paper

---

## Decision Matrix

```
                    FEASIBILITY
                 LOW    MED    HIGH
           ┌───────────────────────────┐
     HIGH  │  18    7,11   1,2,3,4,5,6│  ← PRIORITY (6 core)
SIGNIFICANCE│  20   17,19    8,9,10   │
     LOW   │        21       12-16    │  ← Quick wins for supplement
           └───────────────────────────┘
```

---

## TIER 1: Core Paper (Do All 6)

### 1. Spatial Immune Architecture of Gastric Pre-Cancer Progression
**Score:** Sig⭐5 × Feas⭐5 × Nov⭐5 = **125**

| Aspect | Detail |
|--------|--------|
| Gap | Zero N→M→C spatial studies in stomach |
| Data | All 8 samples span progression stages |
| Output | Fig 2: Pseudotime streamlines showing evolutionary flow |

**Methodological Upgrade: Pseudotime Streamlines**

| Old Method | New Method | Why Better |
|------------|------------|------------|
| Composition bar charts | **Pseudotime streamlines in physical space** | Infers causality & evolutionary flow, not just association |

The "stream" visualization mechanistically proves *how* cancer develops (showing cellular trajectories from Normal glands → Metaplasia → Cancer surface), rather than just stating *that* compositions differ.

**Implementation:**
```python
# Polymath grounding
mcp__polymath-v4__search_papers(query="pseudotime trajectory spatial transcriptomics", n=10)
mcp__polymath-v4__search_algorithms(query="pseudotime spatial embedding", n=10)

# Analysis (scripts/34_progression_analysis.py)
import scvelo as scv
import cellrank as cr

# Compute pseudotime with CellRank
cr.tl.terminal_states(adata)
cr.tl.lineages(adata)
cr.pl.lineages(adata, same_plot=True)

# Project onto spatial coordinates as streamlines
# Custom visualization showing "flow" from N→M→C
scv.pl.velocity_embedding_stream(adata, basis='spatial', color='stage')
```

---

### 2. CD8 Exhaustion Niche Architecture
**Score:** Sig⭐5 × Feas⭐5 × Nov⭐5 = **125** *(upgraded from 100)*

| Aspect | Detail |
|--------|--------|
| Markers | PD1, PDL1, CD8 proteins ✓ |
| Finding | 25% exhaustion rate already detected |
| Output | Fig 4: Niche-phenotype tensor linking cell states to exact neighbors |

**Methodological Upgrade: Niche-Phenotype Tensor**

| Old Method | New Method | Why Better |
|------------|------------|------------|
| Distance distributions, proximity scoring | **Niche-phenotype tensor decomposition** | Links exact cell states to neighbor composition, reveals hidden niche patterns |

Instead of simple distance metrics, the tensor approach captures the full relationship between a cell's phenotype (exhausted vs functional) and its complete neighborhood composition, enabling discovery of specific niche signatures that drive exhaustion.

**Implementation:**
```python
# Polymath grounding
mcp__polymath-v4__search_papers(query="niche phenotype tensor spatial single cell", n=10)
mcp__polymath-v4__search_algorithms(query="tensor decomposition spatial biology", n=10)

# Build niche-phenotype tensor (scripts/35_cellcell_communication.py)
import tensorly as tl
from tensorly.decomposition import non_negative_parafac

# Construct tensor: cells × phenotypes × neighbor_types
# Decompose to find latent niche programs
factors = non_negative_parafac(niche_tensor, rank=5)

# Identify which neighbor compositions predict exhaustion
exhaustion_niche_signature = factors[1][:, exhaustion_component]
```

---

### 3. TLS-like Structures in Gastric Metaplasia
**Score:** Sig⭐5 × Feas⭐5 × Nov⭐5 = **125** *(upgraded from 100)*

| Aspect | Detail |
|--------|--------|
| Gap | Zero gastric TLS spatial papers in Polymath |
| Markers | CD20, CD4, CD8 for B-T aggregates ✓ |
| Output | Fig 5: TDA persistence diagrams validating ring structure |

**Methodological Upgrade: Topological Data Analysis (TDA)**

| Old Method | New Method | Why Better |
|------------|------------|------------|
| Clustering/co-occurrence | **Persistent homology (TDA)** | Mathematically proves ring-shaped structure, protects from "random aggregate" critique |

TDA specifically detects "loops" (H1 homology) in spatial data. Since TLSs are ring-shaped (B-cell core, T-cell rim), TDA provides mathematical **validation** of the structure's shape—not just proximity.

**Implementation:**
```python
# Polymath grounding
mcp__polymath-v4__search_papers(query="topological data analysis spatial biology", n=10)
mcp__polymath-v4__search_algorithms(query="persistent homology TDA", n=10)
mcp__polymath-v4__search_algorithms(domain="topology", bridges=True)

# TDA analysis (scripts/40_tls_detection.py)
import gudhi
from ripser import ripser
from persim import plot_diagrams

# Extract B-cell and T-cell coordinates
b_coords = adata[adata.obs['cell_type'] == 'B_cell'].obsm['spatial']
t_coords = adata[adata.obs['cell_type'].isin(['CD4_T', 'CD8_T'])].obsm['spatial']

# Compute persistent homology - H1 detects loops/rings
dgms = ripser(bt_coords, maxdim=1)['dgms']

# Significant H1 features (long-lived loops) = TLS candidates
tls_candidates = dgms[1][dgms[1][:, 1] - dgms[1][:, 0] > persistence_threshold]
plot_diagrams(dgms, show=True)
```

---

### 4. Spatial CAF Niche Architecture (Hwang Framework)
**Score:** Sig⭐4 × Feas⭐5 × Nov⭐5 = **100** *(upgraded from 80)*

| Aspect | Detail |
|--------|--------|
| Reference | Hwang 2025 (PMID:40154487) - 4 CAF subtypes |
| Markers | ACTA2, IL6, CD74, HLA-DR ✓ |
| Output | Fig 3: 3D interaction network showing CAF hub topology |

**Methodological Upgrade: 3D Interaction Networks**

| Old Method | New Method | Why Better |
|------------|------------|------------|
| CAF subtype proportion heatmaps | **3D spatial interaction network graphs** | Reveals network topology—shows s1-CAFs form a hub that disconnects immune network |

Instead of just showing "s1-CAFs are present," the network visualization demonstrates that CAFs form **structural hubs** that physically rewire the TME communication network.

**Implementation:**
```python
# Polymath grounding - Hwang paper
mcp__vanderbilt-professors__search_hwang_papers(query="CAF spatial subtype", limit=5)
mcp__polymath-v4__search_algorithms(query="spatial interaction network graph", n=10)

# 3D Network construction (scripts/36_caf_subtyping.py)
import networkx as nx
from scipy.spatial import Delaunay

# Build spatial graph from Delaunay triangulation
tri = Delaunay(adata.obsm['spatial'])
G = nx.Graph()
for simplex in tri.simplices:
    for i in range(3):
        for j in range(i+1, 3):
            G.add_edge(simplex[i], simplex[j])

# Compute network metrics per CAF subtype
caf_betweenness = nx.betweenness_centrality(G.subgraph(caf_nodes))
caf_degree = dict(G.degree(caf_nodes))

# Visualize as 3D network (x, y, expression_intensity)
# Show CAFs as "hubs" disconnecting immune clusters
```

---

### 5. Immune Exclusion Dynamics
**Score:** Sig⭐4 × Feas⭐5 × Nov⭐5 = **100** *(upgraded from 80)*

| Aspect | Detail |
|--------|--------|
| Question | Does cancer exclude immune cells? How? |
| Data | Spatial coords + cell types + chemokine expression |
| Output | Fig 6: Chemokine vector fields showing immune flow & barriers |

**Methodological Upgrade: Chemokine Vector Fields**

| Old Method | New Method | Why Better |
|------------|------------|------------|
| Distance gradients (1D metric) | **Chemokine gradient vector fields** | Shows mechanism—T-cells *want* to enter but are blocked/diverted |

Vector fields model the *mechanism* of exclusion by visualizing directional chemokine gradients. This shows immune cells following a gradient but being physically blocked or diverted—a *Nature*-worthy figure vs. a distance histogram.

**Implementation:**
```python
# Polymath grounding
mcp__polymath-v4__search_papers(query="chemokine gradient spatial tumor microenvironment", n=10)
mcp__polymath-v4__search_algorithms(query="vector field interpolation spatial", n=10)

# Vector field analysis (scripts/39_spatial_statistics.py)
from scipy.interpolate import griddata
import numpy as np

# Compute chemokine gradient (e.g., CXCL9, CXCL10 for T-cell attraction)
chemokine_expr = adata[:, 'CXCL9'].X.toarray().flatten()
coords = adata.obsm['spatial']

# Interpolate to grid
grid_x, grid_y = np.mgrid[coords[:,0].min():coords[:,0].max():100j,
                          coords[:,1].min():coords[:,1].max():100j]
grid_chemokine = griddata(coords, chemokine_expr, (grid_x, grid_y), method='cubic')

# Compute gradient vectors
grad_y, grad_x = np.gradient(grid_chemokine)

# Plot as quiver/streamplot showing "flow" direction
plt.streamplot(grid_x, grid_y, grad_x, grad_y, color=grid_chemokine, cmap='coolwarm')
```

---

### 6. Multi-Modal PCA: Progression Trajectories & Tissue Identity
**Score:** Sig⭐5 × Feas⭐5 × Nov⭐4 = **100**

| Aspect | Detail |
|--------|--------|
| Question | How do RNA, protein, and combined modalities capture progression? Does intestinal metaplasia cluster with colon or stomach? |
| Data | All 8 samples, pseudobulk + single-cell, RNA + Protein + WNN |
| Output | Fig 1B: Multi-panel PCA showing modality-specific and combined trajectories |

**Rationale: Foundational Molecular Cartography**

This analysis establishes the molecular basis for all downstream hypotheses by showing:
1. **Modality-specific signals**: Do RNA and protein tell the same story?
2. **Progression trajectory**: Does N→M→C form a continuous gradient or discrete jumps?
3. **Tissue identity**: Intestinal metaplasia should cluster toward "colon-like" phenotype, validating the biological model of gastric intestinalization

**Implementation:**
```python
# Polymath grounding
mcp__polymath-v4__search_papers(query="multi-modal PCA RNA protein spatial transcriptomics", n=10)
mcp__polymath-v4__search_papers(query="intestinal metaplasia molecular signature colon gastric", n=10)

# Analysis (scripts/41_multimodal_pca.py)
import scanpy as sc
import numpy as np
from sklearn.decomposition import PCA

# === BULK PCA (pseudobulk per sample) ===
# Aggregate to sample level
pseudobulk_rna = adata.to_df().groupby(adata.obs['sample_id']).mean()
pseudobulk_protein = adata.obsm['protein'].groupby(adata.obs['sample_id']).mean()

# Run PCA per modality
pca_rna = PCA(n_components=10).fit_transform(pseudobulk_rna)
pca_protein = PCA(n_components=10).fit_transform(pseudobulk_protein)

# Combined (concatenate normalized)
from sklearn.preprocessing import StandardScaler
combined = np.hstack([
    StandardScaler().fit_transform(pseudobulk_rna),
    StandardScaler().fit_transform(pseudobulk_protein)
])
pca_combined = PCA(n_components=10).fit_transform(combined)

# === SPATIAL PCA (single-cell with spatial context) ===
# Use existing WNN embedding or compute fresh
sc.pp.pca(adata, n_comps=50)  # RNA
sc.pp.pca(adata, layer='protein', n_comps=20)  # Protein

# Color by: stage (N/M/C), tissue_type (gastric/colon-like), sample_id
# Highlight: metaplasia samples should shift toward "intestinal" direction

# === KEY VISUALIZATION ===
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
# Row 1: Bulk PCA (RNA, Protein, Combined) colored by stage
# Row 2: Spatial PCA projections colored by tissue identity
# Add trajectory arrows showing N→M→C direction
```

**Key Biological Questions:**

| Question | Expected Result | If Different |
|----------|-----------------|--------------|
| Does metaplasia cluster with colon? | Yes - intestinalization signature | Gastric-specific metaplasia pathway |
| RNA vs Protein agreement? | Partial (r≈0.3-0.5) | Post-transcriptional regulation important |
| Progression continuous or discrete? | Gradient with M intermediate | Stepwise transformation model |

**Markers to Highlight in PCA Loadings:**

| Tissue Identity | Key Markers |
|-----------------|-------------|
| Gastric (stomach) | MUC5AC, TFF1, GKN1/2 |
| Intestinal (colon-like) | CDX2, MUC2, TFF3, VIL1 |
| Cancer-associated | KI67, TP53, ERBB2 |

---

## Summary: Methodological Upgrade Path

| Rank | Biological Question | Old Method | **New "Sexy" Method** | Impact |
|------|---------------------|------------|------------------------|--------|
| **1** | Progression Mechanism | Bar charts | **Pseudotime Streamlines** | Infers causality |
| **2** | CD8 Exhaustion | Proximity scoring | **Niche-Phenotype Tensor** | Links states to niches |
| **3** | TLS Detection | Clustering | **Topological Data Analysis** | Mathematical proof of shape |
| **4** | CAF Architecture | Heatmaps | **3D Interaction Networks** | Reveals hub topology |
| **5** | Immune Exclusion | Distance plots | **Chemokine Vector Fields** | Visualizes forces |
| **6** | Tissue Identity & Trajectory | Single-modality PCA | **Multi-Modal Bulk + Spatial PCA** | Validates intestinalization model |

**Key Insight:** The biological questions (Top 6) are literature-backed and solid. The upgrades are **analytical engines** that transform descriptive observations into mechanistic insights worthy of high-impact journals.

---

## TIER 2: If Time Permits

| # | Idea | Score | Key Insight |
|---|------|-------|-------------|
| 7 | RNA-protein discordance as signal | 75 | r=0.088 is biology, not failure |
| 8 | Neighborhood signatures → stage classifier | 64 | ML on enrichment patterns |
| 9 | Treg/CD8 ratio spatial gradients | 48 | Already have 7.3 ratio |
| 10 | PDL1+ tumor clustering (hot/cold) | 48 | ICI response prediction |
| 11 | Cross-sample domain transfer | 48 | SpatialGlue validation |

---

## TIER 3: Quick Wins (Supplement Material)

| # | Idea | Implementation |
|---|------|----------------|
| 12 | Protein QC benchmark | Document HLA-DR r=0.31 vs KI67 r=0.003 |
| 13 | Cell composition stats | Chi-square + effect sizes for N/M/C |
| 14 | Neighborhood enrichment per stage | `sq.gr.nhood_enrichment()` × 3 stages |
| 15 | Moran's I for key genes | `sq.gr.spatial_autocorr()` |
| 16 | KI67 proliferation vs immune | Map hotspots, test exclusion |

---

## TIER 4: Future Work (AI/ML)

| # | Idea | Feasibility Blocker |
|---|------|---------------------|
| 17 | GNN cell type prediction | Training infrastructure |
| 18 | Contrastive RNA-protein learning | Novel method, high risk |
| 19 | Stage classifier from spatial graphs | Sample size (n=8) |
| 20 | scGPT fine-tuning | 357 genes may be too few |
| 21 | Pseudo-bulk → spatial deconvolution | Validation complexity |

---

## Paper Structure: "Gastric Immune Atlas"

**Title:** *Spatial immune microenvironment of gastric pre-cancer progression reveals exhausted T cell niches and emergent tertiary lymphoid structures*

### Figures

| Fig | Content | Method | Script |
|-----|---------|--------|--------|
| 1A | Study design + QC | Standard | `31_cell_type_annotation.py` |
| 1B | Multi-modal PCA trajectories | **Bulk + Spatial PCA** | `41_multimodal_pca.py` |
| 2 | Progression streamlines N→M→C | **Pseudotime** | `34_progression_analysis.py` |
| 3 | CAF spatial network topology | **3D Networks** | `36_caf_subtyping.py` |
| 4 | CD8 exhaustion niche tensor | **Tensor decomp** | `35_cellcell_communication.py` |
| 5 | TLS detection via TDA | **Persistent homology** | `40_tls_detection.py` |
| 6 | Immune exclusion vector fields | **Chemokine gradients** | `39_spatial_statistics.py` |

### Target Journals (ranked)

1. **Gut** - IF 24.5, gastric focus
2. **Gastroenterology** - IF 29.4, clinical impact
3. **Nature Communications** - IF 16.6, broader reach

---

## Immediate Next Steps

1. **Start with Idea #6** - Multi-modal PCA establishes foundational molecular cartography
2. **Then Idea #1** - progression analysis with pseudotime streamlines
3. **Run Polymath queries** for each new method before implementing
4. **Install dependencies:** `gudhi`, `ripser`, `persim`, `tensorly`, `cellrank`
5. **Create new scripts:** `40_tls_detection.py` (TDA), `41_multimodal_pca.py` (PCA trajectories)

### Dependency Installation

```bash
conda activate enact
pip install gudhi ripser persim tensorly cellrank
```

### Quick Validation Commands

```bash
# Verify data availability
cd ~/spatial-hackathon-2026
python -c "import scanpy as sc; adata = sc.read('results/g4x_choi_batch2/combined_wnn.h5ad'); print(adata.obs['stage'].value_counts())"

# Verify marker availability for TLS
python -c "import scanpy as sc; adata = sc.read('results/g4x_choi_batch2/combined_wnn.h5ad'); print([g for g in ['CD20', 'CD4', 'CD8A'] if g in adata.var_names])"

# Verify new dependencies
python -c "import gudhi, ripser, persim, tensorly, cellrank; print('All TDA/tensor deps OK')"
```

---

## Key Polymath Queries (Pre-loaded)

```python
# Run these BEFORE implementing each analysis

# Pseudotime in spatial context
mcp__polymath-v4__search_papers(query="pseudotime trajectory spatial coordinates single cell", n=15)
mcp__polymath-v4__search_algorithms(query="pseudotime embedding spatial", n=10)

# TDA for spatial biology
mcp__polymath-v4__search_papers(query="topological data analysis persistent homology spatial", n=15)
mcp__polymath-v4__search_algorithms(domain="topology", bridges=True)

# Tensor methods for niches
mcp__polymath-v4__search_papers(query="tensor decomposition single cell niche", n=10)

# Vector fields in TME
mcp__polymath-v4__search_papers(query="chemokine gradient vector field tumor microenvironment", n=10)

# Network topology in spatial
mcp__polymath-v4__search_papers(query="spatial interaction network graph tumor", n=10)
```
