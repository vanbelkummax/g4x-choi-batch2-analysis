# G4X Gastric Cancer Progression Analysis Plan

## Executive Summary

**Dataset:** Choi_GC_preGC_batch2 (Seock-Jin Chung, VUMC)
**Platform:** G4X (Resolve Biosciences) - Multimodal Gastric
**Samples:** 32 total (8 analyzed from Lane 1, 24 remaining)
**Cells:** 522,188 (Lane 1 only)
**Modalities:** 387 genes (RNA) + 17 protein markers

**Primary Question:** What are the spatial and cellular mechanisms of gastric cancer progression from Normal → Metaplasia → Cancer?

---

## Dataset Structure

### Sample Metadata

| Lane | Pos | Block ID | ROI Type | Tissue | Status |
|------|-----|----------|----------|--------|--------|
| 1 | A | WD112575_CRC_ctrl | Control | Colon | Analyzed |
| 1 | B | SNU-103_ROI1N | Normal | Stomach | Analyzed |
| 1 | C | SNU-103_ROI2M | Metaplasia | Stomach | Analyzed |
| 1 | D | SNU-105_ROI1N | Normal | Stomach | Analyzed |
| 1 | E | SNU-105_ROI2M | Metaplasia | Stomach | Analyzed |
| 1 | F | SNU-105_ROI3C | Cancer | Stomach | Analyzed |
| 1 | G | SNU-109_ROI1C | Cancer | Stomach | Analyzed |
| 1 | H | SNU-109_ROI2X | Unknown | Stomach | Analyzed |
| 2 | A | WD112575_CRC_ctrl | Control | Colon | Pending |
| 2 | B-H | SNU-106/107 | Mixed | Stomach | Pending |
| 3 | A | WD112575_CRC_ctrl | Control | Colon | Pending |
| 3 | B-H | SNU-108/110, UTSW-4 | Mixed | Stomach | Pending |
| 4 | A | WD112575_CRC_ctrl | Control | Colon | Pending |
| 4 | B-H | Hu13/15 | Mixed | Stomach | Pending |

### ROI Type Key
- **N** = Normal gastric mucosa
- **M** = Metaplasia (pre-cancerous)
- **C** = Cancer
- **X** = Unknown/Other
- **ctrl** = CRC colon control (technical reference)

### Matched Progression Series (Priority Samples)
1. **SNU-105:** N → M → C (complete progression, Lane 1)
2. **SNU-107:** N → M → C → C (complete progression, Lane 2)
3. **SNU-103:** N → M (partial, Lane 1)

---

## Critical Flaws in Current Analysis (Must Fix)

### 1. Segmentation Error Vulnerability

**Problem:** G4X is imaging-based. Segmentation errors cause transcript "bleed-through" between adjacent cells.

**Impact (Mitchel et al. 2026):**
- 25-30% false positives in DE analysis
- Ligand-receptor inference dominated by artifacts
- "Interaction-changed genes" are often just contamination

**Required Actions:**
1. Implement cellAdmix for admixture detection
2. Add "Membrane Separation Score" validation
3. Filter "False Autocrine" loops before LR analysis
4. Run all claims through segmentation sensitivity analysis

### 2. Suboptimal Multimodal Integration

**Problem:** Current concatenated PCA treats modalities as equally weighted features.

**Why it fails:**
- Arbitrary scaling dominates results
- No per-cell modality weighting
- Loses modality-specific variance

**Required Action:** Replace with Weighted Nearest Neighbors (WNN)

### 3. Hard Gating Cell Typing

**Problem:** Boolean thresholds force binary decisions on noisy data.

**Result:** 26% "Unknown" cells, forced labels on doublets

**Required Action:** Clustering-first annotation with uncertainty scores

### 4. Cell-Level Statistics (Pseudoreplication)

**Problem:** Treating 500K correlated cells as independent samples

**Result:** Meaningless p-values, false discovery

**Required Action:** Sample-level inference (pseudobulk, mixed models)

---

## Revised Analysis Pipeline

### Phase A: Data & Segmentation QC

#### A1. Per-Cell Hard QC Gates
```python
# Filter criteria
min_genes = 10
max_genes = 300  # Adjust based on panel
min_counts = 20
max_counts = 1000
min_cell_area = 50  # um^2
max_cell_area = 500  # um^2

# Flag suspicious cells
flag_high_density = local_cell_density > threshold
flag_boundary_cells = distance_to_roi_edge < 10  # um
flag_abnormal_shape = perimeter / sqrt(area) > threshold
```

#### A2. Segmentation Sensitivity Analysis
Run ALL downstream results three ways:
1. **Baseline:** Original segmentation
2. **QC-strict:** Remove top 10% most suspicious cells
3. **Corrected:** After cellAdmix correction

**Report stability:** If claims flip across (1-3), it's pipeline sensitivity, not biology.

#### A3. cellAdmix Implementation
```r
# R package: kharchenkolab/cellAdmix
library(cellAdmix)

# Detect admixture
admix_results <- detectAdmixture(
  expr_matrix,
  cell_coords,
  neighbor_graph
)

# Flag contaminated cells
cells_to_filter <- admix_results$contamination_score > 0.3
```

#### A4. Membrane Separation Score (for LR analysis)
For any "interaction-upregulated" gene:
- Calculate membrane signal intensity between transcript location and cell centroid
- If membrane signal is high (crossing boundary), discard that molecule
- Apply before differential expression

---

### Phase B: Cell Typing (Hierarchical + Uncertainty)

#### B1. First-Pass Compartments (Protein-Driven)

| Compartment | Positive Markers | Negative Markers |
|-------------|------------------|------------------|
| Epithelial | PanCK | CD45, CD31, aSMA |
| Immune | CD45 | PanCK, CD31 |
| Endothelial | CD31 | CD45, PanCK |
| Stromal | aSMA | CD45, PanCK, CD31 |

#### B2. Second-Pass (RNA + Protein Clustering)

```python
# Within each compartment:
# 1. Leiden clustering on RNA
sc.tl.leiden(adata_compartment, resolution=0.5)

# 2. Visualize protein markers on clusters
sc.pl.dotplot(adata_compartment,
              var_names=['CD3', 'CD4', 'CD8', 'CD20', 'CD68'],
              groupby='leiden')

# 3. Assign labels with confidence
# Output: (label, confidence, runner_up_label)
```

#### B3. Uncertainty Propagation
- Never force labels on ambiguous cells
- Keep "Low Confidence" bucket explicitly
- Track per-cell classification entropy

#### B4. Cell-Bridging for "Unknown" Recovery
Use Mitchel et al. Cell-Bridging method:
- Detect admixture factors skewed toward cell boundaries
- Decompose "Unknown" cells into constituent types
- Potentially recover usable data from 26% waste

---

### Phase C: Multimodal Integration (WNN)

#### C1. Weighted Nearest Neighbors

```python
import muon as mu

# Separate modality processing
mu.pp.neighbors(mdata['rna'], use_rep='X_pca')
mu.pp.neighbors(mdata['protein'], use_rep='X_pca')

# WNN integration
mu.tl.wnn(mdata,
          modalities=['rna', 'protein'],
          n_neighbors=20)

# Cluster on WNN graph
sc.tl.leiden(mdata, neighbors_key='wnn')

# UMAP on WNN
sc.tl.umap(mdata, neighbors_key='wnn')
```

**Key advantage:** Learns per-cell modality weights rather than assuming equal contribution.

#### C2. MOFA+ for Factor Analysis (Optional)

```python
from mofapy2 import run_mofa

# Multi-view, multi-group setup
mofa_model = run_mofa(
    data={'rna': rna_matrix, 'protein': protein_matrix},
    groups=sample_ids,
    n_factors=15
)

# Interpretable factors: "immune activation", "epithelial program", "batch"
```

---

### Phase D: Spatial Analysis

#### D1. Neighborhood Analysis (Per-Sample)

```python
import squidpy as sq

# Compute per-sample, then meta-analyze
for sample in samples:
    adata_s = adata[adata.obs['sample_id'] == sample]

    # Build spatial graph
    sq.gr.spatial_neighbors(adata_s, coord_type='generic', n_neighs=10)

    # Neighborhood enrichment
    sq.gr.nhood_enrichment(adata_s, cluster_key='cell_type')

    # Store results
    results[sample] = adata_s.uns['nhood_enrichment']

# Meta-analyze across samples of same class (N/M/C)
```

#### D2. Cellular Neighborhood Analysis (CNA)

Instead of pairwise enrichment, cluster the neighborhoods themselves:

```python
# For each cell, compute neighborhood composition vector
# Cluster these vectors to find "Neighborhood Types"

# Example: "Neighborhood Type 3" = Macrophage + Fibroblast + Tumor
# Test: Does abundance of NT3 increase from M → C?
```

#### D3. Barrier Analysis (Immune Exclusion)

```python
# Distance from tumor/stroma interface
# Y-axis: Density of CD8+ T cells vs Tregs

# Hypothesis: In Cancer, Tregs accumulate at interface (blocking)
# while CD8s are excluded. In Metaplasia, barrier not formed.

def compute_barrier_plot(adata, interface_cells):
    distances = compute_distance_to_interface(adata, interface_cells)

    for cell_type in ['CD8_T', 'Treg']:
        density = kernel_density(
            distances[adata.obs['cell_type'] == cell_type]
        )
        plot(distances, density, label=cell_type)
```

---

### Phase E: Differential Analysis (Sample-Level)

#### E1. Pseudobulk Approach

```python
# Aggregate to sample level
pseudobulk = adata.obs.groupby(['sample_id', 'cell_type']).apply(
    lambda x: adata[x.index].X.sum(axis=0)
)

# DESeq2 or edgeR on pseudobulk
# Comparison: N vs M vs C (with patient as random effect)
```

#### E2. Mixed Effects Models

```python
import statsmodels.formula.api as smf

# For each gene/metric:
model = smf.mixedlm(
    "expression ~ roi_type + (1|patient_id)",
    data=df,
    groups=df['patient_id']
)
```

---

### Phase F: Ligand-Receptor Analysis (With Caution)

**Prerequisites before running:**
1. cellAdmix correction applied
2. Membrane Separation Score filter applied
3. False Autocrine loops filtered

#### F1. False Autocrine Filter

```python
# If cell type expresses both Ligand and Receptor
# AND this only occurs in high-density regions
# → Treat as segmentation artifact

def filter_false_autocrine(lr_results, cell_density):
    suspicious = (
        (lr_results['source'] == lr_results['target']) &
        (cell_density > density_threshold)
    )
    return lr_results[~suspicious]
```

#### F2. Validated LR Analysis

```python
import liana

# Only after QC filters applied
liana.mt.rank_aggregate(
    adata,
    groupby='cell_type',
    use_raw=False,
    verbose=True
)
```

---

## Priority Hypotheses

### Tier 1: Progression Mechanism (SNU-105, SNU-107)

| ID | Hypothesis | Method | Markers |
|----|------------|--------|---------|
| H1.1 | Immune infiltration pattern changes N→M→C | Cell type proportions + spatial | CD45, CD3, CD68 |
| H1.2 | Epithelial phenotype shifts during progression | Trajectory analysis (PAGA) | PanCK, RNA markers |
| H1.3 | T cell exhaustion increases in Cancer | PD1+ proportion in CD8 T cells | PD1, CD8, CD3 |
| H1.4 | Treg barrier forms at tumor interface | Barrier analysis | FOXP3, CD4, CD3 |

### Tier 2: Spatial Architecture

| ID | Hypothesis | Method | Expected |
|----|------------|--------|----------|
| H2.1 | Spatial entropy increases with progression | Shannon entropy | C > M > N |
| H2.2 | Immune-epithelial distance changes | Mean distance analysis | Changes with stage |
| H2.3 | Specific neighborhood types associate with cancer | CNA clustering | NT enrichment |

### Tier 3: Multimodal Insights

| ID | Hypothesis | Method | Impact |
|----|------------|--------|--------|
| H3.1 | RNA-protein concordance varies by cell state | Per-celltype correlation | Method validation |
| H3.2 | Protein better resolves immune subtypes | Compare RNA vs protein clustering | Modality utility |

---

## Implementation Priority

### Week 1: Foundation

| Day | Task | Output |
|-----|------|--------|
| 1 | Implement WNN integration | `scripts/30_wnn_integration.py` |
| 2 | Add segmentation QC filters | `scripts/31_segmentation_qc.py` |
| 3 | Implement hierarchical cell typing | `scripts/32_hierarchical_annotation.py` |
| 4 | Set up cellAdmix (R) | `scripts/33_celladmix_correction.R` |
| 5 | Run segmentation sensitivity | `results/sensitivity_analysis/` |

### Week 2: Analysis

| Day | Task | Output |
|-----|------|--------|
| 1 | Process Lanes 2-4 | `results/g4x_choi_batch2/all_lanes/` |
| 2 | Neighborhood analysis | `results/spatial/` |
| 3 | Barrier analysis | `results/immune_exclusion/` |
| 4 | Pseudobulk DE | `results/differential/` |
| 5 | Integrate and validate | `results/validated/` |

---

## File Locations

```
spatial-hackathon-2026/
├── G4X_ANALYSIS_PLAN.md          # This document
├── scripts/
│   ├── 20_g4x_choi_batch2_analysis.py   # Original QC (keep)
│   ├── 21_g4x_cell_annotation.py        # REPLACE with hierarchical
│   ├── 22_g4x_multimodal_integration.py # REPLACE with WNN
│   ├── 30_wnn_integration.py            # NEW: WNN
│   ├── 31_segmentation_qc.py            # NEW: QC filters
│   ├── 32_hierarchical_annotation.py    # NEW: Clustering-first
│   ├── 33_celladmix_correction.R        # NEW: Admixture correction
│   └── 34_progression_analysis.py       # NEW: N→M→C comparison
├── results/
│   └── g4x_choi_batch2/
│       ├── annotated/                   # Current (will update)
│       ├── multimodal/                  # Current (will update)
│       ├── sensitivity_analysis/        # NEW
│       ├── spatial/                     # NEW
│       └── validated/                   # NEW
└── figures/
    └── g4x/
```

---

## Key References

1. **Mitchel et al. (2026)** - Impact of Segmentation Errors in Spatial Transcriptomics
   - cellAdmix: https://github.com/kharchenkolab/cellAdmix
   - Cell-Bridging method for Unknown recovery
   - Membrane Separation Score validation

2. **Seurat WNN** - https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis
   - Per-cell modality weighting

3. **muon WNN** - https://muon-tutorials.readthedocs.io/en/latest/cite-seq/2-CITE-seq-PBMC-5k-Weighted-Neighbours.html
   - Python implementation

4. **MOFA+** - https://biofam.github.io/MOFA2/
   - Multi-view factor analysis

5. **sopa** - https://github.com/gustaveroussy/sopa
   - Technology-invariant spatial pipeline

---

## Success Criteria

1. **Segmentation Robustness:** Core claims stable across QC regimes
2. **Unknown Reduction:** <10% Unknown cells (recovered via Cell-Bridging)
3. **Biological Discovery:** At least one novel progression mechanism identified
4. **Statistical Rigor:** All p-values from sample-level inference
5. **Reproducibility:** Full pipeline documented and versioned

---

*Last Updated: 2026-01-21*
*Author: Max Van Belkum (with Claude Opus 4.5)*
