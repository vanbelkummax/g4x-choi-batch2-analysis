# Deep Research Review Request: G4X Gastric Cancer Spatial Analysis

**GitHub:** https://github.com/vanbelkummax/g4x-choi-batch2-analysis
**Local:** `~/g4x-choi-batch2-analysis/` and `~/spatial-hackathon-2026/results/g4x_choi_batch2/`

---

## PROJECT OVERVIEW

**Dataset:** 522K cells, 8 samples, 357 RNA + 17 protein markers (G4X/Resolve platform)
**Goal:** First spatial atlas of gastric pre-cancer progression (Normal → Metaplasia → Cancer)
**Gap:** Zero N→M→C spatial gastric studies exist (validated via Polymath search of 2,300 papers)

---

## WHAT WE HAVE COMPLETED

| Analysis | Script | Status | Output |
|----------|--------|--------|--------|
| QC & Preprocessing | `01_qc_analysis.py` | ✅ Done | 8 samples loaded |
| Cell Type Annotation | `02_cell_type_annotation.py` | ✅ Done | 7 cell types |
| WNN Integration | `03_multimodal_integration.py` | ✅ Done | RNA+Protein combined |
| Basic Spatial | `04_spatial_analysis.py` | ✅ Done | Neighborhood enrichment |
| Marker Validation | (inline) | ✅ Done | 83-100% marker coverage confirmed |

**Data Files:**
- `results/g4x_choi_batch2/annotated_v2/{A-H}01_annotated.h5ad` - Per-sample annotated

---

## WHAT WE HAVE NOT YET ANALYZED (TIER 1 PRIORITIES)

| Analysis | Script | Method | Status |
|----------|--------|--------|--------|
| **Multi-Modal PCA** | `41_multimodal_pca.py` | Bulk + Spatial PCA trajectories | ❌ TODO |
| **Progression Streamlines** | `34_progression_analysis.py` | STORIES (Optimal Transport) | ❌ TODO |
| **TLS Detection** | `40_tls_detection.py` | Persistent Homology (TDA) | ❌ TODO |
| **CD8 Exhaustion Niche** | `35_exhaustion_tensor.py` | Tensor Decomposition | ❌ TODO |
| **CAF Network Topology** | `36_caf_subtyping.py` | Delaunay + Betweenness | ❌ TODO |
| **Immune Exclusion** | `39_spatial_statistics.py` | Chemokine Vector Fields | ❌ TODO |

---

## OUR PLANNED METHODS (CRITIQUE THESE)

### 1. Progression Analysis - STORIES vs CellRank
**Plan:** Use STORIES (Optimal Transport, Nature Methods 2025) instead of CellRank
**Rationale:** G4X lacks splicing info required for RNA velocity
**Question for reviewer:** Is this the right choice? Alternatives?

### 2. TLS Detection - Topological Data Analysis
**Plan:** Persistent homology (H1) to detect ring-shaped B-T cell aggregates
**Tools:** gudhi, ripser, persim
**Question:** What persistence threshold? How to validate TLS vs random aggregates?

### 3. Niche-Phenotype Tensor
**Plan:** Non-negative PARAFAC decomposition of cells × phenotypes × neighbors
**Question:** How to select tensor rank? Is this overcomplicated vs simpler niche scoring?

### 4. Chemokine Vector Fields
**Plan:** Interpolate CXCL9/10/11 gradients, visualize as streamplots
**Markers confirmed:** CXCL9, CXCL10, CXCL11, CCL5 present
**Question:** Is gradient interpolation valid at single-cell resolution?

---

## KEY DOCUMENTS TO REVIEW

1. **Strategy:** `docs/STRATEGIC_COMMAND_BRIEFING.md` - Full analysis plan
2. **Ideas Ranked:** `docs/PROJECT_IDEAS_RANKED.md` - Tier 1-4 with scores
3. **QC Checklist:** `docs/QC_DETAILED_CHECKLIST.md` - 50+ QC items
4. **QC Script:** `scripts/50_comprehensive_qc.py` - Figure generation

---

## QUESTIONS FOR REVIEWERS

1. **Methods:** Are our novel methods (TDA, tensor, OT) appropriate or overcomplicated?
2. **Statistics:** With n=8 samples, what statistical claims can we make?
3. **Biological:** Is the N→M→C progression model well-supported by our markers?
4. **Gaps:** What critical analyses are we missing?
5. **Publication:** Is this Gut-worthy (IF 24.5) or should we aim elsewhere?

---

## HOW TO CRITIQUE

```bash
# Clone and explore
git clone https://github.com/vanbelkummax/g4x-choi-batch2-analysis
cd g4x-choi-batch2-analysis

# Key files
cat docs/STRATEGIC_COMMAND_BRIEFING.md
cat docs/PROJECT_IDEAS_RANKED.md
ls scripts/*.py
```

**Please provide:**
1. Red team critique of methodology
2. Missing analyses we should add
3. Statistical/biological concerns
4. Publication strategy feedback
