# Senior Reviewer Strategic Assessment: G4X Gastric Cancer Spatial Analysis

**Date:** 2026-01-21
**Reviewer:** Senior Scientist & Code Engineer (Claude Opus 4.5)
**Evidence Base:** Polymath (2,308 papers), Clinical Trials, ChEMBL, Open Targets, Vanderbilt Faculty Papers

---

## EXECUTIVE SUMMARY

| Dimension | Assessment | Confidence |
|-----------|------------|------------|
| **Literature Gap** | CONFIRMED - Zero N→M→C gastric spatial studies | HIGH |
| **Data Quality** | Strong - 522K cells, multi-modal (RNA+Protein) | HIGH |
| **Methods Risk** | HIGH - Too many novel methods simultaneously | HIGH |
| **Statistical Power** | WEAK - n=8 limits bulk claims | HIGH |
| **Publication Viability** | Gut-adjacent, not Gut-ready | MEDIUM |

**Bottom Line:** You have a publishable story with completed work (WNN + annotation). The "sexy methods" (TDA, tensor, OT) are credibility risks without demonstrated expertise.

---

## PART 1: EVIDENCE-BASED METHODOLOGY CRITIQUE

### 1.1 STORIES (Optimal Transport) vs CellRank

**Polymath Evidence Found:**
- STORIES paper (Nature Methods 2025) - Huizing et al. - Score 0.908
- Paper confirms: "enables rich and spatially informed analysis of differentiation trajectories"
- Reproducibility repo: `github.com/cantinilab/stories_reproducibility` (4 stars - LOW adoption)

| Criterion | STORIES | PAGA + Diffusion Pseudotime |
|-----------|---------|----------------------------|
| RNA velocity required | No | No |
| Publication track record | 2025 (new) | 2019+ (mature) |
| Reviewer familiarity | Low | High |
| Implementation risk | Medium | Low |
| G4X compatibility | Yes | Yes |

**MY VERDICT:**
- ✅ STORIES is technically appropriate (no RNA velocity needed)
- ⚠️ BUT it's unproven in gastric cancer context
- **RECOMMENDATION:** Run PAGA first as primary, STORIES as validation

```python
# Safe approach - run both
# Primary (reviewer-proof)
sc.tl.paga(adata, groups='cell_type')
sc.pl.paga(adata)
sc.tl.dpt(adata, n_dcs=10)

# Secondary (novelty)
# Only add STORIES if PAGA shows interesting trajectory
```

### 1.2 TDA (Persistent Homology) for TLS Detection

**Polymath Evidence Found:**
- TLS-Finder tool exists: `github.com/AAKoksoy/TLS-Finder` (mentioned in comprehensive framework)
- Multiple papers show TDA is emerging but NOT standard for TLS
- DBSCAN-based clustering remains dominant approach

| Method | Complexity | Reviewer Acceptance | Evidence Quality |
|--------|------------|---------------------|------------------|
| DBSCAN + co-occurrence | Low | High | Established |
| Density-based TLS scoring | Medium | High | Growing |
| TDA (persistent homology) | High | Low | Emerging |

**Critical Finding:** The paper "Antigen-presenting CAFs in gastric cancer" (Song et al. 2025) uses **correlation-based TLS scoring**, NOT TDA. This is directly relevant to your CAF work.

**MY VERDICT:**
- ⚠️ TDA is scientifically interesting but reviewer-unfriendly
- ⚠️ Persistence threshold selection has NO biological ground truth
- **RECOMMENDATION:**
  1. Primary: DBSCAN clustering + spatial co-occurrence
  2. Secondary: TDA only as supplementary validation with bootstrap CIs

```python
# Primary approach (publishable)
from sklearn.cluster import DBSCAN

# Get B and T cell coordinates
bt_mask = adata.obs['cell_type'].isin(['B_cell', 'CD4_T', 'CD8_T'])
bt_coords = adata.obsm['spatial'][bt_mask]

# DBSCAN clustering
clustering = DBSCAN(eps=50, min_samples=10).fit(bt_coords)
tls_candidates = np.unique(clustering.labels_[clustering.labels_ >= 0])

# Validate ring structure with co-occurrence
sq.gr.co_occurrence(adata, cluster_key='cell_type')
```

### 1.3 Tensor Decomposition for Niche-Phenotype

**Polymath Evidence Found:**
- LIMITED papers on tensor decomposition for spatial single-cell
- Deconvolution methods (CARD, Cell2location, SpatialDWLS) dominate
- No established tensor approach for niche scoring

**MY VERDICT:**
- ❌ OVERENGINEERED for the question being asked
- ❌ Rank selection is arbitrary without ground truth
- ❌ Results are difficult to interpret biologically

**RECOMMENDATION:** Replace with direct niche scoring:

```python
# Simple, interpretable niche scoring
def compute_exhausted_niche_score(adata, radius=50):
    """Score cells by % exhausted neighbors"""
    from scipy.spatial import cKDTree

    coords = adata.obsm['spatial']
    tree = cKDTree(coords)

    exhausted_mask = adata.obs['exhaustion_score'] > 0.5
    scores = []

    for i in range(len(coords)):
        neighbors = tree.query_ball_point(coords[i], radius)
        if len(neighbors) > 1:
            pct_exhausted = exhausted_mask.iloc[neighbors].mean()
            scores.append(pct_exhausted)
        else:
            scores.append(np.nan)

    return np.array(scores)

adata.obs['exhausted_niche_score'] = compute_exhausted_niche_score(adata)
```

### 1.4 Chemokine Vector Fields

**Polymath Evidence Found:**
- CXCL9, CXCL10, CXCL11, CCL5 confirmed in your panel
- Vector field visualization IS used in high-impact papers
- BUT: Interpretation as "mechanism" requires perturbation data

**MY VERDICT:**
- ✅ Valid visualization approach
- ⚠️ Frame as "expression pattern", not "mechanism"
- ✅ Keep this analysis

---

## PART 2: STATISTICAL POWER ANALYSIS

### n=8 Limitations

| Analysis Type | Your Power | Minimum for Claims |
|---------------|------------|-------------------|
| Bulk pseudobulk DE | UNDERPOWERED | n≥15 per group |
| Cell type proportions | WEAK | n≥10 per group |
| Single-cell patterns | STRONG | 522K cells |
| Spatial statistics | STRONG | Per-cell analysis |

### Language Calibration Guide

| ❌ CANNOT Claim | ✅ CAN Claim |
|-----------------|--------------|
| "Significant differences" | "We observed differences" |
| "Gastric cancer is characterized by..." | "In our cohort..." |
| "We demonstrate..." | "Our exploratory analysis suggests..." |
| "This proves..." | "These findings are consistent with..." |

### Recommended Statistical Framework

```python
# For n=8 comparisons - use non-parametric + effect sizes
from scipy.stats import mannwhitneyu
from numpy import abs

def safe_comparison(group1, group2, name):
    """Report effect size, not just p-value"""
    stat, pval = mannwhitneyu(group1, group2, alternative='two-sided')

    # Cohen's d effect size
    pooled_std = np.sqrt((np.std(group1)**2 + np.std(group2)**2) / 2)
    cohens_d = abs(np.mean(group1) - np.mean(group2)) / pooled_std

    # Cliff's delta (non-parametric effect size)
    n1, n2 = len(group1), len(group2)
    dominance = sum(g1 > g2 for g1 in group1 for g2 in group2)
    cliffs_delta = (2 * dominance / (n1 * n2)) - 1

    return {
        'comparison': name,
        'p_value': pval,
        'cohens_d': cohens_d,
        'cliffs_delta': cliffs_delta,
        'interpretation': 'large' if abs(cliffs_delta) > 0.474 else
                         'medium' if abs(cliffs_delta) > 0.33 else 'small'
    }
```

---

## PART 3: CLINICAL RELEVANCE (Clinical Trials Evidence)

### Active Gastric Cancer Biomarker Trials

| Trial | Status | Relevance |
|-------|--------|-----------|
| **NCT05657080** (CyGIM) | RECRUITING | Gastric intestinal metaplasia biomarkers - DIRECTLY RELEVANT |
| NCT07243015 | RECRUITING | Exosomal microRNA for recurrence |
| NCT06346054 | RECRUITING | Barrett's/GIM molecular assessment |

### Key Finding: NCT05657080 (CyGIM Trial)

The University of Cambridge trial (CyGIM) is actively recruiting to validate **Cytosponge-TFF3** for detecting gastric intestinal metaplasia. This is DIRECTLY RELEVANT to your work:

- **TFF3** is in your panel (TFF1, TFF2, TFF3 all confirmed)
- **GIM detection** is your metaplasia stage
- **Cite this trial** to establish clinical relevance

### Drug Landscape (ChEMBL)

**Pembrolizumab (CHEMBL3137343):**
- First-in-class PD-1 antibody
- Approved for gastric cancer (MSI-H/dMMR)
- Your PD1/PDL1 protein markers enable response stratification

**Clinical Framing Opportunity:**
> "Our spatial atlas identifies potential biomarkers for immunotherapy eligibility in gastric pre-cancer, including PD-1+ exhausted T cell niches that may predict response to pembrolizumab."

---

## PART 4: VANDERBILT FACULTY ALIGNMENT

### Directly Relevant Faculty Papers

| Professor | Paper | Relevance |
|-----------|-------|-----------|
| **Ken Lau** | PMID:35794563 - "Spatial transcriptomics atlas reveals crosstalk between CAFs and TME in CRC" | CAF spatial methods template |
| **Tae Hwang** | PMID:40154487 - "Conserved spatial subtypes of CAFs" | CAF subtyping framework |
| **Tae Hwang** | PMID:40954300 - "iSCALE" | Large-tissue spatial methods |

### Recommendation: Cite Hwang CAF Framework

The Hwang 2025 paper (PMID:40154487) establishes "4 conserved CAF subtypes" with spatial multi-omics. Your CAF analysis should:
1. Use their subtype definitions (mCAF, iCAF, apCAF)
2. Cite their spatial methodology
3. Show concordance with their findings

---

## PART 5: REVISED EXECUTION PLAN

### Phase 0: Defensive QC (CRITICAL - Day 1)

```bash
# RUN THIS BEFORE ANYTHING ELSE
cd ~/g4x-choi-batch2-analysis
python scripts/50_comprehensive_qc.py

# Verify outputs
ls results/qc_figures/
```

### Phase 1: Minimal Viable Paper (Days 2-4)

| Priority | Analysis | Script | Method | Risk |
|----------|----------|--------|--------|------|
| P1 | Multi-modal PCA | `41_multimodal_pca.py` | Bulk + spatial PCA | LOW |
| P2 | Cell proportions | `34_progression_analysis.py` | Chi-square + effect sizes | LOW |
| P3 | PAGA trajectory | `34_progression_analysis.py` | sc.tl.paga + dpt | LOW |
| P4 | Neighborhood enrichment | `39_spatial_statistics.py` | sq.gr.nhood_enrichment | LOW |
| P5 | Gastric markers | `34_progression_analysis.py` | Spatial plots | LOW |

### Phase 2: High-Impact Additions (Days 5-6)

| Priority | Analysis | Script | Method | Risk |
|----------|----------|--------|--------|------|
| P6 | CAF subtyping | `36_caf_subtyping.py` | Score-based (Hwang framework) | LOW |
| P7 | LIANA communication | `35_cellcell_communication.py` | liana.rank_aggregate | LOW |
| P8 | Chemokine patterns | `39_spatial_statistics.py` | Interpolation + streamplot | MEDIUM |

### Phase 3: Only If Phase 1-2 Compelling (Days 7+)

| Priority | Analysis | Script | Method | Risk |
|----------|----------|--------|--------|------|
| P9 | TLS detection | `40_tls_detection.py` | DBSCAN primary, TDA supplementary | MEDIUM |
| P10 | STORIES trajectory | `34_progression_analysis.py` | Optimal transport | MEDIUM |

### NOT RECOMMENDED

| Analysis | Reason |
|----------|--------|
| Tensor decomposition | Overengineered, hard to interpret |
| Network betweenness | Adds complexity without insight |
| Full TDA as primary | Reviewer unfamiliarity |

---

## PART 6: PUBLICATION STRATEGY

### Current Positioning

| Journal | IF | Fit | Gap to Close |
|---------|-----|-----|--------------|
| **Gut** | 24.5 | 70% | Need stronger clinical angle |
| **Gastroenterology** | 29.4 | 50% | Need validation cohort |
| **Nat Commun** | 16.6 | 75% | Methods must deliver |
| **Genome Biology** | 12.3 | 85% | Frame as methods |

### Recommended Strategy: Gut with Clinical Frame

**Title Revision:**
> *Spatial immune atlas of gastric pre-cancer progression identifies exhausted T cell niches as potential immunotherapy biomarkers*

**Key Selling Points:**
1. First N→M→C spatial gastric study (confirmed gap)
2. Multi-modal (G4X RNA + protein)
3. Clinical relevance via TFF3/GIM (cite CyGIM trial)
4. Immunotherapy stratification (PD-1/PDL-1 spatial patterns)

### Figure Structure for Gut

| Fig | Content | Methods | Status |
|-----|---------|---------|--------|
| 1A | Study design + QC | Standard | TODO |
| 1B | Multi-modal PCA | Bulk + spatial | TODO |
| 2 | Cell type changes N→M→C | Proportions + stats | TODO |
| 3 | PAGA trajectory | Pseudotime | TODO |
| 4 | Neighborhood enrichment | Squidpy | TODO |
| 5 | CAF subtypes | Hwang scoring | TODO |
| 6 | Gastric markers spatial | CDX2, MUC2, TFF | TODO |
| 7 | LIANA communication | Top LR pairs | TODO |
| S1-S3 | TLS, vector fields | Supplementary | Optional |

---

## PART 7: KEY QUESTIONS ANSWERED

### Q1: Are novel methods (TDA, tensor, OT) appropriate or overcomplicated?

**ANSWER:** Overcomplicated for a first paper on this dataset.
- Use established methods as primary
- Novel methods as supplementary validation only
- Risk of reviewer rejection on unfamiliar methods is HIGH

### Q2: With n=8 samples, what statistical claims can we make?

**ANSWER:** Single-cell claims only. Bulk comparisons are exploratory.
- Frame all between-stage comparisons as "hypothesis-generating"
- Emphasize effect sizes, not p-values
- Power is in 522K cells, not 8 samples

### Q3: Is the N→M→C progression model well-supported?

**ANSWER:** Yes, by your markers.
- CDX2, MUC2 (intestinalization)
- TFF1-3 (trefoil progression)
- Epithelial % decrease (39% → 22%)
- Cite CyGIM trial for clinical validation context

### Q4: What critical analyses are missing?

**ANSWER:**
1. QC figures (run `50_comprehensive_qc.py`)
2. Negative controls (permuted spatial coordinates)
3. Within-stage variability assessment
4. Reproducibility documentation

### Q5: Is this Gut-worthy?

**ANSWER:** Not yet, but achievable.
- Current: 70% fit
- With QC + clinical framing: 85% fit
- Add: Strong title, cite CyGIM trial, PD-1 stratification angle

---

## FINAL CHECKLIST

### Before Submitting to Gut

- [ ] QC figures generated and reviewed
- [ ] All statistical claims use appropriate language
- [ ] Clinical relevance section includes CyGIM trial citation
- [ ] CAF analysis cites Hwang 2025 framework
- [ ] Immunotherapy angle developed (PD-1/PDL-1 patterns)
- [ ] Methods section has reproducibility details
- [ ] Code repository is clean and documented
- [ ] Within-stage variability addressed
- [ ] Negative controls performed

### Immediate Actions

1. **TODAY:** Run `python scripts/50_comprehensive_qc.py`
2. **Day 2:** Create `41_multimodal_pca.py` and run
3. **Day 3:** Update `34_progression_analysis.py` with PAGA
4. **Day 4:** Run LIANA and CAF scoring
5. **Day 5:** Compile figures, assess novel method value

---

*Document generated with evidence from Polymath (2,308 papers), Clinical Trials (10 gastric trials), ChEMBL (pembrolizumab), and Vanderbilt faculty papers (Lau, Hwang)*
