# Next Session: G4X Choi Batch 2 - Graphical Analysis & Hypothesis Testing

**Repository:** https://github.com/vanbelkummax/g4x-choi-batch2-analysis
**Data:** `/mnt/x/Choi_Batch_2_Tuesday/` (392GB, 32 samples, 2.3M cells)
**Environment:** `conda activate enact`

---

## Session Objectives

### Phase 1: Polymath Knowledge Base Enrichment (REQUIRED FIRST)

Before analysis, ensure Polymath has sufficient coverage for G4X/multimodal spatial analysis.

#### 1.1 Assess Current Coverage
```bash
cd /home/user/polymath-v4
source ~/miniforge3/etc/profile.d/conda.sh && conda activate enact

# Check G4X-specific papers
python scripts/q.py "G4X spatial transcriptomics" --fast -n 20
python scripts/q.py "imaging-based spatial transcriptomics CosMx Xenium MERSCOPE" --fast -n 20
python scripts/q.py "multimodal spatial proteomics RNA" --fast -n 20

# Check immune panel analysis papers
python scripts/q.py "immune checkpoint PD1 PDL1 spatial" --fast -n 20
python scripts/q.py "tumor microenvironment spatial transcriptomics" --fast -n 20
python scripts/q.py "T cell exhaustion spatial" --fast -n 20

# Check methodology papers
python scripts/q.py "neighborhood enrichment spatial statistics" --fast -n 20
python scripts/q.py "cell type annotation protein markers" --fast -n 20

# Check algorithm registry
python scripts/algo.py "neighborhood" 2>/dev/null | head -20
python scripts/algo.py "spatial statistics" 2>/dev/null | head -20

# Check relevant repos
python scripts/q.py "squidpy scanpy spatial" --repos -n 10
```

#### 1.2 Identify & Fill Knowledge Gaps

**Papers to retrieve if missing:**
- G4X platform methodology papers
- CosMx/Xenium/MERSCOPE comparison papers
- Multimodal (protein + RNA) integration methods
- Immune checkpoint spatial analysis
- TME heterogeneity in cancer spatial studies

**Use these skills to obtain papers:**
```bash
# Use /lit skill for literature search
/lit "G4X spatial transcriptomics multimodal"
/lit "imaging spatial transcriptomics immune checkpoint"

# Use /ingest skill for paper ingestion
/ingest <paper_paths>
```

**Repos to ingest if missing:**
- `nanostring-biostats/CosMx-Analysis-Scratch-Space`
- `Nanostring-Biostats/InSituType`
- `BayraktarLab/cell2location` (if not present)
- `theislab/squidpy` (verify code chunks extracted)

#### 1.3 Generate Domain-Specific Hypotheses from Polymath

Query Polymath to generate testable hypotheses:
```bash
# Cross-domain algorithm discovery
python scripts/algo.py --spatial
python scripts/algo.py --bridges

# Neo4j graph queries for hypothesis generation
python -c "
from neo4j import GraphDatabase
driver = GraphDatabase.driver('bolt://localhost:7687', auth=('neo4j', 'polymathic2026'))
with driver.session() as session:
    # Find methods applicable to multimodal spatial data
    result = session.run('''
        MATCH (a:Algorithm)-[:BELONGS_TO]->(d:AlgorithmDomain)
        WHERE a.spatial_biology_uses IS NOT NULL
          AND (a.name CONTAINS 'multimodal' OR a.name CONTAINS 'integration')
        RETURN a.name, d.name, a.spatial_biology_uses[0]
        LIMIT 20
    ''')
    for r in result:
        print(r)
driver.close()
"
```

---

### Phase 2: Metadata Deep-Dive & Group Identification

#### 2.1 Extract All Available Metadata

**Sample-level metadata to compile:**
```python
# From core metrics files
metrics_df = pd.read_csv('/mnt/x/Choi_Batch_2_Tuesday/choi_preGC_b2_core_metrics.csv')
protein_metrics = pd.read_csv('/mnt/x/Choi_Batch_2_Tuesday/choi_preGC_b2_protein_core_metrics.csv')

# Key variables for grouping:
# - Lane (L001-L004) - technical batch
# - Sample position (A-H, 01-04) - spatial layout on chip
# - Cell count quartiles
# - Transcript detection quartiles
# - Tissue area
# - Empty cell percentage
```

**Check for clinical/experimental metadata:**
```bash
# Look for any annotation files
find /mnt/x/Choi_Batch_2_Tuesday/ -name "*.xlsx" -o -name "*annotation*" -o -name "*clinical*" -o -name "*sample_info*"

# Check QuPath projects for ROI annotations
ls -la /mnt/x/Choi_Batch_2_Tuesday/Qupath_HE/
ls -la /mnt/x/Choi_Batch_2_Tuesday/Qupath_FL/
```

#### 2.2 Define Comparison Groups

Based on current analysis, potential groupings:

| Grouping Variable | Groups | Rationale |
|-------------------|--------|-----------|
| **Lane** | L001, L002, L003, L004 | Technical batch effects |
| **CD8/Treg Ratio** | High (>5), Low (<5) | Immune favorability |
| **PDL1+ Tumor** | High (>30%), Low (<30%) | Checkpoint status |
| **Immune Infiltration** | High (>25%), Low (<15%) | Hot vs cold |
| **Exhaustion Score** | High, Medium, Low | T cell functionality |
| **Cell Density** | High, Medium, Low | Tissue cellularity |

---

### Phase 3: Per-Sample Graphical Analysis

Generate comprehensive visual profiles for EACH of the 32 samples.

#### 3.1 Sample Profile Dashboard (per sample)

For each sample, generate a multi-panel figure including:

1. **Spatial Map** - Cell positions colored by cell type
2. **H&E Context** - Registered H&E image (from QuPath if available)
3. **Protein Expression Heatmap** - All 17 markers across cells
4. **Cell Type Pie Chart** - Composition breakdown
5. **Marker Violin Plots** - Distribution of key markers
6. **Neighborhood Graph** - Local cell-cell interaction network
7. **Spatial Density Maps** - KDE of immune cell distribution
8. **Checkpoint Co-expression** - PD1 vs PDL1 scatter

#### 3.2 Implementation Script Structure

```python
# scripts/06_sample_profiles.py
def generate_sample_profile(sample_id):
    """Generate comprehensive visual profile for one sample."""

    # Load data
    df = load_sample_data(sample_id)

    # Create 4x2 figure
    fig = plt.figure(figsize=(20, 24))

    # Panel 1: Spatial map by cell type
    # Panel 2: H&E image overlay
    # Panel 3: Protein heatmap
    # Panel 4: Cell type composition
    # Panel 5: Marker distributions
    # Panel 6: Spatial density (KDE)
    # Panel 7: Neighborhood network
    # Panel 8: Checkpoint scatter

    # Save
    fig.savefig(f'results/sample_profiles/{sample_id}_profile.png')
```

---

### Phase 4: Between-Group Comparisons

#### 4.1 Statistical Comparisons

For each grouping variable, compare:

| Metric | Test | Visualization |
|--------|------|---------------|
| Cell type proportions | Mann-Whitney U / Kruskal-Wallis | Box plots |
| Marker expression | t-test / Wilcoxon | Violin plots |
| Spatial statistics | Permutation test | Bar plots with CI |
| Neighborhood enrichment | Bootstrap | Heatmap difference |

#### 4.2 Specific Comparisons to Run

```python
# 1. Lane effects (technical batch)
compare_groups(df, group_var='lane', metrics=['all'])

# 2. Immune hot vs cold
hot_samples = df[df['total_immune_pct'] > 25]['sample_id']
cold_samples = df[df['total_immune_pct'] < 15]['sample_id']
compare_groups(hot_samples, cold_samples, metrics=['spatial', 'markers'])

# 3. Checkpoint high vs low
pdl1_high = df[df['pct_pdl1_pos_tumor'] > 30]['sample_id']
pdl1_low = df[df['pct_pdl1_pos_tumor'] < 20]['sample_id']
compare_groups(pdl1_high, pdl1_low, metrics=['immune', 'exhaustion'])

# 4. CD8/Treg ratio comparison
favorable = df[df['cd8_treg_ratio'] > 5]['sample_id']
unfavorable = df[df['cd8_treg_ratio'] < 2]['sample_id']
compare_groups(favorable, unfavorable, metrics=['spatial', 'markers'])
```

---

### Phase 5: Hypothesis Generation & Testing

#### 5.1 Hypotheses from Current Analysis

| # | Hypothesis | Test | Expected |
|---|------------|------|----------|
| H1 | High CD8/Treg samples have distinct spatial organization | Neighborhood enrichment comparison | Different immune clustering |
| H2 | PDL1+ tumors show increased immune exclusion | Spatial proximity analysis | Greater epithelial-immune distance |
| H3 | Lane effects are minimal | ANOVA across lanes | No significant differences |
| H4 | Exhausted T cells localize near PDL1+ cells | Co-localization analysis | Positive correlation |
| H5 | High immune samples have more TLS-like structures | B cell + T cell cluster detection | More aggregates |

#### 5.2 Polymath-Derived Hypotheses

Query Polymath for novel testable hypotheses:
```bash
# Find methods that could reveal new biology
python scripts/q.py "tertiary lymphoid structure spatial" --fast
python scripts/q.py "immune exclusion tumor spatial" --fast
python scripts/q.py "macrophage polarization spatial" --fast

# Algorithm-based hypotheses
python scripts/algo.py "clustering spatial"
python scripts/algo.py "graph community detection"
```

#### 5.3 Hypothesis Testing Framework

```python
# scripts/07_hypothesis_testing.py

def test_hypothesis(h_id, samples_a, samples_b, metric, test='mannwhitneyu'):
    """
    Test a specific hypothesis with proper statistics.

    Returns:
    - Effect size (Cohen's d or rank-biserial)
    - P-value (with multiple testing correction)
    - Confidence intervals (bootstrap)
    - Visualization
    """
    pass
```

---

### Phase 6: Deliverables

#### 6.1 New Figures to Generate

| Figure | Description |
|--------|-------------|
| `sample_profiles/*.png` | 32 individual sample profiles |
| `group_comparisons.png` | Between-group statistical comparisons |
| `hypothesis_results.png` | Hypothesis testing summary |
| `lane_batch_effects.png` | Technical batch analysis |
| `immune_spatial_patterns.png` | Immune organization by group |
| `tls_detection.png` | Tertiary lymphoid structure analysis |

#### 6.2 New Tables

| Table | Contents |
|-------|----------|
| `group_statistics.csv` | All pairwise group comparisons |
| `hypothesis_results.csv` | Hypothesis test results with p-values |
| `sample_metadata_complete.csv` | All extracted metadata |
| `batch_effect_analysis.csv` | Lane/position effects |

---

## Quick Start Commands

```bash
# 1. Activate environment
cd /home/user/g4x-choi-batch2-analysis
source ~/miniforge3/etc/profile.d/conda.sh
conda activate enact

# 2. Check Polymath coverage first
cd /home/user/polymath-v4
python scripts/q.py "G4X multimodal spatial" --fast -n 10

# 3. Return to analysis
cd /home/user/g4x-choi-batch2-analysis

# 4. Load previous results
python -c "
import pandas as pd
hyp = pd.read_csv('results/tables/hypothesis_analysis_results.csv')
print(hyp[['sample_id', 'total_immune_pct', 'cd8_treg_ratio', 'pct_pdl1_pos_tumor']].head(10))
"
```

---

## Current Analysis Results Summary

From completed analysis (`~/g4x-choi-batch2-analysis/`):

| Metric | Value |
|--------|-------|
| Total samples | 32 |
| Total cells | 2,308,968 |
| Mean immune % | 20.0% |
| CD8/Treg ratio | 7.3 ± 9.5 |
| PD1+ T cells | 29.0% ± 13.2% |
| PDL1+ tumor | 25.3% ± 9.0% |
| RNA-protein corr | r = 0.088 |

**Notable samples:**
- High CD8/Treg: B01 (39.3), H03 (38.4), G03 (25.0)
- Low CD8/Treg: A03 (0.9), A04 (1.3), A01 (1.5)
- Low transcript: F01 (34/cell), H04 (25/cell)

---

## Files Reference

| File | Location |
|------|----------|
| Analysis scripts | `~/g4x-choi-batch2-analysis/scripts/` |
| Results tables | `~/g4x-choi-batch2-analysis/results/tables/` |
| Figures | `~/g4x-choi-batch2-analysis/results/figures/` |
| Raw data | `/mnt/x/Choi_Batch_2_Tuesday/` |
| Polymath v4 | `/home/user/polymath-v4/` |

---

## Priority Order

1. **Polymath enrichment** - Ensure sufficient G4X/multimodal/immune papers
2. **Metadata extraction** - Compile all available sample annotations
3. **Group definitions** - Define biologically meaningful comparison groups
4. **Sample profiles** - Generate per-sample visual dashboards
5. **Group comparisons** - Statistical tests between groups
6. **Hypothesis testing** - Test specific biological hypotheses
7. **GitHub update** - Push new analysis to repo
