# G4X Hackathon Priority Analysis Plan

**Date:** 2026-01-22
**Dataset:** G4X Choi Batch 2 - Gastric Cancer Progression (1.83M cells, 29 samples)
**Goal:** Generate publishable results leveraging the within-patient paired design

---

## Strategic Assessment

### Why Current Approaches Failed
| Approach | Result | Reason |
|----------|--------|--------|
| Pseudobulk DE | 0 significant genes | N=29 underpowered, stage-patient confounded |
| PCA stage separation | p=0.578 | Cell type dominates (19.2%), stage only 1.3% |
| Protein-RNA correlation | r=0.095 | Expected for imaging platforms |

### The Opportunity
**SNU-105 and SNU-107 each have complete Normal → Metaplasia → Cancer progression from the SAME patient.** This paired design eliminates patient confounding and enables powerful within-patient comparisons.

| Patient | Normal | Metaplasia | Cancer | Total |
|---------|--------|------------|--------|-------|
| SNU-105 | 189,590 | 127,313 | 241,149 | 558,052 |
| SNU-107 | 198,875 | 220,549 | 307,858 | 727,282 |

---

## Implementation Plan (3 Scripts)

### Script 1: `49_within_patient_composition.py` [HIGHEST PRIORITY]

**Purpose:** Detect cell type proportion shifts during N→M→C progression within each patient.

**Key Functions:**
```python
def subset_to_paired_patients(adata, patients=['SNU-105', 'SNU-107'])
def compute_cell_counts(adata, groupby=['patient', 'stage', 'cell_type'])
def chi_square_composition_test(counts_df, patient, stage1, stage2)
def fisher_exact_per_celltype(counts_df, patient, stage1, stage2, min_cells=10)
def run_paired_composition_analysis(adata) -> Dict with all results
```

**Statistical Methods:**
- Chi-square test for overall composition difference
- Fisher's exact test per cell type (2x2: this_type vs others)
- BH-FDR correction
- Stouffer's method to combine p-values across patients

**Outputs:**
```
results/within_patient_composition/
├── figures/
│   ├── fig1_stacked_bar_composition.png
│   ├── fig2_alluvial_SNU105.png
│   ├── fig2_alluvial_SNU107.png
│   ├── fig3_fold_change_heatmap.png
│   └── fig4_proportion_trajectories.png
├── cell_counts.csv
├── proportions.csv
├── chi_square_results.csv
├── fisher_per_celltype.csv
├── consistent_changes.csv
└── COMPOSITION_REPORT.md
```

**Expected Result:** "Fibroblasts increase 2.3x from N→C in both patients (p<0.001)"

---

### Script 2: `50_protein_spatial_neighborhoods.py`

**Purpose:** Identify cell type co-localization changes during progression using squidpy.

**Key Functions:**
```python
def setup_spatial_coordinates(adata):
    # Creates adata.obsm['spatial'] from cell_x, cell_y
    adata.obsm['spatial'] = np.column_stack([adata.obs['cell_x'], adata.obs['cell_y']])

def build_spatial_graph_per_sample(adata, sample_id, n_neighbors=15):
    sample_adata = adata[adata.obs['sample_id'] == sample_id].copy()
    sq.gr.spatial_neighbors(sample_adata, coord_type='generic', n_neighs=n_neighbors)

def compute_nhood_enrichment_per_sample(adata, cluster_key='cell_type', n_perms=500)
def aggregate_nhood_enrichment_by_stage(adata, stage)
def compare_nhood_enrichment_stages(z_stage1, z_stage2, ...)
```

**Statistical Methods:**
- Squidpy `nhood_enrichment()` with 500 permutations
- Welch's t-test comparing sample-level z-scores between stages
- BH-FDR correction

**Outputs:**
```
results/spatial_neighborhoods/
├── figures/
│   ├── fig1_nhood_enrichment_normal.png
│   ├── fig2_nhood_enrichment_metaplasia.png
│   ├── fig3_nhood_enrichment_cancer.png
│   ├── fig4_nhood_diff_N_to_M.png
│   └── fig5_top_changing_pairs.png
├── nhood_enrichment_per_stage.csv
├── nhood_comparison_N_vs_M.csv
├── significant_changes.csv
└── SPATIAL_NEIGHBORHOOD_REPORT.md
```

**Expected Result:** "CD8-Tumor co-localization decreases 40% from M→C (immune exclusion)"

---

### Script 3: `51_immune_niche_analysis.py`

**Purpose:** Characterize CD8 T cell microenvironments; compare exhausted vs functional CD8 niches.

**Key Functions:**
```python
EXHAUSTION_MARKERS = ['PDCD1', 'LAG3', 'HAVCR2', 'TIGIT', 'CTLA4']

def score_cd8_exhaustion(adata, markers, method='mean_z')
def classify_cd8_exhaustion(scores, threshold_percentile=75)
def compute_niche_composition(adata, center_cells, k=15, cluster_key='cell_type')
def compare_niches(niche_df, group_col, group1, group2)  # Mann-Whitney U
def test_exhausted_vs_functional_niches(adata)
```

**Statistical Methods:**
- Z-score normalized exhaustion scoring
- k=15 nearest neighbor niche composition
- Mann-Whitney U for proportion comparisons
- Cohen's d effect size

**Outputs:**
```
results/immune_niche_analysis/
├── figures/
│   ├── fig1_exhaustion_score_distribution.png
│   ├── fig2_niche_composition_comparison.png
│   ├── fig3_niche_radar_plot.png
│   └── fig4_niche_by_stage.png
├── cd8_exhaustion_scores.csv
├── niche_composition_full.csv
├── niche_comparison.csv
└── IMMUNE_NICHE_REPORT.md
```

**Expected Result:** "Exhausted CD8s have 2.1x more CAF neighbors than functional CD8s (p<0.01)"

---

## Critical Implementation Details

### Data Loading Pattern
```python
BASE_DIR = Path('/home/user/g4x-choi-batch2-analysis')
INPUT_PATH = BASE_DIR / 'results' / 'qc_all_samples' / 'merged' / 'merged_corrected.h5ad'
OUTPUT_DIR = BASE_DIR / 'results' / '{analysis_name}'

adata = sc.read_h5ad(INPUT_PATH)
# 1,835,026 cells × 341 genes
```

### Key Column Names
```python
# Metadata
adata.obs['patient']     # 'SNU-105', 'SNU-107', 'SNU-484', 'ctrl'
adata.obs['stage']       # 'normal', 'metaplasia', 'cancer', 'control'
adata.obs['cell_type']   # 'Epithelial', 'CD8_T', 'Stromal', 'Macrophage', etc.
adata.obs['sample_id']   # 'A01', 'B02', etc.

# Spatial
adata.obs['cell_x']      # X coordinate
adata.obs['cell_y']      # Y coordinate

# Proteins (17 markers)
adata.obs['CD8_intensity_mean']
adata.obs['PD1_intensity_mean']
# etc.
```

### Constants (All Scripts)
```python
STAGE_ORDER = ['normal', 'metaplasia', 'cancer']
STAGE_COLORS = {'normal': '#2ecc71', 'metaplasia': '#f39c12', 'cancer': '#e74c3c'}
PAIRED_PATIENTS = ['SNU-105', 'SNU-107']
```

### Squidpy Spatial Graph (CRITICAL)
```python
# MUST build per-sample to avoid cross-section connections
for sample_id in adata.obs['sample_id'].unique():
    sample_adata = adata[adata.obs['sample_id'] == sample_id].copy()
    sample_adata.obsm['spatial'] = np.column_stack([
        sample_adata.obs['cell_x'],
        sample_adata.obs['cell_y']
    ])
    sq.gr.spatial_neighbors(sample_adata, coord_type='generic', n_neighs=15)
```

---

## Execution Order

```bash
conda activate enact
cd ~/g4x-choi-batch2-analysis

# Step 1 (PRIORITY - run first, ~10 min)
python scripts/49_within_patient_composition.py 2>&1 | tee logs/49_composition.log

# Step 2 (can run in parallel terminal, ~30 min)
python scripts/50_protein_spatial_neighborhoods.py 2>&1 | tee logs/50_spatial.log

# Step 3 (after 1 & 2 complete, ~20 min)
python scripts/51_immune_niche_analysis.py 2>&1 | tee logs/51_niche.log
```

---

## Verification

### Script 49 (Composition)
```bash
# Check outputs exist
ls results/within_patient_composition/figures/*.png
cat results/within_patient_composition/consistent_changes.csv | head

# Verify statistics
grep "Chi-square" results/within_patient_composition/COMPOSITION_REPORT.md
```

### Script 50 (Spatial Neighborhoods)
```bash
# Check nhood enrichment computed
ls results/spatial_neighborhoods/nhood_*.csv
head results/spatial_neighborhoods/significant_changes.csv
```

### Script 51 (Immune Niche)
```bash
# Check CD8 analysis
wc -l results/immune_niche_analysis/cd8_exhaustion_scores.csv  # Should be ~11,131
cat results/immune_niche_analysis/niche_comparison.csv | grep "p_adj"
```

### Integration Test
```python
# Quick validation in Python
import pandas as pd
comp = pd.read_csv('results/within_patient_composition/consistent_changes.csv')
print(f"Consistent changes found: {len(comp[comp['significant']])}")
```

---

## Files to Create

| Script | Lines (est.) | Priority |
|--------|--------------|----------|
| `scripts/49_within_patient_composition.py` | ~600 | 1 (FIRST) |
| `scripts/50_protein_spatial_neighborhoods.py` | ~500 | 2 |
| `scripts/51_immune_niche_analysis.py` | ~550 | 3 |

---

## Edge Cases Handled

1. **Missing cell types in stages**: Pseudocount (0.5) for Fisher's exact
2. **Rare populations (<10 cells)**: Exclude from per-celltype tests, document
3. **Sample-specific spatial graphs**: Build per-sample, never cross-section
4. **CD8 T cell scarcity**: 11,131 total; aggregate across samples within stage
5. **Zero-variance scores**: Guard with `np.std() < 1e-9` checks
