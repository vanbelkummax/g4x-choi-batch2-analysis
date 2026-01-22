# Resolution Optimization Plan: Biology-Guided Clustering

## Problem
Silhouette score optimization leads to **underclustering** (fewer clusters = higher silhouette).
We need a metric that rewards **biological fidelity** - clusters should correspond to real cell types.

---

## Mathematical Framework

### Metric 1: Cluster Purity (Higher = Better)
How "pure" is each cluster in terms of cell type?

```
Purity = (1/N) × Σ max_k(n_ik)

where:
- N = total cells
- n_ik = cells of type k in cluster i
```

**Interpretation:** If every cluster contains only one cell type, purity = 1.0

### Metric 2: Normalized Mutual Information (NMI)
How much information do clusters share with cell type labels?

```
NMI = I(Clusters; CellTypes) / sqrt(H(Clusters) × H(CellTypes))

where:
- I = mutual information
- H = entropy
```

**Range:** 0 (no correspondence) to 1 (perfect correspondence)

### Metric 3: Adjusted Rand Index (ARI)
Chance-corrected measure of cluster-celltype agreement.

```
ARI = (RI - Expected_RI) / (Max_RI - Expected_RI)
```

**Range:** -1 to 1 (1 = perfect, 0 = random, <0 = worse than random)

### Composite Score (Proposed)
```
OptimalScore = 0.4×NMI + 0.3×ARI + 0.3×Purity
```

Maximize across resolutions.

---

## Algorithm

```
For each resolution in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]:
    For each sample in [E02, F02, G02]:
        1. Cluster with Leiden at this resolution
        2. Get cell type labels from marker scoring
        3. Compute Purity, NMI, ARI
        4. Store results

    Compute mean score across samples (for consistency)

Select resolution with highest mean composite score
```

---

## Implementation Steps

### Step 1: Marker-based cell type annotation (already done)
- Uses 19 cell types, 63 marker genes
- Each cell gets a label based on highest marker score

### Step 2: Multi-resolution clustering
- Test resolutions: 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0
- Fixed: all 337 genes, 30 PCs

### Step 3: Compute metrics per resolution per sample
- Purity, NMI, ARI
- Also track: n_clusters, silhouette (for comparison)

### Step 4: Aggregate across samples
- Mean ± std for each metric
- Identify resolution that works best across all 3 samples

### Step 5: Visualize
- Line plot: resolution vs each metric
- Highlight optimal resolution
- Confusion matrices at optimal resolution

---

## Expected Outcome

| Resolution | Clusters | Purity | NMI | ARI | Composite |
|------------|----------|--------|-----|-----|-----------|
| 0.1 | ~2-3 | Low | Low | Low | Low |
| 0.3 | ~4-6 | Med | Med | Med | Med |
| 0.5 | ~6-10 | ? | ? | ? | ? |
| 0.8 | ~10-15 | ? | ? | ? | ? |
| 1.0 | ~12-20 | Med | Med | Med | Med |

**Hypothesis:** Optimal will be somewhere in 0.4-0.6 range where clusters capture major cell types without over-splitting.

---

## Output Files

1. `resolution_optimization.csv` - All metrics per resolution per sample
2. `resolution_comparison.png` - Line plots of metrics vs resolution
3. `optimal_confusion_matrix.png` - Cluster × CellType at optimal resolution
4. `optimal_recommendation.txt` - Final recommendation with justification

---

## Timeline
- Script: ~10 min to write
- Execution: ~15 min (8 resolutions × 3 samples × clustering + metrics)
- Total: ~25 min
