# Next Session: Cell Type-Specific Differential Expression

## Context
PCA Deep Dive complete. **Key finding: No significant stage separation at sample level.**

- Kruskal-Wallis p = 0.578 (stage groups not different)
- Cell type explains 19.2% of PC1 variance; stage only 1.3%
- Conclusion: Cell composition masks progression signal → need **cell type-stratified analysis**

## Priority: Cell Type-Specific DE Analysis

Since bulk PCA doesn't separate stages, the biological signal likely exists **within** cell types across progression.

### Hypothesis
Epithelial cells, CAFs, and immune cells each have stage-specific transcriptional programs that are masked when analyzing all cells together.

### Proposed Analysis (`47_celltype_specific_de.py`)

1. **Epithelial cells: N→M→C trajectory**
   - Subset to epithelial (PanCK+)
   - DE: Normal vs Metaplasia vs Cancer
   - Key markers: CDX2, MUC2, MUC5AC, TFF1-3, intestinal vs gastric phenotype
   - Question: Does metaplasia show intermediate expression?

2. **CAF subtyping across stages**
   - Subset to stromal (aSMA+ or PDGFRA+)
   - Score mCAF/iCAF/apCAF signatures
   - Question: Do CAF phenotypes shift with progression?

3. **T cell exhaustion trajectory**
   - Subset to CD8+ T cells
   - Score exhaustion (PDCD1, LAG3, HAVCR2, TIGIT)
   - Question: Does exhaustion increase from N→M→C?

4. **Macrophage polarization**
   - Subset to CD68+ macrophages
   - M1 vs M2 scoring
   - Question: Does TME become more immunosuppressive?

### Statistical Approach
- Wilcoxon rank-sum for pairwise comparisons
- Benjamini-Hochberg FDR correction
- Effect size (Cohen's d or log2FC)
- Mixed model accounting for patient as random effect

### Expected Output
```
results/celltype_de/
├── epithelial_de_NvM.csv
├── epithelial_de_MvC.csv
├── epithelial_de_NvC.csv
├── caf_subtype_scores.csv
├── cd8_exhaustion_scores.csv
├── macrophage_polarization.csv
├── figures/
│   ├── epithelial_volcano_plots.png
│   ├── caf_composition_by_stage.png
│   ├── cd8_exhaustion_trajectory.png
│   └── heatmap_top_de_genes.png
└── CELLTYPE_DE_REPORT.md
```

### Data Location
- Input: `results/qc_all_samples/merged/merged_corrected.h5ad` (1.83M cells, 29 samples)
- Cell types in `adata.obs['cell_type']`
- Stage in `adata.obs['stage']` (normal, control, metaplasia, cancer)

### Commands
```bash
conda activate enact
cd ~/g4x-choi-batch2-analysis
python scripts/47_celltype_specific_de.py 2>&1 | tee logs/47_celltype_de.log
```

## Alternative: Spatial Analysis First?
If cell type DE also shows weak signal, consider:
- **Neighborhood-aware DE**: Cells in different spatial contexts may behave differently
- **Niche identification**: Define spatial niches, then compare niche composition across stages
- **Spatial gradients**: Look for immune exclusion zones, tumor-stroma boundaries

## Files
- PCA results: `results/pca_deep_dive/`
- Merged data: `results/qc_all_samples/merged/merged_corrected.h5ad`
- Analysis plan: `TIERED_ANALYSIS_PLAN.md`
