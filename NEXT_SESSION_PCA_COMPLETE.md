# PCA Deep Dive H&E Fix Complete (2026-01-22)

## Issue Resolved
QuPath `thumbnail.jpg` files were mislabeled and did not match sample IDs.

## Solution Implemented
Updated `get_he_thumbnail()` in `scripts/46_pca_deep_dive.py` to use G4X source images:
- **OLD (wrong)**: `/mnt/x/Choi_Batch_2_Tuesday/Qupath_HE/data/{folder}/thumbnail.jpg`
- **NEW (correct)**: `/mnt/x/Choi_Batch_2_Tuesday/g4-028-083-FC1-L00{lane}_*//{sample_id}/h_and_e/h_and_e_thumbnail.jpg`

## Verification
See `results/pca_deep_dive/figures/old_vs_new_comparison.png` - B04 shows clear mismatch between old and new sources.

## Generated Outputs
```
results/pca_deep_dive/
├── figures/
│   ├── fig0_batch_validation.png      # LISI/Silhouette metrics
│   ├── fig1_variance_decomposition.png # Variance by factor
│   ├── fig2_sample_overview.png       # ✅ CORRECTED H&E thumbnails
│   ├── fig3_pseudobulk_pca.png        # Sample-level PCA
│   ├── fig4_cell_pca.png              # Cell-level PCA by stage/patient/cell type
│   ├── fig5_loading_analysis.png      # PC gene loadings
│   ├── fig6_multimodal_comparison.png # RNA vs Protein vs WNN
│   ├── fig7_summary_grid.png          # All sample panels
│   ├── fig8_modality_pcas.png         # RNA/Protein/Combined PCA
│   ├── fig9_spatial_pca.png           # Spatially-integrated PCA
│   └── sample_panels/{sample}_panel.png # 29 individual panels
├── data/
│   ├── variance_partition.csv
│   ├── pc_loadings.csv
│   ├── pseudobulk_pca.csv
│   └── progression_stats.json
└── PCA_DEEP_DIVE_REPORT.md
```

## Key Findings

### Variance Decomposition (PC1)
| Factor | R² | Category |
|--------|-----|----------|
| cell_type | 0.192 | Biological |
| lane | 0.056 | Technical |
| stage | 0.013 | Biological |
| patient | 0.012 | Technical |

**Interpretation**: Cell type is the dominant driver of PC1 variance (19.2%), while stage contributes only 1.3%. This suggests cell type composition differences dominate over stage-related transcriptional changes.

### Batch Correction
| Metric | Pre | Post |
|--------|-----|------|
| LISI | 2.462 | 2.662 |
| Silhouette | -0.0013 | -0.0060 |

Harmony correction improved batch mixing (higher LISI) while maintaining biological signal.

### Progression Statistics (Pseudobulk)
- Kruskal-Wallis p = 0.578 (not significant)
- Spearman ρ = -0.08 (p = 0.673)

**Note**: Limited statistical power at sample level (n=4-10 per stage). Cell-level analyses provide more robust inference.

## Next Steps
1. **PC Loading Deep Dive**: Examine which genes drive each PC
2. **Cell Type-Specific DE**: Differential expression within cell types across stages
3. **Trajectory Analysis**: CellRank on Harmony-corrected data
4. **Spatial Statistics**: Neighborhood enrichment, Ripley's K

## Code Change Summary
```python
# OLD CONFIG (removed)
HE_BASE = Path('/mnt/x/Choi_Batch_2_Tuesday/Qupath_HE/data')
SAMPLE_TO_QUPATH = {'A01': '1', ...}

# NEW CONFIG
G4X_BASE = Path('/mnt/x/Choi_Batch_2_Tuesday')
LANE_DIRS = {
    '01': 'g4-028-083-FC1-L001_5XGaAe5DB2dm7sRK',
    '02': 'g4-028-083-FC1-L002_5XGaAe5DB2dm7sRK',
    ...
}

def get_he_thumbnail(sample_id: str, cache_dir: Path):
    lane_suffix = sample_id[-2:]  # 'A01' -> '01'
    lane_dir = LANE_DIRS[lane_suffix]
    return G4X_BASE / lane_dir / sample_id / 'h_and_e' / 'h_and_e_thumbnail.jpg'
```
