# Continue G4X: Spatial Autocorrelation Methods Comparison

## Resume
```bash
cd /home/user/g4x-choi-batch2-analysis
conda activate enact
```

## Task: Compare Spatial Autocorrelation Methods on Cancer Composite

Create a multi-panel figure comparing different spatial statistics methods for detecting cancer clustering:

### Methods to Compare (6 rows × 3 columns)
| Row | Method | What It Shows |
|-----|--------|---------------|
| 1 | **Raw Expression** | cancer_composite score (continuous) |
| 2 | **Moran's I** | Global spatial autocorrelation |
| 3 | **Local Moran's I (LISA)** | Local clusters (HH, LL, HL, LH) |
| 4 | **Getis-Ord Gi*** | Hot/cold spots |
| 5 | **Ripley's K** | Point pattern clustering at multiple scales |
| 6 | **Neighborhood Enrichment** | Cell type co-localization z-scores |

### Implementation Notes
```python
# For LISA and Getis-Ord, use:
from esda.moran import Moran_Local
from esda.getisord import G_Local
from libpysal.weights import KNN

# Or squidpy (if GPU available):
import squidpy as sq
sq.gr.spatial_autocorr(adata, mode='moran')  # Global
sq.gr.spatial_autocorr(adata, mode='geary')  # Geary's C
```

### Key Files
| File | Description |
|------|-------------|
| `output/overnight_annotated_v2.h5ad` | 134K cells with cancer_composite |
| `scripts/overnight_04_generate_figures.py` | Reference for spatial plotting |

### Figure Output
Save to: `C:\Users\User\Desktop\G4X_Overnight_Results\fig_spatial_autocorr_comparison.png`

## What's Done ✅
- 134,467 cells annotated with 25 cell types
- Cancer composite = mean(z(KI67) + z(Hwang) + z(Intestinal) + z(Proliferating))
- Cancer_Diffuse NOT in composite (weakens signal, markers are backwards)
- 10 figures generated for overnight analysis
- Validation: Cancer 6.3% → 8.0% → 34.7% (correct progression)

## Figures Already Created
1. fig_celltype_proportions_new.png
2. fig_spatial_cancer_new.png
3. fig_umap_new.png
4. fig_foxp3_cd8_spatial.png
5. fig_old_vs_new_comparison.png
6. fig_morans_i_heatmap.png
7. fig_neighborhood_enrichment.png
8. fig_marker_panel_table.png (25 cell types, 95 markers)
9. fig_annotation_methods.png (4-step pipeline)
10. fig_rna_vs_protein_markers.png
11. fig_cancer_components_spatial.png (6 rows: KI67, Hwang, Intestinal, Diffuse, Proliferating, Composite)
12. fig_cancer_overlay_spatial.png (multi-positive overlay)
13. fig_cancer_markers_table.png (composite formula + color key)

## Output Location
`C:\Users\User\Desktop\G4X_Overnight_Results\` (13 figures)
