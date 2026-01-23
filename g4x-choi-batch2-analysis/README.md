# G4X RCTD + IM Analysis Pipeline

## Results Summary
- **134,467 cells** analyzed from 3 samples (E02 Normal, F02 Metaplasia, G02 Cancer)
- **87.7% RCTD annotated** (117,887 cells) with 18 cell types
- **Hwang Cancer Score**: Cancer (0.101) > Normal (-0.052) - validates progression
- **IM Ratio**: Cancer (0.469) > Normal (0.388) - increased intestinalization

## Scripts → Outputs

| Script | Output |
|--------|--------|
| `07_rctd_deconvolution.R` | `rctd_results.rds`, `rctd_reference.rds` |
| `07b_extract_rctd_results.R` | `g4x_rctd_annotations.csv`, `g4x_rctd_weights.csv` |
| `08_im_scoring.py` | `g4x_im_scores.csv` |
| `09_merge_results.py` | `pilot_final_annotated.h5ad` |
| `10_publication_figures.py` | `fig_integrated_panel.png`, `fig_qc_panel.png`, `fig_celltype_focus.png` |

## Top Cell Types
1. Enterocyte: 26,691 (22.6%)
2. Cancer: 15,082 (12.8%)
3. Fibroblast: 12,296 (10.4%)

## Next Steps
1. Spatial cell-cell interaction analysis by RCTD types
2. Validate with CARD/Tangram alternative methods
3. Integrate with G4X protein markers
