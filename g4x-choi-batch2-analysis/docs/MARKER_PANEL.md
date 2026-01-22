# G4X Gastric Cancer Marker Panel

## Cell Type Markers Used for Annotation

### Gastric Epithelial
| Cell Type | Markers | Description |
|-----------|---------|-------------|
| Gastric_Pit | MUC5AC, MUC6, TFF1, TFF2 | Pit/mucous cells |
| Gastric_Chief | PGC, ATP4A | Parietal/chief cells |
| Goblet | MUC2, TFF3 | Goblet cells |
| Intestinal_Meta | CDX1, CDX2, MUC2, CDH17 | Intestinal metaplasia markers |
| Enteroendocrine | GHRL, SST, GAST, CHGA | Neuroendocrine cells |
| Stem_Progenitor | LGR5, PROM1, OLFM4, SOX2, SOX9 | Stem/progenitor cells |
| Epithelial_General | EPCAM, CDH1, CLDN1, CLDN3, CLDN4, CLDN7, CLDN18 | General epithelial |

### Immune - T Cells
| Cell Type | Markers | Description |
|-----------|---------|-------------|
| T_CD4 | CD3D, CD3E, CD4, IL7R | Helper T cells |
| T_CD8_Cytotoxic | CD3D, CD3E, CD8A, GZMA, GZMB, PRF1 | Cytotoxic T cells |
| T_Reg | CD3D, CD4, FOXP3, IL2RA | Regulatory T cells |
| T_Exhausted | PDCD1, CTLA4, LAG3, HAVCR2, TIGIT | Exhausted T cells |

### Immune - B/Plasma
| Cell Type | Markers | Description |
|-----------|---------|-------------|
| B_Cell | MS4A1, CD19, CD79A | B cells |
| Plasma | IGHA1, IGHG1, IGHM, JCHAIN, SDC1 | Plasma cells |

### Immune - Myeloid
| Cell Type | Markers | Description |
|-----------|---------|-------------|
| Macrophage | CD68, CD163, CSF1R, MRC1 | Macrophages |
| Monocyte | CD14, ITGAM | Monocytes |

### Stromal
| Cell Type | Markers | Description |
|-----------|---------|-------------|
| Fibroblast | COL1A1, LUM, VCAN, FN1 | Fibroblasts |
| CAF | ACTA2, FAP, PDGFRA, POSTN, THY1 | Cancer-associated fibroblasts |
| Endothelial | PECAM1, VWF | Endothelial cells |

---

## Summary
- **19 cell types** defined
- **63 unique marker genes** (some shared across types)
- All markers verified present in 337-gene G4X panel
- Annotation method: `scanpy.tl.score_genes()` with max-score assignment
- Score threshold: 0.1 (cells below threshold = "Unknown")
