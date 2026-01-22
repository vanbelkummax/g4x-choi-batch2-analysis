# G4X Pilot Analysis - Code Review Context

## Project Overview

**Goal:** Analyze 3 gastric cancer progression samples (Normal → Metaplasia → Cancer) from G4X platform with optimized clustering, hierarchical cell annotation, spatial visualization with H&E overlay, and exploratory DEG analysis.

**Repository:** `/home/user/g4x-choi-batch2-analysis/`
**Branch:** `pilot-clean`
**Environment:** `conda activate enact`

---

## Files to Review

### Existing Scripts (UPDATE NEEDED)
| File | Lines | Action |
|------|-------|--------|
| `scripts/01_hvg_optimize.py` | 75 | Minor updates |
| `scripts/02_annotate.py` | 97 | **REWRITE** - uses non-existent markers |
| `scripts/03_spatial_viz.py` | 85 | Add H&E overlay |

### New Scripts (TO CREATE)
| File | Lines Est. | Purpose |
|------|------------|---------|
| `scripts/00_load_pilot_samples.py` | ~150 | Load raw data with count preservation |
| `scripts/04_exploratory_deg.py` | ~150 | Cell-type specific DEGs |

---

## Raw Data Locations

### Sample Data (3 samples, same patient, L002)
```
/mnt/x/Choi_Batch_2_Tuesday/g4-028-083-FC1-L002_5XGaAe5DB2dm7sRK/
├── E02/  (Normal, 31,790 cells)
├── F02/  (Metaplasia, 33,857 cells)
└── G02/  (Cancer, 79,634 cells)

Each sample contains:
├── single_cell_data/
│   ├── cell_by_transcript.csv.gz    # RNA counts (337 genes)
│   ├── cell_by_protein.csv.gz       # Protein IF data
│   ├── cell_metadata.csv.gz         # Coordinates + cell metrics
│   └── feature_matrix.h5            # Alternative format
├── h_and_e/
│   ├── h_and_e.jp2                  # Full H&E image (~197MB)
│   ├── h_and_e_thumbnail.jpg        # Thumbnail (3047x3840)
│   ├── eosin.jp2, nuclear.jp2       # Separate stains
│   └── *_thumbnail.png
└── g4x_viewer/
    ├── E02_HE.ome.tiff              # H&E with metadata
    └── E02_run_metadata.json        # Platform metrics
```

### Coordinate System
- **Cell coordinates:** X: 36-19,107 | Y: 29-15,199 (pixels)
- **H&E thumbnail:** 3,047 × 3,840 pixels
- **Scale factor:** NOT uniform (0.16 X, 0.25 Y) - needs investigation

---

## Gene Panel (337 genes)

**Full list:** `/home/user/g4x-choi-batch2-analysis/g4x_gene_panel.txt`

### Key Markers VERIFIED Present (95/95)
```python
MARKERS = {
    # Gastric Epithelial
    'Gastric_Pit': ['MUC5AC', 'MUC6', 'TFF1', 'TFF2'],
    'Gastric_Chief': ['PGC', 'ATP4A'],
    'Goblet': ['MUC2', 'TFF3'],
    'Intestinal_Meta': ['CDX1', 'CDX2', 'MUC2', 'CDH17'],
    'Enteroendocrine': ['GHRL', 'SST', 'GAST'],
    'Stem_Progenitor': ['LGR5', 'PROM1', 'OLFM4', 'SOX2', 'SOX9'],
    'Epithelial_General': ['EPCAM', 'CDH1', 'CLDN1', 'CLDN3', 'CLDN4', 'CLDN7', 'CLDN18', 'TACSTD2'],

    # Immune - T Cells
    'T_CD4': ['CD3D', 'CD3E', 'CD4', 'IL7R'],
    'T_CD8_Cytotoxic': ['CD3D', 'CD3E', 'CD8A', 'GZMA', 'GZMB', 'GZMH', 'GZMK', 'PRF1', 'IFNG'],
    'T_Reg': ['CD3D', 'CD4', 'FOXP3', 'IL2RA'],
    'T_Exhausted': ['PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT'],

    # Immune - B/Plasma
    'B_Cell': ['MS4A1', 'CD19', 'CD79A'],
    'Plasma': ['IGHA1', 'IGHG1', 'IGHM', 'JCHAIN', 'SDC1'],

    # Immune - Myeloid
    'Macrophage': ['CD68', 'CD163', 'CSF1R', 'MRC1'],
    'Monocyte': ['CD14', 'ITGAM'],

    # Stromal
    'Fibroblast': ['COL1A1', 'LUM', 'VCAN', 'FN1'],
    'CAF_Myofibroblast': ['ACTA2', 'TAGLN', 'PDGFRA'],
    'CAF_Inflammatory': ['FAP', 'PDPN', 'THY1', 'POSTN'],
    'Pericyte': ['PDGFRB', 'ACTA2'],
    'Endothelial': ['PECAM1', 'VWF'],

    # Cytokines
    'Inflammatory_Cytokines': ['CCL2', 'CCL3', 'CCL4', 'CCL5', 'CXCL9', 'CXCL10', 'CXCL11', 'CXCL13', 'IL1B', 'IL6', 'IL10', 'TGFB1', 'TNF'],
}
```

### Markers NOT in Panel (from existing script)
```
KRT8, KRT18, GKN1, VIL1, COL1A2, DCN, CDH5, CD2
```
**Existing `02_annotate.py` uses these - will fail silently!**

### Additional Useful Markers Available (not in current plan)
```
VIM, CHGA, MKI67, CD44, CD47, CD276, STAT3, STAT1, S100A9,
HLA-DRA, HLA-A, IDO1, NKG7, GNLY, KLRD1, CD40, CD80, CD86,
PRDM1, TBX21, GATA3, EOMES, KIT, CD34, LYVE1, RGS5, DES, MYH11
```

---

## Protein IF Markers (from metadata)
```
CD45 (pan-immune), CD3 (T), CD4, CD8, CD68 (macro), CD20 (B),
FOXP3 (Treg), PD1, PDL1, PanCK (epithelial), aSMA (stroma),
HLA-DR, KI67, CD11c, CD31, ATPase, Isotype (control)
```

---

## Critical Issues Identified

### 🔴 HIGH PRIORITY
1. **`02_annotate.py` uses 7 non-existent markers** - script will score incorrectly
2. **CellTypist not installed** - `pip install celltypist` needed
3. **H&E coordinate mismatch** - non-uniform scaling requires transform

### 🟡 MEDIUM PRIORITY
4. No "Unknown" threshold in annotation (all cells assigned)
5. DEG uses Wilcoxon but cells aren't independent (pseudoreplication)
6. HVG selection may be inappropriate for 337-gene panel
7. Fixed Leiden resolution (0.5) is arbitrary
8. Protein IF data loaded but not used for validation

### 🟢 LOW PRIORITY
9. Missing proliferation signature (MKI67, PCNA, TOP2A)
10. Missing lymphatic endothelium (LYVE1)
11. Consider MAST instead of Wilcoxon for sparse data
12. Add CHGA for enteroendocrine validation

---

## Verification Commands

```bash
# Check marker availability
zcat /mnt/x/.../E02/single_cell_data/cell_by_transcript.csv.gz | head -1 | tr ',' '\n' | grep -i "EPCAM\|CD3D\|FAP"

# Check CellTypist
python -c "import celltypist; print('OK')"

# Check H&E metadata
python -c "import tifffile; print(tifffile.TiffFile('/mnt/x/.../E02_HE.ome.tiff').ome_metadata[:500])"

# Run tests
cd /home/user/g4x-choi-batch2-analysis
python scripts/01_hvg_optimize.py  # Should work
python scripts/02_annotate.py      # Will have issues with missing markers
```

---

## Expected Outputs

```
results/pilot/
├── E02_raw.h5ad, F02_raw.h5ad, G02_raw.h5ad
├── E02_annotated.h5ad, F02_annotated.h5ad, G02_annotated.h5ad
├── hvg_results.csv
├── deg_by_celltype.csv
├── figures/
│   ├── hvg_silhouette_comparison.png
│   ├── umap_{sample}.png
│   ├── spatial_{sample}.png
│   ├── spatial_{sample}_HE_overlay.png
│   ├── progression_3panel.png
│   ├── celltype_proportions.png
│   └── deg_volcano_*.png
└── PILOT_ANALYSIS_REPORT.md
```

---

## Statistical Caveats (Must Document)

```markdown
⚠️ EXPLORATORY ANALYSIS - NOT FOR PUBLICATION WITHOUT VALIDATION

- N = 1 patient (all 3 samples from same individual)
- Cells are technical replicates, NOT biological replicates
- P-values are MEANINGLESS for group comparisons
- Rankings by fold-change only, not statistical significance
- Any findings require validation in independent cohorts
```

---

## Related Files

| File | Purpose |
|------|---------|
| `/home/user/CLAUDE.md` | Project configuration |
| `/home/user/g4x-choi-batch2-analysis/PLAN.md` | High-level plan |
| `/home/user/g4x-choi-batch2-analysis/NEXT_SESSION.md` | QC data source |
| `/tmp/g4x_genes.txt` | Full gene panel list |

---

## Quick Start

```bash
conda activate enact
cd /home/user/g4x-choi-batch2-analysis

# Review existing scripts
cat scripts/01_hvg_optimize.py
cat scripts/02_annotate.py
cat scripts/03_spatial_viz.py

# Check raw data
ls -la /mnt/x/Choi_Batch_2_Tuesday/g4-028-083-FC1-L002_5XGaAe5DB2dm7sRK/E02/

# Gene panel
zcat /mnt/x/.../E02/single_cell_data/cell_by_transcript.csv.gz | head -1 | tr ',' '\n' | wc -l
```
