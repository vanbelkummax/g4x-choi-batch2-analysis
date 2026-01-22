# Next Session: RCTD Cell Type Annotation for G4X

## Goal
Annotate G4X gastric samples (E02 Normal, F02 Metaplasia, G02 Cancer) using RCTD with GSE134520 reference.

## Files

**G4X Query Data:**
- `/home/user/g4x-choi-batch2-analysis/output/pilot_preprocessed.h5ad`
- 134,467 cells (E02: 27,775 / F02: 30,229 / G02: 76,463)

**Reference (Zhang et al. 2019):**
- `/home/user/g4x-choi-batch2-analysis/output/GSE134520_reference.h5ad`
- 44,012 cells, 18 cell types (Enterocyte, Goblet, PMC, Cancer, etc.)

**Marker genes:**
- `/home/user/g4x-choi-batch2-analysis/reference/GSE134520_markers.json`

## Workflow

1. **Convert to Seurat** - h5ad → Seurat objects (SeuratDisk)
2. **Build RCTD reference** - from GSE134520 Seurat object
3. **Run RCTD** - spacexr::run_RCTD() on each G4X sample
4. **Extract results** - cell_id, cell_type columns
5. **Merge back** - add annotations to AnnData
6. **Publication figures** - cell type composition, spatial distribution

## Expected Output
- `output/g4x_cell_annotations.csv` (cell_id, sample, cell_type, confidence)
- `output/pilot_annotated.h5ad` (with cell_type in obs)
- `output/fig_celltype_*.png/pdf` (publication-quality figures)

## Environment
```bash
conda activate enact
# R packages: Seurat 5.4.0, spacexr (RCTD), SeuratDisk
```

## Quick Start
```r
# In R
library(Seurat)
library(spacexr)
library(SeuratDisk)

# Load reference
Convert("GSE134520_reference.h5ad", dest="h5seurat")
ref <- LoadH5Seurat("GSE134520_reference.h5seurat")

# Build RCTD reference
# Run RCTD on G4X samples
# Export annotations
```
