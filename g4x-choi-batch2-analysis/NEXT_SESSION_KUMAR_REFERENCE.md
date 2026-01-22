# Continue Session: Kumar et al. True Reference Investigation

**Resume:** `cd /home/user/g4x-choi-batch2-analysis && conda activate enact`

## ISSUE
Current "reference annotation" is a **pseudo-reference**:
- Cells from Kumar et al. GSE183904 (real paper)
- BUT cell types assigned by OUR marker scoring, not paper's expert annotations
- Need TRUE reference labels from original paper

## INVESTIGATE

### 1. Check GEO for metadata
```bash
# Check if Kumar et al. provided cell type annotations
ls -la reference/GSE183904/
# Look for metadata files with cell type labels
find reference/GSE183904/ -name "*.txt" -o -name "*metadata*" -o -name "*annotation*"
```

### 2. Paper details
- **Paper:** Kumar et al. Cancer Discovery 2022 (PMID: 34642171)
- **GEO:** GSE183904
- **Title:** "Single-Cell Atlas of Lineage States, Tumor Microenvironment..."
- Check supplementary tables for cell type annotations

### 3. If annotations exist
Update `scripts/06_reference_annotation.py` to use paper's labels instead of marker scoring.

### 4. If no annotations
Options:
- Use CellTypist with gastric model
- Use scType automatic annotation
- Accept pseudo-reference (trends still valid)

## COMPLETED (2026-01-22)
- ✅ Professional figures (2a, 2b, 3, 4, 5, 6) in `~/Desktop/G4X_Figures/`
- ✅ HVG optimization (best: all genes)
- ✅ Resolution optimization: **OPTIMAL = 0.6** (NMI=0.259, ARI=0.182)
- ✅ GitHub pushed: `addf3af5ff`

## RESOLUTION RESULTS
| Sample | Best Res | Composite |
|--------|----------|-----------|
| E02 | 0.6 | 0.263 |
| F02 | 0.4 | 0.252 |
| G02 | 0.6 | 0.283 |

Results in: `results/pilot/resolution_optimization.csv`

## KEY FINDING (VALIDATED)
Goblet expansion: 10.9% → 33.5% (3x increase Normal→Cancer)
Both marker-based AND pseudo-reference confirm this trend.

## FILES
- `reference/kumar_gastric_normal_reference.h5ad` - 20,814 cells, NO celltype column
- `scripts/06_reference_annotation.py` - current pseudo-reference approach
- `results/pilot/*_reference_annotated.h5ad` - annotated with pseudo-reference
