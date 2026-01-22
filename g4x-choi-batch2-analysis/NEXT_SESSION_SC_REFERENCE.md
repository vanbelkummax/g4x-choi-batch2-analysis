# Next Session: SC Reference Annotation & Resolution Optimization

**Resume:** `cd /home/user/g4x-choi-batch2-analysis && cat NEXT_SESSION_SC_REFERENCE.md`

---

## Goal
Use Kumar et al. GSE183904 (gastric cancer scRNA-seq atlas) to:
1. Annotate G4X cells with reference-based labels
2. Find optimal clustering resolution using NMI/ARI/Purity
3. Validate Goblet expansion finding (6% → 48%)

---

## Reference Dataset
- **GEO:** GSE183904
- **Paper:** Kumar et al. Cancer Discovery 2022 (PMID: 34642171)
- **Data:** 31 gastric tumors, 34 cell lineage states
- **Download:** `wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183904/suppl/GSE183904_RAW.tar`

---

## Implementation Plan

### Step 1: Download Reference (~5 min)
```bash
mkdir -p ~/data/kumar_gastric_reference
cd ~/data/kumar_gastric_reference
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183904/suppl/GSE183904_RAW.tar"
tar -xf GSE183904_RAW.tar
```

### Step 2: Process Reference + Label Transfer (~30 min)
- Find overlapping genes (G4X 337 ∩ Kumar ~20K)
- Use scANVI or kNN-based label transfer
- Assign Kumar cell types to G4X cells

### Step 3: Resolution Optimization (~20 min)
```python
resolutions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
# For each: cluster → compute NMI, ARI, Purity vs reference labels
# Composite = 0.4*NMI + 0.3*ARI + 0.3*Purity
# Pick resolution with max composite
```

### Step 4: Validate & Compare
- Compare marker-based vs reference-based Goblet %
- Generate comparison figures

---

## Current State (Completed)
- ✅ 134,467 cells loaded (E02/F02/G02)
- ✅ Manual marker annotation (19 types, 63 genes)
- ✅ HVG optimization (best: all genes, res=0.3)
- ✅ Key finding: Goblet 6%→48% in cancer
- ✅ Export package: `Desktop/G4X_Marker_Analysis.zip`

---

## Key Files
```
/home/user/g4x-choi-batch2-analysis/
├── results/pilot/*_annotated.h5ad    # Current data
├── output/data/celltype_proportions.csv
├── export_marker_analysis/           # Marker-based results
└── PLAN_SC_REFERENCE_OPTIMIZATION.md # Detailed plan
```

---

## Quick Start Commands
```bash
conda activate enact
cd /home/user/g4x-choi-batch2-analysis

# Check current marker-based proportions
head output/data/celltype_proportions.csv

# View detailed plan
cat PLAN_SC_REFERENCE_OPTIMIZATION.md
```
