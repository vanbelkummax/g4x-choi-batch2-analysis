# Next Session: Add QC Data

## Goal
Add detailed QC metrics for pilot samples (E02, F02, G02) to the repo.

## Data Source
Raw transcripts from Lane 2:
```
/mnt/x/Choi_Batch_2_Tuesday/g4-028-083-FC1-L002_5XGaAe5DB2dm7sRK/
├── E02/resolved/E02_transcripts.csv
├── F02/resolved/F02_transcripts.csv
└── G02/resolved/G02_transcripts.csv
```

## Tasks

### 1. Create QC script
```bash
# scripts/00_qc_summary.py
# Compute per-sample:
# - n_cells, n_genes
# - median/mean counts per cell
# - median/mean genes per cell
# - % zero-count cells
# - protein detection rate
# - nuclei_area distribution
```

### 2. Generate QC report
Output: `results/qc_summary.csv`

| Metric | E02 (Normal) | F02 (Meta) | G02 (Cancer) |
|--------|--------------|------------|--------------|
| n_cells | ~32K | ~34K | ~80K |
| median_counts | ~51 | ~56 | ~89 |
| median_genes | ~33 | ~35 | ~50 |
| pct_zero | <1% | <1% | <0.2% |
| protein_positive | 100% | 100% | 100% |

### 3. Add QC figures
- `figures/qc_violin_counts.png` - counts distribution
- `figures/qc_violin_genes.png` - genes per cell
- `figures/qc_scatter_counts_genes.png` - counts vs genes

### 4. Update README
Add QC section with table and figure references.

## Quick Start
```bash
conda activate enact
cd ~/g4x-choi-batch2-analysis
python scripts/00_qc_summary.py
```

## Expected Output
```
results/
├── qc_summary.csv
└── figures/
    ├── qc_violin_counts.png
    ├── qc_violin_genes.png
    └── qc_scatter_counts_genes.png
```
