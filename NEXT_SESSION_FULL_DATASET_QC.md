# NEXT SESSION: Full 32-Sample Analysis & Ruthless QC

**Date Created:** 2026-01-21
**Priority:** CRITICAL - We've been analyzing only 25% of the data!

---

## CRITICAL DISCOVERY

| What We Have | What We Analyzed | Gap |
|--------------|------------------|-----|
| **32 samples** | 8 samples (Lane 1 only) | **24 samples unprocessed** |
| **2.3M cells** | ~570K cells | **~1.7M cells missing** |
| **6 Normal** | 2 Normal | 4 Normal missing |
| **7 Metaplasia** | 2 Metaplasia | 5 Metaplasia missing |
| **10 Cancer** | 2 Cancer | 8 Cancer missing |

**Statistical Impact:** With all 32 samples, we have proper power for stage comparisons!

---

## SESSION OBJECTIVES

### 1. Load ALL 32 Samples (Not Just Lane 1)
### 2. Ruthless Visual & Numerical QC
### 3. Identify and Document Outliers
### 4. Create QC Report with Pass/Fail for Each Sample
### 5. Re-run Pipeline on QC-Passing Samples

---

## KEY REFERENCES (Read These First)

```bash
# Project instructions
cat ~/CLAUDE.md

# Full dataset documentation
cat ~/g4x-choi-batch2-analysis/README.md

# Sample metadata (32 samples)
cat /mnt/x/Choi_Batch_2_Tuesday/choi_preGC_b2_core_metrics.csv

# What's been analyzed (only Lane 1!)
ls ~/spatial-hackathon-2026/results/g4x_choi_batch2/annotated_v2/

# Strategic assessment with methods recommendations
cat ~/g4x-choi-batch2-analysis/docs/SENIOR_REVIEWER_STRATEGIC_ASSESSMENT.md
```

---

## DATA PATHS

### Raw Data (32 samples across 4 lanes)
```
/mnt/x/Choi_Batch_2_Tuesday/
├── g4-028-083-FC1-L001_5XGaAe5DB2dm7sRK/  # Lane 1: A01-H01 (ANALYZED)
│   ├── A01/, B01/, C01/, D01/, E01/, F01/, G01/, H01/
├── g4-028-083-FC1-L002_5XGaAe5DB2dm7sRK/  # Lane 2: A02-H02 (NOT ANALYZED)
│   ├── A02/, B02/, C02/, D02/, E02/, F02/, G02/, H02/
├── g4-028-083-FC1-L003_5XGaAe5DB2dm7sRK/  # Lane 3: A03-H03 (NOT ANALYZED)
│   ├── A03/, B03/, C03/, D03/, E03/, F03/, G03/, H03/
├── g4-028-083-FC1-L004_5XGaAe5DB2dm7sRK/  # Lane 4: A04-H04 (NOT ANALYZED)
│   ├── A04/, B04/, C04/, D04/, E04/, F04/, G04/, H04/
├── choi_preGC_b2_core_metrics.csv          # Sample-level metrics
├── choi_preGC_b2_protein_core_metrics.csv  # Protein SNR data
└── choi_preGC_b2_transcript_core_metrics.csv
```

### Each Sample Contains
```
{sample}/
├── single_cell_data/      # Cell-level data (key!)
│   └── *.parquet or *.csv
├── rna/                   # RNA count matrices
├── protein/               # Protein intensities
├── segmentation/          # Cell boundaries
├── h_and_e/              # H&E images
├── metrics/              # QC metrics
└── summary_{sample}.html  # Resolve's QC report
```

### Current Results (Lane 1 only)
```
~/spatial-hackathon-2026/results/g4x_choi_batch2/
├── annotated_v2/         # A01-H01 only!
├── wnn_integrated/
└── segmentation_qc/
```

---

## SAMPLE INVENTORY WITH QC FLAGS

| Sample | Lane | Patient | Stage | Cells | Trans/Cell | Genes/Cell | **QC Flag** |
|--------|:----:|---------|:-----:|------:|:----------:|:----------:|-------------|
| A01 | 1 | CRC ctrl | Ctrl | 127,792 | 42 | 28 | ⚠️ Low depth |
| B01 | 1 | SNU-103 | N | 30,194 | 69 | 36 | ✅ |
| C01 | 1 | SNU-103 | M | 38,588 | 47 | 29 | ✅ |
| D01 | 1 | SNU-105 | N | 31,711 | 64 | 35 | ✅ |
| E01 | 1 | SNU-105 | M | 34,353 | 60 | 35 | ✅ |
| F01 | 1 | SNU-105 | C | 62,875 | 34 | 26 | ⚠️ Low depth |
| G01 | 1 | SNU-109 | C | 158,727 | 81 | 46 | ✅ |
| H01 | 1 | SNU-109 | X | 85,842 | 84 | 42 | ✅ |
| A02 | 2 | CRC ctrl | Ctrl | 120,260 | 48 | 32 | ❓ Not analyzed |
| B02 | 2 | SNU-106 | N | 37,511 | 64 | 38 | ❓ Not analyzed |
| C02 | 2 | SNU-106 | C | 52,588 | 55 | 35 | ❓ Not analyzed |
| D02 | 2 | SNU-106 | N | 56,035 | 48 | 32 | ❓ Not analyzed |
| E02 | 2 | SNU-107 | N | 31,790 | 51 | 33 | ❓ Not analyzed |
| F02 | 2 | SNU-107 | M | 33,857 | 56 | 35 | ❓ Not analyzed |
| G02 | 2 | SNU-107 | C | 79,634 | 89 | 50 | ❓ Not analyzed |
| H02 | 2 | SNU-107 | C | 87,215 | 109 | 55 | ❓ Not analyzed |
| A03 | 3 | CRC ctrl | Ctrl | 145,070 | 48 | 31 | ❓ Not analyzed |
| B03 | 3 | SNU-108 | M | 38,095 | 57 | 33 | ❓ Not analyzed |
| C03 | 3 | SNU-108 | C | 21,231 | 47 | 30 | ⚠️ Low cell count |
| D03 | 3 | SNU-108 | M | 33,227 | 48 | 31 | ❓ Not analyzed |
| E03 | 3 | SNU-110 | C | 86,003 | 63 | 41 | ❓ Not analyzed |
| F03 | 3 | SNU-110 | X | 82,704 | 50 | 33 | ❓ Not analyzed |
| G03 | 3 | UTSW-4 | body | 45,583 | 67 | 40 | ❓ Not analyzed |
| H03 | 3 | UTSW-4 | fundus | 35,585 | 126 | 36 | ⚠️ Very high depth |
| A04 | 4 | CRC ctrl | Ctrl | 130,891 | 65 | 38 | ❓ Not analyzed |
| B04 | 4 | Hu13 | C | 119,211 | 67 | 29 | ❓ Not analyzed |
| C04 | 4 | Hu13 | M | 89,387 | 74 | 33 | ❓ Not analyzed |
| D04 | 4 | Hu13 | C | 169,683 | 75 | 32 | ❓ Not analyzed |
| E04 | 4 | Hu15 | X | 68,014 | 67 | 34 | ❓ Not analyzed |
| F04 | 4 | Hu15 | X | 78,831 | 62 | 30 | ❓ Not analyzed |
| G04 | 4 | Hu15 | X | 64,147 | 53 | 29 | ❓ Not analyzed |
| **H04** | 4 | Hu15 | X | 32,334 | **25** | **11** | **❌ LIKELY FAIL** |

### Pre-Identified Issues
1. **H04**: Only 25 transcripts/cell, 11 genes/cell - likely QC failure
2. **C03**: Lowest cell count (21K) - check tissue quality
3. **F01, A01**: Low transcript depth (34, 42) - may need different normalization
4. **Lane 4**: Smaller cell areas (~55-70 µm² vs ~85-125 µm²) - different segmentation?

---

## PATIENT-MATCHED PROGRESSION SERIES

### Complete N→M→C Series (MOST VALUABLE)
| Patient | Normal | Metaplasia | Cancer | All Analyzed? |
|---------|--------|------------|--------|---------------|
| **SNU-105** | D01 | E01 | F01 | ✅ Yes (Lane 1) |
| **SNU-107** | E02 | F02 | G02, H02 | ❌ No (Lane 2) |

### Partial Series
| Patient | Normal | Metaplasia | Cancer | Notes |
|---------|--------|------------|--------|-------|
| SNU-103 | B01 | C01 | - | N→M only |
| SNU-106 | B02, D02 | - | C02 | N→C, no M |
| SNU-108 | - | B03, D03 | C03 | M→C only |
| Hu13 | - | C04 | B04, D04 | M→C only |

**Key Insight:** We need Lane 2 (SNU-107) for a second complete N→M→C progression!

---

## RUTHLESS QC PROTOCOL

### Phase 1: Numerical QC (Automated)

```python
# QC thresholds to apply
QC_THRESHOLDS = {
    'min_cells': 20_000,           # Minimum cells per sample
    'min_transcripts_per_cell': 30, # Minimum median transcripts
    'min_genes_per_cell': 20,       # Minimum median genes
    'max_pct_empty': 5.0,           # Maximum % empty cells
    'min_pct_in_cells': 80.0,       # Minimum % transcripts in cells
}

# Automatic flags
def flag_sample(row):
    flags = []
    if row['number_cells'] < QC_THRESHOLDS['min_cells']:
        flags.append('LOW_CELLS')
    if row['median_transcripts_per_cell'] < QC_THRESHOLDS['min_transcripts_per_cell']:
        flags.append('LOW_DEPTH')
    if row['median_unique_genes_per_cell'] < QC_THRESHOLDS['min_genes_per_cell']:
        flags.append('LOW_GENES')
    if row['pct_empty_cells'] > QC_THRESHOLDS['max_pct_empty']:
        flags.append('HIGH_EMPTY')
    if row['pct_transcripts_in_cells'] < QC_THRESHOLDS['min_pct_in_cells']:
        flags.append('LOW_ASSIGNMENT')
    return flags if flags else ['PASS']
```

### Phase 2: Visual QC (Per Sample)

Generate these figures for EVERY sample:

1. **Spatial scatter**: Total counts heatmap
2. **Violin plots**: total_counts, n_genes, pct_mt (if available)
3. **Scatter**: total_counts vs n_genes
4. **Histogram**: genes per cell distribution
5. **Protein SNR bar**: All 17 markers
6. **Spatial cell type map**: After annotation

### Phase 3: Cross-Sample Comparison

1. **Ridge plot**: total_counts across all 32 samples
2. **PCA**: Pseudobulk per sample, color by lane (batch effect?)
3. **UMAP**: Single-cell, color by sample_id (check mixing)
4. **Heatmap**: Sample-sample correlation matrix

### Phase 4: Lane/Batch Effect Assessment

1. **Compare Lane 1 vs 2 vs 3 vs 4** on same patient (A01-A04 are all CRC control)
2. **Check cell area distributions** per lane (Lane 4 is different!)
3. **Protein SNR** per lane

---

## EXECUTION PLAN

### Step 1: Create Full-Dataset Loader Script

```python
# scripts/60_load_all_samples.py
"""Load all 32 samples from G4X Choi Batch 2"""

import os
import pandas as pd
import anndata as ad
from pathlib import Path

DATA_ROOT = Path("/mnt/x/Choi_Batch_2_Tuesday")
LANES = {
    'L001': 'g4-028-083-FC1-L001_5XGaAe5DB2dm7sRK',
    'L002': 'g4-028-083-FC1-L002_5XGaAe5DB2dm7sRK',
    'L003': 'g4-028-083-FC1-L003_5XGaAe5DB2dm7sRK',
    'L004': 'g4-028-083-FC1-L004_5XGaAe5DB2dm7sRK',
}
WELLS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

def get_all_sample_paths():
    """Return dict of sample_id -> path"""
    samples = {}
    for lane_num, lane_dir in enumerate(LANES.values(), 1):
        for well in WELLS:
            sample_id = f"{well}0{lane_num}"
            sample_path = DATA_ROOT / lane_dir / sample_id
            if sample_path.exists():
                samples[sample_id] = sample_path
    return samples

def load_sample(sample_path, sample_id):
    """Load a single G4X sample into AnnData"""
    # Find single_cell_data
    sc_dir = sample_path / 'single_cell_data'
    # ... implementation
    pass

# Load all
all_samples = get_all_sample_paths()
print(f"Found {len(all_samples)} samples")  # Should be 32
```

### Step 2: Create Comprehensive QC Script

```bash
# Create and run
python scripts/61_comprehensive_qc_all32.py
```

### Step 3: Generate QC Report

Output: `results/qc_all_samples/QC_REPORT.md`

```markdown
# G4X Choi Batch 2 - Full QC Report

## Summary
- Total samples: 32
- Passed QC: XX
- Failed QC: XX
- Flagged for review: XX

## Failed Samples
| Sample | Reason | Action |
|--------|--------|--------|
| H04 | 11 genes/cell | EXCLUDE |
| ... | ... | ... |

## Per-Sample QC Figures
[Links to individual sample QC figures]
```

### Step 4: Re-Run Pipeline on QC-Passing Samples

```bash
# Only process samples that pass QC
python scripts/62_process_qc_passing.py
```

---

## EXPECTED OUTCOME

After this session:

1. **QC Report** for all 32 samples
2. **Processed AnnData** for QC-passing samples (~28-30 expected)
3. **Identified outliers** with documented exclusion reasons
4. **Lane/batch effect assessment**
5. **Updated sample inventory** with QC status

### Final Statistical Power

| Stage | Samples (Expected) | Cells (Expected) |
|-------|:------------------:|:----------------:|
| Normal | 5-6 | ~180K |
| Metaplasia | 6-7 | ~260K |
| Cancer | 8-10 | ~800K |
| Control/Other | 8-9 | ~750K |
| **TOTAL** | ~28-30 | ~2.0M |

---

## ENVIRONMENT

```bash
conda activate enact

# Verify
python -c "import scanpy, squidpy, anndata; print('Ready')"
```

---

## CHECKPOINTS

- [ ] All 32 sample paths verified accessible
- [ ] Sample loader script created and tested
- [ ] Numerical QC applied to all samples
- [ ] Visual QC generated for all samples
- [ ] Cross-sample comparison figures generated
- [ ] Lane/batch effects assessed
- [ ] H04 and other outliers documented
- [ ] QC report compiled
- [ ] Pipeline re-run on QC-passing samples

---

## NOTES FOR CLAUDE

1. **READ ~/CLAUDE.md first** - contains environment, tools, skills
2. **Use scanpy-analysis skill** for preprocessing
3. **Use squidpy-spatial skill** for spatial stats
4. **Do NOT skip QC** - this is the most important step
5. **Document everything** in QC report
6. **Create TodoWrite todos** for each phase
7. **Save figures to** `results/qc_all_samples/`

---

*Document created by senior reviewer strategic assessment session*
