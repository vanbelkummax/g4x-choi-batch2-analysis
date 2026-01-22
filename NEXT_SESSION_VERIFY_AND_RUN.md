# G4X Full Dataset QC - Verification and Execution

**Created:** 2026-01-21
**Purpose:** Verify pipeline plan before executing full 32-sample QC

---

## Task

Review the G4X full dataset QC pipeline, verify the plan is sound, and execute if approved.

---

## File Locations

### Raw Data (32 samples across 4 lanes)

```
/mnt/x/Choi_Batch_2_Tuesday/
├── g4-028-083-FC1-L001_5XGaAe5DB2dm7sRK/   # Lane 1
│   ├── A01/single_cell_data/feature_matrix.h5
│   ├── B01/single_cell_data/feature_matrix.h5
│   ├── C01/single_cell_data/feature_matrix.h5
│   ├── D01/single_cell_data/feature_matrix.h5
│   ├── E01/single_cell_data/feature_matrix.h5
│   ├── F01/single_cell_data/feature_matrix.h5
│   ├── G01/single_cell_data/feature_matrix.h5
│   └── H01/single_cell_data/feature_matrix.h5
├── g4-028-083-FC1-L002_5XGaAe5DB2dm7sRK/   # Lane 2 (A02-H02)
├── g4-028-083-FC1-L003_5XGaAe5DB2dm7sRK/   # Lane 3 (A03-H03)
├── g4-028-083-FC1-L004_5XGaAe5DB2dm7sRK/   # Lane 4 (A04-H04)
├── choi_preGC_b2_core_metrics.csv          # Resolve baseline QC
└── choi_preGC_b2_protein_core_metrics.csv  # Protein metrics
```

### GitHub Repository

```
~/g4x-choi-batch2-analysis/
├── scripts/
│   ├── 60_load_all_samples.py      # Load 32 samples → h5ad
│   ├── 61_comprehensive_qc.py      # Sample-level QC + batch effects
│   ├── 62_process_qc_passing.py    # Cell-level QC + WNN integration
│   └── preflight_check.py          # Pre-execution verification
├── docs/
│   └── FULL_DATASET_QC_IMPLEMENTATION_PLAN.md  # Detailed plan
├── logs/                           # Pipeline output logs
└── results/qc_all_samples/         # Output directory
    ├── raw/                        # Per-sample h5ad files
    ├── figures/                    # QC visualizations
    └── final_processed/            # Integrated results
```

### Skills

```
~/claude-skills/g4x-spatial-qc.md   # QC workflow reference
```

---

## Sample Metadata

| Sample | Stage | Tissue | Patient | Notes |
|--------|-------|--------|---------|-------|
| A0x | Control | CRC | ctrl | Positive control |
| B0x | Normal | Adjacent normal | SNU-105 | |
| C0x | Metaplasia | IM | SNU-105 | |
| D0x | Cancer | GC | SNU-105 | D01 has 11.9% empty (WARN) |
| E0x | Normal | Adjacent normal | SNU-107 | |
| F0x | Metaplasia | IM | SNU-107 | |
| G0x | Cancer | GC | SNU-107 | |
| H0x | Cancer | GC | SNU-484 | **H04 = FAIL** (25 trans/cell) |

---

## Bug Fixes Applied (2026-01-21)

| Issue | Fix |
|-------|-----|
| Cell-level QC not applied | Added filtering for <10 counts, <5 genes, high admixture BEFORE integration |
| **--harmony removed** | Flag was no-op; batch correction only makes sense post-merge |
| **--keep-admixed added** | Flags but doesn't remove admixed cells for sensitivity analysis |
| LISI threshold | Changed 2.0 → 3.0 |
| **min_pct_in_cells removed** | Metric unavailable in data; Resolve baseline validates capture instead |
| **Lane variance now uses R²** | Explicit PASS/FAIL at 20% threshold (not arbitrary F-stat cutoff) |
| Protein QC | Added min_median_protein_counts, min_pct_protein_positive |
| WNN oversimplified | Documented as variance-weighted concatenation |

---

## Verification Checklist

Please verify before execution:

1. [x] **Data exists:** 32/32 samples found
2. [x] **Scripts valid:** All syntax checks pass
3. [x] **Thresholds appropriate:**
   - Sample: min_cells=20K, min_median_trans=30, max_empty=5%
   - Cell: min_counts=10, min_genes=5, filter_admixed=True
   - Batch: LISI>3.0, silhouette<0.3, PC1_lane_R²<20%
4. [x] **Known exclusions:** H04 expected to FAIL
5. [x] **Cell QC order:** Admixture → Filter → Normalize → PCA → WNN
6. [x] **WNN limitation acknowledged:** Simplified, not true Seurat WNN
7. [x] **Critiques addressed:** See fixes above

---

## Execution Commands

```bash
conda activate enact
cd ~/g4x-choi-batch2-analysis

# Pre-flight (already passed)
python scripts/preflight_check.py

# Step 1: Load all samples (~20 min)
python scripts/60_load_all_samples.py 2>&1 | tee logs/60_loading.log

# Step 2: Comprehensive QC (~30 min)
python scripts/61_comprehensive_qc.py 2>&1 | tee logs/61_qc.log

# Step 3: Review QC report
cat results/qc_all_samples/QC_REPORT.md

# Step 4: Process passing samples (~3 hours, 8 workers)
# Default: filters admixed cells
python scripts/62_process_qc_passing.py --parallel 8 2>&1 | tee logs/62_processing.log

# Alternative: Keep admixed cells flagged for sensitivity analysis
# python scripts/62_process_qc_passing.py --parallel 8 --keep-admixed 2>&1 | tee logs/62_processing_sensitivity.log
```

---

## Expected Outcomes

| Metric | Expected |
|--------|----------|
| Samples loaded | 32 |
| Samples passing QC | 31 (H04 fails) |
| Total cells (raw) | ~2.3M |
| Total cells (post-QC) | ~2.0M |
| Complete N→M→C series | 2 patients (SNU-105, SNU-107) |

---

## Open Questions (RESOLVED)

1. ~~Should admixed cells be flagged instead of removed?~~ **RESOLVED:** Added `--keep-admixed` flag for sensitivity analysis
2. ~~Add Harmony post-processing script if LISI < 3.0?~~ **RESOLVED:** Will add post-merge if needed; per-sample Harmony doesn't make sense
3. ~~Switch to true muon WNN for production?~~ **RESOLVED:** Keep simplified WNN; 17 markers don't warrant complexity of true WNN

---

## Quick Verification Commands

```bash
# Check data exists
ls /mnt/x/Choi_Batch_2_Tuesday/g4-028-083-FC1-L001_5XGaAe5DB2dm7sRK/

# Check sample count
find /mnt/x/Choi_Batch_2_Tuesday -name "feature_matrix.h5" | wc -l
# Expected: 32

# Check Resolve baseline
head -3 /mnt/x/Choi_Batch_2_Tuesday/choi_preGC_b2_core_metrics.csv

# Run preflight
conda activate enact && python ~/g4x-choi-batch2-analysis/scripts/preflight_check.py
```
