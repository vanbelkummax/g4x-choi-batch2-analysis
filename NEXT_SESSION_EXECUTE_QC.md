# G4X Full Dataset QC - Execute Pipeline

**Created:** 2026-01-21
**Status:** VERIFIED AND READY TO EXECUTE

---

## Task

Execute the verified G4X QC pipeline for all 32 samples. The plan has been reviewed, critiques addressed, and preflight checks passed.

---

## File Locations

### Raw Data (32 samples across 4 lanes)

```
/mnt/x/Choi_Batch_2_Tuesday/
├── g4-028-083-FC1-L001_5XGaAe5DB2dm7sRK/   # Lane 1 (A01-H01)
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
│   ├── 60_load_all_samples.py      # Step 1: Load 32 samples → h5ad
│   ├── 61_comprehensive_qc.py      # Step 2: Sample-level QC + batch effects
│   ├── 62_process_qc_passing.py    # Step 3: Cell-level QC + WNN integration
│   └── preflight_check.py          # Pre-execution verification (PASSED)
├── docs/
│   └── FULL_DATASET_QC_IMPLEMENTATION_PLAN.md
├── logs/                           # Pipeline output logs (created)
├── results/qc_all_samples/         # Output directory (created)
│   ├── raw/                        # Per-sample h5ad files
│   ├── figures/                    # QC visualizations
│   │   ├── per_sample/
│   │   ├── cross_sample/
│   │   └── batch_effects/
│   └── final_processed/            # Integrated results
├── NEXT_SESSION_VERIFY_AND_RUN.md  # Previous session doc (updated)
└── NEXT_SESSION_EXECUTE_QC.md      # THIS FILE
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

## Pipeline Overview

### Script 60: Load All Samples
- **Input:** 32 feature_matrix.h5 files
- **Output:** `results/qc_all_samples/raw/{sample}_raw.h5ad` (32 files)
- **Actions:**
  - Load RNA + Protein from H5
  - Separate into RNA (X) and Protein (obsm['protein'])
  - Add metadata (sample_id, lane, stage, patient)
  - Compute basic QC metrics (n_counts, n_genes, n_proteins)

### Script 61: Comprehensive QC
- **Input:** 32 raw h5ad files
- **Output:**
  - `sample_qc_summary.csv` - Pass/Warn/Fail for each sample
  - `QC_REPORT.md` - Comprehensive report
  - `figures/` - QC visualizations
- **Thresholds:**
  | Metric | Value | Action |
  |--------|-------|--------|
  | min_cells | 20,000 | WARN if below |
  | min_median_transcripts | 30 | FAIL if below |
  | min_median_genes | 20 | FAIL if below |
  | max_pct_empty | 5% | WARN if above |
  | LISI (lane) | > 3.0 | CONCERN if below |
  | Silhouette (lane) | < 0.3 | CONCERN if above |
  | PC1 Lane R² | < 20% | CONCERN if above |

### Script 62: Process QC-Passing Samples
- **Input:** QC-passing samples from step 61
- **Output:** `results/qc_all_samples/final_processed/{sample}_final.h5ad`
- **Pipeline per sample:**
  1. Compute admixture scores (conflicting marker co-expression)
  2. Apply cell-level QC (min_counts=10, min_genes=5, filter admixed)
  3. Normalize RNA (CPM + log1p)
  4. PCA on RNA
  5. WNN integration (variance-weighted RNA + Protein)
  6. Lineage annotation (epithelial/immune/stromal/endothelial)
  7. Cell type annotation (T cells, B cells, macrophages, etc.)
  8. Leiden clustering + UMAP

---

## QC Thresholds Summary

### Sample-Level (Script 61)
```python
QC_THRESHOLDS = {
    'min_cells': 20_000,
    'min_median_transcripts_per_cell': 30,
    'min_median_genes_per_cell': 20,
    'max_pct_empty': 5.0,
    'min_median_protein_counts': 5.0,
    'min_pct_protein_positive': 80.0,
}

BATCH_THRESHOLDS = {
    'max_silhouette_batch': 0.3,
    'min_lisi': 3.0,
    'max_pc1_lane_variance': 0.20,  # R² < 20%
}
```

### Cell-Level (Script 62)
```python
CELL_QC_THRESHOLDS = {
    'min_counts': 10,
    'min_genes': 5,
    'filter_admixed': True,  # Can override with --keep-admixed
}

ADMIX_THRESHOLD = 0.7   # 70th percentile for "high" expression
ADMIX_CUTOFF = 0.3      # Flag cells above this score
```

---

## Critiques Addressed (2026-01-21)

| Issue | Severity | Resolution |
|-------|----------|------------|
| Harmony flag no-op | HIGH | **Removed** - batch correction only makes sense post-merge |
| min_pct_in_cells unused | MEDIUM | **Removed** - metric unavailable; Resolve baseline validates instead |
| Lane variance arbitrary F-stat | MEDIUM | **Fixed** - Now uses R² with 20% threshold |
| Admixed cells hard-dropped | MEDIUM | **Fixed** - Added `--keep-admixed` flag |
| WNN simplified | MEDIUM | **Kept** - documented; 17 markers don't need true WNN |

---

## Execution Commands

```bash
# Activate environment
conda activate enact
cd ~/g4x-choi-batch2-analysis

# Verify ready (already passed)
python scripts/preflight_check.py

# Step 1: Load all samples (~20 min)
python scripts/60_load_all_samples.py 2>&1 | tee logs/60_loading.log

# Step 2: Comprehensive QC (~30 min)
python scripts/61_comprehensive_qc.py 2>&1 | tee logs/61_qc.log

# Step 2b (if interrupted): Resume from checkpoint
# python scripts/61_comprehensive_qc.py --resume 2>&1 | tee logs/61_qc_resume.log

# Step 3: Review QC report before proceeding
cat results/qc_all_samples/QC_REPORT.md

# Step 4: Process passing samples (~3 hours with 8 workers)
python scripts/62_process_qc_passing.py --parallel 8 2>&1 | tee logs/62_processing.log

# Alternative: Sensitivity analysis (keep admixed cells flagged but not removed)
# python scripts/62_process_qc_passing.py --parallel 8 --keep-admixed 2>&1 | tee logs/62_processing_sensitivity.log

# Step 5: Merge and batch correct (if needed)
python scripts/63_merge_and_batch_correct.py 2>&1 | tee logs/63_merge.log

# Alternative: Force Harmony even if metrics pass
# python scripts/63_merge_and_batch_correct.py --force 2>&1 | tee logs/63_merge.log

# Alternative: Use scVI for more sophisticated correction
# python scripts/63_merge_and_batch_correct.py --method scvi 2>&1 | tee logs/63_merge_scvi.log

# Alternative: Use true WNN with muon
# python scripts/63_merge_and_batch_correct.py --wnn-muon 2>&1 | tee logs/63_merge_wnn.log
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

### Known Issues
- **H04** will FAIL (25 transcripts/cell, 11 genes/cell)
- **D01** will WARN (11.9% empty cells)
- **C03** may WARN (21K cells - borderline)
- **Lane 4** has smaller cell areas (54-70 µm² vs 80-120 µm²) - monitor for confounding

---

## Checkpoints

| After Step | Verify |
|------------|--------|
| Script 60 | `ls results/qc_all_samples/raw/*.h5ad | wc -l` → 32 files |
| Script 61 | `grep "FAIL" results/qc_all_samples/sample_qc_summary.csv` → H04 only |
| Script 61 | `cat results/qc_all_samples/QC_REPORT.md` → Review batch effects |
| Script 62 | `ls results/qc_all_samples/final_processed/*.h5ad | wc -l` → 31 files |
| Script 63 | `cat results/qc_all_samples/merged/MERGE_REPORT.md` → Review batch metrics |
| Script 63 | `ls results/qc_all_samples/merged/merged_corrected.h5ad` → Final merged file |

---

## Post-Processing (After QC Complete)

Script 63 handles merge and batch correction automatically. After completion:

1. ✅ **Merge all processed samples** - Done by script 63
2. ✅ **Batch correction** - Harmony applied automatically if LISI < 3.0 or silhouette > 0.3
3. **Continue with ROADMAP.md** advanced analyses:
   - Progression analysis with pseudotime streamlines
   - CD8 exhaustion niche tensor
   - CAF network topology
   - TLS detection via TDA
   - Immune exclusion vector fields

### Batch Correction Options

| Method | When to Use | Command |
|--------|-------------|---------|
| Harmony (default) | Standard batch correction | `--method harmony` |
| scVI | Stronger correction needed | `--method scvi` |
| None | Skip correction | `--method none` |
| Force | Apply even if metrics pass | `--force` |
| True WNN | Multi-modal integration | `--wnn-muon` |

---

## Quick Reference

```bash
# Check sample count
find /mnt/x/Choi_Batch_2_Tuesday -name "feature_matrix.h5" | wc -l
# Expected: 32

# Check Resolve baseline
head -3 /mnt/x/Choi_Batch_2_Tuesday/choi_preGC_b2_core_metrics.csv

# View QC report after step 2
cat ~/g4x-choi-batch2-analysis/results/qc_all_samples/QC_REPORT.md

# Check processing results after step 4
cat ~/g4x-choi-batch2-analysis/results/qc_all_samples/final_processed/processing_results.csv
```

---

## Environment

```bash
conda activate enact
# Key packages: scanpy 1.10.3, squidpy 1.6.1, anndata 0.10.9
```

---

## Decision Log

| Question | Decision | Rationale |
|----------|----------|-----------|
| Enforce min_pct_in_cells? | No | Metric unavailable in data |
| Lane effect: F-stat or R²? | R² with 20% threshold | More interpretable |
| Keep admixed cells? | Optional via flag | Default filters; `--keep-admixed` for sensitivity |
| Implement Harmony? | Post-merge only | Per-sample doesn't make sense |
| True WNN vs simplified? | Keep simplified | 17 markers don't warrant complexity |
