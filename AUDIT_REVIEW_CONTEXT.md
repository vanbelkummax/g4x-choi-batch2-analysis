# G4X QC Pipeline - Full Review Audit

**Generated:** 2026-01-21
**Project:** g4x-choi-batch2-analysis
**Branch:** master
**Audit Baseline Commit:** `c237dde` (this document)
**Environment:** `conda activate enact` (see `~/work/polymax/docs/SPATIAL_ANALYSIS_ENVIRONMENT.md`)

> **Note:** Line numbers reference commit `6c201a9` (QC code). If code changes, verify line numbers.

---

## 1. Project Overview

| Attribute | Value |
|-----------|-------|
| **Study** | Choi_GC_preGC_batch2 (Seock-Jin Chung, VUMC) |
| **Platform** | G4X (Resolve Biosciences) |
| **Total Samples** | 32 (4 lanes × 8 positions) |
| **RNA Panel** | 387 genes |
| **Protein Panel** | 17 markers |
| **Primary Question** | Spatial mechanisms of gastric cancer progression (N→M→C) |

---

## 2. Directory Structure

```
~/g4x-choi-batch2-analysis/
├── AUDIT_REVIEW_CONTEXT.md          # THIS FILE
├── README.md                        # Project overview
├── ROADMAP.md                       # Development milestones
├── PROVENANCE.md                    # Data lineage
│
├── G4X_ANALYSIS_PLAN.md             # MASTER ANALYSIS PLAN (primary)
├── FULL_DATASET_QC_PLAN.md          # QC strategy for 32 samples
├── VALIDATION.md                    # Validation criteria
├── METHOD_COMPARISON.md             # RNA vs Protein vs WNN comparison
│
├── config/
│   ├── qc_thresholds.yaml           # Threshold definitions
│   └── samples.yaml                 # Sample metadata
│
├── docs/
│   ├── DATA_DOCUMENTATION.md        # Data format documentation
│   ├── FULL_DATASET_QC_IMPLEMENTATION_PLAN.md  # Detailed QC implementation
│   ├── QC_DETAILED_CHECKLIST.md     # Step-by-step QC checklist
│   ├── G4X_Project_Ideas_Evidence_Based.md    # Analysis ideas
│   ├── PROJECT_IDEAS_RANKED.md      # Prioritized project list
│   ├── REVIEWER_CRITIQUE_PROMPT.md  # External review prompt
│   ├── SENIOR_REVIEWER_STRATEGIC_ASSESSMENT.md  # Strategic assessment
│   └── STRATEGIC_COMMAND_BRIEFING.md  # Executive summary
│
├── scripts/                         # ALL ANALYSIS SCRIPTS
│   ├── 60_load_all_samples.py       # Step 1: Load 32 samples
│   ├── 61_comprehensive_qc.py       # Step 2: Sample-level QC
│   ├── 62_process_qc_passing.py     # Step 3: Cell-level QC + WNN
│   ├── 63_merge_and_batch_correct.py # Step 4: Merge + batch correction
│   └── [01-59]_*.py                 # Legacy/pilot analysis scripts
│
├── results/
│   ├── qc_all_samples/              # FULL DATASET QC OUTPUT
│   │   ├── raw/                     # 32 raw h5ad files
│   │   ├── final_processed/         # QC-passing processed files
│   │   ├── figures/                 # QC figures
│   │   └── sample_qc_summary.csv    # QC metrics
│   ├── g4x_choi_batch2/             # Pilot analysis results (Lane 1)
│   │   ├── wnn_integrated/          # WNN integration output
│   │   └── validation/              # Validation results
│   ├── figures/                     # Generated figures
│   └── tables/                      # Generated tables
│
├── outputs/tables/                  # Analysis output tables
│
└── NEXT_SESSION_*.md                # Session continuation prompts
```

---

## 3. QC Pipeline Scripts (Audit Focus)

### Script 60: Load All Samples
**File:** `scripts/60_load_all_samples.py`
**Purpose:** Load all 32 G4X samples from raw feature matrices

| Input | Output |
|-------|--------|
| `/mnt/x/Choi_Batch_2_Tuesday/*/feature_matrix.h5` | `results/qc_all_samples/raw/{sample}_raw.h5ad` |
| | `results/qc_all_samples/raw/loading_manifest.csv` |

**Key Functions:**
- `load_sample()`: Parse Resolve HDF5 format
- `extract_metadata()`: Parse sample ID, lane, stage from filename
- Stores RNA in `.X`, protein in `.obsm['protein']`

---

### Script 61: Comprehensive QC
**File:** `scripts/61_comprehensive_qc.py`
**Purpose:** Sample-level QC with batch effect assessment

| Input | Output |
|-------|--------|
| `results/qc_all_samples/raw/*.h5ad` | `results/qc_all_samples/sample_qc_summary.csv` |
| `loading_manifest.csv` | `results/qc_all_samples/QC_REPORT.md` |
| | `results/qc_all_samples/figures/per_sample/*.png` |
| | `results/qc_all_samples/figures/cross_sample/*.png` |
| | `results/qc_all_samples/figures/batch_effects/*.png` |

**QC Thresholds (Line 104-112) with Provenance:**
```python
QC_THRESHOLDS = {
    'min_cells': 20_000,                      # G4X platform spec: typical 50K-200K cells/ROI
    'min_median_transcripts_per_cell': 30,    # Resolve validation: median 50-150 for gastric
    'min_median_genes_per_cell': 20,          # 387-gene panel: expect 30-80 genes/cell
    'max_pct_empty': 5.0,                     # Industry standard: <5% empty acceptable
    'min_median_protein_counts': 5.0,         # 17-marker panel: expect 10-50 counts
    'min_pct_protein_positive': 80.0,         # Platform QC: >80% cells should have signal
}
```

**Expected Value Ranges (for verification):**
| Metric | Expected Range | Red Flag If |
|--------|----------------|-------------|
| Cells per sample | 50,000 - 200,000 | <20,000 or >300,000 |
| Median transcripts/cell | 50 - 150 | <30 |
| Median genes/cell | 30 - 80 | <20 |
| Transcript density | 0.3 - 0.8 trans/µm² | >2 SD from lane mean |
| LISI (batch mixing) | >3.0 | <2.0 |
| Silhouette (lane) | <0.3 | >0.5 |

**Key Metrics Computed:**
- Cell counts, transcript counts, gene counts per sample
- **Transcript density (transcripts/µm²)** - NEW in commit `6c201a9`
- LISI (Local Inverse Simpson Index) for batch mixing
- Silhouette score for lane separation
- PC variance explained by lane (R²)

**Critical Code Review Points:**
1. `compute_sample_qc_metrics()` (Line 126-170) - Metric computation
2. `apply_qc_thresholds()` (Line 172-215) - Pass/Warn/Fail logic
3. `compute_lisi()` (Line 218-254) - Batch mixing assessment
4. `plot_cross_sample_comparison()` (Line 359-470) - Density plots
5. `generate_qc_report()` (Line 595-750) - Report generation

---

### Script 62: Process QC-Passing Samples
**File:** `scripts/62_process_qc_passing.py`
**Purpose:** Cell-level QC + WNN integration + annotation

| Input | Output |
|-------|--------|
| `results/qc_all_samples/raw/*.h5ad` | `results/qc_all_samples/final_processed/{sample}_final.h5ad` |
| `sample_qc_summary.csv` | `results/qc_all_samples/final_processed/processing_results.csv` |

**Cell-Level QC Thresholds (Line 167-171):**
```python
CELL_QC_THRESHOLDS = {
    'min_counts': 10,
    'min_genes': 5,
    'filter_admixed': True,
}
```

**Conflict Pairs for Admixture Detection (Line 152-164):**
```python
CONFLICT_PAIRS = [
    ('PanCK', 'CD45'),      # Epithelial + Immune
    ('PanCK', 'aSMA'),      # Epithelial + Stromal (Tumor-Stroma)
    ('CD45', 'aSMA'),       # Immune + Stromal
    ('CD31', 'CD45'),       # Endothelial + Immune
    ('CD3', 'CD68'),        # T cell + Macrophage
    ('CD4', 'CD8'),         # CD4 + CD8
    ('CD3', 'PanCK'),       # T cell + Epithelial (Tumor-Immune) - NEW
]
```

**Key Functions:**
1. `apply_cell_level_qc()` (Line 178-234) - Filter low-quality cells
2. `compute_admixture_score()` (Line 310-356) - Detect doublets/multiplets
3. `compute_wnn_integration()` (Line 237-307) - Simplified WNN
4. `annotate_lineage()` (Line 359-393) - Assign lineages
5. `annotate_cell_types()` (Line 396-432) - Refine cell types

**WNN Implementation Note (Line 243-258):**
This is a SIMPLIFIED WNN using variance-weighted concatenation, NOT true per-cell WNN. Acceptable given 17 proteins but documented trade-off.

---

### Script 63: Merge and Batch Correct
**File:** `scripts/63_merge_and_batch_correct.py`
**Purpose:** Merge processed samples + batch correction

| Input | Output |
|-------|--------|
| `results/qc_all_samples/final_processed/*.h5ad` | `results/qc_all_samples/merged/merged_counts.h5ad` |
| | `results/qc_all_samples/merged/merged_normalized.h5ad` |
| | `results/qc_all_samples/merged/merged_corrected.h5ad` |
| | `results/qc_all_samples/merged/MERGE_REPORT.md` |

**Key Fix in Commit `00d47c2`:**
- Always saves `merged_counts.h5ad` with TRUE raw counts before normalization
- scVI uses `.layers['counts']` from raw, not normalized data

---

## 4. Code Review Findings Addressed

### Finding 1: Lane 4 Cell Size Anomaly ✅ DETECTION ADDED
**Commit:** `6c201a9`
**Issue:** Lane 4 cells (54-70 µm²) smaller than Lanes 1-3 (80-120 µm²)
**Risk:** Segmentation artifact could confound biological conclusions

**Fix Applied:**
- Added `transcript_density` metric (transcripts/µm²) in Script 61
- Added 3-panel density comparison plots
- Added density analysis section to QC report with z-score assessment

**Note:** This adds DETECTION, not automatic correction. If density is constant across lanes despite area differences, downstream analysis should normalize by area or use density-based metrics.

**How to Verify:**
```python
# In sample_qc_summary.csv, check:
# - median_transcript_density should be similar across lanes
# - If constant density but different area → segmentation artifact (normalize by area)
# - If density differs → potential staining issue (investigate)
```

---

### Finding 2: Admixture Discrepancy ✅ FIXED
**Commit:** `6c201a9`
**Issue:** Missing Tumor-Immune conflict pair

**Fix Applied:**
- Added `('CD3', 'PanCK')` to CONFLICT_PAIRS in Script 62
- Added documentation about heuristic vs cellAdmix approach

---

### Finding 3: Statistical Power ⚠️ DOCUMENTED
**Issue:** Only 2 patients (SNU-105, SNU-107) have complete N→M→C series
**Status:** Not a code fix - downstream analysis consideration
**Mitigation:** Use all partial series in mixed effects models, not just matched

---

### Finding 4: WNN "Lite" ⚠️ DOCUMENTED
**Issue:** Simplified WNN uses global weights, not per-cell weights
**Status:** Acceptable trade-off given 17 proteins
**Location:** Script 62, Line 243-258, documented in docstring

---

## 5. Execution Checklist

### Pre-Execution Verification
- [ ] Raw data exists: `/mnt/x/Choi_Batch_2_Tuesday/`
- [ ] Environment active: `conda activate enact`
- [ ] Disk space >5GB: `df -h ~`
- [ ] Previous outputs backed up (if re-running)

### Execution Commands
```bash
cd ~/g4x-choi-batch2-analysis
conda activate enact

# Step 1: Load all samples (~20 min)
python scripts/60_load_all_samples.py 2>&1 | tee logs/60_loading.log

# Step 2: Comprehensive QC (~30 min)
python scripts/61_comprehensive_qc.py 2>&1 | tee logs/61_qc.log

# Step 3: Process passing samples (~3 hours)
python scripts/62_process_qc_passing.py --parallel 8 2>&1 | tee logs/62_processing.log

# Step 4: Merge and batch correct (~1 hour)
python scripts/63_merge_and_batch_correct.py 2>&1 | tee logs/63_merge.log
```

### Post-Execution Verification
- [ ] `sample_qc_summary.csv` has 32 rows
- [ ] `QC_REPORT.md` generated with density analysis
- [ ] H04 sample marked FAIL (expected: zero transcripts)
- [ ] Total FAIL count ≤ 2 (see decision tree below)
- [ ] `merged_counts.h5ad` exists and has raw counts
- [ ] `merged_corrected.h5ad` exists with batch correction
- [ ] Transcript density within expected range (0.3-0.8 trans/µm²)

### Failure Decision Tree

```
How many samples FAIL QC?
│
├── 0-1 failures → PROCEED (expected: H04 only)
│
├── 2-3 failures → INVESTIGATE
│   ├── Same lane? → Lane-level technical issue
│   ├── Same stage? → Biological concern
│   └── Random? → Review thresholds, may proceed
│
├── 4-5 failures → PAUSE
│   ├── Check raw data integrity
│   ├── Verify loading script
│   └── Consider relaxing thresholds with justification
│
└── 6+ failures → STOP
    ├── Major technical problem
    ├── Do NOT proceed without PI consultation
    └── Document and escalate
```

### Rollback Procedure

If pipeline fails mid-execution:
```bash
# Step 60 failed: No cleanup needed, re-run
# Step 61 failed: Remove partial outputs
rm -rf results/qc_all_samples/figures/per_sample/*
rm -f results/qc_all_samples/sample_qc_summary.csv

# Step 62 failed: Remove processed files
rm -rf results/qc_all_samples/final_processed/*

# Step 63 failed: Remove merged files
rm -rf results/qc_all_samples/merged/*

# Then re-run from failed step (--resume flag available for 61)
```

---

## 6. Key Files for Audit Review

### Must-Read Documents
| File | Purpose | Priority |
|------|---------|----------|
| `G4X_ANALYSIS_PLAN.md` | Master analysis strategy | 1 |
| `docs/FULL_DATASET_QC_IMPLEMENTATION_PLAN.md` | QC implementation details | 1 |
| `docs/QC_DETAILED_CHECKLIST.md` | Step-by-step checklist | 2 |
| `VALIDATION.md` | Validation criteria | 2 |

### Must-Review Scripts
| Script | Lines to Focus | Key Logic |
|--------|----------------|-----------|
| `61_comprehensive_qc.py` | 126-170, 359-470, 692-725 | Density metrics, plots, report |
| `62_process_qc_passing.py` | 152-164, 237-307, 310-356 | Conflict pairs, WNN, admixture |
| `63_merge_and_batch_correct.py` | (entire) | Batch correction, count preservation |

### Must-Check Outputs
| Output | Verify |
|--------|--------|
| `sample_qc_summary.csv` | 32 rows, transcript_density column exists |
| `transcript_density_by_lane.png` | 3-panel plot showing density vs area |
| `QC_REPORT.md` | Density analysis section with lane comparison |
| `processing_results.csv` | n_admixed_raw, n_cells columns |

---

## 7. Known Limitations & Decisions

| Issue | Decision | Rationale |
|-------|----------|-----------|
| H04 sample dead | Exclude | Zero transcripts |
| Simplified WNN | Keep | Only 17 proteins; global weights acceptable |
| Python admixture vs cellAdmix | Python heuristic | Faster; covers known conflicts; note in docs |
| Lane 4 small cells | Monitor via density | Don't exclude unless density anomaly |
| N=2 matched series | Use all partial series | Maximize statistical power |

---

## 8. Git History (Recent Commits)

```
6c201a9 fix(qc): Add transcript density metrics and expanded conflict detection
00d47c2 fix(qc): Always store raw counts regardless of magnitude
5e62a1c docs: Add env check override docs, update output file checkpoints
1fe5cc8 fix(qc): Critical fixes for scVI, file naming, and env check flexibility
282e072 fix(qc): Cap HVGs for targeted panels, relax env check
1a0da30 feat(qc): Add environment checks, resume capability, and batch correction
d00dee3 fix: Address critical QC pipeline issues before full dataset run
8665f60 feat: Add full dataset QC pipeline for all 32 samples
```

---

## 9. Audit Sign-Off Checklist

### Code Quality
- [ ] All functions have docstrings
- [ ] Thresholds are evidence-based and documented
- [ ] Error handling for edge cases (empty files, missing columns)
- [ ] Logging captures key metrics

### Scientific Validity
- [ ] QC thresholds appropriate for G4X platform
- [ ] Batch effect metrics (LISI, silhouette) correctly implemented
- [ ] Admixture detection covers major lineage conflicts
- [ ] WNN limitations documented

### Reproducibility
- [ ] Random seeds set where applicable
- [ ] Input/output paths configurable
- [ ] Resume capability works
- [ ] Disk space monitoring prevents data loss

### Documentation
- [ ] CLAUDE.md updated with current status
- [ ] Execution commands copy-paste ready
- [ ] Post-QC analysis checklist exists

---

## 10. Post-QC Handoff

After successful QC pipeline completion, proceed to advanced analysis phases:

| Phase | Script | Key Outputs |
|-------|--------|-------------|
| **1. Progression Analysis** | `34_progression_analysis.py` | Cell type proportions N→M→C, CellRank streamlines |
| **2. Cell-Cell Communication** | `35_cellcell_communication.py` | LIANA+ L-R, niche-phenotype tensor |
| **3. CAF Network Topology** | `36_caf_subtyping.py` | mCAF/iCAF/apCAF scores, Delaunay centrality |
| **4. TLS Detection** | `40_tls_detection.py` | Persistent homology H1, TLS density |
| **5. Spatial Statistics** | `39_spatial_statistics.py` | Chemokine gradients, Ripley's K |
| **6. Multi-Modal PCA** | `41_multimodal_pca.py` | Pseudobulk trajectories |

**Reference:** See `G4X_ANALYSIS_PLAN.md` for detailed methodology.

**Input for Advanced Analysis:**
- Primary: `results/qc_all_samples/merged/merged_corrected.h5ad`
- Fallback: `results/qc_all_samples/merged/merged_normalized.h5ad`

---

**Reviewer Notes:**
_Use this section to add review comments_

