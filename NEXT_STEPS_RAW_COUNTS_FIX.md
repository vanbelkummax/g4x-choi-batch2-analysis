# Next Steps: Raw Counts Fix Execution Plan

**Status:** Ready to Execute
**Date:** 2026-01-22
**Estimated Time:** 3-4 hours total

---

## Pre-Execution Checklist

- [ ] Read `CRITICAL_AUDIT_RAW_COUNTS_BUG.md` to understand the issue
- [ ] Verify disk space (need ~20GB free for new merged files)
- [ ] Ensure `enact` conda environment is active

---

## Execution Steps

### Step 1: Validate Current State (5 min)

First, run validation to confirm the bug exists:

```bash
cd ~/g4x-choi-batch2-analysis
conda activate enact

# Validate current state - should show the bug
python scripts/66_validate_raw_counts.py --verbose
```

**Expected output:** Should confirm raw files are valid but processed/merged files have the bug.

---

### Step 2: Fix Per-Sample Processing (30-60 min)

Re-process all samples with raw counts preservation:

```bash
# Process all samples (parallel mode)
python scripts/64_fix_raw_counts_preservation.py --parallel 8 2>&1 | tee logs/64_fix.log

# Or process a single sample first to test
python scripts/64_fix_raw_counts_preservation.py --sample A01
```

**Output:** Creates `results/qc_all_samples/final_processed_fixed/` with properly preserved counts.

---

### Step 3: Rebuild Merged Files (2-3 hours)

Rebuild all merged files from fixed processed files:

```bash
# Rebuild with Harmony batch correction
python scripts/65_rebuild_merged_from_raw.py --method harmony 2>&1 | tee logs/65_rebuild.log

# OR rebuild directly from raw files (alternative path)
python scripts/65_rebuild_merged_from_raw.py --from-raw --method harmony 2>&1 | tee logs/65_rebuild_raw.log
```

**Output:** Creates `results/qc_all_samples/merged_fixed/` with:
- `merged_counts.h5ad` - TRUE raw counts in X
- `merged_normalized.h5ad` - normalized X, raw in layers['counts']
- `merged_corrected.h5ad` - batch-corrected, raw in layers['counts']

---

### Step 4: Validate Fix (5 min)

Verify the fix worked:

```bash
python scripts/66_validate_raw_counts.py --fix-check --verbose
```

**Expected output:** All validations should PASS.

---

### Step 5: Update Analysis to Use Fixed Files

After validation passes, update downstream scripts to use the fixed merged files:

```python
# In analysis scripts, change from:
adata = sc.read_h5ad("results/qc_all_samples/merged/merged_corrected.h5ad")

# To:
adata = sc.read_h5ad("results/qc_all_samples/merged_fixed/merged_corrected.h5ad")
```

---

## Quick Validation Commands

```bash
# Quick check of merged_counts.h5ad
python -c "
import anndata as ad
import numpy as np
a = ad.read_h5ad('results/qc_all_samples/merged_fixed/merged_counts.h5ad', backed='r')
X = a.X[:1000, :100]
if hasattr(X, 'toarray'): X = X.toarray()
print(f'Max: {X.max():.0f} (should be >50)')
print(f'Integer-like: {np.allclose(X, np.round(X))} (should be True)')
"

# Quick check of merged_corrected.h5ad layers
python -c "
import anndata as ad
import numpy as np
a = ad.read_h5ad('results/qc_all_samples/merged_fixed/merged_corrected.h5ad', backed='r')
L = a.layers['counts'][:1000, :100]
if hasattr(L, 'toarray'): L = L.toarray()
print(f'layers[counts] max: {L.max():.0f} (should be >50)')
print(f'Integer-like: {np.allclose(L, np.round(L))} (should be True)')
"
```

---

## Troubleshooting

### Issue: Script 64 fails on a sample
```bash
# Check the specific sample
python scripts/64_fix_raw_counts_preservation.py --sample SAMPLE_ID --verbose
```

### Issue: Not enough disk space
```bash
# Check disk space
df -h .

# Remove old files if needed (AFTER validation)
rm -rf results/qc_all_samples/final_processed/  # Old buggy files
```

### Issue: scVI fails
The scripts default to Harmony which doesn't require raw counts. If you need scVI:
```bash
# Verify raw counts are available
python scripts/66_validate_raw_counts.py --fix-check

# Then run with scVI (after fix)
python scripts/65_rebuild_merged_from_raw.py --method scvi
```

---

## Files Created by Fix

```
results/qc_all_samples/
├── final_processed_fixed/     # NEW: Fixed per-sample files
│   ├── A01_final_fixed.h5ad
│   ├── B01_final_fixed.h5ad
│   └── ...
├── merged_fixed/              # NEW: Fixed merged files
│   ├── merged_counts.h5ad     # TRUE raw counts
│   ├── merged_normalized.h5ad # Normalized + raw in layers
│   ├── merged_corrected.h5ad  # Batch-corrected + raw in layers
│   └── VALIDATION_REPORT.md
└── RAW_COUNTS_VALIDATION_REPORT.md
```

---

## Post-Fix: Within-Patient Analysis

After fixing raw counts, proceed with the recommended within-patient analysis:

```bash
# SNU-105 and SNU-107 each have Normal → Metaplasia → Cancer
# This eliminates patient confounding completely

python scripts/49_within_patient_progression.py 2>&1 | tee logs/49_within_patient.log
```

---

## Contact

For issues with this fix, check:
1. `logs/64_fix.log` - per-sample processing logs
2. `logs/65_rebuild.log` - merge logs
3. `CRITICAL_AUDIT_RAW_COUNTS_BUG.md` - original audit report
