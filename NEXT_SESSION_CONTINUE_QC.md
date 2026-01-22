# G4X QC Pipeline - Continue Session

**Last Updated:** 2026-01-21 20:44

## Current State

**Step 4 (Merge & Batch Correct) is RUNNING in background**

Task ID: `b5280a8`
Log: `/home/user/g4x-choi-batch2-analysis/logs/63_merge_v2.log`

### Pipeline Progress

| Step | Script | Status | Notes |
|------|--------|--------|-------|
| 1 | 60_load_all_samples.py | ✅ Complete | 32 samples, 2.3M cells, protein fixed |
| 2 | 61_comprehensive_qc.py | ✅ Complete | 30 PASS, 1 WARN (D01), 1 FAIL (H04) |
| 3 | 62_process_qc_passing.py | ✅ Complete | 29 samples processed (C02, H02 failed NaN) |
| 4 | 63_merge_and_batch_correct.py | 🔄 Running | PCA + Harmony on 1.8M cells |

### Key Fixes Applied This Session

1. **Protein loading bug** (script 60): Changed from feature_matrix.h5 to cell_by_protein.csv.gz
2. **Python 3.9 compatibility** (script 63): Changed `dict | None` to `Optional[dict]`
3. **PyTorch tensor handling** (script 63): Fixed Harmony output conversion for CUDA

### To Check If Step 4 Completed

```bash
# Check if process is still running
ps aux | grep 63_merge | grep python

# Check output
tail -50 /home/user/g4x-choi-batch2-analysis/logs/63_merge_v2.log

# Check output files
ls -la /home/user/g4x-choi-batch2-analysis/results/qc_all_samples/merged/
```

### Expected Output Files

```
results/qc_all_samples/merged/
├── merged_counts.h5ad        # Raw counts (6.4GB) - CREATED
├── merged_normalized.h5ad    # Pre-correction - PENDING
├── merged_corrected.h5ad     # Post-Harmony - PENDING
├── batch_assessment.csv      # LISI/silhouette - PENDING
└── MERGE_REPORT.md           # Summary - PENDING
```

### If Step 4 Failed

If Harmony still has issues, re-run with:
```bash
conda activate enact
python scripts/63_merge_and_batch_correct.py --method none 2>&1 | tee logs/63_merge_skip_correction.log
```

This skips batch correction and just merges.

### If Step 4 Succeeded

Commit the final changes:
```bash
git add -A && git commit -m "feat: Complete QC pipeline with batch correction"
git push
```

### Post-QC Analysis Checklist

After QC pipeline completes:

1. **34_progression_analysis.py** - Cell type proportions N→M→C
2. **35_cellcell_communication.py** - LIANA+ ligand-receptor analysis
3. **36_caf_subtyping.py** - mCAF/iCAF/apCAF scoring
4. **40_tls_detection.py** - TLS via persistent homology
5. **39_spatial_statistics.py** - Ripley's K, gradient fields

### Summary Stats

- **Total cells loaded:** 2,308,968
- **Cells after cell-level QC:** 1,835,026
- **Samples processed:** 29/32
- **Failed samples:** H04 (low quality), C02 & H02 (NaN during PCA)
- **Batch correction:** Harmony (LISI=2.46 triggered correction)

### Commits This Session

- `f1f083c` - fix: Load protein data from cell_by_protein.csv.gz
