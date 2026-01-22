# Next Session: Add QC Data from Master

## Quick Task
Pull QC data from existing report - **don't recompute**.

**Source:** `https://github.com/vanbelkummax/g4x-choi-batch2-analysis/blob/master/reports/G4X_QC_REPORT.md`

## Pilot Sample QC (from master report)

| Sample | Stage | Cells | Median Trans | Median Genes | Status |
|--------|-------|-------|--------------|--------------|--------|
| **E02** | Normal | 31,790 | 51 | 33 | PASS |
| **F02** | Metaplasia | 33,857 | 56 | 35 | PASS |
| **G02** | Cancer | 79,634 | 89 | 50 | PASS |

**Total: 145,281 cells** (same patient, same lane L002)

## To Do
1. `git show origin/master:reports/G4X_QC_REPORT.md > QC.md`
2. Extract E02/F02/G02 sections
3. Add to this branch
4. Done

## Branch
`pilot-clean` on https://github.com/vanbelkummax/g4x-choi-batch2-analysis
