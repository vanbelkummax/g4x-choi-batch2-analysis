# Continue Session: G4X Hypothesis Testing

**Resume:** `cd ~/g4x-choi-batch2-analysis`

## Session Summary (2026-01-20)

### Completed
1. ✅ **Polymath Assessment** - G4X coverage: 1 paper, CosMx: 25, Xenium: 83
2. ✅ **32 Sample Profiles** - 8-panel dashboards in `results/figures/sample_profiles/`
3. ✅ **Group Comparisons** - Statistical tests completed

### Key Statistical Findings

**Lane Effects (Batch):**
- Only PDL1+ tumor shows lane effect (p=0.008) - monitor but not blocking

**CD8/Treg Ratio Groups (n=16 each):**
| Metric | p-value | Effect |
|--------|---------|--------|
| Treg % | 0.002 | 0.65 |
| PD1+ T cells | 0.019 | 0.49 |
| PDL1+ tumor | 0.007 | -0.56 |

**PDL1 Status Groups:**
| Metric | p-value | Effect |
|--------|---------|--------|
| CD8/Treg ratio | 0.002 | -0.66 |
| PD1+ T cells | 0.009 | 0.55 |
| Treg % | 0.014 | 0.52 |

**Exhaustion Groups (most significant!):**
- n_cells: p=0.0009, effect=-0.70
- Treg %: p=0.0002, effect=-0.79
- CD8/Treg: p=0.004, effect=0.59

## Next Steps

### 1. Hypothesis Testing Script (scripts/08_hypothesis_testing.py)
Test the 5 specific hypotheses:
- H1: CD8/Treg affects spatial organization (neighborhood enrichment)
- H2: PDL1+ tumors show immune exclusion (spatial proximity)
- H3: Lane effects minimal ✅ (mostly confirmed)
- H4: Exhausted T cells near PDL1+ cells (co-localization)
- H5: High immune = TLS structures (cluster detection)

### 2. Polymath Enrichment
Papers to retrieve for context:
- PDAC spatial papers (925 OA available)
- Immune checkpoint spatial papers (have 10+)
- TLS detection methods papers

### 3. GitHub Push
```bash
cd ~/g4x-choi-batch2-analysis
git add -A
git commit -m "Add sample profiles and group comparisons"
git push
```

## Key Files

| File | Description |
|------|-------------|
| `results/figures/sample_profiles/*.png` | 32 sample dashboards |
| `results/figures/group_comparisons.png` | Group stat comparisons |
| `results/figures/correlation_matrix.png` | Metric correlations |
| `results/tables/group_comparison_results.csv` | All p-values |
| `results/tables/samples_with_groups.csv` | Samples with group labels |

## Quick Start

```bash
cd ~/g4x-choi-batch2-analysis
source ~/miniforge3/etc/profile.d/conda.sh
conda activate enact

# View results
ls results/figures/sample_profiles/ | head
cat results/tables/group_comparison_results.csv | head -20

# Continue with hypothesis testing
# Create scripts/08_hypothesis_testing.py
```

## Data Summary

- **32 samples**, 2.3M cells total
- **17 protein markers** including PD1, PDL1, CD8, FOXP3
- **Groups defined:** CD8/Treg, PDL1, Exhaustion (balanced 16/16)
- **Strongest signal:** Exhaustion groups show most differences
