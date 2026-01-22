# G4X Block5 Pilot Analysis

**3 samples from gastric cancer progression (same patient)**

| Sample | Stage | Cells |
|--------|-------|-------|
| E02 | Normal | 31,790 |
| F02 | Metaplasia | 33,857 |
| G02 | Cancer | 79,634 |

## Run

```bash
conda activate enact
cd ~/g4x-choi-batch2-analysis

python scripts/01_hvg_optimize.py    # ~5 min
python scripts/02_annotate.py        # ~10 min
python scripts/03_spatial_viz.py     # ~2 min
```

## Data

- **Input:** `results/qc_all_samples/merged/merged_corrected.h5ad`
- **Output:** `results/pilot/`
