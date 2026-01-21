# Code Provenance

Tracking the origin and authorship of code in this repository.

## Original Code (This Repository)

| File | Description | Author |
|------|-------------|--------|
| scripts/data_loading.py | Unified data loaders for Visium and G4X | Max Van Belkum |
| scripts/preprocess.py | QC filtering pipeline | Max Van Belkum |
| scripts/run_qc_all.py | Batch QC runner | Max Van Belkum |
| scripts/utils/qc_metrics.py | QC metric calculations | Max Van Belkum |
| scripts/utils/plotting.py | QC visualization functions | Max Van Belkum |
| config/samples.yaml | Sample registry | Max Van Belkum |
| config/qc_thresholds.yaml | Platform-specific QC thresholds | Max Van Belkum |

## Adapted Code

| File | Source | Modifications |
|------|--------|---------------|
| scripts/data_loading.py | ~/data/hackathon/notebooks/01_data_loading.py | Extended to all 10 samples, added clinical metadata, added H5 loader for G4X |

## Patterns from Polymath Repos

| Pattern | Source Repo | Usage |
|---------|-------------|-------|
| scanpy_workflow | omicverse/scSLAT | Preprocessing pipeline structure |
| QC thresholds | harpy, stlearn | Platform-specific filtering defaults |
| H5 loading | squidpy | Sparse matrix construction from H5 |

## External Dependencies

- scanpy: Core single-cell analysis
- squidpy: Spatial analysis
- anndata: Data structure
- pandas, numpy: Data manipulation
- matplotlib: Visualization

## Data Sources

| Dataset | Source | Description |
|---------|--------|-------------|
| PDAC Visium | Lab data | 8 samples, treatment response study |
| G4X Multimodal | Lab data | 2 samples, gastric cancer progression |

---
*Last updated: 2026-01-20*
