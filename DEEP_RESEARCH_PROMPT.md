# Deep Research Prompt: G4X Spatial Transcriptomics Analysis

**Purpose:** Comprehensive literature and code resource discovery for Nanostring G4X (GeoMx/CosMx successor) multimodal spatial transcriptomics analysis.

---

## RESEARCH QUERY

I am analyzing a **Nanostring G4X spatial transcriptomics dataset** with the following characteristics:
- **32 tissue samples** from a single study
- **2.3 million cells** with single-cell resolution
- **Multimodal data**: 17 protein markers + full transcriptome RNA
- **Protein panel includes**: CD3, CD4, CD8, FOXP3, PD1, PDL1, CD20, CD68, CD11c, HLA-DR, PanCK, aSMA, CD31, CD45, KI67
- **Data format**: Cell-level CSVs with XY coordinates, protein intensities, transcript counts

I need comprehensive resources covering **BOTH basic workflows AND advanced analyses**.

---

## SECTION 1: Platform-Specific Resources

### 1.1 Nanostring G4X / CosMx / GeoMx Documentation
Find:
- Official Nanostring analysis pipelines and best practices
- G4X vs CosMx vs GeoMx differences and compatibility
- QC metrics specific to imaging-based spatial transcriptomics
- Recommended normalization methods for protein + RNA multimodal data
- Cell segmentation validation approaches

### 1.2 GitHub Repositories - Nanostring Official
Search for:
- `nanostring-biostats` organization repos
- `CosMx-Analysis-Scratch-Space`
- `InSituType` (cell typing)
- `SpatialDecon` (deconvolution)
- `GeomxTools`
- Any G4X-specific analysis repos

### 1.3 Community Analysis Pipelines
Find GitHub repos for:
- End-to-end CosMx/G4X analysis workflows
- Multimodal (protein + RNA) integration pipelines
- Benchmarking imaging spatial transcriptomics platforms

---

## SECTION 2: Foundational Analysis Methods

### 2.1 Cell Type Annotation
Find papers and code for:
- Protein-based cell type gating strategies (like flow cytometry)
- Marker panel optimization for immune cell identification
- Hierarchical gating: CD45+ → CD3+ → CD4/CD8 → FOXP3
- Handling ambiguous/double-positive cells
- Confidence scoring for cell type assignments
- **Specific question**: Best practices for annotating cells with 17 protein markers in spatial context

### 2.2 Quality Control
Find resources on:
- Cell segmentation quality assessment
- Transcript/protein detection thresholds
- Empty cell identification and filtering
- Batch effect detection across lanes/FOVs
- Signal-to-noise ratio calculations for protein markers

### 2.3 Normalization
Find methods for:
- Protein intensity normalization (CLR, z-score, quantile)
- RNA count normalization for imaging-based ST
- Joint normalization of protein + RNA
- Batch correction methods applicable to spatial data (Harmony, ComBat, scVI)

---

## SECTION 3: Spatial Statistics & Analysis

### 3.1 Neighborhood Analysis
Find papers and repos for:
- Neighborhood enrichment analysis (squidpy method)
- Cell-cell interaction inference from spatial proximity
- K-nearest neighbor graph construction for spatial data
- Permutation testing for spatial statistics
- **Tools**: squidpy, Giotto, Seurat v5 spatial

### 3.2 Spatial Variable Gene Detection
Find methods for:
- SpatialDE, SPARK, trendsceek
- Moran's I and Geary's C for spatial autocorrelation
- Hotspot detection algorithms
- Comparing SVG methods on imaging-based ST

### 3.3 Cell-Cell Communication
Find resources on:
- Ligand-receptor analysis in spatial context (CellChat, CellPhoneDB, LIANA)
- Spatial constraints on communication inference
- Receptor-ligand co-localization analysis

---

## SECTION 4: Immune Microenvironment Analysis

### 4.1 Immune Checkpoint Analysis
Find papers on:
- PD1/PDL1 spatial co-expression patterns
- T cell exhaustion signatures in spatial data
- Immune exclusion vs infiltration quantification
- Checkpoint inhibitor response prediction from spatial data

### 4.2 T Cell Subpopulations
Find methods for:
- CD8+ T cell exhaustion scoring
- Treg identification and spatial distribution
- CD8/Treg ratio as prognostic marker
- Spatial organization of T cell subsets

### 4.3 Tertiary Lymphoid Structures (TLS)
Find papers and code for:
- TLS detection algorithms in spatial data
- B cell + T cell aggregation scoring
- TLS maturity classification
- TLS as biomarker for immunotherapy response

### 4.4 Macrophage Polarization
Find resources on:
- M1/M2 classification from protein markers
- Spatial distribution of macrophage subsets
- Macrophage-tumor cell interactions

---

## SECTION 5: Advanced Multimodal Analysis

### 5.1 Protein-RNA Integration
Find methods for:
- Correlating protein and RNA expression at single-cell level
- Multi-omics factor analysis (MOFA, MOFA+)
- Canonical correlation analysis for multimodal data
- Discordance analysis (when protein ≠ RNA)

### 5.2 Deep Learning Approaches
Find papers and repos on:
- Graph neural networks for spatial transcriptomics
- Attention mechanisms for spatial context
- Foundation models for spatial biology (Nicheformer, etc.)
- Self-supervised learning on spatial data

### 5.3 Deconvolution & Cell Type Mapping
Find tools for:
- Cell2location, RCTD, Tangram, CARD
- Reference-based vs reference-free deconvolution
- Deconvolution validation strategies

---

## SECTION 6: Visualization & Interpretation

### 6.1 Spatial Visualization
Find tools and examples for:
- Publication-quality spatial plots
- Interactive visualization (Vitessce, napari, TissUUmaps)
- Multi-panel figure generation
- H&E overlay with molecular data

### 6.2 Dimensionality Reduction
Find methods for:
- UMAP/tSNE with spatial constraints
- Spatial-aware clustering (Leiden with spatial neighbors)
- Trajectory analysis in spatial context

---

## SECTION 7: Reproducibility & Best Practices

### 7.1 Benchmarking Studies
Find papers comparing:
- Imaging ST platforms (CosMx vs Xenium vs MERSCOPE)
- Cell segmentation methods
- Normalization approaches
- Spatial statistics methods

### 7.2 Reporting Standards
Find guidelines on:
- Minimum information for spatial transcriptomics (MIST)
- Reproducible analysis workflows
- Data sharing standards

---

## SECTION 8: Specific Analysis Questions

Please find resources addressing these specific challenges:

1. **Cell segmentation errors**: How do segmentation errors (cell admixture) affect downstream analysis? Methods to detect and correct?

2. **Rare cell types**: How to reliably identify rare populations (<1%) in imaging-based ST?

3. **Spatial heterogeneity**: Methods to quantify tumor microenvironment heterogeneity across samples?

4. **Batch effects**: Best practices for multi-sample/multi-lane normalization?

5. **Statistical power**: How many cells/samples needed for robust spatial statistics?

---

## OUTPUT FORMAT REQUESTED

Please organize findings into:

### A. Key Papers (with PMIDs/DOIs)
| Paper Title | Year | Key Contribution | PMID/DOI |
|-------------|------|------------------|----------|

### B. GitHub Repositories
| Repo | Stars | Language | Use Case |
|------|-------|----------|----------|

### C. Software Tools
| Tool | Platform | Purpose | Documentation URL |
|------|----------|---------|-------------------|

### D. Tutorials & Workflows
| Resource | Type | Platform Coverage |
|----------|------|-------------------|

### E. Benchmark Datasets
| Dataset | Platform | Size | Access |
|---------|----------|------|--------|

---

## PRIORITY AREAS

**HIGH PRIORITY** (need immediately):
1. Nanostring G4X/CosMx official analysis pipelines
2. Multimodal protein+RNA integration methods
3. Immune checkpoint spatial analysis papers
4. Cell segmentation error detection/correction

**MEDIUM PRIORITY**:
5. TLS detection algorithms
6. Spatial statistics benchmarking
7. Deep learning for spatial biology

**LOWER PRIORITY** (nice to have):
8. Interactive visualization tools
9. Foundation models for spatial data
10. Advanced trajectory analysis

---

## SEARCH TERMS TO USE

```
"G4X spatial transcriptomics"
"CosMx SMI analysis"
"Nanostring spatial proteomics"
"imaging-based spatial transcriptomics"
"multimodal spatial omics protein RNA"
"cell segmentation spatial transcriptomics"
"neighborhood enrichment spatial"
"immune checkpoint spatial analysis"
"PD1 PDL1 spatial co-expression"
"tertiary lymphoid structure detection"
"T cell exhaustion spatial transcriptomics"
"squidpy tutorial"
"Giotto spatial analysis"
"cell2location deconvolution"
"spatial variable genes detection"
```

---

## DATABASES TO SEARCH

- **PubMed/PMC**: Recent papers (2022-2026)
- **bioRxiv/medRxiv**: Preprints
- **GitHub**: Code repositories
- **Zenodo**: Datasets and workflows
- **protocols.io**: Analysis protocols
- **YouTube**: Tutorial videos
- **Nature Methods/Protocols**: Methods papers

---

*This prompt is designed for use with deep research tools like Perplexity Pro, ChatGPT with browsing, Claude with web search, or similar AI research assistants.*
