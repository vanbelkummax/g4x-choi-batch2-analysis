#!/usr/bin/env python3
"""
Day 3: Polymathic Analysis - Cross-Domain Algorithms
=====================================================

Apply algorithms from topology, graph theory, and information theory
to spatial transcriptomics data, as identified through Polymath KB.

Methods:
1. Graph Centrality (from Graph Theory) - Hub cell identification
2. Persistent Homology (from Topology) - Tissue architecture
3. Mutual Information (from Information Theory) - Gene selection

Author: Max Van Belkum
Date: 2026-01-19
"""

import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import warnings
import logging

warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

# Check for optional TDA imports
try:
    from gtda.homology import VietorisRipsPersistence
    from gtda.diagrams import BettiCurve
    HAS_TDA = True
except ImportError:
    HAS_TDA = False
    logger.warning("giotto-tda not available. Install with: pip install giotto-tda")

# =============================================================================
# Configuration
# =============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
ADATA_DIR = OUTPUT_DIR / "adata" / "annotated"
POLYMATHIC_DIR = OUTPUT_DIR / "adata" / "polymathic"
FIG_DIR = OUTPUT_DIR / "figures" / "polymathic"
TABLE_DIR = OUTPUT_DIR / "tables"

POLYMATHIC_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Sample metadata (PDAC only)
PDAC_METADATA = {
    "YP03A": {"patient": "YP03", "timepoint": "Pre", "response": "NR"},
    "YP03C": {"patient": "YP03", "timepoint": "Post", "response": "NR"},
    "YP04C": {"patient": "YP04", "timepoint": "Post", "response": "NR"},  # YP04A excluded
    "YP12A": {"patient": "YP12", "timepoint": "Pre", "response": "R"},
    "YP12C": {"patient": "YP12", "timepoint": "Post", "response": "R"},
    "YP15A": {"patient": "YP15", "timepoint": "Pre", "response": "R"},
    "YP15C": {"patient": "YP15", "timepoint": "Post", "response": "R"},
}


# =============================================================================
# Graph Centrality Analysis (from Graph Theory)
# =============================================================================

def compute_graph_centrality(adata, n_neighs: int = 6) -> Dict[str, float]:
    """
    Compute hub cell statistics using graph centrality measures.

    Cross-domain application: Graph centrality algorithms from social network
    analysis applied to spatial cell interaction networks.

    Returns dict with centrality statistics.
    """
    logger.info("Computing spatial graph centrality...")

    # Ensure spatial graph exists
    if 'spatial_connectivities' not in adata.obsp:
        sq.gr.spatial_neighbors(adata, n_neighs=n_neighs)

    # Convert to networkx
    G = nx.from_scipy_sparse_array(adata.obsp['spatial_connectivities'])

    # Compute centrality measures
    logger.info("  - Betweenness centrality...")
    betweenness = nx.betweenness_centrality(G)

    logger.info("  - PageRank...")
    pagerank = nx.pagerank(G)

    logger.info("  - Degree centrality...")
    degree = nx.degree_centrality(G)

    # Store in adata
    adata.obs['betweenness'] = [betweenness.get(i, 0) for i in range(adata.n_obs)]
    adata.obs['pagerank'] = [pagerank.get(i, 0) for i in range(adata.n_obs)]
    adata.obs['degree_centrality'] = [degree.get(i, 0) for i in range(adata.n_obs)]

    # Identify hub cells (top 5% by pagerank)
    threshold = np.percentile(adata.obs['pagerank'], 95)
    adata.obs['is_hub_cell'] = adata.obs['pagerank'] > threshold

    # Calculate statistics
    stats = {
        'n_cells': adata.n_obs,
        'n_hub_cells': adata.obs['is_hub_cell'].sum(),
        'hub_fraction': adata.obs['is_hub_cell'].mean(),
        'mean_betweenness': adata.obs['betweenness'].mean(),
        'max_betweenness': adata.obs['betweenness'].max(),
        'mean_pagerank': adata.obs['pagerank'].mean(),
        'max_pagerank': adata.obs['pagerank'].max(),
        'mean_degree': adata.obs['degree_centrality'].mean(),
    }

    # Hub cell type distribution
    if 'cell_type' in adata.obs.columns:
        hub_cells = adata[adata.obs['is_hub_cell']]
        hub_type_counts = hub_cells.obs['cell_type'].value_counts()
        stats['dominant_hub_type'] = hub_type_counts.index[0] if len(hub_type_counts) > 0 else 'Unknown'

    return stats


# =============================================================================
# Persistent Homology Analysis (from Topology)
# =============================================================================

def maxmin_subsample(coords: np.ndarray, n_landmarks: int, seed: int = 42) -> np.ndarray:
    """
    MaxMin (farthest point) landmark selection for topology-preserving subsampling.

    Preserves topological features by selecting points that maximize coverage
    of the point cloud shape, rather than random selection which can break loops.

    Algorithm:
    1. Start with random point
    2. Repeatedly add the point farthest from all current landmarks
    3. Continue until n_landmarks reached

    References:
    - De Silva & Carlsson (2004) "Topological estimation using witness complexes"
    - Polymath KB: Search "landmark selection topology"

    Args:
        coords: Array of shape (n_points, n_dims) with coordinates
        n_landmarks: Number of landmark points to select
        seed: Random seed for reproducibility

    Returns:
        Array of indices for selected landmark points
    """
    rng = np.random.default_rng(seed)
    n = len(coords)

    if n <= n_landmarks:
        return np.arange(n)

    # Start with random point
    landmarks = [rng.integers(n)]

    # Distance to nearest landmark for each point
    min_dists = np.full(n, np.inf)

    for _ in range(n_landmarks - 1):
        # Update distances to nearest landmark
        last_landmark = coords[landmarks[-1]]
        dists_to_last = np.linalg.norm(coords - last_landmark, axis=1)
        min_dists = np.minimum(min_dists, dists_to_last)

        # Add farthest point
        next_landmark = np.argmax(min_dists)
        landmarks.append(next_landmark)

    return np.array(landmarks)


def compute_topology_features(adata, max_cells: int = 5000) -> Optional[Dict[str, float]]:
    """
    Compute persistent homology features for tissue architecture.

    Cross-domain application: Topological data analysis from algebraic topology
    applied to quantify tissue structure complexity.

    H0 = connected components (tissue clusters)
    H1 = holes/voids (tissue gaps/boundaries)

    Uses MaxMin landmark selection (not random) to preserve topological features.
    """
    if not HAS_TDA:
        logger.warning("Skipping topology - giotto-tda not installed")
        return None

    logger.info("Computing persistent homology...")

    coords = adata.obsm['spatial']

    # Subsample if too large (TDA is O(n^3))
    # Use MaxMin landmark selection to preserve topology
    if len(coords) > max_cells:
        idx = maxmin_subsample(coords, max_cells, seed=42)
        coords = coords[idx]
        logger.info(f"  - MaxMin subsampled to {max_cells} cells (preserves topology)")

    # Normalize coordinates
    coords = (coords - coords.min(axis=0)) / (coords.max(axis=0) - coords.min(axis=0) + 1e-6)

    # Compute Vietoris-Rips persistence
    VR = VietorisRipsPersistence(
        metric='euclidean',
        homology_dimensions=[0, 1],
        n_jobs=-1
    )
    diagrams = VR.fit_transform([coords])

    # Compute Betti curves
    betti = BettiCurve()
    features = betti.fit_transform(diagrams)

    stats = {
        'betti_0_mean': features[0, :, 0].mean(),  # Connected components
        'betti_0_max': features[0, :, 0].max(),
        'betti_1_mean': features[0, :, 1].mean(),  # Holes/voids
        'betti_1_max': features[0, :, 1].max(),
        'topology_complexity': features[0, :, 0].mean() + features[0, :, 1].mean(),
    }

    return stats


# =============================================================================
# Mutual Information Analysis (from Information Theory)
# =============================================================================

def compute_cross_sample_mi_response_pseudobulk(all_adatas: Dict, metadata: Dict, top_n: int = 50) -> pd.DataFrame:
    """
    Compute MI at SAMPLE level (correct statistical unit) to avoid pseudoreplication.

    CRITICAL: The original cell-level approach treated 32K cells as independent,
    but cells from the same sample share biological and technical variation.
    This pseudobulk approach aggregates to sample means first, giving n=7 observations.

    Args:
        all_adatas: Dict of sample_name -> AnnData objects
        metadata: Dict of sample_name -> {response: 'R' or 'NR', ...}
        top_n: Number of top genes to return

    Returns:
        DataFrame with genes sorted by MI against response, with Welch's t-test validation
    """
    from sklearn.feature_selection import mutual_info_classif
    from scipy.stats import ttest_ind
    import anndata as ad

    logger.info("Computing cross-sample MI using PSEUDOBULK aggregation (correct statistical unit)...")

    # Step 1: Find common genes across all samples
    common_genes = None
    for sample, adata in all_adatas.items():
        if sample not in metadata:
            continue
        if common_genes is None:
            common_genes = set(adata.var_names)
        else:
            common_genes = common_genes.intersection(set(adata.var_names))

    if not common_genes:
        logger.warning("  - No common genes found")
        return pd.DataFrame()

    gene_names = sorted(list(common_genes))
    logger.info(f"  - Found {len(gene_names)} common genes across samples")

    # Step 2: Aggregate expression by sample (pseudobulk = mean expression per gene)
    sample_expr = {}
    sample_labels = {}

    for sample, adata in all_adatas.items():
        if sample not in metadata:
            continue
        # Subset to common genes and compute mean expression per gene
        adata_subset = adata[:, gene_names]
        X = adata_subset.X.toarray() if hasattr(adata_subset.X, 'toarray') else adata_subset.X
        sample_expr[sample] = X.mean(axis=0)  # Mean across all cells in sample
        sample_labels[sample] = 1 if metadata[sample]['response'] == 'R' else 0

    # Step 3: Build sample-level matrices (n=7 samples x n_genes)
    samples = sorted(sample_expr.keys())
    X = np.array([sample_expr[s] for s in samples])  # (7, n_genes)
    y = np.array([sample_labels[s] for s in samples])  # (7,)

    logger.info(f"  - Pseudobulk matrix: {X.shape[0]} samples x {X.shape[1]} genes")
    logger.info(f"  - Response labels: {sum(y)} R, {len(y) - sum(y)} NR")

    # Step 4: Limit genes for MI computation (top variable genes)
    n_genes_limit = min(5000, X.shape[1])
    if X.shape[1] > n_genes_limit:
        gene_var = np.var(X, axis=0)
        top_var_idx = np.argsort(gene_var)[-n_genes_limit:]
        X = X[:, top_var_idx]
        gene_names = [gene_names[i] for i in top_var_idx]
        logger.info(f"  - Using top {n_genes_limit} variable genes")

    # Step 5: Compute MI (now on 7 INDEPENDENT samples)
    # Use n_neighbors=2 since we have very few samples
    logger.info(f"  - Computing MI for {len(gene_names)} genes vs response (sample-level)...")
    mi_scores = mutual_info_classif(X, y, random_state=42, n_neighbors=2)

    # Step 6: Create results dataframe with validation via Welch's t-test
    results = pd.DataFrame({
        'gene': gene_names,
        'mi_vs_response': mi_scores
    }).sort_values('mi_vs_response', ascending=False)

    # Step 7: Add Welch's t-test p-values for top candidates (statistical validation)
    logger.info("  - Validating top candidates with Welch's t-test...")
    welch_pvals = []
    r_means = []
    nr_means = []

    for gene in results['gene']:
        gene_idx = gene_names.index(gene)
        r_vals = [sample_expr[s][gene_idx] for s in samples if sample_labels[s] == 1]
        nr_vals = [sample_expr[s][gene_idx] for s in samples if sample_labels[s] == 0]

        if len(r_vals) >= 2 and len(nr_vals) >= 2:
            _, welch_p = ttest_ind(r_vals, nr_vals, equal_var=False)
            welch_pvals.append(welch_p)
            r_means.append(np.mean(r_vals))
            nr_means.append(np.mean(nr_vals))
        else:
            welch_pvals.append(np.nan)
            r_means.append(np.nan)
            nr_means.append(np.nan)

    results['welch_p'] = welch_pvals
    results['r_mean'] = r_means
    results['nr_mean'] = nr_means
    results['fold_change'] = results['r_mean'] / (results['nr_mean'] + 1e-10)

    logger.info(f"  - Top 5 by MI (with Welch validation):")
    for _, row in results.head(5).iterrows():
        logger.info(f"    {row['gene']:15s} MI={row['mi_vs_response']:.4f}, Welch p={row['welch_p']:.4f}")

    return results.head(top_n)


def compute_cross_sample_mi_response(all_adatas: Dict, metadata: Dict, top_n: int = 50) -> pd.DataFrame:
    """
    DEPRECATED: Use compute_cross_sample_mi_response_pseudobulk instead.

    This cell-level approach has pseudoreplication issues (treating 32K cells as independent).
    Kept for reference/comparison only.
    """
    logger.warning("DEPRECATED: Using cell-level MI (pseudoreplication issue). Use pseudobulk version.")
    from sklearn.feature_selection import mutual_info_classif
    import anndata as ad

    logger.info("Computing cross-sample MI against TREATMENT RESPONSE (DEPRECATED cell-level)...")

    # Merge all samples with response labels
    adatas_list = []
    for sample, adata in all_adatas.items():
        if sample not in metadata:
            continue
        adata_copy = adata.copy()
        adata_copy.obs['sample'] = sample
        adata_copy.obs['response'] = metadata[sample]['response']
        adatas_list.append(adata_copy)

    if not adatas_list:
        logger.warning("  - No samples to merge")
        return pd.DataFrame()

    # Concatenate (inner join on genes)
    merged = ad.concat(adatas_list, join='inner')
    logger.info(f"  - Merged {len(adatas_list)} samples: {merged.n_obs} cells, {merged.n_vars} genes")

    # Get expression and response labels
    X = merged.X.toarray() if hasattr(merged.X, 'toarray') else merged.X
    y = (merged.obs['response'] == 'R').astype(int).values  # 1=R, 0=NR

    # Limit genes for speed
    n_genes = min(5000, X.shape[1])
    if X.shape[1] > n_genes:
        gene_var = np.var(X, axis=0)
        top_var_idx = np.argsort(gene_var)[-n_genes:]
        X = X[:, top_var_idx]
        gene_names = merged.var_names[top_var_idx]
    else:
        gene_names = merged.var_names

    logger.info(f"  - Computing MI for {len(gene_names)} genes vs response...")
    mi_scores = mutual_info_classif(X, y, random_state=42)

    mi_df = pd.DataFrame({
        'gene': gene_names,
        'mi_vs_response': mi_scores
    }).sort_values('mi_vs_response', ascending=False)

    return mi_df.head(top_n)


def compute_mutual_information_genes(adata, target_col: str = 'cell_type', top_n: int = 50) -> pd.DataFrame:
    """
    Compute mutual information between genes and a target variable.

    Cross-domain application: Information theory from signal processing/statistics
    applied to identify non-linearly associated genes.

    MI captures dependencies that correlation misses.
    """
    from sklearn.feature_selection import mutual_info_classif

    logger.info(f"Computing mutual information for genes vs {target_col}...")

    if target_col not in adata.obs.columns:
        logger.warning(f"  - {target_col} not found, skipping MI analysis")
        return pd.DataFrame()

    # Get expression matrix and target
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    y = adata.obs[target_col].astype('category').cat.codes.values

    # Remove cells with unknown labels
    mask = y >= 0
    X = X[mask]
    y = y[mask]

    # Compute MI (this can be slow for many genes)
    n_genes = min(5000, X.shape[1])  # Limit for speed
    if X.shape[1] > n_genes:
        # Use most variable genes
        gene_var = np.var(X, axis=0)
        top_var_idx = np.argsort(gene_var)[-n_genes:]
        X = X[:, top_var_idx]
        gene_names = adata.var_names[top_var_idx]
    else:
        gene_names = adata.var_names

    logger.info(f"  - Computing MI for {len(gene_names)} genes...")
    mi_scores = mutual_info_classif(X, y, random_state=42)

    # Create results DataFrame
    mi_df = pd.DataFrame({
        'gene': gene_names,
        'mutual_information': mi_scores
    }).sort_values('mutual_information', ascending=False)

    return mi_df.head(top_n)


# =============================================================================
# Main Analysis
# =============================================================================

def analyze_sample(sample_name: str, adata_path: Path) -> Dict:
    """Run full polymathic analysis on a single sample."""
    logger.info(f"\n{'='*60}")
    logger.info(f"Analyzing {sample_name}")
    logger.info(f"{'='*60}")

    adata = sc.read_h5ad(adata_path)
    logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # Get metadata
    meta = PDAC_METADATA.get(sample_name, {})

    results = {
        'sample': sample_name,
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        **meta
    }

    # 1. Graph Centrality
    centrality_stats = compute_graph_centrality(adata)
    results.update({f'centrality_{k}': v for k, v in centrality_stats.items()})

    # 2. Topology Features
    topology_stats = compute_topology_features(adata)
    if topology_stats:
        results.update({f'topology_{k}': v for k, v in topology_stats.items()})

    # 3. Mutual Information (only for cell type)
    mi_df = compute_mutual_information_genes(adata, 'cell_type')
    if len(mi_df) > 0:
        results['top_mi_gene'] = mi_df.iloc[0]['gene']
        results['top_mi_score'] = mi_df.iloc[0]['mutual_information']

        # Save MI results
        mi_path = TABLE_DIR / f"{sample_name}_mi_genes.csv"
        mi_df.to_csv(mi_path, index=False)

    # Save updated adata
    out_path = POLYMATHIC_DIR / f"{sample_name}_polymathic.h5ad"
    adata.write_h5ad(out_path)
    logger.info(f"Saved to {out_path}")

    return results


def main():
    """Run polymathic analysis on all PDAC samples."""
    logger.info("="*60)
    logger.info("DAY 3: POLYMATHIC ANALYSIS")
    logger.info("Cross-domain algorithms for spatial transcriptomics")
    logger.info("="*60)

    all_results = []
    all_adatas = {}  # Keep track for cross-sample MI

    # Analyze PDAC samples
    for sample_name in PDAC_METADATA.keys():
        adata_path = ADATA_DIR / f"{sample_name}_annotated.h5ad"

        if not adata_path.exists():
            logger.warning(f"Skipping {sample_name} - file not found")
            continue

        try:
            results = analyze_sample(sample_name, adata_path)
            all_results.append(results)
            # Load for cross-sample MI
            all_adatas[sample_name] = sc.read_h5ad(adata_path)
        except Exception as e:
            logger.error(f"Error analyzing {sample_name}: {e}")
            continue

    # Compute cross-sample MI against treatment response using PSEUDOBULK (correct statistical unit)
    if len(all_adatas) >= 2:
        logger.info("\n" + "="*60)
        logger.info("CROSS-SAMPLE MI: Treatment Response Biomarkers (Pseudobulk)")
        logger.info("="*60)
        mi_response_df = compute_cross_sample_mi_response_pseudobulk(all_adatas, PDAC_METADATA, top_n=50)
        if len(mi_response_df) > 0:
            mi_response_path = TABLE_DIR / "mi_vs_response_biomarkers.csv"
            mi_response_df.to_csv(mi_response_path, index=False)
            logger.info(f"Saved MI vs Response (Pseudobulk) to {mi_response_path}")
            logger.info(f"Top 10 Treatment Response Biomarkers (by MI, with Welch validation):")
            for _, row in mi_response_df.head(10).iterrows():
                welch_str = f", Welch p={row['welch_p']:.4f}" if 'welch_p' in row else ""
                logger.info(f"  {row['gene']:15s} MI={row['mi_vs_response']:.4f}{welch_str}")

    # Combine results
    if all_results:
        results_df = pd.DataFrame(all_results)
        results_path = TABLE_DIR / "polymathic_analysis_results.csv"
        results_df.to_csv(results_path, index=False)
        logger.info(f"\nSaved combined results to {results_path}")

        # Print summary comparison
        logger.info("\n" + "="*60)
        logger.info("R vs NR COMPARISON (Polymathic Features)")
        logger.info("="*60)

        r_samples = results_df[results_df['response'] == 'R']
        nr_samples = results_df[results_df['response'] == 'NR']

        metrics = ['centrality_hub_fraction', 'centrality_mean_pagerank',
                   'centrality_mean_betweenness']

        if HAS_TDA:
            metrics.extend(['topology_betti_0_mean', 'topology_betti_1_mean',
                           'topology_complexity'])

        for metric in metrics:
            if metric in results_df.columns:
                r_mean = r_samples[metric].mean()
                nr_mean = nr_samples[metric].mean()
                diff_pct = ((r_mean - nr_mean) / (nr_mean + 1e-6)) * 100
                logger.info(f"  {metric:35s}: R={r_mean:.4f}, NR={nr_mean:.4f} ({diff_pct:+.1f}%)")

    logger.info("\n" + "="*60)
    logger.info("ANALYSIS COMPLETE")
    logger.info("="*60)


if __name__ == '__main__':
    main()
