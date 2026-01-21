#!/usr/bin/env python3
"""
G4X WNN Integration with Segmentation Artifact Control
=======================================================
Spatial Hackathon 2026 - Max Van Belkum

Weighted integration of RNA and Protein modalities + admixture scoring.
Simplified implementation that works directly on AnnData (no MuData).

Usage:
    python 30_wnn_integration.py [--sample SAMPLE_ID] [--all] [--parallel N]
"""

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from scipy import sparse, stats
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import warnings
import logging
import gc

warnings.filterwarnings('ignore')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# =============================================================================
# Configuration
# =============================================================================

INPUT_DIR = Path("/home/user/spatial-hackathon-2026/results/g4x_choi_batch2")
OUTPUT_DIR = Path("/home/user/spatial-hackathon-2026/results/g4x_choi_batch2/wnn_integrated")
FIGURE_DIR = Path("/home/user/spatial-hackathon-2026/figures/g4x/wnn")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# Marker conflict definitions for admixture detection
CONFLICT_PAIRS = [
    ('PanCK', 'CD45'),      # Epithelial + Immune
    ('PanCK', 'aSMA'),      # Epithelial + Stromal
    ('CD45', 'aSMA'),       # Immune + Stromal
    ('CD31', 'CD45'),       # Endothelial + Immune
    ('CD3', 'CD68'),        # T cell + Macrophage
    ('CD4', 'CD8'),         # CD4 + CD8
]

LINEAGE_MARKERS = {
    'epithelial': ['PanCK'],
    'immune': ['CD45'],
    'stromal': ['aSMA'],
    'endothelial': ['CD31'],
}

ADMIX_THRESHOLD = 0.7  # Percentile for "high" expression
ADMIX_CUTOFF = 0.3     # Flag cells above this score


# =============================================================================
# Weighted Nearest Neighbors (Direct Implementation)
# =============================================================================

def compute_weighted_integration(adata, n_pcs_rna=30, n_pcs_prot=15, n_neighbors=20):
    """
    Compute weighted integration of RNA and protein modalities.

    Returns AnnData with:
    - obsm['X_wnn']: Combined weighted representation
    - obs['rna_weight']: Per-modality weight (global)
    - obs['protein_weight']: Per-modality weight (global)
    """
    logger.info("Computing weighted multimodal integration...")

    # Get RNA PCA
    if 'X_pca' not in adata.obsm:
        logger.info("  Computing RNA PCA...")
        sc.pp.pca(adata, n_comps=min(n_pcs_rna, adata.n_vars - 1))

    rna_pca = adata.obsm['X_pca'][:, :min(n_pcs_rna, adata.obsm['X_pca'].shape[1])].copy()

    # Get protein data and compute PCA
    logger.info("  Processing protein modality...")
    protein_names = adata.uns.get('protein_names', [])
    protein_data = adata.obsm['protein'].copy()

    if sparse.issparse(protein_data):
        protein_data = protein_data.toarray()

    # Z-score normalize protein
    scaler = StandardScaler()
    protein_scaled = scaler.fit_transform(protein_data)

    # PCA on protein
    from sklearn.decomposition import PCA
    n_prot_comps = min(n_pcs_prot, protein_data.shape[1] - 1)
    pca_prot = PCA(n_components=n_prot_comps)
    protein_pca = pca_prot.fit_transform(protein_scaled)
    logger.info(f"  Protein PCA: {protein_pca.shape[1]} components, {pca_prot.explained_variance_ratio_.sum():.1%} variance")

    # Store protein PCA
    adata.obsm['X_pca_protein'] = protein_pca

    # Compute variance-based weights
    rna_var = float(np.var(rna_pca, axis=0).sum())
    prot_var = float(np.var(protein_pca, axis=0).sum())

    total_var = rna_var + prot_var + 1e-10
    rna_weight = prot_var / total_var  # Inverse variance weighting
    prot_weight = rna_var / total_var

    logger.info(f"  Modality weights - RNA: {rna_weight:.3f}, Protein: {prot_weight:.3f}")

    # Z-score normalize each modality
    rna_norm = StandardScaler().fit_transform(rna_pca)
    prot_norm = StandardScaler().fit_transform(protein_pca)

    # Weighted concatenation
    combined = np.hstack([
        rna_norm * np.sqrt(rna_weight),
        prot_norm * np.sqrt(prot_weight)
    ])

    # Store results
    adata.obsm['X_wnn'] = combined
    adata.obs['rna_weight'] = rna_weight
    adata.obs['protein_weight'] = prot_weight

    logger.info(f"  Combined representation: {combined.shape}")

    return adata


def compute_wnn_clustering(adata, n_neighbors=20, resolutions=[0.5, 1.0]):
    """
    Compute neighbors, clustering, and UMAP on WNN representation.
    """
    logger.info("Computing WNN clustering and embedding...")

    # Build neighbor graph on combined representation
    logger.info("  Building neighbors...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_wnn', key_added='wnn')

    # Leiden clustering at multiple resolutions
    for res in resolutions:
        key = f'leiden_wnn_{str(res).replace(".", "_")}'
        logger.info(f"  Leiden clustering (res={res})...")
        sc.tl.leiden(adata, resolution=res, key_added=key, neighbors_key='wnn')
        n_clusters = adata.obs[key].nunique()
        logger.info(f"    Found {n_clusters} clusters")

    # UMAP
    logger.info("  Computing UMAP...")
    sc.tl.umap(adata, neighbors_key='wnn')
    adata.obsm['X_wnn_umap'] = adata.obsm['X_umap'].copy()

    return adata


# =============================================================================
# Admixture/Contamination Scoring
# =============================================================================

def compute_lineage_scores(adata):
    """Compute lineage probability scores from protein markers."""
    protein_names = list(adata.uns.get('protein_names', []))
    protein_data = adata.obsm['protein']

    if sparse.issparse(protein_data):
        protein_data = protein_data.toarray()

    name_to_idx = {p: i for i, p in enumerate(protein_names)}

    scores = {}
    for lineage, markers in LINEAGE_MARKERS.items():
        indices = [name_to_idx[m] for m in markers if m in name_to_idx]
        if indices:
            scores[lineage] = np.mean(protein_data[:, indices], axis=1)
        else:
            scores[lineage] = np.zeros(adata.n_obs)

    # Normalize to 0-1
    df = pd.DataFrame(scores, index=adata.obs_names)
    for col in df.columns:
        vmin, vmax = np.percentile(df[col], [1, 99])
        df[col] = np.clip((df[col] - vmin) / (vmax - vmin + 1e-10), 0, 1)

    # Softmax for probabilities
    exp_scores = np.exp(df.values * 3)
    probs = exp_scores / exp_scores.sum(axis=1, keepdims=True)

    return pd.DataFrame(probs, columns=[f'p_{c}' for c in df.columns], index=adata.obs_names)


def compute_marker_conflicts(adata):
    """Detect marker conflicts indicating admixture."""
    protein_names = list(adata.uns.get('protein_names', []))
    protein_data = adata.obsm['protein']

    if sparse.issparse(protein_data):
        protein_data = protein_data.toarray()

    name_to_idx = {p: i for i, p in enumerate(protein_names)}

    # Compute thresholds
    thresholds = {}
    for name, idx in name_to_idx.items():
        thresholds[name] = np.percentile(protein_data[:, idx], ADMIX_THRESHOLD * 100)

    # Count conflicts
    n_conflicts = np.zeros(adata.n_obs)
    conflict_scores = np.zeros(adata.n_obs)

    for m1, m2 in CONFLICT_PAIRS:
        if m1 not in name_to_idx or m2 not in name_to_idx:
            continue

        expr1 = protein_data[:, name_to_idx[m1]]
        expr2 = protein_data[:, name_to_idx[m2]]

        # Binary conflict
        conflict = (expr1 > thresholds[m1]) & (expr2 > thresholds[m2])
        n_conflicts += conflict.astype(int)

        # Continuous score
        norm1 = np.clip((expr1 - thresholds[m1]) / (expr1.max() - thresholds[m1] + 1e-10), 0, 1)
        norm2 = np.clip((expr2 - thresholds[m2]) / (expr2.max() - thresholds[m2] + 1e-10), 0, 1)
        conflict_scores += np.minimum(norm1, norm2)

    return pd.DataFrame({
        'n_conflicts': n_conflicts,
        'conflict_score': conflict_scores / max(len(CONFLICT_PAIRS), 1)
    }, index=adata.obs_names)


def compute_neighborhood_consistency(adata, n_neighbors=10):
    """Compute spatial neighborhood consistency."""
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    elif 'cell_x' in adata.obs and 'cell_y' in adata.obs:
        coords = adata.obs[['cell_x', 'cell_y']].values
    else:
        return pd.DataFrame({'neighborhood_consistency': np.ones(adata.n_obs)}, index=adata.obs_names)

    # Get cell types
    if 'cell_type' in adata.obs:
        cell_types = adata.obs['cell_type'].values
    else:
        # Use dominant lineage
        probs = compute_lineage_scores(adata)
        lineages = ['epithelial', 'immune', 'stromal', 'endothelial']
        cell_types = np.array([lineages[i] for i in np.argmax(probs.values, axis=1)])

    # Build spatial neighbors
    nn = NearestNeighbors(n_neighbors=min(n_neighbors + 1, adata.n_obs))
    nn.fit(coords)
    _, indices = nn.kneighbors(coords)

    # Consistency score
    consistency = np.zeros(adata.n_obs)
    for i in range(adata.n_obs):
        neighbor_idx = indices[i, 1:]  # Exclude self
        consistency[i] = np.mean([cell_types[j] == cell_types[i] for j in neighbor_idx])

    return pd.DataFrame({'neighborhood_consistency': consistency}, index=adata.obs_names)


def compute_admixture_score(adata):
    """Compute composite admixture/contamination score."""
    logger.info("Computing admixture scores...")

    conflicts = compute_marker_conflicts(adata)
    consistency = compute_neighborhood_consistency(adata)
    probs = compute_lineage_scores(adata)

    # Entropy from lineage probabilities
    probs_arr = np.clip(probs.values, 1e-10, 1)
    entropy = -np.sum(probs_arr * np.log(probs_arr), axis=1)
    max_entropy = -np.log(1 / probs_arr.shape[1])
    norm_entropy = entropy / max_entropy

    # Composite score
    admix_score = (
        0.5 * conflicts['conflict_score'].values +
        0.3 * (1 - consistency['neighborhood_consistency'].values) +
        0.2 * norm_entropy
    )

    result = pd.DataFrame({
        'admixture_score': admix_score,
        'admixture_flag': admix_score > ADMIX_CUTOFF,
        'conflict_score': conflicts['conflict_score'].values,
        'neighborhood_consistency': consistency['neighborhood_consistency'].values,
        'lineage_entropy': norm_entropy,
    }, index=adata.obs_names)

    # Add lineage probabilities
    for col in probs.columns:
        result[col] = probs[col].values

    n_flagged = result['admixture_flag'].sum()
    pct_flagged = 100 * n_flagged / len(result)
    logger.info(f"  Flagged {n_flagged} cells ({pct_flagged:.1f}%) as potential admixture")

    return result


# =============================================================================
# Visualization
# =============================================================================

def plot_wnn_comparison(adata, output_path):
    """Compare RNA-only, WNN UMAP, and spatial admixture."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # RNA-only UMAP
    if 'X_umap' in adata.obsm and 'leiden' in adata.obs:
        sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, title='RNA-only')
    else:
        axes[0].set_title('RNA-only (not available)')

    # WNN UMAP
    if 'X_wnn_umap' in adata.obsm:
        adata_temp = adata.copy()
        adata_temp.obsm['X_umap'] = adata.obsm['X_wnn_umap']
        color_key = 'leiden_wnn_0_5' if 'leiden_wnn_0_5' in adata.obs else 'leiden'
        sc.pl.umap(adata_temp, color=color_key, ax=axes[1], show=False, title='WNN Integrated')

    # Admixture score spatial
    if 'spatial' in adata.obsm and 'admixture_score' in adata.obs:
        coords = adata.obsm['spatial']
        n = min(20000, adata.n_obs)
        idx = np.random.choice(adata.n_obs, n, replace=False)
        scatter = axes[2].scatter(
            coords[idx, 0], coords[idx, 1],
            c=adata.obs['admixture_score'].values[idx],
            cmap='RdYlBu_r', s=1, alpha=0.5
        )
        axes[2].set_title('Spatial Admixture Score')
        axes[2].set_aspect('equal')
        plt.colorbar(scatter, ax=axes[2])

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_admixture_summary(adata, output_path):
    """Plot admixture score distribution."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    admix = adata.obs['admixture_score']

    # Distribution
    axes[0, 0].hist(admix, bins=50, edgecolor='black', alpha=0.7)
    axes[0, 0].axvline(ADMIX_CUTOFF, color='red', linestyle='--', label=f'Cutoff ({ADMIX_CUTOFF})')
    axes[0, 0].set_xlabel('Admixture Score')
    axes[0, 0].set_ylabel('Count')
    axes[0, 0].set_title('Admixture Score Distribution')
    axes[0, 0].legend()

    # Conflict vs consistency
    scatter = axes[0, 1].scatter(
        adata.obs['conflict_score'],
        adata.obs['neighborhood_consistency'],
        c=admix, cmap='RdYlBu_r', alpha=0.3, s=1
    )
    axes[0, 1].set_xlabel('Conflict Score')
    axes[0, 1].set_ylabel('Neighborhood Consistency')
    axes[0, 1].set_title('Conflict vs Consistency')
    plt.colorbar(scatter, ax=axes[0, 1])

    # Clean vs flagged
    flagged = adata.obs['admixture_flag']
    axes[1, 0].bar(['Clean', 'Flagged'], [sum(~flagged), sum(flagged)],
                   color=['steelblue', 'salmon'])
    axes[1, 0].set_ylabel('Number of Cells')
    axes[1, 0].set_title(f'Clean vs Flagged ({100*flagged.mean():.1f}% flagged)')

    # Lineage distribution in flagged
    lineage_cols = [c for c in adata.obs.columns if c.startswith('p_')]
    if lineage_cols and flagged.any():
        means = adata.obs.loc[flagged, lineage_cols].mean()
        axes[1, 1].bar(means.index, means.values, color='salmon')
        axes[1, 1].set_ylabel('Mean Probability')
        axes[1, 1].set_title('Lineage Distribution (Flagged Cells)')
        axes[1, 1].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


# =============================================================================
# Main Pipeline
# =============================================================================

def process_sample(sample_path):
    """Process single sample through WNN + admixture scoring."""
    sample_id = sample_path.stem.replace('_analyzed', '')

    try:
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing: {sample_id}")
        logger.info(f"{'='*60}")

        # Load
        adata = sc.read_h5ad(sample_path)
        logger.info(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

        if 'protein' not in adata.obsm:
            logger.error(f"  No protein data, skipping")
            return None

        # WNN integration
        adata = compute_weighted_integration(adata)
        adata = compute_wnn_clustering(adata)

        # Admixture scoring
        admix_df = compute_admixture_score(adata)
        for col in admix_df.columns:
            adata.obs[col] = admix_df[col].values

        # Save
        out_path = OUTPUT_DIR / f"{sample_id}_wnn.h5ad"
        adata.write(out_path)
        logger.info(f"  Saved: {out_path}")

        # Figures
        plot_wnn_comparison(adata, FIGURE_DIR / f"{sample_id}_wnn_comparison.png")
        plot_admixture_summary(adata, FIGURE_DIR / f"{sample_id}_admixture.png")

        # Summary
        result = {
            'sample_id': sample_id,
            'n_cells': adata.n_obs,
            'n_flagged': int(admix_df['admixture_flag'].sum()),
            'pct_flagged': float(100 * admix_df['admixture_flag'].mean()),
            'mean_admix': float(admix_df['admixture_score'].mean()),
            'rna_weight': float(adata.obs['rna_weight'].iloc[0]),
            'protein_weight': float(adata.obs['protein_weight'].iloc[0]),
            'n_clusters': int(adata.obs['leiden_wnn_0_5'].nunique()) if 'leiden_wnn_0_5' in adata.obs else 0,
        }

        # Clean up
        del adata
        gc.collect()

        return result

    except Exception as e:
        logger.error(f"  Error processing {sample_id}: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    parser = argparse.ArgumentParser(description='G4X WNN Integration')
    parser.add_argument('--sample', type=str, help='Process specific sample')
    parser.add_argument('--all', action='store_true', help='Process all')
    parser.add_argument('--parallel', type=int, default=1, help='Parallel workers')
    args = parser.parse_args()

    print("=" * 70)
    print("G4X WNN Integration + Admixture Scoring")
    print("=" * 70)

    sample_files = sorted(INPUT_DIR.glob("*_analyzed.h5ad"))
    logger.info(f"Found {len(sample_files)} samples")

    if args.sample:
        sample_files = [f for f in sample_files if args.sample in f.stem]

    results = []

    if args.parallel > 1 and len(sample_files) > 1:
        logger.info(f"Processing in parallel with {args.parallel} workers")
        with ProcessPoolExecutor(max_workers=args.parallel) as executor:
            futures = {executor.submit(process_sample, f): f for f in sample_files}
            for future in tqdm(as_completed(futures), total=len(futures)):
                result = future.result()
                if result:
                    results.append(result)
    else:
        for f in tqdm(sample_files, desc="Processing"):
            result = process_sample(f)
            if result:
                results.append(result)

    # Summary
    if results:
        df = pd.DataFrame(results)
        df.to_csv(OUTPUT_DIR / "wnn_summary.csv", index=False)
        print("\n" + "=" * 70)
        print("Summary")
        print("=" * 70)
        print(df.to_string(index=False))
        print(f"\nMean flagged: {df['pct_flagged'].mean():.1f}%")

    print("\nDone!")


if __name__ == "__main__":
    main()
