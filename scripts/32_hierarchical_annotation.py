#!/usr/bin/env python3
"""
G4X Hierarchical Cell Type Annotation
======================================
Spatial Hackathon 2026 - Max Van Belkum

Clustering-first annotation with uncertainty propagation.
Replaces hard gating to reduce "Unknown" cells.

Strategy:
1. First-pass: Protein-based lineage compartments (Epithelial/Immune/Stromal/Endothelial)
2. Second-pass: Leiden clustering within each compartment
3. Third-pass: Annotate clusters using marker profiles
4. Output: Cell type + confidence score + entropy

Usage:
    python 32_hierarchical_annotation.py [--sample SAMPLE_ID] [--all]
"""

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import sparse, stats
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
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

INPUT_DIR = Path("/home/user/spatial-hackathon-2026/results/g4x_choi_batch2/wnn_integrated")
OUTPUT_DIR = Path("/home/user/spatial-hackathon-2026/results/g4x_choi_batch2/annotated_v2")
FIGURE_DIR = Path("/home/user/spatial-hackathon-2026/figures/g4x/annotation")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# Marker Definitions
# =============================================================================

# Lineage markers (protein-based)
LINEAGE_MARKERS = {
    'Epithelial': {'positive': ['PanCK'], 'negative': ['CD45', 'aSMA', 'CD31']},
    'Immune': {'positive': ['CD45'], 'negative': ['PanCK']},
    'Stromal': {'positive': ['aSMA'], 'negative': ['CD45', 'PanCK']},
    'Endothelial': {'positive': ['CD31'], 'negative': ['CD45', 'PanCK']},
}

# Immune subtype markers
IMMUNE_MARKERS = {
    'T_cell': {'positive': ['CD3'], 'negative': ['CD20', 'CD68']},
    'CD4_T': {'positive': ['CD3', 'CD4'], 'negative': ['CD8']},
    'CD8_T': {'positive': ['CD3', 'CD8'], 'negative': ['CD4']},
    'Treg': {'positive': ['CD3', 'CD4', 'FOXP3'], 'negative': ['CD8']},
    'B_cell': {'positive': ['CD20'], 'negative': ['CD3']},
    'Macrophage': {'positive': ['CD68'], 'negative': ['CD3', 'CD20']},
    'Dendritic': {'positive': ['CD11c', 'HLA-DR'], 'negative': ['CD68']},
}

# Functional markers
FUNCTIONAL_MARKERS = {
    'Exhausted': ['PD1'],
    'Proliferating': ['KI67'],
    'Checkpoint_ligand': ['PDL1'],
}

# Confidence thresholds
CONFIDENCE_THRESHOLD = 0.6  # Below this = "Low_confidence"


# =============================================================================
# Annotation Functions
# =============================================================================

def get_protein_data(adata):
    """Extract protein expression data and names."""
    protein_names = list(adata.uns.get('protein_names', []))
    protein_data = adata.obsm['protein'].copy()

    if sparse.issparse(protein_data):
        protein_data = protein_data.toarray()

    name_to_idx = {p: i for i, p in enumerate(protein_names)}

    return protein_data, protein_names, name_to_idx


def compute_marker_scores(protein_data, name_to_idx, marker_dict, n_cells):
    """Compute positive/negative marker scores."""
    scores = {}

    for cell_type, markers in marker_dict.items():
        pos_markers = markers.get('positive', [])
        neg_markers = markers.get('negative', [])

        # Positive marker score (mean)
        pos_idx = [name_to_idx[m] for m in pos_markers if m in name_to_idx]
        if pos_idx:
            pos_score = np.mean(protein_data[:, pos_idx], axis=1)
        else:
            pos_score = np.zeros(n_cells)

        # Negative marker score (should be low)
        neg_idx = [name_to_idx[m] for m in neg_markers if m in name_to_idx]
        if neg_idx:
            neg_score = np.mean(protein_data[:, neg_idx], axis=1)
        else:
            neg_score = np.zeros(n_cells)

        # Combined score: high positive, low negative
        scores[cell_type] = pos_score - 0.5 * neg_score

    return scores


def assign_lineages(adata):
    """
    Assign cells to lineage compartments using probabilistic scoring.
    """
    logger.info("Assigning lineage compartments...")

    protein_data, protein_names, name_to_idx = get_protein_data(adata)
    n_cells = adata.n_obs

    # Compute lineage scores
    scores = compute_marker_scores(protein_data, name_to_idx, LINEAGE_MARKERS, n_cells)

    # Convert to DataFrame and normalize
    score_df = pd.DataFrame(scores, index=adata.obs_names)

    # Percentile normalization per lineage
    for col in score_df.columns:
        vmin, vmax = np.percentile(score_df[col], [5, 95])
        score_df[col] = np.clip((score_df[col] - vmin) / (vmax - vmin + 1e-10), 0, 1)

    # Softmax for probabilities
    temperature = 2.0  # Lower = more confident assignments
    exp_scores = np.exp(score_df.values / temperature)
    probs = exp_scores / exp_scores.sum(axis=1, keepdims=True)
    prob_df = pd.DataFrame(probs, columns=[f'p_{c}' for c in score_df.columns], index=adata.obs_names)

    # Assign lineage based on max probability
    lineages = score_df.columns[np.argmax(score_df.values, axis=1)]
    confidences = np.max(probs, axis=1)

    # Compute entropy (uncertainty)
    entropy = -np.sum(probs * np.log(probs + 1e-10), axis=1)
    max_entropy = -np.log(1 / probs.shape[1])
    normalized_entropy = entropy / max_entropy

    # Store results
    adata.obs['lineage'] = pd.Categorical(lineages)
    adata.obs['lineage_confidence'] = confidences
    adata.obs['lineage_entropy'] = normalized_entropy

    for col in prob_df.columns:
        adata.obs[col] = prob_df[col].values

    # Summary
    for lin in score_df.columns:
        count = np.sum(lineages == lin)
        mean_conf = confidences[lineages == lin].mean()
        logger.info(f"  {lin}: {count} cells ({100*count/n_cells:.1f}%), mean confidence: {mean_conf:.2f}")

    return adata


def cluster_within_lineage(adata, lineage, resolution=0.5, min_cells=100):
    """
    Perform Leiden clustering within a lineage compartment.
    """
    mask = adata.obs['lineage'] == lineage
    n_cells = np.sum(mask)

    if n_cells < min_cells:
        logger.warning(f"  {lineage}: Too few cells ({n_cells}), skipping subclustering")
        return None, None

    adata_sub = adata[mask].copy()

    # Use WNN representation if available
    if 'X_wnn' in adata_sub.obsm:
        use_rep = 'X_wnn'
    elif 'X_pca' in adata_sub.obsm:
        use_rep = 'X_pca'
    else:
        logger.warning(f"  {lineage}: No suitable representation found")
        return None, None

    # Build neighbors and cluster
    try:
        sc.pp.neighbors(adata_sub, n_neighbors=15, use_rep=use_rep)
        sc.tl.leiden(adata_sub, resolution=resolution, key_added='subcluster')

        n_clusters = adata_sub.obs['subcluster'].nunique()
        logger.info(f"  {lineage}: {n_cells} cells → {n_clusters} subclusters")

        return adata_sub, mask

    except Exception as e:
        logger.error(f"  {lineage} clustering failed: {e}")
        return None, None


def annotate_immune_clusters(adata_sub):
    """
    Annotate immune subclusters based on marker profiles.
    """
    protein_data, protein_names, name_to_idx = get_protein_data(adata_sub)

    # Compute immune subtype scores
    scores = compute_marker_scores(protein_data, name_to_idx, IMMUNE_MARKERS, adata_sub.n_obs)
    score_df = pd.DataFrame(scores, index=adata_sub.obs_names)

    # Normalize
    for col in score_df.columns:
        vmin, vmax = np.percentile(score_df[col], [5, 95])
        score_df[col] = np.clip((score_df[col] - vmin) / (vmax - vmin + 1e-10), 0, 1)

    # Softmax probabilities
    exp_scores = np.exp(score_df.values * 2)
    probs = exp_scores / exp_scores.sum(axis=1, keepdims=True)

    # Per-cluster assignment (majority voting with confidence)
    cluster_annotations = {}
    for cluster in adata_sub.obs['subcluster'].unique():
        mask = adata_sub.obs['subcluster'] == cluster
        cluster_probs = probs[mask].mean(axis=0)
        best_type = score_df.columns[np.argmax(cluster_probs)]
        confidence = np.max(cluster_probs)

        # Add functional modifiers
        modifiers = []
        cluster_protein = protein_data[mask].mean(axis=0)

        for modifier, markers in FUNCTIONAL_MARKERS.items():
            marker_idx = [name_to_idx[m] for m in markers if m in name_to_idx]
            if marker_idx:
                modifier_expr = np.mean(cluster_protein[marker_idx])
                threshold = np.percentile(protein_data[:, marker_idx].flatten(), 70)
                if modifier_expr > threshold:
                    modifiers.append(modifier)

        # Build annotation
        if confidence < CONFIDENCE_THRESHOLD:
            annotation = f"Low_confidence_{best_type}"
        else:
            annotation = best_type

        if modifiers:
            annotation = f"{annotation}_{'+'.join(modifiers)}"

        cluster_annotations[cluster] = (annotation, confidence)

    return cluster_annotations, score_df


def annotate_epithelial_clusters(adata_sub):
    """
    Annotate epithelial subclusters (simpler - mainly KI67/proliferation status).
    """
    protein_data, protein_names, name_to_idx = get_protein_data(adata_sub)

    cluster_annotations = {}
    for cluster in adata_sub.obs['subcluster'].unique():
        mask = adata_sub.obs['subcluster'] == cluster
        cluster_protein = protein_data[mask].mean(axis=0)

        # Check proliferation
        if 'KI67' in name_to_idx:
            ki67_expr = cluster_protein[name_to_idx['KI67']]
            ki67_threshold = np.percentile(protein_data[:, name_to_idx['KI67']], 70)
            if ki67_expr > ki67_threshold:
                annotation = "Epithelial_proliferating"
            else:
                annotation = "Epithelial"
        else:
            annotation = "Epithelial"

        cluster_annotations[cluster] = (annotation, 0.8)  # High confidence for epithelial

    return cluster_annotations


def annotate_stromal_clusters(adata_sub):
    """Annotate stromal subclusters."""
    cluster_annotations = {}
    for cluster in adata_sub.obs['subcluster'].unique():
        cluster_annotations[cluster] = ("Fibroblast", 0.7)
    return cluster_annotations


def annotate_endothelial_clusters(adata_sub):
    """Annotate endothelial subclusters."""
    cluster_annotations = {}
    for cluster in adata_sub.obs['subcluster'].unique():
        cluster_annotations[cluster] = ("Endothelial", 0.7)
    return cluster_annotations


def hierarchical_annotation(adata):
    """
    Main hierarchical annotation pipeline.
    """
    logger.info("Starting hierarchical annotation...")

    # Step 1: Assign lineages
    adata = assign_lineages(adata)

    # Initialize cell type columns
    adata.obs['cell_type'] = 'Unknown'
    adata.obs['cell_type_confidence'] = 0.0
    adata.obs['cell_type_detailed'] = 'Unknown'

    # Step 2: Cluster and annotate within each lineage
    annotation_funcs = {
        'Immune': annotate_immune_clusters,
        'Epithelial': annotate_epithelial_clusters,
        'Stromal': annotate_stromal_clusters,
        'Endothelial': annotate_endothelial_clusters,
    }

    for lineage in adata.obs['lineage'].cat.categories:
        logger.info(f"\nProcessing {lineage}...")

        adata_sub, mask = cluster_within_lineage(adata, lineage)

        if adata_sub is None:
            # Assign lineage as cell type for small compartments
            adata.obs.loc[adata.obs['lineage'] == lineage, 'cell_type'] = lineage
            adata.obs.loc[adata.obs['lineage'] == lineage, 'cell_type_detailed'] = lineage
            adata.obs.loc[adata.obs['lineage'] == lineage, 'cell_type_confidence'] = \
                adata.obs.loc[adata.obs['lineage'] == lineage, 'lineage_confidence']
            continue

        # Get annotation function
        annotate_func = annotation_funcs.get(lineage)

        if annotate_func is not None:
            if lineage == 'Immune':
                cluster_annotations, score_df = annotate_func(adata_sub)
            else:
                cluster_annotations = annotate_func(adata_sub)
        else:
            # Default: use lineage as type
            cluster_annotations = {c: (lineage, 0.7) for c in adata_sub.obs['subcluster'].unique()}

        # Map annotations back to original adata
        cell_indices = adata_sub.obs_names
        for i, cell_id in enumerate(cell_indices):
            cluster = adata_sub.obs.loc[cell_id, 'subcluster']
            annotation, confidence = cluster_annotations.get(cluster, (lineage, 0.5))

            # Combine lineage and cluster confidence
            lineage_conf = adata.obs.loc[cell_id, 'lineage_confidence']
            final_conf = 0.7 * confidence + 0.3 * lineage_conf

            adata.obs.loc[cell_id, 'cell_type_detailed'] = annotation
            adata.obs.loc[cell_id, 'cell_type_confidence'] = final_conf

            # Simplified cell type (without modifiers)
            base_type = annotation.split('_')[0] if '_' in annotation else annotation
            if base_type.startswith('Low'):
                base_type = annotation.split('_')[-1]  # Extract actual type from Low_confidence_X
            adata.obs.loc[cell_id, 'cell_type'] = base_type

    # Summary
    logger.info("\nCell Type Summary:")
    for ct in adata.obs['cell_type'].value_counts().index[:15]:
        count = np.sum(adata.obs['cell_type'] == ct)
        mean_conf = adata.obs.loc[adata.obs['cell_type'] == ct, 'cell_type_confidence'].mean()
        logger.info(f"  {ct}: {count} ({100*count/adata.n_obs:.1f}%), confidence: {mean_conf:.2f}")

    # Count low confidence
    low_conf = np.sum(adata.obs['cell_type_confidence'] < CONFIDENCE_THRESHOLD)
    logger.info(f"\nLow confidence cells: {low_conf} ({100*low_conf/adata.n_obs:.1f}%)")

    return adata


# =============================================================================
# Visualization
# =============================================================================

def plot_annotation_summary(adata, output_path):
    """Plot cell type annotation summary."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Cell type proportions
    ct_counts = adata.obs['cell_type'].value_counts()
    axes[0, 0].barh(ct_counts.index[:15], ct_counts.values[:15])
    axes[0, 0].set_xlabel('Number of Cells')
    axes[0, 0].set_title('Cell Type Distribution')
    axes[0, 0].invert_yaxis()

    # Confidence distribution
    axes[0, 1].hist(adata.obs['cell_type_confidence'], bins=50, edgecolor='black', alpha=0.7)
    axes[0, 1].axvline(CONFIDENCE_THRESHOLD, color='red', linestyle='--', label=f'Threshold ({CONFIDENCE_THRESHOLD})')
    axes[0, 1].set_xlabel('Cell Type Confidence')
    axes[0, 1].set_ylabel('Count')
    axes[0, 1].set_title('Confidence Distribution')
    axes[0, 1].legend()

    # Lineage breakdown
    lineage_counts = adata.obs['lineage'].value_counts()
    axes[1, 0].pie(lineage_counts.values, labels=lineage_counts.index, autopct='%1.1f%%')
    axes[1, 0].set_title('Lineage Compartments')

    # Entropy distribution
    axes[1, 1].hist(adata.obs['lineage_entropy'], bins=50, edgecolor='black', alpha=0.7, color='orange')
    axes[1, 1].set_xlabel('Lineage Entropy')
    axes[1, 1].set_ylabel('Count')
    axes[1, 1].set_title('Classification Uncertainty')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_spatial_annotation(adata, output_path):
    """Plot spatial distribution of cell types."""
    if 'spatial' not in adata.obsm:
        return

    coords = adata.obsm['spatial']

    # Get top cell types
    top_types = adata.obs['cell_type'].value_counts().index[:8]
    n_types = len(top_types)

    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    # Subsample for large datasets
    n = min(30000, adata.n_obs)
    idx = np.random.choice(adata.n_obs, n, replace=False)

    for i, ct in enumerate(top_types):
        if i >= len(axes):
            break

        mask = adata.obs['cell_type'].values == ct

        # Background
        axes[i].scatter(coords[idx, 0], coords[idx, 1], c='lightgray', s=0.5, alpha=0.3)

        # Highlighted type
        type_idx = np.where(mask)[0]
        type_idx = np.intersect1d(type_idx, idx)
        if len(type_idx) > 0:
            axes[i].scatter(coords[type_idx, 0], coords[type_idx, 1], c='red', s=1, alpha=0.7)

        axes[i].set_title(f'{ct}\n({np.sum(mask)} cells)')
        axes[i].set_aspect('equal')
        axes[i].axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


# =============================================================================
# Main Pipeline
# =============================================================================

def process_sample(sample_path):
    """Process single sample through hierarchical annotation."""
    sample_id = sample_path.stem.replace('_wnn', '')

    try:
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing: {sample_id}")
        logger.info(f"{'='*60}")

        # Load WNN-integrated data
        adata = sc.read_h5ad(sample_path)
        logger.info(f"  Loaded: {adata.n_obs} cells")

        # Run hierarchical annotation
        adata = hierarchical_annotation(adata)

        # Save
        out_path = OUTPUT_DIR / f"{sample_id}_annotated.h5ad"
        adata.write(out_path)
        logger.info(f"  Saved: {out_path}")

        # Figures
        plot_annotation_summary(adata, FIGURE_DIR / f"{sample_id}_annotation_summary.png")
        plot_spatial_annotation(adata, FIGURE_DIR / f"{sample_id}_spatial_annotation.png")

        # Summary
        result = {
            'sample_id': sample_id,
            'n_cells': adata.n_obs,
            'n_cell_types': adata.obs['cell_type'].nunique(),
            'mean_confidence': float(adata.obs['cell_type_confidence'].mean()),
            'pct_low_confidence': float(100 * np.mean(adata.obs['cell_type_confidence'] < CONFIDENCE_THRESHOLD)),
            'pct_epithelial': float(100 * np.mean(adata.obs['lineage'] == 'Epithelial')),
            'pct_immune': float(100 * np.mean(adata.obs['lineage'] == 'Immune')),
            'pct_stromal': float(100 * np.mean(adata.obs['lineage'] == 'Stromal')),
            'pct_endothelial': float(100 * np.mean(adata.obs['lineage'] == 'Endothelial')),
        }

        del adata
        gc.collect()

        return result

    except Exception as e:
        logger.error(f"  Error processing {sample_id}: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    parser = argparse.ArgumentParser(description='G4X Hierarchical Annotation')
    parser.add_argument('--sample', type=str, help='Process specific sample')
    parser.add_argument('--all', action='store_true', help='Process all')
    args = parser.parse_args()

    print("=" * 70)
    print("G4X Hierarchical Cell Type Annotation")
    print("=" * 70)

    sample_files = sorted(INPUT_DIR.glob("*_wnn.h5ad"))
    logger.info(f"Found {len(sample_files)} WNN-integrated samples")

    if args.sample:
        sample_files = [f for f in sample_files if args.sample in f.stem]

    results = []
    for f in tqdm(sample_files, desc="Processing"):
        result = process_sample(f)
        if result:
            results.append(result)

    # Summary
    if results:
        df = pd.DataFrame(results)
        df.to_csv(OUTPUT_DIR / "annotation_summary.csv", index=False)
        print("\n" + "=" * 70)
        print("Summary")
        print("=" * 70)
        print(df.to_string(index=False))

    print("\nDone!")


if __name__ == "__main__":
    main()
