#!/usr/bin/env python
"""
05_spatial_statistics.py - Spatial neighborhood analysis for G4X pilot samples

Analyses:
1. Neighborhood enrichment (which cell types co-localize?)
2. Ripley's L function (spatial clustering)
3. Cell type co-occurrence patterns across disease stages

Author: Claude Code
Date: 2026-01-22
"""

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
RESULTS_DIR = Path("results/pilot")
FIGURES_DIR = RESULTS_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

def load_annotated_samples():
    """Load annotated h5ad files."""
    samples = {}
    for sample_id in ['E02', 'F02', 'G02']:
        path = RESULTS_DIR / f"{sample_id}_annotated.h5ad"
        if path.exists():
            samples[sample_id] = sc.read_h5ad(path)
            print(f"Loaded {sample_id}: {samples[sample_id].n_obs} cells")
        else:
            print(f"WARNING: {path} not found")
    return samples

def build_spatial_graph(adata, n_neighbors=15, coord_type='generic'):
    """Build spatial connectivity graph."""
    # Check for spatial coordinates
    if 'spatial' not in adata.obsm:
        if 'X_spatial' in adata.obsm:
            adata.obsm['spatial'] = adata.obsm['X_spatial']
        elif 'x_centroid' in adata.obs and 'y_centroid' in adata.obs:
            adata.obsm['spatial'] = np.column_stack([
                adata.obs['x_centroid'].values,
                adata.obs['y_centroid'].values
            ])
        else:
            print("No spatial coordinates found!")
            return None

    # Build spatial neighbors graph
    sq.gr.spatial_neighbors(
        adata,
        n_neighs=n_neighbors,
        coord_type=coord_type
    )
    print(f"  Built spatial graph: {n_neighbors} neighbors")
    return adata

def neighborhood_enrichment(adata, cluster_key='cell_type_marker'):
    """Calculate neighborhood enrichment."""
    if 'spatial_connectivities' not in adata.obsp:
        print("  No spatial graph - building...")
        adata = build_spatial_graph(adata)
        if adata is None:
            return None

    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
    return adata

def analyze_sample(adata, sample_id, stage):
    """Run full spatial analysis for one sample."""
    print(f"\n{'='*60}")
    print(f"Analyzing {sample_id} ({stage})")
    print(f"{'='*60}")

    # Build spatial graph
    adata = build_spatial_graph(adata)
    if adata is None:
        return None

    # Get cluster key (try different possibilities)
    cluster_key = None
    for key in ['cell_type_marker', 'cell_type', 'leiden', 'cluster']:
        if key in adata.obs.columns:
            cluster_key = key
            break

    if cluster_key is None:
        print("  No cluster annotation found!")
        return None

    print(f"  Using cluster key: {cluster_key}")
    n_clusters = adata.obs[cluster_key].nunique()
    print(f"  Number of cell types: {n_clusters}")

    # Neighborhood enrichment
    print("  Computing neighborhood enrichment...")
    try:
        adata = neighborhood_enrichment(adata, cluster_key=cluster_key)
    except Exception as e:
        print(f"  Enrichment failed: {e}")

    # Ripley's L function (if not too many cells)
    if adata.n_obs < 50000:
        print("  Computing Ripley's L function...")
        try:
            sq.gr.ripley(adata, cluster_key=cluster_key, mode='L')
        except Exception as e:
            print(f"  Ripley's L failed: {e}")
    else:
        print(f"  Skipping Ripley's L (n={adata.n_obs} > 50000)")

    # Co-occurrence
    print("  Computing co-occurrence...")
    try:
        sq.gr.co_occurrence(
            adata,
            cluster_key=cluster_key,
            spatial_key='spatial',
            n_steps=50
        )
    except Exception as e:
        print(f"  Co-occurrence failed: {e}")

    return adata

def plot_neighborhood_enrichment(samples, cluster_key='cell_type_marker'):
    """Plot neighborhood enrichment heatmaps for all samples."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    stages = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}

    for idx, (sample_id, stage) in enumerate(stages.items()):
        adata = samples.get(sample_id)
        if adata is None or f'{cluster_key}_nhood_enrichment' not in adata.uns:
            axes[idx].text(0.5, 0.5, f'{sample_id}\nNo data', ha='center', va='center')
            axes[idx].set_title(f"{sample_id} ({stage})")
            continue

        # Get enrichment matrix
        zscore = adata.uns[f'{cluster_key}_nhood_enrichment']['zscore']
        labels = adata.obs[cluster_key].cat.categories

        # Plot heatmap
        im = axes[idx].imshow(zscore, cmap='RdBu_r', vmin=-5, vmax=5, aspect='auto')
        axes[idx].set_xticks(range(len(labels)))
        axes[idx].set_yticks(range(len(labels)))
        axes[idx].set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        axes[idx].set_yticklabels(labels, fontsize=8)
        axes[idx].set_title(f"{sample_id} ({stage})")

    plt.colorbar(im, ax=axes, label='Z-score', shrink=0.6)
    plt.suptitle("Neighborhood Enrichment (cell type co-localization)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "spatial_neighborhood_enrichment.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {FIGURES_DIR / 'spatial_neighborhood_enrichment.png'}")

def plot_co_occurrence(samples, cluster_key='cell_type_marker'):
    """Plot co-occurrence analysis."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    stages = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}

    for idx, (sample_id, stage) in enumerate(stages.items()):
        adata = samples.get(sample_id)
        if adata is None:
            continue

        try:
            sq.pl.co_occurrence(
                adata,
                cluster_key=cluster_key,
                clusters=adata.obs[cluster_key].cat.categories[:5].tolist(),  # Top 5
                ax=axes[idx]
            )
            axes[idx].set_title(f"{sample_id} ({stage})")
        except Exception as e:
            axes[idx].text(0.5, 0.5, f"Error: {e}", ha='center', va='center')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "spatial_co_occurrence.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {FIGURES_DIR / 'spatial_co_occurrence.png'}")

def summarize_spatial_patterns(samples, cluster_key='cell_type_marker'):
    """Extract key spatial patterns comparing stages."""
    results = []

    stages = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}

    for sample_id, stage in stages.items():
        adata = samples.get(sample_id)
        if adata is None:
            continue

        if f'{cluster_key}_nhood_enrichment' not in adata.uns:
            continue

        zscore = adata.uns[f'{cluster_key}_nhood_enrichment']['zscore']
        labels = list(adata.obs[cluster_key].cat.categories)

        # Find strongest co-localizations (off-diagonal)
        n = len(labels)
        for i in range(n):
            for j in range(i+1, n):
                results.append({
                    'sample': sample_id,
                    'stage': stage,
                    'cell_type_1': labels[i],
                    'cell_type_2': labels[j],
                    'zscore': zscore[i, j]
                })

    df = pd.DataFrame(results)

    # Get top co-localizations per stage
    print("\n" + "="*60)
    print("TOP CELL TYPE CO-LOCALIZATIONS BY STAGE")
    print("="*60)

    for stage in ['Normal', 'Metaplasia', 'Cancer']:
        stage_df = df[df['stage'] == stage].nlargest(5, 'zscore')
        print(f"\n{stage}:")
        for _, row in stage_df.iterrows():
            print(f"  {row['cell_type_1']} <-> {row['cell_type_2']}: z={row['zscore']:.2f}")

    # Save full results
    df.to_csv(RESULTS_DIR / "spatial_colocalization_summary.csv", index=False)
    print(f"\nSaved: {RESULTS_DIR / 'spatial_colocalization_summary.csv'}")

    return df

def main():
    print("="*60)
    print("G4X PILOT - SPATIAL STATISTICS ANALYSIS")
    print("="*60)

    # Load samples
    samples = load_annotated_samples()
    if not samples:
        print("ERROR: No annotated samples found!")
        return

    # Analyze each sample
    stages = {'E02': 'Normal', 'F02': 'Metaplasia', 'G02': 'Cancer'}
    analyzed_samples = {}

    for sample_id, stage in stages.items():
        if sample_id in samples:
            result = analyze_sample(samples[sample_id], sample_id, stage)
            if result is not None:
                analyzed_samples[sample_id] = result
                # Save updated h5ad with spatial stats
                result.write_h5ad(RESULTS_DIR / f"{sample_id}_spatial.h5ad")

    # Determine cluster key
    cluster_key = 'cell_type_marker'
    if analyzed_samples:
        first_sample = list(analyzed_samples.values())[0]
        for key in ['cell_type_marker', 'cell_type', 'leiden']:
            if key in first_sample.obs.columns:
                cluster_key = key
                break

    # Generate plots
    print("\n" + "="*60)
    print("GENERATING SPATIAL FIGURES")
    print("="*60)

    if analyzed_samples:
        try:
            plot_neighborhood_enrichment(analyzed_samples, cluster_key)
        except Exception as e:
            print(f"Enrichment plot failed: {e}")

        try:
            plot_co_occurrence(analyzed_samples, cluster_key)
        except Exception as e:
            print(f"Co-occurrence plot failed: {e}")

        # Summarize patterns
        summarize_spatial_patterns(analyzed_samples, cluster_key)

    print("\n" + "="*60)
    print("SPATIAL ANALYSIS COMPLETE")
    print("="*60)

if __name__ == "__main__":
    main()
