#!/usr/bin/env python3
"""
G4X Cell Type Annotation - Spatial Hackathon 2026
=================================================
Uses 17 protein markers to annotate cell types in G4X data.

Protein Panel:
- Immune: CD3, CD4, CD8, CD11c, CD20, CD45, CD68, FOXP3, HLA-DR
- Checkpoint: PD1, PDL1
- Proliferation: KI67
- Structural: PanCK, aSMA, CD31, ATPase
- Control: Isotype
"""

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Configuration
RESULTS_DIR = Path("/home/user/spatial-hackathon-2026/results/g4x_choi_batch2")
OUTPUT_DIR = Path("/home/user/spatial-hackathon-2026/results/g4x_choi_batch2/annotated")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR = Path("/home/user/spatial-hackathon-2026/figures/g4x")
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# Cell type definitions based on protein markers
CELL_TYPE_MARKERS = {
    'T cell (CD4+)': {'pos': ['CD3', 'CD4', 'CD45'], 'neg': ['CD8', 'CD20', 'CD68']},
    'T cell (CD8+)': {'pos': ['CD3', 'CD8', 'CD45'], 'neg': ['CD4', 'CD20', 'CD68']},
    'Treg': {'pos': ['CD3', 'CD4', 'FOXP3', 'CD45'], 'neg': ['CD8']},
    'B cell': {'pos': ['CD20', 'CD45'], 'neg': ['CD3', 'CD68']},
    'Macrophage': {'pos': ['CD68', 'CD45'], 'neg': ['CD3', 'CD20']},
    'Dendritic cell': {'pos': ['CD11c', 'HLA-DR', 'CD45'], 'neg': ['CD3', 'CD20', 'CD68']},
    'Epithelial': {'pos': ['PanCK'], 'neg': ['CD45', 'aSMA', 'CD31']},
    'Fibroblast': {'pos': ['aSMA'], 'neg': ['CD45', 'PanCK', 'CD31']},
    'Endothelial': {'pos': ['CD31'], 'neg': ['CD45', 'PanCK', 'aSMA']},
    'Immune (other)': {'pos': ['CD45'], 'neg': ['CD3', 'CD20', 'CD68', 'CD11c']},
}

def normalize_protein_names(adata):
    """Normalize protein names to match our definitions."""
    # Get protein names from uns
    if 'protein_names' not in adata.uns:
        print("  No protein names found in adata.uns")
        return None

    protein_names = adata.uns['protein_names']

    # Create mapping from common variations
    name_mapping = {
        'cd3': 'CD3',
        'cd4': 'CD4',
        'cd8': 'CD8',
        'cd11c': 'CD11c',
        'cd20': 'CD20',
        'cd31': 'CD31',
        'cd45': 'CD45',
        'cd68': 'CD68',
        'foxp3': 'FOXP3',
        'hla-dr': 'HLA-DR',
        'hladr': 'HLA-DR',
        'pd1': 'PD1',
        'pd-1': 'PD1',
        'pdl1': 'PDL1',
        'pd-l1': 'PDL1',
        'ki67': 'KI67',
        'panck': 'PanCK',
        'pan-ck': 'PanCK',
        'asma': 'aSMA',
        'a-sma': 'aSMA',
        'alpha-sma': 'aSMA',
        'atpase': 'ATPase',
        'isotype': 'Isotype',
    }

    normalized = []
    for name in protein_names:
        lower_name = name.lower().replace('_', '-').replace(' ', '')
        normalized.append(name_mapping.get(lower_name, name))

    return normalized

def score_cell_types(adata, protein_data, protein_names):
    """Score cell types based on protein expression."""
    # Z-score normalize protein data
    protein_z = (protein_data - protein_data.mean(axis=0)) / (protein_data.std(axis=0) + 1e-8)

    scores = pd.DataFrame(index=adata.obs_names)

    for cell_type, markers in CELL_TYPE_MARKERS.items():
        pos_markers = [m for m in markers['pos'] if m in protein_names]
        neg_markers = [m for m in markers['neg'] if m in protein_names]

        if not pos_markers:
            continue

        # Positive score
        pos_idx = [protein_names.index(m) for m in pos_markers]
        pos_score = protein_z[:, pos_idx].mean(axis=1)

        # Negative penalty
        neg_score = 0
        if neg_markers:
            neg_idx = [protein_names.index(m) for m in neg_markers]
            neg_score = protein_z[:, neg_idx].mean(axis=1)

        scores[cell_type] = pos_score - 0.5 * neg_score

    return scores

def assign_cell_types(scores, min_score=0.5):
    """Assign cell types based on highest score."""
    if scores.empty:
        return pd.Series(['Unknown'] * len(scores.index), index=scores.index)

    # Get best scoring cell type
    best_type = scores.idxmax(axis=1)
    best_score = scores.max(axis=1)

    # Mark low confidence assignments
    best_type = best_type.where(best_score >= min_score, 'Unknown')

    return best_type

def annotate_sample(sample_path):
    """Annotate a single sample."""
    sample_id = sample_path.stem.replace('_analyzed', '')
    print(f"\n--- Annotating {sample_id} ---")

    # Load data
    adata = sc.read_h5ad(sample_path)
    print(f"  Loaded: {adata.n_obs} cells")

    # Check for protein data
    if 'protein' not in adata.obsm:
        print("  No protein data, skipping")
        return None

    protein_data = adata.obsm['protein']
    protein_names = normalize_protein_names(adata)

    if protein_names is None:
        print("  Could not normalize protein names, skipping")
        return None

    print(f"  Proteins: {protein_names}")

    # Score and assign cell types
    scores = score_cell_types(adata, protein_data, protein_names)
    cell_types = assign_cell_types(scores)

    # Add to adata
    adata.obs['cell_type'] = cell_types
    for col in scores.columns:
        adata.obs[f'score_{col}'] = scores[col].values

    # Summary
    type_counts = cell_types.value_counts()
    print(f"  Cell types:")
    for ct, count in type_counts.items():
        pct = 100 * count / len(cell_types)
        print(f"    {ct}: {count:,} ({pct:.1f}%)")

    # Save annotated file
    out_path = OUTPUT_DIR / f"{sample_id}_annotated.h5ad"
    adata.write(out_path)
    print(f"  Saved: {out_path}")

    return {
        'sample_id': sample_id,
        'n_cells': adata.n_obs,
        'cell_type_counts': type_counts.to_dict()
    }

def plot_cell_type_composition(results, output_path):
    """Plot cell type composition across samples."""
    # Prepare data
    all_types = set()
    for r in results:
        if r and 'cell_type_counts' in r:
            all_types.update(r['cell_type_counts'].keys())

    all_types = sorted(all_types)

    # Create matrix
    samples = [r['sample_id'] for r in results if r]
    data = []
    for r in results:
        if r and 'cell_type_counts' in r:
            row = [r['cell_type_counts'].get(ct, 0) for ct in all_types]
            data.append(row)

    df = pd.DataFrame(data, columns=all_types, index=samples)

    # Normalize to percentages
    df_pct = df.div(df.sum(axis=1), axis=0) * 100

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    df_pct.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
    ax.set_ylabel('Percentage')
    ax.set_xlabel('Sample')
    ax.set_title('Cell Type Composition Across G4X Samples')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")

    return df, df_pct

def main():
    """Main annotation pipeline."""
    print("=" * 60)
    print("G4X Cell Type Annotation")
    print("=" * 60)

    # Get all analyzed samples
    sample_files = sorted(RESULTS_DIR.glob("*_analyzed.h5ad"))
    print(f"Found {len(sample_files)} analyzed samples")

    # Annotate each sample
    results = []
    for sample_path in tqdm(sample_files, desc="Annotating samples"):
        result = annotate_sample(sample_path)
        results.append(result)

    # Generate summary
    valid_results = [r for r in results if r]
    if valid_results:
        # Save summary
        summary = []
        for r in valid_results:
            row = {'sample_id': r['sample_id'], 'n_cells': r['n_cells']}
            row.update(r['cell_type_counts'])
            summary.append(row)

        summary_df = pd.DataFrame(summary)
        summary_df.to_csv(OUTPUT_DIR / "cell_type_summary.csv", index=False)
        print(f"\nSummary saved to {OUTPUT_DIR / 'cell_type_summary.csv'}")

        # Plot composition
        plot_cell_type_composition(valid_results, FIGURE_DIR / "cell_type_composition.png")

    print("\nAnnotation complete!")
    return results

if __name__ == "__main__":
    main()
