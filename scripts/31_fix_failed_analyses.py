#!/usr/bin/env python3
"""
Fix Failed Analyses - Pathway Activity & Niche Identification
=============================================================

Fixes:
1. fig16 - Pathway Activity: Convert pandas nullable Float32 to numpy float32
2. fig19 - Niche Identification: Store and plot niche results properly

Author: Max Van Belkum
Date: 2026-01-20
"""

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
import decoupler as dc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
ADATA_DIR = OUTPUT_DIR / "adata" / "polymathic"
FIG_DIR = OUTPUT_DIR / "figures" / "advanced"

plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

RESPONSE_COLORS = {'R': '#2ecc71', 'NR': '#e74c3c'}
SAMPLE_ORDER = ['YP12A', 'YP12C', 'YP15A', 'YP15C', 'YP03A', 'YP03C', 'YP04C']

SAMPLES = {
    "YP03A": {"response": "NR", "timepoint": "Pre"},
    "YP03C": {"response": "NR", "timepoint": "Post"},
    "YP04C": {"response": "NR", "timepoint": "Post"},
    "YP12A": {"response": "R", "timepoint": "Pre"},
    "YP12C": {"response": "R", "timepoint": "Post"},
    "YP15A": {"response": "R", "timepoint": "Pre"},
    "YP15C": {"response": "R", "timepoint": "Post"},
}


def load_all_samples():
    """Load all polymathic h5ad files."""
    adatas = {}
    for sample in SAMPLES.keys():
        path = ADATA_DIR / f"{sample}_polymathic.h5ad"
        if path.exists():
            adata = sc.read_h5ad(path)
            adata.obs['sample'] = sample
            adata.obs['response'] = SAMPLES[sample]['response']
            adata.obs['timepoint'] = SAMPLES[sample]['timepoint']
            adatas[sample] = adata
            print(f"  Loaded {sample}: {adata.n_obs} spots")
    return adatas


def fix_pathway_analysis(adatas):
    """Run pathway analysis with fixed dtype conversion."""
    print("\n" + "="*60)
    print("FIX 1: Pathway Activity (decoupler)")
    print("="*60)

    # Get models and FIX: convert nullable dtypes
    print("\n  Loading pathway models...")
    progeny = dc.get_progeny(organism='human', top=500)
    progeny['weight'] = progeny['weight'].astype('float32')
    progeny['p_value'] = progeny['p_value'].astype('float32')
    print(f"    PROGENy: {len(progeny)} gene-pathway associations")

    dorothea = dc.get_dorothea(organism='human', levels=['A', 'B', 'C'])
    dorothea['weight'] = dorothea['weight'].astype('float64').astype('float32')
    print(f"    DoRothEA: {len(dorothea)} gene-TF associations")

    pathway_results = {}

    for sample, adata in adatas.items():
        print(f"\n  Processing {sample}...")

        # Run PROGENy
        try:
            dc.run_mlm(
                mat=adata,
                net=progeny,
                source='source',
                target='target',
                weight='weight',
                verbose=False,
                use_raw=False
            )
            adata.obsm['pathway_activities'] = adata.obsm.pop('mlm_estimate')
            adata.obsm['pathway_pvals'] = adata.obsm.pop('mlm_pvals')
            print(f"    PROGENy: {adata.obsm['pathway_activities'].shape[1]} pathways")
            pathway_results[sample] = {
                'progeny': adata.obsm['pathway_activities'].mean(axis=0).to_dict()
            }
        except Exception as e:
            print(f"    PROGENy failed: {e}")

        # Run DoRothEA
        try:
            dc.run_mlm(
                mat=adata,
                net=dorothea,
                source='source',
                target='target',
                weight='weight',
                verbose=False,
                use_raw=False
            )
            adata.obsm['tf_activities'] = adata.obsm.pop('mlm_estimate')
            adata.obsm['tf_pvals'] = adata.obsm.pop('mlm_pvals')
            print(f"    DoRothEA: {adata.obsm['tf_activities'].shape[1]} TFs")
        except Exception as e:
            print(f"    DoRothEA failed: {e}")

    # Generate figure
    print("\n  Generating pathway figure...")

    # Collect pathway activities
    pathway_data = []
    for sample, adata in adatas.items():
        if 'pathway_activities' in adata.obsm:
            mean_activities = adata.obsm['pathway_activities'].mean(axis=0)
            for pathway, activity in mean_activities.items():
                pathway_data.append({
                    'sample': sample,
                    'response': SAMPLES[sample]['response'],
                    'timepoint': SAMPLES[sample]['timepoint'],
                    'pathway': pathway,
                    'activity': activity
                })

    if pathway_data:
        df = pd.DataFrame(pathway_data)

        fig, axes = plt.subplots(2, 2, figsize=(14, 12))

        # A) Heatmap of pathway activities
        ax = axes[0, 0]
        pivot = df.pivot_table(index='pathway', columns='sample', values='activity')
        pivot = pivot[[s for s in SAMPLE_ORDER if s in pivot.columns]]
        sns.heatmap(pivot, cmap='RdBu_r', center=0, ax=ax, cbar_kws={'label': 'Activity'})
        ax.set_title('A) Pathway Activities by Sample')
        ax.set_xlabel('')

        # Add response annotation
        colors = [RESPONSE_COLORS[SAMPLES[s]['response']] for s in pivot.columns]
        for i, (col, color) in enumerate(zip(pivot.columns, colors)):
            ax.add_patch(plt.Rectangle((i, -0.5), 1, 0.5, color=color, clip_on=False))

        # B) R vs NR comparison
        ax = axes[0, 1]
        r_means = df[df['response'] == 'R'].groupby('pathway')['activity'].mean()
        nr_means = df[df['response'] == 'NR'].groupby('pathway')['activity'].mean()

        pathways = r_means.index.tolist()
        x = np.arange(len(pathways))
        width = 0.35

        ax.barh(x - width/2, r_means.values, width, label='Responders', color=RESPONSE_COLORS['R'])
        ax.barh(x + width/2, nr_means.values, width, label='Non-Responders', color=RESPONSE_COLORS['NR'])
        ax.set_yticks(x)
        ax.set_yticklabels(pathways)
        ax.set_xlabel('Mean Activity')
        ax.set_title('B) R vs NR Pathway Activities')
        ax.legend()
        ax.axvline(0, color='black', linewidth=0.5)

        # C) Pre vs Post (Responders)
        ax = axes[1, 0]
        r_pre = df[(df['response'] == 'R') & (df['timepoint'] == 'Pre')].groupby('pathway')['activity'].mean()
        r_post = df[(df['response'] == 'R') & (df['timepoint'] == 'Post')].groupby('pathway')['activity'].mean()

        diff = r_post - r_pre
        diff = diff.sort_values()

        colors = ['#e74c3c' if v < 0 else '#2ecc71' for v in diff.values]
        ax.barh(range(len(diff)), diff.values, color=colors)
        ax.set_yticks(range(len(diff)))
        ax.set_yticklabels(diff.index)
        ax.set_xlabel('Activity Change (Post - Pre)')
        ax.set_title('C) Responders: Pre→Post Change')
        ax.axvline(0, color='black', linewidth=0.5)

        # D) Top differentially active pathways
        ax = axes[1, 1]

        # Calculate effect sizes
        from scipy.stats import ttest_ind
        effects = []
        for pathway in pathways:
            r_vals = df[(df['pathway'] == pathway) & (df['response'] == 'R')]['activity'].values
            nr_vals = df[(df['pathway'] == pathway) & (df['response'] == 'NR')]['activity'].values

            if len(r_vals) > 1 and len(nr_vals) > 1:
                pooled_std = np.sqrt((np.var(r_vals) + np.var(nr_vals)) / 2)
                if pooled_std > 0:
                    d = (np.mean(r_vals) - np.mean(nr_vals)) / pooled_std
                else:
                    d = 0
                _, p = ttest_ind(r_vals, nr_vals, equal_var=False)
                effects.append({'pathway': pathway, 'effect_size': d, 'p_value': p})

        eff_df = pd.DataFrame(effects).sort_values('effect_size')
        colors = ['#e74c3c' if v < 0 else '#2ecc71' for v in eff_df['effect_size'].values]

        ax.barh(range(len(eff_df)), eff_df['effect_size'].values, color=colors)
        ax.set_yticks(range(len(eff_df)))
        ax.set_yticklabels(eff_df['pathway'].values)
        ax.set_xlabel("Cohen's d (R - NR)")
        ax.set_title("D) Effect Sizes (R vs NR)")
        ax.axvline(0, color='black', linewidth=0.5)

        # Add significance markers
        for i, (_, row) in enumerate(eff_df.iterrows()):
            if row['p_value'] < 0.05:
                ax.text(row['effect_size'] + 0.1 * np.sign(row['effect_size']), i, '*',
                       ha='center', va='center', fontsize=14)

        plt.suptitle('Figure 16: Pathway Activity Analysis (PROGENy)', fontsize=14, y=1.02)
        plt.tight_layout()
        plt.savefig(FIG_DIR / 'fig16_pathway_activity.png', bbox_inches='tight', dpi=300)
        plt.close()
        print("    Saved fig16_pathway_activity.png")
    else:
        print("    No pathway data available")

    return pathway_results


def fix_niche_identification(adatas):
    """Run niche identification with proper result storage."""
    print("\n" + "="*60)
    print("FIX 2: Niche Identification")
    print("="*60)

    niche_results = {}

    for sample, adata in adatas.items():
        print(f"\n  Processing {sample}...")

        # Check for cell type column
        ct_col = None
        for col in ['cell_type', 'cell_type_broad', 'leiden']:
            if col in adata.obs.columns:
                ct_col = col
                break

        if ct_col is None:
            print(f"    No cell type column found")
            continue

        # Build spatial neighbors if needed
        if 'spatial_neighbors' not in adata.uns:
            sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

        # Run neighborhood enrichment
        try:
            sq.gr.nhood_enrichment(adata, cluster_key=ct_col)

            # Store results
            zscore = adata.uns[f'{ct_col}_nhood_enrichment']['zscore']
            cell_types = adata.obs[ct_col].cat.categories.tolist()

            niche_results[sample] = {
                'zscore': zscore,
                'cell_types': cell_types,
                'response': SAMPLES[sample]['response']
            }
            print(f"    Computed enrichment for {len(cell_types)} cell types")

        except Exception as e:
            print(f"    Failed: {e}")

    # Generate figure
    print("\n  Generating niche figure...")

    if niche_results:
        n_samples = len(niche_results)
        fig, axes = plt.subplots(2, 4, figsize=(16, 10))
        axes = axes.flatten()

        # Plot individual sample heatmaps
        for idx, (sample, result) in enumerate(niche_results.items()):
            if idx >= 7:
                break
            ax = axes[idx]

            zscore = result['zscore']
            cell_types = result['cell_types']
            response = result['response']

            # Handle inf values
            zscore_clean = np.clip(zscore, -10, 10)

            im = ax.imshow(zscore_clean, cmap='RdBu_r', vmin=-5, vmax=5)
            ax.set_xticks(range(len(cell_types)))
            ax.set_yticks(range(len(cell_types)))
            ax.set_xticklabels(cell_types, rotation=45, ha='right', fontsize=6)
            ax.set_yticklabels(cell_types, fontsize=6)

            title_color = RESPONSE_COLORS[response]
            ax.set_title(f'{sample} ({response})', color=title_color, fontweight='bold')

        # Summary panel - average by response
        ax = axes[7]

        # Get common cell types
        all_ct = set()
        for r in niche_results.values():
            all_ct.update(r['cell_types'])
        common_ct = list(all_ct)[:10]  # Top 10

        # Average z-scores by response
        r_zscores = []
        nr_zscores = []

        for sample, result in niche_results.items():
            if result['response'] == 'R':
                r_zscores.append(np.nanmean(result['zscore']))
            else:
                nr_zscores.append(np.nanmean(result['zscore']))

        # Bar plot of mean enrichment
        x = ['Responders', 'Non-Responders']
        y = [np.mean(r_zscores) if r_zscores else 0,
             np.mean(nr_zscores) if nr_zscores else 0]
        yerr = [np.std(r_zscores) if len(r_zscores) > 1 else 0,
                np.std(nr_zscores) if len(nr_zscores) > 1 else 0]

        bars = ax.bar(x, y, yerr=yerr, capsize=5,
                     color=[RESPONSE_COLORS['R'], RESPONSE_COLORS['NR']])
        ax.set_ylabel('Mean Nhood Enrichment Z-score')
        ax.set_title('H) Average Spatial Organization')
        ax.axhline(0, color='black', linewidth=0.5)

        # Add text box
        textstr = '\n'.join([
            'NICHE IDENTIFICATION',
            '=' * 25,
            'Method: squidpy nhood_enrichment',
            f'Samples analyzed: {len(niche_results)}',
            '',
            'Interpretation:',
            '• Positive z-score: cell types',
            '  co-localize more than expected',
            '• Negative z-score: spatial',
            '  avoidance/segregation',
            '',
            'R vs NR patterns may reveal',
            'tumor microenvironment differences'
        ])
        ax.text(1.3, 0.5, textstr, transform=ax.transAxes, fontsize=8,
               verticalalignment='center', fontfamily='monospace',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.suptitle('Figure 19: Niche Identification (Cell Type Co-localization)',
                    fontsize=14, y=1.02)
        plt.tight_layout()
        plt.savefig(FIG_DIR / 'fig19_niche_identification.png', bbox_inches='tight', dpi=300)
        plt.close()
        print("    Saved fig19_niche_identification.png")
    else:
        print("    No niche results to plot")

    return niche_results


def main():
    print("="*60)
    print("FIXING FAILED ANALYSES")
    print("="*60)

    print("\nLoading samples...")
    adatas = load_all_samples()

    if not adatas:
        print("ERROR: No samples found!")
        return

    # Fix 1: Pathway analysis
    pathway_results = fix_pathway_analysis(adatas)

    # Fix 2: Niche identification
    niche_results = fix_niche_identification(adatas)

    print("\n" + "="*60)
    print("FIXES COMPLETE")
    print("="*60)

    # List generated figures
    print("\nGenerated figures:")
    for f in sorted(FIG_DIR.glob('fig1[69]*.png')):
        print(f"  {f.name}")


if __name__ == '__main__':
    main()
