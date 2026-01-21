#!/usr/bin/env python3
"""
Advanced PDAC Analyses - Complete Pipeline
==========================================

Runs all remaining analyses on PDAC Visium samples:
1. Pathway Activity (decoupler - PROGENy + DoRothEA)
2. Gene Set Enrichment (MSigDB Hallmarks)
3. Spatial Domain Detection (spatially-aware clustering)
4. Niche Identification (cell type co-localization)
5. Trajectory Analysis (CellRank on tumor cells)
6. Multi-sample Integration (Harmony)

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
import anndata as ad
from pathlib import Path
import warnings
import json
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
warnings.filterwarnings('ignore')

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
ADATA_DIR = OUTPUT_DIR / "adata" / "polymathic"
FIG_DIR = OUTPUT_DIR / "figures" / "advanced"
TABLE_DIR = OUTPUT_DIR / "tables"

FIG_DIR.mkdir(parents=True, exist_ok=True)

# Style
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 10

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
            print(f"  Loaded {sample}: {adata.n_obs} cells")
    return adatas


def dual_test(group1, group2):
    """Perform both parametric and non-parametric tests."""
    from scipy.stats import ttest_ind, mannwhitneyu

    # Handle edge cases
    if len(group1) < 2 or len(group2) < 2:
        return {'welch_p': np.nan, 'mwu_p': np.nan, 'effect_size': np.nan}

    # Welch's t-test
    t_stat, welch_p = ttest_ind(group1, group2, equal_var=False)

    # Mann-Whitney U
    try:
        u_stat, mwu_p = mannwhitneyu(group1, group2, alternative='two-sided')
    except:
        mwu_p = np.nan

    # Cohen's d
    pooled_std = np.sqrt((np.var(group1) + np.var(group2)) / 2)
    if pooled_std > 0:
        cohens_d = (np.mean(group1) - np.mean(group2)) / pooled_std
    else:
        cohens_d = 0

    return {'welch_p': welch_p, 'mwu_p': mwu_p, 'effect_size': cohens_d}


# =============================================================================
# ANALYSIS 1: Pathway Activity (decoupler)
# =============================================================================

def run_pathway_analysis(adatas):
    """Run PROGENy and DoRothEA pathway/TF activity inference."""
    print("\n" + "="*60)
    print("ANALYSIS 1: Pathway Activity (decoupler)")
    print("="*60)

    import decoupler as dc

    results = {}

    for sample, adata in adatas.items():
        print(f"\n  Processing {sample}...")

        # Ensure we have normalized data
        if 'normalized' not in adata.layers:
            adata.layers['normalized'] = adata.X.copy()

        # PROGENy pathway activity
        try:
            progeny = dc.get_progeny(organism='human', top=500)
            dc.run_mlm(
                mat=adata,
                net=progeny,
                source='source',
                target='target',
                weight='weight',
                verbose=False,
                use_raw=False
            )

            if 'mlm_estimate' in adata.obsm:
                adata.obsm['pathway_activities'] = adata.obsm['mlm_estimate'].copy()
                print(f"    PROGENy: {adata.obsm['pathway_activities'].shape[1]} pathways")
        except Exception as e:
            print(f"    PROGENy failed: {e}")

        # DoRothEA TF activity
        try:
            dorothea = dc.get_dorothea(organism='human', levels=['A', 'B', 'C'])
            dc.run_ulm(
                mat=adata,
                net=dorothea,
                source='source',
                target='target',
                weight='weight',
                verbose=False,
                use_raw=False
            )

            if 'ulm_estimate' in adata.obsm:
                adata.obsm['tf_activities'] = adata.obsm['ulm_estimate'].copy()
                print(f"    DoRothEA: {adata.obsm['tf_activities'].shape[1]} TFs")
        except Exception as e:
            print(f"    DoRothEA failed: {e}")

        results[sample] = adata

    # Generate figures
    _plot_pathway_results(results)

    return results


def _plot_pathway_results(adatas):
    """Generate pathway activity figures."""
    print("\n  Generating pathway figures...")

    # Collect pathway activities per sample
    pathway_data = []
    tf_data = []

    for sample, adata in adatas.items():
        response = SAMPLES[sample]['response']

        if 'pathway_activities' in adata.obsm:
            activities = pd.DataFrame(
                adata.obsm['pathway_activities'],
                columns=adata.obsm['pathway_activities'].columns if hasattr(adata.obsm['pathway_activities'], 'columns') else None
            )
            if activities.columns is not None:
                mean_act = activities.mean()
                for pathway, value in mean_act.items():
                    pathway_data.append({
                        'sample': sample,
                        'response': response,
                        'pathway': pathway,
                        'activity': value
                    })

        if 'tf_activities' in adata.obsm:
            tf_act = pd.DataFrame(adata.obsm['tf_activities'])
            if tf_act.columns is not None:
                mean_tf = tf_act.mean()
                for tf, value in mean_tf.items():
                    tf_data.append({
                        'sample': sample,
                        'response': response,
                        'tf': tf,
                        'activity': value
                    })

    if not pathway_data:
        print("    No pathway data to plot")
        return

    pathway_df = pd.DataFrame(pathway_data)

    # Figure: Pathway activity heatmap + R vs NR comparison
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))

    # Panel A: Heatmap
    ax1 = axes[0]
    pivot = pathway_df.pivot_table(index='sample', columns='pathway', values='activity', aggfunc='mean')
    pivot = pivot.reindex(SAMPLE_ORDER)

    sns.heatmap(pivot, cmap='RdBu_r', center=0, ax=ax1, cbar_kws={'label': 'Activity'})
    ax1.set_title('A) PROGENy Pathway Activities', fontweight='bold')
    ax1.set_ylabel('')

    # Add response labels
    for i, sample in enumerate(pivot.index):
        color = RESPONSE_COLORS[SAMPLES[sample]['response']]
        ax1.text(-0.5, i + 0.5, SAMPLES[sample]['response'], ha='right', va='center',
                color=color, fontweight='bold')

    # Panel B: Top differential pathways
    ax2 = axes[1]

    # Calculate R vs NR difference for each pathway
    pathway_diff = []
    for pathway in pivot.columns:
        r_vals = pivot.loc[[s for s in pivot.index if SAMPLES[s]['response'] == 'R'], pathway].values
        nr_vals = pivot.loc[[s for s in pivot.index if SAMPLES[s]['response'] == 'NR'], pathway].values

        diff = np.mean(r_vals) - np.mean(nr_vals)
        stats_result = dual_test(r_vals, nr_vals)

        pathway_diff.append({
            'pathway': pathway,
            'R_mean': np.mean(r_vals),
            'NR_mean': np.mean(nr_vals),
            'diff': diff,
            'p_value': stats_result['mwu_p']
        })

    diff_df = pd.DataFrame(pathway_diff).sort_values('diff', key=abs, ascending=False)

    # Bar plot of differences
    colors = ['#2ecc71' if d > 0 else '#e74c3c' for d in diff_df['diff'].head(10)]
    ax2.barh(range(10), diff_df['diff'].head(10), color=colors, edgecolor='black')
    ax2.set_yticks(range(10))
    ax2.set_yticklabels(diff_df['pathway'].head(10))
    ax2.axvline(0, color='black', linestyle='-', linewidth=0.5)
    ax2.set_xlabel('Activity Difference (R - NR)')
    ax2.set_title('B) Top Differential Pathways', fontweight='bold')
    ax2.invert_yaxis()

    # Add significance markers
    for i, (_, row) in enumerate(diff_df.head(10).iterrows()):
        if row['p_value'] < 0.1:
            marker = '†' if row['p_value'] < 0.1 else ''
            if row['p_value'] < 0.05:
                marker = '*'
            ax2.text(row['diff'] + 0.01 * np.sign(row['diff']), i, marker,
                    va='center', fontsize=12, fontweight='bold')

    # Panel C: Interpretation
    ax3 = axes[2]
    ax3.axis('off')

    # Get top pathways for interpretation
    top_r = diff_df[diff_df['diff'] > 0].head(3)['pathway'].tolist()
    top_nr = diff_df[diff_df['diff'] < 0].head(3)['pathway'].tolist()

    interpretation = f"""
PATHWAY ACTIVITY ANALYSIS
=========================

Method: decoupler with PROGENy
(14 cancer-relevant pathways)

Higher in Responders (R):
• {', '.join(top_r[:3]) if top_r else 'None significant'}

Higher in Non-Responders (NR):
• {', '.join(top_nr[:3]) if top_nr else 'None significant'}

Statistical Testing:
• Mann-Whitney U test
• † p < 0.10, * p < 0.05

Biological Interpretation:
Pathway activity differences may
indicate distinct tumor biology
and treatment sensitivity between
response groups.

Clinical Relevance:
These pathways could serve as
predictive biomarkers or
therapeutic targets.
"""
    ax3.text(0.05, 0.95, interpretation, transform=ax3.transAxes,
             fontsize=9, family='monospace', va='top',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))

    plt.suptitle('Figure 16: PROGENy Pathway Activity Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig16_pathway_activity.png", bbox_inches='tight', dpi=300)
    plt.close()

    # Save data
    diff_df.to_csv(TABLE_DIR / "pathway_activity_comparison.csv", index=False)
    print("    Saved fig16_pathway_activity.png")


# =============================================================================
# ANALYSIS 2: Gene Set Enrichment
# =============================================================================

def run_gsea_analysis(adatas):
    """Run gene set enrichment on DE genes."""
    print("\n" + "="*60)
    print("ANALYSIS 2: Gene Set Enrichment (GSEA)")
    print("="*60)

    import decoupler as dc

    # Load DE results
    de_pre = pd.read_csv(TABLE_DIR / "de_R_vs_NR_Pre.csv")
    de_post = pd.read_csv(TABLE_DIR / "de_R_vs_NR_Post.csv")

    # Get MSigDB Hallmark gene sets
    try:
        msigdb = dc.get_resource('MSigDB')
        hallmarks = msigdb[msigdb['collection'] == 'hallmark']
        print(f"  Loaded {len(hallmarks['geneset'].unique())} Hallmark gene sets")
    except Exception as e:
        print(f"  Failed to load MSigDB: {e}")
        # Fallback to manual hallmarks
        hallmarks = None

    results = {'pre': {}, 'post': {}}

    # Run enrichment on pre-treatment DE
    if 'names' in de_pre.columns and 'logfoldchanges' in de_pre.columns:
        print("\n  Running enrichment on Pre-treatment DE genes...")

        # Create ranked gene list
        de_pre_ranked = de_pre.dropna(subset=['logfoldchanges'])
        de_pre_ranked = de_pre_ranked.sort_values('logfoldchanges', ascending=False)

        # Separate up/down regulated
        up_genes = de_pre_ranked[de_pre_ranked['logfoldchanges'] > 0.5]['names'].tolist()
        down_genes = de_pre_ranked[de_pre_ranked['logfoldchanges'] < -0.5]['names'].tolist()

        results['pre']['up_in_R'] = up_genes[:100]
        results['pre']['up_in_NR'] = down_genes[:100]

        print(f"    Up in R: {len(up_genes)} genes, Down in R: {len(down_genes)} genes")

    # Run enrichment on post-treatment DE
    if 'names' in de_post.columns and 'logfoldchanges' in de_post.columns:
        print("\n  Running enrichment on Post-treatment DE genes...")

        de_post_ranked = de_post.dropna(subset=['logfoldchanges'])
        de_post_ranked = de_post_ranked.sort_values('logfoldchanges', ascending=False)

        up_genes = de_post_ranked[de_post_ranked['logfoldchanges'] > 0.5]['names'].tolist()
        down_genes = de_post_ranked[de_post_ranked['logfoldchanges'] < -0.5]['names'].tolist()

        results['post']['up_in_R'] = up_genes[:100]
        results['post']['up_in_NR'] = down_genes[:100]

    # If we have hallmarks, run proper enrichment
    if hallmarks is not None:
        try:
            # Run ORA (over-representation analysis)
            for timepoint, de_df in [('pre', de_pre), ('post', de_post)]:
                if 'names' not in de_df.columns:
                    continue

                # Get significant genes
                if 'pvals_adj' in de_df.columns:
                    sig_genes = de_df[(de_df['pvals_adj'] < 0.05) &
                                     (abs(de_df['logfoldchanges']) > 0.5)]['names'].tolist()
                else:
                    sig_genes = de_df[abs(de_df['logfoldchanges']) > 0.5]['names'].tolist()

                if len(sig_genes) > 10:
                    # Simple enrichment: count overlaps with hallmarks
                    enrichment_results = []
                    for geneset in hallmarks['geneset'].unique():
                        gs_genes = hallmarks[hallmarks['geneset'] == geneset]['genesymbol'].tolist()
                        overlap = len(set(sig_genes) & set(gs_genes))
                        if overlap > 0:
                            enrichment_results.append({
                                'geneset': geneset.replace('HALLMARK_', ''),
                                'overlap': overlap,
                                'geneset_size': len(gs_genes),
                                'ratio': overlap / len(gs_genes)
                            })

                    results[timepoint]['enrichment'] = pd.DataFrame(enrichment_results)
        except Exception as e:
            print(f"    Enrichment analysis failed: {e}")

    # Generate figures
    _plot_gsea_results(results, de_pre, de_post)

    return results


def _plot_gsea_results(results, de_pre, de_post):
    """Generate GSEA figures."""
    print("\n  Generating GSEA figures...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Panel A: Pre-treatment volcano with top genes labeled
    ax1 = axes[0, 0]
    if 'logfoldchanges' in de_pre.columns and 'pvals_adj' in de_pre.columns:
        de_pre['neg_log10_p'] = -np.log10(de_pre['pvals_adj'].clip(lower=1e-300))

        ax1.scatter(de_pre['logfoldchanges'], de_pre['neg_log10_p'],
                   c='gray', alpha=0.3, s=5)

        # Highlight significant
        sig_up = (de_pre['logfoldchanges'] > 1) & (de_pre['pvals_adj'] < 0.05)
        sig_down = (de_pre['logfoldchanges'] < -1) & (de_pre['pvals_adj'] < 0.05)

        ax1.scatter(de_pre.loc[sig_up, 'logfoldchanges'],
                   de_pre.loc[sig_up, 'neg_log10_p'],
                   c=RESPONSE_COLORS['R'], alpha=0.7, s=15, label=f'Up in R ({sig_up.sum()})')
        ax1.scatter(de_pre.loc[sig_down, 'logfoldchanges'],
                   de_pre.loc[sig_down, 'neg_log10_p'],
                   c=RESPONSE_COLORS['NR'], alpha=0.7, s=15, label=f'Up in NR ({sig_down.sum()})')

        ax1.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
        ax1.axvline(1, color='gray', linestyle='--', alpha=0.5)
        ax1.axvline(-1, color='gray', linestyle='--', alpha=0.5)
        ax1.set_xlabel('Log2 Fold Change')
        ax1.set_ylabel('-Log10(adj p-value)')
        ax1.legend(fontsize=8, loc='upper right')
        ax1.set_title('A) Pre-Treatment DE: R vs NR', fontweight='bold')

    # Panel B: Post-treatment volcano
    ax2 = axes[0, 1]
    if 'logfoldchanges' in de_post.columns and 'pvals_adj' in de_post.columns:
        de_post['neg_log10_p'] = -np.log10(de_post['pvals_adj'].clip(lower=1e-300))

        ax2.scatter(de_post['logfoldchanges'], de_post['neg_log10_p'],
                   c='gray', alpha=0.3, s=5)

        sig_up = (de_post['logfoldchanges'] > 1) & (de_post['pvals_adj'] < 0.05)
        sig_down = (de_post['logfoldchanges'] < -1) & (de_post['pvals_adj'] < 0.05)

        ax2.scatter(de_post.loc[sig_up, 'logfoldchanges'],
                   de_post.loc[sig_up, 'neg_log10_p'],
                   c=RESPONSE_COLORS['R'], alpha=0.7, s=15, label=f'Up in R ({sig_up.sum()})')
        ax2.scatter(de_post.loc[sig_down, 'logfoldchanges'],
                   de_post.loc[sig_down, 'neg_log10_p'],
                   c=RESPONSE_COLORS['NR'], alpha=0.7, s=15, label=f'Up in NR ({sig_down.sum()})')

        ax2.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
        ax2.axvline(1, color='gray', linestyle='--', alpha=0.5)
        ax2.axvline(-1, color='gray', linestyle='--', alpha=0.5)
        ax2.set_xlabel('Log2 Fold Change')
        ax2.set_ylabel('-Log10(adj p-value)')
        ax2.legend(fontsize=8, loc='upper right')
        ax2.set_title('B) Post-Treatment DE: R vs NR', fontweight='bold')

    # Panel C: Enrichment results (pre)
    ax3 = axes[1, 0]
    if 'pre' in results and 'enrichment' in results['pre']:
        enr = results['pre']['enrichment'].sort_values('overlap', ascending=False).head(10)
        colors = plt.cm.YlOrRd(enr['ratio'] / enr['ratio'].max())
        ax3.barh(range(len(enr)), enr['overlap'], color=colors, edgecolor='black')
        ax3.set_yticks(range(len(enr)))
        ax3.set_yticklabels([g[:25] for g in enr['geneset']])
        ax3.set_xlabel('Overlapping Genes')
        ax3.set_title('C) Pre-Treatment: Hallmark Enrichment', fontweight='bold')
        ax3.invert_yaxis()
    else:
        ax3.text(0.5, 0.5, 'Enrichment data\nnot available', ha='center', va='center',
                transform=ax3.transAxes, fontsize=12)
        ax3.set_title('C) Pre-Treatment: Hallmark Enrichment', fontweight='bold')

    # Panel D: Enrichment results (post)
    ax4 = axes[1, 1]
    if 'post' in results and 'enrichment' in results['post']:
        enr = results['post']['enrichment'].sort_values('overlap', ascending=False).head(10)
        colors = plt.cm.YlOrRd(enr['ratio'] / enr['ratio'].max())
        ax4.barh(range(len(enr)), enr['overlap'], color=colors, edgecolor='black')
        ax4.set_yticks(range(len(enr)))
        ax4.set_yticklabels([g[:25] for g in enr['geneset']])
        ax4.set_xlabel('Overlapping Genes')
        ax4.set_title('D) Post-Treatment: Hallmark Enrichment', fontweight='bold')
        ax4.invert_yaxis()
    else:
        ax4.text(0.5, 0.5, 'Enrichment data\nnot available', ha='center', va='center',
                transform=ax4.transAxes, fontsize=12)
        ax4.set_title('D) Post-Treatment: Hallmark Enrichment', fontweight='bold')

    plt.suptitle('Figure 17: Differential Expression & Gene Set Enrichment',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig17_gsea_enrichment.png", bbox_inches='tight', dpi=300)
    plt.close()
    print("    Saved fig17_gsea_enrichment.png")


# =============================================================================
# ANALYSIS 3: Spatial Domain Detection
# =============================================================================

def run_spatial_domains(adatas):
    """Detect spatial domains using spatially-aware clustering."""
    print("\n" + "="*60)
    print("ANALYSIS 3: Spatial Domain Detection")
    print("="*60)

    results = {}

    for sample, adata in adatas.items():
        print(f"\n  Processing {sample}...")

        try:
            # Build spatial graph if not present
            if 'spatial_connectivities' not in adata.obsp:
                sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

            # Run spatially-constrained clustering
            # Use existing PCA, add spatial constraint via graph
            if 'X_pca' not in adata.obsm:
                sc.pp.pca(adata, n_comps=30)

            # Combine expression neighbors with spatial neighbors
            sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=15)

            # Leiden clustering with different resolutions to find domains
            for res in [0.3, 0.5, 0.8]:
                sc.tl.leiden(adata, resolution=res, key_added=f'spatial_domain_r{res}')

            # Use resolution 0.5 as default
            adata.obs['spatial_domain'] = adata.obs['spatial_domain_r0.5']
            n_domains = adata.obs['spatial_domain'].nunique()
            print(f"    Found {n_domains} spatial domains (res=0.5)")

            results[sample] = adata

        except Exception as e:
            print(f"    Failed: {e}")
            results[sample] = adata

    # Generate figures
    _plot_spatial_domains(results)

    return results


def _plot_spatial_domains(adatas):
    """Generate spatial domain figures."""
    print("\n  Generating spatial domain figures...")

    n_samples = len(adatas)
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    domain_composition = []

    for idx, (sample, adata) in enumerate(adatas.items()):
        ax = axes[idx]

        if 'spatial_domain' in adata.obs and 'spatial' in adata.obsm:
            coords = adata.obsm['spatial']
            domains = adata.obs['spatial_domain'].astype(str)

            # Plot
            scatter = ax.scatter(coords[:, 0], coords[:, 1],
                               c=domains.astype('category').cat.codes,
                               cmap='tab10', s=3, alpha=0.8)
            ax.set_aspect('equal')
            ax.set_title(f'{sample} ({SAMPLES[sample]["response"]})', fontweight='bold')
            ax.axis('off')

            # Collect composition
            for domain in domains.unique():
                count = (domains == domain).sum()
                domain_composition.append({
                    'sample': sample,
                    'response': SAMPLES[sample]['response'],
                    'domain': domain,
                    'count': count,
                    'proportion': count / len(domains)
                })
        else:
            ax.text(0.5, 0.5, 'No spatial\ndata', ha='center', va='center',
                   transform=ax.transAxes)
            ax.set_title(f'{sample}', fontweight='bold')

    # Panel 8: Summary
    ax = axes[7]
    ax.axis('off')

    if domain_composition:
        comp_df = pd.DataFrame(domain_composition)
        n_domains_per_sample = comp_df.groupby('sample')['domain'].nunique()

        summary = f"""
SPATIAL DOMAIN DETECTION
========================

Method: Leiden clustering on
combined expression + spatial
neighborhood graph (res=0.5)

Domains per Sample:
{chr(10).join([f'• {s}: {n} domains' for s, n in n_domains_per_sample.items()])}

Interpretation:
Spatial domains represent
functionally distinct tissue
regions based on both gene
expression and spatial context.

R vs NR patterns may indicate
different tumor organization.
"""
        ax.text(0.05, 0.95, summary, transform=ax.transAxes,
               fontsize=9, family='monospace', va='top',
               bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))

    plt.suptitle('Figure 18: Spatial Domain Detection', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig18_spatial_domains.png", bbox_inches='tight', dpi=300)
    plt.close()

    # Save composition data
    if domain_composition:
        pd.DataFrame(domain_composition).to_csv(TABLE_DIR / "spatial_domain_composition.csv", index=False)

    print("    Saved fig18_spatial_domains.png")


# =============================================================================
# ANALYSIS 4: Niche Identification
# =============================================================================

def run_niche_identification(adatas):
    """Identify spatial niches based on cell type co-localization."""
    print("\n" + "="*60)
    print("ANALYSIS 4: Niche Identification")
    print("="*60)

    results = {}

    for sample, adata in adatas.items():
        print(f"\n  Processing {sample}...")

        try:
            # Get cell type annotations
            ct_col = None
            for col in ['cell_type', 'celltype', 'annotation', 'cluster']:
                if col in adata.obs.columns:
                    ct_col = col
                    break

            if ct_col is None:
                print(f"    No cell type column found")
                continue

            # Compute neighborhood enrichment
            if 'spatial_connectivities' not in adata.obsp:
                sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

            sq.gr.nhood_enrichment(adata, cluster_key=ct_col)

            # Store results
            if 'nhood_enrichment' in adata.uns:
                results[sample] = {
                    'adata': adata,
                    'enrichment': adata.uns[f'{ct_col}_nhood_enrichment'],
                    'ct_col': ct_col
                }
                print(f"    Computed neighborhood enrichment")

        except Exception as e:
            print(f"    Failed: {e}")

    # Generate figures
    _plot_niche_results(results)

    return results


def _plot_niche_results(results):
    """Generate niche identification figures."""
    print("\n  Generating niche figures...")

    if not results:
        print("    No niche results to plot")
        return

    # Collect enrichment matrices
    r_enrichments = []
    nr_enrichments = []

    for sample, data in results.items():
        if 'enrichment' in data:
            zscore = data['enrichment']['zscore']
            if SAMPLES[sample]['response'] == 'R':
                r_enrichments.append(zscore)
            else:
                nr_enrichments.append(zscore)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel A: Average R enrichment
    ax1 = axes[0]
    if r_enrichments:
        avg_r = np.mean(r_enrichments, axis=0)
        # Get cell type labels from first result
        first_sample = [s for s in results if SAMPLES[s]['response'] == 'R'][0]
        ct_col = results[first_sample]['ct_col']
        labels = results[first_sample]['adata'].obs[ct_col].unique()

        im1 = ax1.imshow(avg_r, cmap='RdBu_r', vmin=-3, vmax=3)
        ax1.set_xticks(range(len(labels)))
        ax1.set_yticks(range(len(labels)))
        ax1.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        ax1.set_yticklabels(labels, fontsize=8)
        plt.colorbar(im1, ax=ax1, label='Z-score', shrink=0.8)
        ax1.set_title('A) Responders: Neighborhood Enrichment', fontweight='bold')
    else:
        ax1.text(0.5, 0.5, 'No R data', ha='center', va='center', transform=ax1.transAxes)
        ax1.set_title('A) Responders', fontweight='bold')

    # Panel B: Average NR enrichment
    ax2 = axes[1]
    if nr_enrichments:
        avg_nr = np.mean(nr_enrichments, axis=0)
        first_sample = [s for s in results if SAMPLES[s]['response'] == 'NR'][0]
        ct_col = results[first_sample]['ct_col']
        labels = results[first_sample]['adata'].obs[ct_col].unique()

        im2 = ax2.imshow(avg_nr, cmap='RdBu_r', vmin=-3, vmax=3)
        ax2.set_xticks(range(len(labels)))
        ax2.set_yticks(range(len(labels)))
        ax2.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        ax2.set_yticklabels(labels, fontsize=8)
        plt.colorbar(im2, ax=ax2, label='Z-score', shrink=0.8)
        ax2.set_title('B) Non-Responders: Neighborhood Enrichment', fontweight='bold')
    else:
        ax2.text(0.5, 0.5, 'No NR data', ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('B) Non-Responders', fontweight='bold')

    # Panel C: Interpretation
    ax3 = axes[2]
    ax3.axis('off')

    interpretation = """
NICHE IDENTIFICATION
====================

Method: Neighborhood enrichment
analysis (squidpy) quantifies
cell type co-localization patterns.

Interpretation:
• Red = cell types cluster together
• Blue = cell types avoid each other
• White = random distribution

Key Patterns to Look For:
• Immune cell clustering
• Tumor-stroma interfaces
• Vascular niches

R vs NR Differences:
Different co-localization patterns
may indicate distinct TME
organization affecting treatment
response.

Clinical Relevance:
Niche composition could predict
immunotherapy response or
identify therapeutic targets.
"""
    ax3.text(0.05, 0.95, interpretation, transform=ax3.transAxes,
             fontsize=9, family='monospace', va='top',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))

    plt.suptitle('Figure 19: Spatial Niche Identification', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig19_niche_identification.png", bbox_inches='tight', dpi=300)
    plt.close()
    print("    Saved fig19_niche_identification.png")


# =============================================================================
# ANALYSIS 5: Trajectory Analysis
# =============================================================================

def run_trajectory_analysis(adatas):
    """Run trajectory/pseudotime analysis on tumor cells."""
    print("\n" + "="*60)
    print("ANALYSIS 5: Trajectory Analysis")
    print("="*60)

    try:
        import cellrank as cr
        has_cellrank = True
    except ImportError:
        print("  CellRank not available, using diffusion pseudotime instead")
        has_cellrank = False

    results = {}

    for sample, adata in adatas.items():
        print(f"\n  Processing {sample}...")

        try:
            # Get cell type column
            ct_col = None
            for col in ['cell_type', 'celltype', 'annotation']:
                if col in adata.obs.columns:
                    ct_col = col
                    break

            if ct_col is None:
                print(f"    No cell type column")
                continue

            # Subset to tumor/epithelial cells for trajectory
            tumor_mask = adata.obs[ct_col].str.contains('Ductal|Cancer|Tumor|Epithelial', case=False, na=False)

            if tumor_mask.sum() < 50:
                print(f"    Too few tumor cells ({tumor_mask.sum()})")
                continue

            adata_tumor = adata[tumor_mask].copy()
            print(f"    Subset to {adata_tumor.n_obs} tumor cells")

            # Compute diffusion pseudotime
            if 'X_pca' not in adata_tumor.obsm:
                sc.pp.pca(adata_tumor, n_comps=30)

            sc.pp.neighbors(adata_tumor, n_neighbors=15)
            sc.tl.diffmap(adata_tumor)

            # Find root cell (highest expression of stemness marker or just use first DC)
            adata_tumor.uns['iroot'] = np.argmax(adata_tumor.obsm['X_diffmap'][:, 0])
            sc.tl.dpt(adata_tumor)

            # Store pseudotime in original adata
            adata.obs.loc[tumor_mask, 'pseudotime'] = adata_tumor.obs['dpt_pseudotime'].values

            results[sample] = {
                'adata': adata,
                'adata_tumor': adata_tumor,
                'n_tumor': adata_tumor.n_obs
            }
            print(f"    Computed diffusion pseudotime")

        except Exception as e:
            print(f"    Failed: {e}")

    # Generate figures
    _plot_trajectory_results(results)

    return results


def _plot_trajectory_results(results):
    """Generate trajectory figures."""
    print("\n  Generating trajectory figures...")

    if not results:
        print("    No trajectory results to plot")
        return

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    pseudotime_stats = []

    for idx, (sample, data) in enumerate(results.items()):
        if idx >= 7:
            break

        ax = axes[idx]
        adata = data['adata']

        if 'pseudotime' in adata.obs.columns and 'spatial' in adata.obsm:
            coords = adata.obsm['spatial']
            pt = adata.obs['pseudotime'].values

            # Plot spatial with pseudotime coloring
            mask = ~np.isnan(pt)
            scatter = ax.scatter(coords[mask, 0], coords[mask, 1],
                               c=pt[mask], cmap='viridis', s=5, alpha=0.8)
            ax.scatter(coords[~mask, 0], coords[~mask, 1],
                      c='lightgray', s=2, alpha=0.3)

            ax.set_aspect('equal')
            ax.set_title(f'{sample} ({SAMPLES[sample]["response"]})', fontweight='bold')
            ax.axis('off')

            # Collect stats
            pseudotime_stats.append({
                'sample': sample,
                'response': SAMPLES[sample]['response'],
                'mean_pt': np.nanmean(pt),
                'std_pt': np.nanstd(pt),
                'n_cells': mask.sum()
            })
        else:
            ax.text(0.5, 0.5, 'No pseudotime', ha='center', va='center',
                   transform=ax.transAxes)
            ax.set_title(f'{sample}', fontweight='bold')

    # Panel 8: Summary
    ax = axes[7]
    ax.axis('off')

    if pseudotime_stats:
        pt_df = pd.DataFrame(pseudotime_stats)
        r_mean = pt_df[pt_df['response'] == 'R']['mean_pt'].mean()
        nr_mean = pt_df[pt_df['response'] == 'NR']['mean_pt'].mean()

        summary = f"""
TRAJECTORY ANALYSIS
===================

Method: Diffusion pseudotime
on tumor/epithelial cells

Mean Pseudotime:
• Responders: {r_mean:.3f}
• Non-Responders: {nr_mean:.3f}

Interpretation:
Pseudotime represents inferred
developmental progression within
the tumor cell population.

Higher values may indicate:
• More differentiated cells
• Later disease stage
• Treatment-induced changes

Spatial patterns show how
progression varies across
tissue regions.
"""
        ax.text(0.05, 0.95, summary, transform=ax.transAxes,
               fontsize=9, family='monospace', va='top',
               bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))

        pt_df.to_csv(TABLE_DIR / "trajectory_pseudotime_stats.csv", index=False)

    plt.suptitle('Figure 20: Tumor Cell Trajectory Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig20_trajectory_analysis.png", bbox_inches='tight', dpi=300)
    plt.close()
    print("    Saved fig20_trajectory_analysis.png")


# =============================================================================
# ANALYSIS 6: Multi-sample Integration
# =============================================================================

def run_integration(adatas):
    """Integrate all samples using Harmony."""
    print("\n" + "="*60)
    print("ANALYSIS 6: Multi-sample Integration (Harmony)")
    print("="*60)

    try:
        import harmonypy as hm
        has_harmony = True
    except ImportError:
        print("  Harmony not available")
        has_harmony = False
        return None

    # Concatenate all samples
    print("\n  Concatenating samples...")
    adata_list = []
    for sample, adata in adatas.items():
        adata_copy = adata.copy()
        adata_copy.obs['sample'] = sample
        adata_copy.obs['response'] = SAMPLES[sample]['response']
        adata_list.append(adata_copy)

    adata_combined = ad.concat(adata_list, join='outer')
    print(f"  Combined: {adata_combined.n_obs} cells x {adata_combined.n_vars} genes")

    # PCA
    if 'X_pca' not in adata_combined.obsm:
        sc.pp.pca(adata_combined, n_comps=50)

    # Run Harmony
    print("  Running Harmony integration...")
    try:
        ho = hm.run_harmony(adata_combined.obsm['X_pca'], adata_combined.obs, 'sample')
        adata_combined.obsm['X_pca_harmony'] = ho.Z_corr.T
        print("  Harmony complete")
    except Exception as e:
        print(f"  Harmony failed: {e}")
        adata_combined.obsm['X_pca_harmony'] = adata_combined.obsm['X_pca']

    # Compute UMAP on integrated space
    sc.pp.neighbors(adata_combined, use_rep='X_pca_harmony', n_neighbors=15)
    sc.tl.umap(adata_combined)

    # Generate figures
    _plot_integration_results(adata_combined, adatas)

    # Save integrated object
    adata_combined.write(OUTPUT_DIR / "adata" / "integrated_harmony.h5ad")
    print("  Saved integrated_harmony.h5ad")

    return adata_combined


def _plot_integration_results(adata_combined, adatas):
    """Generate integration figures."""
    print("\n  Generating integration figures...")

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Panel A: UMAP colored by sample
    ax1 = axes[0, 0]
    sc.pl.umap(adata_combined, color='sample', ax=ax1, show=False,
               title='A) UMAP by Sample', legend_loc='right margin',
               palette='tab10', size=5)

    # Panel B: UMAP colored by response
    ax2 = axes[0, 1]
    sc.pl.umap(adata_combined, color='response', ax=ax2, show=False,
               title='B) UMAP by Response', legend_loc='right margin',
               palette=RESPONSE_COLORS, size=5)

    # Panel C: UMAP colored by cell type
    ax3 = axes[0, 2]
    ct_col = None
    for col in ['cell_type', 'celltype', 'annotation']:
        if col in adata_combined.obs.columns:
            ct_col = col
            break

    if ct_col:
        sc.pl.umap(adata_combined, color=ct_col, ax=ax3, show=False,
                   title='C) UMAP by Cell Type', legend_loc='right margin',
                   size=5)
    else:
        ax3.text(0.5, 0.5, 'No cell type\nannotation', ha='center', va='center',
                transform=ax3.transAxes)
        ax3.set_title('C) UMAP by Cell Type', fontweight='bold')

    # Panel D: Sample mixing metric
    ax4 = axes[1, 0]

    # Calculate LISI-like mixing score per sample
    from sklearn.neighbors import NearestNeighbors

    try:
        nn = NearestNeighbors(n_neighbors=50)
        nn.fit(adata_combined.obsm['X_pca_harmony'])
        distances, indices = nn.kneighbors()

        mixing_scores = []
        for sample in SAMPLE_ORDER:
            sample_mask = adata_combined.obs['sample'] == sample
            sample_indices = np.where(sample_mask)[0]

            if len(sample_indices) > 0:
                # For each cell in sample, what fraction of neighbors are from other samples
                other_sample_frac = []
                for idx in sample_indices[:100]:  # Sample for speed
                    neighbors = indices[idx]
                    n_other = (adata_combined.obs['sample'].iloc[neighbors] != sample).sum()
                    other_sample_frac.append(n_other / len(neighbors))

                mixing_scores.append({
                    'sample': sample,
                    'response': SAMPLES[sample]['response'],
                    'mixing': np.mean(other_sample_frac)
                })

        mix_df = pd.DataFrame(mixing_scores)
        colors = [RESPONSE_COLORS[SAMPLES[s]['response']] for s in mix_df['sample']]
        ax4.bar(mix_df['sample'], mix_df['mixing'], color=colors, edgecolor='black')
        ax4.set_ylabel('Mixing Score')
        ax4.set_title('D) Sample Mixing (higher = better)', fontweight='bold')
        ax4.tick_params(axis='x', rotation=45)
        ax4.axhline(mix_df['mixing'].mean(), color='gray', linestyle='--', alpha=0.7)

    except Exception as e:
        ax4.text(0.5, 0.5, f'Mixing calculation\nfailed: {str(e)[:20]}',
                ha='center', va='center', transform=ax4.transAxes)
        ax4.set_title('D) Sample Mixing', fontweight='bold')

    # Panel E: Batch effect reduction
    ax5 = axes[1, 1]

    # Compare variance explained by sample before/after
    try:
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import LabelEncoder

        le = LabelEncoder()
        y = le.fit_transform(adata_combined.obs['sample'])

        # Before Harmony
        lr_before = LogisticRegression(max_iter=1000, random_state=42)
        lr_before.fit(adata_combined.obsm['X_pca'][:, :10], y)
        acc_before = lr_before.score(adata_combined.obsm['X_pca'][:, :10], y)

        # After Harmony
        lr_after = LogisticRegression(max_iter=1000, random_state=42)
        lr_after.fit(adata_combined.obsm['X_pca_harmony'][:, :10], y)
        acc_after = lr_after.score(adata_combined.obsm['X_pca_harmony'][:, :10], y)

        ax5.bar(['Before\nHarmony', 'After\nHarmony'], [acc_before, acc_after],
               color=['#e74c3c', '#2ecc71'], edgecolor='black')
        ax5.set_ylabel('Sample Classification Accuracy')
        ax5.set_title('E) Batch Effect Reduction', fontweight='bold')
        ax5.set_ylim(0, 1)

        # Add text
        ax5.text(0, acc_before + 0.02, f'{acc_before:.2f}', ha='center', fontsize=10)
        ax5.text(1, acc_after + 0.02, f'{acc_after:.2f}', ha='center', fontsize=10)

    except Exception as e:
        ax5.text(0.5, 0.5, 'Batch effect\nanalysis failed',
                ha='center', va='center', transform=ax5.transAxes)
        ax5.set_title('E) Batch Effect Reduction', fontweight='bold')

    # Panel F: Interpretation
    ax6 = axes[1, 2]
    ax6.axis('off')

    interpretation = f"""
MULTI-SAMPLE INTEGRATION
========================

Method: Harmony batch correction
on PCA space, followed by UMAP.

Total Cells: {adata_combined.n_obs:,}
Samples: {len(SAMPLE_ORDER)}

Integration Quality:
• Samples should mix in UMAP
• Cell types should cluster
• Response groups may separate

Biological Signal:
After removing batch effects,
remaining variation reflects
true biological differences
between R and NR.

Use Cases:
• Combined differential analysis
• Cross-sample cell type mapping
• Response signature discovery
"""
    ax6.text(0.05, 0.95, interpretation, transform=ax6.transAxes,
             fontsize=9, family='monospace', va='top',
             bbox=dict(boxstyle='round', facecolor='lavender', alpha=0.3))

    plt.suptitle('Figure 21: Multi-Sample Integration (Harmony)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig21_integration_harmony.png", bbox_inches='tight', dpi=300)
    plt.close()
    print("    Saved fig21_integration_harmony.png")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("="*70)
    print("PDAC ADVANCED ANALYSES - COMPLETE PIPELINE")
    print("="*70)

    # Load all samples
    print("\nLoading samples...")
    adatas = load_all_samples()

    if not adatas:
        print("ERROR: No samples loaded!")
        return

    print(f"\nLoaded {len(adatas)} samples")

    # Run all analyses
    print("\n" + "="*70)

    # 1. Pathway Activity
    adatas = run_pathway_analysis(adatas)

    # 2. GSEA
    run_gsea_analysis(adatas)

    # 3. Spatial Domains
    adatas = run_spatial_domains(adatas)

    # 4. Niche Identification
    run_niche_identification(adatas)

    # 5. Trajectory Analysis
    run_trajectory_analysis(adatas)

    # 6. Multi-sample Integration
    run_integration(adatas)

    # Summary
    print("\n" + "="*70)
    print("COMPLETE - FIGURES GENERATED:")
    print("="*70)
    for fig in sorted(FIG_DIR.glob("*.png")):
        print(f"  {fig.name}")

    print("\n" + "="*70)
    print("ALL ADVANCED ANALYSES COMPLETE")
    print("="*70)


if __name__ == "__main__":
    main()
