#!/usr/bin/env python3
"""
Create Comprehensive Excel Data Files for Manuscript
=====================================================

Generates Excel files with:
1. Sample metadata and patient structure
2. Cell type proportions (all types, pre/post labeled)
3. DE genes (extended, with stats)
4. Pathway activity (all pathways, all samples)
5. Spatial metrics (entropy, topology, etc.)
6. Niche enrichment z-scores
7. Additional R vs NR discriminators

Author: Max Van Belkum
Date: 2026-01-20
"""

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import decoupler as dc
from scipy import stats
from scipy.spatial.distance import pdist
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
ADATA_DIR = OUTPUT_DIR / "adata" / "polymathic"
TABLE_DIR = OUTPUT_DIR / "tables"
TABLE_DIR.mkdir(parents=True, exist_ok=True)

# Sample metadata with PATIENT structure
SAMPLES = {
    "YP03A": {"patient": "YP03", "response": "NR", "timepoint": "Pre", "has_pair": True},
    "YP03C": {"patient": "YP03", "response": "NR", "timepoint": "Post", "has_pair": True},
    "YP04C": {"patient": "YP04", "response": "NR", "timepoint": "Post", "has_pair": False},  # NO PRE!
    "YP12A": {"patient": "YP12", "response": "R", "timepoint": "Pre", "has_pair": True},
    "YP12C": {"patient": "YP12", "response": "R", "timepoint": "Post", "has_pair": True},
    "YP15A": {"patient": "YP15", "response": "R", "timepoint": "Pre", "has_pair": True},
    "YP15C": {"patient": "YP15", "response": "R", "timepoint": "Post", "has_pair": True},
}


def load_all_samples():
    """Load all polymathic h5ad files."""
    adatas = {}
    for sample in SAMPLES.keys():
        path = ADATA_DIR / f"{sample}_polymathic.h5ad"
        if path.exists():
            adata = sc.read_h5ad(path)
            for key, value in SAMPLES[sample].items():
                adata.obs[key] = value
            adata.obs['sample'] = sample
            adatas[sample] = adata
            print(f"  Loaded {sample}: {adata.n_obs} spots")
    return adatas


def create_sample_metadata_excel(adatas):
    """Create sample metadata Excel with patient structure."""
    print("\n" + "="*60)
    print("1. SAMPLE METADATA")
    print("="*60)

    rows = []
    for sample, adata in adatas.items():
        info = SAMPLES[sample]

        # Basic stats
        n_spots = adata.n_obs
        n_genes = adata.n_vars
        total_counts = adata.obs['total_counts'].sum() if 'total_counts' in adata.obs else np.nan
        mean_counts = adata.obs['total_counts'].mean() if 'total_counts' in adata.obs else np.nan
        median_genes = adata.obs['n_genes_by_counts'].median() if 'n_genes_by_counts' in adata.obs else np.nan

        rows.append({
            'sample_id': sample,
            'patient_id': info['patient'],
            'response': info['response'],
            'response_full': 'Responder' if info['response'] == 'R' else 'Non-Responder',
            'timepoint': info['timepoint'],
            'timepoint_full': 'Pre-treatment' if info['timepoint'] == 'Pre' else 'Post-treatment',
            'has_paired_sample': info['has_pair'],
            'n_spots': n_spots,
            'n_genes_detected': n_genes,
            'total_UMI_counts': total_counts,
            'mean_counts_per_spot': mean_counts,
            'median_genes_per_spot': median_genes,
        })

    df = pd.DataFrame(rows)

    # Add patient-level summary
    patient_summary = []
    for patient in ['YP03', 'YP04', 'YP12', 'YP15']:
        patient_samples = [s for s in SAMPLES if SAMPLES[s]['patient'] == patient]
        patient_summary.append({
            'patient_id': patient,
            'response': SAMPLES[patient_samples[0]]['response'],
            'n_samples': len(patient_samples),
            'has_pre': any(SAMPLES[s]['timepoint'] == 'Pre' for s in patient_samples),
            'has_post': any(SAMPLES[s]['timepoint'] == 'Post' for s in patient_samples),
            'samples': ', '.join(patient_samples)
        })

    patient_df = pd.DataFrame(patient_summary)

    # Save to Excel with multiple sheets
    excel_path = TABLE_DIR / 'S1_sample_metadata.xlsx'
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Sample_Details', index=False)
        patient_df.to_excel(writer, sheet_name='Patient_Summary', index=False)

        # Add notes sheet
        notes = pd.DataFrame({
            'Note': [
                'IMPORTANT: YP04 only has post-treatment sample (no pre-treatment pair)',
                'Total patients: 4 (2 Responders, 2 Non-Responders)',
                'Total samples: 7 (4 Responders, 3 Non-Responders)',
                'Paired analysis (Pre vs Post) limited to: YP03, YP12, YP15',
                'Sample naming: A suffix = Pre-treatment, C suffix = Post-treatment'
            ]
        })
        notes.to_excel(writer, sheet_name='Notes', index=False)

    print(f"  Saved: {excel_path}")
    return df


def create_cell_type_excel(adatas):
    """Create cell type proportions Excel."""
    print("\n" + "="*60)
    print("2. CELL TYPE PROPORTIONS")
    print("="*60)

    # Collect proportions
    all_props = []

    for sample, adata in adatas.items():
        info = SAMPLES[sample]

        if 'cell_type' in adata.obs.columns:
            ct_col = 'cell_type'
        elif 'cell_type_broad' in adata.obs.columns:
            ct_col = 'cell_type_broad'
        else:
            continue

        counts = adata.obs[ct_col].value_counts()
        total = counts.sum()

        for ct, count in counts.items():
            all_props.append({
                'sample_id': sample,
                'patient_id': info['patient'],
                'response': info['response'],
                'timepoint': info['timepoint'],
                'cell_type': ct,
                'count': count,
                'proportion': count / total,
                'percentage': 100 * count / total
            })

    df = pd.DataFrame(all_props)

    # Create pivot table for easy viewing
    pivot = df.pivot_table(
        index='cell_type',
        columns='sample_id',
        values='percentage',
        aggfunc='first'
    )

    # Reorder columns: R samples first, then NR
    r_samples = [s for s in SAMPLES if SAMPLES[s]['response'] == 'R']
    nr_samples = [s for s in SAMPLES if SAMPLES[s]['response'] == 'NR']
    col_order = [s for s in r_samples + nr_samples if s in pivot.columns]
    pivot = pivot[col_order]

    # Calculate R vs NR statistics
    ct_stats = []
    for ct in df['cell_type'].unique():
        ct_data = df[df['cell_type'] == ct]
        r_vals = ct_data[ct_data['response'] == 'R']['percentage'].values
        nr_vals = ct_data[ct_data['response'] == 'NR']['percentage'].values

        if len(r_vals) >= 2 and len(nr_vals) >= 2:
            t_stat, p_val = stats.ttest_ind(r_vals, nr_vals, equal_var=False)
            u_stat, mwu_p = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')

            pooled_std = np.sqrt((np.var(r_vals) + np.var(nr_vals)) / 2)
            cohens_d = (np.mean(r_vals) - np.mean(nr_vals)) / pooled_std if pooled_std > 0 else 0
        else:
            t_stat, p_val, u_stat, mwu_p, cohens_d = np.nan, np.nan, np.nan, np.nan, np.nan

        ct_stats.append({
            'cell_type': ct,
            'R_mean': np.mean(r_vals),
            'R_std': np.std(r_vals),
            'R_n': len(r_vals),
            'NR_mean': np.mean(nr_vals),
            'NR_std': np.std(nr_vals),
            'NR_n': len(nr_vals),
            'difference': np.mean(r_vals) - np.mean(nr_vals),
            'cohens_d': cohens_d,
            'welch_t': t_stat,
            'welch_p': p_val,
            'mwu_U': u_stat,
            'mwu_p': mwu_p,
            'significant_0.05': p_val < 0.05 if not np.isnan(p_val) else False
        })

    stats_df = pd.DataFrame(ct_stats).sort_values('cohens_d', ascending=False)

    # Pre vs Post analysis (paired samples only)
    paired_stats = []
    for ct in df['cell_type'].unique():
        for patient in ['YP03', 'YP12', 'YP15']:  # Exclude YP04 (no pre)
            ct_patient = df[(df['cell_type'] == ct) & (df['patient_id'] == patient)]
            pre = ct_patient[ct_patient['timepoint'] == 'Pre']['percentage'].values
            post = ct_patient[ct_patient['timepoint'] == 'Post']['percentage'].values

            if len(pre) == 1 and len(post) == 1:
                paired_stats.append({
                    'cell_type': ct,
                    'patient_id': patient,
                    'response': SAMPLES[f'{patient}A']['response'],
                    'pre_pct': pre[0],
                    'post_pct': post[0],
                    'change': post[0] - pre[0],
                    'fold_change': post[0] / pre[0] if pre[0] > 0 else np.nan
                })

    paired_df = pd.DataFrame(paired_stats)

    # Save
    excel_path = TABLE_DIR / 'S2_cell_type_proportions.xlsx'
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Raw_Data', index=False)
        pivot.to_excel(writer, sheet_name='Proportions_Matrix')
        stats_df.to_excel(writer, sheet_name='R_vs_NR_Statistics', index=False)
        paired_df.to_excel(writer, sheet_name='Pre_vs_Post_Paired', index=False)

    print(f"  Saved: {excel_path}")
    print(f"  Top discriminating cell types:")
    for _, row in stats_df.head(5).iterrows():
        sig = '*' if row['significant_0.05'] else ''
        print(f"    {row['cell_type']}: d={row['cohens_d']:.2f}, p={row['welch_p']:.3f}{sig}")

    return df, stats_df


def create_de_genes_excel(adatas):
    """Create DE genes Excel with extended results."""
    print("\n" + "="*60)
    print("3. DIFFERENTIAL EXPRESSION GENES")
    print("="*60)

    # Concatenate all samples
    adata_list = []
    for sample, adata in adatas.items():
        adata_copy = adata.copy()
        adata_list.append(adata_copy)

    combined = sc.concat(adata_list, join='outer')

    # Run DE: R vs NR (all samples)
    print("  Computing R vs NR DE (all samples)...")
    sc.tl.rank_genes_groups(combined, groupby='response', method='wilcoxon',
                           key_added='de_response')

    # Extract results for both groups
    de_results = []
    for group in ['R', 'NR']:
        names = combined.uns['de_response']['names'][group]
        scores = combined.uns['de_response']['scores'][group]
        pvals = combined.uns['de_response']['pvals'][group]
        pvals_adj = combined.uns['de_response']['pvals_adj'][group]
        logfc = combined.uns['de_response']['logfoldchanges'][group]

        for i in range(min(500, len(names))):  # Top 500 genes per group
            de_results.append({
                'gene': names[i],
                'upregulated_in': group,
                'score': scores[i],
                'log2FC': logfc[i],
                'pval': pvals[i],
                'pval_adj': pvals_adj[i],
                'significant_0.05': pvals_adj[i] < 0.05,
                'significant_0.01': pvals_adj[i] < 0.01,
            })

    de_df = pd.DataFrame(de_results)

    # Pre vs Post DE (responders only)
    print("  Computing Pre vs Post DE (Responders)...")
    r_samples = [s for s, info in SAMPLES.items() if info['response'] == 'R']
    r_combined = sc.concat([adatas[s] for s in r_samples if s in adatas], join='outer')

    sc.tl.rank_genes_groups(r_combined, groupby='timepoint', method='wilcoxon',
                           key_added='de_timepoint')

    pre_post_results = []
    for group in ['Pre', 'Post']:
        names = r_combined.uns['de_timepoint']['names'][group]
        scores = r_combined.uns['de_timepoint']['scores'][group]
        pvals = r_combined.uns['de_timepoint']['pvals'][group]
        pvals_adj = r_combined.uns['de_timepoint']['pvals_adj'][group]
        logfc = r_combined.uns['de_timepoint']['logfoldchanges'][group]

        for i in range(min(200, len(names))):
            pre_post_results.append({
                'gene': names[i],
                'upregulated_in': group,
                'score': scores[i],
                'log2FC': logfc[i],
                'pval': pvals[i],
                'pval_adj': pvals_adj[i],
            })

    prepost_df = pd.DataFrame(pre_post_results)

    # Per-sample mean expression of top DE genes
    print("  Computing per-sample expression of top genes...")
    top_genes = de_df.nlargest(100, 'score')['gene'].tolist()

    expr_data = []
    for sample, adata in adatas.items():
        info = SAMPLES[sample]
        for gene in top_genes:
            if gene in adata.var_names:
                expr = adata[:, gene].X.mean()
                expr_data.append({
                    'gene': gene,
                    'sample_id': sample,
                    'patient_id': info['patient'],
                    'response': info['response'],
                    'timepoint': info['timepoint'],
                    'mean_expression': float(expr)
                })

    expr_df = pd.DataFrame(expr_data)

    # Save
    excel_path = TABLE_DIR / 'S3_differential_expression.xlsx'
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        de_df.to_excel(writer, sheet_name='R_vs_NR_Top500', index=False)
        prepost_df.to_excel(writer, sheet_name='Pre_vs_Post_Responders', index=False)
        expr_df.to_excel(writer, sheet_name='Top_Gene_Expression', index=False)

        # Summary
        summary = pd.DataFrame({
            'Comparison': ['R vs NR', 'R vs NR', 'Pre vs Post (R only)', 'Pre vs Post (R only)'],
            'Group': ['Up in R', 'Up in NR', 'Up in Pre', 'Up in Post'],
            'Significant_genes_0.05': [
                len(de_df[(de_df['upregulated_in'] == 'R') & (de_df['significant_0.05'])]),
                len(de_df[(de_df['upregulated_in'] == 'NR') & (de_df['significant_0.05'])]),
                len(prepost_df[(prepost_df['upregulated_in'] == 'Pre') & (prepost_df['pval_adj'] < 0.05)]),
                len(prepost_df[(prepost_df['upregulated_in'] == 'Post') & (prepost_df['pval_adj'] < 0.05)]),
            ]
        })
        summary.to_excel(writer, sheet_name='Summary', index=False)

    print(f"  Saved: {excel_path}")
    print(f"  Significant genes (adj.p<0.05): Up in R={len(de_df[(de_df['upregulated_in']=='R') & de_df['significant_0.05']])}, Up in NR={len(de_df[(de_df['upregulated_in']=='NR') & de_df['significant_0.05']])}")

    return de_df


def create_pathway_excel(adatas):
    """Create pathway activity Excel."""
    print("\n" + "="*60)
    print("4. PATHWAY ACTIVITY")
    print("="*60)

    # Get PROGENy model
    progeny = dc.get_progeny(organism='human', top=500)
    progeny['weight'] = progeny['weight'].astype('float32')
    progeny['p_value'] = progeny['p_value'].astype('float32')

    # Run for all samples
    pathway_data = []

    for sample, adata in adatas.items():
        info = SAMPLES[sample]
        print(f"  Processing {sample}...")

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

            activities = adata.obsm['mlm_estimate']
            pvals = adata.obsm['mlm_pvals']

            for pathway in activities.columns:
                pathway_data.append({
                    'sample_id': sample,
                    'patient_id': info['patient'],
                    'response': info['response'],
                    'timepoint': info['timepoint'],
                    'pathway': pathway,
                    'mean_activity': activities[pathway].mean(),
                    'std_activity': activities[pathway].std(),
                    'median_activity': activities[pathway].median(),
                    'mean_pval': pvals[pathway].mean(),
                })
        except Exception as e:
            print(f"    Failed: {e}")

    df = pd.DataFrame(pathway_data)

    # Calculate R vs NR statistics
    pathway_stats = []
    for pathway in df['pathway'].unique():
        pw_data = df[df['pathway'] == pathway]
        r_vals = pw_data[pw_data['response'] == 'R']['mean_activity'].values
        nr_vals = pw_data[pw_data['response'] == 'NR']['mean_activity'].values

        if len(r_vals) >= 2 and len(nr_vals) >= 2:
            t_stat, p_val = stats.ttest_ind(r_vals, nr_vals, equal_var=False)
            u_stat, mwu_p = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')

            pooled_std = np.sqrt((np.var(r_vals) + np.var(nr_vals)) / 2)
            cohens_d = (np.mean(r_vals) - np.mean(nr_vals)) / pooled_std if pooled_std > 0 else 0
        else:
            t_stat, p_val, u_stat, mwu_p, cohens_d = np.nan, np.nan, np.nan, np.nan, np.nan

        pathway_stats.append({
            'pathway': pathway,
            'R_mean': np.mean(r_vals),
            'R_std': np.std(r_vals),
            'NR_mean': np.mean(nr_vals),
            'NR_std': np.std(nr_vals),
            'difference': np.mean(r_vals) - np.mean(nr_vals),
            'cohens_d': cohens_d,
            'welch_t': t_stat,
            'welch_p': p_val,
            'mwu_p': mwu_p,
            'significant_0.05': p_val < 0.05 if not np.isnan(p_val) else False
        })

    stats_df = pd.DataFrame(pathway_stats).sort_values('cohens_d', key=abs, ascending=False)

    # Pivot table
    pivot = df.pivot_table(
        index='pathway',
        columns='sample_id',
        values='mean_activity'
    )

    # Save
    excel_path = TABLE_DIR / 'S4_pathway_activity.xlsx'
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Raw_Data', index=False)
        pivot.to_excel(writer, sheet_name='Activity_Matrix')
        stats_df.to_excel(writer, sheet_name='R_vs_NR_Statistics', index=False)

    print(f"  Saved: {excel_path}")
    print(f"  Significant pathways (p<0.05):")
    for _, row in stats_df[stats_df['significant_0.05']].iterrows():
        print(f"    {row['pathway']}: d={row['cohens_d']:.2f}, p={row['welch_p']:.3f}")

    return df, stats_df


def create_spatial_metrics_excel(adatas):
    """Create spatial metrics Excel."""
    print("\n" + "="*60)
    print("5. SPATIAL METRICS")
    print("="*60)

    metrics_data = []

    for sample, adata in adatas.items():
        info = SAMPLES[sample]
        print(f"  Processing {sample}...")

        row = {
            'sample_id': sample,
            'patient_id': info['patient'],
            'response': info['response'],
            'timepoint': info['timepoint'],
            'n_spots': adata.n_obs,
        }

        # Spatial entropy (if computed)
        if 'spatial_entropy' in adata.obs.columns:
            row['spatial_entropy_mean'] = adata.obs['spatial_entropy'].mean()
            row['spatial_entropy_std'] = adata.obs['spatial_entropy'].std()

        # Compute new spatial metrics
        if 'spatial' in adata.obsm:
            coords = adata.obsm['spatial']

            # Tissue area (convex hull approximation)
            from scipy.spatial import ConvexHull
            try:
                hull = ConvexHull(coords)
                row['tissue_area'] = hull.volume  # 2D area
                row['tissue_perimeter'] = hull.area  # 2D perimeter
            except:
                pass

            # Spot density
            if 'tissue_area' in row and row['tissue_area'] > 0:
                row['spot_density'] = adata.n_obs / row['tissue_area']

            # Mean nearest neighbor distance
            from scipy.spatial.distance import cdist
            dists = cdist(coords, coords)
            np.fill_diagonal(dists, np.inf)
            nn_dists = dists.min(axis=1)
            row['mean_nn_distance'] = np.mean(nn_dists)
            row['std_nn_distance'] = np.std(nn_dists)

        # Morans I for top genes (spatial autocorrelation)
        if 'spatial_neighbors' in adata.uns or 'spatial' in adata.obsm:
            try:
                if 'spatial_neighbors' not in adata.uns:
                    sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

                # Calculate Moran's I for a few marker genes
                test_genes = ['KRT19', 'EPCAM', 'CD3D', 'CD68', 'COL1A1']
                for gene in test_genes:
                    if gene in adata.var_names:
                        sq.gr.spatial_autocorr(adata, genes=[gene], mode='moran')
                        if 'moranI' in adata.uns:
                            moran_val = adata.uns['moranI'].loc[gene, 'I']
                            row[f'moranI_{gene}'] = moran_val
            except Exception as e:
                print(f"    Moran's I failed: {e}")

        # Cell type diversity (Shannon)
        if 'cell_type' in adata.obs.columns:
            ct_counts = adata.obs['cell_type'].value_counts()
            ct_props = ct_counts / ct_counts.sum()
            row['ct_diversity_shannon'] = -np.sum(ct_props * np.log(ct_props + 1e-10))
            row['ct_richness'] = len(ct_counts)

        metrics_data.append(row)

    df = pd.DataFrame(metrics_data)

    # Calculate R vs NR statistics for each metric
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    metric_stats = []

    for col in numeric_cols:
        if col in ['n_spots']:
            continue

        r_vals = df[df['response'] == 'R'][col].dropna().values
        nr_vals = df[df['response'] == 'NR'][col].dropna().values

        if len(r_vals) >= 2 and len(nr_vals) >= 2:
            t_stat, p_val = stats.ttest_ind(r_vals, nr_vals, equal_var=False)
            pooled_std = np.sqrt((np.var(r_vals) + np.var(nr_vals)) / 2)
            cohens_d = (np.mean(r_vals) - np.mean(nr_vals)) / pooled_std if pooled_std > 0 else 0

            metric_stats.append({
                'metric': col,
                'R_mean': np.mean(r_vals),
                'R_std': np.std(r_vals),
                'NR_mean': np.mean(nr_vals),
                'NR_std': np.std(nr_vals),
                'cohens_d': cohens_d,
                'welch_p': p_val,
                'significant_0.05': p_val < 0.05
            })

    stats_df = pd.DataFrame(metric_stats).sort_values('cohens_d', key=abs, ascending=False)

    # Save
    excel_path = TABLE_DIR / 'S5_spatial_metrics.xlsx'
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='All_Metrics', index=False)
        stats_df.to_excel(writer, sheet_name='R_vs_NR_Statistics', index=False)

    print(f"  Saved: {excel_path}")
    return df, stats_df


def create_niche_excel(adatas):
    """Create niche enrichment Excel."""
    print("\n" + "="*60)
    print("6. NICHE ENRICHMENT")
    print("="*60)

    niche_data = []

    for sample, adata in adatas.items():
        info = SAMPLES[sample]
        print(f"  Processing {sample}...")

        ct_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

        if 'spatial_neighbors' not in adata.uns:
            sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

        try:
            sq.gr.nhood_enrichment(adata, cluster_key=ct_col)

            zscore = adata.uns[f'{ct_col}_nhood_enrichment']['zscore']
            cell_types = adata.obs[ct_col].cat.categories.tolist()

            for i, ct1 in enumerate(cell_types):
                for j, ct2 in enumerate(cell_types):
                    niche_data.append({
                        'sample_id': sample,
                        'patient_id': info['patient'],
                        'response': info['response'],
                        'timepoint': info['timepoint'],
                        'cell_type_1': ct1,
                        'cell_type_2': ct2,
                        'zscore': zscore[i, j],
                        'interaction_type': 'self' if i == j else 'hetero'
                    })
        except Exception as e:
            print(f"    Failed: {e}")

    df = pd.DataFrame(niche_data)

    # Calculate mean z-scores by interaction
    interaction_summary = df.groupby(['cell_type_1', 'cell_type_2', 'response']).agg({
        'zscore': ['mean', 'std', 'count']
    }).reset_index()
    interaction_summary.columns = ['cell_type_1', 'cell_type_2', 'response',
                                   'mean_zscore', 'std_zscore', 'n_samples']

    # R vs NR statistics for each interaction
    interaction_stats = []
    for (ct1, ct2), group in df.groupby(['cell_type_1', 'cell_type_2']):
        r_vals = group[group['response'] == 'R']['zscore'].dropna().values
        nr_vals = group[group['response'] == 'NR']['zscore'].dropna().values

        if len(r_vals) >= 2 and len(nr_vals) >= 2:
            r_vals = r_vals[~np.isinf(r_vals)]
            nr_vals = nr_vals[~np.isinf(nr_vals)]

            if len(r_vals) >= 2 and len(nr_vals) >= 2:
                t_stat, p_val = stats.ttest_ind(r_vals, nr_vals, equal_var=False)
                pooled_std = np.sqrt((np.var(r_vals) + np.var(nr_vals)) / 2)
                cohens_d = (np.mean(r_vals) - np.mean(nr_vals)) / pooled_std if pooled_std > 0 else 0

                interaction_stats.append({
                    'cell_type_1': ct1,
                    'cell_type_2': ct2,
                    'R_mean': np.mean(r_vals),
                    'NR_mean': np.mean(nr_vals),
                    'cohens_d': cohens_d,
                    'welch_p': p_val,
                    'significant_0.05': p_val < 0.05
                })

    stats_df = pd.DataFrame(interaction_stats).sort_values('cohens_d', key=abs, ascending=False)

    # Save
    excel_path = TABLE_DIR / 'S6_niche_enrichment.xlsx'
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Raw_Data', index=False)
        interaction_summary.to_excel(writer, sheet_name='Interaction_Summary', index=False)
        stats_df.to_excel(writer, sheet_name='R_vs_NR_Statistics', index=False)

    print(f"  Saved: {excel_path}")
    print(f"  Significant interactions (p<0.05): {len(stats_df[stats_df['significant_0.05']])}")

    return df, stats_df


def run_additional_discriminators(adatas):
    """Run additional analyses to discriminate R vs NR."""
    print("\n" + "="*60)
    print("7. ADDITIONAL R vs NR DISCRIMINATORS")
    print("="*60)

    results = []

    for sample, adata in adatas.items():
        info = SAMPLES[sample]
        print(f"  Processing {sample}...")

        row = {
            'sample_id': sample,
            'patient_id': info['patient'],
            'response': info['response'],
            'timepoint': info['timepoint'],
        }

        # 1. Immune infiltration score
        immune_markers = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'PTPRC']
        immune_genes = [g for g in immune_markers if g in adata.var_names]
        if immune_genes:
            row['immune_score'] = adata[:, immune_genes].X.mean()

        # 2. Tumor cell proportion (epithelial markers)
        tumor_markers = ['KRT19', 'EPCAM', 'KRT7', 'KRT8', 'MUC1']
        tumor_genes = [g for g in tumor_markers if g in adata.var_names]
        if tumor_genes:
            row['tumor_score'] = adata[:, tumor_genes].X.mean()

        # 3. Stromal score
        stromal_markers = ['COL1A1', 'COL1A2', 'ACTA2', 'FAP', 'PDGFRA', 'PDGFRB']
        stromal_genes = [g for g in stromal_markers if g in adata.var_names]
        if stromal_genes:
            row['stromal_score'] = adata[:, stromal_genes].X.mean()

        # 4. Angiogenesis score
        angio_markers = ['VEGFA', 'VEGFB', 'VEGFC', 'KDR', 'FLT1', 'PECAM1', 'CDH5']
        angio_genes = [g for g in angio_markers if g in adata.var_names]
        if angio_genes:
            row['angiogenesis_score'] = adata[:, angio_genes].X.mean()

        # 5. Proliferation score
        prolif_markers = ['MKI67', 'TOP2A', 'PCNA', 'MCM2', 'MCM6']
        prolif_genes = [g for g in prolif_markers if g in adata.var_names]
        if prolif_genes:
            row['proliferation_score'] = adata[:, prolif_genes].X.mean()

        # 6. Hypoxia score
        hypoxia_markers = ['HIF1A', 'VEGFA', 'SLC2A1', 'LDHA', 'PGK1', 'ENO1']
        hypoxia_genes = [g for g in hypoxia_markers if g in adata.var_names]
        if hypoxia_genes:
            row['hypoxia_score'] = adata[:, hypoxia_genes].X.mean()

        # 7. KRAS signature (important for PDAC)
        kras_markers = ['KRAS', 'HRAS', 'NRAS', 'BRAF', 'RAF1', 'MAP2K1', 'MAPK1']
        kras_genes = [g for g in kras_markers if g in adata.var_names]
        if kras_genes:
            row['kras_pathway_score'] = adata[:, kras_genes].X.mean()

        # 8. TGF-beta score
        tgfb_markers = ['TGFB1', 'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2', 'SMAD2', 'SMAD3']
        tgfb_genes = [g for g in tgfb_markers if g in adata.var_names]
        if tgfb_genes:
            row['tgfb_score'] = adata[:, tgfb_genes].X.mean()

        # 9. Interferon response
        ifn_markers = ['ISG15', 'MX1', 'OAS1', 'IFIT1', 'IFIT2', 'IRF7', 'STAT1']
        ifn_genes = [g for g in ifn_markers if g in adata.var_names]
        if ifn_genes:
            row['ifn_response_score'] = adata[:, ifn_genes].X.mean()

        # 10. Cytotoxicity score
        cyto_markers = ['GZMA', 'GZMB', 'GZMK', 'PRF1', 'NKG7', 'GNLY']
        cyto_genes = [g for g in cyto_markers if g in adata.var_names]
        if cyto_genes:
            row['cytotoxicity_score'] = adata[:, cyto_genes].X.mean()

        results.append(row)

    df = pd.DataFrame(results)

    # Calculate R vs NR statistics
    score_cols = [c for c in df.columns if c.endswith('_score')]
    score_stats = []

    for col in score_cols:
        r_vals = df[df['response'] == 'R'][col].dropna().values.astype(float)
        nr_vals = df[df['response'] == 'NR'][col].dropna().values.astype(float)

        if len(r_vals) >= 2 and len(nr_vals) >= 2:
            t_stat, p_val = stats.ttest_ind(r_vals, nr_vals, equal_var=False)
            pooled_std = np.sqrt((np.var(r_vals) + np.var(nr_vals)) / 2)
            cohens_d = (np.mean(r_vals) - np.mean(nr_vals)) / pooled_std if pooled_std > 0 else 0

            score_stats.append({
                'score': col,
                'R_mean': np.mean(r_vals),
                'R_std': np.std(r_vals),
                'NR_mean': np.mean(nr_vals),
                'NR_std': np.std(nr_vals),
                'cohens_d': cohens_d,
                'welch_p': p_val,
                'discriminates_R': cohens_d > 0,  # True if higher in R
                'significant_0.05': p_val < 0.05
            })

    stats_df = pd.DataFrame(score_stats).sort_values('cohens_d', key=abs, ascending=False)

    # Save
    excel_path = TABLE_DIR / 'S7_discriminator_scores.xlsx'
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Raw_Scores', index=False)
        stats_df.to_excel(writer, sheet_name='R_vs_NR_Statistics', index=False)

        # Add gene lists
        gene_lists = pd.DataFrame({
            'Score': ['immune', 'tumor', 'stromal', 'angiogenesis', 'proliferation',
                     'hypoxia', 'kras_pathway', 'tgfb', 'ifn_response', 'cytotoxicity'],
            'Genes': [
                ', '.join(immune_markers),
                ', '.join(tumor_markers),
                ', '.join(stromal_markers),
                ', '.join(angio_markers),
                ', '.join(prolif_markers),
                ', '.join(hypoxia_markers),
                ', '.join(kras_markers),
                ', '.join(tgfb_markers),
                ', '.join(ifn_markers),
                ', '.join(cyto_markers)
            ]
        })
        gene_lists.to_excel(writer, sheet_name='Gene_Lists', index=False)

    print(f"  Saved: {excel_path}")
    print(f"\n  TOP DISCRIMINATORS (by effect size):")
    for _, row in stats_df.head(10).iterrows():
        direction = "↑R" if row['discriminates_R'] else "↑NR"
        sig = "*" if row['significant_0.05'] else ""
        print(f"    {row['score']}: d={row['cohens_d']:.2f} ({direction}), p={row['welch_p']:.3f}{sig}")

    return df, stats_df


def main():
    print("="*60)
    print("CREATING COMPREHENSIVE EXCEL DATA FILES")
    print("="*60)

    print("\nLoading samples...")
    adatas = load_all_samples()

    if not adatas:
        print("ERROR: No samples found!")
        return

    # 1. Sample metadata
    create_sample_metadata_excel(adatas)

    # 2. Cell type proportions
    create_cell_type_excel(adatas)

    # 3. DE genes
    create_de_genes_excel(adatas)

    # 4. Pathway activity
    create_pathway_excel(adatas)

    # 5. Spatial metrics
    create_spatial_metrics_excel(adatas)

    # 6. Niche enrichment
    create_niche_excel(adatas)

    # 7. Additional discriminators
    run_additional_discriminators(adatas)

    print("\n" + "="*60)
    print("ALL EXCEL FILES CREATED")
    print("="*60)
    print(f"\nFiles saved to: {TABLE_DIR}")
    for f in sorted(TABLE_DIR.glob('S*.xlsx')):
        print(f"  {f.name}")


if __name__ == '__main__':
    main()
