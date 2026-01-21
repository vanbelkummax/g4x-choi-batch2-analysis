#!/usr/bin/env python3
"""
Deep-Dive Statistical Analysis for Spatial Hackathon 2026
=========================================================

Performs comprehensive statistical analysis for p < 0.1 findings:
- Bootstrap confidence intervals (1000 iterations)
- Leave-one-out sensitivity analysis
- Dual statistical testing (MWU + Welch)

Metrics analyzed:
1. Betweenness centrality (from polymathic data)
2. Betti-1 AUC (topology - recomputed)
3. Acinar cell proportion (from cell type data)
4. Spatial entropy (computed from cell type distributions)

Author: Max Van Belkum
Date: 2026-01-20
"""

import matplotlib
matplotlib.use('Agg')

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import entropy
import warnings
warnings.filterwarnings('ignore')

# Import our statistical framework
from statistical_framework import (
    dual_test, bootstrap_ci, leave_one_out_sensitivity,
    format_stats_annotation, format_effect_annotation,
    interpret_significance, apply_fdr_correction
)

# Try TDA imports
try:
    from gtda.homology import VietorisRipsPersistence
    from gtda.diagrams import BettiCurve
    HAS_TDA = True
except ImportError:
    HAS_TDA = False
    print("Warning: giotto-tda not available, skipping topology analysis")

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
ADATA_DIR = OUTPUT_DIR / "adata" / "polymathic"
TABLE_DIR = OUTPUT_DIR / "tables"
FIG_DIR = OUTPUT_DIR / "figures" / "showcase_v2"

# Sample metadata
SAMPLES = {
    "YP03A": {"response": "NR", "timepoint": "Pre"},
    "YP03C": {"response": "NR", "timepoint": "Post"},
    "YP04C": {"response": "NR", "timepoint": "Post"},
    "YP12A": {"response": "R", "timepoint": "Pre"},
    "YP12C": {"response": "R", "timepoint": "Post"},
    "YP15A": {"response": "R", "timepoint": "Pre"},
    "YP15C": {"response": "R", "timepoint": "Post"},
}


def compute_spatial_entropy(adata, cell_type_col='cell_type'):
    """Compute spatial entropy based on cell type distribution."""
    if cell_type_col not in adata.obs.columns:
        return np.nan

    # Get cell type proportions
    ct_counts = adata.obs[cell_type_col].value_counts(normalize=True)

    # Shannon entropy (higher = more diverse)
    return entropy(ct_counts.values, base=2)


def compute_betti1_auc(coords, max_points=2000, seed=42):
    """Compute Betti-1 AUC from spatial coordinates using TDA."""
    if not HAS_TDA:
        return np.nan

    rng = np.random.default_rng(seed)

    # Subsample if needed (MaxMin for topology preservation)
    if len(coords) > max_points:
        # Simple random for speed in this analysis
        idx = rng.choice(len(coords), max_points, replace=False)
        coords = coords[idx]

    try:
        # Compute persistence
        VR = VietorisRipsPersistence(
            metric='euclidean',
            homology_dimensions=[0, 1],
            n_jobs=1
        )
        diagrams = VR.fit_transform([coords])

        # Compute Betti curves
        betti = BettiCurve(n_bins=100)
        curves = betti.fit_transform(diagrams)

        # AUC of Betti-1 curve (holes/loops)
        betti1_curve = curves[0, :, 1]  # Dimension 1
        auc = np.trapz(betti1_curve)

        return auc
    except Exception as e:
        print(f"  TDA error: {e}")
        return np.nan


def extract_all_metrics():
    """Extract all metrics for each sample."""
    print("Extracting metrics from all samples...")

    # Load cell type proportions
    props_df = pd.read_csv(TABLE_DIR / "cell_type_proportions.csv")

    results = []

    for sample, meta in SAMPLES.items():
        print(f"  Processing {sample}...")

        # Load adata
        adata_path = ADATA_DIR / f"{sample}_polymathic.h5ad"
        if not adata_path.exists():
            print(f"    Warning: {sample} not found")
            continue

        adata = sc.read_h5ad(adata_path)

        # 1. Betweenness centrality (mean per sample)
        betweenness = adata.obs['betweenness'].mean() if 'betweenness' in adata.obs else np.nan

        # 2. Betti-1 AUC
        coords = adata.obsm['spatial']
        betti1_auc = compute_betti1_auc(coords)

        # 3. Acinar proportion
        props_row = props_df[props_df['sample'] == sample]
        acinar = props_row['Acinar'].values[0] if len(props_row) > 0 and 'Acinar' in props_row.columns else np.nan

        # 4. Spatial entropy
        spatial_entropy = compute_spatial_entropy(adata)

        results.append({
            'sample': sample,
            'response': meta['response'],
            'timepoint': meta['timepoint'],
            'n_cells': adata.n_obs,
            'betweenness': betweenness,
            'betti1_auc': betti1_auc,
            'acinar_prop': acinar,
            'spatial_entropy': spatial_entropy,
        })

    return pd.DataFrame(results)


def run_deep_dive_analysis(metrics_df):
    """Run comprehensive statistical analysis on all metrics."""
    print("\nRunning deep-dive statistical analysis...")

    # Separate by response
    R_df = metrics_df[metrics_df['response'] == 'R']
    NR_df = metrics_df[metrics_df['response'] == 'NR']

    metrics_to_analyze = ['betweenness', 'betti1_auc', 'acinar_prop', 'spatial_entropy']
    metric_names = {
        'betweenness': 'Betweenness Centrality',
        'betti1_auc': 'Betti-1 AUC',
        'acinar_prop': 'Acinar Proportion',
        'spatial_entropy': 'Spatial Entropy'
    }

    results = []

    for metric in metrics_to_analyze:
        print(f"\n  Analyzing: {metric_names[metric]}")

        R_vals = R_df[metric].dropna().values
        NR_vals = NR_df[metric].dropna().values

        if len(R_vals) < 2 or len(NR_vals) < 2:
            print(f"    Insufficient data (R={len(R_vals)}, NR={len(NR_vals)})")
            continue

        # Dual test
        dual_result = dual_test(R_vals, NR_vals)

        # Bootstrap CI
        ci_result = bootstrap_ci(R_vals, NR_vals, n_bootstrap=1000)

        # LOO sensitivity
        R_labels = R_df[R_df[metric].notna()]['sample'].tolist()
        NR_labels = NR_df[NR_df[metric].notna()]['sample'].tolist()
        loo_result = leave_one_out_sensitivity(R_vals, NR_vals, R_labels, NR_labels)

        # Extract LOO p-value range
        loo_welch_ps = [r['welch_p'] for r in loo_result['loo_g1'] + loo_result['loo_g2'] if not np.isnan(r['welch_p'])]
        loo_mwu_ps = [r['mwu_p'] for r in loo_result['loo_g1'] + loo_result['loo_g2'] if not np.isnan(r['mwu_p'])]

        result = {
            'metric': metric,
            'metric_name': metric_names[metric],
            'R_mean': dual_result['g1_mean'],
            'R_std': dual_result['g1_std'],
            'R_n': dual_result['g1_n'],
            'NR_mean': dual_result['g2_mean'],
            'NR_std': dual_result['g2_std'],
            'NR_n': dual_result['g2_n'],
            'welch_p': dual_result['welch_p'],
            'mwu_p': dual_result['mwu_p'],
            'effect_size': dual_result['effect_size'],
            'fold_change': dual_result['fold_change'],
            'effect_ci_low': ci_result['effect_size_ci'][0],
            'effect_ci_high': ci_result['effect_size_ci'][1],
            'fc_ci_low': ci_result['fold_change_ci'][0],
            'fc_ci_high': ci_result['fold_change_ci'][1],
            'loo_welch_min': min(loo_welch_ps) if loo_welch_ps else np.nan,
            'loo_welch_max': max(loo_welch_ps) if loo_welch_ps else np.nan,
            'loo_mwu_min': min(loo_mwu_ps) if loo_mwu_ps else np.nan,
            'loo_mwu_max': max(loo_mwu_ps) if loo_mwu_ps else np.nan,
        }

        # Interpretation
        result['interpretation'] = interpret_significance(
            dual_result['welch_p'], dual_result['mwu_p'],
            dual_result['g1_n'], dual_result['g2_n'],
            dual_result['effect_size']
        )

        results.append(result)

        # Print summary
        print(f"    R: {result['R_mean']:.4f} ± {result['R_std']:.4f} (n={result['R_n']})")
        print(f"    NR: {result['NR_mean']:.4f} ± {result['NR_std']:.4f} (n={result['NR_n']})")
        print(f"    Welch p={result['welch_p']:.4f}, MWU p={result['mwu_p']:.4f}")
        print(f"    Effect size d={result['effect_size']:.3f} [{result['effect_ci_low']:.3f}, {result['effect_ci_high']:.3f}]")
        print(f"    LOO stability: Welch p=[{result['loo_welch_min']:.3f}, {result['loo_welch_max']:.3f}]")

    return pd.DataFrame(results)


def generate_deep_dive_figure(metrics_df, stats_df):
    """Generate comprehensive deep-dive visualization."""
    print("\nGenerating deep-dive figure...")

    fig = plt.figure(figsize=(16, 12))

    # Create grid
    gs = fig.add_gridspec(3, 4, hspace=0.35, wspace=0.3)

    metrics = ['betweenness', 'betti1_auc', 'acinar_prop', 'spatial_entropy']
    titles = ['Betweenness Centrality', 'Betti-1 AUC\n(Topological Loops)',
              'Acinar Cell Proportion', 'Spatial Entropy\n(Cell Type Diversity)']

    colors = {'R': '#2ecc71', 'NR': '#e74c3c'}

    for idx, (metric, title) in enumerate(zip(metrics, titles)):
        row = idx // 2
        col = (idx % 2) * 2

        # Get stats for this metric
        stat_row = stats_df[stats_df['metric'] == metric]
        if len(stat_row) == 0:
            continue
        stat_row = stat_row.iloc[0]

        # Main comparison plot
        ax1 = fig.add_subplot(gs[row, col])

        R_vals = metrics_df[metrics_df['response'] == 'R'][metric].dropna()
        NR_vals = metrics_df[metrics_df['response'] == 'NR'][metric].dropna()

        # Box + strip plot
        data_for_plot = []
        for v in R_vals:
            data_for_plot.append({'Response': 'R', 'value': v})
        for v in NR_vals:
            data_for_plot.append({'Response': 'NR', 'value': v})
        plot_df = pd.DataFrame(data_for_plot)

        sns.boxplot(data=plot_df, x='Response', y='value', palette=colors, ax=ax1, width=0.5)
        sns.stripplot(data=plot_df, x='Response', y='value', color='black', size=8, ax=ax1, alpha=0.7)

        ax1.set_title(title, fontsize=11, fontweight='bold')
        ax1.set_xlabel('')
        ax1.set_ylabel(metric.replace('_', ' ').title())

        # Add stats annotation
        welch_p = stat_row['welch_p']
        mwu_p = stat_row['mwu_p']
        d = stat_row['effect_size']

        sig_marker = ''
        if min(welch_p, mwu_p) < 0.05:
            sig_marker = '*'
        elif min(welch_p, mwu_p) < 0.1:
            sig_marker = '†'

        stats_text = f"Welch p={welch_p:.3f}\nMWU p={mwu_p:.3f}\nd={d:.2f} {sig_marker}"
        ax1.text(0.98, 0.98, stats_text, transform=ax1.transAxes,
                fontsize=9, va='top', ha='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        # Effect size with CI plot
        ax2 = fig.add_subplot(gs[row, col+1])

        d_low = stat_row['effect_ci_low']
        d_high = stat_row['effect_ci_high']

        ax2.errorbar([0], [d], yerr=[[d - d_low], [d_high - d]],
                    fmt='o', markersize=12, color='steelblue', capsize=10, capthick=2)
        ax2.axhline(0, color='gray', linestyle='--', alpha=0.5)
        ax2.axhline(0.8, color='green', linestyle=':', alpha=0.3, label='Large effect')
        ax2.axhline(-0.8, color='green', linestyle=':', alpha=0.3)
        ax2.axhline(0.5, color='orange', linestyle=':', alpha=0.3, label='Medium effect')
        ax2.axhline(-0.5, color='orange', linestyle=':', alpha=0.3)

        ax2.set_xlim(-0.5, 0.5)
        ax2.set_xticks([])
        ax2.set_ylabel("Cohen's d (95% CI)")
        ax2.set_title("Effect Size", fontsize=10)

        # LOO range annotation
        loo_text = f"LOO range:\np=[{stat_row['loo_welch_min']:.2f}, {stat_row['loo_welch_max']:.2f}]"
        ax2.text(0.5, 0.02, loo_text, transform=ax2.transAxes,
                fontsize=8, ha='center', va='bottom',
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))

    # Summary panel
    ax_summary = fig.add_subplot(gs[2, :])
    ax_summary.axis('off')

    # Create summary table
    summary_text = "DEEP-DIVE STATISTICAL SUMMARY\n" + "="*60 + "\n\n"
    summary_text += f"{'Metric':<25} {'Welch p':<12} {'MWU p':<12} {'d [95% CI]':<25} {'LOO Stable?':<12}\n"
    summary_text += "-"*86 + "\n"

    for _, row in stats_df.iterrows():
        ci_str = f"{row['effect_size']:.2f} [{row['effect_ci_low']:.2f}, {row['effect_ci_high']:.2f}]"
        loo_stable = "Yes" if abs(row['loo_welch_max'] - row['loo_welch_min']) < 0.2 else "Variable"
        summary_text += f"{row['metric_name']:<25} {row['welch_p']:<12.4f} {row['mwu_p']:<12.4f} {ci_str:<25} {loo_stable:<12}\n"

    summary_text += "\n" + "-"*86 + "\n"
    summary_text += "† p < 0.1 (trending)  * p < 0.05 (significant)\n"
    summary_text += "Note: n=4 Responders, n=3 Non-Responders (limited power, effect sizes more informative)"

    ax_summary.text(0.5, 0.5, summary_text, transform=ax_summary.transAxes,
                   fontsize=10, family='monospace', ha='center', va='center',
                   bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray'))

    plt.suptitle('Deep-Dive Analysis: Bootstrap CI & LOO Sensitivity\nR vs NR Treatment Response',
                fontsize=14, fontweight='bold', y=0.98)

    # Save
    plt.savefig(FIG_DIR / "fig_deepdive_analysis.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig_deepdive_analysis.pdf", bbox_inches='tight')
    plt.close()
    print("  Saved fig_deepdive_analysis.png/pdf")


def main():
    print("="*60)
    print("DEEP-DIVE STATISTICAL ANALYSIS")
    print("="*60)

    # Extract metrics
    metrics_df = extract_all_metrics()
    print(f"\nExtracted {len(metrics_df)} samples")
    print(metrics_df.to_string())

    # Save metrics
    metrics_df.to_csv(TABLE_DIR / "deep_dive_metrics.csv", index=False)
    print(f"\nSaved metrics to {TABLE_DIR / 'deep_dive_metrics.csv'}")

    # Run analysis
    stats_df = run_deep_dive_analysis(metrics_df)

    # Save stats
    stats_df.to_csv(TABLE_DIR / "deep_dive_stats.csv", index=False)
    print(f"\nSaved statistics to {TABLE_DIR / 'deep_dive_stats.csv'}")

    # Generate figure
    generate_deep_dive_figure(metrics_df, stats_df)

    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)


if __name__ == "__main__":
    main()
