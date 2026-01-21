#!/usr/bin/env python3
"""
Figure 0: Sample Exclusion Rationale
=====================================

Documents the exclusion of YP04A from analysis due to low cell count after QC.
This figure provides transparency about data quality decisions.

Author: Max Van Belkum
Date: 2026-01-20
"""

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
from pathlib import Path
import seaborn as sns

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
TABLE_DIR = OUTPUT_DIR / "tables"
FIG_DIR = OUTPUT_DIR / "figures" / "showcase_v2"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Visual settings
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 10

# Sample metadata
PDAC_SAMPLES_ALL = {
    "YP03A": {"patient": "YP03", "timepoint": "Pre", "response": "NR"},
    "YP03C": {"patient": "YP03", "timepoint": "Post", "response": "NR"},
    "YP04A": {"patient": "YP04", "timepoint": "Pre", "response": "NR", "excluded": True},
    "YP04C": {"patient": "YP04", "timepoint": "Post", "response": "NR"},
    "YP12A": {"patient": "YP12", "timepoint": "Pre", "response": "R"},
    "YP12C": {"patient": "YP12", "timepoint": "Post", "response": "R"},
    "YP15A": {"patient": "YP15", "timepoint": "Pre", "response": "R"},
    "YP15C": {"patient": "YP15", "timepoint": "Post", "response": "R"},
}

RESPONSE_COLORS = {'R': '#2ecc71', 'NR': '#e74c3c'}
EXCLUDED_COLOR = '#95a5a6'
MIN_CELL_THRESHOLD = 500


def fig0_sample_exclusion():
    """Create figure explaining sample exclusion rationale."""
    print("Generating Figure 0: Sample Exclusion Rationale")

    # Load filtering stats
    filtering_df = pd.read_csv(TABLE_DIR / "filtering_stats.csv")

    # Filter to PDAC samples only
    pdac_df = filtering_df[filtering_df['platform'] == 'Visium'].copy()
    pdac_df = pdac_df.sort_values('sample')

    # Create figure with 4 panels
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

    # Panel A: Cell count comparison (before vs after QC)
    ax_a = fig.add_subplot(gs[0, 0])

    samples = pdac_df['sample'].values
    cells_before = pdac_df['cells_before'].values
    cells_after = pdac_df['cells_after'].values

    x = np.arange(len(samples))
    width = 0.35

    # Color bars by exclusion status
    colors_after = []
    for s in samples:
        if s == 'YP04A':
            colors_after.append(EXCLUDED_COLOR)
        else:
            colors_after.append(RESPONSE_COLORS[PDAC_SAMPLES_ALL[s]['response']])

    bars1 = ax_a.bar(x - width/2, cells_before, width, label='Before QC', color='#3498db', alpha=0.7)
    bars2 = ax_a.bar(x + width/2, cells_after, width, label='After QC', color=colors_after, edgecolor='black')

    # Add threshold line
    ax_a.axhline(y=MIN_CELL_THRESHOLD, color='red', linestyle='--', linewidth=2, label=f'Min threshold ({MIN_CELL_THRESHOLD})')

    # Highlight excluded sample
    yp04a_idx = list(samples).index('YP04A')
    ax_a.annotate('EXCLUDED', xy=(yp04a_idx + width/2, cells_after[yp04a_idx]),
                  xytext=(yp04a_idx + 0.8, cells_after[yp04a_idx] + 200),
                  fontsize=10, fontweight='bold', color='red',
                  arrowprops=dict(arrowstyle='->', color='red'))

    ax_a.set_xlabel('Sample')
    ax_a.set_ylabel('Number of Cells')
    ax_a.set_title('A) Cell Count Before and After QC', fontweight='bold')
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(samples, rotation=45, ha='right')
    ax_a.legend(loc='upper right')
    ax_a.set_ylim(0, max(cells_before) * 1.15)

    # Add value labels on bars
    for i, (before, after) in enumerate(zip(cells_before, cells_after)):
        ax_a.text(i - width/2, before + 50, str(int(before)), ha='center', va='bottom', fontsize=8)
        ax_a.text(i + width/2, after + 50, str(int(after)), ha='center', va='bottom', fontsize=8,
                  color='red' if samples[i] == 'YP04A' else 'black')

    # Panel B: Retention rate comparison
    ax_b = fig.add_subplot(gs[0, 1])

    retention = pdac_df['cells_retained_pct'].values

    colors_ret = []
    for s in samples:
        if s == 'YP04A':
            colors_ret.append(EXCLUDED_COLOR)
        else:
            colors_ret.append(RESPONSE_COLORS[PDAC_SAMPLES_ALL[s]['response']])

    bars = ax_b.bar(x, retention, color=colors_ret, edgecolor='black')

    # Add threshold indication
    ax_b.axhline(y=retention[yp04a_idx], color='red', linestyle=':', linewidth=1.5, alpha=0.7)

    ax_b.set_xlabel('Sample')
    ax_b.set_ylabel('Retention Rate (%)')
    ax_b.set_title('B) Cell Retention Rate After QC', fontweight='bold')
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(samples, rotation=45, ha='right')
    ax_b.set_ylim(0, 100)

    # Add value labels
    for i, ret in enumerate(retention):
        color = 'red' if samples[i] == 'YP04A' else 'black'
        weight = 'bold' if samples[i] == 'YP04A' else 'normal'
        ax_b.text(i, ret + 2, f'{ret:.1f}%', ha='center', va='bottom', fontsize=9,
                  color=color, fontweight=weight)

    # Panel C: Gene retention vs cell retention scatter
    ax_c = fig.add_subplot(gs[1, 0])

    gene_retention = pdac_df['genes_retained_pct'].values

    for i, s in enumerate(samples):
        color = EXCLUDED_COLOR if s == 'YP04A' else RESPONSE_COLORS[PDAC_SAMPLES_ALL[s]['response']]
        marker = 'X' if s == 'YP04A' else 'o'
        size = 200 if s == 'YP04A' else 120
        ax_c.scatter(retention[i], gene_retention[i], c=color, s=size, marker=marker,
                     edgecolors='black', linewidth=1.5, zorder=3)

        # Label points
        offset = (-15, 10) if s != 'YP04A' else (10, -20)
        ax_c.annotate(s, (retention[i], gene_retention[i]),
                      textcoords='offset points', xytext=offset, fontsize=9,
                      fontweight='bold' if s == 'YP04A' else 'normal',
                      color='red' if s == 'YP04A' else 'black')

    ax_c.set_xlabel('Cell Retention Rate (%)')
    ax_c.set_ylabel('Gene Retention Rate (%)')
    ax_c.set_title('C) QC Quality: Cell vs Gene Retention', fontweight='bold')
    ax_c.set_xlim(0, 100)
    ax_c.set_ylim(50, 105)
    ax_c.axvline(x=MIN_CELL_THRESHOLD/max(cells_before)*100, color='red', linestyle='--', alpha=0.5)

    # Add legend
    legend_elements = [
        plt.scatter([], [], c=RESPONSE_COLORS['R'], s=100, label='Responder', edgecolors='black'),
        plt.scatter([], [], c=RESPONSE_COLORS['NR'], s=100, label='Non-Responder', edgecolors='black'),
        plt.scatter([], [], c=EXCLUDED_COLOR, s=150, marker='X', label='Excluded (YP04A)', edgecolors='black'),
    ]
    ax_c.legend(handles=legend_elements, loc='lower right')

    # Panel D: Rationale text box
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.axis('off')

    # Calculate statistics
    yp04a_row = pdac_df[pdac_df['sample'] == 'YP04A'].iloc[0]
    other_retention = retention[samples != 'YP04A']

    rationale_text = f"""
    SAMPLE EXCLUSION RATIONALE
    {'='*40}

    EXCLUDED: YP04A (Patient YP04, Pre-treatment, Non-Responder)

    STATISTICS:
    - Cells before QC: {int(yp04a_row['cells_before']):,}
    - Cells after QC:  {int(yp04a_row['cells_after']):,}
    - Retention rate:  {yp04a_row['cells_retained_pct']:.1f}%
    - Gene retention:  {yp04a_row['genes_retained_pct']:.1f}%

    COMPARISON TO OTHER SAMPLES:
    - Other samples range: {int(min(cells_after[samples != 'YP04A'])):,} - {int(max(cells_after[samples != 'YP04A'])):,} cells
    - Other retention rates: {min(other_retention):.1f}% - {max(other_retention):.1f}%
    - YP04A is {min(other_retention)/yp04a_row['cells_retained_pct']:.1f}x below the next lowest

    EXCLUSION CRITERIA:
    - Minimum cells required: {MIN_CELL_THRESHOLD}
    - YP04A has only {int(yp04a_row['cells_after'])} cells ({int(yp04a_row['cells_after'])/MIN_CELL_THRESHOLD*100:.0f}% of threshold)

    DECISION:
    YP04A was excluded due to insufficient cells for reliable
    statistical analysis. With only 134 cells, spatial statistics
    and cell type proportion analyses would be underpowered and
    potentially misleading.

    NOTE: YP04C (Post-treatment from same patient) was retained
    with 1,914 cells, allowing partial representation of patient YP04.
    """

    ax_d.text(0.05, 0.95, rationale_text, transform=ax_d.transAxes,
              fontsize=10, verticalalignment='top', fontfamily='monospace',
              bbox=dict(boxstyle='round', facecolor='#f8f9fa', edgecolor='#dee2e6', alpha=0.9))

    ax_d.set_title('D) Exclusion Rationale', fontweight='bold')

    # Main title
    fig.suptitle('Figure 0: Sample Exclusion Documentation\n'
                 '(YP04A excluded due to insufficient cells after QC)',
                 fontsize=14, fontweight='bold', y=0.98)

    # Save figure
    plt.savefig(FIG_DIR / "fig0_sample_exclusion.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig0_sample_exclusion.pdf", bbox_inches='tight')
    plt.close()

    print(f"  Saved to {FIG_DIR / 'fig0_sample_exclusion.png'}")
    print(f"  YP04A: {int(yp04a_row['cells_before'])} -> {int(yp04a_row['cells_after'])} cells ({yp04a_row['cells_retained_pct']:.1f}% retention)")


if __name__ == "__main__":
    fig0_sample_exclusion()
