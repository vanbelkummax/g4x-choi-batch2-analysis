#!/usr/bin/env python3
"""
Create publication-ready marker table figure.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from pathlib import Path

output_dir = Path('/home/user/g4x-choi-batch2-analysis/output/figures')

# Marker data
markers = [
    # Gastric Epithelial
    ("Gastric Epithelial", "Gastric_Pit", "MUC5AC, MUC6, TFF1, TFF2", "Pit/mucous cells"),
    ("Gastric Epithelial", "Gastric_Chief", "PGC, ATP4A", "Parietal/chief cells"),
    ("Gastric Epithelial", "Goblet", "MUC2, TFF3", "Goblet cells"),
    ("Gastric Epithelial", "Intestinal_Meta", "CDX1, CDX2, MUC2, CDH17", "Intestinal metaplasia"),
    ("Gastric Epithelial", "Enteroendocrine", "GHRL, SST, GAST, CHGA", "Neuroendocrine cells"),
    ("Gastric Epithelial", "Stem_Progenitor", "LGR5, PROM1, OLFM4, SOX2, SOX9", "Stem/progenitor cells"),
    ("Gastric Epithelial", "Epithelial_General", "EPCAM, CDH1, CLDN1, CLDN3, CLDN4, CLDN7, CLDN18", "General epithelial"),
    # Immune - T Cells
    ("Immune - T Cells", "T_CD4", "CD3D, CD3E, CD4, IL7R", "Helper T cells"),
    ("Immune - T Cells", "T_CD8_Cytotoxic", "CD3D, CD3E, CD8A, GZMA, GZMB, PRF1", "Cytotoxic T cells"),
    ("Immune - T Cells", "T_Reg", "CD3D, CD4, FOXP3, IL2RA", "Regulatory T cells"),
    ("Immune - T Cells", "T_Exhausted", "PDCD1, CTLA4, LAG3, HAVCR2, TIGIT", "Exhausted T cells"),
    # Immune - B/Plasma
    ("Immune - B/Plasma", "B_Cell", "MS4A1, CD19, CD79A", "B cells"),
    ("Immune - B/Plasma", "Plasma", "IGHA1, IGHG1, IGHM, JCHAIN, SDC1", "Plasma cells"),
    # Immune - Myeloid
    ("Immune - Myeloid", "Macrophage", "CD68, CD163, CSF1R, MRC1", "Macrophages"),
    ("Immune - Myeloid", "Monocyte", "CD14, ITGAM", "Monocytes"),
    # Stromal
    ("Stromal", "Fibroblast", "COL1A1, LUM, VCAN, FN1", "Fibroblasts"),
    ("Stromal", "CAF", "ACTA2, FAP, PDGFRA, POSTN, THY1", "Cancer-associated fibroblasts"),
    ("Stromal", "Endothelial", "PECAM1, VWF", "Endothelial cells"),
]

# Category colors
category_colors = {
    "Gastric Epithelial": "#2E86AB",  # Blue
    "Immune - T Cells": "#A23B72",     # Magenta
    "Immune - B/Plasma": "#F18F01",    # Orange
    "Immune - Myeloid": "#C73E1D",     # Red
    "Stromal": "#3B8C3B",              # Green
}

# Create figure
fig, ax = plt.subplots(figsize=(14, 10))
ax.axis('off')

# Table parameters
row_height = 0.045
header_height = 0.05
col_widths = [0.18, 0.18, 0.44, 0.20]  # Category, Cell Type, Markers, Description
start_y = 0.92
start_x = 0.02

# Draw header
header_y = start_y + header_height
headers = ["Category", "Cell Type", "Markers", "Description"]
x = start_x
for i, (header, width) in enumerate(zip(headers, col_widths)):
    rect = mpatches.FancyBboxPatch(
        (x, header_y - header_height), width - 0.005, header_height,
        boxstyle="round,pad=0.01", facecolor='#2C3E50', edgecolor='white', linewidth=1
    )
    ax.add_patch(rect)
    ax.text(x + width/2, header_y - header_height/2, header,
            ha='center', va='center', fontsize=11, fontweight='bold', color='white')
    x += width

# Draw rows
y = start_y
prev_category = None
for category, cell_type, marker_genes, description in markers:
    x = start_x

    # Category column (merge cells for same category)
    if category != prev_category:
        # Count rows in this category
        cat_count = sum(1 for m in markers if m[0] == category)
        cat_height = row_height * cat_count

        rect = mpatches.FancyBboxPatch(
            (x, y - cat_height), col_widths[0] - 0.005, cat_height,
            boxstyle="round,pad=0.01", facecolor=category_colors[category],
            edgecolor='white', linewidth=1, alpha=0.8
        )
        ax.add_patch(rect)
        ax.text(x + col_widths[0]/2, y - cat_height/2, category.replace(" - ", "\n"),
                ha='center', va='center', fontsize=9, fontweight='bold', color='white')
        prev_category = category

    x += col_widths[0]

    # Cell Type column
    rect = mpatches.FancyBboxPatch(
        (x, y - row_height), col_widths[1] - 0.005, row_height,
        boxstyle="round,pad=0.01", facecolor='#ECF0F1', edgecolor='#BDC3C7', linewidth=0.5
    )
    ax.add_patch(rect)
    ax.text(x + col_widths[1]/2, y - row_height/2, cell_type.replace("_", " "),
            ha='center', va='center', fontsize=9, fontweight='bold', color='#2C3E50')
    x += col_widths[1]

    # Markers column
    rect = mpatches.FancyBboxPatch(
        (x, y - row_height), col_widths[2] - 0.005, row_height,
        boxstyle="round,pad=0.01", facecolor='white', edgecolor='#BDC3C7', linewidth=0.5
    )
    ax.add_patch(rect)
    # Style markers - italicize gene names
    ax.text(x + 0.01, y - row_height/2, marker_genes,
            ha='left', va='center', fontsize=8, fontstyle='italic', color='#34495E')
    x += col_widths[2]

    # Description column
    rect = mpatches.FancyBboxPatch(
        (x, y - row_height), col_widths[3] - 0.005, row_height,
        boxstyle="round,pad=0.01", facecolor='#F8F9FA', edgecolor='#BDC3C7', linewidth=0.5
    )
    ax.add_patch(rect)
    ax.text(x + 0.01, y - row_height/2, description,
            ha='left', va='center', fontsize=8, color='#5D6D7E')

    y -= row_height

# Title removed per user request
# ax.text(0.5, 0.99, 'Cell Type Markers Used for Manual Annotation',
#         ha='center', va='top', fontsize=14, fontweight='bold', transform=ax.transAxes)

# Summary stats
total_markers = sum(len(m[2].split(", ")) for m in markers)
ax.text(0.5, 0.03, f'18 Cell Types | {total_markers} Total Markers | G4X 337-Gene Panel',
        ha='center', va='bottom', fontsize=10, color='#7F8C8D', transform=ax.transAxes)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

plt.tight_layout()
fig.savefig(output_dir / 'table_marker_panel.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(output_dir / 'table_marker_panel.pdf', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print(f"Saved: {output_dir / 'table_marker_panel.png'}")

# Also create a simpler version as a proper matplotlib table
fig, ax = plt.subplots(figsize=(12, 8))
ax.axis('off')

# Create DataFrame for table
df = pd.DataFrame(markers, columns=['Category', 'Cell Type', 'Markers', 'Description'])

# Create table
table = ax.table(
    cellText=df.values,
    colLabels=df.columns,
    cellLoc='left',
    loc='center',
    colColours=['#2C3E50']*4,
)

# Style the table
table.auto_set_font_size(False)
table.set_fontsize(8)
table.scale(1.2, 1.5)

# Color header
for j in range(4):
    table[(0, j)].set_text_props(color='white', fontweight='bold')
    table[(0, j)].set_facecolor('#2C3E50')

# Color category cells
for i, (category, _, _, _) in enumerate(markers, start=1):
    table[(i, 0)].set_facecolor(category_colors[category])
    table[(i, 0)].set_text_props(color='white', fontweight='bold')

ax.set_title('Cell Type Markers for G4X Annotation', fontsize=14, fontweight='bold', pad=20)

plt.tight_layout()
fig.savefig(output_dir / 'table_marker_panel_simple.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print(f"Saved: {output_dir / 'table_marker_panel_simple.png'}")
print("Done!")
