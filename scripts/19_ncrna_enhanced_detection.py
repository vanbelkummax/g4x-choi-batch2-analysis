#!/usr/bin/env python3
"""
Script 19: Enhanced ncRNA Detection
===================================
Three-pronged approach:
1. mygene.info API (NCBI biotypes)
2. Pattern matching (LINC, MIR, etc.)
3. Expanded patterns from GENCODE conventions

Expected: 200-500 ncRNAs (vs current ~60)
"""

import sys
sys.path.insert(0, '/home/user/spatial-hackathon-2026')

import mygene
import scanpy as sc
import pandas as pd
import numpy as np
import re
from pathlib import Path
from scripts.config import POLYMATHIC_DIR, PDAC_METADATA

# Derive SAMPLES and METADATA from config
SAMPLES = list(PDAC_METADATA.keys())
METADATA = PDAC_METADATA


def detect_ncrna_comprehensive(adata, verbose=True):
    """
    Three-pronged ncRNA detection:
    1. mygene.info API (NCBI biotypes)
    2. Pattern matching (LINC, MIR, etc.)
    3. Extended patterns for GENCODE conventions

    Returns:
        tuple: (set of ncRNA genes, annotated adata)
    """
    genes = adata.var_names.tolist()

    if verbose:
        print(f"Analyzing {len(genes)} genes...")

    # Method 1: mygene API query
    if verbose:
        print("Method 1: Querying mygene.info...")
    mg = mygene.MyGeneInfo()

    # Query in batches to avoid timeout
    batch_size = 1000
    mygene_ncrna = set()
    mygene_pseudogene = set()

    ncrna_types = ['ncRNA', 'lncRNA', 'snoRNA', 'snRNA', 'miRNA', 'rRNA', 'tRNA', 'scRNA']
    pseudo_types = ['pseudo', 'pseudogene']

    for i in range(0, len(genes), batch_size):
        batch = genes[i:i+batch_size]
        try:
            results = mg.querymany(batch, scopes='symbol',
                                   fields='type_of_gene,ensembl.type_of_gene',
                                   species='human', verbose=False)
            for r in results:
                if isinstance(r, dict):
                    gene_type = r.get('type_of_gene', '')
                    if gene_type in ncrna_types:
                        mygene_ncrna.add(r.get('query'))
                    elif gene_type in pseudo_types:
                        mygene_pseudogene.add(r.get('query'))
        except Exception as e:
            if verbose:
                print(f"  Warning: mygene batch {i//batch_size} failed: {e}")

    if verbose:
        print(f"  mygene ncRNA: {len(mygene_ncrna)}")
        print(f"  mygene pseudogenes: {len(mygene_pseudogene)}")

    # Method 2: Pattern matching (expanded)
    if verbose:
        print("Method 2: Pattern matching...")

    patterns = [
        # Long non-coding RNAs
        (r'^LINC\d+', 'lincRNA'),           # Long intergenic ncRNA
        (r'^LINC-', 'lincRNA'),
        (r'^LOC\d{6,}', 'lincRNA'),          # Uncharacterized (many are lncRNA)

        # microRNAs
        (r'^MIR\d+', 'miRNA'),               # microRNA
        (r'^MIRLET\d+', 'miRNA'),
        (r'^MIR\d+HG', 'miRNA_host'),        # miRNA host genes

        # Small nucleolar RNAs
        (r'^SNORD\d+', 'snoRNA'),            # C/D box snoRNA
        (r'^SNORA\d+', 'snoRNA'),            # H/ACA box snoRNA
        (r'^SNHG\d+', 'snoRNA_host'),        # snoRNA host genes

        # Antisense transcripts
        (r'-AS\d*$', 'antisense'),           # Antisense
        (r'-IT\d*$', 'intronic'),            # Intronic transcript
        (r'-OT\d*$', 'overlapping'),         # Overlapping transcript

        # GENCODE lncRNA (clone-based)
        (r'^AC\d{6}\.\d+', 'gencode_lncRNA'),
        (r'^AL\d{6}\.\d+', 'gencode_lncRNA'),
        (r'^AP\d{6}\.\d+', 'gencode_lncRNA'),
        (r'^BX\d{6}\.\d+', 'gencode_lncRNA'),
        (r'^CR\d{6}\.\d+', 'gencode_lncRNA'),
        (r'^CTD-\d+', 'gencode_lncRNA'),
        (r'^CTB-\d+', 'gencode_lncRNA'),
        (r'^CTC-\d+', 'gencode_lncRNA'),
        (r'^RP\d+-\d+', 'gencode_lncRNA'),
        (r'^RP\d+\.\d+', 'gencode_lncRNA'),
        (r'^KB-\d+', 'gencode_lncRNA'),
        (r'^AC\d{6}', 'gencode_lncRNA'),
        (r'^AL\d{6}', 'gencode_lncRNA'),

        # Well-known lncRNAs
        (r'^NEAT\d+', 'known_lncRNA'),       # Nuclear enriched
        (r'^MALAT\d*', 'known_lncRNA'),      # Metastasis associated
        (r'^XIST$', 'known_lncRNA'),         # X-inactive specific
        (r'^HOTAIR$', 'known_lncRNA'),       # HOX antisense
        (r'^MEG\d+', 'known_lncRNA'),        # Maternally expressed
        (r'^H19$', 'known_lncRNA'),          # Imprinted lncRNA
        (r'^KCNQ1OT\d*', 'known_lncRNA'),
        (r'^PVT1$', 'known_lncRNA'),
        (r'^GAS5$', 'known_lncRNA'),
        (r'^DANCR$', 'known_lncRNA'),
        (r'^TUG1$', 'known_lncRNA'),
        (r'^ANRIL$', 'known_lncRNA'),

        # Other RNA types
        (r'^RN7S', 'scRNA'),                 # 7SL RNA
        (r'^RNU\d+', 'snRNA'),               # U snRNA
        (r'^RNVU\d+', 'snRNA'),
        (r'^RNA5S', 'rRNA'),                 # 5S rRNA
        (r'^RNA5-', 'rRNA'),
        (r'^MT-R', 'mt_rRNA'),               # Mitochondrial rRNA
        (r'^MT-T', 'mt_tRNA'),               # Mitochondrial tRNA

        # Vault RNAs
        (r'^VTRNA\d+', 'vault_RNA'),

        # Y RNAs
        (r'^RNY\d+', 'Y_RNA'),
    ]

    pattern_ncrna = {}  # gene -> category
    for gene in genes:
        for pattern, category in patterns:
            if re.search(pattern, gene, re.IGNORECASE):
                pattern_ncrna[gene] = category
                break

    if verbose:
        print(f"  Pattern-matched ncRNA: {len(pattern_ncrna)}")
        # Show breakdown by category
        from collections import Counter
        cat_counts = Counter(pattern_ncrna.values())
        for cat, cnt in sorted(cat_counts.items(), key=lambda x: -x[1])[:10]:
            print(f"    {cat}: {cnt}")

    # Combine all sources
    all_ncrna = mygene_ncrna | set(pattern_ncrna.keys())
    all_pseudogene = mygene_pseudogene

    # Annotate adata
    adata.var['is_ncrna'] = adata.var_names.isin(all_ncrna)
    adata.var['is_pseudogene'] = adata.var_names.isin(all_pseudogene)
    adata.var['ncrna_source'] = 'protein_coding'

    # Assign source with priority: mygene > pattern
    for gene in mygene_ncrna:
        if gene in adata.var_names:
            adata.var.loc[gene, 'ncrna_source'] = 'mygene'

    for gene in pattern_ncrna:
        if gene in adata.var_names and adata.var.loc[gene, 'ncrna_source'] == 'protein_coding':
            adata.var.loc[gene, 'ncrna_source'] = f'pattern_{pattern_ncrna[gene]}'

    for gene in mygene_pseudogene:
        if gene in adata.var_names:
            adata.var.loc[gene, 'ncrna_source'] = 'pseudogene'

    # Add category column
    adata.var['ncrna_category'] = 'protein_coding'
    for gene, cat in pattern_ncrna.items():
        if gene in adata.var_names:
            adata.var.loc[gene, 'ncrna_category'] = cat

    if verbose:
        print(f"\nTotal unique ncRNAs: {len(all_ncrna)}")
        print(f"Total pseudogenes: {len(all_pseudogene)}")
        print(f"Combined non-coding: {len(all_ncrna | all_pseudogene)}")

    return all_ncrna, all_pseudogene, adata


def analyze_ncrna_expression(adata, ncrnas, group_col='response'):
    """
    Analyze ncRNA expression patterns between groups.
    """
    ncrna_list = list(ncrnas & set(adata.var_names))

    if len(ncrna_list) == 0:
        print("No ncRNAs found in data")
        return None

    # Subset to ncRNAs
    adata_ncrna = adata[:, ncrna_list].copy()

    # Calculate mean expression per gene
    mean_expr = np.array(adata_ncrna.X.mean(axis=0)).flatten()

    # Filter to expressed ncRNAs (mean > 0.1)
    expressed_mask = mean_expr > 0.1
    expressed_ncrna = [g for g, m in zip(ncrna_list, expressed_mask) if m]

    print(f"Expressed ncRNAs (mean > 0.1): {len(expressed_ncrna)}/{len(ncrna_list)}")

    # Top expressed ncRNAs
    expr_df = pd.DataFrame({
        'gene': ncrna_list,
        'mean_expression': mean_expr,
        'category': [adata.var.loc[g, 'ncrna_category'] if g in adata.var_names else 'unknown'
                     for g in ncrna_list]
    })
    expr_df = expr_df.sort_values('mean_expression', ascending=False)

    print("\nTop 20 expressed ncRNAs:")
    print(expr_df.head(20).to_string())

    return expr_df


def main():
    """Run enhanced ncRNA detection on all samples."""
    print("=" * 60)
    print("Enhanced ncRNA Detection (mygene + patterns)")
    print("=" * 60)

    # Collect all genes from all samples first
    all_genes = set()
    sample_files = list(POLYMATHIC_DIR.glob("*_polymathic.h5ad"))

    print(f"\nLoading genes from {len(sample_files)} samples...")
    for f in sample_files:
        adata = sc.read_h5ad(f)
        all_genes.update(adata.var_names.tolist())

    print(f"Total unique genes across all samples: {len(all_genes)}")

    # Create dummy adata for gene detection
    dummy_adata = sc.AnnData(
        X=np.zeros((1, len(all_genes))),
        var=pd.DataFrame(index=list(all_genes))
    )

    # Run detection
    ncrnas, pseudogenes, dummy_adata = detect_ncrna_comprehensive(dummy_adata)

    # Save results
    output_dir = Path("/home/user/spatial-hackathon-2026/outputs")
    output_dir.mkdir(exist_ok=True)

    # Create summary dataframe
    ncrna_summary = dummy_adata.var[dummy_adata.var['is_ncrna'] | dummy_adata.var['is_pseudogene']].copy()
    ncrna_summary.to_csv(output_dir / "ncrna_detection_results.csv")
    print(f"\nSaved ncRNA annotations to: {output_dir / 'ncrna_detection_results.csv'}")

    # Now annotate each sample
    print("\n" + "=" * 60)
    print("Annotating individual samples...")
    print("=" * 60)

    for f in sample_files:
        sample_name = f.stem.replace("_polymathic", "")
        print(f"\n{sample_name}:")

        adata = sc.read_h5ad(f)

        # Transfer annotations
        adata.var['is_ncrna'] = adata.var_names.isin(ncrnas)
        adata.var['is_pseudogene'] = adata.var_names.isin(pseudogenes)
        adata.var['ncrna_category'] = 'protein_coding'

        for gene in adata.var_names:
            if gene in dummy_adata.var_names:
                adata.var.loc[gene, 'ncrna_category'] = dummy_adata.var.loc[gene, 'ncrna_category']

        sample_ncrna = sum(adata.var['is_ncrna'])
        sample_pseudo = sum(adata.var['is_pseudogene'])
        print(f"  ncRNAs: {sample_ncrna}")
        print(f"  Pseudogenes: {sample_pseudo}")

        # Save updated adata
        adata.write_h5ad(f)

        # Quick expression analysis
        if sample_ncrna > 0:
            ncrna_genes = adata.var_names[adata.var['is_ncrna']].tolist()
            mean_expr = np.array(adata[:, ncrna_genes].X.mean(axis=0)).flatten()
            top_idx = np.argsort(mean_expr)[-5:]
            print(f"  Top expressed ncRNAs: {[ncrna_genes[i] for i in top_idx]}")

    # Summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total genes analyzed: {len(all_genes)}")
    print(f"ncRNAs detected: {len(ncrnas)} ({100*len(ncrnas)/len(all_genes):.1f}%)")
    print(f"Pseudogenes detected: {len(pseudogenes)} ({100*len(pseudogenes)/len(all_genes):.1f}%)")
    print(f"Protein-coding (assumed): {len(all_genes) - len(ncrnas) - len(pseudogenes)}")

    # Category breakdown
    print("\nCategory breakdown:")
    cat_counts = ncrna_summary['ncrna_category'].value_counts()
    for cat, cnt in cat_counts.items():
        print(f"  {cat}: {cnt}")

    return ncrnas, pseudogenes


if __name__ == "__main__":
    ncrnas, pseudogenes = main()
