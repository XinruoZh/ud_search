#!/usr/bin/env python3
"""
Map cluster_id and interpro from clustered_proteins.tsv back to downstream_genes.tsv.

This script enriches downstream_genes.tsv with cluster information by joining
on (strain, locus_tag) extracted from the origins field in clustered_proteins.tsv.
"""

import argparse
import pandas as pd
import sys
import os


def parse_origins(origins_str, protein_id, gene_names, interpro):
    """
    Parse the origins field and return a list of (strain, locus_tag) -> info mappings.

    Origins format: "strain:locus_tag;strain:locus_tag;..."
    Example: "85_B26_S45:FOOMLCIH_00200;181_GCF_002096775_1_ASM209677v1:CFGCKAJJ_00855"
    """
    mappings = []
    if pd.isna(origins_str) or not origins_str:
        return mappings

    for origin in origins_str.split(';'):
        origin = origin.strip()
        if ':' in origin:
            parts = origin.split(':', 1)
            if len(parts) == 2:
                strain, locus_tag = parts
                mappings.append({
                    'strain': strain.strip(),
                    'locus_tag': locus_tag.strip(),
                    'cluster_id': protein_id,
                    'cluster_gene_names': gene_names if pd.notna(gene_names) else '',
                    'cluster_interpro': interpro if pd.notna(interpro) else ''
                })
    return mappings


def build_lookup_table(clusters_df):
    """
    Build a lookup table from clustered_proteins.tsv.
    Returns dict: (strain, locus_tag) -> {cluster_id, cluster_gene_names, cluster_interpro}
    """
    lookup = {}

    for _, row in clusters_df.iterrows():
        protein_id = row['protein_id']
        gene_names = row.get('gene_names', '')
        interpro = row.get('interpro', '')
        origins = row.get('origins', '')

        mappings = parse_origins(origins, protein_id, gene_names, interpro)
        for m in mappings:
            key = (m['strain'], m['locus_tag'])
            lookup[key] = {
                'cluster_id': m['cluster_id'],
                'cluster_gene_names': m['cluster_gene_names'],
                'cluster_interpro': m['cluster_interpro']
            }

    return lookup


def enrich_downstream_genes(downstream_df, lookup):
    """
    Add cluster information to downstream_genes.tsv.
    """
    cluster_ids = []
    cluster_gene_names = []
    cluster_interpros = []

    for _, row in downstream_df.iterrows():
        strain = row['strain']
        locus_tag = row['locus_tag']
        key = (strain, locus_tag)

        if key in lookup:
            info = lookup[key]
            cluster_ids.append(info['cluster_id'])
            cluster_gene_names.append(info['cluster_gene_names'])
            cluster_interpros.append(info['cluster_interpro'])
        else:
            cluster_ids.append('')
            cluster_gene_names.append('')
            cluster_interpros.append('')

    downstream_df['cluster_id'] = cluster_ids
    downstream_df['cluster_gene_names'] = cluster_gene_names
    downstream_df['cluster_interpro'] = cluster_interpros

    return downstream_df


def main():
    parser = argparse.ArgumentParser(
        description="Map cluster_id and interpro from clustered_proteins.tsv to downstream_genes.tsv"
    )
    parser.add_argument(
        "--downstream", "-d",
        required=True,
        help="Path to downstream_genes.tsv"
    )
    parser.add_argument(
        "--clusters", "-c",
        required=True,
        help="Path to clustered_proteins.tsv"
    )
    parser.add_argument(
        "--output", "-o",
        default="output/downstream_genes_enriched.tsv",
        help="Output path (default: output/downstream_genes_enriched.tsv)"
    )

    args = parser.parse_args()

    # Check input files exist
    if not os.path.exists(args.downstream):
        print(f"Error: {args.downstream} not found", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.clusters):
        print(f"Error: {args.clusters} not found", file=sys.stderr)
        sys.exit(1)

    print(f"Loading {args.clusters}...", file=sys.stderr)
    clusters_df = pd.read_csv(args.clusters, sep='\t')
    print(f"  Loaded {len(clusters_df)} cluster records", file=sys.stderr)

    print("Building lookup table from origins...", file=sys.stderr)
    lookup = build_lookup_table(clusters_df)
    print(f"  Created {len(lookup)} (strain, locus_tag) -> cluster mappings", file=sys.stderr)

    print(f"Loading {args.downstream}...", file=sys.stderr)
    downstream_df = pd.read_csv(args.downstream, sep='\t')
    print(f"  Loaded {len(downstream_df)} downstream gene records", file=sys.stderr)

    print("Enriching with cluster information...", file=sys.stderr)
    enriched_df = enrich_downstream_genes(downstream_df, lookup)

    # Count how many rows got mapped
    mapped_count = (enriched_df['cluster_id'] != '').sum()
    print(f"  Mapped {mapped_count}/{len(enriched_df)} rows ({100*mapped_count/len(enriched_df):.1f}%)", file=sys.stderr)

    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"Saving to {args.output}...", file=sys.stderr)
    enriched_df.to_csv(args.output, sep='\t', index=False)
    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
