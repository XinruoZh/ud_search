#!/usr/bin/env python3
"""
Update unique_proteins.tsv with InterProScan results

This script:
1. Reads InterProScan TSV output
2. Parses InterPro family memberships and GO terms
3. Updates the unique_proteins.tsv with new annotations

InterProScan TSV format (columns):
0: Protein accession (e.g., Drgg_0001)
1: Sequence MD5 digest
2: Sequence length
3: Analysis (e.g., Pfam, TIGRFAM, Gene3D)
4: Signature accession (e.g., PF00001)
5: Signature description
6: Start location
7: Stop location
8: Score (e.g., e-value)
9: Status (T = match)
10: Date
11: InterPro accession (e.g., IPR000001)
12: InterPro description
13: GO annotations (e.g., GO:0003674|GO:0008150)
14: Pathways (optional)

Author: Generated for Rgg144 downstream analysis pipeline
"""

import argparse
import os
import sys
from collections import defaultdict
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Update protein TSV with InterProScan results"
    )
    parser.add_argument(
        "--interpro",
        required=True,
        help="InterProScan TSV output file"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input protein TSV file (e.g., unique_proteins.tsv or clustered_proteins.tsv)"
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output TSV file (default: overwrite input)"
    )
    return parser.parse_args()


def parse_interpro_tsv(interpro_path):
    """
    Parse InterProScan TSV output

    Returns dict: protein_id -> {
        'interpro_entries': set of (IPR_accession, description),
        'go_terms': set of GO terms,
        'pfam_new': set of Pfam accessions (for verification)
    }
    """
    results = defaultdict(lambda: {
        'interpro_entries': set(),
        'go_terms': set(),
        'pfam_new': set()
    })

    if not os.path.exists(interpro_path):
        print(f"Error: InterPro file not found: {interpro_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Parsing InterPro results from {interpro_path}...", file=sys.stderr)

    with open(interpro_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue

            protein_id = parts[0]
            analysis = parts[3] if len(parts) > 3 else ''
            sig_accession = parts[4] if len(parts) > 4 else ''
            interpro_acc = parts[11] if len(parts) > 11 else ''
            interpro_desc = parts[12] if len(parts) > 12 else ''
            go_terms = parts[13] if len(parts) > 13 else ''

            # Store InterPro entries
            if interpro_acc and interpro_acc != '-':
                results[protein_id]['interpro_entries'].add(
                    (interpro_acc, interpro_desc)
                )

            # Store GO terms
            if go_terms and go_terms != '-':
                for go in go_terms.split('|'):
                    go = go.strip()
                    if go.startswith('GO:'):
                        results[protein_id]['go_terms'].add(go)

            # Track new Pfam hits for verification
            if analysis == 'Pfam' and sig_accession:
                results[protein_id]['pfam_new'].add(sig_accession)

    return results


def main():
    args = parse_args()

    # Parse InterPro results
    interpro_data = parse_interpro_tsv(args.interpro)
    print(f"  Found annotations for {len(interpro_data)} proteins", file=sys.stderr)

    # Load existing proteins TSV
    print(f"Loading {args.input}...", file=sys.stderr)
    df = pd.read_csv(args.input, sep='\t')
    print(f"  Loaded {len(df)} proteins", file=sys.stderr)

    # Update columns
    updated_count = 0

    for idx, row in df.iterrows():
        protein_id = row['protein_id']

        if protein_id in interpro_data:
            data = interpro_data[protein_id]
            updated_count += 1

            # Update interpro column
            if data['interpro_entries']:
                interpro_str = ';'.join(
                    f"{acc}:{desc}" for acc, desc in sorted(data['interpro_entries'])
                )
                df.at[idx, 'interpro'] = interpro_str

                # If gene_names is empty/NaN or includes "hypothetical", use first interpro description
                current_gene_names = row.get('gene_names', '')
                if pd.isna(current_gene_names) or not str(current_gene_names).strip() or 'hypothetical' in str(current_gene_names).lower():
                    # Extract just the description (not IPR ID) from first entry
                    first_entry = sorted(data['interpro_entries'])[0]  # (acc, desc)
                    first_desc = first_entry[1] if first_entry[1] else first_entry[0]
                    df.at[idx, 'gene_names'] = first_desc
                    df.at[idx, 'gene_name_origin'] = 'interpro'

            # Update go_terms column
            if data['go_terms']:
                existing_go = set()
                if pd.notna(row.get('go_terms')) and row['go_terms']:
                    existing_go = set(row['go_terms'].split(','))
                all_go = existing_go | data['go_terms']
                df.at[idx, 'go_terms'] = ','.join(sorted(all_go))

            # Optionally verify/supplement Pfam
            if data['pfam_new']:
                existing_pfam = set()
                if pd.notna(row.get('pfam')) and row['pfam']:
                    existing_pfam = set(row['pfam'].split(','))

                new_pfam = data['pfam_new'] - existing_pfam
                if new_pfam:
                    all_pfam = existing_pfam | data['pfam_new']
                    df.at[idx, 'pfam'] = ','.join(sorted(all_pfam))

    print(f"  Updated {updated_count} proteins with InterPro annotations", file=sys.stderr)

    # Save output
    output_path = args.output if args.output else args.input
    df.to_csv(output_path, sep='\t', index=False)
    print(f"  Saved to {output_path}", file=sys.stderr)

    # Summary statistics
    print(f"\nSummary:", file=sys.stderr)
    has_interpro = df['interpro'].notna() & (df['interpro'] != '')
    has_go = df['go_terms'].notna() & (df['go_terms'] != '')
    print(f"  Proteins with InterPro: {has_interpro.sum()}/{len(df)}", file=sys.stderr)
    print(f"  Proteins with GO terms: {has_go.sum()}/{len(df)}", file=sys.stderr)

    print("\nDone!", file=sys.stderr)


if __name__ == "__main__":
    main()
