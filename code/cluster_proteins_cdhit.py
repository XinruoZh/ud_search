#!/usr/bin/env python3
"""
Cluster protein sequences using CD-HIT at 98% identity

This script:
1. Reads downstream_genes.tsv and filters by rank/direction
2. Looks up protein sequences from FAA files
3. Clusters sequences using CD-HIT at 98% identity (-c 0.98 -n 5)
4. Outputs FASTA (representatives) and TSV with cluster metadata

Author: Generated for Rgg144 downstream analysis pipeline
"""

import argparse
import os
import sys
import subprocess
import tempfile
import gc
from collections import defaultdict
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster proteins using CD-HIT at 98% identity"
    )
    parser.add_argument(
        "--input", "-i",
        default="output/downstream_genes.tsv",
        help="Input TSV file from extract_neighbors.py"
    )
    parser.add_argument(
        "--faa_base",
        default="mid/annotation",
        help="Base directory containing strain FAA files"
    )
    parser.add_argument(
        "--max_rank", "-r",
        type=int,
        default=5,
        help="Maximum Neighbor_Rank to include (0 = query gene itself)"
    )
    parser.add_argument(
        "--direction", "-d",
        choices=["upstream", "downstream", "both"],
        default="downstream",
        help="Filter by search direction"
    )
    parser.add_argument(
        "--prefix", "-p",
        default="D",
        help="Prefix for output protein IDs"
    )
    parser.add_argument(
        "--output_dir", "-o",
        default="output",
        help="Output directory"
    )
    parser.add_argument(
        "--batch_size", "-b",
        type=int,
        default=10,
        help="Number of strains to process at once (memory management)"
    )
    parser.add_argument(
        "--identity", "-c",
        type=float,
        default=0.98,
        help="CD-HIT sequence identity threshold (default: 0.98 = 98%%)"
    )
    parser.add_argument(
        "--word_size", "-n",
        type=int,
        default=5,
        help="CD-HIT word size (default: 5)"
    )
    parser.add_argument(
        "--threads", "-T",
        type=int,
        default=4,
        help="Number of threads for CD-HIT (default: 4)"
    )
    parser.add_argument(
        "--length-diff", "-s",
        type=float,
        default=0.8,
        dest="length_diff",
        help="CD-HIT -s: min shorter/longer length ratio to cluster together (default: 0.8)"
    )
    parser.add_argument(
        "--aL",
        type=float,
        default=0.8,
        help="CD-HIT -aL: min alignment coverage of the longer sequence (default: 0.8)"
    )
    return parser.parse_args()


def load_faa(faa_path):
    """Load FAA file into dict: locus_tag -> sequence"""
    sequences = {}
    current_tag = None
    current_seq = []

    if not os.path.exists(faa_path):
        print(f"  Warning: FAA file not found: {faa_path}", file=sys.stderr)
        return sequences

    with open(faa_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_tag and current_seq:
                    sequences[current_tag] = ''.join(current_seq)
                header = line[1:].split()[0]
                current_tag = header
                current_seq = []
            else:
                current_seq.append(line)

        if current_tag and current_seq:
            sequences[current_tag] = ''.join(current_seq)

    return sequences


def find_strain_dir(strain, base_path):
    """Find the annotation directory for a strain."""
    exact = os.path.join(base_path, strain)
    if os.path.isdir(exact):
        return exact

    try:
        for d in os.listdir(base_path):
            dir_path = os.path.join(base_path, d)
            if os.path.isdir(dir_path):
                if strain in d or d in strain:
                    return dir_path
    except OSError:
        pass

    return None


def escape_tsv_value(value):
    """Escape special characters for TSV output"""
    if pd.isna(value) or value is None:
        return ""
    value = str(value)
    value = value.replace('\t', ' ').replace('\n', ' ').replace('\r', '')
    return value


def run_cdhit(input_fasta, output_fasta, identity=0.98, word_size=5, threads=4, length_diff=0.8, aL=0.8):
    """
    Run CD-HIT on input FASTA file.

    Args:
        input_fasta: Path to input FASTA
        output_fasta: Path to output FASTA (representatives)
        identity: Sequence identity threshold (0-1)
        word_size: Word size (-n parameter)
        threads: Number of threads
        length_diff: Min shorter/longer length ratio (-s)
        aL: Min alignment coverage of the longer sequence (-aL)

    Returns:
        Path to cluster file (.clstr)
    """
    clstr_file = output_fasta + ".clstr"

    cmd = [
        "cd-hit",
        "-i", input_fasta,
        "-o", output_fasta,
        "-c", str(identity),
        "-n", str(word_size),
        "-T", str(threads),
        "-M", "0",  # Unlimited memory
        "-d", "0",  # Full sequence name in output
        "-s", str(length_diff),   # Min shorter/longer length ratio
        "-aL", str(aL),           # Min alignment coverage of longer sequence
    ]

    print(f"Running CD-HIT: {' '.join(cmd)}", file=sys.stderr)

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(result.stderr, file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"CD-HIT failed: {e.stderr}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("Error: cd-hit not found. Please install CD-HIT.", file=sys.stderr)
        sys.exit(1)

    return clstr_file


def parse_clstr_file(clstr_path):
    """
    Parse CD-HIT .clstr file.

    Returns:
        dict: cluster_id -> {
            'representative': seq_id,
            'members': [seq_id, ...],
            'size': int
        }
    """
    clusters = {}
    current_cluster = None
    current_members = []
    representative = None

    with open(clstr_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>Cluster'):
                # Save previous cluster
                if current_cluster is not None:
                    clusters[current_cluster] = {
                        'representative': representative,
                        'members': current_members,
                        'size': len(current_members)
                    }
                # Start new cluster
                current_cluster = int(line.split()[1])
                current_members = []
                representative = None
            else:
                # Parse member line: "0	100aa, >seq_id... *" or "0	100aa, >seq_id... at 98.00%"
                if '>' in line:
                    # Extract sequence ID
                    start = line.index('>') + 1
                    end = line.index('...') if '...' in line else len(line)
                    seq_id = line[start:end].strip()
                    current_members.append(seq_id)

                    # Check if this is the representative (marked with *)
                    if line.endswith('*'):
                        representative = seq_id

    # Don't forget last cluster
    if current_cluster is not None:
        clusters[current_cluster] = {
            'representative': representative,
            'members': current_members,
            'size': len(current_members)
        }

    return clusters


def main():
    args = parse_args()

    # Validate input file
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Load downstream genes data
    print(f"Loading {args.input}...", file=sys.stderr)
    df = pd.read_csv(args.input, sep='\t')
    print(f"  Loaded {len(df)} records", file=sys.stderr)

    # Filter by direction
    if args.direction != "both":
        df = df[df['search_direction'] == args.direction]
        print(f"  After direction filter ({args.direction}): {len(df)} records", file=sys.stderr)

    # Filter by rank
    df = df[df['Neighbor_Rank'].abs() <= args.max_rank]
    print(f"  After rank filter (|rank| <= {args.max_rank}): {len(df)} records", file=sys.stderr)

    # Build strain -> locus_tags mapping
    strain_locus_map = defaultdict(list)
    for _, row in df.iterrows():
        strain = row['strain']
        locus_tag = row['locus_tag']
        if pd.notna(locus_tag) and locus_tag:
            strain_locus_map[strain].append({
                'locus_tag': locus_tag.strip(),
                'eggnog_id': escape_tsv_value(row.get('eggnog_id', '')),
                'pfam': escape_tsv_value(row.get('Pfam', '')),
                'gene_name': escape_tsv_value(row.get('gene_name', '')),
                'product': escape_tsv_value(row.get('product', ''))
            })

    strains = list(strain_locus_map.keys())
    print(f"  Found {len(strains)} unique strains", file=sys.stderr)

    # Collect ALL sequences (no deduplication yet)
    all_sequences = {}  # unique_id -> {sequence, strain, locus_tag, metadata}

    stats = {
        'strains_processed': 0,
        'strains_missing': 0,
        'sequences_found': 0,
        'sequences_missing': 0
    }

    batch_size = args.batch_size
    total_batches = (len(strains) + batch_size - 1) // batch_size

    for batch_idx in range(total_batches):
        batch_start = batch_idx * batch_size
        batch_end = min(batch_start + batch_size, len(strains))
        batch_strains = strains[batch_start:batch_end]

        print(f"Processing batch {batch_idx + 1}/{total_batches} ({len(batch_strains)} strains)...", file=sys.stderr)

        for strain in batch_strains:
            strain_dir = find_strain_dir(strain, args.faa_base)

            if strain_dir is None:
                print(f"  Warning: No directory found for strain: {strain}", file=sys.stderr)
                stats['strains_missing'] += 1
                continue

            faa_path = os.path.join(strain_dir, f"{os.path.basename(strain_dir)}.faa")
            if not os.path.exists(faa_path):
                faa_path = os.path.join(strain_dir, f"{strain}.faa")

            faa_dict = load_faa(faa_path)

            if not faa_dict:
                stats['strains_missing'] += 1
                continue

            stats['strains_processed'] += 1

            for entry in strain_locus_map[strain]:
                locus_tag = entry['locus_tag']
                seq = faa_dict.get(locus_tag)

                if seq is None:
                    for key in faa_dict:
                        if key.startswith(locus_tag) or locus_tag.startswith(key):
                            seq = faa_dict[key]
                            break

                if seq is None:
                    stats['sequences_missing'] += 1
                    continue

                stats['sequences_found'] += 1

                # Create unique ID for this sequence occurrence
                unique_id = f"{strain}_{locus_tag}"
                all_sequences[unique_id] = {
                    'sequence': seq,
                    'strain': strain,
                    'locus_tag': locus_tag,
                    'eggnog_id': entry['eggnog_id'],
                    'pfam': entry['pfam'],
                    'gene_name': entry['gene_name'],
                    'product': entry['product']
                }

            del faa_dict

        gc.collect()

    print(f"\nStatistics:", file=sys.stderr)
    print(f"  Strains processed: {stats['strains_processed']}", file=sys.stderr)
    print(f"  Strains missing: {stats['strains_missing']}", file=sys.stderr)
    print(f"  Sequences found: {stats['sequences_found']}", file=sys.stderr)
    print(f"  Sequences missing: {stats['sequences_missing']}", file=sys.stderr)
    print(f"  Total sequences before clustering: {len(all_sequences)}", file=sys.stderr)

    # Write temporary FASTA for CD-HIT
    print(f"\nWriting sequences for CD-HIT...", file=sys.stderr)

    temp_fasta = os.path.join(args.output_dir, "temp_all_proteins.fasta")
    with open(temp_fasta, 'w') as f:
        for seq_id, data in all_sequences.items():
            f.write(f">{seq_id}\n")
            seq = data['sequence']
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')

    # Run CD-HIT
    print(f"\nRunning CD-HIT clustering at {args.identity*100:.0f}% identity...", file=sys.stderr)

    cdhit_output = os.path.join(args.output_dir, "cdhit_clusters")
    clstr_file = run_cdhit(temp_fasta, cdhit_output, args.identity, args.word_size, args.threads,
                           args.length_diff, args.aL)

    # Parse cluster results
    print(f"\nParsing cluster results...", file=sys.stderr)
    clusters = parse_clstr_file(clstr_file)
    print(f"  Found {len(clusters)} clusters", file=sys.stderr)

    # Build cluster metadata
    cluster_data = []
    for cluster_id, cluster_info in clusters.items():
        rep_id = cluster_info['representative']
        members = cluster_info['members']

        # Collect metadata from all members
        origins = []
        eggnog_ids = set()
        pfams = set()
        gene_names = set()
        products = set()

        for member_id in members:
            if member_id in all_sequences:
                data = all_sequences[member_id]
                origins.append(f"{data['strain']}:{data['locus_tag']}")
                if data['eggnog_id']:
                    eggnog_ids.add(data['eggnog_id'])
                if data['pfam']:
                    for p in data['pfam'].split(','):
                        if p.strip():
                            pfams.add(p.strip())
                if data['gene_name']:
                    gene_names.add(data['gene_name'])
                if data['product']:
                    products.add(data['product'])

        # Get representative sequence
        rep_seq = all_sequences[rep_id]['sequence'] if rep_id in all_sequences else ''

        cluster_data.append({
            'cluster_id': cluster_id,
            'representative': rep_id,
            'sequence': rep_seq,
            'size': len(members),
            'origins': origins,
            'eggnog_ids': eggnog_ids,
            'pfams': pfams,
            'gene_names': gene_names,
            'products': products
        })

    # Sort by cluster size (most conserved first)
    cluster_data.sort(key=lambda x: x['size'], reverse=True)

    # Write output files
    print(f"\nWriting output files...", file=sys.stderr)

    fasta_path = os.path.join(args.output_dir, "clustered_proteins.fasta")
    tsv_path = os.path.join(args.output_dir, "clustered_proteins.tsv")

    tsv_rows = []

    with open(fasta_path, 'w') as fasta_out:
        for idx, data in enumerate(cluster_data, 1):
            protein_id = f"{args.prefix}_{idx:04d}"

            # Write FASTA
            fasta_out.write(f">{protein_id}\n")
            seq = data['sequence']
            for i in range(0, len(seq), 60):
                fasta_out.write(seq[i:i+60] + '\n')

            # Determine gene_names and origin
            if data['gene_names']:
                gene_names_str = ','.join(sorted(data['gene_names']))
                gene_name_origin = 'annotation'
            else:
                non_hypo_products = [p for p in data['products'] if 'hypothetical' not in p.lower()]
                if non_hypo_products:
                    gene_names_str = ';'.join(sorted(non_hypo_products))
                    gene_name_origin = 'anno_product'
                else:
                    gene_names_str = ''
                    gene_name_origin = ''

            tsv_rows.append({
                'protein_id': protein_id,
                'sequence_length': len(seq),
                'gene_names': gene_names_str,
                'gene_name_origin': gene_name_origin,
                'eggnog_ids': ','.join(sorted(data['eggnog_ids'])) if data['eggnog_ids'] else '',
                'pfam': ','.join(sorted(data['pfams'])) if data['pfams'] else '',
                'products': ';'.join(sorted(data['products'])) if data['products'] else '',
                'interpro': '',
                'go_terms': '',
                'cluster_size': data['size'],
                'representative': data['representative'],
                'origin_count': len(data['origins']),
                'origins': ';'.join(data['origins'])
            })

    # Write TSV
    tsv_df = pd.DataFrame(tsv_rows)
    tsv_df.to_csv(tsv_path, sep='\t', index=False)

    # Clean up temporary files
    if os.path.exists(temp_fasta):
        os.remove(temp_fasta)

    print(f"\nOutput files:", file=sys.stderr)
    print(f"  FASTA: {fasta_path}", file=sys.stderr)
    print(f"  TSV:   {tsv_path}", file=sys.stderr)
    print(f"  CD-HIT clusters: {clstr_file}", file=sys.stderr)
    print(f"\nClustering: {len(all_sequences)} sequences -> {len(clusters)} clusters ({args.identity*100:.0f}% identity)", file=sys.stderr)
    print(f"\nDone!", file=sys.stderr)


if __name__ == "__main__":
    main()
