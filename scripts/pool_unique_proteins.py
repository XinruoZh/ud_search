#!/usr/bin/env python3
"""
pool_unique_proteins.py — Pool neighborhood FAAs across strains, deduplicate by
sequence, and write a mapping table that tracks every origin for each unique protein.

Inputs:
  --subset_dir   Directory containing per-strain neighborhood FAA files
                 ({subset_dir}/{strain}/{strain}_neighborhood.faa)
  --gff_base     Directory containing per-strain Prokka GFF files
                 ({gff_base}/{strain}/{strain}.gff)

Outputs:
  --output_faa      unique_proteins.faa  (one record per unique sequence, ID = uniq_NNNN)
  --output_origins  protein_origins.tsv  (one row per unique_id × strain occurrence)

protein_origins.tsv columns:
  unique_id   strain   locus_tag   contig   start   end   strand   original_header

Usage:
    python pool_unique_proteins.py \
        --subset_dir  path/to/eggnog_subset_faa \
        --gff_base    path/to/annotation \
        --output_faa  path/to/eggnog/unique_proteins.faa \
        --output_origins path/to/eggnog/protein_origins.tsv
"""

import argparse
import hashlib
import os
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# FAA parser
# ---------------------------------------------------------------------------

def read_faa(faa_path):
    """Return list of (locus_tag, header_line, sequence_str).
    locus_tag = first word of header (no '>').
    """
    records = []
    tag = header = seq_parts = None
    with open(faa_path) as f:
        for line in f:
            if line.startswith('>'):
                if tag is not None:
                    records.append((tag, header, ''.join(seq_parts)))
                header = line.rstrip('\n')
                tag = header[1:].split()[0]
                seq_parts = []
            else:
                if seq_parts is not None:
                    seq_parts.append(line.strip())
    if tag is not None:
        records.append((tag, header, ''.join(seq_parts)))
    return records


# ---------------------------------------------------------------------------
# GFF parser — build locus_tag -> (contig, start, end, strand)
# ---------------------------------------------------------------------------

def parse_gff_locations(gff_path):
    """Return dict: locus_tag -> (contig, start, end, strand)."""
    locations = {}
    if not os.path.isfile(gff_path):
        return locations
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) != 9 or parts[2] != 'CDS':
                continue
            contig = parts[0]
            start = parts[3]
            end = parts[4]
            strand = parts[6]
            attrs = parts[8]
            locus_tag = None
            for chunk in attrs.split(';'):
                chunk = chunk.strip()
                if chunk.startswith('locus_tag='):
                    locus_tag = chunk.split('=', 1)[1].strip()
                    break
                if chunk.startswith('ID='):
                    locus_tag = chunk.split('=', 1)[1].strip()
            if locus_tag:
                locations[locus_tag] = (contig, start, end, strand)
    return locations


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--subset_dir', required=True,
                        help='Dir with per-strain neighborhood FAAs: {subset_dir}/{strain}/{strain}_neighborhood.faa')
    parser.add_argument('--gff_base', required=True,
                        help='Dir with per-strain Prokka GFFs: {gff_base}/{strain}/{strain}.gff')
    parser.add_argument('--output_faa', required=True,
                        help='Output FASTA of deduplicated unique proteins')
    parser.add_argument('--output_origins', required=True,
                        help='Output TSV mapping unique_id to strain/locus/location')
    args = parser.parse_args()

    subset_dir = Path(args.subset_dir)
    gff_base   = Path(args.gff_base)

    # Collect all per-strain neighborhood FAA files
    faa_files = sorted(subset_dir.glob('*/*_neighborhood.faa'))
    if not faa_files:
        print(f'[ERROR] No *_neighborhood.faa files found under {subset_dir}', file=sys.stderr)
        sys.exit(1)

    print(f'[INFO] Found {len(faa_files)} neighborhood FAA files', file=sys.stderr)

    # Pass 1: hash every sequence, assign stable unique IDs
    # seq_hash -> unique_id
    hash_to_uid = {}
    # uid -> sequence string (for writing FAA)
    uid_to_seq = {}
    uid_counter = 0

    # Collect all origins: list of (unique_id, strain, locus_tag, header_line)
    origins = []

    for faa_path in faa_files:
        # Derive strain name from parent directory name
        strain = faa_path.parent.name
        records = read_faa(str(faa_path))

        for locus_tag, header_line, seq in records:
            if not seq:
                continue
            seq_upper = seq.upper()
            seq_hash = hashlib.md5(seq_upper.encode()).hexdigest()

            if seq_hash not in hash_to_uid:
                uid_counter += 1
                uid = f'uniq_{uid_counter:04d}'
                hash_to_uid[seq_hash] = uid
                uid_to_seq[uid] = seq_upper
            else:
                uid = hash_to_uid[seq_hash]

            origins.append((uid, strain, locus_tag, header_line))

    print(f'[INFO] {len(origins)} total protein records across all strains', file=sys.stderr)
    print(f'[INFO] {len(uid_to_seq)} unique protein sequences after deduplication', file=sys.stderr)

    # Pass 2: look up GFF locations per strain (cached per strain)
    gff_cache = {}

    def get_location(strain, locus_tag):
        if strain not in gff_cache:
            gff_path = gff_base / strain / f'{strain}.gff'
            gff_cache[strain] = parse_gff_locations(str(gff_path))
            if not gff_cache[strain]:
                print(f'[WARN] No GFF locations loaded for strain {strain} ({gff_path})',
                      file=sys.stderr)
        loc = gff_cache[strain].get(locus_tag)
        if loc is None:
            return ('', '', '', '')
        return loc

    # Write outputs
    os.makedirs(os.path.dirname(args.output_faa) or '.', exist_ok=True)

    # unique_proteins.faa — sorted by uid for reproducibility
    with open(args.output_faa, 'w') as f:
        for uid in sorted(uid_to_seq, key=lambda u: int(u.split('_')[1])):
            f.write(f'>{uid}\n')
            seq = uid_to_seq[uid]
            # wrap at 60 chars
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')
    print(f'[INFO] Wrote {len(uid_to_seq)} unique sequences -> {args.output_faa}', file=sys.stderr)

    # protein_origins.tsv
    header = ['unique_id', 'strain', 'locus_tag', 'contig', 'start', 'end', 'strand', 'original_header']
    with open(args.output_origins, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for uid, strain, locus_tag, orig_header in origins:
            contig, start, end, strand = get_location(strain, locus_tag)
            # Strip leading '>' from original_header for the TSV field
            orig_header_clean = orig_header.lstrip('>').strip()
            f.write('\t'.join([uid, strain, locus_tag, contig, start, end, strand, orig_header_clean]) + '\n')
    print(f'[INFO] Wrote {len(origins)} origin rows -> {args.output_origins}', file=sys.stderr)


if __name__ == '__main__':
    main()
