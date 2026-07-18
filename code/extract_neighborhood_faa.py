#!/usr/bin/env python3
"""
Extract CDS from both sides of a Rgg144 hit from a Prokka GFF+FAA pair.

--max_genes N extends N genes on EACH SIDE of Rgg144 (max total = 2*N).

Usage:
    python extract_neighborhood_faa.py \
        --rgg144_fasta data/Rgg144_standard_dna.fasta \
        --gff mid/annotation/{sample}/{sample}.gff \
        --faa mid/annotation/{sample}/{sample}.faa \
        --sample_name {sample} \
        --max_genes 20 \
        --output mid/eggnog_neighborhood/{sample}/{sample}_neighborhood.faa
"""

import argparse
import re
import sys


def parse_rgg144_fasta(fasta_file, sample_name, gene_name='Rgg144'):
    """
    Find the target gene hit for sample_name in the FASTA.
    Header format: >{gene_name}_vs_{sample}|contig={contig}|...|SubjectStart={n}|SubjectEnd={m}|...
    Returns (contig, midpoint) or None if not found.
    """
    pattern = re.compile(
        rf'>{re.escape(gene_name)}_vs_([^|]+)\|contig=([^|]+)\|.*?SubjectStart=(\d+)\|SubjectEnd=(\d+)'
    )
    with open(fasta_file) as f:
        for line in f:
            if not line.startswith('>'):
                continue
            m = pattern.search(line)
            if not m:
                continue
            fasta_sample = m.group(1).strip()
            if fasta_sample == sample_name:
                contig = m.group(2).strip()
                start = int(m.group(3))
                end = int(m.group(4))
                midpoint = (start + end) // 2
                return contig, midpoint
    return None


def parse_gff_cds(gff_file, contig):
    """
    Parse CDS features from a Prokka GFF on the given contig.
    Returns list of (locus_tag, midpoint).
    """
    cds_list = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) != 9:
                continue
            if parts[0] != contig or parts[2] != 'CDS':
                continue
            start = int(parts[3])
            end = int(parts[4])
            midpoint = (start + end) // 2
            attrs = parts[8]
            locus_tag = None
            for chunk in attrs.split(';'):
                if chunk.startswith('locus_tag='):
                    locus_tag = chunk.split('=', 1)[1].strip()
                    break
                if chunk.startswith('ID='):
                    locus_tag = chunk.split('=', 1)[1].strip()
            if locus_tag:
                cds_list.append((locus_tag, midpoint))
    return cds_list


def read_faa(faa_file):
    """
    Read a FASTA protein file. Returns dict: locus_tag -> full record string.
    Prokka headers: >locus_tag description
    """
    proteins = {}
    current_tag = None
    current_lines = []
    with open(faa_file) as f:
        for line in f:
            if line.startswith('>'):
                if current_tag:
                    proteins[current_tag] = ''.join(current_lines)
                current_tag = line.split()[0][1:]  # strip '>' and take first word
                current_lines = [line]
            else:
                current_lines.append(line)
    if current_tag:
        proteins[current_tag] = ''.join(current_lines)
    return proteins


def main():
    parser = argparse.ArgumentParser(
        description='Extract <=N CDS nearest to Rgg144 from a Prokka GFF+FAA pair.'
    )
    parser.add_argument('--rgg144_fasta', required=True,
                        help='Target gene FASTA with hit coordinates')
    parser.add_argument('--gff', required=True,
                        help='Prokka GFF file for this sample')
    parser.add_argument('--faa', required=True,
                        help='Prokka FAA file for this sample')
    parser.add_argument('--sample_name', required=True,
                        help='Sample name (matches {gene_name}_vs_{sample} in FASTA header)')
    parser.add_argument('--gene_name', default='Rgg144',
                        help='Target gene name prefix in FASTA header (default: Rgg144)')
    parser.add_argument('--max_genes', type=int, default=20,
                        help='Max CDS per side of target gene; total output <= 2*max_genes (default: 20)')
    parser.add_argument('--output', required=True,
                        help='Output subset FAA file')
    args = parser.parse_args()

    # 1. Find target gene hit position for this sample
    hit = parse_rgg144_fasta(args.rgg144_fasta, args.sample_name, args.gene_name)
    if hit is None:
        print(f'[INFO] {args.sample_name}: no {args.gene_name} hit found — writing empty file.', file=sys.stderr)
        open(args.output, 'w').close()
        return

    contig, rgg144_mid = hit

    # 2. Parse CDS features from GFF on the matching contig
    cds_list = parse_gff_cds(args.gff, contig)
    if not cds_list:
        print(f'[WARN] {args.sample_name}: no CDS found on contig {contig} — writing empty file.', file=sys.stderr)
        open(args.output, 'w').close()
        return

    # 3. Split into upstream (left) and downstream (right) of Rgg144, take <=max_genes per side
    upstream   = sorted([c for c in cds_list if c[1] <= rgg144_mid], key=lambda x: -x[1])
    downstream = sorted([c for c in cds_list if c[1] >  rgg144_mid], key=lambda x:  x[1])
    selected = upstream[:args.max_genes] + downstream[:args.max_genes]
    # Write in genomic order (ascending position)
    selected.sort(key=lambda x: x[1])

    # 4. Read FAA and write only selected sequences
    proteins = read_faa(args.faa)
    written = 0
    with open(args.output, 'w') as out:
        for tag, _ in selected:
            if tag in proteins:
                out.write(proteins[tag])
                written += 1
            else:
                print(f'[WARN] {args.sample_name}: locus_tag {tag} not found in FAA.',
                      file=sys.stderr)

    print(f'[INFO] {args.sample_name}: wrote {written} CDS '
          f'({len(upstream[:args.max_genes])} upstream + {len(downstream[:args.max_genes])} downstream '
          f'of {args.gene_name} on {contig}:{rgg144_mid}) -> {args.output}',
          file=sys.stderr)


if __name__ == '__main__':
    main()
