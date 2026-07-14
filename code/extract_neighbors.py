import argparse
import re
import os
import csv
import sys
import urllib.parse

def parse_fasta_headers(fasta_file):
    """
    Parses FASTA headers to extract Sample, Contig, Coordinates, and Strand.
    """
    targets = []
    # Regex designed for BLAST/Diamond style headers
    # Accepts SubjectStart/SubjectEnd separated by either '_' or '|'
    header_pattern = re.compile(r'>([^|]+)\|contig=([^|]+)\|.*SubjectStart=(\d+)[_|]SubjectEnd=(\d+)')

    print(f"Reading FASTA targets from {fasta_file}...", file=sys.stderr)

    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                match = header_pattern.search(line)
                if match:
                    sample_id = match.group(1).strip()
                    # Strip query_vs_ prefix if present (e.g. "Rgg144_vs_SampleID" -> "SampleID")
                    if '_vs_' in sample_id:
                        sample_id = sample_id.split('_vs_', 1)[1]
                    contig_id = match.group(2).strip()
                    raw_start = int(match.group(3))
                    raw_end = int(match.group(4))
                    
                    # Determine Strand
                    if raw_start < raw_end:
                        strand = '+'
                        p_start = raw_start
                        p_end = raw_end
                    else:
                        strand = '-'
                        p_start = raw_end 
                        p_end = raw_start 
                    
                    targets.append({
                        'sample': sample_id,
                        'contig': contig_id,
                        'strand': strand,
                        'p_start': p_start,
                        'p_end': p_end,
                        'original_header': line.strip()
                    })
    return targets

def parse_gff_line(line):
    """
    Parses a GFF line. Unquotes percent-encoded attributes.
    """
    parts = line.strip().split('\t')
    if len(parts) != 9: return None
    
    attr_str = parts[8]
    attributes = {}
    for chunk in attr_str.split(';'):
        if '=' in chunk:
            k, v = chunk.split('=', 1)
            attributes[k] = urllib.parse.unquote(v)
            
    return {
        'contig': parts[0],
        'type': parts[2],
        'start': int(parts[3]),
        'end': int(parts[4]),
        'strand': parts[6],
        'attributes': attributes
    }

def get_search_region(target, window_bp, direction):
    """
    Calculates search boundaries.
    Always includes the target gene itself, plus the window in the specified direction.
    """
    # Always include the target gene
    search_min = target['p_start']
    search_max = target['p_end']

    if target['strand'] == '+':
        if direction == 'upstream':
            search_min = target['p_start'] - window_bp
            # search_max stays at p_end to include target
        elif direction == 'downstream':
            # search_min stays at p_start to include target
            search_max = target['p_end'] + window_bp
        else: # both
            search_min = target['p_start'] - window_bp
            search_max = target['p_end'] + window_bp
    else:
        # Negative Strand: Upstream is Right (Higher Coords), Downstream is Left (Lower Coords)
        if direction == 'upstream':
            # search_min stays at p_start to include target
            search_max = target['p_end'] + window_bp
        elif direction == 'downstream':
            search_min = target['p_start'] - window_bp
            # search_max stays at p_end to include target
        else: # both
            search_min = target['p_start'] - window_bp
            search_max = target['p_end'] + window_bp

    return search_min, search_max

def resolve_contig(gff_path, p_start, p_end):
    """
    Fallback: find which GFF contig contains the query position by coordinate overlap.
    Used when the FASTA contig name doesn't match GFF contig names (e.g. Prokka renamed them).
    """
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) != 9: continue
            if parts[2] != 'CDS': continue
            if int(parts[3]) <= p_end and int(parts[4]) >= p_start:
                return parts[0]
    return None

def find_neighbors(target, gff_path, window_bp, direction):
    neighbors = []
    region_min, region_max = get_search_region(target, window_bp, direction)

    if not os.path.exists(gff_path):
        print(f"Warning: GFF not found for {target['sample']}", file=sys.stderr)
        return neighbors

    # Resolve the GFF contig name: try name match first, fall back to position scan
    fasta_contig = target['contig']
    def contig_matches(gff_contig):
        return (gff_contig == fasta_contig or
                gff_contig.startswith(fasta_contig + '_'))

    # Quick check: does the FASTA contig name appear anywhere in the GFF?
    gff_contig = None
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            col = line.split('\t', 1)[0]
            if contig_matches(col):
                gff_contig = col
                break

    if gff_contig is None:
        # Prokka renamed the contig — find it by position
        gff_contig = resolve_contig(gff_path, target['p_start'], target['p_end'])
        if gff_contig is None:
            return neighbors

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            if gff_contig not in line: continue

            feature = parse_gff_line(line)
            if not feature: continue
            if feature['type'] != 'CDS': continue

            if feature['contig'] == gff_contig:
                # Overlap Check
                if feature['start'] < region_max and feature['end'] > region_min:
                    
                    n_center = (feature['start'] + feature['end']) // 2
                    distance_bp = n_center - target['p_end']
                    
                    attrs = feature['attributes']
                    
                    neighbors.append({
                        'strain': target['sample'],
                        'contig': feature['contig'],
                        'rgg_strand': target['strand'],
                        'search_direction': direction,
                        # Rank will be calculated in the next step
                        'Neighbor_Rank': 0,
                        'distance_from_end': distance_bp,
                        'locus_tag': attrs.get('locus_tag', ''),
                        'gene_name': attrs.get('gene', ''),
                        'eggnog_id': attrs.get('seed_ortholog', ''),
                        'product': attrs.get('product', ''),
                        'start': feature['start'],
                        'end': feature['end'],
                        'strand': feature['strand'],
                        'EC_number': attrs.get('EC_number', ''),
                        'Pfam': attrs.get('Pfam', ''),
                        'CAZy': attrs.get('CAZy', ''),
                        'KEGG_Pathway': attrs.get('KEGG_Pathway', ''),
                        'GO_terms': attrs.get('Ontology_term', ''),
                        'Dbxref': attrs.get('Dbxref', '')
                    })
    return neighbors

def assign_ranks(hits, target_strand, query_len=None):
    """
    Robust Ranking Logic:
    1. Split into Left (negative distance) and Right (positive distance) buckets.
    2. Sort buckets by proximity to Rgg144 start (absolute distance).
    3. Assign -1/-2... or +1/+2... based on strand orientation.

    The rank-0 gene (query itself) is identified as the annotated CDS whose center
    falls within half the query hit length from the BLAST hit end.  This handles
    cases where the BLAST hit end and the annotated CDS center differ by more than
    a few bp (e.g. partial hits, annotation boundary differences).
    Falls back to the single closest gene if nothing qualifies.
    """
    if not hits: return []

    # Expected distance_from_end for rank-0: Rgg144's center is ~-query_len/2 from its end
    rgg_half_len = (query_len // 2) if query_len else 0
    expected_rank0_dist = -rgg_half_len
    center_threshold = max(10, rgg_half_len)

    # Buckets
    center_hits = [] # Distance ~ -rgg_half_len from Rgg144 end (The query itself)
    left_hits = []   # Physical Left (Distance < 0)
    right_hits = []  # Physical Right (Distance > 0)

    for h in hits:
        dist = h['distance_from_end']
        if abs(dist - expected_rank0_dist) <= center_threshold:
            center_hits.append(h)
        elif dist < 0:
            left_hits.append(h)
        else:
            right_hits.append(h)

    # If nothing landed in the center bucket, promote the single closest gene
    if not center_hits:
        all_sorted = sorted(hits, key=lambda x: abs(x['distance_from_end'] - expected_rank0_dist))
        closest = all_sorted[0]
        # Only promote if it's plausibly the query (within one gene length)
        if query_len and abs(closest['distance_from_end'] - expected_rank0_dist) <= query_len:
            center_hits.append(closest)
            if closest in left_hits:
                left_hits.remove(closest)
            elif closest in right_hits:
                right_hits.remove(closest)
            
    # --- Sorting ---
    # Left bucket: Sort Descending (Highest coord is closest to center)
    left_hits.sort(key=lambda x: x['start'], reverse=True)
    
    # Right bucket: Sort Ascending (Lowest coord is closest to center)
    right_hits.sort(key=lambda x: x['start'])
    
    # --- Rank Assignment ---
    # 0 for the Query Gene
    for h in center_hits:
        h['Neighbor_Rank'] = 0

    if target_strand == '+':
        # (+) Strand: Left is Upstream (-), Right is Downstream (+)
        for i, h in enumerate(left_hits):
            h['Neighbor_Rank'] = -1 * (i + 1)
        for i, h in enumerate(right_hits):
            h['Neighbor_Rank'] = 1 * (i + 1)
            
    else:
        # (-) Strand: Left is Downstream (+), Right is Upstream (-)
        
        # Left (Physical) -> Downstream (Biological)
        for i, h in enumerate(left_hits):
            h['Neighbor_Rank'] = 1 * (i + 1)
            
        # Right (Physical) -> Upstream (Biological)
        for i, h in enumerate(right_hits):
            h['Neighbor_Rank'] = -1 * (i + 1)
            
    # Combine back
    return center_hits + left_hits + right_hits

def main():
    parser = argparse.ArgumentParser(description="Extract genomic neighbors of Rgg144 hits")
    parser.add_argument("--fasta", required=True, help="FASTA file with Rgg144 hits")
    parser.add_argument("--gff_base", required=True, help="Base folder containing enriched GFFs")
    parser.add_argument("--window", type=int, default=5000, help="Window size in bp")
    parser.add_argument("--output", required=True, help="Output CSV file")
    parser.add_argument("--direction", choices=['upstream', 'downstream', 'both'], default='both', help="Search direction")
    
    args = parser.parse_args()
    
    targets = parse_fasta_headers(args.fasta)
    print(f"Found {len(targets)} targets. Mode: {args.direction}", file=sys.stderr)
    
    all_neighbors = []
    
    for t in targets:
        gff_path = os.path.join(args.gff_base, t['sample'], f"{t['sample']}.gff")
        
        hits = find_neighbors(t, gff_path, args.window, args.direction)
        query_len = t['p_end'] - t['p_start']
        ranked_hits = assign_ranks(hits, t['strand'], query_len=query_len)
        
        # Final Sort for output readability (Rank -5 to +5)
        ranked_hits.sort(key=lambda x: x['Neighbor_Rank'])
        
        all_neighbors.extend(ranked_hits)
        
    if all_neighbors:
        keys = all_neighbors[0].keys()
        with open(args.output, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=keys, delimiter='\t')
            writer.writeheader()
            writer.writerows(all_neighbors)
        print(f"Success! Wrote {len(all_neighbors)} neighbors to {args.output}", file=sys.stderr)
    else:
        print("No neighbors found.", file=sys.stderr)

if __name__ == "__main__":
    main()