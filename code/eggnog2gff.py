import argparse
import csv
import sys
import urllib.parse

def _parse_emapper_file(eggnog_file):
    """
    Low-level parser: reads an emapper.annotations file.
    Returns: dict { query_id: annotation_dict }
    query_id is whatever key the file uses (locus_tag OR uniq_NNNN).
    """
    annotations = {}

    with open(eggnog_file, 'r') as f:
        line = f.readline()
        while line.startswith('##'):
            line = f.readline()

        header = line.strip().lstrip('#').split('\t')
        reader = csv.DictReader(f, fieldnames=header, delimiter='\t')

        for row in reader:
            query_id = row['query']

            data = {
                'preferred_name': (row.get('Preferred_name') or '').strip(),
                'description': (row.get('Description') or '').strip(),
                'GOs': (row.get('GOs') or '').strip(),
                'KEGG': (row.get('KEGG_ko') or '').strip(),
                'COG': (row.get('COG_category') or '').strip(),
                'PFAMs': (row.get('PFAMs') or '').strip(),
                'EC': (row.get('EC') or '').strip(),
                'KEGG_Pathway': (row.get('KEGG_Pathway') or '').strip(),
                'CAZy': (row.get('CAZy') or '').strip(),
                'seed_ortholog': (row.get('seed_ortholog') or '').strip()
            }

            for k, v in data.items():
                if v == '-':
                    data[k] = ''

            annotations[query_id] = data

    return annotations


def parse_eggnog_annotations(eggnog_file, origins_file=None, sample_name=None):
    """
    Parses EggNOG annotation file into a dictionary keyed by locus_tag.

    Old mode (origins_file=None):
      eggnog_file is a per-strain file already keyed by locus_tag.

    New mode (origins_file + sample_name provided):
      eggnog_file is the global unique-protein file keyed by uniq_NNNN.
      origins_file (protein_origins.tsv) maps unique_id → strain/locus_tag.
      Builds { locus_tag: annotation } for the given sample_name.
    """
    if origins_file and sample_name:
        # New global mode: join via protein_origins.tsv
        uid_annotations = _parse_emapper_file(eggnog_file)
        print(f"Loaded {len(uid_annotations)} global annotations from EggNOG.", file=sys.stderr)

        annotations = {}
        with open(origins_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row['strain'] != sample_name:
                    continue
                locus_tag = row['locus_tag']
                uid = row['unique_id']
                if uid in uid_annotations:
                    annotations[locus_tag] = uid_annotations[uid]

        print(f"Resolved {len(annotations)} locus-tag annotations for sample '{sample_name}'.", file=sys.stderr)
    else:
        # Old per-strain mode: file is already keyed by locus_tag
        annotations = _parse_emapper_file(eggnog_file)
        print(f"Loaded {len(annotations)} annotations from EggNOG.", file=sys.stderr)

    return annotations

def parse_gff_attributes(attr_string):
    """
    Parses a GFF attribute string (key=value;key=value) into a dictionary.
    """
    attributes = {}
    if not attr_string or attr_string == '.':
        return attributes
        
    parts = attr_string.strip().split(';')
    for part in parts:
        if '=' in part:
            key, value = part.split('=', 1)
            attributes[key] = value
            
    return attributes

def reconstruct_gff_attributes(attributes):
    """
    Reconstructs dictionary back to GFF attribute string.
    """
    parts = []
    for k, v in attributes.items():
        parts.append(f"{k}={v}")
    return ";".join(parts)

def escape_gff_value(value):
    """
    Escapes only critical GFF3 characters (;, =, \n, \t) to preserve readability.
    """
    if not value:
        return ""
    
    trans = str.maketrans({
        ';': '%3B',
        '=': '%3D',
        '\n': '%0A',
        '\r': '%0D',
        '\t': '%09',
        ',': '%2C' # Keep comma escaped to prevent list-splitting in single values
    })
    
    return value.translate(trans)

def process_gff(gff_file, annotations, output_file):
    with open(gff_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            
            parts = line.strip().split('\t')
            if len(parts) != 9:
                fout.write(line)
                continue
            
            feature_type = parts[2]
            attributes_str = parts[8]
            
            # Only enrich CDS or Genes or mRNA
            if feature_type not in ['CDS', 'gene', 'mRNA']:
                fout.write(line)
                continue

            attributes = parse_gff_attributes(attributes_str)
            
            # Identify the Locus Tag (The link to EggNOG)
            locus_tag = attributes.get('locus_tag', attributes.get('ID', ''))
            
            # If the locus tag matches our EggNOG data, inject info
            if locus_tag in annotations:
                anno = annotations[locus_tag]
                
                # 1. Update Gene Name
                if anno['preferred_name'] and 'gene' not in attributes:
                    attributes['gene'] = anno['preferred_name']
                
                # 2. Update Product/Description
                current_product = attributes.get('product', '')
                new_desc = anno['description']
                
                if new_desc:
                    safe_desc = escape_gff_value(new_desc)
                    if "hypothetical" in current_product.lower() and new_desc:
                        attributes['product'] = safe_desc
                    else:
                        existing_note = attributes.get('Note', '')
                        if existing_note:
                            attributes['Note'] = existing_note + escape_gff_value(", EggNOG:" + new_desc)
                        else:
                            attributes['Note'] = escape_gff_value("EggNOG:" + new_desc)

                # 3. Add GO Terms
                if anno['GOs']:
                    go_terms = anno['GOs']
                    if 'Ontology_term' in attributes:
                        attributes['Ontology_term'] += ',' + go_terms
                    else:
                        attributes['Ontology_term'] = go_terms
                        
                # 4. Add Dbxref (KEGG, COG)
                refs = []
                if anno['KEGG']:
                    kegg_id = anno['KEGG'].replace('ko:', '')
                    refs.append(f"KEGG:{kegg_id}")
                if anno['COG']:
                     refs.append(f"COG:{anno['COG']}")
                
                if refs:
                    ref_str = ",".join(refs)
                    if 'Dbxref' in attributes:
                        attributes['Dbxref'] += ',' + ref_str
                    else:
                        attributes['Dbxref'] = ref_str

                # 5. Add EC Number
                if anno['EC']:
                    attributes['EC_number'] = anno['EC']

                # 6. Add PFAMs
                if anno['PFAMs']:
                    attributes['Pfam'] = anno['PFAMs']

                # 7. Add KEGG Pathways
                if anno['KEGG_Pathway']:
                    attributes['KEGG_Pathway'] = anno['KEGG_Pathway']

                # 8. Add CAZy
                if anno['CAZy']:
                    attributes['CAZy'] = anno['CAZy']

                # 9. NEW: Add Seed Ortholog
                if anno['seed_ortholog']:
                    attributes['seed_ortholog'] = anno['seed_ortholog']

            # Reconstruct line
            parts[8] = reconstruct_gff_attributes(attributes)
            fout.write('\t'.join(parts) + '\n')
            
    print(f"Enrichment complete. Saved to {output_file}", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge EggNOG annotations into Prokka GFF (Rich Version)")
    parser.add_argument("--gff", required=True, help="Input GFF file")
    parser.add_argument("--eggnog", required=True, help="Input EggNOG annotation file")
    parser.add_argument("--output", required=True, help="Output enriched GFF file")
    parser.add_argument("--origins", default=None,
                        help="protein_origins.tsv from pool_unique_proteins.py (global mode)")
    parser.add_argument("--sample", default=None,
                        help="Strain name to filter protein_origins.tsv (required with --origins)")

    args = parser.parse_args()

    if bool(args.origins) != bool(args.sample):
        parser.error("--origins and --sample must be provided together")

    annos = parse_eggnog_annotations(args.eggnog, args.origins, args.sample)
    process_gff(args.gff, annos, args.output)