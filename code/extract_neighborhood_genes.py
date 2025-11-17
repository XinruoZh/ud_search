#!/usr/bin/env python3
import os
import re
import csv
import argparse
from collections import defaultdict

# Regex to parse your specific FASTA header format
# >239_P28_S15|contig=NODE_3_...|...SubjectStart=72226_SubjectEnd=73089...
HEADER_REGEX = re.compile(
    r">([^|]+)\|contig=([^|]+)\|.*?SubjectStart=(\d+)_SubjectEnd=(\d+)"
)

def parse_goi_fasta(fasta_file):
    """
    Parses the input FASTA file to get the list of genes of interest (GOI).
    Returns a list of dictionaries.
    """
    genes_of_interest = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                match = HEADER_REGEX.search(line)
                if match:
                    strain, contig, start, end = match.groups()
                    # Ensure start is always less than end
                    start_pos = min(int(start), int(end))
                    end_pos = max(int(start), int(end))
                    
                    genes_of_interest.append({
                        "strain": strain,
                        "contig": contig,
                        "start": start_pos,
                        "end": end_pos,
                        "header": line.strip()
                    })
                else:
                    print(f"Warning: Could not parse header: {line.strip()}")
                    
    print(f"Found {len(genes_of_interest)} genes of interest from FASTA.")
    return genes_of_interest

def parse_gff(gff_file):
    """
    Parses a single GFF file.
    Returns a dictionary mapping contig IDs to a list of gene features.
    """
    contig_db = defaultdict(list)
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            # We only care about gene or CDS features for location
            feature_type = parts[2]
            if feature_type not in ("gene", "CDS"):
                continue

            contig = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes_str = parts[8]
            
            # Simple attribute parsing to get a unique ID or product
            gene_id = "N/A"
            product = "N/A"
            
            if "ID=" in attributes_str:
                gene_id = attributes_str.split("ID=")[1].split(";")[0]
            if "product=" in attributes_str:
                product = attributes_str.split("product=")[1].split(";")[0]

            contig_db[contig].append({
                "type": feature_type,
                "start": start,
                "end": end,
                "strand": strand,
                "id": gene_id,
                "product": product
            })
    return contig_db

def load_all_annotations(annotation_dir):
    """
    Walks the annotation directory, finds all .gff files,
    and loads them into a master database.
    
    Assumes directory name is the strain ID (e.g., '1_103U_S92')
    """
    print(f"Scanning for .gff files in {annotation_dir}...")
    annotations_db = {}
    
    for root, dirs, files in os.walk(annotation_dir):
        for file in files:
            if file.endswith(".gff"):
                # Assumes the strain ID is the name of the parent directory
                strain_id = os.path.basename(root)
                gff_path = os.path.join(root, file)
                
                print(f"  Loading annotations for strain: {strain_id}")
                annotations_db[strain_id] = parse_gff(gff_path)
                
    print(f"Loaded annotations for {len(annotations_db)} strains.")
    return annotations_db

def find_neighbors(genes_of_interest, annotations_db, window_size):
    """
    Main function to find all neighbors for all GOIs.
    """
    results = []
    
    for goi in genes_of_interest:
        strain = goi["strain"]
        contig = goi["contig"]
        goi_start = goi["start"]
        goi_end = goi["end"]
        
        # Check if we have annotations for this strain
        if strain not in annotations_db:
            print(f"Warning: No annotation found for strain {strain}. Skipping.")
            continue
            
        strain_annotations = annotations_db[strain]
        
        # Check if the contig exists in this strain's annotations
        if contig not in strain_annotations:
            print(f"Warning: Contig {contig} not found in annotations for strain {strain}. Skipping.")
            continue
            
        # Define the search window
        goi_center = (goi_start + goi_end) / 2
        window_start = goi_center - window_size
        window_end = goi_center + window_size
        
        # Get all genes on this contig
        all_genes_on_contig = strain_annotations[contig]
        
        for gene in all_genes_on_contig:
            gene_start = gene["start"]
            gene_end = gene["end"]
            
            # Skip the GOI itself
            if gene_start == goi_start and gene_end == goi_end:
                continue
                
            # Check if the gene is within the window
            # We check if the gene's center falls within the window
            gene_center = (gene_start + gene_end) / 2
            
            if window_start <= gene_center <= window_end:
                # This is a neighbor!
                results.append({
                    "goi_strain": strain,
                    "goi_contig": contig,
                    "goi_start": goi_start,
                    "goi_end": goi_end,
                    "goi_header": goi["header"],
                    "neighbor_id": gene["id"],
                    "neighbor_type": gene["type"],
                    "neighbor_start": gene_start,
                    "neighbor_end": gene_end,
                    "neighbor_strand": gene["strand"],
                    "neighbor_product": gene["product"]
                })
                
    return results

def write_output_csv(results, output_file):
    """
    Writes the final list of neighbors to a CSV file.
    """
    if not results:
        print("No neighbors found. Output file will be empty.")
        return

    # Define CSV headers
    headers = [
        "goi_strain", "goi_contig", "goi_start", "goi_end", 
        "neighbor_id", "neighbor_type", "neighbor_start", "neighbor_end", 
        "neighbor_strand", "neighbor_product", "goi_header"
    ]
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(results)
        
    print(f"Successfully wrote {len(results)} neighbors to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Find gene neighbors from a FASTA list and GFF annotations."
    )
    parser.add_argument(
        "-f", "--fasta", 
        required=True, 
        help="Input FASTA file with genes of interest."
    )
    parser.add_argument(
        "-a", "--annotations", 
        required=True, 
        help="Main directory containing strain annotation subdirectories (e.g., 'mid/annotation')."
    )
    parser.add_argument(
        "-w", "--window", 
        required=True, 
        type=int, 
        help="Window size (length) to search upstream and downstream (e.g., 10000 for 10kb)."
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Name of the output CSV file (e.g., 'gene_neighborhoods.csv')."
    )
    
    args = parser.parse_args()
    
    # 1. Parse GOIs from FASTA
    genes_of_interest = parse_goi_fasta(args.fasta)
    
    # 2. Load all GFF annotations
    annotations_db = load_all_annotations(args.annotations)
    
    # 3. Find neighbors
    neighbors = find_neighbors(genes_of_interest, annotations_db, args.window)
    
    # 4. Write output
    write_output_csv(neighbors, args.output)

if __name__ == "__main__":
    main()