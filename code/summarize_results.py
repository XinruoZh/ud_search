#!/usr/bin/env python3

import csv
import argparse
from collections import defaultdict

def summarize_neighbors(input_file, output_strains, output_genes):
    """
    Reads the neighbor CSV and creates two summary reports:
    1. A strain-centric report (strains as rows)
    2. A gene-centric report (genes as rows)
    """
    
    # Use defaultdict(set) to automatically handle unique items
    strain_data = defaultdict(set)
    gene_data = defaultdict(set)

    print(f"Reading and processing {input_file}...")
    
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                strain = row.get('goi_strain')
                neighbor = row.get('neighbor_id')
                
                if strain and neighbor:
                    # Add neighbor to this strain's set
                    strain_data[strain].add(neighbor)
                    # Add this strain to this neighbor's set
                    gene_data[neighbor].add(strain)

    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        return
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    if not strain_data:
        print("No valid data found in the input file.")
        return

    # --- 1. Write the Strain-centric CSV ---
    
    print(f"Writing strain summary to {output_strains}...")
    
    # Find the max number of neighbors for column padding
    max_neighbors = max(len(genes) for genes in strain_data.values())
    
    # Create headers
    strain_headers = ['Strain', 'Neighbor_Count'] + \
                     [f'Neighbor_{i+1}' for i in range(max_neighbors)]
    
    with open(output_strains, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(strain_headers)
        
        # Sort by strain name for consistent output
        for strain, neighbors_set in sorted(strain_data.items()):
            neighbor_list = sorted(list(neighbors_set))
            count = len(neighbor_list)
            
            # Create the row and pad with empty strings
            row = [strain, count] + neighbor_list
            row.extend([''] * (max_neighbors - count))
            writer.writerow(row)

    print(f"Successfully wrote {len(strain_data)} strains.")

    # --- 2. Write the Gene-centric CSV ---
    
    print(f"Writing gene summary to {output_genes}...")
    
    # Find the max number of strains for column padding
    max_strains = max(len(strains) for strains in gene_data.values())
    
    # Create headers
    gene_headers = ['Gene_Name', 'Strain_Count'] + \
                   [f'Strain_{i+1}' for i in range(max_strains)]
                   
    with open(output_genes, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(gene_headers)
        
        # Sort by gene name for consistent output
        for gene, strains_set in sorted(gene_data.items()):
            strains_list = sorted(list(strains_set))
            count = len(strains_list)
            
            # Create the row and pad with empty strings
            row = [gene, count] + strains_list
            row.extend([''] * (max_strains - count))
            writer.writerow(row)
            
    print(f"Successfully wrote {len(gene_data)} genes.")


def main():
    parser = argparse.ArgumentParser(
        description="Summarize gene neighbor CSV into two reports."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input CSV file (the one you provided)."
    )
    parser.add_argument(
        "-o1", "--output_strains",
        required=True,
        help="Output file for the strain-centric summary (e.g., 'strains_summary.csv')."
    )
    parser.add_argument(
        "-o2", "--output_genes",
        required=True,
        help="Output file for the gene-centric summary (e.g., 'genes_summary.csv')."
    )
    
    args = parser.parse_args()
    summarize_neighbors(args.input, args.output_strains, args.output_genes)

if __name__ == "__main__":
    main()