#!/bin/bash

source /home/hiller-lab/anaconda3/etc/profile.d/conda.sh
conda activate ud_search

# Set paths
BASE="/home/hiller-lab/xinruoz/ud_search"
DATA="/home/hiller-lab/database"
GENOME="$DATA/SMitis_Database_pubMLST_102225"
OUT="$BASE/output"
MID="$BASE/mid"
HMM="$DATA/hmm/eggNOG_strep_combined.hmm"
ANNO="$MID/annotation"
LOG="$BASE/log/run"

# Create output directory if it doesn't exist
mkdir -p "$MID/annotation"
mkdir -p "$LOG"
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")

# Loop through all fasta files in the genome directory
for fasta_file in "$GENOME"/*.fasta; do

    # Get the base filename without extension
    filename=$(basename "$fasta_file" .fasta)
    
    # Create output directory for this specific genome
    output_dir="$ANNO/$filename"
    
    echo "Annotating: $filename"
    
    # Run prokka
    prokka \
        --cpus 20 \
        --outdir "$output_dir" \
        --prefix "$filename" \
        --hmms "$HMM" \
        --quiet \
        "$fasta_file" >>"$LOG/${timestamp}_prokka_master.log" 2>&1
    
    echo "Completed annotation for: $filename"
    echo "----------------------------------------"
done

echo "All annotations completed!"