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

# TODO actrually start from here

# -------------------------------------------------------------------
# ⚙️ USER SETTINGS: Please update these paths
# -------------------------------------------------------------------
# This section contains variables you will likely need to change
# when reusing this script for a different project.

# 1. Project base directory
BASE="/home/hiller-lab/xinruoz/ud_search"

# 2. Main database directory (where all your large database files are stored)
DATA="/home/hiller-lab/database"

# 3. Path to your local eggnog-mapper database folder
# (This is the folder containing files like 'eggnog.db', 'eggnog.taxa.db', etc.)
EGGNOG_DB="$DATA/eggnog_data/emapperdb-5.0.2"

# 4. Conda environment name to activate
CONDA_ENV="ud_search"

# 5. Full path to your main conda activation script
CONDA_PATH="/home/hiller-lab/anaconda3/etc/profile.d/conda.sh"

# 6. Number of CPUs/threads to use
CPUS=20

# -------------------------------------------------------------------
# 🗂️ DERIVED PATHS: (Generally no need to change)
# -------------------------------------------------------------------
# This section defines paths based on the settings above.
# As long as your directory structure is consistent, you shouldn't
# need to edit this part.

# Input directory (should contain all the prokka output folders like '1_103U_S92')
ANNO="$BASE/mid/annotation"
EGGNOG_OUT="$BASE/mid/eggnog"

# Then in the loop:


# Log directory (where the master log file will be saved)
LOG="$BASE/log/run"

# -------------------------------------------------------------------
# 📜 SCRIPT START
# -------------------------------------------------------------------

# --- Activate Conda Environment ---
echo "Activating conda environment: $CONDA_ENV"
#TODO source "$CONDA_PATH"
#TODO conda activate "$CONDA_ENV"

# --- Setup Logging ---
mkdir -p "$LOG"
mkdir -p "$EGGNOG_OUT"
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")
LOG_FILE="$LOG/${timestamp}_emapper_master.log"

echo "[$(date)]Starting EggNOG-mapper analysis..." | tee -a "$LOG_FILE"
echo "Using database: $EGGNOG_DB" | tee -a "$LOG_FILE"
echo "Annotation input directory: $ANNO" | tee -a "$LOG_FILE"
echo "----------------------------------------" | tee -a "$LOG_FILE"

# --- Main Loop ---
# Loop through all subdirectories in the ANNO directory (e.g., "1_103U_S92/")
# The "*/" ensures we only match directories

for sample_dir in "$ANNO"/*; do
    filename=$(basename "$sample_dir")

    # Define input/output paths
    input_faa="$sample_dir/$filename.faa"
    input_gff="$sample_dir/$filename.gff"
    output_dir="$EGGNOG_OUT/$filename"
    mkdir -p "$output_dir"

    # Check if the input .faa file actually exists
    if [ -f "$input_faa" ]; then
        echo "[$(date)] Running eggNOG-mapper for: $filename" | tee -a "$LOG_FILE"

        emapper.py \
            -i "$input_faa" \
            --itype proteins \
            -o "$filename" \
            --output_dir "$output_dir" \
            --data_dir "$EGGNOG_DB" \
            --cpu "$CPUS" \
            -m diamond \
            --tax_scope auto \
            --go_evidence non-electronic \
            --target_orthologs all \
            --report_orthologs \
            --excel >>"$LOG_FILE" 2>&1

        echo "[$(date)] Completed: $filename" | tee -a "$LOG_FILE"
        echo "----------------------------------------" | tee -a "$LOG_FILE"
    else
        echo "WARNING: $input_faa not found. Skipping $filename." | tee -a "$LOG_FILE"
    fi
done


echo "All EggNOG-mapper runs completed!" | tee -a "$LOG_FILE"