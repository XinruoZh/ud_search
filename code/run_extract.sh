#!/bin/bash

source /home/hiller-lab/anaconda3/etc/profile.d/conda.sh
conda activate ud_search

#Change
PRE="mitis_Rgg144"
BASE="/home/hiller-lab/xinruoz/ud_search"
ALLELE="/home/hiller-lab/xinruoz/ud_search/data/mitis_Rgg144_standard.fasta"
GENOME="$DATA/SMitis_Database_pubMLST_102225"

# Set paths
OUT="$BASE/output"
MID="$BASE/mid"
CODE="$BASE/code"
HMM="$DATA/hmm/eggNOG_strep_combined.hmm"
ANNO="$MID/annotation"
LOG="$BASE/log/run"
PRE="mitis_Rgg144"


# Create output directory if it doesn't exist
mkdir -p "$MID/annotation"
mkdir -p "$LOG"
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")

python "$CODE/extract_neighborhood_genes.py" \
    --fasta "$ALLELE" \
    --annotations "$ANNO" \
    --window 10000 \
    --output "${OUT}/${PRE}_neighborhood_genes.csv"