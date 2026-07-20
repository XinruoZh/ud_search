#!/usr/bin/env bash
# ============================================================
# run_annotate_prokka.sh — Neighborhood FAA extraction using central Prokka annotations
#
# Prokka is run ONCE per genome and stored centrally.
# This script reads from prokka_base_dir/{prokka_species_name}/{strain}/
# and extracts neighborhood FAA per strain for each target gene query.
#
# Run run_annotate_eggnog.sh afterwards to pool proteins and run EggNOG-mapper.
#
# Usage:
#   bash scripts/run_annotate_prokka.sh [config.yaml] [--test]
#   Default config: scripts/annotate_prokka.yaml
#   --test: process only 10 strains per set; output to *_test/ dirs
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG="${1:-${SCRIPT_DIR}/annotate_prokka.yaml}"
TEST_MODE=false
TEST_N=10
for _arg in "${@:2}"; do
    [[ "$_arg" == "--test" ]] && TEST_MODE=true
done

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

# ============================================================
# YAML parsers
# ============================================================

parse_yaml_scalar() {
    # Usage: parse_yaml_scalar KEY
    python3 - "$CONFIG" "$1" <<'PYEOF'
import sys
try:
    import yaml
    with open(sys.argv[1]) as f:
        cfg = yaml.safe_load(f)
    val = cfg.get(sys.argv[2], '')
    print('' if val is None else val)
except ImportError:
    with open(sys.argv[1]) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or ':' not in line:
                continue
            k, _, v = line.partition(':')
            if k.strip() == sys.argv[2]:
                v = v.strip()
                if ' #' in v:
                    v = v[:v.index(' #')].strip()
                print(v)
                sys.exit(0)
    sys.exit(1)
PYEOF
}

parse_genome_sets() {
    # Emits one tab-separated line per genome set:
    # genome_name\tprokka_species_name\ttarget_fasta
    python3 - "$CONFIG" <<'PYEOF'
import sys
try:
    import yaml
except ImportError:
    print("ERROR: PyYAML required. Install with: pip install pyyaml", file=sys.stderr)
    sys.exit(1)

with open(sys.argv[1]) as f:
    cfg = yaml.safe_load(f)

for entry in cfg.get('genome_sets', []):
    name    = str(entry.get('genome_name', '')).strip()
    species = str(entry.get('prokka_species_name', '')).strip()
    tfasta  = str(entry.get('target_fasta', '')).strip()
    if name and species and tfasta:
        print(f"{name}\t{species}\t{tfasta}")
PYEOF
}

# ============================================================
# Load scalar config values
# ============================================================
TARGET_GENE=$(parse_yaml_scalar target_gene)
OUTPUT_DIR=$(parse_yaml_scalar output_dir)
MAX_GENES=$(parse_yaml_scalar max_genes)
PROKKA_BASE_DIR=$(parse_yaml_scalar prokka_base_dir)
CONDA_ENV=$(parse_yaml_scalar conda_env)
PARALLEL_JOBS=$(parse_yaml_scalar parallel_jobs)

EXTRACT_SCRIPT="${BASE_DIR}/code/extract_neighborhood_faa.py"

# ============================================================
# Activate conda
# ============================================================
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"
set -u

# ============================================================
# Per-strain function (exported for GNU parallel)
# ============================================================

run_extract_neighborhood() {
    local name=$1
    local anno_dir=$2
    local subset_dir=$3
    local target_fasta=$4
    local target_gene=$5
    local extract_script=$6
    local max_genes=$7
    local log=$8

    local faa="$anno_dir/$name/$name.faa"
    local gff="$anno_dir/$name/$name.gff"
    local out_dir="$subset_dir/$name"
    local subset_faa="$out_dir/${name}_neighborhood.faa"

    if [[ -f "$subset_faa" ]]; then
        echo "[SKIP] Extract: $name" | tee -a "$log"
        return 0
    fi

    if [[ ! -f "$faa" || ! -f "$gff" ]]; then
        echo "[ERROR] Extract: FAA or GFF not found for $name" | tee -a "$log"
        return 1
    fi

    mkdir -p "$out_dir"
    python3 "$extract_script" \
        --rgg144_fasta "$target_fasta" \
        --gff          "$gff" \
        --faa          "$faa" \
        --sample_name  "$name" \
        --gene_name    "$target_gene" \
        --max_genes    "$max_genes" \
        --output       "$subset_faa" \
        >> "$log" 2>&1

    if [[ -s "$subset_faa" ]]; then
        echo "[DONE]  Extract: $name" | tee -a "$log"
    else
        echo "[INFO]  Extract: $name — no hit, empty file written" | tee -a "$log"
    fi
}

export -f run_extract_neighborhood

# ============================================================
# Main loop over genome sets
# ============================================================
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")

echo "============================================"
echo "Neighborhood FAA extraction pipeline"
echo "Target gene:     $TARGET_GENE"
echo "Output dir:      $OUTPUT_DIR"
echo "Prokka base dir: $PROKKA_BASE_DIR"
echo "Config:          $CONFIG"
echo "Started:         $(date)"
[[ "$TEST_MODE" == true ]] && echo "*** TEST MODE: max $TEST_N strains per set ***"
echo "============================================"

while IFS=$'\t' read -r GENOME_NAME PROKKA_SPECIES TARGET_FASTA; do

    ANNO_DIR="${PROKKA_BASE_DIR}/${PROKKA_SPECIES}"

    if [[ ! -d "$ANNO_DIR" ]]; then
        echo "[ERROR] Central Prokka dir not found: $ANNO_DIR" >&2
        continue
    fi

    if [[ "$TEST_MODE" == true ]]; then
        SET_DIR="${OUTPUT_DIR}/${GENOME_NAME}_${TARGET_GENE}_test"
    else
        SET_DIR="${OUTPUT_DIR}/${GENOME_NAME}_${TARGET_GENE}"
    fi
    SUBSET_DIR="${SET_DIR}/eggnog_subset_faa"
    LOG_DIR="${SET_DIR}/log"

    mkdir -p "$SUBSET_DIR" "$LOG_DIR"
    LOG="${LOG_DIR}/${timestamp}_annotate_prokka.log"

    n_strains=$(find "$ANNO_DIR" -maxdepth 1 -mindepth 1 -type d | wc -l)

    echo "" | tee -a "$LOG"
    echo "============================================" | tee -a "$LOG"
    if [[ "$TEST_MODE" == true ]]; then
        echo "Genome set: $GENOME_NAME  (TEST: $TEST_N of $n_strains strains)" | tee -a "$LOG"
    else
        echo "Genome set: $GENOME_NAME  ($n_strains strains)" | tee -a "$LOG"
    fi
    echo "Prokka dir:   $ANNO_DIR" | tee -a "$LOG"
    echo "Target FASTA: $TARGET_FASTA" | tee -a "$LOG"
    echo "Output:       $SET_DIR" | tee -a "$LOG"
    echo "============================================" | tee -a "$LOG"

    # ----------------------------------------------------------
    # Extract neighborhood FAA (parallel)
    # ----------------------------------------------------------
    echo "" | tee -a "$LOG"
    echo "--- Extract neighborhood FAA (parallel -j ${PARALLEL_JOBS}) ---" | tee -a "$LOG"

    if [[ "$TEST_MODE" == true ]]; then
        find "$ANNO_DIR" -maxdepth 1 -mindepth 1 -type d -exec basename {} \; | head -"$TEST_N"
    else
        find "$ANNO_DIR" -maxdepth 1 -mindepth 1 -type d -exec basename {} \;
    fi | parallel -j "$PARALLEL_JOBS" --line-buffer \
        run_extract_neighborhood {} "$ANNO_DIR" "$SUBSET_DIR" \
        "$TARGET_FASTA" "$TARGET_GENE" "$EXTRACT_SCRIPT" "$MAX_GENES" "$LOG" || true

    n_subset=$(find "$SUBSET_DIR" -name "*_neighborhood.faa" | wc -l)
    echo "Complete: $n_subset neighborhood FAA files written" | tee -a "$LOG"

    echo "" | tee -a "$LOG"
    echo "============================================" | tee -a "$LOG"
    echo "Genome set $GENOME_NAME complete: $(date)" | tee -a "$LOG"
    echo "Key outputs:" | tee -a "$LOG"
    echo "  Neighborhood FAA: $SUBSET_DIR/" | tee -a "$LOG"
    echo "  Log:              $LOG" | tee -a "$LOG"
    echo "Next: run run_annotate_eggnog.sh to pool proteins and run EggNOG-mapper." | tee -a "$LOG"
    echo "============================================" | tee -a "$LOG"

done < <(parse_genome_sets)

echo ""
echo "============================================"
echo "All genome sets complete: $(date)"
echo "============================================"
