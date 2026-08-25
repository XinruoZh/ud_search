#!/usr/bin/env bash
# ============================================================
# run_annotate_eggnog.sh — Pool neighborhood FAAs + EggNOG-mapper
#
# Expects run_annotate_prokka.sh to have been run first.
#
# For each genome set defined in the YAML:
#   Step 4: Pool + deduplicate neighborhood proteins (Python)
#   Step 5: Run EggNOG-mapper once on unique proteins
#
# Usage:
#   bash scripts/run_annotate_eggnog.sh [config.yaml]
#   Default config: scripts/annotate_eggnog.yaml
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONFIG="${1:-${SCRIPT_DIR}/annotate_eggnog.yaml}"

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
    # Emits one tab-separated line per genome set: genome_name\tprokka_dir
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
    name      = str(entry.get('genome_name', '')).strip()
    prokka_dir = str(entry.get('prokka_dir', '')).strip()
    if name:
        print(f"{name}\t{prokka_dir}")
PYEOF
}

# ============================================================
# Load scalar config values
# ============================================================
TARGET_GENE=$(parse_yaml_scalar target_gene)
OUTPUT_DIR=$(parse_yaml_scalar output_dir)
EGGNOG_DB=$(parse_yaml_scalar eggnog_db)
CONDA_ENV=$(parse_yaml_scalar conda_env)
EGGNOG_CPUS=$(parse_yaml_scalar eggnog_cpus)

POOL_SCRIPT="${SCRIPT_DIR}/pool_unique_proteins.py"

# ============================================================
# Activate conda
# ============================================================
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"
set -u

# ============================================================
# Main loop over genome sets
# ============================================================
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")

echo "============================================"
echo "EggNOG-mapper pipeline"
echo "Target gene: $TARGET_GENE"
echo "Output dir:  $OUTPUT_DIR"
echo "Config:      $CONFIG"
echo "Started:     $(date)"
echo "============================================"

while IFS=$'\t' read -r GENOME_NAME PROKKA_DIR; do

    SET_DIR="${OUTPUT_DIR}/${GENOME_NAME}_${TARGET_GENE}"
    ANNO_DIR="$PROKKA_DIR"
    SUBSET_DIR="${SET_DIR}/eggnog_subset_faa"
    EGGNOG_DIR="${SET_DIR}/eggnog"
    LOG_DIR="${SET_DIR}/log"

    mkdir -p "$EGGNOG_DIR" "$LOG_DIR"
    LOG="${LOG_DIR}/${timestamp}_annotate_eggnog.log"

    if [[ ! -d "$SUBSET_DIR" ]]; then
        echo "[ERROR] Neighborhood FAA dir not found: $SUBSET_DIR — run run_annotate_prokka.sh first" | tee -a "$LOG" >&2
        continue
    fi

    echo "" | tee -a "$LOG"
    echo "============================================" | tee -a "$LOG"
    echo "Genome set: $GENOME_NAME" | tee -a "$LOG"
    echo "Output: $SET_DIR" | tee -a "$LOG"
    echo "============================================" | tee -a "$LOG"

    # ----------------------------------------------------------
    # Step 4: Pool + deduplicate proteins
    # ----------------------------------------------------------
    echo "" | tee -a "$LOG"
    echo "--- Step 4/5: Pool + deduplicate proteins ---" | tee -a "$LOG"

    UNIQUE_FAA="${EGGNOG_DIR}/unique_proteins.faa"
    ORIGINS_TSV="${EGGNOG_DIR}/protein_origins.tsv"

    n_faa=$(find "$SUBSET_DIR" -name "*_neighborhood.faa" -size +0 | wc -l)
    if [[ $n_faa -eq 0 ]]; then
        echo "[SKIP] No non-empty neighborhood FAA files in $SUBSET_DIR — skipping EggNOG for $GENOME_NAME" | tee -a "$LOG"
        continue
    fi

    if [[ -f "$UNIQUE_FAA" && -f "$ORIGINS_TSV" ]]; then
        echo "[SKIP] unique_proteins.faa and protein_origins.tsv already exist" | tee -a "$LOG"
    else
        python3 "$POOL_SCRIPT" \
            --subset_dir     "$SUBSET_DIR" \
            --gff_base       "$ANNO_DIR" \
            --output_faa     "$UNIQUE_FAA" \
            --output_origins "$ORIGINS_TSV" \
            2>&1 | tee -a "$LOG"
    fi

    n_unique=$(grep -c '^>' "$UNIQUE_FAA" 2>/dev/null || echo 0)
    n_origins=$(( $(wc -l < "$ORIGINS_TSV") - 1 ))
    echo "Step 4 complete: $n_unique unique proteins, $n_origins origin rows" | tee -a "$LOG"

    # ----------------------------------------------------------
    # Step 5: EggNOG-mapper on unique proteins (single run)
    # ----------------------------------------------------------
    echo "" | tee -a "$LOG"
    echo "--- Step 5/5: EggNOG-mapper on $n_unique unique proteins ---" | tee -a "$LOG"

    EMAPPER_PREFIX="${GENOME_NAME}_${TARGET_GENE}"
    EMAPPER_ANNOTATIONS="${EGGNOG_DIR}/${EMAPPER_PREFIX}.emapper.annotations"

    if [[ -f "$EMAPPER_ANNOTATIONS" ]]; then
        echo "[SKIP] EggNOG annotations already exist: $EMAPPER_ANNOTATIONS" | tee -a "$LOG"
    else
        echo "[START] EggNOG-mapper" | tee -a "$LOG"
        emapper.py \
            -i                  "$UNIQUE_FAA" \
            --itype             proteins \
            --output            "$EMAPPER_PREFIX" \
            --output_dir        "$EGGNOG_DIR" \
            --data_dir          "$EGGNOG_DB" \
            --cpu               "$EGGNOG_CPUS" \
            -m                  diamond \
            --tax_scope         auto \
            --target_orthologs  one2one \
            --go_evidence       non-electronic \
            --override \
            >> "$LOG" 2>&1

        if [[ -f "$EMAPPER_ANNOTATIONS" ]]; then
            echo "[DONE]  EggNOG-mapper: $EMAPPER_ANNOTATIONS" | tee -a "$LOG"
        else
            echo "[FAIL]  EggNOG-mapper did not produce output — check $LOG" | tee -a "$LOG"
        fi
    fi

    echo "" | tee -a "$LOG"
    echo "============================================" | tee -a "$LOG"
    echo "Genome set $GENOME_NAME complete: $(date)" | tee -a "$LOG"
    echo "Key outputs:" | tee -a "$LOG"
    echo "  Unique proteins:  $UNIQUE_FAA" | tee -a "$LOG"
    echo "  Origins table:    $ORIGINS_TSV" | tee -a "$LOG"
    echo "  EggNOG result:    $EMAPPER_ANNOTATIONS" | tee -a "$LOG"
    echo "  Log:              $LOG" | tee -a "$LOG"
    echo "============================================" | tee -a "$LOG"

done < <(parse_genome_sets)

echo ""
echo "============================================"
echo "All genome sets complete: $(date)"
echo "============================================"
