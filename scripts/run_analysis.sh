#!/usr/bin/env bash
# ============================================================
# Post-annotation analysis — runs after run_annotate_contigs.sh
# Stages 3-9: eggnog2gff → extract_neighbors → cluster_cdhit
#             → (interpro) → update_interpro → map_cluster → analyze_chains
#
# Usage: bash scripts/run_analysis.sh [config.yaml]
# Default config: scripts/analysis.yaml
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONFIG="${1:-${SCRIPT_DIR}/analysis.yaml}"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

# --- Parse YAML with Python ---
parse_yaml() {
    python3 - "$CONFIG" "$1" <<'PYEOF'
import sys, pathlib
try:
    import yaml
except ImportError:
    # Minimal fallback parser for simple key: value yaml
    config_file, key = sys.argv[1], sys.argv[2]
    with open(config_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or ':' not in line:
                continue
            k, _, v = line.partition(':')
            if k.strip() == key:
                v = v.strip()
                if ' #' in v:
                    v = v[:v.index(' #')].strip()
                print(v)
                sys.exit(0)
    sys.exit(1)

config_file, key = sys.argv[1], sys.argv[2]
with open(config_file) as f:
    cfg = yaml.safe_load(f)
val = cfg.get(key, '')
print('' if val is None else val)
PYEOF
}

# Load all parameters
BASE_DIR=$(parse_yaml base_dir)
_anno_dir=$(parse_yaml annotation_dir)
if [[ "$_anno_dir" == /* ]]; then
    ANNO_DIR="$_anno_dir"
else
    ANNO_DIR="${BASE_DIR}/${_anno_dir}"
fi
EGGNOG_ANNOTATIONS="${BASE_DIR}/$(parse_yaml eggnog_annotations)"
PROTEIN_ORIGINS="${BASE_DIR}/$(parse_yaml protein_origins)"
QUERY_FASTA="$(parse_yaml query_fasta)"
ENRICHED_DIR="${BASE_DIR}/$(parse_yaml annotation_with_eggnog_dir)"
OUTPUT_DIR="${BASE_DIR}/$(parse_yaml output_dir)"
LOG_DIR="${BASE_DIR}/$(parse_yaml log_dir)"
WINDOW=$(parse_yaml window)
DIRECTION=$(parse_yaml direction)
MAX_RANK=$(parse_yaml max_rank)
CDHIT_IDENTITY=$(parse_yaml cdhit_identity)
CDHIT_THREADS=$(parse_yaml cdhit_threads)
CLUSTER_PREFIX=$(parse_yaml cluster_prefix)
CHAIN_MAX_RANK=$(parse_yaml chain_max_rank)
CONDA_ENV=$(parse_yaml conda_env)
BLAST_NAMES=$(parse_yaml balst_protein_name 2>/dev/null || true)
END_AT_GENE=$(parse_yaml end_at_gene 2>/dev/null || true)
STOP_BEFORE_GENE=$(parse_yaml stop_before_gene 2>/dev/null || true)
MASTER_COLOR_MAP=$(parse_yaml master_color_map 2>/dev/null || true)
INTERPROSCAN_CPUS=$(parse_yaml interproscan_cpus 2>/dev/null || true)
[[ -z "$INTERPROSCAN_CPUS" ]] && INTERPROSCAN_CPUS=4

CODE_DIR="$(cd "${SCRIPT_DIR}/../code" && pwd)"

# --- Activate conda ---
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"
set -u

# --- Setup ---
mkdir -p "$ENRICHED_DIR" "$OUTPUT_DIR" "${OUTPUT_DIR}/cluster" "$LOG_DIR"
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOG="${LOG_DIR}/${TIMESTAMP}_analysis.log"

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

log "============================================"
log "Analysis pipeline — config: $CONFIG"
log "Base dir:   $BASE_DIR"
log "Output dir: $OUTPUT_DIR"
log "Direction:  $DIRECTION  Window: ${WINDOW}bp  Max_rank: $MAX_RANK"
log "============================================"

# ============================================================
# Stage 3: eggnog2gff — merge EggNOG annotations into GFF
# ============================================================
log ""
log "--- Stage 3: eggnog2gff ---"
count=0; skipped=0

if [[ ! -f "$EGGNOG_ANNOTATIONS" || ! -f "$PROTEIN_ORIGINS" ]]; then
    log "[ERROR] Missing global EggNOG files:"
    log "  annotations: $EGGNOG_ANNOTATIONS"
    log "  origins:     $PROTEIN_ORIGINS"
    log "  Run run_annotate_contigs.sh first."
    exit 1
fi

for sample_path in "$ANNO_DIR"/*/; do
    [[ -d "$sample_path" ]] || continue
    name=$(basename "$sample_path")
    gff="${sample_path}${name}.gff"
    out_dir="${ENRICHED_DIR}/${name}"
    out_gff="${out_dir}/${name}.gff"

    if [[ -f "$out_gff" ]]; then
        skipped=$(( skipped + 1 )); continue
    fi
    if [[ ! -f "$gff" ]]; then
        log "[SKIP] $name: missing Prokka GFF"
        skipped=$(( skipped + 1 )); continue
    fi

    mkdir -p "$out_dir"
    python3 "${CODE_DIR}/eggnog2gff.py" \
        --gff     "$gff" \
        --eggnog  "$EGGNOG_ANNOTATIONS" \
        --origins "$PROTEIN_ORIGINS" \
        --sample  "$name" \
        --output  "$out_gff" \
        >> "$LOG" 2>&1
    count=$(( count + 1 ))
done

log "Stage 3 complete: processed=$count skipped=$skipped"

# ============================================================
# Stage 4: extract_neighbors — find genes near Rgg144 hits
# ============================================================
log ""
log "--- Stage 4: extract_neighbors ---"
DOWNSTREAM_TSV="${OUTPUT_DIR}/downstream_genes.tsv"

if [[ -f "$DOWNSTREAM_TSV" ]]; then
    log "[SKIP] $DOWNSTREAM_TSV already exists"
else
    python3 "${CODE_DIR}/extract_neighbors.py" \
        --fasta     "$QUERY_FASTA" \
        --gff_base  "$ENRICHED_DIR" \
        --output    "$DOWNSTREAM_TSV" \
        --direction "$DIRECTION" \
        --window    "$WINDOW" \
        2>&1 | tee -a "$LOG"
fi

log "Stage 4 complete: $DOWNSTREAM_TSV"

# ============================================================
# Stage 5: cluster_proteins_cdhit — CD-HIT clustering
# ============================================================
log ""
log "--- Stage 5: cluster_proteins_cdhit ---"
CLUSTER_TSV="${OUTPUT_DIR}/cluster/clustered_proteins.tsv"
CLUSTER_FASTA="${OUTPUT_DIR}/cluster/clustered_proteins.fasta"

if [[ -f "$CLUSTER_TSV" ]]; then
    log "[SKIP] $CLUSTER_TSV already exists"
else
    python3 "${CODE_DIR}/cluster_proteins_cdhit.py" \
        --input      "$DOWNSTREAM_TSV" \
        --faa_base   "$ANNO_DIR" \
        --max_rank   "$MAX_RANK" \
        --direction  "$DIRECTION" \
        --prefix     "$CLUSTER_PREFIX" \
        --identity   "$CDHIT_IDENTITY" \
        --threads    "$CDHIT_THREADS" \
        --output_dir "${OUTPUT_DIR}/cluster" \
        2>&1 | tee -a "$LOG"
fi

log "Stage 5 complete: $CLUSTER_TSV"

# ============================================================
# Stage 6: InterProScan on clustered proteins
# ============================================================
log ""
log "--- Stage 6: interproscan ---"
INTERPRO_TSV="${OUTPUT_DIR}/interpro_results/clustered_proteins.fasta.tsv"

if [[ -f "$INTERPRO_TSV" ]]; then
    log "[SKIP] InterProScan TSV already exists: $INTERPRO_TSV"
elif [[ ! -s "$CLUSTER_FASTA" ]]; then
    log "[SKIP] Stage 6: $CLUSTER_FASTA not found or empty"
else
    mkdir -p "${OUTPUT_DIR}/interpro_results"
    log "[START] interproscan -i $CLUSTER_FASTA --cpu $INTERPROSCAN_CPUS"
    interproscan.sh \
        -i   "$CLUSTER_FASTA" \
        -f   TSV \
        -o   "$INTERPRO_TSV" \
        --cpu "$INTERPROSCAN_CPUS" \
        --goterms \
        --iprlookup \
        -dp \
        >> "$LOG" 2>&1 \
        && log "[DONE]  Stage 6: $INTERPRO_TSV" \
        || log "[WARN]  Stage 6: interproscan failed — stage 7 will be skipped"
fi

# ============================================================
# Stage 7: update_tsv_with_interpro
# ============================================================
UPDATED_CLUSTER_TSV="${OUTPUT_DIR}/cluster/clustered_proteins_updated.tsv"

if [[ -f "$INTERPRO_TSV" ]]; then
    log ""
    log "--- Stage 7: update_tsv_with_interpro ---"
    if [[ -f "$UPDATED_CLUSTER_TSV" ]]; then
        log "[SKIP] $UPDATED_CLUSTER_TSV already exists"
    else
        python3 "${CODE_DIR}/update_tsv_with_interpro.py" \
            --interpro "$INTERPRO_TSV" \
            --input    "$CLUSTER_TSV" \
            --output   "$UPDATED_CLUSTER_TSV" \
            2>&1 | tee -a "$LOG"
        log "Stage 7 complete: $UPDATED_CLUSTER_TSV"
    fi
    CLUSTER_TSV_FOR_MAP="$UPDATED_CLUSTER_TSV"
else
    log ""
    log "[INFO] No InterProScan results (stage 6 failed or was skipped) — skipping stage 7"
    CLUSTER_TSV_FOR_MAP="$CLUSTER_TSV"
fi

# ============================================================
# Stage 8: map_cluster_to_downstream
# ============================================================
log ""
log "--- Stage 8: map_cluster_to_downstream ---"
ENRICHED_TSV="${OUTPUT_DIR}/downstream_genes_enriched.tsv"

if [[ -f "$ENRICHED_TSV" ]]; then
    log "[SKIP] $ENRICHED_TSV already exists"
else
    python3 "${CODE_DIR}/map_cluster_to_downstream.py" \
        --clusters   "$CLUSTER_TSV_FOR_MAP" \
        --downstream "$DOWNSTREAM_TSV" \
        --output     "$ENRICHED_TSV" \
        2>&1 | tee -a "$LOG"
fi

log "Stage 8 complete: $ENRICHED_TSV"

# ============================================================
# Stage 9: analyze_gene_chains
# ============================================================
log ""
log "--- Stage 9: analyze_gene_chains ---"
CHAIN_DIR="${OUTPUT_DIR}/gene_chain_${CHAIN_MAX_RANK}"

CHAIN_EXTRA=""
[[ -n "$BLAST_NAMES" ]] && CHAIN_EXTRA="$CHAIN_EXTRA --blast-names ${BASE_DIR}/${BLAST_NAMES}"
[[ -n "$END_AT_GENE" ]] && CHAIN_EXTRA="$CHAIN_EXTRA --end-at ${END_AT_GENE} --same-strand-only"
[[ -n "$STOP_BEFORE_GENE" ]] && CHAIN_EXTRA="$CHAIN_EXTRA --stop-before ${STOP_BEFORE_GENE}"
[[ -n "$MASTER_COLOR_MAP" ]] && CHAIN_EXTRA="$CHAIN_EXTRA --master-colors ${BASE_DIR}/${MASTER_COLOR_MAP}"

python3 "${CODE_DIR}/analyze_gene_chains.py" \
    --input      "$ENRICHED_TSV" \
    --output-dir "$CHAIN_DIR" \
    --max-rank   "$CHAIN_MAX_RANK" \
    $CHAIN_EXTRA \
    2>&1 | tee -a "$LOG"

log "Stage 9 complete: $CHAIN_DIR"

# ============================================================
# Stage 9b: unique_patterns_with_strains — add all strains per pattern
# ============================================================
log ""
log "--- Stage 9b: unique_patterns_with_strains ---"
UNIQUE_PATTERNS_TSV="${CHAIN_DIR}/unique_patterns.tsv"
CHAIN_PATTERNS_TSV="${CHAIN_DIR}/gene_chain_patterns.tsv"
STRAINS_TSV="${CHAIN_DIR}/unique_patterns_with_strains.tsv"

if [[ ! -f "$UNIQUE_PATTERNS_TSV" || ! -f "$CHAIN_PATTERNS_TSV" ]]; then
    log "[SKIP] Missing input files for stage 9b (unique_patterns.tsv or gene_chain_patterns.tsv)"
else
    python3 - "$UNIQUE_PATTERNS_TSV" "$CHAIN_PATTERNS_TSV" "$STRAINS_TSV" <<'PYEOF'
import csv, sys
from collections import defaultdict

unique_tsv, chain_tsv, out_tsv = sys.argv[1], sys.argv[2], sys.argv[3]

pattern_strains = defaultdict(list)
with open(chain_tsv) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        pattern_strains[row["pattern_type"]].append(row["strain"])

with open(unique_tsv) as fin, open(out_tsv, "w", newline="") as fout:
    reader = csv.DictReader(fin, delimiter="\t")
    writer = csv.DictWriter(fout, fieldnames=reader.fieldnames + ["all_strains"],
                            delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in reader:
        row["all_strains"] = ";".join(pattern_strains.get(row["pattern_type"], []))
        writer.writerow(row)
PYEOF
    log "Stage 9b complete: $STRAINS_TSV"
fi

log ""
log "============================================"
log "Analysis complete. Log: $LOG"
log "============================================"
