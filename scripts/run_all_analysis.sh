#!/usr/bin/env bash
# ============================================================
# run_all_analysis.sh
# Runs post-annotation analysis (stages 3-9) for every
# query × genome-set combination produced by
# run_all_sequence_annotation.sh.
#
# Queries:     rgg144, rgg939, rgg1518, TprA, TprB, TprC
# Genome sets: S_mitis, S_pneumo, S_oralis, S_therm, Portugal
#
# For each combination it generates a per-combo analysis config
# in scripts/configs/ and calls run_analysis.sh.
# Combinations whose output dir does not yet exist are skipped.
#
# Usage:
#   bash scripts/run_all_analysis.sh [options]
#
# Options:
#   --queries  q1,q2,...   comma-separated subset of queries  (default: all)
#   --sets     S1,S2,...   comma-separated subset of genome sets (default: all)
#   --dry-run              print what would run, without running it
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG_DIR="${SCRIPT_DIR}/configs"

# ============================================================
# Paths — must match run_all_sequence_annotation.sh
# ============================================================
ALLELE_SEARCH_RESULTS="/mnt/extra_space/xinruoz/allele_search/results"
DOWNSTREAM_ALLELE_RESULTS="/mnt/extra_space/xinruoz/downstream_analysis_paper/allele_results"
PROKKA_BASE_DIR="/mnt/extra_space/xinruoz/ud_search/annotation_prokka"
OUTPUT_DIR="/mnt/extra_space/xinruoz/ud_search/annotation_results"

# Species subfolder inside PROKKA_BASE_DIR
declare -A SPECIES_FOLDER=([S_mitis]=mitis [S_pneumo]=pneumo [S_oralis]=oralis [S_therm]=therm [Portugal]=port)

# Allele-search result prefix per genome set
declare -A ALLELE_PREFIX=([S_mitis]=sm [S_pneumo]=sp [S_oralis]=so [S_therm]=sth [Portugal]=port)

# Per-genome-set allele results base directory (S_therm and Portugal live in a different root)
declare -A ALLELE_BASE=([S_mitis]="${ALLELE_SEARCH_RESULTS}" [S_pneumo]="${ALLELE_SEARCH_RESULTS}" [S_oralis]="${ALLELE_SEARCH_RESULTS}" [S_therm]="${DOWNSTREAM_ALLELE_RESULTS}" [Portugal]="${DOWNSTREAM_ALLELE_RESULTS}")

# ============================================================
# Analysis parameters (shared across all combos)
# ============================================================
CONDA_ENV="ud_search"
WINDOW=20000
MAX_RANK=10
CDHIT_IDENTITY=0.80
CDHIT_THREADS=4
CHAIN_MAX_RANK=10
_ncpu=$(nproc 2>/dev/null || echo 4)
INTERPROSCAN_CPUS=$(( _ncpu < 8 ? _ncpu : 8 ))

# ============================================================
# Defaults
# ============================================================
# Per-query search direction:
#   downstream = regulated genes are 3' of the query (rgg144, TprA)
#   upstream   = divergent regulation — genes on the opposite/upstream side (rgg939, rgg1518)
#   both       = unknown direction, search both sides (TprB, TprC)
declare -A QUERY_DIRECTION=(
    [rgg144]=downstream
    [rgg939]=upstream
    [rgg1518]=upstream
    [TprA]=downstream
    [TprB]=both
    [TprC]=both
    [TprA2]=downstream
    [RtgR]=upstream
)

ALL_QUERIES=(rgg144 rgg939 rgg1518 TprA TprB TprC TprA2 RtgR)
ALL_SETS=(S_mitis S_pneumo S_oralis S_therm Portugal)
DRY_RUN=false

# ============================================================
# Argument parsing
# ============================================================
SELECTED_QUERIES=()
SELECTED_SETS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --queries)
            IFS=',' read -ra SELECTED_QUERIES <<< "$2"; shift 2 ;;
        --sets)
            IFS=',' read -ra SELECTED_SETS <<< "$2"; shift 2 ;;
        --dry-run)
            DRY_RUN=true; shift ;;
        *)
            echo "Unknown option: $1"; exit 1 ;;
    esac
done

[[ ${#SELECTED_QUERIES[@]} -eq 0 ]] && SELECTED_QUERIES=("${ALL_QUERIES[@]}")
[[ ${#SELECTED_SETS[@]}    -eq 0 ]] && SELECTED_SETS=("${ALL_SETS[@]}")

mkdir -p "$CONFIG_DIR"

echo "============================================"
echo "run_all_analysis.sh"
echo "Queries:    ${SELECTED_QUERIES[*]}"
echo "Sets:       ${SELECTED_SETS[*]}"
echo "Output dir: $OUTPUT_DIR"
echo "Started:    $(date)"
[[ "$DRY_RUN" == true ]] && echo "*** DRY RUN ***"
echo "============================================"

# ============================================================
# Main loop
# ============================================================
n_done=0
n_skip=0

for query in "${SELECTED_QUERIES[@]}"; do
    for gset in "${SELECTED_SETS[@]}"; do

        SET_DIR="${OUTPUT_DIR}/${gset}_${query}"

        if [[ ! -d "$SET_DIR" ]]; then
            echo "[SKIP] $gset/$query: output dir not found ($SET_DIR)"
            n_skip=$(( n_skip + 1 ))
            continue
        fi

        EGGNOG_FILE="${SET_DIR}/eggnog/${gset}_${query}.emapper.annotations"
        if [[ ! -f "$EGGNOG_FILE" ]]; then
            echo "[SKIP] $gset/$query: EggNOG annotations not found ($EGGNOG_FILE) — run annotation first"
            n_skip=$(( n_skip + 1 ))
            continue
        fi

        DIRECTION="${QUERY_DIRECTION[$query]}"
        prefix="${ALLELE_PREFIX[$gset]}"
        _gene_short=$(echo "$query" | sed 's/^[Rr][Gg][Gg]//')
        [[ "$_gene_short" == "$query" ]] && _gene_short=$(echo "$query" | tr '[:upper:]' '[:lower:]')
        CLUSTER_PREFIX="D${_gene_short}_${prefix}"
        species="${SPECIES_FOLDER[$gset]}"
        ANNO_DIR="${PROKKA_BASE_DIR}/${species}"
        QUERY_FASTA="${ALLELE_BASE[$gset]}/${query}/${prefix}_${query}/standard_dna_seq/${query}.fasta"

        ANALYSIS_CONFIG="${CONFIG_DIR}/${gset}_${query}_analysis.yaml"

        cat > "$ANALYSIS_CONFIG" << YAML
# Auto-generated by run_all_analysis.sh — $(date)
# Genome set: ${gset}   Query: ${query}

base_dir: ${SET_DIR}

annotation_dir: ${ANNO_DIR}

eggnog_annotations: eggnog/${gset}_${query}.emapper.annotations
protein_origins:    eggnog/protein_origins.tsv

query_fasta: ${QUERY_FASTA}

annotation_with_eggnog_dir: annotation_with_eggnog
output_dir:                 output
log_dir:                    log

window:    ${WINDOW}
direction: ${DIRECTION}

max_rank:       ${MAX_RANK}
cdhit_identity: ${CDHIT_IDENTITY}
cdhit_threads:  ${CDHIT_THREADS}
cluster_prefix: ${CLUSTER_PREFIX}

chain_max_rank: ${CHAIN_MAX_RANK}

interproscan_cpus: ${INTERPROSCAN_CPUS}

conda_env: ${CONDA_ENV}
YAML

        echo ""
        echo "============================================"
        echo "Analyzing: $gset / $query  ($(date))"
        echo "Config:    $ANALYSIS_CONFIG"
        echo "============================================"

        if [[ "$DRY_RUN" == true ]]; then
            echo "[DRY RUN] Would run: bash ${SCRIPT_DIR}/run_analysis.sh $ANALYSIS_CONFIG"
        else
            bash "${SCRIPT_DIR}/run_analysis.sh" "$ANALYSIS_CONFIG"
        fi

        n_done=$(( n_done + 1 ))

    done
done

echo ""
echo "============================================"
echo "All done: $(date)"
echo "Processed: $n_done    Skipped: $n_skip"
echo "Configs:   $CONFIG_DIR"
echo "============================================"
