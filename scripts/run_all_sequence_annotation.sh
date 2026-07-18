#!/usr/bin/env bash
# ============================================================
# run_all_sequence_annotation.sh
# Runs Prokka + EggNOG annotation for all queries from allele_search.sh
#
# Queries:     rgg144, rgg939, rgg1518, TprA
# Genome sets: S_mitis (sm), S_pneumo (sp), S_oralis (so)
#
# For each query:
#   Step A: Prokka annotate all genomes + extract neighborhood FAA
#   Step B: Pool unique proteins + run EggNOG-mapper
#
# Generated per-query YAML configs are written to scripts/configs/
#
# Usage:
#   bash scripts/run_all_sequence_annotation.sh [--test]
#   --test: passes --test to run_annotate_prokka.sh (10 genomes per set)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG_DIR="${SCRIPT_DIR}/configs"

TEST_MODE=false
for _arg in "$@"; do
    [[ "$_arg" == "--test" ]] && TEST_MODE=true
done

mkdir -p "$CONFIG_DIR"

# ============================================================
# Paths — derived from allele_search.sh
# ============================================================
ALLELE_SEARCH_RESULTS="/mnt/extra_space/xinruoz/allele_search/results"

SM_GENOME_DIR="/mnt/extra_space/database/SMitis_Database_pubMLST_102225"
SP_GENOME_DIR="/mnt/extra_space/database/GoldenSet_7548Genomes_2023"
SO_GENOME_DIR="/mnt/extra_space/database/OralisSetX"

OUTPUT_DIR="${BASE_DIR}/annotation_results"

# ============================================================
# Stable settings
# ============================================================
HMM_FILE="/mnt/extra_space/database/hmm/eggNOG_strep_combined.hmm"
EGGNOG_DB="/mnt/extra_space/database/eggnog_data/emapperdb-5.0.2"
CONDA_ENV="ud_search"
PROKKA_CPUS=4
PROKKA_JOBS=5
EGGNOG_CPUS=20
MAX_GENES=20

# ============================================================
# Queries (order matches allele_search.sh)
# ============================================================
QUERIES=(rgg144 rgg939 rgg1518 TprA)

echo "============================================"
echo "run_all_sequence_annotation.sh"
echo "Queries:    ${QUERIES[*]}"
echo "Output dir: $OUTPUT_DIR"
echo "Started:    $(date)"
[[ "$TEST_MODE" == true ]] && echo "*** TEST MODE ***"
echo "============================================"

for query in "${QUERIES[@]}"; do

    echo ""
    echo "============================================"
    echo "Query: $query  ($(date))"
    echo "============================================"

    PROKKA_CONFIG="${CONFIG_DIR}/${query}_prokka.yaml"
    EGGNOG_CONFIG="${CONFIG_DIR}/${query}_eggnog.yaml"

    # ----------------------------------------------------------
    # Generate Prokka config for this query
    # ----------------------------------------------------------
    cat > "$PROKKA_CONFIG" << YAML
target_gene: ${query}
output_dir: ${OUTPUT_DIR}
max_genes: ${MAX_GENES}
genome_sets:
  - genome_name: S_mitis
    genome_dir:  ${SM_GENOME_DIR}
    target_fasta: ${ALLELE_SEARCH_RESULTS}/${query}/sm_${query}/standard_dna_seq/${query}.fasta
  - genome_name: S_pneumo
    genome_dir:  ${SP_GENOME_DIR}
    target_fasta: ${ALLELE_SEARCH_RESULTS}/${query}/sp_${query}/standard_dna_seq/${query}.fasta
  - genome_name: S_oralis
    genome_dir:  ${SO_GENOME_DIR}
    target_fasta: ${ALLELE_SEARCH_RESULTS}/${query}/so_${query}/standard_dna_seq/${query}.fasta
hmm_file: ${HMM_FILE}
conda_env: ${CONDA_ENV}
prokka_cpus: ${PROKKA_CPUS}
prokka_parallel_jobs: ${PROKKA_JOBS}
YAML

    # ----------------------------------------------------------
    # Generate EggNOG config for this query
    # ----------------------------------------------------------
    cat > "$EGGNOG_CONFIG" << YAML
target_gene: ${query}
output_dir: ${OUTPUT_DIR}
genome_sets:
  - genome_name: S_mitis
  - genome_name: S_pneumo
  - genome_name: S_oralis
eggnog_db: ${EGGNOG_DB}
conda_env: ${CONDA_ENV}
eggnog_cpus: ${EGGNOG_CPUS}
YAML

    # ----------------------------------------------------------
    # Step A: Prokka + neighborhood FAA extraction
    # ----------------------------------------------------------
    echo ""
    echo "--- Step A: Prokka ($query) ---"
    if [[ "$TEST_MODE" == true ]]; then
        bash "$SCRIPT_DIR/run_annotate_prokka.sh" "$PROKKA_CONFIG" --test
    else
        bash "$SCRIPT_DIR/run_annotate_prokka.sh" "$PROKKA_CONFIG"
    fi

    # ----------------------------------------------------------
    # Step B: Pool proteins + EggNOG-mapper
    # ----------------------------------------------------------
    echo ""
    echo "--- Step B: EggNOG ($query) ---"
    bash "$SCRIPT_DIR/run_annotate_eggnog.sh" "$EGGNOG_CONFIG"

    echo ""
    echo "[DONE] Query $query complete: $(date)"

done

echo ""
echo "============================================"
echo "All queries complete: $(date)"
echo "Configs written to: $CONFIG_DIR"
echo "Results under:      $OUTPUT_DIR"
echo "============================================"
