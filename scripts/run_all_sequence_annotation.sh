#!/usr/bin/env bash
# ============================================================
# run_all_sequence_annotation.sh
# Runs Prokka + EggNOG annotation for all queries from allele_search.sh
#
# Queries:     rgg144, rgg939, rgg1518, TprA, TprB, TprC
# Genome sets: S_mitis (sm), S_pneumo (sp), S_oralis (so), S_therm (sth), Portugal (port)
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
DOWNSTREAM_ALLELE_RESULTS="/mnt/extra_space/xinruoz/downstream_analysis_paper/allele_results"

# Per-genome-set allele results base directory (S_therm and Portugal live in a different root)
declare -A ALLELE_BASE=([S_mitis]="${ALLELE_SEARCH_RESULTS}" [S_pneumo]="${ALLELE_SEARCH_RESULTS}" [S_oralis]="${ALLELE_SEARCH_RESULTS}" [S_therm]="${DOWNSTREAM_ALLELE_RESULTS}" [Portugal]="${DOWNSTREAM_ALLELE_RESULTS}")

# Central Prokka annotation directory — created on first run, shared across all queries.
# Structure: {PROKKA_BASE_DIR}/{species_folder}/{strain}/{strain}.gff/.faa
PROKKA_BASE_DIR="/mnt/extra_space/xinruoz/ud_search/annotation_prokka"

# Mapping: genome set name -> species subfolder in PROKKA_BASE_DIR
declare -A SPECIES_FOLDER=([S_mitis]=mitis [S_pneumo]=pneumo [S_oralis]=oralis [S_therm]=therm [Portugal]=port)

# Genome FASTA directories (used to re-run Prokka on strains with missing GFF)
SM_GENOME_DIR="/mnt/extra_space/database/SMitis_Database_pubMLST_102225"
SP_GENOME_DIR="/mnt/extra_space/database/GoldenSet_7548Genomes_2023"
SO_GENOME_DIR="/mnt/extra_space/database/OralisSetX"
STH_GENOME_DIR="/mnt/extra_space/database/ThermSetX"
PORT_GENOME_DIR="/mnt/extra_space/database/PortugalStrains"

OUTPUT_DIR="/mnt/extra_space/xinruoz/ud_search/annotation_results"

# ============================================================
# Stable settings
# ============================================================
HMM_FILE="/mnt/extra_space/database/hmm/eggNOG_strep_combined.hmm"
EGGNOG_DB="/mnt/extra_space/database/eggnog_data/emapperdb-5.0.2"
CONDA_ENV="ud_search"
PROKKA_CPUS=4
_ncpu=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
if [[ $_ncpu -le 8 ]]; then
    PARALLEL_JOBS=$(( _ncpu - 2 ))
else
    PARALLEL_JOBS=$(( _ncpu - 4 ))
fi
[[ $PARALLEL_JOBS -lt 1 ]] && PARALLEL_JOBS=1
EGGNOG_CPUS=20
MAX_GENES=20

# ============================================================
# Queries (order matches allele_search.sh)
# ============================================================
QUERIES=(rgg144 rgg939 rgg1518 TprA TprB TprC TprA2 RtgR)

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
prokka_base_dir: ${PROKKA_BASE_DIR}
hmm_file: ${HMM_FILE}
prokka_cpus: ${PROKKA_CPUS}
conda_env: ${CONDA_ENV}
parallel_jobs: ${PARALLEL_JOBS}
genome_sets:
  - genome_name: S_mitis
    prokka_species_name: ${SPECIES_FOLDER[S_mitis]}
    genome_dir: ${SM_GENOME_DIR}
    target_fasta: ${ALLELE_SEARCH_RESULTS}/${query}/sm_${query}/standard_dna_seq/${query}.fasta
  - genome_name: S_pneumo
    prokka_species_name: ${SPECIES_FOLDER[S_pneumo]}
    genome_dir: ${SP_GENOME_DIR}
    target_fasta: ${ALLELE_SEARCH_RESULTS}/${query}/sp_${query}/standard_dna_seq/${query}.fasta
  - genome_name: S_oralis
    prokka_species_name: ${SPECIES_FOLDER[S_oralis]}
    genome_dir: ${SO_GENOME_DIR}
    target_fasta: ${ALLELE_SEARCH_RESULTS}/${query}/so_${query}/standard_dna_seq/${query}.fasta
  - genome_name: S_therm
    prokka_species_name: ${SPECIES_FOLDER[S_therm]}
    genome_dir: ${STH_GENOME_DIR}
    target_fasta: ${ALLELE_BASE[S_therm]}/${query}/sth_${query}/standard_dna_seq/${query}.fasta
  - genome_name: Portugal
    prokka_species_name: ${SPECIES_FOLDER[Portugal]}
    genome_dir: ${PORT_GENOME_DIR}
    target_fasta: ${ALLELE_BASE[Portugal]}/${query}/port_${query}/standard_dna_seq/${query}.fasta
YAML

    # ----------------------------------------------------------
    # Generate EggNOG config for this query
    # ----------------------------------------------------------
    cat > "$EGGNOG_CONFIG" << YAML
target_gene: ${query}
output_dir: ${OUTPUT_DIR}
genome_sets:
  - genome_name: S_mitis
    prokka_dir: ${PROKKA_BASE_DIR}/${SPECIES_FOLDER[S_mitis]}
  - genome_name: S_pneumo
    prokka_dir: ${PROKKA_BASE_DIR}/${SPECIES_FOLDER[S_pneumo]}
  - genome_name: S_oralis
    prokka_dir: ${PROKKA_BASE_DIR}/${SPECIES_FOLDER[S_oralis]}
  - genome_name: S_therm
    prokka_dir: ${PROKKA_BASE_DIR}/${SPECIES_FOLDER[S_therm]}
  - genome_name: Portugal
    prokka_dir: ${PROKKA_BASE_DIR}/${SPECIES_FOLDER[Portugal]}
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
