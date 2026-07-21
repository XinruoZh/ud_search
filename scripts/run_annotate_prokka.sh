#!/usr/bin/env bash
# ============================================================
# run_annotate_prokka.sh — Neighborhood FAA extraction using central Prokka annotations
#
# Prokka is run ONCE per genome and stored centrally.
# This script reads from prokka_base_dir/{prokka_species_name}/{strain}/
# and extracts neighborhood FAA per strain for each target gene query.
#
# If a strain directory exists but is missing a .gff (Prokka previously failed
# due to contig names >37 chars), this script automatically re-runs Prokka
# with shortened contig names and restores original names in the output GFF.
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
    # genome_name\tprokka_species_name\ttarget_fasta\tgenome_dir
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
    name       = str(entry.get('genome_name', '')).strip()
    species    = str(entry.get('prokka_species_name', '')).strip()
    tfasta     = str(entry.get('target_fasta', '')).strip()
    genome_dir = str(entry.get('genome_dir', '')).strip()
    if name and species and tfasta:
        print(f"{name}\t{species}\t{tfasta}\t{genome_dir}")
PYEOF
}

# ============================================================
# Load scalar config values
# ============================================================
TARGET_GENE=$(parse_yaml_scalar target_gene)
OUTPUT_DIR=$(parse_yaml_scalar output_dir)
MAX_GENES=$(parse_yaml_scalar max_genes)
PROKKA_BASE_DIR=$(parse_yaml_scalar prokka_base_dir)
HMM_FILE=$(parse_yaml_scalar hmm_file)
CONDA_ENV=$(parse_yaml_scalar conda_env)
PROKKA_CPUS=$(parse_yaml_scalar prokka_cpus)
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
# Python snippets for contig-name repair (exported as variables)
# ============================================================

# Reads genome FASTA, shortens any contig name >37 chars to rctg{N:05d},
# writes temp FASTA and mapping file.
# Exit 0 = some names were shortened; exit 2 = nothing to shorten.
SHORTEN_PY='
import sys

fasta_in  = sys.argv[1]
fasta_out = sys.argv[2]
map_out   = sys.argv[3]

mapping = {}
idx = 0
lines_out = []

with open(fasta_in) as f:
    for line in f:
        if line.startswith(">"):
            original = line[1:].split()[0]
            rest = line[1 + len(original):]
            if len(original) > 37:
                idx += 1
                short = f"rctg{idx:05d}"
                mapping[short] = original
                lines_out.append(f">{short}{rest}")
            else:
                lines_out.append(line)
        else:
            lines_out.append(line)

if not mapping:
    sys.exit(2)  # no long names found — Prokka failed for another reason

with open(fasta_out, "w") as f:
    f.writelines(lines_out)
with open(map_out, "w") as f:
    for short, orig in mapping.items():
        f.write(f"{short}\t{orig}\n")
'

# Reads Prokka GFF with shortened contig names, restores original names.
RESTORE_PY='
import sys

gff_in   = sys.argv[1]
gff_out  = sys.argv[2]
map_file = sys.argv[3]

mapping = {}
with open(map_file) as f:
    for line in f:
        short, original = line.rstrip("\n").split("\t", 1)
        mapping[short] = original

in_fasta = False
with open(gff_in) as fin, open(gff_out, "w") as fout:
    for line in fin:
        if line.startswith("##FASTA"):
            in_fasta = True
            fout.write(line)
            continue
        if in_fasta:
            if line.startswith(">"):
                name = line[1:].rstrip("\n")
                fout.write(f">{mapping.get(name, name)}\n")
            else:
                fout.write(line)
        elif line.startswith("#"):
            fout.write(line)
        elif "\t" in line:
            parts = line.split("\t")
            parts[0] = mapping.get(parts[0], parts[0])
            fout.write("\t".join(parts))
        else:
            fout.write(line)
'

export SHORTEN_PY RESTORE_PY

# ============================================================
# Prokka repair function (exported for GNU parallel)
# ============================================================

run_prokka_repair() {
    local name=$1
    local anno_dir=$2
    local genome_dir=$3
    local hmm_file=$4
    local prokka_cpus=$5
    local log=$6

    local genome_fasta="${genome_dir}/${name}.fasta"
    local tmpdir="/tmp/prokka_repair_${name}"
    local tmp_fasta="${tmpdir}/input.fasta"
    local map_file="${tmpdir}/contig_map.tsv"
    local prokka_out="${tmpdir}/out"
    local raw_gff="${prokka_out}/${name}.gff"
    local fixed_gff="${tmpdir}/${name}_fixed.gff"
    local strain_dir="${anno_dir}/${name}"

    if [[ ! -f "$genome_fasta" ]]; then
        echo "[WARN] Repair $name: genome FASTA not found at $genome_fasta — skipping" | tee -a "$log"
        return 0
    fi

    echo "[REPAIR] $name: re-running Prokka with contig-name fix" | tee -a "$log"
    rm -rf "$tmpdir"
    mkdir -p "$tmpdir" "$prokka_out"

    # Step 1: shorten contig names
    python3 -c "$SHORTEN_PY" "$genome_fasta" "$tmp_fasta" "$map_file"
    local rc=$?
    if [[ $rc -eq 2 ]]; then
        echo "[WARN] Repair $name: no contig names >37 chars — Prokka failed for another reason, skipping" | tee -a "$log"
        rm -rf "$tmpdir"
        return 0
    fi

    # Step 2: run Prokka on shortened FASTA
    prokka \
        --cpus   "$prokka_cpus" \
        --outdir "$prokka_out" \
        --prefix "$name" \
        --hmms   "$hmm_file" \
        --quiet \
        --force \
        "$tmp_fasta" \
        >> "$log" 2>&1 || {
        echo "[ERROR] Repair $name: Prokka failed — tmpdir kept for debugging: $tmpdir" | tee -a "$log"
        return 1
    }

    # Step 3: restore original contig names in GFF
    python3 -c "$RESTORE_PY" "$raw_gff" "$fixed_gff" "$map_file"

    # Step 4: copy outputs to central strain dir
    cp "$fixed_gff" "${strain_dir}/${name}.gff"
    cp "${prokka_out}/${name}.faa" "${strain_dir}/${name}.faa"
    for ext in err ffn fna gbk log sqn tbl tsv txt; do
        [[ -f "${prokka_out}/${name}.${ext}" ]] && \
            cp "${prokka_out}/${name}.${ext}" "${strain_dir}/${name}.${ext}"
    done

    rm -rf "$tmpdir"
    echo "[REPAIRED] $name: GFF restored with original contig names" | tee -a "$log"
}

export -f run_prokka_repair

# ============================================================
# Neighborhood extraction function (exported for GNU parallel)
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

while IFS=$'\t' read -r GENOME_NAME PROKKA_SPECIES TARGET_FASTA GENOME_DIR; do

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
    # Step 0: Repair strains with missing GFF (long contig names)
    # ----------------------------------------------------------
    if [[ -n "$GENOME_DIR" && -d "$GENOME_DIR" ]]; then

        # Collect strains missing a GFF
        MISSING_GFF=()
        while IFS= read -r d; do
            name=$(basename "$d")
            [[ ! -f "${d}/${name}.gff" ]] && MISSING_GFF+=("$name")
        done < <(find "$ANNO_DIR" -maxdepth 1 -mindepth 1 -type d | sort)

        if [[ ${#MISSING_GFF[@]} -gt 0 ]]; then
            echo "" | tee -a "$LOG"
            echo "--- Step 0: Prokka repair for ${#MISSING_GFF[@]} strains missing GFF ---" | tee -a "$LOG"

            if [[ "$TEST_MODE" == true ]]; then
                # In test mode limit to a few repairs so it doesn't take too long
                MISSING_GFF=("${MISSING_GFF[@]:0:3}")
            fi

            printf '%s\n' "${MISSING_GFF[@]}" | \
                parallel -j "$PARALLEL_JOBS" --line-buffer \
                    run_prokka_repair {} "$ANNO_DIR" "$GENOME_DIR" \
                    "$HMM_FILE" "$PROKKA_CPUS" "$LOG" || true

            n_repaired=$(for n in "${MISSING_GFF[@]}"; do
                [[ -f "${ANNO_DIR}/${n}/${n}.gff" ]] && echo ok; done | wc -l)
            echo "Step 0 complete: repaired=${n_repaired} of ${#MISSING_GFF[@]}" | tee -a "$LOG"
        else
            echo "" | tee -a "$LOG"
            echo "--- Step 0: No missing GFF files — skipping Prokka repair ---" | tee -a "$LOG"
        fi

    else
        echo "" | tee -a "$LOG"
        echo "--- Step 0: genome_dir not set or not found — skipping Prokka repair ---" | tee -a "$LOG"
    fi

    # ----------------------------------------------------------
    # Step 1: Extract neighborhood FAA (parallel)
    # ----------------------------------------------------------
    echo "" | tee -a "$LOG"
    echo "--- Step 1: Extract neighborhood FAA (parallel -j ${PARALLEL_JOBS}) ---" | tee -a "$LOG"

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
