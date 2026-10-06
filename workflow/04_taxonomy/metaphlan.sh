#!/bin/bash
#SBATCH -A lu2026-2-117
#SBATCH -J metaphlan
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH -t 08:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#SBATCH --mail-user=ch2068va-s@student.lu.se
#SBATCH --mail-type=END,FAIL
# =============================================================================
# 04_taxonomy/metaphlan.sh — MetaPhlAn 4 taxonomic profile per sample
#
# Input : CLEAN_FASTQ_DIR (KneadData output; paired + unmatched reads)
# Output: METAPHLAN_DIR/profiles/<S>_profile.txt
#         METAPHLAN_DIR/mapout/<S>.mapout.bz2   (allows fast re-profiling)
#
# Usage:
#   ./submit.sh --array-samples workflow/04_taxonomy/metaphlan.sh
# Then merge into abundance tables with merge_metaphlan.sh.
# =============================================================================
set -euo pipefail
# Find config.sh: submit from the repository folder (SLURM_SUBMIT_DIR).
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${PIPELINE_DIR}/config.sh" ]] || { echo "config.sh not found in ${PIPELINE_DIR}: submit from the repository folder" >&2; exit 1; }
source "${PIPELINE_DIR}/config.sh"

SAMPLE="$(resolve_item "${SELECTED_SAMPLES_FILE}" "${1:-}")"
T="$(threads)"
require_vars METAPHLAN_DB_DIR CONDA_MODULE
require_paths METAPHLAN_DB_DIR

PROF_DIR="${METAPHLAN_DIR}/profiles"
MAP_DIR="${METAPHLAN_DIR}/mapout"
make_out_dirs "${PROF_DIR}" "${MAP_DIR}"

if [[ -s "${PROF_DIR}/${SAMPLE}_profile.txt" ]]; then
    log "[skip] ${SAMPLE} already profiled"; exit 0
fi

# ---- Inputs: every non-empty cleaned read file of the sample -------------------
INPUTS=()
for f in "${CLEAN_FASTQ_DIR}/${SAMPLE}_paired_1.fastq.gz" \
         "${CLEAN_FASTQ_DIR}/${SAMPLE}_paired_2.fastq.gz" \
         "${CLEAN_FASTQ_DIR}/${SAMPLE}_unmatched_1.fastq.gz" \
         "${CLEAN_FASTQ_DIR}/${SAMPLE}_unmatched_2.fastq.gz" \
         "${CLEAN_FASTQ_DIR}/${SAMPLE}.fastq.gz"; do
    if nonempty_gz "$f"; then INPUTS+=("$f"); fi
done
(( ${#INPUTS[@]} > 0 )) || die "No cleaned reads for ${SAMPLE} in ${CLEAN_FASTQ_DIR} (run 03_preprocessing)"

make_scratch metaphlan "${SAMPLE}"
trap 'rm -rf "${WORK}"' EXIT
LOCAL=()
for f in "${INPUTS[@]}"; do
    cp "$f" "${WORK}/in/"
    LOCAL+=("${WORK}/in/$(basename "$f")")
done
IN_LIST="$(IFS=,; echo "${LOCAL[*]}")"   # MetaPhlAn takes a comma-separated list

load_conda_env "${ENV_METAPHLAN}"
metaphlan --version 2>&1 | head -n1 >&2 || true

# MetaPhlAn >= 4.2 renamed --bowtie2db/--bowtie2out to --db_dir/--mapout.
if metaphlan --help 2>&1 | grep -q -- '--mapout'; then
    DB_FLAG=(--db_dir "${METAPHLAN_DB_DIR}")
    MAP_FLAG=(--mapout "${WORK}/out/${SAMPLE}.mapout.bz2")
else
    DB_FLAG=(--bowtie2db "${METAPHLAN_DB_DIR}")
    MAP_FLAG=(--bowtie2out "${WORK}/out/${SAMPLE}.mapout.bz2")
fi

metaphlan "${IN_LIST}" \
    --input_type fastq \
    "${DB_FLAG[@]}" \
    --index "${METAPHLAN_INDEX}" \
    --offline \
    --nproc "${T}" \
    --tmp_dir "${WORK}/tmp" \
    "${MAP_FLAG[@]}" \
    -o "${WORK}/out/${SAMPLE}_profile.txt"

[[ -s "${WORK}/out/${SAMPLE}_profile.txt" ]] || die "MetaPhlAn produced no profile for ${SAMPLE}"
cp "${WORK}/out/${SAMPLE}.mapout.bz2" "${MAP_DIR}/" 2>/dev/null || log "WARNING: no mapout file"
cp "${WORK}/out/${SAMPLE}_profile.txt" "${PROF_DIR}/"   # copied last = completion marker
log "Finished ${SAMPLE}"
