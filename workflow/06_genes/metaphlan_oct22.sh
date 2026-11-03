#!/bin/bash
#SBATCH -A lu2026-2-117
#SBATCH -J mpa_oct22
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH -t 08:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#SBATCH --mail-user=ch2068va-s@student.lu.se
#SBATCH --mail-type=END,FAIL
# =============================================================================
# 06_genes/metaphlan_oct22.sh — MetaPhlAn profiles with the Oct22 index,
#                               used ONLY as input for HUMAnN 4
#
# HUMAnN 4 needs profiles made with mpa_vOct22_CHOCOPhlAnSGB_202403 (as in
# CRASH's humann_meta.sh). Uses all non-empty cleaned read files of a sample,
# like 04_taxonomy. The Jun23 profiles in 04_taxonomy stay the taxonomy result.
#
#   export DATASET=<name>
#   ./submit.sh --array-samples workflow/06_genes/metaphlan_oct22.sh
# Output: <dataset>/06_genes/metaphlan_oct22/profiles/<S>_profile.txt
# =============================================================================
set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${PIPELINE_DIR}/config.sh" ]] || { echo "config.sh not found in ${PIPELINE_DIR}: submit from the repository folder" >&2; exit 1; }
source "${PIPELINE_DIR}/config.sh"

SAMPLE="$(resolve_item "${SELECTED_SAMPLES_FILE}" "${1:-}")"
T="$(threads)"
require_vars ENV_HUMANN MPA_OCT22_DB_DIR
[[ -s "${MPA_OCT22_DB_DIR}/${MPA_OCT22_INDEX}.pkl" ]] || die "Oct22 index not found: ${MPA_OCT22_DB_DIR}/${MPA_OCT22_INDEX}.pkl"
PROF_DIR="${GENES_DIR}/metaphlan_oct22/profiles"
make_out_dirs "${PROF_DIR}"

if [[ -s "${PROF_DIR}/${SAMPLE}_profile.txt" ]]; then
    log "[skip] ${SAMPLE} already profiled"; exit 0
fi

INPUTS=()
for f in "${CLEAN_FASTQ_DIR}/${SAMPLE}_paired_1.fastq.gz" "${CLEAN_FASTQ_DIR}/${SAMPLE}_paired_2.fastq.gz" \
         "${CLEAN_FASTQ_DIR}/${SAMPLE}_unmatched_1.fastq.gz" "${CLEAN_FASTQ_DIR}/${SAMPLE}_unmatched_2.fastq.gz" \
         "${CLEAN_FASTQ_DIR}/${SAMPLE}.fastq.gz"; do
    if nonempty_gz "$f"; then INPUTS+=("$f"); fi
done
(( ${#INPUTS[@]} > 0 )) || die "No cleaned reads for ${SAMPLE} in ${CLEAN_FASTQ_DIR}"

make_scratch mpa_oct22 "${SAMPLE}"
trap 'rm -rf "${WORK}"' EXIT
LOCAL=()
for f in "${INPUTS[@]}"; do cp "$f" "${WORK}/in/"; LOCAL+=("${WORK}/in/$(basename "$f")"); done
IN_LIST="$(IFS=,; echo "${LOCAL[*]}")"

load_conda_env "${ENV_HUMANN}"
metaphlan --version 2>&1 | head -n 1 >&2 || true

# MetaPhlAn 4.1.1 option names (--bowtie2db / -x), settings as in CRASH.
metaphlan "${IN_LIST}" \
    --input_type fastq \
    --bowtie2db "${MPA_OCT22_DB_DIR}" \
    -x "${MPA_OCT22_INDEX}" \
    --offline \
    --no_map \
    --nproc "${T}" \
    --tmp_dir "${WORK}/tmp" \
    -t rel_ab_w_read_stats \
    -o "${WORK}/out/${SAMPLE}_profile.txt"

[[ -s "${WORK}/out/${SAMPLE}_profile.txt" ]] || die "No Oct22 profile produced for ${SAMPLE}"
cp "${WORK}/out/${SAMPLE}_profile.txt" "${PROF_DIR}/"
log "Finished ${SAMPLE}"
