#!/bin/bash
#SBATCH -A lu2026-2-117
#SBATCH -J humann
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -t 72:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#SBATCH --mail-user=ch2068va-s@student.lu.se
#SBATCH --mail-type=END,FAIL
# =============================================================================
# 06_genes/humann.sh — HUMAnN 4 functional profile per sample
#
# Settings as in CRASH's humann4.sh (24 threads per sample, prescreen 1.0),
# except that ALL non-empty cleaned read files are used (paired + unpaired),
# consistent with the MetaPhlAn runs. HUMAnN does not use read pairing, so the
# files are simply concatenated (gzip files can be joined with cat).
#
#   export DATASET=<name>
#   ./submit.sh --array-samples workflow/06_genes/humann.sh
#   (run metaphlan_oct22.sh first)
# Output: <dataset>/06_genes/humann/per_sample/<S>/<S>_*.tsv  (gene families,
#         reactions, pathways) and <S>.humann.log
# =============================================================================
set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${PIPELINE_DIR}/config.sh" ]] || { echo "config.sh not found in ${PIPELINE_DIR}: submit from the repository folder" >&2; exit 1; }
source "${PIPELINE_DIR}/config.sh"

SAMPLE="$(resolve_item "${SELECTED_SAMPLES_FILE}" "${1:-}")"
T="$(threads)"
require_vars ENV_HUMANN HUMANN_NUC_DB HUMANN_PROT_DB HUMANN_UTILITY_DB
require_paths HUMANN_NUC_DB HUMANN_PROT_DB HUMANN_UTILITY_DB
OUT_S="${HUMANN_DIR}/per_sample/${SAMPLE}"
make_out_dirs "${HUMANN_DIR}/per_sample"

if is_done "${HUMANN_DIR}" "${SAMPLE}"; then
    log "[skip] ${SAMPLE} already done"; exit 0
fi
PROFILE="${GENES_DIR}/metaphlan_oct22/profiles/${SAMPLE}_profile.txt"
[[ -s "${PROFILE}" ]] || die "Oct22 profile missing for ${SAMPLE} (run metaphlan_oct22.sh first): ${PROFILE}"

INPUTS=()
for f in "${CLEAN_FASTQ_DIR}/${SAMPLE}_paired_1.fastq.gz" "${CLEAN_FASTQ_DIR}/${SAMPLE}_paired_2.fastq.gz" \
         "${CLEAN_FASTQ_DIR}/${SAMPLE}_unmatched_1.fastq.gz" "${CLEAN_FASTQ_DIR}/${SAMPLE}_unmatched_2.fastq.gz" \
         "${CLEAN_FASTQ_DIR}/${SAMPLE}.fastq.gz"; do
    if nonempty_gz "$f"; then INPUTS+=("$f"); fi
done
(( ${#INPUTS[@]} > 0 )) || die "No cleaned reads for ${SAMPLE} in ${CLEAN_FASTQ_DIR}"

make_scratch humann "${SAMPLE}"
trap 'rm -rf "${WORK}"' EXIT
cat "${INPUTS[@]}" > "${WORK}/in/${SAMPLE}.fastq.gz"
log "${SAMPLE}: ${#INPUTS[@]} read files joined; ${T} threads on $(hostname)"

load_conda_env "${ENV_HUMANN}"
humann --version >&2 || true

humann \
    --input "${WORK}/in/${SAMPLE}.fastq.gz" \
    --taxonomic-profile "${PROFILE}" \
    --output "${WORK}/out" \
    --output-basename "${SAMPLE}" \
    --threads "${T}" \
    --nucleotide-database "${HUMANN_NUC_DB}" \
    --protein-database "${HUMANN_PROT_DB}" \
    --utility-database "${HUMANN_UTILITY_DB}" \
    --prescreen-threshold "${HUMANN_PRESCREEN}" \
    --remove-temp-output \
    2>&1 | tee "${WORK}/out/${SAMPLE}.humann.log"

shopt -s nullglob
gf=( "${WORK}/out/${SAMPLE}"_*genefamilies.tsv )
(( ${#gf[@]} > 0 )) || die "HUMAnN produced no gene-family table for ${SAMPLE}"
mkdir -p "${OUT_S}"
cp "${WORK}/out/${SAMPLE}"_*.tsv "${WORK}/out/${SAMPLE}.humann.log" "${OUT_S}/"
mark_done "${HUMANN_DIR}" "${SAMPLE}"
log "Finished ${SAMPLE}: $(find "${OUT_S}" -maxdepth 1 -type f | wc -l) files in ${OUT_S}"
