#!/bin/bash
#SBATCH -A lu2026-2-117
#SBATCH -J kneaddata_summary
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH -t 01:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=ch2068va-s@student.lu.se
#SBATCH --mail-type=END,FAIL
# =============================================================================
# 03_preprocessing/kneaddata_summary.sh — read-count table + MultiQC report
#
# Run after all kneaddata.sh tasks finish:
#   ./submit.sh workflow/03_preprocessing/kneaddata_summary.sh
# Outputs in KNEADDATA_DIR:
#   kneaddata_read_counts.tsv  raw / trimmed / host-removed reads per sample
#   multiqc/                   combined FastQC report
# =============================================================================
set -euo pipefail
# Find config.sh: submit from the repository folder (SLURM_SUBMIT_DIR).
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${PIPELINE_DIR}/config.sh" ]] || { echo "config.sh not found in ${PIPELINE_DIR}: submit from the repository folder" >&2; exit 1; }
source "${PIPELINE_DIR}/config.sh"

ls "${KNEADDATA_LOG_DIR}"/*.log >/dev/null 2>&1 || die "No KneadData logs in ${KNEADDATA_LOG_DIR}"
load_conda_env "${ENV_KNEADDATA}"

kneaddata_read_count_table \
    --input "${KNEADDATA_LOG_DIR}" \
    --output "${KNEADDATA_DIR}/kneaddata_read_counts.tsv"
log "Read counts -> ${KNEADDATA_DIR}/kneaddata_read_counts.tsv"

if command -v multiqc >/dev/null; then
    multiqc --force --outdir "${KNEADDATA_DIR}/multiqc" "${KNEADDATA_QC_DIR}"
    log "MultiQC -> ${KNEADDATA_DIR}/multiqc"
fi

# Expected sample count vs. completed
n_exp="$(count_items "${SELECTED_SAMPLES_FILE}")"
n_done="$(find "${KNEADDATA_DIR}/.done" -name '*.done' 2>/dev/null | wc -l)"
log "Samples completed: ${n_done}/${n_exp}"
