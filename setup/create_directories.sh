#!/usr/bin/env bash
# =============================================================================
# setup/create_directories.sh — create the output folder tree on COSMOS
#
# Creates PROJECT_DIR and all stage folders under ANALYSIS with `mkdir -p`
# (existing folders and files are left untouched). Never writes to DATA.
# Safe to run on the login node:  bash setup/create_directories.sh
# =============================================================================
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PIPELINE_DIR
source "${PIPELINE_DIR}/config.sh"

[[ -d "${ANALYSIS_ROOT}" ]] || die "ANALYSIS_ROOT not found: ${ANALYSIS_ROOT} (are you on COSMOS?)"
[[ -w "${ANALYSIS_ROOT}" ]] || die "ANALYSIS_ROOT is not writable for you: ${ANALYSIS_ROOT}"

DIRS=(
    "${PROJECT_DIR}"
    "${METADATA_RAW_DIR}" "${METADATA_SELECTED_DIR}"
    "${SRA_DIR}" "${RAW_FASTQ_DIR}"
    "${CLEAN_FASTQ_DIR}" "${KNEADDATA_LOG_DIR}" "${KNEADDATA_QC_DIR}"
    "${METAPHLAN_DIR}"
    "${ASSEMBLY_DIR}" "${NXF_CACHE_DIR}"
    "${GENES_DIR}" "${MAG_DIR}" "${MAG_FILTERED_DIR}"
    "${ENZYME_DIR}" "${INTEGRATION_DIR}"
)
mkdir -p "${PIPELINE_DIR}/logs"     # SLURM .out/.err files (see #SBATCH -o lines)

make_out_dirs "${DIRS[@]}"
for d in "${DIRS[@]}"; do log "ok  ${d}"; done
log "Output tree ready under ${PROJECT_DIR}"
