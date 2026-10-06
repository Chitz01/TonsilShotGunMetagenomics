#!/usr/bin/env bash
# =============================================================================
# 01_metadata/run_metadata.sh — find the tonsil samples in PRJNA48479
#
# Downloads METADATA ONLY (a few MB) from NCBI and ENA, then selects runs whose
# body-site metadata matches TARGET_BODY_SITE_REGEX. No sequencing data.
#
# Light and network-bound: run it directly on the COSMOS login node
# (compute nodes may not have internet access):
#     bash workflow/01_metadata/run_metadata.sh            # fetch + select
#     bash workflow/01_metadata/run_metadata.sh --select   # re-select only
#
# Afterwards, READ ${METADATA_SELECTED_DIR}/<label>_selection_report.txt and
# confirm the matched values are what you want before downloading anything.
# =============================================================================
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export PIPELINE_DIR
source "${PIPELINE_DIR}/config.sh"

HERE="${PIPELINE_DIR}/workflow/01_metadata"
command -v python3 >/dev/null || die "python3 not found (try: module load a Python module, or activate any conda env)"

make_out_dirs "${METADATA_RAW_DIR}" "${METADATA_SELECTED_DIR}"

RUNINFO="${METADATA_RAW_DIR}/${BIOPROJECT}_sra_runinfo.csv"
BIOSAMPLE="${METADATA_RAW_DIR}/${BIOPROJECT}_biosample_attributes.tsv"
ENA="${METADATA_RAW_DIR}/${BIOPROJECT}_ena_read_run.tsv"

if [[ "${1:-}" != "--select" ]]; then
    EMAIL_ARG=()
    is_placeholder "${SLURM_EMAIL}" || EMAIL_ARG=(--email "${SLURM_EMAIL}")
    log "Fetching metadata for ${BIOPROJECT} -> ${METADATA_RAW_DIR}"
    python3 "${HERE}/fetch_metadata.py" \
        --bioproject "${BIOPROJECT}" \
        --outdir "${METADATA_RAW_DIR}" \
        "${EMAIL_ARG[@]}"
fi

[[ -s "${RUNINFO}" ]] || die "Missing ${RUNINFO}; run without --select first"

ONE_ARG=()
[[ "${ONE_SAMPLE_PER_SUBJECT}" == "true" ]] && ONE_ARG=(--one-per-subject)

log "Selecting runs matching /${TARGET_BODY_SITE_REGEX}/ (${LIBRARY_STRATEGY_KEEP})"
python3 "${HERE}/select_samples.py" \
    --runinfo   "${RUNINFO}" \
    --biosample "${BIOSAMPLE}" \
    --ena       "${ENA}" \
    --regex     "${TARGET_BODY_SITE_REGEX}" \
    --label     "${TARGET_LABEL}" \
    --strategy  "${LIBRARY_STRATEGY_KEEP}" \
    --outdir    "${METADATA_SELECTED_DIR}" \
    "${ONE_ARG[@]}"

log "Done. Review before downloading:"
log "  ${METADATA_SELECTED_DIR}/${TARGET_LABEL}_selection_report.txt"
log "  ${METADATA_SELECTED_DIR}/body_site_summary.tsv"
log "Then copy the small tables into the git repo's metadata/ folder for version control."
