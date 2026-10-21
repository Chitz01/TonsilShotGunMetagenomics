#!/usr/bin/env bash
# =============================================================================
# 01_metadata/run_metadata.sh — find the samples of the chosen DATASET
#
# Downloads METADATA ONLY (a few MB) from NCBI and ENA, then selects runs whose
# body-site metadata matches TARGET_BODY_SITE_REGEX (datasets/<name>.sh), with
# the optional filters SUBJECTS_FILE, PLATFORM_REGEX, MIN_BASES and
# ONE_SAMPLE_PER_SUBJECT (true = per person, site = per person and site,
# false = keep all). No sequencing data are downloaded.
#
# Run on the login node, from the repository folder:
#     export DATASET=<name>
#     bash workflow/01_metadata/run_metadata.sh            # fetch + select
#     bash workflow/01_metadata/run_metadata.sh --select   # re-select only
# Then READ <label>_selection_report.txt before downloading anything.
# =============================================================================
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export PIPELINE_DIR
source "${PIPELINE_DIR}/config.sh"
log "Dataset: ${DATASET} (${BIOPROJECT}) -> ${PROJECT_DIR}"

HERE="${PIPELINE_DIR}/workflow/01_metadata"
command -v python3 >/dev/null || die "python3 not found"
make_out_dirs "${METADATA_RAW_DIR}" "${METADATA_SELECTED_DIR}"

RUNINFO="${METADATA_RAW_DIR}/${BIOPROJECT}_sra_runinfo.csv"
BIOSAMPLE="${METADATA_RAW_DIR}/${BIOPROJECT}_biosample_attributes.tsv"
ENA="${METADATA_RAW_DIR}/${BIOPROJECT}_ena_read_run.tsv"

if [[ "${1:-}" != "--select" ]]; then
    EMAIL_ARG=()
    is_placeholder "${SLURM_EMAIL}" || EMAIL_ARG=(--email "${SLURM_EMAIL}")
    log "Fetching metadata for ${BIOPROJECT} -> ${METADATA_RAW_DIR}"
    python3 "${HERE}/fetch_metadata.py" --bioproject "${BIOPROJECT}" \
        --outdir "${METADATA_RAW_DIR}" "${EMAIL_ARG[@]}"
fi
[[ -s "${RUNINFO}" ]] || die "Missing ${RUNINFO}; run without --select first"

# Optional filters (defined only in some datasets/<name>.sh files)
case "${ONE_SAMPLE_PER_SUBJECT:-false}" in
    true)  ONE_MODE="subject" ;;
    site)  ONE_MODE="site" ;;
    *)     ONE_MODE="none" ;;
esac
EXTRA=(--one-per-subject "${ONE_MODE}")
if [[ -n "${SUBJECTS_FILE:-}" ]]; then
    [[ -s "${SUBJECTS_FILE}" ]] || die "SUBJECTS_FILE not found or empty: ${SUBJECTS_FILE}"
    EXTRA+=(--subjects "${SUBJECTS_FILE}")
fi
[[ -n "${PLATFORM_REGEX:-}" ]] && EXTRA+=(--platform "${PLATFORM_REGEX}")
if [[ -n "${MATCH_TO_FILE:-}" ]]; then
    [[ -s "${MATCH_TO_FILE}" ]] || die "MATCH_TO_FILE not found or empty: ${MATCH_TO_FILE}"
    EXTRA+=(--match-to "${MATCH_TO_FILE}")
fi
[[ -n "${MIN_BASES:-}" ]] && EXTRA+=(--min-bases "${MIN_BASES}")

log "Selecting runs matching /${TARGET_BODY_SITE_REGEX}/ (${LIBRARY_STRATEGY_KEEP}; ${EXTRA[*]})"
python3 "${HERE}/select_samples.py" \
    --runinfo   "${RUNINFO}" \
    --biosample "${BIOSAMPLE}" \
    --ena       "${ENA}" \
    --regex     "${TARGET_BODY_SITE_REGEX}" \
    --label     "${TARGET_LABEL}" \
    --strategy  "${LIBRARY_STRATEGY_KEEP}" \
    --outdir    "${METADATA_SELECTED_DIR}" \
    "${EXTRA[@]}"

log "Done. Review before downloading:"
log "  ${METADATA_SELECTED_DIR}/${TARGET_LABEL}_selection_report.txt"
