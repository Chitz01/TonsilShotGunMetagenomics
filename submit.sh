#!/usr/bin/env bash
# =============================================================================
# submit.sh — convenience wrapper around sbatch, mainly for ARRAY jobs
#
# Account, e-mail, time and log paths are in each script's #SBATCH header
# (CRASH / supervisor style), so single jobs can also be run with plain sbatch
# from this folder:   sbatch workflow/04_taxonomy/merge_metaphlan.sh
#
# What this wrapper adds:
#   --array-runs       one array task per run in SELECTED_RUNS_FILE
#   --array-samples    one array task per sample in SELECTED_SAMPLES_FILE
#   --array-file F     one array task per line of file F
#   --sbatch "OPTS"    extra sbatch options, e.g. --sbatch "--mem=200G"
#   --dry-run          print the sbatch command instead of submitting (safe)
# The array size is counted from the list, so you never type it by hand.
#
# Examples (run from the repository folder):
#   ./submit.sh --dry-run --array-runs workflow/02_download/download_runs.sh
#   ./submit.sh --array-samples        workflow/03_preprocessing/kneaddata.sh
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR
# shellcheck source=config.sh
source "${PIPELINE_DIR}/config.sh"
cd "${PIPELINE_DIR}"     # logs/ in the #SBATCH lines is relative to this folder

usage() { sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

DRY_RUN=false
ARRAY_FILE=""
EXTRA=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)       DRY_RUN=true; shift ;;
        --array-runs)    ARRAY_FILE="${SELECTED_RUNS_FILE}"; shift ;;
        --array-samples) ARRAY_FILE="${SELECTED_SAMPLES_FILE}"; shift ;;
        --array-file)    ARRAY_FILE="$2"; shift 2 ;;
        --sbatch)        read -r -a more <<< "$2"; EXTRA+=("${more[@]}"); shift 2 ;;
        -h|--help)       usage; exit 0 ;;
        --)              shift; break ;;
        -*)              die "Unknown option: $1 (see --help)" ;;
        *)               break ;;
    esac
done
[[ $# -ge 1 ]] || { usage; exit 1; }
SCRIPT="$1"; shift
[[ -f "$SCRIPT" ]] || die "Script not found: $SCRIPT (run from the repository folder)"

ARGS=(--export="ALL,PIPELINE_DIR=${PIPELINE_DIR}")
if [[ -n "$ARRAY_FILE" ]]; then
    [[ -s "$ARRAY_FILE" ]] || die "Array list missing or empty: $ARRAY_FILE (run 01_metadata first)"
    N="$(count_items "$ARRAY_FILE")"
    (( N > 0 )) || die "No entries in $ARRAY_FILE"
    ARGS+=(--array="1-${N}%${SLURM_MAX_PARALLEL}")
    log "Array of ${N} tasks from ${ARRAY_FILE}"
fi

CMD=(sbatch "${ARGS[@]}" "${EXTRA[@]}" "$SCRIPT" "$@")
if $DRY_RUN; then
    printf 'DRY RUN — would execute:\n  '; printf '%q ' "${CMD[@]}"; printf '\n'
    exit 0
fi
mkdir -p logs
"${CMD[@]}"
