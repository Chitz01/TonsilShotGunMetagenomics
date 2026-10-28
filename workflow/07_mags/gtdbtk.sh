#!/bin/bash
#SBATCH -A lu2026-2-117
#SBATCH -J gtdbtk
#SBATCH -N 1
#SBATCH --ntasks-per-node=25
#SBATCH -t 160:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=ch2068va-s@student.lu.se
#SBATCH --mail-type=END,FAIL
# =============================================================================
# 07_mags/gtdbtk.sh — GTDB taxonomy of the quality-filtered MAGs
#
# Follows the supervisor's GTDB script (GTDB-Tk classify_wf, GTDB release 226,
# --skip_ani_screen, 25 cores). Run from the repository folder:
#     sbatch workflow/07_mags/gtdbtk.sh
# Input : ${MAG_FILTERED_DIR}/*.fa   (from filter_mags.sh)
# Output: ${MAG_DIR}/gtdbtk/          incl. gtdbtk.bac120.summary.tsv
# =============================================================================
set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${PIPELINE_DIR}/config.sh" ]] || { echo "config.sh not found in ${PIPELINE_DIR}: submit from the repository folder" >&2; exit 1; }
source "${PIPELINE_DIR}/config.sh"

T="$(threads)"
require_vars ENV_GTDBTK GTDBTK_DATA_PATH
require_paths GTDBTK_DATA_PATH
OUT="${MAG_DIR}/gtdbtk"
make_out_dirs "${OUT}"
# Short temporary path. GTDB-Tk's helper processes talk through socket files
# whose full path must stay below ~108 characters; under the long project path
# this fails with "OSError: AF_UNIX path too long". Node-local scratch
# ($SNIC_TMP) has a short path; /tmp is the fallback.
if [[ "${SNIC_TMP:-}" == /* && -d "${SNIC_TMP}" ]]; then GT_BASE="${SNIC_TMP}"; else GT_BASE=/tmp; fi
GT_TMP="${GT_BASE}/gtdbtk_${SLURM_JOB_ID:-$$}"
mkdir -p "${GT_TMP}"
export TMPDIR="${GT_TMP}"
trap 'rm -rf "${GT_TMP}"' EXIT
log "Temporary folder: ${GT_TMP}"

shopt -s nullglob
MAGS=( "${MAG_FILTERED_DIR}"/*.fa )
(( ${#MAGS[@]} > 0 )) || die "No MAGs in ${MAG_FILTERED_DIR} (run filter_mags.sh)"
log "Classifying ${#MAGS[@]} MAGs with ${T} cores"

load_conda_env "${ENV_GTDBTK}"
export GTDBTK_DATA_PATH
gtdbtk --version >&2 || true

# GTDB-Tk >= 2.4 requires either --mash_db or --skip_ani_screen.
ANI_OPT=()
gtdbtk classify_wf -h 2>&1 | grep -q -- '--skip_ani_screen' && ANI_OPT=(--skip_ani_screen)

gtdbtk classify_wf \
    --genome_dir "${MAG_FILTERED_DIR}" \
    --out_dir "${OUT}" \
    --extension fa \
    --cpus "${T}" \
    --tmpdir "${GT_TMP}" \
    "${ANI_OPT[@]}"

ls "${OUT}"/gtdbtk.*.summary.tsv >&2 || log "WARNING: no summary table produced"
log "GTDB-Tk results in ${OUT}"
