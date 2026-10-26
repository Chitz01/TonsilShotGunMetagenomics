#!/bin/bash
#SBATCH -A lu2026-2-117
#SBATCH -J bakta
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH -t 160:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=ch2068va-s@student.lu.se
#SBATCH --mail-type=END,FAIL
# =============================================================================
# 07_mags/bakta.sh — gene/function annotation of each MAG with Bakta
#
# Follows the supervisor's Bakta script (one MAG at a time, --meta, 20 cores,
# MAGs already annotated are skipped, so the job can simply be resubmitted).
#     sbatch workflow/07_mags/bakta.sh
# Input : ${MAG_FILTERED_DIR}/*.fa
# Output: ${MAG_DIR}/bakta/<MAG>/<MAG>.gff3, .faa, .tsv, ...
# The protein files (.faa) are the input for 08_enzyme_discovery later.
# =============================================================================
set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${PIPELINE_DIR}/config.sh" ]] || { echo "config.sh not found in ${PIPELINE_DIR}: submit from the repository folder" >&2; exit 1; }
source "${PIPELINE_DIR}/config.sh"

T="$(threads)"
# Optional splitting: "bakta.sh CHUNK NCHUNKS" annotates only MAGs number
# CHUNK, CHUNK+NCHUNKS, ... (0-based), so several jobs can run side by side,
# each with its own temporary folder. Without arguments: all MAGs, one job.
CHUNK="${1:-0}"; NCHUNKS="${2:-1}"
if ! [[ "${CHUNK}" =~ ^[0-9]+$ && "${NCHUNKS}" =~ ^[0-9]+$ ]] || (( NCHUNKS < 1 || CHUNK >= NCHUNKS )); then
    die "Usage: bakta.sh [CHUNK NCHUNKS], e.g. bakta.sh 0 8 (CHUNK from 0 to NCHUNKS-1)"
fi
require_vars ENV_BAKTA BAKTA_DB
require_paths BAKTA_DB
BAKDIR="${MAG_DIR}/bakta"
make_out_dirs "${BAKDIR}" "${MAG_DIR}/bakta_tmp_${CHUNK}"

shopt -s nullglob
MAGS=( "${MAG_FILTERED_DIR}"/*.fa )
(( ${#MAGS[@]} > 0 )) || die "No MAGs in ${MAG_FILTERED_DIR} (run filter_mags.sh)"

load_conda_env "${ENV_BAKTA}"
bakta --version >&2 || true

done_n=0; new_n=0; i=-1
for mag in "${MAGS[@]}"; do
    i=$((i + 1))
    (( i % NCHUNKS == CHUNK )) || continue
    f="$(basename "${mag}" .fa)"
    if [[ -s "${BAKDIR}/${f}/${f}.gff3" ]]; then
        done_n=$((done_n+1)); continue
    fi
    log "Annotating ${f}"
    bakta "${mag}" \
        --db "${BAKTA_DB}" \
        --threads "${T}" \
        --prefix "${f}" \
        --output "${BAKDIR}/${f}" \
        --meta \
        --force \
        --tmp-dir "${MAG_DIR}/bakta_tmp_${CHUNK}"
    new_n=$((new_n+1))
done

rm -rf "${MAG_DIR}/bakta_tmp_${CHUNK}"
log "Bakta chunk ${CHUNK}/${NCHUNKS}: ${new_n} annotated now, ${done_n} already done (${#MAGS[@]} MAGs in total)"
