#!/bin/bash
#SBATCH -A lu2026-2-117
#SBATCH -J download
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 06:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#SBATCH --mail-user=ch2068va-s@student.lu.se
#SBATCH --mail-type=END,FAIL
# =============================================================================
# 02_download/download_runs.sh — download the SELECTED runs only
#
# Downloads one run per array task from SELECTED_RUNS_FILE (made by 01_metadata).
# Output: RAW_FASTQ_DIR/<RUN>_1.fastq.gz, <RUN>_2.fastq.gz (paired) or
#         <RUN>.fastq.gz (single-end / orphan reads).
# Method from config.sh DOWNLOAD_METHOD:
#   sra : prefetch + vdb-validate + fasterq-dump (conda env ENV_SRA)
#   ena : direct FASTQ from ENA with MD5 check (needs only curl)
#
# Usage:
#   ./submit.sh --array-runs workflow/02_download/download_runs.sh   # all runs
#   bash workflow/02_download/download_runs.sh SRR0000000            # one run
#
# Check first whether COSMOS compute nodes have outbound internet. If they do
# not, use DOWNLOAD_METHOD=ena and run single runs from the login node, or ask
# LUNARC support for the recommended transfer route.
# Finished runs are skipped on resubmission (marker in RAW_FASTQ_DIR/.done/).
# =============================================================================
set -euo pipefail
# Find config.sh: submit from the repository folder (SLURM_SUBMIT_DIR).
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${PIPELINE_DIR}/config.sh" ]] || { echo "config.sh not found in ${PIPELINE_DIR}: submit from the repository folder" >&2; exit 1; }
source "${PIPELINE_DIR}/config.sh"

RUN="$(resolve_item "${SELECTED_RUNS_FILE}" "${1:-}")"
T="$(threads)"
make_out_dirs "${RAW_FASTQ_DIR}" "${SRA_DIR}"

if is_done "${RAW_FASTQ_DIR}" "${RUN}"; then
    log "[skip] ${RUN} already downloaded"
    exit 0
fi

make_scratch download "${RUN}"
trap 'rm -rf "${WORK}"' EXIT
log "Downloading ${RUN} with method '${DOWNLOAD_METHOD}' (${T} threads) on $(hostname)"

case "${DOWNLOAD_METHOD}" in
# -----------------------------------------------------------------------------
sra)
    load_conda_env "${ENV_SRA}"
    # prefetch stores RUN/RUN.sra inside SRA_DIR; fasterq-dump finds it when
    # run from that directory.
    prefetch "${RUN}" --output-directory "${SRA_DIR}" --max-size 200G
    ( cd "${SRA_DIR}" && vdb-validate "${RUN}" )
    ( cd "${SRA_DIR}" && fasterq-dump "${RUN}" \
        --split-3 \
        --threads "${T}" \
        --temp "${WORK}/tmp" \
        --outdir "${WORK}/out" )
    shopt -s nullglob
    fq=( "${WORK}/out/${RUN}"*.fastq )
    (( ${#fq[@]} > 0 )) || die "fasterq-dump produced no FASTQ for ${RUN}"
    pigz -p "${T}" "${fq[@]}"
    ;;
# -----------------------------------------------------------------------------
ena)
    [[ -s "${SELECTED_RUNS_TABLE}" ]] || die "Missing ${SELECTED_RUNS_TABLE} (run 01_metadata)"
    # Columns of <label>_runs_full.tsv: find ena_fastq_ftp / ena_fastq_md5 by name.
    read -r URLS MD5S < <(awk -F'\t' -v r="${RUN}" '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        $col["run"] == r { print $col["ena_fastq_ftp"], $col["ena_fastq_md5"] }' \
        "${SELECTED_RUNS_TABLE}") || true
    [[ -n "${URLS:-}" ]] || die "No ENA FASTQ URLs for ${RUN}; use DOWNLOAD_METHOD=sra"
    IFS=';' read -r -a URL_ARR <<< "${URLS}"
    IFS=';' read -r -a MD5_ARR <<< "${MD5S:-}"
    for i in "${!URL_ARR[@]}"; do
        url="${URL_ARR[$i]}"; f="${WORK}/out/$(basename "${url}")"
        log "  fetching ${url}"
        curl -fL --retry 5 --retry-delay 20 -o "${f}" "https://${url#ftp://}"
        if [[ -n "${MD5_ARR[$i]:-}" ]]; then
            echo "${MD5_ARR[$i]}  ${f}" | md5sum -c - >/dev/null \
                || die "MD5 mismatch for $(basename "${f}")"
            log "  md5 ok: $(basename "${f}")"
        fi
    done
    ;;
# -----------------------------------------------------------------------------
*) die "Unknown DOWNLOAD_METHOD='${DOWNLOAD_METHOD}' (use sra or ena)" ;;
esac

# ---- Sanity check and move to project storage -------------------------------
shopt -s nullglob
out=( "${WORK}/out/${RUN}"*.fastq.gz )
(( ${#out[@]} > 0 )) || die "No FASTQ files for ${RUN}"
for f in "${out[@]}"; do
    gzip -t "${f}" || die "Corrupt gzip: ${f}"
done
if [[ -f "${WORK}/out/${RUN}_1.fastq.gz" && ! -f "${WORK}/out/${RUN}_2.fastq.gz" ]]; then
    die "${RUN}: found _1 but not _2"
fi
mv "${out[@]}" "${RAW_FASTQ_DIR}/"
ls -lh "${RAW_FASTQ_DIR}/${RUN}"*.fastq.gz >&2

if [[ "${DOWNLOAD_METHOD}" == "sra" && "${KEEP_SRA_FILES}" != "true" ]]; then
    rm -rf "${SRA_DIR:?}/${RUN}"
    log "Removed cached .sra for ${RUN} (KEEP_SRA_FILES=false)"
fi

mark_done "${RAW_FASTQ_DIR}" "${RUN}"
log "Finished ${RUN}"
