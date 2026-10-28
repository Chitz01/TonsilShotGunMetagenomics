#!/bin/bash
# =============================================================================
# 03_preprocessing/cleanup_raw.sh — delete raw FASTQ files that are no longer needed
#
# A sample's raw FASTQ files (all its runs) are deleted ONLY if
#   1. KneadData finished for that sample (marker in KNEADDATA_DIR/.done), and
#   2. its cleaned reads exist and are not empty in CLEAN_FASTQ_DIR.
# Everything downstream (MetaPhlAn, nf-core/mag) uses the cleaned reads only.
# Raw data are public and can be downloaded again from the run accessions in
# metadata/ (remove RAW_FASTQ_DIR/.done/<RUN>.done first so 02_download fetches it).
#
# Every deleted file is recorded in RAW_FASTQ_DIR/deleted_raw_fastq.tsv.
#
# Usage (from the repository folder, after `export DATASET=<name>`):
#   bash workflow/03_preprocessing/cleanup_raw.sh --dry-run     # show only
#   bash workflow/03_preprocessing/cleanup_raw.sh               # all finished samples
#   bash workflow/03_preprocessing/cleanup_raw.sh SRS013947     # one sample
# kneaddata.sh calls it automatically per sample when
# DELETE_RAW_AFTER_KNEADDATA=true (config.sh).
# =============================================================================
set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
export PIPELINE_DIR
source "${PIPELINE_DIR}/config.sh"

DRY=false
if [[ "${1:-}" == "--dry-run" ]]; then DRY=true; shift; fi
if (( $# > 0 )); then
    SAMPLES=("$@")
else
    [[ -s "${SELECTED_SAMPLES_FILE}" ]] || die "Missing ${SELECTED_SAMPLES_FILE}"
    mapfile -t SAMPLES < <(awk 'NF && $1 !~ /^#/ {print $1}' "${SELECTED_SAMPLES_FILE}")
fi
assert_writable_target "${RAW_FASTQ_DIR}"
log "Dataset: ${DATASET} -> ${RAW_FASTQ_DIR}$($DRY && echo '   (DRY RUN: nothing is deleted)')"

MANIFEST="${RAW_FASTQ_DIR}/deleted_raw_fastq.tsv"
if ! $DRY && [[ ! -f "${MANIFEST}" ]]; then
    printf 'date\tdataset\tsample\trun\tfile\tbytes\n' > "${MANIFEST}"
fi

n_files=0; bytes=0; skipped=()
for S in "${SAMPLES[@]}"; do
    if ! is_done "${KNEADDATA_DIR}" "${S}"; then
        skipped+=("${S}: KneadData not finished"); continue
    fi
    if [[ -s "${CLEAN_FASTQ_DIR}/${S}_paired_1.fastq.gz" && -s "${CLEAN_FASTQ_DIR}/${S}_paired_2.fastq.gz" ]] \
       || [[ -s "${CLEAN_FASTQ_DIR}/${S}.fastq.gz" ]]; then
        :
    else
        skipped+=("${S}: cleaned reads missing"); continue
    fi
    mapfile -t RUNS < <(runs_for_sample "${S}")
    for RUN in "${RUNS[@]}"; do
        for f in "${RAW_FASTQ_DIR}/${RUN}_1.fastq.gz" "${RAW_FASTQ_DIR}/${RUN}_2.fastq.gz" "${RAW_FASTQ_DIR}/${RUN}.fastq.gz"; do
            [[ -f "${f}" ]] || continue
            size=$(stat -c %s "${f}")
            if $DRY; then
                printf '  would delete %s (%s MB)\n' "$(basename "${f}")" $(( size / 1048576 )) >&2
            else
                rm -f -- "${f}"
                printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(date '+%F %T')" "${DATASET}" "${S}" "${RUN}" "$(basename "${f}")" "${size}" >> "${MANIFEST}"
            fi
            n_files=$((n_files + 1)); bytes=$((bytes + size))
        done
    done
done

gb=$(awk -v b="${bytes}" 'BEGIN {printf "%.1f", b / 1e9}')
if $DRY; then
    log "DRY RUN: would delete ${n_files} raw FASTQ files (${gb} GB) from ${#SAMPLES[@]} samples"
else
    log "Deleted ${n_files} raw FASTQ files (${gb} GB freed); record: ${MANIFEST}"
fi
if (( ${#skipped[@]} > 0 )); then
    log "Kept raw files of ${#skipped[@]} sample(s):"
    printf '  %s\n' "${skipped[@]}" >&2
fi
