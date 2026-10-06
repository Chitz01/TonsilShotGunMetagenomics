#!/bin/bash
# =============================================================================
# 05_assembly/make_samplesheet.sh — nf-core/mag input sheet from KneadData output
#
# Light; run on the login node from the repository folder:
#     bash workflow/05_assembly/make_samplesheet.sh
# Writes NFMAG_SAMPLESHEET (CSV) with one row per selected sample:
#     sample,group,short_reads_1,short_reads_2,long_reads
# All samples get group 0, so with NFMAG_BINNING_MAP_MODE=group every sample's
# reads are mapped to every assembly to compute co-abundance for binning.
# Only the repaired, properly paired reads are used (nf-core/mag accepts one
# pair of files per sample).
# =============================================================================
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export PIPELINE_DIR
source "${PIPELINE_DIR}/config.sh"

[[ -s "${SELECTED_SAMPLES_FILE}" ]] || die "Missing ${SELECTED_SAMPLES_FILE} (run 01_metadata)"
make_out_dirs "${ASSEMBLY_DIR}"

TMP="$(mktemp)"
echo "sample,group,short_reads_1,short_reads_2,long_reads" > "${TMP}"
n=0; missing=()
while read -r S; do
    [[ -z "$S" || "$S" == \#* ]] && continue
    R1="${CLEAN_FASTQ_DIR}/${S}_paired_1.fastq.gz"
    R2="${CLEAN_FASTQ_DIR}/${S}_paired_2.fastq.gz"
    if [[ -s "$R1" && -s "$R2" ]]; then
        echo "${S},0,${R1},${R2}," >> "${TMP}"
        n=$((n+1))
    else
        missing+=("$S")
    fi
done < "${SELECTED_SAMPLES_FILE}"

if (( ${#missing[@]} > 0 )); then
    rm -f "${TMP}"
    die "Cleaned reads missing for ${#missing[@]} sample(s): ${missing[*]} (finish 03_preprocessing first)"
fi
mv "${TMP}" "${NFMAG_SAMPLESHEET}"
log "Wrote ${n} samples -> ${NFMAG_SAMPLESHEET}"
head -n 3 "${NFMAG_SAMPLESHEET}" >&2
