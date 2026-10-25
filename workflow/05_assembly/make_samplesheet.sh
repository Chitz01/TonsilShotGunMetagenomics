#!/bin/bash
# =============================================================================
# 05_assembly/make_samplesheet.sh — nf-core/mag input sheet from KneadData output
#
# Light; run on the login node from the repository folder:
#     bash workflow/05_assembly/make_samplesheet.sh
# Writes NFMAG_SAMPLESHEET (CSV) with one row per selected sample:
#     sample,group,short_reads_1,short_reads_2,long_reads
# Group: 0 for all samples, or one number per body site if the dataset sets
# SAMPLESHEET_GROUP_BY=site. With NFMAG_BINNING_MAP_MODE=group, each sample's
# reads are mapped to every assembly of its group to compute co-abundance.
# Only the repaired, properly paired reads are used (nf-core/mag accepts one
# pair of files per sample).
# =============================================================================
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export PIPELINE_DIR
source "${PIPELINE_DIR}/config.sh"

[[ -s "${SELECTED_SAMPLES_FILE}" ]] || die "Missing ${SELECTED_SAMPLES_FILE} (run 01_metadata)"
make_out_dirs "${ASSEMBLY_DIR}"

# Groups: 0 for all samples (tonsil), or one number per body site when the
# dataset sets SAMPLESHEET_GROUP_BY=site (binning then uses co-abundance
# within each site). The number-to-site key is written next to the sheet.
GROUP_KEY="${ASSEMBLY_DIR}/samplesheet_groups.tsv"
sample_group() { echo 0; }
if [[ "${SAMPLESHEET_GROUP_BY:-none}" == "site" ]]; then
    [[ -s "${SELECTED_RUNS_TABLE}" ]] || die "Missing ${SELECTED_RUNS_TABLE} (run 01_metadata)"
    awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} {print $h["sample_id"] "\t" $h["body_site"]}' \
        "${SELECTED_RUNS_TABLE}" | sort -u > "${GROUP_KEY}.tmp"
    cut -f2 "${GROUP_KEY}.tmp" | sort -u | awk 'BEGIN {OFS = "\t"} {print $0, NR}' > "${GROUP_KEY}.sites"
    { printf 'sample_id\tgroup\tbody_site\n'
      awk -F'\t' 'BEGIN {OFS = "\t"} NR == FNR {g[$1] = $2; next} {print $1, g[$2], $2}' "${GROUP_KEY}.sites" "${GROUP_KEY}.tmp"
    } > "${GROUP_KEY}"
    rm -f "${GROUP_KEY}.tmp" "${GROUP_KEY}.sites"
    sample_group() { awk -F'\t' -v s="$1" 'NR > 1 && $1 == s {print $2; exit}' "${GROUP_KEY}"; }
    log "Groups by body site (key -> ${GROUP_KEY}):"
    awk -F'\t' 'NR > 1 {n[$2 "\t" $3]++} END {for (k in n) print "  group " k ": " n[k] " samples"}' "${GROUP_KEY}" | sort -n -k2 >&2
fi

TMP="$(mktemp)"
echo "sample,group,short_reads_1,short_reads_2,long_reads" > "${TMP}"
n=0; missing=(); excluded=()
if [[ -n "${MAG_EXCLUDE_FILE:-}" ]]; then
    [[ -s "${MAG_EXCLUDE_FILE}" ]] || die "MAG_EXCLUDE_FILE not found or empty: ${MAG_EXCLUDE_FILE}"
fi
while read -r S; do
    [[ -z "$S" || "$S" == \#* ]] && continue
    # samples listed in MAG_EXCLUDE_FILE (first column) stay out of the MAG branch
    if [[ -n "${MAG_EXCLUDE_FILE:-}" ]] && awk -F'\t' -v s="$S" 'NR > 1 && $1 == s {f = 1} END {exit !f}' "${MAG_EXCLUDE_FILE}"; then
        excluded+=("$S"); continue
    fi
    R1="${CLEAN_FASTQ_DIR}/${S}_paired_1.fastq.gz"
    R2="${CLEAN_FASTQ_DIR}/${S}_paired_2.fastq.gz"
    if [[ -s "$R1" && -s "$R2" ]]; then
        g="$(sample_group "${S}")"
        [[ -n "${g}" ]] || die "No group found for ${S} in ${GROUP_KEY}"
        echo "${S},${g},${R1},${R2}," >> "${TMP}"
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
if (( ${#excluded[@]} > 0 )); then
    log "Left out of the MAG branch (listed in ${MAG_EXCLUDE_FILE}): ${excluded[*]}"
fi
head -n 3 "${NFMAG_SAMPLESHEET}" >&2
