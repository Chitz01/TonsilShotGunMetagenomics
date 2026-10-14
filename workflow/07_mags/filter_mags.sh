#!/bin/bash
# =============================================================================
# 07_mags/filter_mags.sh — keep MAGs that pass the quality thresholds
#
# Light; run on the login node from the repository folder:
#     bash workflow/07_mags/filter_mags.sh
# Input : nf-core/mag DAS Tool refined bins + CheckM summary
#         ${NFMAG_OUTDIR}/GenomeBinning/DASTool/bins/*.fa
#         ${NFMAG_OUTDIR}/GenomeBinning/QC/checkm_summary.tsv
# Output: ${MAG_FILTERED_DIR}/*.fa            MAGs with completeness >= 50 %
#                                             and contamination < 10 %
#         ${MAG_DIR}/mag_quality_summary.tsv  all refined bins + pass/fail
# Rule (supervisor, Oct 2026): completeness > 50 %, contamination <= 5 %,
# and DAS Tool 'unbinned' contig sets are never kept.
# =============================================================================
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export PIPELINE_DIR
source "${PIPELINE_DIR}/config.sh"

QC="${NFMAG_OUTDIR}/GenomeBinning/QC/checkm_summary.tsv"
BIN_DIR="${NFMAG_OUTDIR}/GenomeBinning/DASTool/bins"
[[ -s "${QC}" ]]     || die "Missing ${QC} (has nf-core/mag finished?)"
[[ -d "${BIN_DIR}" ]] || die "Missing ${BIN_DIR} (DAS Tool produced no bins?)"
make_out_dirs "${MAG_DIR}" "${MAG_FILTERED_DIR}"

SUMMARY="${MAG_DIR}/mag_quality_summary.tsv"
# Look up columns by name: "Bin Id", "Completeness", "Contamination".
awk -F'\t' -v c="${MAG_MIN_COMPLETENESS}" -v x="${MAG_MAX_CONTAMINATION}" '
    NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i
              if (!("Bin Id" in h) || !("Completeness" in h) || !("Contamination" in h)) {
                  print "unexpected CheckM header" > "/dev/stderr"; exit 1 }
              print "bin\tcompleteness\tcontamination\tpass"; next }
    { pass = ($h["Completeness"] > c && $h["Contamination"] <= x && $h["Bin Id"] !~ /unbinned/) ? "yes" : "no"
      print $h["Bin Id"] "\t" $h["Completeness"] "\t" $h["Contamination"] "\t" pass }' \
    "${QC}" > "${SUMMARY}"

rm -f "${MAG_FILTERED_DIR}"/*.fa
n=0; nf=0
while IFS=$'\t' read -r bin comp cont pass; do
    [[ "${bin}" == "bin" || "${pass}" != "yes" ]] && continue
    if [[ -f "${BIN_DIR}/${bin}.fa" ]]; then
        cp "${BIN_DIR}/${bin}.fa" "${MAG_FILTERED_DIR}/"
    elif [[ -f "${BIN_DIR}/${bin}.fa.gz" ]]; then
        gzip -cd "${BIN_DIR}/${bin}.fa.gz" > "${MAG_FILTERED_DIR}/${bin}.fa"
    else
        log "WARNING: passing bin ${bin} not found in ${BIN_DIR}"; nf=$((nf+1)); continue
    fi
    n=$((n+1))
done < "${SUMMARY}"

total=$(( $(wc -l < "${SUMMARY}") - 1 ))
log "${n} of ${total} CheckM-assessed bins pass (completeness > ${MAG_MIN_COMPLETENESS}%, contamination <= ${MAG_MAX_CONTAMINATION}%, unbinned sets excluded)"
(( nf == 0 )) || log "WARNING: ${nf} passing bins could not be located"
log "MAGs -> ${MAG_FILTERED_DIR}   summary -> ${SUMMARY}"
