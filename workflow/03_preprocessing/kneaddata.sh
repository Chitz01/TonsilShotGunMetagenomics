#!/bin/bash
#SBATCH -A lu2026-2-117
#SBATCH -J kneaddata
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH -t 12:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#SBATCH --mail-user=ch2068va-s@student.lu.se
#SBATCH --mail-type=END,FAIL
# =============================================================================
# 03_preprocessing/kneaddata.sh — quality trimming + human host-read removal
#
# One array task per SAMPLE (SELECTED_SAMPLES_FILE). HMP samples often have
# several sequencing runs; all runs of a sample are concatenated first.
#
# Steps: FastQC (before) -> Trimmomatic -> Bowtie2 against the human reference
#        (KNEADDATA_DB) -> FastQC (after) -> BBTools repair.sh, which puts the
#        two read files back in matching order (as in the supervisor's
#        assembly script) so nf-core/mag gets properly paired reads.
# Kept   : CLEAN_FASTQ_DIR/<S>_paired_1/2.fastq.gz, <S>_unmatched_1/2.fastq.gz
#          (or <S>.fastq.gz for single-end), log, FastQC reports.
# Deleted: the human (contaminant) reads stay on node-local scratch and are
#          removed with it — they are never copied to project storage.
#
# Note: HMP deposited reads in SRA after removing human DNA, so expect very
# few host reads here. This step is still useful for trimming, for a uniform
# QC record, and as an independent check of residual human sequence.
#
# Usage:
#   ./submit.sh --array-samples workflow/03_preprocessing/kneaddata.sh
#   ./submit.sh workflow/03_preprocessing/kneaddata.sh SRS000000   # one sample
# =============================================================================
set -euo pipefail
# Find config.sh: submit from the repository folder (SLURM_SUBMIT_DIR).
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${PIPELINE_DIR}/config.sh" ]] || { echo "config.sh not found in ${PIPELINE_DIR}: submit from the repository folder" >&2; exit 1; }
source "${PIPELINE_DIR}/config.sh"

SAMPLE="$(resolve_item "${SELECTED_SAMPLES_FILE}" "${1:-}")"
T="$(threads)"

require_vars KNEADDATA_DB CONDA_MODULE
require_paths KNEADDATA_DB
make_out_dirs "${CLEAN_FASTQ_DIR}" "${KNEADDATA_LOG_DIR}" "${KNEADDATA_QC_DIR}"

if is_done "${KNEADDATA_DIR}" "${SAMPLE}"; then
    log "[skip] ${SAMPLE} already processed"
    exit 0
fi

# ---- Gather the runs of this sample -------------------------------------------
mapfile -t RUNS < <(runs_for_sample "${SAMPLE}")
(( ${#RUNS[@]} > 0 )) || die "No runs listed for ${SAMPLE} in ${SAMPLE_RUN_MAP}"
log "${SAMPLE}: ${#RUNS[@]} run(s): ${RUNS[*]}"

PAIRED=(); SINGLE=()
for RUN in "${RUNS[@]}"; do
    if [[ -f "${RAW_FASTQ_DIR}/${RUN}_1.fastq.gz" && -f "${RAW_FASTQ_DIR}/${RUN}_2.fastq.gz" ]]; then
        PAIRED+=("${RUN}")
        [[ -f "${RAW_FASTQ_DIR}/${RUN}.fastq.gz" ]] && \
            log "  note: ${RUN}.fastq.gz (orphan reads from --split-3) is not used"
    elif [[ -f "${RAW_FASTQ_DIR}/${RUN}.fastq.gz" ]]; then
        SINGLE+=("${RUN}")
    else
        die "Raw FASTQ for ${RUN} not found in ${RAW_FASTQ_DIR} (run 02_download)"
    fi
done
if (( ${#PAIRED[@]} > 0 && ${#SINGLE[@]} > 0 )); then
    die "${SAMPLE} mixes paired (${PAIRED[*]}) and single-end (${SINGLE[*]}) runs; decide how to handle it manually"
fi

make_scratch kneaddata "${SAMPLE}"
trap 'rm -rf "${WORK}"' EXIT
cd "${WORK}"

if (( ${#PAIRED[@]} > 0 )); then
    MODE="paired"
    for RUN in "${PAIRED[@]}"; do
        cat "${RAW_FASTQ_DIR}/${RUN}_1.fastq.gz" >> "in/${SAMPLE}_1.fastq.gz"
        cat "${RAW_FASTQ_DIR}/${RUN}_2.fastq.gz" >> "in/${SAMPLE}_2.fastq.gz"
    done
    INPUT_ARGS=(--input1 "in/${SAMPLE}_1.fastq.gz" --input2 "in/${SAMPLE}_2.fastq.gz")
else
    MODE="single"
    for RUN in "${SINGLE[@]}"; do
        cat "${RAW_FASTQ_DIR}/${RUN}.fastq.gz" >> "in/${SAMPLE}.fastq.gz"
    done
    INPUT_ARGS=(--unpaired "in/${SAMPLE}.fastq.gz")
fi
log "${SAMPLE}: ${MODE}-end input prepared on $(hostname), ${T} threads"

# ---- KneadData -------------------------------------------------------------------
load_conda_env "${ENV_KNEADDATA}"
TRIM_JAR="$(find "${CONDA_PREFIX}/share" -name 'trimmomatic*.jar' -print -quit)"
[[ -n "${TRIM_JAR}" ]] || die "Trimmomatic jar not found in ${CONDA_PREFIX}/share"
TRIM_DIR="$(dirname "${TRIM_JAR}")"
kneaddata --version 2>&1 | head -n1 >&2 || true

read -r -a EXTRA_OPTS <<< "${KNEADDATA_EXTRA_OPTS}"

kneaddata \
    "${INPUT_ARGS[@]}" \
    --reference-db "${KNEADDATA_DB}" \
    --output "out" \
    --output-prefix "${SAMPLE}" \
    --threads "${T}" \
    --processes 1 \
    --trimmomatic "${TRIM_DIR}" \
    --sequencer-source "${KNEADDATA_SEQUENCER_SOURCE}" \
    --max-memory "${KNEADDATA_MAX_MEMORY}" \
    --run-fastqc-start \
    --run-fastqc-end \
    --remove-intermediate-output \
    --log "out/${SAMPLE}.log" \
    "${EXTRA_OPTS[@]}"

# ---- Re-synchronise read pairs (paired mode) -------------------------------------
if [[ "${MODE}" == "paired" ]]; then
    repair.sh \
        in1="out/${SAMPLE}_paired_1.fastq" in2="out/${SAMPLE}_paired_2.fastq" \
        out1="out/${SAMPLE}_repaired_1.fastq" out2="out/${SAMPLE}_repaired_2.fastq" \
        outs="out/${SAMPLE}_repair_singletons.fastq" \
        repair=t "-Xmx${REPAIR_MEMORY}" 2> "out/${SAMPLE}.repair.txt"
    mv "out/${SAMPLE}_repaired_1.fastq" "out/${SAMPLE}_paired_1.fastq"
    mv "out/${SAMPLE}_repaired_2.fastq" "out/${SAMPLE}_paired_2.fastq"
    # reads whose mate was lost go with the other unpaired reads
    [[ -s "out/${SAMPLE}_repair_singletons.fastq" ]] && \
        cat "out/${SAMPLE}_repair_singletons.fastq" >> "out/${SAMPLE}_unmatched_1.fastq"
    tail -n 5 "out/${SAMPLE}.repair.txt" >&2 || true
fi

# ---- Keep only the cleaned reads (host reads are discarded) ----------------------
if [[ "${MODE}" == "paired" ]]; then
    KEEP=( "out/${SAMPLE}_paired_1.fastq" "out/${SAMPLE}_paired_2.fastq"
           "out/${SAMPLE}_unmatched_1.fastq" "out/${SAMPLE}_unmatched_2.fastq" )
    for f in "${KEEP[@]:0:2}"; do [[ -s "$f" ]] || die "Expected output missing/empty: $f"; done
    for f in "${KEEP[@]:2:2}"; do [[ -e "$f" ]] || : > "$f"; done   # may legitimately be empty
else
    KEEP=( "out/${SAMPLE}.fastq" )
    [[ -s "${KEEP[0]}" ]] || die "Expected output missing/empty: ${KEEP[0]}"
fi

pigz -p "${T}" "${KEEP[@]}"
for f in "${KEEP[@]}"; do cp "${f}.gz" "${CLEAN_FASTQ_DIR}/"; done
cp "out/${SAMPLE}.log" "${KNEADDATA_LOG_DIR}/"
# (.txt, not .log, so kneaddata_read_count_table ignores it)
[[ -f "out/${SAMPLE}.repair.txt" ]] && cp "out/${SAMPLE}.repair.txt" "${KNEADDATA_LOG_DIR}/"
mkdir -p "${KNEADDATA_QC_DIR}/${SAMPLE}"
find "out" -type f \( -name '*fastqc.html' -o -name '*fastqc.zip' \) \
    -exec cp {} "${KNEADDATA_QC_DIR}/${SAMPLE}/" \;

ls -lh "${CLEAN_FASTQ_DIR}/${SAMPLE}"*.fastq.gz >&2
mark_done "${KNEADDATA_DIR}" "${SAMPLE}"

# Delete this sample's raw FASTQ now that its cleaned reads are saved
# (supervisor's request; set DELETE_RAW_AFTER_KNEADDATA=false in config.sh to keep them).
if [[ "${DELETE_RAW_AFTER_KNEADDATA:-false}" == "true" ]]; then
    bash "${PIPELINE_DIR}/workflow/03_preprocessing/cleanup_raw.sh" "${SAMPLE}" \
        || log "WARNING: raw FASTQ cleanup failed for ${SAMPLE} (results are fine; run cleanup_raw.sh later)"
fi
log "Finished ${SAMPLE} (host reads discarded with scratch)"
