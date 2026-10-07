#!/usr/bin/env bash
# =============================================================================
# config.sh — central configuration for TonsilShotGunMetagenomics
#             (BINP52 tonsil/oral microbiome, tire-particle enzymes; LUNARC COSMOS)
#
# Every pipeline script sources this file (same idea as CRASH/scripts/config.sh).
# Submit jobs FROM THE REPOSITORY FOLDER, either way works:
#       sbatch workflow/04_taxonomy/merge_metaphlan.sh          # plain sbatch
#       ./submit.sh --array-samples workflow/04_taxonomy/metaphlan.sh   # arrays
#
# This file only DEFINES variables and helper functions. Sourcing it never
# creates directories, submits jobs or downloads anything.
#
# Values set to __SET_ME__ are not yet known and MUST be filled in on COSMOS.
# Run `bash setup/check_config.sh` to list everything that is still missing.
# See README.md -> "Values you still need to obtain/check on COSMOS".
# =============================================================================

# Marker used for every unknown value. Helper functions refuse to run a step
# while a variable it needs still contains this string.
PLACEHOLDER="__SET_ME__"

# -----------------------------------------------------------------------------
# 0. Location of the pipeline code (this git repository)
# -----------------------------------------------------------------------------
# Resolved automatically from where config.sh lives. Job scripts find it via
# SLURM_SUBMIT_DIR, so always submit from the repository folder.
export PIPELINE_DIR="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# -----------------------------------------------------------------------------
# 1. Storage on COSMOS
# -----------------------------------------------------------------------------
# Your LUNARC project storage area (given by you).
export LUNARC_PROJECT_ROOT="/lunarc/nobackup/projects/lu2025-12-46/CHETANA"

# Writable analysis area. ALL pipeline outputs are written under here.
export ANALYSIS_ROOT="${LUNARC_PROJECT_ROOT}/ANALYSIS"

# DATA area. Treated as READ-ONLY: the pipeline never writes here, and
# assert_writable_target() refuses any output path inside it. Kept only so
# you can point to pre-existing input files later if needed.
export DATA_ROOT="${LUNARC_PROJECT_ROOT}/DATA"

# Project output directory: the stage folders (01_metadata ... 09_integration)
# live directly inside ANALYSIS, matching the folders already created on COSMOS.
# The pipeline code itself lives in ANALYSIS/scripts/.
export PROJECT_DIR="${ANALYSIS_ROOT}"

# -----------------------------------------------------------------------------
# 2. Dataset selection
# -----------------------------------------------------------------------------
# Public BioProject: HMP WGS Production Phase (healthy human subjects).
export BIOPROJECT="PRJNA48479"

# Short label for the body site being analysed; used in output file names.
export TARGET_LABEL="tonsil"

# Case-insensitive regular expression matched against SRA/BioSample body-site
# metadata to pick the samples. Check the report from 01_metadata to confirm it
# only matches what you intend (e.g. "Palatine Tonsils").
# Later, for oral-niche comparisons, change TARGET_LABEL and this regex
# (e.g. "tongue", "buccal", "supragingival") and rerun 01_metadata.
export TARGET_BODY_SITE_REGEX="tonsil"

# Keep only shotgun metagenomes (excludes any amplicon runs).
export LIBRARY_STRATEGY_KEEP="WGS"

# Keep ONE sample per person (HMP re-sampled some people at a later visit).
# Rule, fixed before any results are seen: keep the earliest visit if a visit
# number is in the BioSample metadata, otherwise the earliest-registered sample
# (lowest SRS accession). Removed samples are listed with the reason in
# 01_metadata/selected/<label>_excluded_samples.tsv. Set "false" to keep all.
export ONE_SAMPLE_PER_SUBJECT="true"

# -----------------------------------------------------------------------------
# 3. Stage directories (all under PROJECT_DIR, i.e. inside ANALYSIS)
# -----------------------------------------------------------------------------
# 01 — metadata tables downloaded from NCBI/ENA and the selected sample lists
export METADATA_DIR="${PROJECT_DIR}/01_metadata"
export METADATA_RAW_DIR="${METADATA_DIR}/raw"            # unfiltered project tables
export METADATA_SELECTED_DIR="${METADATA_DIR}/selected"  # filtered tonsil tables

# Files produced by 01_metadata and consumed by later stages:
export SELECTED_RUNS_FILE="${METADATA_SELECTED_DIR}/${TARGET_LABEL}_run_accessions.txt"   # one SRR per line
export SELECTED_SAMPLES_FILE="${METADATA_SELECTED_DIR}/${TARGET_LABEL}_sample_ids.txt"    # one SRS per line
export SAMPLE_RUN_MAP="${METADATA_SELECTED_DIR}/${TARGET_LABEL}_sample_runs.tsv"          # sample -> run(s)
export SELECTED_RUNS_TABLE="${METADATA_SELECTED_DIR}/${TARGET_LABEL}_runs_full.tsv"       # full metadata

# 02 — downloads
export DOWNLOAD_DIR="${PROJECT_DIR}/02_download"
export SRA_DIR="${DOWNLOAD_DIR}/sra"             # SRA Toolkit prefetch cache (.sra files)
export RAW_FASTQ_DIR="${DOWNLOAD_DIR}/raw_fastq" # raw per-run FASTQ (RUN_1/_2.fastq.gz)

# 03 — preprocessing (KneadData: trimming + host-read removal)
export PREPROCESS_DIR="${PROJECT_DIR}/03_preprocessing"
export KNEADDATA_DIR="${PREPROCESS_DIR}/kneaddata"         # KneadData output root
export CLEAN_FASTQ_DIR="${KNEADDATA_DIR}/clean_fastq"      # host-depleted, trimmed reads
export KNEADDATA_LOG_DIR="${KNEADDATA_DIR}/logs"           # per-sample KneadData logs
export KNEADDATA_QC_DIR="${KNEADDATA_DIR}/fastqc"          # FastQC before/after

# 04 — taxonomic profiling
export TAXONOMY_DIR="${PROJECT_DIR}/04_taxonomy"
export METAPHLAN_DIR="${TAXONOMY_DIR}/metaphlan"

# ---- MAG branch (runs separately from the taxonomy branch) -------------------
# 05 — assembly + binning + bin refinement + bin QC with nf-core/mag
export ASSEMBLY_DIR="${PROJECT_DIR}/05_assembly"
export NFMAG_SAMPLESHEET="${ASSEMBLY_DIR}/samplesheet.csv"   # made by make_samplesheet.sh
export NFMAG_OUTDIR="${ASSEMBLY_DIR}/nfcore_mag"              # nf-core/mag results
export NFMAG_WORKDIR="${ASSEMBLY_DIR}/nf_work"                # Nextflow work dir (large, temporary)
export NXF_CACHE_DIR="${ASSEMBLY_DIR}/nf_cache"               # Nextflow home + container images

# 06 — genes (planned)
export GENES_DIR="${PROJECT_DIR}/06_genes"

# 07 — MAG quality filtering, GTDB-Tk taxonomy and Bakta annotation
export MAG_DIR="${PROJECT_DIR}/07_mags"
export MAG_FILTERED_DIR="${MAG_DIR}/quality_filtered"   # MAGs passing the thresholds below
export ENZYME_DIR="${PROJECT_DIR}/08_enzyme_discovery"
export INTEGRATION_DIR="${PROJECT_DIR}/09_integration"

# SLURM stdout/stderr logs go to logs/ inside the repository folder (set by the
# "#SBATCH -o logs/..." lines, as in CRASH); they are ignored by git.

# Fallback temporary directory, used only when node-local scratch is unavailable
# (e.g. a light step run on the login node).
export FALLBACK_TMP_DIR="${PROJECT_DIR}/tmp"

# -----------------------------------------------------------------------------
# 4. SLURM settings (applied by submit.sh, because #SBATCH lines cannot read
#    shell variables)
# -----------------------------------------------------------------------------
# Compute allocation (given by supervisor). This is also written in every
# script's "#SBATCH -A" line; change both if it ever changes (see README).
export SLURM_ACCOUNT="lu2026-2-117"

# E-mail for job notifications (also in each script's "#SBATCH --mail-user").
export SLURM_EMAIL="ch2068va-s@student.lu.se"
export SLURM_MAIL_TYPE="END,FAIL"

# Partition. Empty = cluster default. Only set this if LUNARC tells you to.
export SLURM_PARTITION=""

# Maximum number of array tasks running at the same time.
export SLURM_MAX_PARALLEL="10"

# -----------------------------------------------------------------------------
# 5. Software
# -----------------------------------------------------------------------------
# Environment module that provides conda/mamba on COSMOS (checked with
# `module spider Miniforge3`; same version CRASH used).
export CONDA_MODULE="Miniforge3/25.3.0-3"

# Optional: directory holding your conda environments, if you create them with
# `conda env create -p <dir>/<name>` instead of the default location (useful
# because environments in $HOME can exceed the home quota). Empty = activate
# environments by name.
export CONDA_ENV_DIR=""

# Conda environment names (created from envs/*.yml; see README).
export ENV_SRA="binp52_sra"              # sra-tools (prefetch, fasterq-dump)
export ENV_KNEADDATA="binp52_kneaddata"  # kneaddata, trimmomatic, bowtie2, fastqc, multiqc
export ENV_METAPHLAN="binp52_metaphlan"  # MetaPhlAn 4
export ENV_NEXTFLOW="binp52_nextflow"    # Nextflow (runs nf-core/mag; tools come in containers)

# GTDB-Tk + Bakta: the supervisor's scripts activate a shared group environment
# by path. Set the path once confirmed (candidate seen in his scripts:
# /lunarc/nobackup/projects/lu2025-12-46/util/gtdbtk2). If it is not usable,
# create your own from envs/gtdbtk.yml / envs/bakta.yml and put the env NAME here.
export ENV_GTDBTK="/lunarc/nobackup/projects/lu2025-12-46/util/gtdbtk2"
export ENV_BAKTA="/lunarc/nobackup/projects/lu2025-12-46/util/gtdbtk2"

# nf-core/mag version (as advised) and the Nextflow version used to run it.
# 2.5.1 is an older release (2023); a Nextflow version from that period is
# pinned so the pipeline code is interpreted as its authors intended.
export NFMAG_VERSION="2.5.1"
export NXF_VER="23.10.1"

# -----------------------------------------------------------------------------
# 6. Reference databases (paths unknown until you check/install them)
# -----------------------------------------------------------------------------
# First check whether LUNARC already provides shared copies before downloading
# your own (they are large). If you install them yourself, put them under
# ANALYSIS (not DATA, not $HOME).

# KneadData host reference: directory containing the human Bowtie2 index
# (e.g. from `kneaddata_database --download human_genome bowtie2 <dir>`).
# The supervisor's script uses $HOME/util/KD_HOMO in his home folder; ask
# whether a copy readable by you exists (e.g. under lu2025-12-46/util).
export KNEADDATA_DB="/lunarc/nobackup/projects/lu2025-12-46/Vik/Util/kneaddata_db"

# MetaPhlAn 4 database directory and index name.
export METAPHLAN_DB_DIR="/lunarc/nobackup/projects/lu2025-12-46/util/util_12-47/metaphlan_db_202403"
export METAPHLAN_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"  # change if you install another version

# GTDB-Tk reference data directory. The supervisor uses GTDB release 226
# (needs GTDB-Tk >= 2.4.1). Candidate seen in his scripts (check access):
# /lunarc/nobackup/projects/lu2024-12-47/util/GTDB/release226
export GTDBTK_DATA_PATH="/lunarc/nobackup/projects/lu2025-12-46/util/util_12-47/GTDB/release226"

# Bakta database directory. Candidate seen in his scripts (check access):
# /lunarc/nobackup/projects/lu2025-12-46/util/util_12-47/BAKTA/db
export BAKTA_DB="/lunarc/nobackup/projects/lu2025-12-46/util/util_12-47/BAKTA/db"

# -----------------------------------------------------------------------------
# 7. Tool parameters (scientific choices; adjust and record in your thesis)
# -----------------------------------------------------------------------------
# Download method: "sra" (prefetch + fasterq-dump) or "ena" (direct FASTQ from
# ENA with MD5 check). "ena" only needs curl and is light enough to run on a
# login node if compute nodes have no internet access.
export DOWNLOAD_METHOD="sra"
export KEEP_SRA_FILES="false"       # delete .sra files after successful conversion

# Trimmomatic adapter set used by KneadData: NexteraPE | TruSeq2 | TruSeq3 | none.
# HMP Phase 1 used Illumina TruSeq-type libraries; confirm with FastQC adapter
# content on a few raw files before the full run.
export KNEADDATA_SEQUENCER_SOURCE="TruSeq3"
export KNEADDATA_MAX_MEMORY="16g"   # memory given to Trimmomatic inside KneadData
export KNEADDATA_EXTRA_OPTS=""      # any extra kneaddata options, space separated

# BBTools repair.sh re-synchronises read pairs after KneadData (as in the
# supervisor's assembly script) so assemblers and mappers get matched pairs.
export REPAIR_MEMORY="16g"

# ---- nf-core/mag settings (written into the params file by run_nfcore_mag.sh)
# Assembler: "megahit", "spades" or "both". One assembler avoids duplicate MAGs
# of the same genome; MEGAHIT is faster and lighter. Confirm with supervisor.
export NFMAG_ASSEMBLER="megahit"
# Binning coverage: "group" = map every sample to every assembly in the same
# group (all samples are put in group 0) — usually more and better MAGs;
# "own" = each sample only to its own assembly (like the supervisor's script).
export NFMAG_BINNING_MAP_MODE="group"
# Largest single task nf-core may request. Check node sizes with:
#   sinfo -o "%P %c %m"     (cores, memory in MB per node)
export NFMAG_MAX_CPUS="48"
export NFMAG_MAX_MEMORY="240.GB"
export NFMAG_MAX_TIME="160.h"               # same walltime the supervisor uses
export NFMAG_QUEUE_SIZE="20"                # max nf-core jobs in the queue at once

# MAG quality thresholds (same as the supervisor's metaWRAP refinement -c 50 -x 10)
export MAG_MIN_COMPLETENESS="50"
export MAG_MAX_CONTAMINATION="10"

# =============================================================================
# 8. Helper functions (used by every script)
# =============================================================================

# Timestamped message to stderr.
log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }

# Print an error and stop.
die() { log "ERROR: $*"; exit 1; }

# True if a value is empty or still contains the placeholder.
is_placeholder() { [[ -z "${1:-}" || "$1" == *"${PLACEHOLDER}"* ]]; }

# Stop unless every named variable has a real value.
#   require_vars SLURM_ACCOUNT KNEADDATA_DB
require_vars() {
    local v missing=()
    for v in "$@"; do
        is_placeholder "${!v:-}" && missing+=("$v")
    done
    (( ${#missing[@]} == 0 )) || die "Please set these in config.sh first: ${missing[*]}"
}

# Stop unless every named variable points to an existing file/directory.
#   require_paths KNEADDATA_DB METAPHLAN_DB_DIR
require_paths() {
    local v
    for v in "$@"; do
        [[ -e "${!v}" ]] || die "${v} does not exist: ${!v}"
    done
}

# Refuse to write into DATA (treated as read-only) or outside ANALYSIS.
assert_writable_target() {
    local p="$1"
    case "$p" in
        "${DATA_ROOT}"|"${DATA_ROOT}"/*) die "Refusing to write inside DATA (read-only by design): $p" ;;
        "${ANALYSIS_ROOT}"|"${ANALYSIS_ROOT}"/*) ;;
        *) die "Output path is outside ANALYSIS_ROOT (${ANALYSIS_ROOT}): $p" ;;
    esac
}

# Create output directories after checking they are allowed.
make_out_dirs() {
    local d
    for d in "$@"; do
        assert_writable_target "$d"
        mkdir -p "$d"
    done
}

# Load the conda module and activate an environment.
# Lmod and conda scripts reference unset variables, so nounset is paused.
load_conda_env() {
    local env="$1"
    require_vars CONDA_MODULE
    is_placeholder "${env}" && die "Conda environment not set in config.sh"
    set +u
    module purge
    module load "${CONDA_MODULE}"
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    if [[ "${env}" == /* ]]; then
        conda activate "${env}"                      # environment given as a path
    elif [[ -n "${CONDA_ENV_DIR}" ]]; then
        conda activate "${CONDA_ENV_DIR}/${env}"
    else
        conda activate "${env}"
    fi
    set -u
    log "Activated conda environment: ${env} (${CONDA_PREFIX})"
}

# Create a per-task working directory and set $WORK and $TMPDIR.
# Uses node-local scratch ($SNIC_TMP inside LUNARC jobs) when available.
make_scratch() {
    local tool="$1" id="$2" base
    if [[ "${SNIC_TMP:-}" == /* && -d "${SNIC_TMP}" ]]; then
        base="${SNIC_TMP}"
    else
        log "WARNING: SNIC_TMP not a usable directory; using ${FALLBACK_TMP_DIR} (slower shared storage)"
        make_out_dirs "${FALLBACK_TMP_DIR}"
        base="${FALLBACK_TMP_DIR}"
    fi
    WORK="${base}/${tool}_${USER:-user}_${SLURM_JOB_ID:-local$$}_${id}"
    mkdir -p "${WORK}"/{in,tmp,out}
    TMPDIR="${WORK}/tmp"
    export WORK TMPDIR
}

# Number of threads for this task. Works with --cpus-per-task (CRASH style)
# and with --ntasks-per-node (supervisor style).
threads() { echo "${SLURM_CPUS_PER_TASK:-${SLURM_NTASKS_PER_NODE:-${SLURM_CPUS_ON_NODE:-4}}}"; }

# Count non-empty, non-comment lines in a list file.
count_items() { awk 'NF && $1 !~ /^#/ {n++} END {print n+0}' "$1"; }

# Return item N (1-based) of a list file; N defaults to SLURM_ARRAY_TASK_ID.
array_item() {
    local file="$1" idx="${2:-${SLURM_ARRAY_TASK_ID:-}}" item
    [[ -n "$idx" ]] || die "No index given and SLURM_ARRAY_TASK_ID is not set"
    [[ -s "$file" ]] || die "List file missing or empty: $file (run 01_metadata first)"
    item="$(awk 'NF && $1 !~ /^#/ {print $1}' "$file" | sed -n "${idx}p")"
    [[ -n "$item" ]] || die "No item number ${idx} in ${file}"
    printf '%s\n' "$item"
}

# Resolve the ID a script should process:
#   - an explicit ID argument (e.g. SRR123456 / SRS012345) is used as-is,
#   - a number is treated as a line number in the list file,
#   - no argument uses SLURM_ARRAY_TASK_ID.
resolve_item() {
    local file="$1" arg="${2:-}"
    if [[ -z "$arg" ]]; then
        array_item "$file"
    elif [[ "$arg" =~ ^[0-9]+$ ]]; then
        array_item "$file" "$arg"
    else
        printf '%s\n' "$arg"
    fi
}

# List the run accessions belonging to one sample (from SAMPLE_RUN_MAP).
runs_for_sample() {
    [[ -s "${SAMPLE_RUN_MAP}" ]] || die "Missing ${SAMPLE_RUN_MAP} (run 01_metadata first)"
    awk -F'\t' -v s="$1" 'NR > 1 && $1 == s {print $2}' "${SAMPLE_RUN_MAP}"
}

# True if a gzip file exists and contains at least one byte of data.
nonempty_gz() { [[ -s "$1" ]] && [[ -n "$(gzip -cd "$1" 2>/dev/null | head -c 1)" ]]; }

# Marker files record finished work so resubmitted jobs skip completed items.
done_marker() { echo "$1/.done/$2.done"; }
mark_done() { mkdir -p "$1/.done"; date '+%F %T' > "$(done_marker "$1" "$2")"; }
is_done() { [[ -f "$(done_marker "$1" "$2")" ]]; }
