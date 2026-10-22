# shellcheck shell=bash
# datasets/hmp_oral.sh — HMP oral sites of the SAME 23 people as the tonsil set
# (DATASET=hmp_oral). Supervisor's decision: same participants, all their oral
# sites, same selection rule and pipeline settings as the tonsils.
export BIOPROJECT="PRJNA48479"
export TARGET_LABEL="hmp_oral"
# All HMP oral sites except tonsils (hard palate has none of these people).
export TARGET_BODY_SITE_REGEX="Tongue dorsum|Supragingival plaque|Buccal mucosa|Subgingival plaque|Throat|Keratinized gingiva|Saliva|Hard palate"
export LIBRARY_STRATEGY_KEEP="WGS"
# One sample per person AND site; earliest-registered sample if repeated.
export ONE_SAMPLE_PER_SUBJECT="site"
# Only the 23 tonsil participants (made from tonsil_runs_full.tsv).
export SUBJECTS_FILE="${PIPELINE_DIR}/metadata/hmp_oral/tonsil_subjects.txt"
# Illumina only (excludes the single 454 run).
export PLATFORM_REGEX="Illumina"
export PROJECT_DIR="${ANALYSIS_ROOT}/hmp_oral"
# Per person and site, keep the sample closest to the person's tonsil sample
# (same visit whenever one exists; supervisor's decision, Oct 2026).
export MATCH_TO_FILE="${PIPELINE_DIR}/metadata/tonsil_runs_full.tsv"
# nf-core/mag: one group per body site (binning co-abundance within each site)
export SAMPLESHEET_GROUP_BY="site"
