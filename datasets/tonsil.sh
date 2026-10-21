# shellcheck shell=bash
# datasets/tonsil.sh — HMP tonsil samples (sourced by config.sh, DATASET=tonsil)
export BIOPROJECT="PRJNA48479"
export TARGET_LABEL="tonsil"
export TARGET_BODY_SITE_REGEX="tonsil"      # all matches were "G_DNA_Palatine Tonsils"
export LIBRARY_STRATEGY_KEEP="WGS"
export ONE_SAMPLE_PER_SUBJECT="true"
export PROJECT_DIR="${ANALYSIS_ROOT}/hmp_tonsil"   # tonsil results (moved Oct 2026)
