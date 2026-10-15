# shellcheck shell=bash
# datasets/oral.sh — oral microbiome dataset ENA PRJEB36291 (DATASET=oral)
export BIOPROJECT="PRJEB36291"
export TARGET_LABEL="oral"
# EXPLORATION: "." matches any body-site value, so the first metadata run lists
# ALL shotgun runs. Replace with the exact site wording after reviewing the report.
export TARGET_BODY_SITE_REGEX="."
export LIBRARY_STRATEGY_KEEP="WGS"
export ONE_SAMPLE_PER_SUBJECT="true"
export PROJECT_DIR="${ANALYSIS_ROOT}/oral_PRJEB36291"   # own output tree
