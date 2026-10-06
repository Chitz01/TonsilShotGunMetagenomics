#!/usr/bin/env bash
# =============================================================================
# setup/check_config.sh — report which config.sh values are still missing
#
# Read-only: prints every important variable with a status and checks that
# paths exist. Run it on the COSMOS login node:  bash setup/check_config.sh
# Exit code is 1 while any required value is still a placeholder.
# =============================================================================
set -uo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PIPELINE_DIR
source "${PIPELINE_DIR}/config.sh"

problems=0
row() { printf '%-26s %-14s %s\n' "$1" "$2" "$3"; }

check_value() {   # variable must not be a placeholder
    local v="$1"
    if is_placeholder "${!v:-}"; then row "$v" "PLACEHOLDER" "${!v:-<empty>}"; problems=$((problems+1))
    else row "$v" "set" "${!v}"; fi
}
check_path() {    # variable must be set and exist
    local v="$1"
    if is_placeholder "${!v:-}"; then row "$v" "PLACEHOLDER" "${!v:-<empty>}"; problems=$((problems+1))
    elif [[ -e "${!v}" ]]; then row "$v" "exists" "${!v}"
    else row "$v" "NOT FOUND" "${!v}"; problems=$((problems+1)); fi
}

echo "== Storage =================================================================="
check_path  LUNARC_PROJECT_ROOT
check_path  ANALYSIS_ROOT
if [[ -d "${ANALYSIS_ROOT}" && -w "${ANALYSIS_ROOT}" ]]; then row "ANALYSIS writable" "yes" ""
else row "ANALYSIS writable" "NO" "pipeline outputs need this"; problems=$((problems+1)); fi
if [[ -d "${DATA_ROOT}" ]]; then
    row "DATA_ROOT" "read-only use" "${DATA_ROOT} (pipeline never writes here)"
fi
row "PROJECT_DIR" "$( [[ -d ${PROJECT_DIR} ]] && echo exists || echo 'not created')" "${PROJECT_DIR}"

echo "== Dataset =================================================================="
row BIOPROJECT "set" "${BIOPROJECT}"
row TARGET_BODY_SITE_REGEX "set" "${TARGET_BODY_SITE_REGEX}"

echo "== SLURM ===================================================================="
check_value SLURM_ACCOUNT
if is_placeholder "${SLURM_EMAIL}"; then row SLURM_EMAIL "optional" "not set -> no e-mails"
else row SLURM_EMAIL "set" "${SLURM_EMAIL}"; fi
row SLURM_PARTITION "$( [[ -n ${SLURM_PARTITION} ]] && echo set || echo default)" "${SLURM_PARTITION:-<cluster default>}"

echo "== Software ================================================================="
check_value CONDA_MODULE
[[ -n "${CONDA_ENV_DIR}" ]] && check_path CONDA_ENV_DIR
for c in sbatch module; do
    if command -v "$c" >/dev/null 2>&1 || type "$c" >/dev/null 2>&1; then row "command: $c" "found" ""
    else row "command: $c" "missing" "expected on COSMOS"; fi
done

echo "== Software environments =============================================="
row ENV_NEXTFLOW "name" "${ENV_NEXTFLOW} (create from envs/nextflow.yml)"
check_value ENV_GTDBTK
check_value ENV_BAKTA

echo "== Databases (needed from the stage noted) =================================="
check_path KNEADDATA_DB        # 03_preprocessing (shared)
check_path METAPHLAN_DB_DIR    # 04_taxonomy (taxonomy branch)
row METAPHLAN_INDEX "set" "${METAPHLAN_INDEX}"
check_path GTDBTK_DATA_PATH    # 07_mags (MAG branch)
check_path BAKTA_DB            # 07_mags (MAG branch)

echo "== nf-core/mag (MAG branch) ================================================="
row NFMAG_VERSION "set" "${NFMAG_VERSION} (Nextflow ${NXF_VER})"
row NFMAG_ASSEMBLER "set" "${NFMAG_ASSEMBLER}"
check_value NFMAG_MAX_CPUS
check_value NFMAG_MAX_MEMORY
for c in apptainer singularity; do
    if command -v "$c" >/dev/null 2>&1; then row "command: $c" "found" ""; else row "command: $c" "missing" ""; fi
done

echo "== Scratch =================================================================="
if [[ -n "${SNIC_TMP:-}" ]]; then row SNIC_TMP "set" "${SNIC_TMP}"
else row SNIC_TMP "not set" "normal on login node; check inside a job"; fi

echo
if (( problems > 0 )); then
    echo "${problems} value(s) still need attention. See README: 'Values you still need to obtain/check on COSMOS'."
    exit 1
fi
echo "All checked values are set."
