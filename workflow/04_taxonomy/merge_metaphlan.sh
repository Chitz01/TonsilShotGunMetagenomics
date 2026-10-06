#!/bin/bash
#SBATCH -A lu2026-2-117
#SBATCH -J metaphlan_merge
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -t 00:30:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=ch2068va-s@student.lu.se
#SBATCH --mail-type=END,FAIL
# =============================================================================
# 04_taxonomy/merge_metaphlan.sh — merge per-sample profiles into tables
#
#   ./submit.sh workflow/04_taxonomy/merge_metaphlan.sh
# Outputs in METAPHLAN_DIR:
#   merged_abundance_table.tsv   all ranks (MetaPhlAn format, relative abundance)
#   species_abundance.tsv        species rows only (your taxonomic table)
#   genus_abundance.tsv          genus rows only
# =============================================================================
set -euo pipefail
# Find config.sh: submit from the repository folder (SLURM_SUBMIT_DIR).
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${PIPELINE_DIR}/config.sh" ]] || { echo "config.sh not found in ${PIPELINE_DIR}: submit from the repository folder" >&2; exit 1; }
source "${PIPELINE_DIR}/config.sh"

shopt -s nullglob
PROFILES=( "${METAPHLAN_DIR}/profiles/"*_profile.txt )
(( ${#PROFILES[@]} > 0 )) || die "No profiles in ${METAPHLAN_DIR}/profiles"
log "Merging ${#PROFILES[@]} profiles (expected $(count_items "${SELECTED_SAMPLES_FILE}"))"

load_conda_env "${ENV_METAPHLAN}"
MERGED="${METAPHLAN_DIR}/merged_abundance_table.tsv"
merge_metaphlan_tables.py "${PROFILES[@]}" > "${MERGED}"

# Keep header lines ('#...' and 'clade_name ...') plus rows of one rank.
awk '/^#|^clade_name/ {print; next} /\|s__/ && !/\|t__/ {print}' \
    "${MERGED}" > "${METAPHLAN_DIR}/species_abundance.tsv"
awk '/^#|^clade_name/ {print; next} /\|g__/ && !/\|s__/ {print}' \
    "${MERGED}" > "${METAPHLAN_DIR}/genus_abundance.tsv"

log "Tables written to ${METAPHLAN_DIR}"
