#!/bin/bash
#SBATCH -A lu2026-2-117
#SBATCH -J humann_post
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=192G
#SBATCH -t 24:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=ch2068va-s@student.lu.se
#SBATCH --mail-type=END,FAIL
# =============================================================================
# 06_genes/humann_postprocess.sh — combine and normalise HUMAnN 4 results
#
# As CRASH's humann_utility.sh:
#   1. regroup each sample's gene families (UniRef90) to KEGG Orthologs (KO),
#      EC numbers (level-4 EC) and EggNOG groups (mapping files in
#      HUMANN_UTILITY_DB);
#   2. join the samples into one table per output type;
#   3. normalise: gene families, KO, EC, EggNOG, reactions -> CPM;
#      pathway abundance -> relative abundance;
#   4. split each normalised table into stratified (per species) and
#      unstratified (community total) versions.
# Stops if any selected sample is unfinished (add "--partial" to override).
#
#   export DATASET=<name>
#   sbatch workflow/06_genes/humann_postprocess.sh
# Output: <dataset>/06_genes/humann/tables/{joined,normalized,split}/
# =============================================================================
set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${PIPELINE_DIR}/config.sh" ]] || { echo "config.sh not found in ${PIPELINE_DIR}: submit from the repository folder" >&2; exit 1; }
source "${PIPELINE_DIR}/config.sh"

require_vars ENV_HUMANN HUMANN_UTILITY_DB
require_paths HUMANN_UTILITY_DB
P="${HUMANN_DIR}/per_sample"
RG="${HUMANN_DIR}/regrouped"
TB="${HUMANN_DIR}/tables"
FLAT="${TB}/flat"
make_out_dirs "${RG}" "${TB}/joined" "${TB}/normalized" "${TB}/split"

n_exp="$(count_items "${SELECTED_SAMPLES_FILE}")"
n_done="$(find "${HUMANN_DIR}/.done" -type f -name '*.done' 2>/dev/null | wc -l)"
log "HUMAnN finished for ${n_done} of ${n_exp} selected samples"
if (( n_done < n_exp )) && [[ "${1:-}" != "--partial" ]]; then
    die "Not all samples are finished; rerun humann.sh first (or pass --partial)"
fi

load_conda_env "${ENV_HUMANN}"

# ---- 1. flat folders per output type (links to the per-sample files) --------
rm -rf "${FLAT}"
for type in genefamilies reactions pathabundance pathcoverage; do
    mkdir -p "${FLAT}/${type}"
    find "${P}" -mindepth 2 -maxdepth 2 -name "*_${type}.tsv" -exec ln -s {} "${FLAT}/${type}/" \;
    [[ -n "$(ls -A "${FLAT}/${type}")" ]] || rmdir "${FLAT}/${type}"
done
[[ -d "${FLAT}/genefamilies" ]] || die "No gene-family tables found under ${P}"

# ---- 2. regroup gene families (kept between runs; finished files skipped) ---
for g in ko level4ec eggnog; do
    map="${HUMANN_UTILITY_DB}/map_${g}_uniref90.txt.gz"
    if [[ ! -s "${map}" ]]; then log "No mapping file for ${g} (${map}); skipping"; continue; fi
    mkdir -p "${RG}/${g}" "${FLAT}/${g}"
    for gf in "${FLAT}/genefamilies"/*.tsv; do
        out="${RG}/${g}/$(basename "${gf}" .tsv)_${g}.tsv"
        [[ -s "${out}" ]] || humann_regroup_table -i "${gf}" -c "${map}" -o "${out}"
        ln -s "${out}" "${FLAT}/${g}/"
    done
    log "Regrouped to ${g}: $(find "${FLAT}/${g}" -maxdepth 1 -name '*.tsv' | wc -l) samples"
done

# ---- 3. join, normalise, split ------------------------------------------------
for dir in "${FLAT}"/*/; do
    name="$(basename "${dir}")"
    joined="${TB}/joined/${name}.tsv"
    humann_join_tables -i "${dir}" -o "${joined}"
    units="cpm"; [[ "${name}" == pathabundance ]] && units="relab"
    if [[ "${name}" == pathcoverage ]]; then
        cp "${joined}" "${TB}/normalized/${name}.tsv"; norm="${TB}/normalized/${name}.tsv"
    else
        norm="${TB}/normalized/${name}_${units}.tsv"
        humann_renorm_table -i "${joined}" -o "${norm}" --units "${units}"
    fi
    humann_split_stratified_table -i "${norm}" -o "${TB}/split"
    log "${name}: joined, normalised (${units}) and split"
done

log "Tables in ${TB}:"
ls -lh "${TB}/normalized" "${TB}/split" >&2
