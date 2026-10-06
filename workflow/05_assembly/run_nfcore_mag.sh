#!/bin/bash
#SBATCH -A lu2026-2-117
#SBATCH -J nfcore_mag
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH -t 160:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=ch2068va-s@student.lu.se
#SBATCH --mail-type=END,FAIL
# =============================================================================
# 05_assembly/run_nfcore_mag.sh — MAG branch: nf-core/mag (assembly -> binning
#                                 -> DAS Tool refinement -> CheckM bin QC)
#
# This job is only the Nextflow "head": it uses 2 cores and submits every
# pipeline step as its own SLURM job (account from config.sh).
#
#   sbatch workflow/05_assembly/run_nfcore_mag.sh test   # 1) small nf-core test dataset
#   sbatch workflow/05_assembly/run_nfcore_mag.sh        # 2) your samples
#
# Resubmitting continues where it stopped (-resume).
# Settings that differ from nf-core defaults (all in the params file):
#   skip_clipping       reads were already trimmed by KneadData (Trimmomatic)
#   (no host removal)   KneadData already removed human reads
#   assembler           NFMAG_ASSEMBLER (default megahit only)
#   skip_concoct        MetaBAT2 + MaxBin2, as in the supervisor's metaWRAP script
#   refine_bins_dastool combine MetaBAT2 + MaxBin2 bins (like metaWRAP bin_refinement)
#   binqc_tool checkm   CheckM, as in the supervisor's script
#   skip_gtdbtk         run separately with GTDB release 226 (07_mags/gtdbtk.sh);
#                       nf-core/mag 2.5.1 ships GTDB-Tk 2.3.2 (release 214 only)
#   skip_prokka         annotation is done with Bakta instead (07_mags/bakta.sh)
# =============================================================================
set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${PIPELINE_DIR}/config.sh" ]] || { echo "config.sh not found in ${PIPELINE_DIR}: submit from the repository folder" >&2; exit 1; }
source "${PIPELINE_DIR}/config.sh"

MODE="${1:-run}"
require_vars SLURM_ACCOUNT CONDA_MODULE

# ---- Nextflow environment -------------------------------------------------------
make_out_dirs "${ASSEMBLY_DIR}" "${NXF_CACHE_DIR}/home" "${NXF_CACHE_DIR}/singularity"
export NXF_VER                                   # pinned Nextflow version (config.sh)
export NXF_HOME="${NXF_CACHE_DIR}/home"          # pipeline code + Nextflow jars
export NXF_SINGULARITY_CACHEDIR="${NXF_CACHE_DIR}/singularity"
export NXF_OPTS="-Xms1g -Xmx4g"
export NF_SLURM_ACCOUNT="${SLURM_ACCOUNT}"
export NFMAG_QUEUE_SIZE
load_conda_env "${ENV_NEXTFLOW}"
nextflow -version >&2

CONF="${PIPELINE_DIR}/workflow/05_assembly/cosmos.config"

# ---- Test mode: tiny public dataset, checks containers + SLURM setup -------------
if [[ "${MODE}" == "test" ]]; then
    OUT="${ASSEMBLY_DIR}/nfcore_mag_test"
    make_out_dirs "${OUT}"
    cd "${ASSEMBLY_DIR}"
    nextflow run nf-core/mag -r "${NFMAG_VERSION}" \
        -profile test,singularity -c "${CONF}" \
        -work-dir "${ASSEMBLY_DIR}/nf_work_test" \
        --outdir "${OUT}" -resume
    log "Test finished: ${OUT}"
    exit 0
fi

# ---- Real run --------------------------------------------------------------------
require_vars NFMAG_MAX_CPUS NFMAG_MAX_MEMORY
[[ -s "${NFMAG_SAMPLESHEET}" ]] || die "Missing ${NFMAG_SAMPLESHEET} (run make_samplesheet.sh)"
make_out_dirs "${NFMAG_OUTDIR}" "${NFMAG_WORKDIR}"

case "${NFMAG_ASSEMBLER}" in
    megahit) SKIP_MEGAHIT=false; SKIP_SPADES=true  ;;
    spades)  SKIP_MEGAHIT=true;  SKIP_SPADES=false ;;
    both)    SKIP_MEGAHIT=false; SKIP_SPADES=false ;;
    *) die "NFMAG_ASSEMBLER must be megahit, spades or both" ;;
esac

PARAMS="${ASSEMBLY_DIR}/params_nfcore_mag.yaml"
cat > "${PARAMS}" <<YAML
# Written by run_nfcore_mag.sh on $(date '+%F %T') — keep with the results.
input: "${NFMAG_SAMPLESHEET}"
outdir: "${NFMAG_OUTDIR}"
skip_clipping: true
skip_megahit: ${SKIP_MEGAHIT}
skip_spades: ${SKIP_SPADES}
skip_spadeshybrid: true
binning_map_mode: "${NFMAG_BINNING_MAP_MODE}"
skip_concoct: true
refine_bins_dastool: true
postbinning_input: "refined_bins_only"
binqc_tool: "checkm"
save_checkm_data: true
skip_gtdbtk: true
skip_prokka: true
skip_krona: true
max_cpus: ${NFMAG_MAX_CPUS}
max_memory: "${NFMAG_MAX_MEMORY}"
max_time: "${NFMAG_MAX_TIME}"
YAML
log "Parameters -> ${PARAMS}"; cat "${PARAMS}" >&2

cd "${ASSEMBLY_DIR}"      # .nextflow.log and .nextflow/ (needed for -resume) live here
nextflow run nf-core/mag -r "${NFMAG_VERSION}" \
    -profile singularity -c "${CONF}" \
    -params-file "${PARAMS}" \
    -work-dir "${NFMAG_WORKDIR}" \
    -resume

log "nf-core/mag finished: ${NFMAG_OUTDIR}"
log "Next: bash workflow/07_mags/filter_mags.sh"
log "When all later steps are done, the work dir can be deleted: ${NFMAG_WORKDIR}"
