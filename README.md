# TonsilShotGunMetagenomics

Mining the human tonsil and oral microbiomes for microbial enzymes with
potential to degrade tire-derived micro- and nanoplastics.
BINP52 MSc project, Lund University — Chetana Varada (supervisor: Ghjuvan Grimaud).
Runs on LUNARC COSMOS with SLURM.

**Data:** HMP shotgun metagenomes, BioProject PRJNA48479. Only tonsil samples
are used: 25 samples from 23 people matched; one sample per person is kept
(23 samples), see [Sample selection](#sample-selection).

## Workflow

```text
            01_metadata ──> 02_download ──> 03_preprocessing (KneadData + pair repair)
                                                     │
                     ┌───────────────────────────────┴──────────────────────────────┐
     TAXONOMY BRANCH │                                                   MAG BRANCH │
                     ▼                                                              ▼
            04_taxonomy (MetaPhlAn 4)                      05_assembly (nf-core/mag 2.5.1:
            -> species / genus tables                       MEGAHIT, MetaBAT2 + MaxBin2,
                                                            DAS Tool, CheckM)
                                                                    │
                                                                    ▼
                                                           07_mags: filter 50/10
                                                           -> GTDB-Tk (R226) -> Bakta
                     └───────────────────────┬──────────────────────────────────────┘
                                             ▼
                      06_genes, 08_enzyme_discovery, 09_integration (planned)
```

The two branches are independent after 03: either can be run, rerun or
changed without touching the other.

## Repository layout

```text
TonsilShotGunMetagenomics/
├── config.sh                  all paths, account, databases, parameters + helpers
├── submit.sh                  sbatch wrapper for ARRAY jobs (counts the samples for you)
├── setup/                     check_config.sh, create_directories.sh
├── envs/                      conda environments (sra, kneaddata, metaphlan, nextflow, gtdbtk, bakta)
├── workflow/
│   ├── 01_metadata/           find tonsil samples (+ one per person)
│   ├── 02_download/           download selected runs
│   ├── 03_preprocessing/      KneadData + repair.sh, summary
│   ├── 04_taxonomy/           MetaPhlAn 4 + merged tables          (taxonomy branch)
│   ├── 05_assembly/           nf-core/mag launcher + COSMOS config (MAG branch)
│   ├── 07_mags/               filter_mags, gtdbtk, bakta           (MAG branch)
│   └── 06_genes/ 08_enzyme_discovery/ 09_integration/   (planned)
├── metadata/                  reviewed sample tables (committed)
├── docs/                      notes, exported environment versions
└── logs/                      SLURM .out/.err (not committed)
```

On COSMOS the code lives in `ANALYSIS/scripts/TonsilShotGunMetagenomics/`, and
all results go to the numbered folders directly under
`/lunarc/nobackup/projects/lu2025-12-46/CHETANA/ANALYSIS/`. `DATA/` is never
written to.

## How to submit jobs

Always work from the repository folder. Account (`lu2026-2-117`), e-mail,
walltime and log paths are in each script's `#SBATCH` header, as in CRASH and
the supervisor's scripts.

```bash
cd /lunarc/nobackup/projects/lu2025-12-46/CHETANA/ANALYSIS/scripts/TonsilShotGunMetagenomics
sbatch workflow/04_taxonomy/merge_metaphlan.sh                           # single job
./submit.sh --array-samples workflow/03_preprocessing/kneaddata.sh       # one task per sample
./submit.sh --dry-run --array-runs workflow/02_download/download_runs.sh # show, don't submit
```

If the account ever changes, update all headers and `config.sh` at once:

```bash
grep -rl "lu2026-2-117" config.sh workflow | xargs sed -i 's/lu2026-2-117/NEW-ACCOUNT/g'
```

## Running the pipeline

```bash
# 01  metadata + selection (login node)
bash workflow/01_metadata/run_metadata.sh            # or --select to redo only the selection
# 02  download
./submit.sh --array-runs workflow/02_download/download_runs.sh
# 03  KneadData (trim, remove human reads, re-pair reads)
./submit.sh --array-samples workflow/03_preprocessing/kneaddata.sh
sbatch workflow/03_preprocessing/kneaddata_summary.sh

# TAXONOMY BRANCH
./submit.sh --array-samples workflow/04_taxonomy/metaphlan.sh
sbatch workflow/04_taxonomy/merge_metaphlan.sh

# MAG BRANCH
bash   workflow/05_assembly/make_samplesheet.sh
sbatch workflow/05_assembly/run_nfcore_mag.sh test    # first: nf-core test dataset
sbatch workflow/05_assembly/run_nfcore_mag.sh         # then: your samples
bash   workflow/07_mags/filter_mags.sh
sbatch workflow/07_mags/gtdbtk.sh
sbatch workflow/07_mags/bakta.sh
```

Every step skips finished samples, so a failed job can simply be resubmitted;
nf-core/mag continues with `-resume`.

## Sample selection

`01_metadata` searches SRA RunInfo and BioSample attributes for "tonsil" and
keeps WGS runs. All matches had `isolation_source = G_DNA_Palatine Tonsils`.

**One sample per person** (`ONE_SAMPLE_PER_SUBJECT=true`). HMP re-sampled some
people at a later visit; two subjects had two tonsil samples. Keeping both
would count those people twice, so one is removed. The rule was fixed before
looking at any results, so it cannot bias them:

1. keep the earliest visit, if a visit number is in the BioSample metadata;
2. otherwise keep the earliest-registered sample (lowest SRS accession).

The two pairs have almost identical depth (about 12 Gbases each), so the
choice does not trade away data. Removed samples and the reason are written to
`01_metadata/selected/tonsil_excluded_samples.tsv` and the selection report.

## MAG branch: nf-core/mag settings

nf-core/mag 2.5.1, run with Nextflow 23.10.1 (pinned in `config.sh`) and
Singularity/Apptainer containers. Non-default choices, written to
`05_assembly/params_nfcore_mag.yaml` at each run:

| Setting | Why |
| --- | --- |
| `skip_clipping`, no host removal | already done by KneadData |
| MEGAHIT only (`NFMAG_ASSEMBLER`) | one assembler avoids duplicate MAGs; faster; to confirm with supervisor |
| `binning_map_mode: group`, all samples group 0 | co-abundance from all samples improves binning |
| MetaBAT2 + MaxBin2, `skip_concoct` | same binners as the supervisor's metaWRAP script |
| `refine_bins_dastool`, `postbinning_input: refined_bins_only` | bin refinement, like metaWRAP `bin_refinement` |
| `binqc_tool: checkm` | CheckM, as in the supervisor's script |
| `skip_gtdbtk` | nf-core/mag 2.5.1 contains GTDB-Tk 2.3.2 (GTDB R214 only); GTDB-Tk runs separately with R226 |
| `skip_prokka` | annotation with Bakta instead |

After nf-core, `filter_mags.sh` keeps MAGs with ≥ 50 % completeness and
< 10 % contamination (as `metawrap bin_refinement -c 50 -x 10`).

### Mapping to the supervisor's scripts

| Supervisor script | Here |
| --- | --- |
| `Assembly_CRA007433.sh` (KneadData, repair.sh, metaSPAdes/MEGAHIT, metaWRAP binning + refinement, CheckM) | `03_preprocessing/kneaddata.sh` (KneadData + repair.sh), `05_assembly/run_nfcore_mag.sh` (nf-core/mag), `07_mags/filter_mags.sh` |
| `GTDB_CRA009581.sh` | `07_mags/gtdbtk.sh` (same command, `--skip_ani_screen`, 25 cores, R226) |
| `Bakta_CRA007433.sh` | `07_mags/bakta.sh` (same options, skips finished MAGs) |

## Comparison with CRASH

Adapted from CRASH (https://github.com/vikkxtar123/CRASH; Natarajan J., MSc
thesis, Lund University).

| Aspect | CRASH | This repository |
| --- | --- | --- |
| Central `config.sh` sourced by every script | yes | yes |
| `#SBATCH -A`, mail, `logs/%x_%j` in headers; plain `sbatch` | yes | yes |
| One conda env per tool (YAML files) | `Environments/` | `envs/` |
| KneadData → MetaPhlAn 4 (`mpa_vJun23_CHOCOPhlAnSGB_202403`, `--offline`) | yes | yes (MetaPhlAn ≥ 4.2 flag names detected) |
| Node-local scratch `$SNIC_TMP`, verify outputs before copying back | yes | yes |
| Parallelism | one node, worker farm | SLURM job arrays (one task per sample) |
| Sample selection from metadata | — | yes |
| Assembly / MAGs | — | nf-core/mag + GTDB-Tk + Bakta |

## Values still to obtain/check on COSMOS

Run `bash setup/check_config.sh` to see what is missing.

| Variable | What | Where to look |
| --- | --- | --- |
| `KNEADDATA_DB` | human Bowtie2 index | ask supervisor for a readable copy (his script uses one in his home), or `kneaddata_database --download human_genome bowtie2 <dir>` |
| `METAPHLAN_DB_DIR` | MetaPhlAn 4 database | shared group copy, or `metaphlan --install` |
| `ENV_GTDBTK`, `ENV_BAKTA` | GTDB-Tk / Bakta environment | supervisor's `/lunarc/nobackup/projects/lu2025-12-46/util/gtdbtk2` (check access), or `envs/gtdbtk.yml`, `envs/bakta.yml` |
| `GTDBTK_DATA_PATH` | GTDB release 226 | supervisor's `.../lu2024-12-47/util/GTDB/release226` (check access) |
| `BAKTA_DB` | Bakta database | supervisor's `.../lu2025-12-46/util/util_12-47/BAKTA/db` (check access) |
| `NFMAG_MAX_CPUS`, `NFMAG_MAX_MEMORY` | largest single nf-core task | `sinfo -o "%P %c %m"` |

Already set: account `lu2026-2-117`, e-mail, `Miniforge3/25.3.0-3`.

## Notes

- HMP removed human DNA before uploading to SRA; KneadData is still run for
  trimming, uniform QC and as a check. Human reads are never stored.
- An HMP sample (SRS) can have several runs (SRR); runs are merged per sample in 03.
- Three samples are shallow (< 2 Gbases): fine for taxonomy, few MAGs expected.
- Record exact tool versions: `conda env export -n <env> > docs/<env>.lock.yml`.
