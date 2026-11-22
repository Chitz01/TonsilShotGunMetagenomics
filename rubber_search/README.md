# Rubber-degrading enzyme search in tonsil and oral MAGs

Searches all proteins of the 1,253 quality-filtered MAGs (180 tonsil, 1,073 oral;
2,119,874 proteins) for homologs of characterised natural-rubber oxygenases
(Lcp, RoxA, RoxB): by Bakta annotation, by sequence (DIAMOND) and by protein family
(HMMER, Pfam), with a positive control, relaxed-threshold sensitivity searches and a
phylogenetic placement of the Rox-family hits.
Code drafted with help from Claude (Anthropic); declared per course rules.
Detailed parameters and reasoning: `METHODS_NOTES.md`.

## Files

| File | Where it runs | What it does |
|---|---|---|
| `reference_accessions.tsv` | input | UniProt IDs of the 6 reference enzymes |
| `pfam_models.txt` | input | Pfam models `MPAB_Lcp_cat` (PF09995), `RoxA-like_Cyt-c` (PF21419) |
| `01_get_references.sh` | login node | Downloads reference proteins (UniProt) and Pfam-A; extracts the 2 models |
| `01b_add_roxb_from_uniparc.sh` | login node | Adds RoxB 35Y from UniParc UPI0009C31814 (UniProt entry deleted) |
| `02_collect_bakta.sh` | SLURM | Combines `MAG.faa` of the quality-filtered MAGs (keep list = filtered GTDB-Tk summaries) into one FASTA (`>MAG\|locus_tag`); scans Bakta `.tsv` for rubber terms |
| `03_search_diamond_hmmer.sh` | SLURM | DIAMOND blastp (`--very-sensitive`, e ≤ 1e-5) and HMMER (`--cut_ga`) |
| `03b_sensitivity_search.sh` | SLURM | Relaxed searches: DIAMOND `--ultra-sensitive` e ≤ 10; HMMER `-E 10`, no GA |
| `04_summarise_hits.R` | R | Joins hits with GTDB-Tk taxonomy, MAG quality and body site (needed when there are DIAMOND hits) |
| `05_get_relatives.sh` | login node | Downloads Swiss-Prot cytochrome c peroxidases and MauG (tree comparison group) |
| `06_rox_tree.sh` | SLURM | Rox-model hits ≥ 300 aa + stringent hits → SeqKit dedup → MAFFT L-INS-i → trimAl `-gt 0.5` → IQ-TREE (`-m MFP -B 1000 -alrt 1000`) |
| `07_rox_tree_summary.R` | R (local) | Tree plot, nearest characterised reference per MAG protein, Rox-clade test (MauG or midpoint root) |
| `references/reference_proteins.faa` | data | The 6 reference sequences actually searched |

## Setup (once, on COSMOS)

```bash
module load Anaconda3/2024.02-1
conda create -p /lunarc/nobackup/projects/lu2025-12-46/CHETANA/envs/rubber_search \
    -c conda-forge -c bioconda diamond hmmer mafft trimal iqtree seqkit -y
```

R packages for script 07: `install.packages(c("ape", "phangorn"))`.

## Run order

1. `bash 01_get_references.sh`, then `bash 01b_add_roxb_from_uniparc.sh` (login node)
2. `sbatch 02_collect_bakta.sh` (expect: 1,253 MAGs, 2,119,874 proteins, 0 duplicated headers)
3. `sbatch 03_search_diamond_hmmer.sh`
4. `sbatch 03b_sensitivity_search.sh`
5. `bash 05_get_relatives.sh` (login node)
6. `sbatch 06_rox_tree.sh`
7. Copy `rubber_search/tree/` locally, then in R: `source("07_rox_tree_summary.R")`

Working folder on COSMOS: `/lunarc/nobackup/projects/lu2025-12-46/CHETANA/ANALYSIS/rubber_search`.
SLURM account: `lu2026-2-117`.

## Results so far (20–22 Nov 2026)

- Annotation scan: 74 proteins named cutinase; none named Lcp, Rox or Dsz.
- DIAMOND (e ≤ 1e-5): 0 hits.
- HMMER (Pfam GA): 1 hit, `MEGAHIT-MetaBAT2Refined-SRS144380.47|MNDINJ_01890` (126 aa, RoxA-like_Cyt-c, score 28.5 vs GA 26.7, model 1–31/39).
- Positive control: Lcp references 50.3–52.6 bits (≈ 90 % of model), Rox references 70.9–73.1 bits (full model).
- Relaxed: 60 DIAMOND hits (e 0.4–3.6), 213 HMMER proteins (156 Rox model, 57 Lcp model).
- Tree: 120 tips (106 unique MAG sequences = 133 proteins, 3 Rox, 7 CCP, 4 MauG); summary pending.
