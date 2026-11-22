#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -J rubber_search
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -A lu2026-2-117
## 03_search_diamond_hmmer.sh
## Search all MAG proteins for rubber-degrading enzyme homologs:
##   (a) DIAMOND blastp: reference proteins (Lcp, RoxA, RoxB) vs all MAG proteins
##   (b) HMMER hmmsearch: Pfam models (MPAB_Lcp_cat, RoxA-like_Cyt-c) vs all MAG proteins
## Thresholds here are deliberately permissive; filtering is done afterwards in R
## (04_summarise_hits.R) so the cut-offs can be decided with the supervisor.
## Code drafted with help from Claude (Anthropic); declare per course rules.
set -euo pipefail

## ---- edit these -------------------------------------------------------------
WORK="/lunarc/nobackup/projects/lu2025-12-46/CHETANA/ANALYSIS/rubber_search"
ENV="/lunarc/nobackup/projects/lu2025-12-46/CHETANA/envs/rubber_search"   # diamond + hmmer
THREADS=20
EVALUE=1e-5
## -----------------------------------------------------------------------------

module load Anaconda3/2024.02-1
source activate "$ENV"

PROT="$WORK/proteins/all_MAG_proteins.faa"
REF="$WORK/references/reference_proteins.faa"
HMM="$WORK/hmm/rubber_models.hmm"
mkdir -p "$WORK/diamond" "$WORK/hmmer"

## (a) DIAMOND: database = all MAG proteins (large), query = references (few)
if [ ! -f "$WORK/diamond/all_MAG_proteins.dmnd" ]; then
    diamond makedb --in "$PROT" -d "$WORK/diamond/all_MAG_proteins" --threads "$THREADS"
fi
diamond blastp \
    -q "$REF" -d "$WORK/diamond/all_MAG_proteins" \
    -o "$WORK/diamond/reference_vs_MAGs.tsv" \
    --very-sensitive --evalue "$EVALUE" --max-target-seqs 100000 --threads "$THREADS" \
    --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp scovhsp
echo "DIAMOND hits: $(wc -l < "$WORK/diamond/reference_vs_MAGs.tsv")"

## (b) HMMER: Pfam gathering thresholds (--cut_ga) = Pfam's own curated family cut-offs
hmmsearch --cpu "$THREADS" --cut_ga \
    --domtblout "$WORK/hmmer/rubber_models_vs_MAGs.domtbl" \
    --tblout    "$WORK/hmmer/rubber_models_vs_MAGs.tbl" \
    "$HMM" "$PROT" > "$WORK/hmmer/hmmsearch.log"
echo "HMMER protein hits: $(grep -vc '^#' "$WORK/hmmer/rubber_models_vs_MAGs.tbl" || true)"
