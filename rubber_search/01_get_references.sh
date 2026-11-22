#!/bin/bash
## 01_get_references.sh
## Download reference protein sequences (UniProt) and Pfam HMMs for the rubber-enzyme search.
## Run on the COSMOS LOGIN node (compute nodes usually have no internet):
##     bash 01_get_references.sh
## Code drafted with help from Claude (Anthropic); declare per course rules.
set -euo pipefail

## ---- edit these -------------------------------------------------------------
WORK="/lunarc/nobackup/projects/lu2025-12-46/CHETANA/ANALYSIS/rubber_search"
ENV="/lunarc/nobackup/projects/lu2025-12-46/CHETANA/envs/rubber_search"   # conda env with hmmer
SCRIPTS="$(cd "$(dirname "$0")" && pwd)"   # folder containing this script and the .tsv/.txt lists
## -----------------------------------------------------------------------------

module load Anaconda3/2024.02-1
source activate "$ENV"

mkdir -p "$WORK/references" "$WORK/hmm"
cd "$WORK/references"

## 1. Reference proteins from UniProt (one FASTA per accession, then combined)
: > reference_proteins.faa
tail -n +2 "$SCRIPTS/reference_accessions.tsv" | while IFS=$'\t' read -r acc name family organism; do
    echo "Downloading $acc ($name)"
    curl -sf "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "${acc}.fasta" \
        || { echo "  FAILED: $acc"; continue; }
    # Rename header to: >NAME|FAMILY|ACCESSION  (keeps the search output readable)
    awk -v n="$name" -v f="$family" -v a="$acc" 'NR==1{print ">" n "|" f "|" a; next} {print}' \
        "${acc}.fasta" >> reference_proteins.faa
done
echo "Reference sequences: $(grep -c '>' reference_proteins.faa)"

## 2. Pfam HMMs: download Pfam-A once, then extract the models we need by name
cd "$WORK/hmm"
if [ ! -f Pfam-A.hmm ]; then
    curl -sf -O https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip Pfam-A.hmm.gz
fi
[ -f Pfam-A.hmm.ssi ] || hmmfetch --index Pfam-A.hmm
hmmfetch -f Pfam-A.hmm "$SCRIPTS/pfam_models.txt" > rubber_models.hmm
echo "HMM models extracted: $(grep -c '^NAME' rubber_models.hmm) (expected $(wc -l < "$SCRIPTS/pfam_models.txt"))"
grep '^NAME\|^ACC\|^GA' rubber_models.hmm
