#!/bin/bash
## 05_get_relatives.sh
## Download curated (Swiss-Prot, reviewed) bacterial relatives of the Rox enzymes:
## di-heme cytochrome c peroxidases (CCP) and MauG. They are the comparison group in the
## phylogenetic tree: a MAG protein that groups with RoxA/RoxB rather than with these
## relatives is a stronger Rox candidate.
## Run on the COSMOS LOGIN node (needs internet):   bash 05_get_relatives.sh
## Code drafted with help from Claude (Anthropic); declare per course rules.
set -euo pipefail

W="/lunarc/nobackup/projects/lu2025-12-46/CHETANA/ANALYSIS/rubber_search"
mkdir -p "$W/tree"
cd "$W/tree"

# fetch NAME QUERY : UniProt REST query -> FASTA with headers >REF_NAME|ENTRYNAME|ACCESSION
fetch() {
  curl -sfG "https://rest.uniprot.org/uniprotkb/stream" \
       --data-urlencode "query=$2" --data-urlencode "format=fasta" > "raw_$1.faa" \
    || { echo "FAILED to download $1" >&2; return 1; }
  echo "$1: $(grep -c '>' "raw_$1.faa") sequences" >&2
  awk -v t="$1" '/^>/{split($0,a,"|"); split(a[3],b," "); print ">REF_" t "|" b[1] "|" a[2]; next}{print}' "raw_$1.faa"
}

fetch CCP  '(reviewed:true) AND (taxonomy_id:2) AND (protein_name:cytochrome) AND (protein_name:peroxidase)' >  relatives.faa
fetch MauG '(reviewed:true) AND (taxonomy_id:2) AND (gene:mauG)'                                                >> relatives.faa

echo "Total relatives: $(grep -c '>' relatives.faa)"
echo "Downloaded entries (check that they are cytochrome c peroxidases / MauG):"
grep -h '>' raw_CCP.faa raw_MauG.faa | cut -c1-110
