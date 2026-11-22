#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH -J rox_tree
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -A lu2026-2-117
## 06_rox_tree.sh
## Phylogenetic placement of MAG proteins that matched the RoxA-like_Cyt-c Pfam model:
##   candidates (MAG proteins) + RoxA/RoxB references + Swiss-Prot CCP and MauG relatives
##   -> remove identical duplicates -> MAFFT alignment -> trimAl trimming -> IQ-TREE ML tree
## Code drafted with help from Claude (Anthropic); declare per course rules.
set -euo pipefail

## ---- settings ------------------------------------------------------------------
W="/lunarc/nobackup/projects/lu2025-12-46/CHETANA/ANALYSIS/rubber_search"
ENV="/lunarc/nobackup/projects/lu2025-12-46/CHETANA/envs/rubber_search"
THREADS=20
MINLEN=300          # candidates shorter than this (fragments) are left out, except stringent hits
UFBOOT=1000         # ultrafast bootstrap replicates
## ---------------------------------------------------------------------------------

module load Anaconda3/2024.02-1
source activate "$ENV"
IQ=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
T="$W/tree"; mkdir -p "$T"; cd "$T"

# 1. Candidate IDs: relaxed Rox-model hits >= MINLEN aa, plus every stringent (GA) hit
grep -v '^#' "$W/sensitivity/hmmer_noGA_E10.domtbl" | \
  awk -v m="$MINLEN" '$4=="RoxA-like_Cyt-c" && $3>=m {print $1}' | sort -u > candidate_ids.txt
grep -v '^#' "$W/hmmer/rubber_models_vs_MAGs.tbl" | awk '$3=="RoxA-like_Cyt-c"{print $1}' >> candidate_ids.txt
sort -u -o candidate_ids.txt candidate_ids.txt
echo "Candidate MAG proteins: $(wc -l < candidate_ids.txt)"

# 2. Extract sequences; collapse identical sequences (same protein recovered in several MAGs)
seqkit grep -f candidate_ids.txt "$W/proteins/all_MAG_proteins.faa" > candidates.faa
seqkit rmdup -s -D duplicates.txt candidates.faa > candidates_unique.faa 2> rmdup.log
echo "Unique candidate sequences: $(grep -c '>' candidates_unique.faa)"

# 3. Combine with Rox references and curated relatives
seqkit grep -r -p '\|Rox' "$W/references/reference_proteins.faa" > rox_references.faa
[ -s relatives.faa ] || { echo "ERROR: relatives.faa missing - run 05_get_relatives.sh first"; exit 1; }
cat rox_references.faa relatives.faa candidates_unique.faa | sed '/^>/ s/|/__/g' > tree_input.faa
N=$(grep -c '>' tree_input.faa)
echo "Sequences in tree: $N (Rox refs $(grep -c '>' rox_references.faa), relatives $(grep -c '>' relatives.faa))"

# 4. Alignment: L-INS-i (most accurate) for <= 200 sequences, otherwise MAFFT --auto
if [ "$N" -le 200 ]; then
  mafft --localpair --maxiterate 1000 --thread "$THREADS" tree_input.faa > tree_input.aln 2> mafft.log
else
  mafft --auto --thread "$THREADS" tree_input.faa > tree_input.aln 2> mafft.log
fi

# 5. Trim gappy columns: keep columns where at least 50 % of sequences have a residue
#    (-gt 0.5). Candidates differ in length (fragments to ~1000 aa), so a fixed gap
#    threshold is more predictable than trimAl's automatic heuristics.
trimal -in tree_input.aln -out tree_input.trim.aln -gt 0.5
echo "Alignment columns: $(seqkit fx2tab -n -l tree_input.aln | head -1 | cut -f2) -> after trimming: $(seqkit fx2tab -n -l tree_input.trim.aln | head -1 | cut -f2)"

# 6. Maximum-likelihood tree: ModelFinder (-m MFP), ultrafast bootstrap + SH-aLRT support
"$IQ" -s tree_input.trim.aln -m MFP -B "$UFBOOT" -alrt 1000 -T AUTO -ntmax "$THREADS" --prefix rox_tree -redo > iqtree_stdout.log
echo "Best model: $(grep 'Best-fit model' rox_tree.iqtree | head -1)"

# 7. Annotation table for the R summary (label, type, dataset, genus, product, length)
ANN="$W/sensitivity/hmmer_noGA_hits_annotated.tsv"
{
  printf "label\ttype\tdataset\tgenus\tproduct\tlength\n"
  seqkit fx2tab -n -l tree_input.faa | while IFS=$'\t' read -r lab len; do
    lab=$(echo "$lab" | awk '{print $1}')
    case "$lab" in
      Rox*)      printf "%s\tRox_reference\t-\t-\t-\t%s\n" "$lab" "$len" ;;
      REF_CCP*)  printf "%s\tCCP_SwissProt\t-\t-\t-\t%s\n" "$lab" "$len" ;;
      REF_MauG*) printf "%s\tMauG_SwissProt\t-\t-\t-\t%s\n" "$lab" "$len" ;;
      *) orig=$(echo "$lab" | sed 's/__/|/g')
         row=$(awk -F'\t' -v p="$orig" '$1==p' "$ANN" | head -1)
         ds=$(echo "$row" | cut -f2); genus=$(echo "$row" | cut -f6); prod=$(echo "$row" | cut -f7)
         printf "%s\tMAG\t%s\t%s\t%s\t%s\n" "$lab" "${ds:--}" "${genus:--}" "${prod:--}" "$len" ;;
    esac
  done
} > tree_annotation.tsv
echo "Done. Tree: $T/rox_tree.treefile   Annotation: $T/tree_annotation.tsv"
