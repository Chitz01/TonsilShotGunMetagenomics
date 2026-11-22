#!/bin/bash
#SBATCH -t 02:00:00
#SBATCH -J collect_bakta
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -A lu2026-2-117
## 02_collect_bakta.sh
## (1) Combine all Bakta protein files (.faa) of the tonsil + oral MAGs into one FASTA,
##     with every header prefixed by its MAG name:  >MAGNAME|LOCUS_TAG
## (2) First-pass annotation scan of the Bakta .tsv files for rubber/tyre-related terms.
## Code drafted with help from Claude (Anthropic); declare per course rules.
set -euo pipefail

## ---- edit these -------------------------------------------------------------
WORK="/lunarc/nobackup/projects/lu2025-12-46/CHETANA/ANALYSIS/rubber_search"
# Bakta output folders: one subfolder per MAG, each containing MAG.faa and MAG.tsv
BAKTA_DIRS=(
  "/lunarc/nobackup/projects/lu2025-12-46/CHETANA/ANALYSIS/hmp_tonsil/07_mags/bakta"
  "/lunarc/nobackup/projects/lu2025-12-46/CHETANA/ANALYSIS/hmp_oral/07_mags/bakta"
)
# Case-insensitive terms searched in the Bakta product names and database cross-references
SCAN_TERMS="EC:1\.13\.11\.85|EC:1\.13\.11\.87|latex.clearing|rubber|polyisoprene|DszA|DszB|DszC|dibenzothiophene|cutinase|PET hydrolase"
# Only MAGs that passed the quality filter (completeness > 50 %, contamination <= 5 %)
# are searched: their names are taken from the filtered GTDB-Tk summaries.
KEEP_FILES=(
  "/lunarc/nobackup/projects/lu2025-12-46/CHETANA/ANALYSIS/hmp_tonsil/07_mags/gtdbtk/gtdbtk.bac120.summary.filtered.tsv"
  "/lunarc/nobackup/projects/lu2025-12-46/CHETANA/ANALYSIS/hmp_oral/07_mags/gtdbtk/gtdbtk.bac120.summary.filtered.tsv"
)
## -----------------------------------------------------------------------------

mkdir -p "$WORK/proteins" "$WORK/annotation_scan"
KEEP="$WORK/proteins/MAGs_kept.txt"
tail -q -n +2 "${KEEP_FILES[@]}" | cut -f1 | sort -u > "$KEEP"
echo "MAGs in filtered GTDB-Tk lists: $(wc -l < "$KEEP")"
OUT="$WORK/proteins/all_MAG_proteins.faa"
: > "$OUT"
printf "MAG\tdataset\tn_proteins\n" > "$WORK/proteins/MAG_protein_counts.tsv"
printf "MAG\tdataset\tSequence_Id\tType\tStart\tStop\tStrand\tLocus_Tag\tGene\tProduct\tDbXrefs\n" \
    > "$WORK/annotation_scan/bakta_term_hits.tsv"

for d in "${BAKTA_DIRS[@]}"; do
    dataset=$(echo "$d" | grep -o 'hmp_[a-z]*' || echo "unknown")
    # One subfolder per MAG; use only MAG/MAG.faa (Bakta also writes MAG.hypotheticals.faa,
    # which repeats a subset of the same proteins and must NOT be included)
    n_found=$(for m in "$d"/*/; do b=$(basename "$m"); [ -f "$m/$b.faa" ] && echo x; done | wc -l)
    echo "$d: $n_found MAG protein files (MAG/MAG.faa)"
    [ "$n_found" -gt 0 ] || { echo "ERROR: no MAG/MAG.faa files in $d (check BAKTA_DIRS)"; exit 1; }

    for magdir in "$d"/*/; do
        mag=$(basename "$magdir")
        grep -qxF "$mag" "$KEEP" || continue      # not in the quality-filtered set
        faa="$magdir/$mag.faa"
        [ -f "$faa" ] || { echo "WARNING: no $mag.faa in $magdir, skipped"; continue; }
        # Prefix headers with MAG name; keep only the first word (locus tag) of the original header
        awk -v m="$mag" '/^>/{sub(/^>/,""); split($0,a," "); print ">" m "|" a[1]; next} {print}' "$faa" >> "$OUT"
        printf "%s\t%s\t%s\n" "$mag" "$dataset" "$(grep -c '^>' "$faa")" >> "$WORK/proteins/MAG_protein_counts.tsv"

        tsv="$magdir/$mag.tsv"
        if [ -f "$tsv" ]; then
            grep -v '^#' "$tsv" | grep -i -E "$SCAN_TERMS" \
              | awk -v m="$mag" -v ds="$dataset" 'BEGIN{OFS="\t"} {print m, ds, $0}' \
              >> "$WORK/annotation_scan/bakta_term_hits.tsv" || true
        fi
    done
done

echo "MAGs:      $(($(wc -l < "$WORK/proteins/MAG_protein_counts.tsv") - 1))"
echo "Proteins:  $(grep -c '^>' "$OUT")"
echo "Duplicated headers (should be 0): $(grep '^>' "$OUT" | sort | uniq -d | wc -l)"
echo "Annotation-scan hits: $(($(wc -l < "$WORK/annotation_scan/bakta_term_hits.tsv") - 1))"
