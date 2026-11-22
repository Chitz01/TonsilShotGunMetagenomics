## 04_summarise_hits.R
## Combine DIAMOND and HMMER hits with MAG taxonomy (GTDB-Tk), MAG quality and
## sample metadata (body site, subject). Writes all hits plus a filtered set.
## Run in R (locally after copying the files, or on COSMOS):
##     source("04_summarise_hits.R")
## Code drafted with help from Claude (Anthropic); declare per course rules.

## ---- edit these -------------------------------------------------------------
work <- "rubber_search"            # folder with diamond/, hmmer/ (copied from COSMOS)
gtdb_files <- c(                   # GTDB-Tk summaries, tonsil + oral
  "hmp_tonsil/07_mags/gtdbtk/gtdbtk.bac120.summary.tsv",
  "hmp_oral/07_mags/gtdbtk/gtdbtk.bac120.summary.tsv")
quality_files <- c(                # MAG quality summaries, tonsil + oral
  "hmp_tonsil/07_mags/mag_quality_summary.tsv",
  "hmp_oral/07_mags/mag_quality_summary.tsv")
quality_id_col   <- "bin"          # column with the MAG name in the quality files (check!)
completeness_col <- "completeness" # (check the exact column names in your files)
contamination_col <- "contamination"
metadata_file <- "metadata/master_sample_table.tsv"

# Filters for the "filtered" output. These are choices, not facts:
# decide them with your supervisor and report them in the Methods.
min_pident   <- 30    # DIAMOND: % identity to the reference
min_qcov     <- 70    # DIAMOND: % of the reference protein covered by the alignment
max_evalue   <- 1e-5
## -----------------------------------------------------------------------------

split_hit <- function(x) {          # "MAGNAME|LOCUS" -> MAG, locus
  data.frame(MAG = sub("\\|.*$", "", x), locus_tag = sub("^[^|]*\\|", "", x))
}

## DIAMOND
dcols <- c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend",
           "sstart", "send", "evalue", "bitscore", "qcovhsp", "scovhsp")
dfile <- file.path(work, "diamond", "reference_vs_MAGs.tsv")
dia <- if (file.exists(dfile) && file.size(dfile) > 0) {
  read.delim(dfile, header = FALSE, col.names = dcols)
} else data.frame(matrix(ncol = length(dcols), nrow = 0, dimnames = list(NULL, dcols)))
dia <- cbind(dia, split_hit(dia$sseqid))
dia$reference <- sub("\\|.*$", "", dia$qseqid)
dia$enzyme_family <- sapply(strsplit(dia$qseqid, "|", fixed = TRUE), `[`, 2)
dia$method <- "DIAMOND"

## HMMER domain table (whitespace-separated; columns 1-22 are fixed)
hfile <- file.path(work, "hmmer", "rubber_models_vs_MAGs.domtbl")
hl <- if (file.exists(hfile)) grep("^#", readLines(hfile), value = TRUE, invert = TRUE) else character(0)
hmm <- if (length(hl)) {
  f <- strsplit(trimws(hl), "\\s+")
  data.frame(target = sapply(f, `[`, 1), tlen = as.numeric(sapply(f, `[`, 3)),
             model = sapply(f, `[`, 4), full_evalue = as.numeric(sapply(f, `[`, 7)),
             full_score = as.numeric(sapply(f, `[`, 8)), dom_ievalue = as.numeric(sapply(f, `[`, 13)),
             dom_score = as.numeric(sapply(f, `[`, 14)),
             ali_from = as.numeric(sapply(f, `[`, 18)), ali_to = as.numeric(sapply(f, `[`, 19)))
} else data.frame(target = character(0), tlen = numeric(0), model = character(0),
                  full_evalue = numeric(0), full_score = numeric(0), dom_ievalue = numeric(0),
                  dom_score = numeric(0), ali_from = numeric(0), ali_to = numeric(0))
hmm <- cbind(hmm, split_hit(hmm$target))
hmm$method <- "HMMER"

## MAG annotation: taxonomy, quality, sample metadata
gtdb <- do.call(rbind, lapply(gtdb_files[file.exists(gtdb_files)], function(f)
  read.delim(f, check.names = FALSE)[, c("user_genome", "classification")]))
names(gtdb) <- c("MAG", "gtdb_classification")
if (is.null(gtdb)) stop("No GTDB-Tk summary files found: check gtdb_files")
gtdb$genus   <- sub(";.*$", "", sub("^.*g__", "", gtdb$gtdb_classification))
gtdb$species <- sub("^.*s__", "", gtdb$gtdb_classification)

qual <- do.call(rbind, lapply(quality_files[file.exists(quality_files)], function(f) {
  q <- read.delim(f, check.names = FALSE)
  miss <- setdiff(c(quality_id_col, completeness_col, contamination_col), names(q))
  if (length(miss)) stop("Column(s) not found in ", f, ": ", paste(miss, collapse = ", "),
                         "\nAvailable: ", paste(names(q), collapse = ", "))
  data.frame(MAG = sub("\\.(fa|fasta|fna)$", "", q[[quality_id_col]]),
             completeness = q[[completeness_col]], contamination = q[[contamination_col]])
}))

meta <- read.delim(metadata_file, check.names = FALSE)[, c("sample_id", "subject_id", "body_site", "dataset")]

annotate <- function(d) {
  if (!nrow(d)) return(d)
  d$sample_id <- ifelse(grepl("SRS[0-9]+", d$MAG),
                        sub("^.*?(SRS[0-9]+).*$", "\\1", d$MAG, perl = TRUE), NA)
  d <- merge(d, gtdb, by = "MAG", all.x = TRUE)
  if (!is.null(qual)) d <- merge(d, qual, by = "MAG", all.x = TRUE)
  merge(d, meta, by = "sample_id", all.x = TRUE)
}
dia <- annotate(dia)
hmm <- annotate(hmm)

if (nrow(dia) && any(is.na(dia$sample_id)))
  warning("Some MAG names contain no SRS sample ID; body site could not be assigned for them.")

## Filtered sets
dia_f <- dia[dia$pident >= min_pident & dia$qcovhsp >= min_qcov & dia$evalue <= max_evalue, ]
# HMMER already used Pfam gathering thresholds (--cut_ga), so all hits are kept

## Best DIAMOND hit per MAG protein (a protein can match several references)
best <- if (nrow(dia_f)) {
  dia_f <- dia_f[order(dia_f$sseqid, -dia_f$bitscore), ]
  dia_f[!duplicated(dia_f$sseqid), ]
} else dia_f

## Write outputs
out <- file.path(work, "summary")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
write.table(dia,  file.path(out, "diamond_hits_all_annotated.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(best, file.path(out, "diamond_hits_filtered_best_per_protein.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(hmm,  file.path(out, "hmmer_hits_annotated.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## Proteins found by both methods
both <- intersect(best$sseqid, hmm$target)
cat("DIAMOND hits (all):", nrow(dia), "| filtered (best per protein):", nrow(best), "\n")
cat("HMMER protein-domain hits:", nrow(hmm), "| unique proteins:", length(unique(hmm$target)), "\n")
cat("Proteins found by both methods:", length(both), "\n")

## Counts per enzyme family x body site (number of distinct MAGs carrying a hit)
count_mags <- function(d, family_col) {
  if (!nrow(d)) return(data.frame())
  u <- unique(d[, c(family_col, "body_site", "MAG")])
  tab <- as.data.frame(table(u[[family_col]], u$body_site, useNA = "ifany"))
  tab[tab$Freq > 0, ]
}
cnt_dia <- count_mags(best, "enzyme_family")
cnt_hmm <- count_mags(hmm, "model")
if (nrow(cnt_dia)) { names(cnt_dia) <- c("enzyme_family", "body_site", "n_MAGs")
  write.table(cnt_dia, file.path(out, "diamond_MAGs_per_family_and_site.tsv"), sep = "\t", quote = FALSE, row.names = FALSE) }
if (nrow(cnt_hmm)) { names(cnt_hmm) <- c("pfam_model", "body_site", "n_MAGs")
  write.table(cnt_hmm, file.path(out, "hmmer_MAGs_per_model_and_site.tsv"), sep = "\t", quote = FALSE, row.names = FALSE) }
writeLines(capture.output(sessionInfo()), file.path(out, "R_session_info.txt"))
cat("Outputs in:", normalizePath(out), "\n")
