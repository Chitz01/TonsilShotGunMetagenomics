## 00_combine_metaphlan_tonsil_oral.R
## Join the tonsil and oral MetaPhlAn species tables into one table (species
## missing from one dataset get 0), in MetaPhlAn merged format, for
## 03_MetaPhlAn4_phyloseq_tonsil_oral.R.

tonsil_file <- file.path("data", "hmp_tonsil", "species_abundance.tsv")
oral_file   <- file.path("data", "hmp_oral", "species_abundance.tsv")
out_file    <- file.path("data", "hmp_tonsil_oral_species_abundance.tsv")

read_mpa <- function(file) {
  first <- readLines(file, n = 20)
  hdr <- grep("clade_name", first)[1]
  tab <- read.delim(file, skip = hdr - 1, check.names = FALSE, stringsAsFactors = FALSE)
  names(tab)[1] <- "clade_name"
  tab[, !names(tab) %in% c("clade_taxid", "NCBI_tax_id"), drop = FALSE]
}

tonsil <- read_mpa(tonsil_file)
oral   <- read_mpa(oral_file)

shared <- intersect(names(tonsil)[-1], names(oral)[-1])
if (length(shared)) stop("Sample columns present in both tables: ", paste(shared, collapse = ", "))

combined <- merge(tonsil, oral, by = "clade_name", all = TRUE)
combined[is.na(combined)] <- 0

cat("Tonsil:", ncol(tonsil) - 1, "samples,", nrow(tonsil), "rows\n")
cat("Oral:  ", ncol(oral) - 1, "samples,", nrow(oral), "rows\n")
cat("Combined:", ncol(combined) - 1, "samples,", nrow(combined), "rows\n")

# Keep MetaPhlAn's database line on top, as in the original merged tables
mpa_line <- grep("^#mpa", readLines(tonsil_file, n = 5), value = TRUE)[1]
if (file.exists(out_file)) file.remove(out_file)
if (!is.na(mpa_line)) writeLines(mpa_line, out_file)
suppressWarnings(write.table(combined, out_file, sep = "\t", quote = FALSE,
                             row.names = FALSE, append = TRUE))
cat("Written:", out_file, "\n")
