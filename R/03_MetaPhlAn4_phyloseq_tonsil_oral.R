###############################################################################
# METAPHLAN 4 -> PHYLOSEQ -> SINGLE-STUDY MICROBIOME ANALYSIS
#
# Inputs:
#   1. MetaPhlAn *_merged_abundance_table.txt (clade_name + sample columns)
#   2. One matching metadata table (.xlsx, .csv or .tsv)
#
# Outputs: phyloseq RDS, sample/taxonomy audits, alpha and beta diversity,
#          group-level and sample-level abundance profiles, MaAsLin2 results.
#
# For another study, normally only edit SECTION 1 below.
# Features are species-level MetaPhlAn relative abundances, NOT raw read counts.
###############################################################################

# 0. PACKAGES -----------------------------------------------------------------
required_packages <- c(
  "data.table", "readxl", "dplyr", "tidyr", "tibble", "ggplot2",
  "phyloseq", "vegan", "scales", "Maaslin2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install missing packages: ", paste(missing_packages, collapse = ", "),
       "\nFor phyloseq and Maaslin2, use BiocManager::install().")
}
set.seed(42)

# 1. CONFIGURATION -------------------------------------------------------------
# HMP PRJNA48479: tonsil + oral sites from the same people (paired by subject_id).
# Paths are relative to the RStudio project folder.
# [Chetana] Edited from the supervisor's PRJNA1226574 version.
wdir <- "."
rdir <- file.path("results", "R", "hmp_tonsil_oral_SINGLE_STUDY")

# Combined tonsil + oral table, made by R/00_combine_metaphlan_tonsil_oral.R
abundance_file <- file.path(wdir, "data", "hmp_tonsil_oral_species_abundance.tsv")
metadata_file  <- file.path(wdir, "metadata", "master_sample_table.tsv")
metadata_sheet <- 1                      # ignored for TSV

sample_column <- "sample_id"             # SRS IDs, matching MetaPhlAn columns
group_column  <- "body_site"             # tonsil vs each oral site
reference_group <- "Palatine Tonsils"    # must match the spelling in body_site exactly
# Optional ordered vector if you want a specific order; NULL discovers groups.
group_levels <- NULL

# Optional confounders: ONLY add columns that actually exist in the metadata.
# log10_clean_pairs is created from clean_pairs in section 4 (sequencing depth).
# E.g.: covariates <- c("log10_clean_pairs", "sequencer")
covariates <- character(0)
numeric_covariates <- c("log10_clean_pairs", "clean_pairs")
# Same people sampled at several sites -> subject as random effect in MaAsLin2
random_effects <- "subject_id"

# [Chetana] Pairing for PERMANOVA/PERMDISP (permutations within subject) and for
# the alpha-diversity tests (paired Wilcoxon, each site vs reference_group).
# Set to NULL to get the original, unpaired behaviour.
pairing_column <- "subject_id"

BACTERIA_ONLY <- TRUE                    # Restrict all analyses to Bacteria
alpha_presence_threshold <- 1e-5        # Relative fraction; 1e-5 = 0.001%
jaccard_presence_threshold <- 0         # Nonzero bacterial abundance -> present
min_prevalence_maaslin <- 0.10
fdr_cutoff <- 0.05
n_permutations <- 9999
top_n_group <- 20
top_n_plot <- 20                     # Global top species for sample stacked bars
top_n_per_sample <- 10                # Export each sample's top species

# Optional: a SINGLE MetaPhlAn profile can have a relative_abundance column.
# Supply its sample ID here if needed. Merged tables do not require this.
single_sample_id <- NULL

# 2. HELPERS AND OUTPUT FOLDERS ------------------------------------------------
dir.create(rdir, showWarnings = FALSE, recursive = TRUE)
dirs <- setNames(file.path(rdir, c("00_QC", "01_PHYLOSEQ", "02_ALPHA",
                                     "03_BETA", "04_ABUNDANCE", "05_MAASLIN2")),
                 c("QC", "PS", "ALPHA", "BETA", "ABUND", "MAASLIN"))
invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

save_plot <- function(p, name, directory, width = 10, height = 6) {
  ggplot2::ggsave(file.path(directory, paste0(name, ".pdf")), p,
                  width = width, height = height, units = "in", limitsize = FALSE)
  ggplot2::ggsave(file.path(directory, paste0(name, ".png")), p,
                  width = width, height = height, units = "in", dpi = 400,
                  limitsize = FALSE)
}
write_tsv <- function(x, file) {
  data.table::fwrite(as.data.frame(x), file, sep = "\t", na = "NA")
}
base_theme <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey95"),
                   axis.text.x = ggplot2::element_text(colour = "black"),
                   axis.text.y = ggplot2::element_text(colour = "black"))
}
pretty_taxon <- function(z) gsub("_", " ", z, fixed = TRUE)

# 3. READ METAPHLAN AND KEEP NONREDUNDANT SPECIES ------------------------------
if (!file.exists(abundance_file)) stop("MetaPhlAn file not found: ", abundance_file)
if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)

# 'skip = "clade_name"' locates the header despite MetaPhlAn's #mpa_v... line.
mpa <- data.table::fread(abundance_file, sep = "\t", skip = "clade_name",
                         header = TRUE, data.table = FALSE, check.names = FALSE)
names(mpa)[1] <- sub("^#", "", names(mpa)[1])
if (!identical(names(mpa)[1], "clade_name")) {
  stop("Expected first MetaPhlAn column to be 'clade_name'.")
}

# Species rows end with s__; including ancestor rows would double-count taxa.
is_species <- grepl("(^|\\|)s__[^|]+$", mpa$clade_name, perl = TRUE)
sp <- mpa[is_species, , drop = FALSE]
if (!nrow(sp)) stop("No terminal s__ taxa found. Check MetaPhlAn format.")

lineages <- strsplit(sp$clade_name, "\\|", perl = TRUE)
get_rank <- function(lin, prefix) {
  y <- lin[startsWith(lin, prefix)]
  if (!length(y)) return(NA_character_)
  substring(y[1], nchar(prefix) + 1L)
}
rank_prefix <- c(Kingdom = "k__", Phylum = "p__", Class = "c__",
                 Order = "o__", Family = "f__", Genus = "g__", Species = "s__")
taxonomy <- as.data.frame(lapply(rank_prefix, function(prefix) {
  vapply(lineages, get_rank, FUN.VALUE = character(1), prefix = prefix)
}), stringsAsFactors = FALSE, check.names = FALSE)
if (BACTERIA_ONLY) {
  keep_bacteria <- !is.na(taxonomy$Kingdom) & taxonomy$Kingdom == "Bacteria"
  sp <- sp[keep_bacteria, , drop = FALSE]
  taxonomy <- taxonomy[keep_bacteria, , drop = FALSE]
}
if (!nrow(sp)) stop("No species remaining after the bacteria-only filter.")

# Feature IDs must be unambiguous and safe as MaAsLin2 variable names.
feature_id <- make.unique(make.names(taxonomy$Species), sep = "__duplicate_")
rownames(taxonomy) <- feature_id

nontaxon_columns <- c("clade_name", "clade_taxid", "NCBI_tax_id", "tax_id",
                      "additional_species", "estimated_number_of_reads_from_the_clade")
sample_columns <- setdiff(names(sp), nontaxon_columns)
if ("relative_abundance" %in% sample_columns) {
  if (length(sample_columns) != 1L || is.null(single_sample_id)) {
    stop("Found a single-profile 'relative_abundance' column. ",
         "Use a merged MetaPhlAn table or set single_sample_id.")
  }
  names(sp)[names(sp) == "relative_abundance"] <- single_sample_id
  sample_columns <- single_sample_id
}
if (!length(sample_columns)) stop("No sample columns found in MetaPhlAn input.")
if (anyDuplicated(sample_columns)) stop("Duplicated sample columns in MetaPhlAn input.")

x_original <- as.matrix(sp[, sample_columns, drop = FALSE])
suppressWarnings(storage.mode(x_original) <- "double")
rownames(x_original) <- feature_id
colnames(x_original) <- trimws(colnames(x_original))
# [Chetana] Column names come from profile file names: "SRS011086_profile" -> "SRS011086"
colnames(x_original) <- sub("(_profile)?(\\.txt|\\.tsv)?$", "", colnames(x_original))
if (anyNA(x_original) || any(!is.finite(x_original)) || any(x_original < 0)) {
  stop("MetaPhlAn abundances contain missing, nonnumeric, negative or infinite values.")
}
if (anyDuplicated(colnames(x_original))) stop("Duplicated sample IDs after trimming.")

# MetaPhlAn merged tables ordinarily contain percentages (0..100).
# Support already-normalized relative fractions (0..1), too.
raw_sums <- colSums(x_original)
scale_denom <- if (max(raw_sums) > 1.5) 100 else 1
if (max(raw_sums / scale_denom) > 1.01) {
  warning("Some species totals exceed 100%; check input units and hierarchy.")
}
x <- x_original / scale_denom

# 4. READ METADATA; ALIGN SAMPLES ---------------------------------------------
ext <- tolower(tools::file_ext(metadata_file))
metadata <- switch(ext,
  "xlsx" = as.data.frame(readxl::read_excel(metadata_file, sheet = metadata_sheet)),
  "xls"  = as.data.frame(readxl::read_excel(metadata_file, sheet = metadata_sheet)),
  "csv"  = utils::read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE),
  "tsv"  = utils::read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE),
  "txt"  = utils::read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE),
  stop("Metadata must be .xlsx, .xls, .csv, .tsv or tab-delimited .txt"))
# [Chetana] Sequencing depth on log scale, usable as covariate
if ("clean_pairs" %in% names(metadata)) {
  metadata$log10_clean_pairs <- log10(as.numeric(metadata$clean_pairs))
}
needed_columns <- c(sample_column, group_column, covariates, random_effects, pairing_column)
if (!all(needed_columns %in% names(metadata))) {
  miss <- setdiff(needed_columns, names(metadata))
  stop("Missing metadata columns: ", paste(miss, collapse = ", "))
}
metadata[[sample_column]] <- trimws(as.character(metadata[[sample_column]]))
if (anyNA(metadata[[sample_column]]) || any(!nzchar(metadata[[sample_column]]))) {
  stop("Metadata has missing or blank sample IDs.")
}
if (anyDuplicated(metadata[[sample_column]])) {
  stop("Metadata sample IDs are not unique. Resolve duplicate rows before analysis.")
}
rownames(metadata) <- metadata[[sample_column]]

common <- colnames(x)[colnames(x) %in% rownames(metadata)]
audit <- rbind(
  data.frame(Sample = colnames(x), source = "MetaPhlAn",
             matched = colnames(x) %in% rownames(metadata)),
  data.frame(Sample = rownames(metadata), source = "Metadata",
             matched = rownames(metadata) %in% colnames(x))
)
write_tsv(audit, file.path(dirs["QC"], "sample_matching_audit.tsv"))
if (length(common) < 4) {
  stop("Too few matched samples (", length(common), "). ",
       "Check the 00_QC/sample_matching_audit.tsv file and sample ID format.")
}
if (length(common) != ncol(x) || length(common) != nrow(metadata)) {
  warning("Some samples are unmatched; see sample_matching_audit.tsv")
}
x <- x[, common, drop = FALSE]
metadata <- metadata[common, , drop = FALSE]

metadata[[group_column]] <- trimws(as.character(metadata[[group_column]]))
metadata[[group_column]][!nzchar(metadata[[group_column]])] <- NA_character_
keep <- !is.na(metadata[[group_column]]) & colSums(x) > 0
excluded <- data.frame(Sample = colnames(x)[!keep],
                       reason = ifelse(is.na(metadata[[group_column]][!keep]),
                                       "missing group", "zero retained species abundance"))
write_tsv(excluded, file.path(dirs["QC"], "excluded_samples.tsv"))
x <- x[, keep, drop = FALSE]
metadata <- metadata[keep, , drop = FALSE]
if (ncol(x) < 4) stop("Fewer than four usable samples remain.")

seen_levels <- sort(unique(as.character(metadata[[group_column]])))
if (!is.null(group_levels)) {
  if (!setequal(group_levels, seen_levels)) {
    stop("group_levels do not exactly match the observed group labels: ",
         paste(seen_levels, collapse = ", "))
  }
} else {
  group_levels <- seen_levels
  if (!is.null(reference_group)) {
    if (!reference_group %in% seen_levels) {
      stop("Reference group not found: '", reference_group, "'. Observed: ",
           paste(seen_levels, collapse = ", "))
    }
    group_levels <- c(reference_group, setdiff(seen_levels, reference_group))
  }
}
if (!is.null(reference_group)) {
  group_levels <- c(reference_group, setdiff(group_levels, reference_group))
}
metadata[[group_column]] <- factor(metadata[[group_column]], levels = group_levels)
if (nlevels(droplevels(metadata[[group_column]])) < 2) {
  stop("The chosen grouping variable needs at least two observed groups.")
}
if (any(table(metadata[[group_column]]) < 2L)) {
  stop("Beta-diversity comparisons require at least two samples per group.")
}

# Remove the omitted/unclassified fraction and renormalize the retained species
# to sum to 1 PER SAMPLE; this is a conditional bacterial-species composition.
pre_normalization <- data.frame(
  Sample = colnames(x), Group = as.character(metadata[[group_column]]),
  retained_species_percent_of_original = colSums(x) * 100
)
write_tsv(pre_normalization, file.path(dirs["QC"], "retained_species_fraction.tsv"))
if (any(colSums(x) < 0.50)) {
  warning("At least one sample has <50% of its original profile retained ",
          "as chosen species; inspect retained_species_fraction.tsv")
}
x <- sweep(x, 2, colSums(x), FUN = "/")
x <- x[rowSums(x) > 0, , drop = FALSE]
taxonomy <- taxonomy[rownames(x), , drop = FALSE]
stopifnot(identical(colnames(x), rownames(metadata)),
          identical(rownames(x), rownames(taxonomy)),
          all(abs(colSums(x) - 1) < 1e-8))

feature_map <- data.frame(FeatureID = rownames(taxonomy),
                          Clade = sp$clade_name[match(rownames(taxonomy), feature_id)],
                          taxonomy, check.names = FALSE, row.names = NULL)
write_tsv(feature_map, file.path(dirs["QC"], "feature_taxonomy_map.tsv"))
write_tsv(data.frame(Group = names(table(metadata[[group_column]])),
                     N = as.integer(table(metadata[[group_column]]))),
          file.path(dirs["QC"], "samples_per_group.tsv"))

# 5. BUILD AND SAVE PHYLOSEQ ---------------------------------------------------
ps <- phyloseq::phyloseq(
  phyloseq::otu_table(x, taxa_are_rows = TRUE),
  phyloseq::tax_table(as.matrix(taxonomy)),
  phyloseq::sample_data(metadata)
)
saveRDS(ps, file.path(dirs["PS"], "phyloseq_species.rds"))
write_tsv(data.frame(FeatureID = rownames(x), x, check.names = FALSE),
          file.path(dirs["PS"], "species_relative_abundance.tsv"))
aligned_metadata_export <- metadata
aligned_metadata_export[[sample_column]] <- NULL
write_tsv(data.frame(Sample = rownames(metadata), aligned_metadata_export, check.names = FALSE),
          file.path(dirs["PS"], "aligned_metadata.tsv"))
message("Constructed ps: ", phyloseq::nsamples(ps), " samples; ",
        phyloseq::ntaxa(ps), " bacterial species/features")

sample_info <- data.frame(Sample = colnames(x),
                          Group = metadata[[group_column]])
colour_values <- if (length(group_levels) <= 5) {
  c("#4169A1", "#B85651", "#369578", "#9677AF", "#D29C3C")[seq_along(group_levels)]
} else {
  scales::hue_pal()(length(group_levels))
}
colours <- stats::setNames(colour_values, group_levels)

# 6. ALPHA DIVERSITY -----------------------------------------------------------
# Shannon/Simpson operate on relative abundances. Observed uses a configurable
# presence cutoff to reduce the influence of tiny MetaPhlAn assignments.
alpha <- phyloseq::estimate_richness(ps, measures = c("Shannon", "Simpson"))
alpha$Observed <- colSums(x[, rownames(alpha), drop = FALSE] > alpha_presence_threshold)
alpha <- tibble::rownames_to_column(as.data.frame(alpha), var = "Sample")
alpha <- dplyr::left_join(alpha, sample_info, by = "Sample")
write_tsv(alpha, file.path(dirs["ALPHA"], "alpha_diversity_per_sample.tsv"))

metrics <- c("Shannon", "Simpson", "Observed")
alpha_tests <- list()
alpha_pairwise <- list()
for (metric in metrics) {
  tmp <- data.frame(Value = alpha[[metric]], Group = alpha$Group)
  if (!is.null(pairing_column)) {
    # [Chetana] Paired design: Wilcoxon signed-rank, each site vs reference_group,
    # using only subjects sampled at both sites; BH across sites.
    tmp$Subject <- as.character(metadata[alpha$Sample, pairing_column])
    ref <- tmp[tmp$Group == group_levels[1], c("Subject", "Value")]
    tab <- do.call(rbind, lapply(group_levels[-1], function(g) {
      both <- merge(ref, tmp[tmp$Group == g, c("Subject", "Value")],
                    by = "Subject", suffixes = c("_ref", "_site"))
      p <- if (nrow(both) >= 2) {
        stats::wilcox.test(both$Value_site, both$Value_ref, paired = TRUE,
                           exact = FALSE)$p.value
      } else NA_real_
      data.frame(Group1 = group_levels[1], Group2 = g, n_pairs = nrow(both),
                 median_diff_site_minus_ref = stats::median(both$Value_site - both$Value_ref),
                 p = p)
    }))
    tab$p_BH <- stats::p.adjust(tab$p, "BH")
    tab$Metric <- metric
    alpha_pairwise[[metric]] <- tab
    p <- NA_real_
    method <- "Paired Wilcoxon signed-rank vs reference (see alpha_pairwise_BH.tsv)"
  } else if (nlevels(droplevels(tmp$Group)) == 2) {
    p <- stats::wilcox.test(Value ~ Group, data = tmp, exact = FALSE)$p.value
    method <- "Wilcoxon rank-sum"
  } else {
    p <- stats::kruskal.test(Value ~ Group, data = tmp)$p.value
    method <- "Kruskal-Wallis"
    pw <- stats::pairwise.wilcox.test(tmp$Value, tmp$Group,
                                     p.adjust.method = "BH", exact = FALSE)
    tab <- as.data.frame(as.table(pw$p.value), stringsAsFactors = FALSE)
    names(tab) <- c("Group1", "Group2", "p_BH")
    tab <- tab[!is.na(tab$p_BH), , drop = FALSE]
    tab$Metric <- metric
    alpha_pairwise[[metric]] <- tab
  }
  alpha_tests[[metric]] <- data.frame(Metric = metric, Method = method,
                                      pvalue = p)
  plot_data <- data.frame(Sample = alpha$Sample,
                          Group = alpha$Group, Value = alpha[[metric]])
  p_alpha <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = Value,
                                                      fill = Group, colour = Group)) +
    ggplot2::geom_boxplot(width = 0.55, alpha = 0.24, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.13, height = 0, size = 1.8, alpha = 0.65) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::labs(x = NULL, y = metric, title = paste("Alpha diversity —", metric)) +
    base_theme() + ggplot2::theme(legend.position = "none",
                                  axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
  save_plot(p_alpha, paste0("alpha_", metric), dirs["ALPHA"], width = 8, height = 5.5)
}
alpha_tests <- dplyr::bind_rows(alpha_tests)
alpha_tests$qvalue_BH_across_metrics <- stats::p.adjust(alpha_tests$pvalue, "BH")
write_tsv(alpha_tests, file.path(dirs["ALPHA"], "alpha_group_tests.tsv"))
if (length(alpha_pairwise)) {
  write_tsv(dplyr::bind_rows(alpha_pairwise),
            file.path(dirs["ALPHA"], "alpha_pairwise_BH.tsv"))
}

# 7. BETA DIVERSITY ------------------------------------------------------------
# vegan uses sample rows and species columns.
vegan_x <- t(x)
beta_meta <- data.frame(Group = sample_info$Group, row.names = sample_info$Sample)
stopifnot(identical(rownames(vegan_x), rownames(beta_meta)))

# [Chetana] Paired design: permute samples only within the same subject.
perm_design <- if (!is.null(pairing_column)) {
  permute::how(blocks = factor(metadata[rownames(beta_meta), pairing_column]),
               nperm = n_permutations)
} else n_permutations

run_beta <- function(distance, distance_name) {
  # adonis2: differences between group centroids, but also sensitive to dispersion.
  ad <- vegan::adonis2(distance ~ Group, data = beta_meta,
                       permutations = perm_design)
  ad_tbl <- data.frame(Term = rownames(as.data.frame(ad)),
                       as.data.frame(ad), check.names = FALSE, row.names = NULL)
  write_tsv(ad_tbl, file.path(dirs["BETA"],
                              paste0(distance_name, "_PERMANOVA.tsv")))

  # PERMDISP tests homogeneity of multivariate dispersion between groups.
  bd <- vegan::betadisper(distance, beta_meta$Group, type = "median",
                          bias.adjust = TRUE, add = "lingoes")
  pd <- vegan::permutest(bd, permutations = perm_design)
  pd_tbl <- data.frame(Term = rownames(pd$tab), pd$tab,
                       check.names = FALSE, row.names = NULL)
  write_tsv(pd_tbl, file.path(dirs["BETA"],
                              paste0(distance_name, "_PERMDISP.tsv")))
  disp <- data.frame(Sample = rownames(beta_meta), Group = beta_meta$Group,
                     DistanceToCentroid = as.numeric(bd$distances))
  write_tsv(disp, file.path(dirs["BETA"],
                             paste0(distance_name, "_distance_to_centroid.tsv")))

  pdisp <- ggplot2::ggplot(disp, ggplot2::aes(Group, DistanceToCentroid,
                                               colour = Group, fill = Group)) +
    ggplot2::geom_boxplot(alpha = .2, outlier.shape = NA) +
    ggplot2::geom_jitter(width = .1, size = 1.5, alpha = .6) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::labs(x = NULL, y = "Distance to group centroid",
                  title = paste(distance_name, "dispersion")) +
    base_theme() + ggplot2::theme(legend.position = "none")
  save_plot(pdisp, paste0(distance_name, "_PERMDISP_plot"),
            dirs["BETA"], width = 8, height = 5.5)

  # Lingoes correction makes non-Euclidean distances suitable for PCoA.
  ord <- vegan::wcmdscale(distance, eig = TRUE, k = 2, add = "lingoes")
  pts <- as.matrix(ord$points)
  if (ncol(pts) < 2) pts <- cbind(pts[, 1], PC2 = 0)
  eig <- pmax(ord$eig, 0)
  pct <- if (sum(eig) > 0) 100 * eig[1:2] / sum(eig) else c(NA_real_, NA_real_)
  coordinates <- data.frame(Sample = rownames(pts),
                            PC1 = pts[, 1], PC2 = pts[, 2])
  coordinates <- dplyr::left_join(coordinates, sample_info, by = "Sample")
  write_tsv(coordinates, file.path(dirs["BETA"],
                                    paste0(distance_name, "_PCoA_coordinates.tsv")))
  p <- ggplot2::ggplot(coordinates,
                       ggplot2::aes(x = PC1, y = PC2, colour = Group)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::labs(title = paste0(distance_name, " PCoA"),
                  x = sprintf("PCoA1 (%.1f%%)", pct[1]),
                  y = sprintf("PCoA2 (%.1f%%)", pct[2])) +
    base_theme()
  if (all(table(coordinates$Group) >= 3)) {
    p <- p + ggplot2::stat_ellipse(linewidth = .6, alpha = .65)
  }
  save_plot(p, paste0(distance_name, "_PCoA"), dirs["BETA"],
            width = 8, height = 6)
}

bray <- vegan::vegdist(vegan_x, method = "bray")
run_beta(bray, "Bray_Curtis")
pa_x <- (vegan_x > jaccard_presence_threshold) * 1L
if (any(rowSums(pa_x) == 0)) {
  stop("Jaccard threshold leaves at least one sample with zero species; lower it.")
}
jaccard <- vegan::vegdist(pa_x, method = "jaccard", binary = TRUE)
run_beta(jaccard, "Jaccard_binary")

# 8. MOST ABUNDANT SPECIES BY GROUP -------------------------------------------
sp_names <- stats::setNames(taxonomy$Species, rownames(taxonomy))
sp_labels <- pretty_taxon(sp_names)
# Keep readable but distinct labels if MetaPhlAn has duplicate species names.
duplicated_species <- duplicated(sp_labels) | duplicated(sp_labels, fromLast = TRUE)
sp_labels[duplicated_species] <- paste0(
  sp_labels[duplicated_species], " [", names(sp_labels)[duplicated_species], "]"
)
means <- do.call(rbind, lapply(group_levels, function(g) {
  ids <- sample_info$Sample[sample_info$Group == g]
  data.frame(Group = g, FeatureID = rownames(x),
             MeanRelativeAbundance = rowMeans(x[, ids, drop = FALSE]))
}))
means$Species <- sp_labels[means$FeatureID]
means$MeanPercent <- 100 * means$MeanRelativeAbundance
write_tsv(means, file.path(dirs["ABUND"], "all_species_mean_abundance_by_group.tsv"))

top_by_group <- means |>
  dplyr::filter(MeanRelativeAbundance > 0) |>
  dplyr::group_by(Group) |>
  dplyr::slice_max(MeanRelativeAbundance, n = top_n_group, with_ties = FALSE) |>
  dplyr::arrange(dplyr::desc(MeanRelativeAbundance), .by_group = TRUE) |>
  dplyr::mutate(Rank = dplyr::row_number()) |>
  dplyr::ungroup()
write_tsv(top_by_group, file.path(dirs["ABUND"], "top_species_by_group.tsv"))

# Separate labels per facet avoid incorrect species ranking in different groups.
top_plot_data <- top_by_group |>
  dplyr::mutate(FacetSpecies = paste(Species, Group, sep = "___"),
                FacetSpecies = stats::reorder(FacetSpecies, MeanPercent))
p_group <- ggplot2::ggplot(top_plot_data,
                           ggplot2::aes(x = MeanPercent, y = FacetSpecies, fill = Group)) +
  ggplot2::geom_col(width = 0.8) +
  ggplot2::facet_wrap(~ Group, scales = "free_y") +
  ggplot2::scale_y_discrete(labels = function(z) sub("___.*$", "", z)) +
  ggplot2::scale_fill_manual(values = colours) +
  ggplot2::labs(y = NULL, x = "Mean relative abundance (%)",
                title = paste("Top", top_n_group, "species per group")) +
  base_theme() + ggplot2::theme(legend.position = "none")
save_plot(p_group, "top_species_per_group", dirs["ABUND"],
          width = max(10, length(group_levels) * 5), height = 8)

# 9. MOST ABUNDANT SPECIES PER SAMPLE -----------------------------------------
sample_top <- do.call(rbind, lapply(seq_len(ncol(x)), function(j) {
  idx <- head(order(x[, j], decreasing = TRUE), top_n_per_sample)
  idx <- idx[x[idx, j] > 0]
  data.frame(Sample = colnames(x)[j], Group = as.character(sample_info$Group[j]),
             Rank = seq_along(idx), FeatureID = rownames(x)[idx],
             Species = sp_labels[rownames(x)[idx]],
             RelativeAbundance = x[idx, j], RelativePercent = 100 * x[idx, j])
}))
write_tsv(sample_top, file.path(dirs["ABUND"], "top_species_per_sample.tsv"))

# Stacked bars show the same global top N across all samples; all other
# identified species are collapsed into the 'Other' category.
global_top <- head(names(sort(rowMeans(x), decreasing = TRUE)), top_n_plot)
ss <- as.data.frame(t(x[global_top, , drop = FALSE]), check.names = FALSE)
ss$Other <- pmax(0, 1 - rowSums(ss))
ss$Sample <- rownames(ss)
ss <- dplyr::left_join(ss, sample_info, by = "Sample")
ss_long <- tidyr::pivot_longer(ss, cols = tidyselect::all_of(c(global_top, "Other")),
                               names_to = "FeatureID", values_to = "RelativeAbundance")
ss_long$Species <- ifelse(ss_long$FeatureID == "Other", "Other",
                           sp_labels[ss_long$FeatureID])
ss_long$RelativePercent <- 100 * ss_long$RelativeAbundance
sample_order <- sample_info$Sample[order(sample_info$Group, sample_info$Sample)]
ss_long$Sample <- factor(ss_long$Sample, levels = sample_order)
ss_long$Species <- factor(ss_long$Species,
                          levels = c(sp_labels[global_top], "Other"))
write_tsv(ss_long, file.path(dirs["ABUND"], "top_global_species_stacked_per_sample.tsv"))
p_samples <- ggplot2::ggplot(ss_long,
                             ggplot2::aes(Sample, RelativePercent, fill = Species)) +
  ggplot2::geom_col(width = 0.93) +
  ggplot2::facet_grid(~ Group, scales = "free_x", space = "free_x") +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::labs(x = "Sample", y = "Relative abundance (%)",
                title = paste("Top", top_n_plot, "species across samples")) +
  base_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                      vjust = 0.5, size = 6),
                 legend.key.height = grid::unit(0.4, "cm"))
save_plot(p_samples, "sample_composition_stacked", dirs["ABUND"],
          width = max(12, min(28, 6 + ncol(x) * .22)), height = 7)

# 10. MAASLIN2 DIFFERENTIAL ABUNDANCE ----------------------------------------
# Use samples x features (not features x samples).
# Group factor reference is controlled explicitly and documented in QC files.
# Complete-case selection is model-specific (including chosen covariates).
model_terms <- unique(c(group_column, covariates, random_effects))
model_meta <- metadata[, model_terms, drop = FALSE]
model_meta <- as.data.frame(lapply(model_meta, function(z) {
  if (is.list(z)) stop("MaAsLin2 metadata contains list-columns; clean Excel import.")
  if (is.character(z)) trimws(z) else z
}), check.names = FALSE, row.names = rownames(metadata))
model_meta[[group_column]] <- factor(as.character(model_meta[[group_column]]),
                                     levels = group_levels)
for (z in intersect(numeric_covariates, covariates)) {
  model_meta[[z]] <- suppressWarnings(as.numeric(as.character(model_meta[[z]])))
}
for (z in setdiff(covariates, numeric_covariates)) {
  if (is.character(model_meta[[z]]) || is.logical(model_meta[[z]])) {
    model_meta[[z]] <- factor(model_meta[[z]])
  }
}
for (z in random_effects) model_meta[[z]] <- factor(model_meta[[z]])

valid <- stats::complete.cases(model_meta)
if (any(!valid)) {
  message("MaAsLin2: removing ", sum(!valid), " samples with missing model terms")
}
model_meta <- droplevels(model_meta[valid, , drop = FALSE])
if (nrow(model_meta) < 4 || nlevels(model_meta[[group_column]]) < 2) {
  stop("Too few complete MaAsLin2 samples, or <2 groups after filtering.")
}
for (z in setdiff(model_terms, group_column)) {
  if (length(unique(model_meta[[z]])) < 2L) {
    stop("Covariate/random effect has fewer than two observed values: ", z)
  }
}
maaslin_abundance <- as.data.frame(t(x[, rownames(model_meta), drop = FALSE]),
                                   check.names = FALSE)
stopifnot(identical(rownames(maaslin_abundance), rownames(model_meta)))
write_tsv(data.frame(Sample = rownames(model_meta), model_meta, check.names = FALSE),
          file.path(dirs["MAASLIN"], "model_metadata.tsv"))
write_tsv(data.frame(Group = names(table(model_meta[[group_column]])),
                     N = as.integer(table(model_meta[[group_column]]))),
          file.path(dirs["MAASLIN"], "model_samples_per_group.tsv"))

# MaAsLin2 normalization="TSS" is retained for consistency with your prior
# pipelines; here the rows already sum to 1, so TSS is essentially a no-op.
# LOG uses MaAsLin2's zero handling. Coefficients refer to *relative* abundance.
fit <- Maaslin2::Maaslin2(
  input_data = maaslin_abundance,
  input_metadata = model_meta,
  output = file.path(dirs["MAASLIN"], "model"),
  fixed_effects = unique(c(group_column, covariates)),
  random_effects = if (length(random_effects)) random_effects else NULL,
  reference = if (!is.null(reference_group)) {
    paste(group_column, reference_group, sep = ",")
  } else NULL,
  normalization = "TSS", transform = "LOG", analysis_method = "LM",
  min_abundance = 0, min_prevalence = min_prevalence_maaslin,
  min_variance = 0, max_significance = fdr_cutoff,
  correction = "BH", standardize = FALSE,
  plot_heatmap = TRUE, plot_scatter = FALSE, cores = 1
)

# Make an annotated table for group-associated species.
all_file <- file.path(dirs["MAASLIN"], "model", "all_results.tsv")
if (file.exists(all_file)) {
  res <- data.table::fread(all_file, data.table = FALSE)
  if (all(c("feature", "metadata", "qval") %in% names(res))) {
    group_res <- res[res$metadata == group_column, , drop = FALSE]
    group_res$Species <- sp_labels[group_res$feature]
    group_res <- group_res[order(group_res$qval), , drop = FALSE]
    write_tsv(group_res, file.path(dirs["MAASLIN"],
                                     "group_associations_all.tsv"))
    sig <- group_res[!is.na(group_res$qval) & group_res$qval < fdr_cutoff,
                     , drop = FALSE]
    write_tsv(sig, file.path(dirs["MAASLIN"],
                                 "group_associations_significant.tsv"))
    message("MaAsLin2 significant group/species associations (q < ",
            fdr_cutoff, "): ", nrow(sig))
  }
}

# 11. FINISH -------------------------------------------------------------------
writeLines(c(
  "MetaPhlAn -> phyloseq single-study pipeline",
  paste("Data:", abundance_file), paste("Metadata:", metadata_file),
  paste("Grouping:", group_column), paste("Reference:", reference_group),
  paste("Covariates:", paste(covariates, collapse = ", ")),
  paste("Random effects (MaAsLin2):", paste(random_effects, collapse = ", ")),
  paste("Pairing (PERMANOVA blocks, paired alpha tests):",
        if (is.null(pairing_column)) "none" else pairing_column),
  paste("Retained species:", nrow(x)), paste("Retained samples:", ncol(x)),
  "Phyloseq object: 01_PHYLOSEQ/phyloseq_species.rds",
  "Metadata/sample audits: 00_QC/",
  "Alpha diversity: 02_ALPHA/",
  "Beta diversity: 03_BETA/",
  "Group and sample abundances: 04_ABUNDANCE/",
  "MaAsLin2: 05_MAASLIN2/",
  "Caveat: MetaPhlAn profiles are compositional estimates, not sequencing read counts.",
  "Taxa and UNCLASSIFIED not retained are excluded before sample-wise renormalization."
), con = file.path(rdir, "README_RESULTS.txt"))
capture.output(utils::sessionInfo(), file = file.path(dirs["QC"], "R_session_info.txt"))
message("Completed analysis. Outputs in: ", normalizePath(rdir))
