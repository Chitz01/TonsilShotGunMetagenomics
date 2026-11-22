## 07_rox_tree_summary.R
## Summarise the IQ-TREE tree from 06_rox_tree.sh:
##  (1) plot the tree, coloured by sequence type (Rox reference, CCP, MauG, tonsil/oral MAG)
##  (2) for every MAG protein: nearest characterised reference (Rox, CCP or MauG) by
##      patristic (along-the-tree) distance - this does not depend on where the root is
##  (3) clade test: which sequences fall in the smallest clade containing all Rox references,
##      and how well that clade is supported (SH-aLRT / UFBoot)
## Run in RStudio from the folder that contains the copied tree/ folder:
##     source("07_rox_tree_summary.R")
## Needs the 'ape' package: install.packages("ape")
## Code drafted with help from Claude (Anthropic); declare per course rules.

library(ape)

tree_dir <- "tree"   # local copy of rubber_search/tree from COSMOS
tr  <- read.tree(file.path(tree_dir, "rox_tree.treefile"))
ann <- read.delim(file.path(tree_dir, "tree_annotation.tsv"), stringsAsFactors = FALSE)
ann <- ann[match(tr$tip.label, ann$label), ]
stopifnot(!any(is.na(ann$label)))

# How many MAG proteins each tip stands for (identical sequences were collapsed)
n_copies <- setNames(rep(1, nrow(ann)), ann$label)
dup_file <- file.path(tree_dir, "duplicates.txt")
if (file.exists(dup_file) && file.size(dup_file) > 0) {
  for (l in readLines(dup_file)) {
    parts <- strsplit(l, "\t")[[1]]
    ids <- gsub("\\|", "__", trimws(strsplit(parts[2], ",")[[1]]))
    hit <- ids[ids %in% names(n_copies)]
    if (length(hit)) n_copies[hit[1]] <- as.numeric(parts[1])
  }
}
ann$n_MAG_copies <- n_copies[ann$label]

# Tip colours
ann$group <- ifelse(ann$type == "MAG", paste0("MAG_", sub("hmp_", "", ann$dataset)), ann$type)
cols <- c(Rox_reference = "red3", CCP_SwissProt = "royalblue", MauG_SwissProt = "purple",
          MAG_tonsil = "darkorange", MAG_oral = "forestgreen", "MAG_-" = "grey40")
tipcol <- cols[ann$group]; tipcol[is.na(tipcol)] <- "grey40"

# Root on MauG if they form one clade (MauG is a separate, well-known subfamily); else midpoint-free unrooted plot
mau <- ann$label[ann$type == "MauG_SwissProt"]
rooted <- FALSE
if (length(mau) >= 1) {
  if (length(mau) == 1 || is.monophyletic(tr, mau)) {
    tr <- root(tr, outgroup = mau, resolve.root = TRUE); rooted <- TRUE
    ann <- ann[match(tr$tip.label, ann$label), ]
    tipcol <- cols[ann$group]; tipcol[is.na(tipcol)] <- "grey40"
  }
}
cat("Tree rooted on MauG:", rooted, "\n")

# Short display labels: genus + locus for MAG proteins
short <- ifelse(ann$type == "MAG",
                paste0(sub("^g__", "", ann$genus), " ", sub("^.*__", "", ann$label),
                       ifelse(ann$n_MAG_copies > 1, paste0(" (x", ann$n_MAG_copies, ")"), "")),
                ann$label)
tr_plot <- tr; tr_plot$tip.label <- short

pdf(file.path(tree_dir, "rox_tree.pdf"), width = 9, height = max(8, 0.13 * Ntip(tr) + 2))
plot(tr_plot, type = if (rooted) "phylogram" else "unrooted", tip.color = tipcol,
     cex = 0.45, label.offset = 0.01, no.margin = TRUE)
if (rooted) nodelabels(tr$node.label, frame = "none", cex = 0.3, adj = c(1.1, -0.3), col = "grey30")
add.scale.bar(cex = 0.6)
legend("topright", legend = names(cols)[1:5], text.col = cols[1:5], bty = "n", cex = 0.7,
       title = "SH-aLRT/UFBoot at nodes")
dev.off()

# (2) Nearest characterised reference for each MAG protein (patristic distance)
d <- cophenetic.phylo(tr)
refs <- ann$label[ann$type != "MAG"]
mags <- ann$label[ann$type == "MAG"]
near <- do.call(rbind, lapply(mags, function(m) {
  dr <- d[m, refs]
  rox <- dr[ann$type[match(refs, ann$label)] == "Rox_reference"]
  data.frame(label = m, nearest_reference = names(which.min(dr)),
             nearest_type = ann$type[match(names(which.min(dr)), ann$label)],
             dist_nearest = round(min(dr), 3),
             dist_nearest_Rox = round(min(rox), 3))
}))
near <- merge(near, ann[, c("label", "dataset", "genus", "product", "length", "n_MAG_copies")], by = "label")
near <- near[order(near$dist_nearest_Rox), ]
write.table(near, file.path(tree_dir, "nearest_reference_per_MAG_protein.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nNearest characterised reference type for MAG proteins:\n")
print(table(near$nearest_type))

# (3) Clade test (needs a rooted tree)
if (rooted) {
  rox <- ann$label[ann$type == "Rox_reference"]
  node <- getMRCA(tr, rox)
  members <- extract.clade(tr, node)$tip.label
  support <- tr$node.label[node - Ntip(tr)]
  cat("\nSmallest clade containing all Rox references: ", length(members), " tips; support (SH-aLRT/UFBoot): ",
      support, "\n", sep = "")
  print(table(ann$type[match(members, ann$label)]))
  writeLines(members, file.path(tree_dir, "rox_clade_members.txt"))
} else {
  cat("\nMauG relatives not monophyletic/absent: clade test skipped; use the nearest-reference table.\n")
}
cat("\nOutputs in", normalizePath(tree_dir), ": rox_tree.pdf, nearest_reference_per_MAG_protein.tsv",
    if (rooted) ", rox_clade_members.txt", "\n")
