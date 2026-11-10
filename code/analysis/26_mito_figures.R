#!/usr/bin/env Rscript
# Step 26 -- final figure and distance tables for the mitogenome phylogeny.
#
# Reads the partitioned IQ-TREE result (step 25) and the concatenated
# supermatrix (step 24), and writes:
#
#   results/mito/figures/mito_tree.png / .pdf   the tree, rooted on the
#                                               Argasidae reference mitogenomes
#                                               (midpoint if that fails); sample
#                                               tips coloured by tick genus,
#                                               published references in black
#                                               bold; SH-aLRT / ultrafast-
#                                               bootstrap support at nodes
#   results/mito/tables/mito_distance_pdistance.tsv
#       pairwise uncorrected p-distance (%) from the concatenated alignment,
#       PAIRWISE DELETION -- compared only over columns where BOTH sequences
#       have real data, which matters here specifically because a sample
#       missing a gene contributes a whole all-gap block (step 24) that must
#       not be counted as sequence difference against a sample that has it.
#   results/mito/tables/mito_distance_patristic.tsv
#       tree-based (patristic) distance between every pair of tips.
#
# Run with: Rscript code/analysis/26_mito_figures.R
# (the tick_micro env; see env/tick_micro.yml for r-ape / r-phangorn)

suppressMessages({
  library(ape)
  library(phangorn)
  library(tidyverse)
})

.args <- commandArgs(trailingOnly = FALSE)
.this <- .args |> str_subset("^--file=") |> str_remove("^--file=")
source(file.path(dirname(normalizePath(.this)), "common.R"))

params   <- load_params()
mito_dir <- file.path(ROOT, params$mito$out_dir)
fig_dir  <- file.path(mito_dir, "figures")
tab_dir  <- file.path(mito_dir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

CAVEAT <- paste(
  "One mitogenome candidate per sample, chosen by homology (order-wide tblastn/",
  "blastn against 13 PCGs + 2 rRNAs); annotation is homology-based, not MITOS2",
  "(its Zenodo reference data was unreachable from this compute node -- see",
  "README_mito.md). Gene recovery varies by sample (results/mito/mito_selection.tsv);",
  "samples below the configured minimum are excluded (results/mito/concat/excluded_samples.tsv).",
  "Published reference mitogenomes (metadata/mito_ref_genomes.tsv, step 20b) in black bold,",
  "annotated by the same homology method as the samples.")

# ---- tree -------------------------------------------------------------------
tree_file <- file.path(mito_dir, "tree", "mito.treefile")
if (!file.exists(tree_file)) stop("missing ", tree_file, " -- run step 25 first")
tree <- read.tree(tree_file)

refs <- mito_ref_table(params)
tree_rooted <- root_mito_tree(tree, refs, params$mito$tree$outgroup_family)
rooting <- attr(tree_rooted, "rooting")
message("tree ", rooting)
CAVEAT <- paste(CAVEAT, "Tree", rooting, "for display.")

sample_id <- str_split_i(tree_rooted$tip.label, fixed("|"), 1)
species   <- str_split_i(tree_rooted$tip.label, fixed("|"), 2) |>
  replace_na("unknown_species") |>
  str_replace_all("_", " ")
genus     <- if_else(species == "unknown species", "unknown",
                     str_split_i(species, fixed(" "), 1))
is_ref    <- str_starts(sample_id, fixed("REF_"))

# the project-wide CSIRO tick-genus palette (common.R), so a genus is the same
# colour here as in every other figure; anything outside it is "unidentified"
genus <- if_else(genus %in% names(TICK_GENUS_PAL), genus, "unidentified")
genus_levels <- intersect(names(TICK_GENUS_PAL), unique(genus[!is_ref]))
pal <- TICK_GENUS_PAL[genus_levels]
tip_col <- if_else(is_ref, "black", unname(pal[genus]))   # marker fill; label text stays dark
tip_font <- if_else(is_ref, 2, 1)
# references: "Species ACCESSION"; samples keep their full name
tree_plot <- tree_rooted
tree_plot$tip.label <- if_else(is_ref,
                               sprintf("%s  %s", species, str_remove(sample_id, "^REF_")),
                               tree_plot$tip.label)

# Node labels from IQ-TREE are "SH-aLRT/UFBoot"; keep both but drop the root's
# empty label and anything under a low-support floor to avoid a cluttered plot.
node_lab <- tree_rooted$node.label
node_lab <- if_else(is.na(node_lab) | node_lab == "", "", node_lab)

n_tips <- length(tree_rooted$tip.label)
plot_height <- max(6, n_tips * 0.16)

for (dev in list(list(fn = file.path(fig_dir, "mito_tree.png"), open = png,
                      args = list(width = 2000, height = plot_height * 150, res = 150)),
                 list(fn = file.path(fig_dir, "mito_tree.svg"), open = svg_device,
                      args = list(width = 2000 / 150, height = plot_height)),
                 list(fn = file.path(fig_dir, "mito_tree.pdf"), open = pdf,
                      args = list(width = 11, height = plot_height)))) {
  do.call(dev$open, c(list(dev$fn), dev$args))
  par(mar = c(1, 1, 3, 12), xpd = NA)
  off <- 0.012 * max(ape::node.depth.edgelength(tree_plot))
  plot.phylo(tree_plot, tip.color = ifelse(is_ref, "black", "grey20"), cex = 0.55,
            label.offset = off, no.margin = FALSE, font = tip_font)
  tiplabels(pch = 22, bg = tip_col, col = MARK_OUTLINE, cex = 0.75, lwd = 0.4, offset = off / 2)
  nodelabels(node_lab, frame = "none", cex = 0.45, adj = c(1.1, -0.3), col = CSIRO_MUTED)
  title(main = sprintf("Tick mitogenome phylogeny -- %d samples + %d reference mitogenomes, %s partitions (IQ-TREE, ModelFinder + SH-aLRT/UFBoot)",
                       sum(!is_ref), sum(is_ref), nrow(read_tsv(file.path(mito_dir, "concat", "mito_partition_summary.tsv"),
                                                                 show_col_types = FALSE))),
        cex.main = 0.85)
  legend("bottomleft", legend = c(genus_levels, "published reference"),
        pt.bg = c(pal, "black"), col = MARK_OUTLINE, pch = 22, pt.cex = 1.1, bty = "n", cex = 0.6,
        title = "sample genus (as named at step 21)", inset = c(-0.02, 0.03))
  mtext(paste(strwrap(CAVEAT, width = 160), collapse = "\n"), side = 1, line = -1,
       cex = 0.45, col = "grey40", adj = 0)
  dev.off()
}
message("wrote ", file.path(fig_dir, "mito_tree.png"), ", .svg and .pdf")

# ---- distances ---------------------------------------------------------------
concat_file <- file.path(mito_dir, "concat", "mito_concat.fasta")
dna <- read.dna(concat_file, format = "fasta")

pdist <- as.matrix(dist.dna(dna, model = "raw", pairwise.deletion = TRUE)) * 100
write.table(round(pdist, 2), file.path(tab_dir, "mito_distance_pdistance.tsv"),
           sep = "\t", quote = FALSE, col.names = NA)
message("wrote ", file.path(tab_dir, "mito_distance_pdistance.tsv"),
       " (", nrow(pdist), " x ", ncol(pdist), ", pairwise-deletion p-distance, %)")

patristic <- cophenetic.phylo(tree_rooted)
patristic <- patristic[rownames(pdist), colnames(pdist)]
write.table(round(patristic, 5), file.path(tab_dir, "mito_distance_patristic.tsv"),
           sep = "\t", quote = FALSE, col.names = NA)
message("wrote ", file.path(tab_dir, "mito_distance_patristic.tsv"),
       " (tree-based patristic distance, substitutions/site)")

message("\n", CAVEAT)
