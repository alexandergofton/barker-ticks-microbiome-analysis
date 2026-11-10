#!/usr/bin/env Rscript
# Step 19b -- root each marker-gene tree and draw it.
#
# Reads every tree written by step 19 and writes, beside it:
#
#   <target>/<target>__<marker>.rooted.nwk   the rooted tree (the one plotted)
#   <target>/<target>__<marker>.png          the figure
#
# ROOTING. On the outgroup when step 18 found one -- tips tagged OUTGROUP --
# and by midpoint otherwise. Which of the two was used is printed on the figure
# itself, because a midpoint root is a display convention, not a phylogenetic
# claim, and a tree read without that caveat will be over-interpreted. The
# unrooted FastTree output is left untouched in <target>__<marker>.nwk.
#
# READING THE FIGURE
#   red      sequence recovered from one of these ticks (tip = sample | tick sp.)
#   black    reference downloaded from GenBank (tip = organism + accession)
#   blue     outgroup
#   numbers  FastTree SH-like local support, drawn only at or above
#            trees.figures.min_support_shown. These are NOT bootstrap values.
#
# Run with: Rscript code/analysis/19b_tree_figures.R [target ...]
# (the tick_micro env; see env/tick_micro.yml for r-ape / r-phangorn)

suppressMessages({
  library(ape)
  library(phangorn)
  library(tidyverse)
})

.args <- commandArgs(trailingOnly = FALSE)
.this <- .args |> str_subset("^--file=") |> str_remove("^--file=")
source(file.path(dirname(normalizePath(.this)), "common.R"))

argv       <- commandArgs(trailingOnly = TRUE)
params     <- load_params()
trees_dir  <- file.path(ROOT, params$trees$out_dir)
fig        <- params$trees$figures
min_sup    <- as.numeric(fig$min_support_shown)
root_mode  <- params$trees$tree$root

if (!dir.exists(trees_dir)) stop("missing ", trees_dir, " -- run step 19 first")

nwk <- list.files(trees_dir, pattern = "\\.nwk$", recursive = TRUE, full.names = TRUE)
nwk <- str_subset(nwk, "\\.rooted\\.nwk$", negate = TRUE)
if (length(argv)) {
  nwk <- keep(nwk, \(p) basename(dirname(p)) %in% argv)
}
if (!length(nwk)) stop("no .nwk files under ", trees_dir)

# ---- curated tick species ---------------------------------------------------
# The species inside a sample tip label was frozen when step 12 ran, often from
# the id-string guess. For the FIGURE, the curated metadata wins when it exists,
# so the trees show current identifications without re-running steps 12-19.
# The .nwk files keep their original labels untouched.
curated <- character()
meta_f <- file.path(ROOT, "metadata", "sample_metadata.tsv")
if (file.exists(meta_f)) {
  curated <- read_tsv(meta_f, col_types = cols(.default = col_character()),
                      na = c("", "NA"), trim_ws = FALSE) |>
    filter(!is.na(tick_species)) |>
    pull(tick_species, name = sample_id)
}

# ---- tip labels -------------------------------------------------------------
# Both naming schemes put the target at field 3 and the marker at field 4:
#   REF|<acc>|<organism>|<target>|<marker>|<len>bp
#   <sample>|<tick species>|<target>|<marker>|<contig>|<len>bp[|copyN]
parse_tips <- function(labels) {
  # field i of every label, "" where a label has fewer than i fields
  fld <- function(i) replace_na(str_split_i(labels, fixed("|"), i), "")
  kind <- fld(1)
  is_ref <- kind %in% c("REF", "OUTGROUP")
  organism <- str_replace_all(fld(3), "_", " ")
  copy <- if_else(fld(7) != "", paste0(" ", fld(7)), "")
  tick <- if_else(fld(1) %in% names(curated), unname(curated[fld(1)]),
                  str_replace_all(fld(2), "_", " "))
  pretty <- if_else(
    is_ref,
    paste0(organism, "  ", fld(2)),
    paste0(fld(1), "  [", tick, "]", copy))
  tibble(is_ref = is_ref, is_out = kind == "OUTGROUP", pretty = pretty)
}

drawn <- 0L
for (path in nwk) {
  tr <- tryCatch(read.tree(path), error = function(e) NULL)
  if (is.null(tr) || is.null(tr$tip.label) || length(tr$tip.label) < 3) {
    message("skipping ", basename(path), " (unreadable or <3 tips)")
    next
  }
  info <- parse_tips(tr$tip.label)
  outs <- tr$tip.label[info$is_out]

  # edgelabel = TRUE is the correct handling for FastTree support values: they
  # describe the SPLIT below a node, so on re-rooting they must travel with the
  # edge, not stay attached to a node that now means something else.
  rooted <- tr
  how <- "unrooted (as built)"
  if (root_mode != "none") {
    if (length(outs) && root_mode == "outgroup_else_midpoint") {
      rooted <- tryCatch(
        root(tr, outgroup = outs, resolve.root = TRUE, edgelabel = TRUE),
        error = function(e) { message("  outgroup rooting failed: ",
                                      conditionMessage(e)); NULL })
      how <- paste0("rooted on outgroup (", length(outs), " tip",
                    if (length(outs) > 1) "s" else "", ")")
      if (is.null(rooted)) { rooted <- tr; how <- "unrooted (outgroup rooting failed)" }
    } else {
      rooted <- tryCatch(midpoint(tr, node.labels = "support"),
                         error = function(e) { message("  midpoint rooting failed: ",
                                                       conditionMessage(e)); NULL })
      how <- "midpoint-rooted (display only -- no outgroup available)"
      if (is.null(rooted)) { rooted <- tr; how <- "unrooted (midpoint rooting failed)" }
    }
  }

  info <- parse_tips(rooted$tip.label)
  maxch <- 70
  labs <- if_else(str_length(info$pretty) > maxch,
                  paste0(str_sub(info$pretty, 1, maxch - 1), "…"), info$pretty)
  plot_tree <- rooted
  plot_tree$tip.label <- labs
  # CSIRO: outgroup Ocean Blue, references black, samples Orange (blue/orange
  # stay distinct under every common colour-vision deficiency; font also differs)
  col <- case_when(info$is_out ~ CSIRO[["OceanBlue"]],
                   info$is_ref ~ "black",
                   .default = CSIRO[["Orange"]])
  fnt <- if_else(info$is_ref, 3, 2)   # references italic, samples bold

  sup <- suppressWarnings(as.numeric(rooted$node.label))
  sup_txt <- if_else(is.na(sup) | sup < min_sup, "", formatC(sup, format = "f", digits = 2))

  n_tip  <- length(labs)
  n_samp <- sum(!info$is_ref)
  stem   <- str_remove(path, "\\.nwk$")
  locus  <- str_replace_all(basename(stem), "__", "  ")
  height <- max(as.numeric(fig$min_height_in),
                n_tip * as.numeric(fig$height_per_tip_in) + 2)

  write.tree(rooted, file = paste0(stem, ".rooted.nwk"))

  # ape sizes the plot from the tree and then hopes the tip labels fit. With
  # organism names and accessions they routinely do not, so reserve the space
  # explicitly: a character at cex 0.55 is about 0.055 in wide, so the longest
  # label claims that fraction of the panel and the tree gets the rest.
  cex_tip <- 0.55
  width_in <- as.numeric(fig$width_in)
  lab_frac <- min(0.62, max(str_length(labs)) * 0.055 * cex_tip / 0.55 / width_in)
  depth <- max(node.depth.edgelength(plot_tree))
  if (!is.finite(depth) || depth <= 0) depth <- 1

  png_svg(paste0(stem, ".png"), width = width_in, height = height, res = 200, draw = function() {
  op <- par(mar = c(3.5, 1, 3.5, 0.5), xpd = NA)
  plot.phylo(plot_tree, cex = cex_tip, tip.color = col, font = fnt,
             label.offset = depth * 0.01, x.lim = c(0, depth / (1 - lab_frac)),
             no.margin = FALSE, edge.width = 0.8)
  if (any(nzchar(sup_txt)))
    nodelabels(sup_txt, frame = "none", adj = c(1.15, -0.35), cex = 0.45, col = "grey35")
  add.scale.bar(lwd = 1.2, cex = 0.6)
  title(main = locus, line = 2.1, cex.main = 1.05, font.main = 2)
  mtext(sprintf("%d tips: %d from these ticks (orange, bold), %d GenBank references (black; outgroup blue) | %s",
                n_tip, n_samp, n_tip - n_samp, how),
        side = 3, line = 0.8, cex = 0.62, col = "grey25")
  mtext(paste("FastTree GTR+CAT, SH-like local support shown at >=", min_sup,
              "-- a screening tree, not a species identification;",
              "many sample sequences are gene fragments"),
        side = 1, line = 2.2, cex = 0.55, col = "grey40")
  par(op)
  })
  drawn <- drawn + 1L
  message(sprintf("  %-42s %3d tips  %s", basename(stem), n_tip, how))
}

message(sprintf("\n%d figures written under %s/", drawn,
                str_remove(trees_dir, fixed(paste0(ROOT, "/")))))
