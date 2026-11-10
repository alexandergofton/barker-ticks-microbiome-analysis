# Shared helpers for the R side of the tick microbiome analysis.
# source() this at the top of every Rscript in code/analysis/.

# ---- project root ----------------------------------------------------------
# `sys.frame(1)$ofile` is NULL under Rscript, so derive the root from the
# --file= argument instead, falling back to TICK_ROOT then getwd().
.script_path <- function() {
  f <- commandArgs(trailingOnly = FALSE) |>
    stringr::str_subset("^--file=") |>
    stringr::str_remove("^--file=")
  if (length(f)) normalizePath(f[1], mustWork = FALSE) else NA_character_
}

project_root <- function() {
  env <- Sys.getenv("TICK_ROOT", "")
  if (nzchar(env)) {
    return(normalizePath(env, mustWork = TRUE))
  }
  sp <- .script_path()
  if (!is.na(sp)) {
    r <- normalizePath(file.path(dirname(sp), "..", ".."), mustWork = FALSE)
    if (dir.exists(file.path(r, "results"))) return(r)
  }
  d <- normalizePath(getwd(), mustWork = TRUE)
  while (!dir.exists(file.path(d, "results")) && dirname(d) != d) {
    d <- dirname(d)
  }
  if (!dir.exists(file.path(d, "results"))) {
    stop("cannot locate project root; set TICK_ROOT")
  }
  d
}

ROOT <- project_root()

load_params <- function(root = ROOT) {
  yaml::read_yaml(file.path(root, "config", "analysis_params.yaml"))
}

# ---- I/O -------------------------------------------------------------------
# Tables are read with data.table::fread() and written with fwrite(), then
# handled as tibbles. readr would read every whole-number column as double
# (not integer) and format numbers differently on write, changing the outputs.
read_matrix <- function(path) {
  # taxa x samples; first column is the taxon label
  df <- tibble::as_tibble(data.table::fread(path, sep = "\t", header = TRUE))
  m <- as.matrix(dplyr::select(df, -1))
  rownames(m) <- dplyr::pull(df, 1)
  m
}

fig_path <- function(name, root = ROOT) {
  d <- file.path(root, "results", "figures")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.path(d, name)
}

# Every figure is written twice: the .png it is named as, and an .svg beside it.
# SVGs use svglite (text stays editable <text>). svglite is not in the pinned
# conda envs; it lives in a separate R library per env, built from CRAN source
# (see env/README_svglite.md), appended LAST to .libPaths() so it can never
# shadow an env package. If it is missing, SVGs fall back to Cairo
# grDevices::svg() (text as outlined glyphs).
local({
  env <- if (grepl("tick_report", R.home(), fixed = TRUE)) "r_lib_svglite_report" else "r_lib_svglite"
  lib <- file.path(ROOT, "env", env)
  if (dir.exists(lib)) .libPaths(c(.libPaths(), lib))
})
HAVE_SVGLITE <- requireNamespace("svglite", quietly = TRUE)
svg_device <- if (HAVE_SVGLITE) function(filename, width, height, ...)
  svglite::svglite(filename, width = width, height = height, ...) else grDevices::svg
# save_plot() replaces ggsave() for ggplot figures; png_svg() runs a
# base-graphics drawing function once per device (sizes in inches).
svg_path <- function(png_file) stringr::str_replace(png_file, "\\.png$", ".svg")
save_plot <- function(filename, plot = ggplot2::last_plot(), ...) {
  ggplot2::ggsave(filename, plot, ...)
  ggplot2::ggsave(svg_path(filename), plot, device = svg_device, ...)
  invisible(filename)
}
png_svg <- function(png_file, width, height, res, draw) {
  grDevices::png(png_file, width = width, height = height, units = "in", res = res)
  draw(); grDevices::dev.off()
  svg_device(svg_path(png_file), width = width, height = height)
  draw(); grDevices::dev.off()
  invisible(png_file)
}

tab_path <- function(name, root = ROOT) {
  d <- file.path(root, "results", "tables")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.path(d, name)
}

# ---- framing ---------------------------------------------------------------
# Reused verbatim on every figure. This study has no extraction blanks or
# no-template controls, so composition is descriptive and contaminant taxa are
# annotated rather than removed. Saying so once per figure is the point.
CAVEAT <- paste(
  "Descriptive only: no extraction blanks or no-template controls were",
  "sequenced, so absolute microbial composition cannot be validated.",
  "Reagent-contaminant taxa are annotated, not removed; unidentifiable taxa",
  "and excluded contaminant genera (config taxon_removal) are removed."
)

say <- function(...) cat(sprintf(...), "\n", sep = "")

# ---- colours and theme: CSIRO, colour-blind checked -------------------------
# The CSIRO brand colours live in CSIRO_colours.R at the project root (one
# object per colour). They are read into a named vector here, never as loose
# globals, so a script's own `Orange` or `Black` cannot clash with them.
CSIRO <- local({
  e <- new.env()
  sys.source(file.path(ROOT, "CSIRO_colours.R"), envir = e)
  unlist(mget(ls(e), envir = e))
})

# Categorical order. Nine CSIRO hues chosen and ordered by computation (the
# dataviz skill's validator, Machado-2009 protan/deutan simulation, OKLab
# Delta E x100): every PAIR of the nine stays >= 8.7 apart under both CVD types
# and >= 15.8 under normal vision, and each prefix is the best-separated set of
# that size (the first 3 >= 25, first 5 >= 18.6). Midday/Midnight Blue, Gold
# and Light Mint sit outside the validator's lightness band and Gold/Light
# Mint have < 3:1 contrast on white -- CSIRO brand values cannot be re-stepped,
# so filled marks always carry a dark outline (MARK_OUTLINE) and identity is
# never carried by colour alone (shape, facet or direct label as well).
CSIRO_CAT <- unname(CSIRO[c(
  "Blueberry",
  "Mint",
  "Gold",
  "MiddayBlue",
  "MidnightBlue",
  "Orange",
  "LightMint",
  "Fuschia",
  "Plum"
)])
CSIRO_OTHER <- unname(CSIRO["Mist"]) # "other" / "unidentified" / not detected
CSIRO_MUTED <- unname(CSIRO["Steel"]) # reference lines, secondary ink
MARK_OUTLINE <- unname(CSIRO["MidnightBlue"])
# Sequential (magnitude): one hue, light -> dark, white to Midnight Blue
# through Midday and Ocean Blue. Lightness is monotone, so it reads under any
# colour-vision type.
CSIRO_SEQ <- c(
  "#FFFFFF",
  "#CCEEF5",
  unname(CSIRO[c("MiddayBlue", "OceanBlue", "MidnightBlue")])
)
csiro_seq <- function(n) grDevices::colorRampPalette(CSIRO_SEQ)(n)
# First n categorical colours; past nine, fold into "other" rather than invent hues.
csiro_cat <- function(n) {
  if (n > length(CSIRO_CAT)) {
    stop(
      "more than ",
      length(CSIRO_CAT),
      " categories: fold the rest into 'other'"
    )
  }
  CSIRO_CAT[seq_len(n)]
}

# Every ggplot: theme_classic with light-grey major gridlines both ways.
theme_tick <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size)
}
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::theme_set(theme_tick())
}

# ---- tick taxonomy: fixed colours and shapes ----------------------------------
# One colour per tick genus, fixed to the ENTITY (never to rank or count), shared
# by every figure and both maps so a genus reads the same everywhere. Largest
# genera take the best-separated slots. All nine genera together still pass
# all-pairs CVD separation (above), and every figure that uses this palette
# also carries a second encoding (facets, labels, or the host_type shape).
TICK_GENUS_PAL <- c(
  Amblyomma = CSIRO_CAT[2], # Mint
  Ixodes = CSIRO_CAT[1], # Blueberry
  Haemaphysalis = CSIRO_CAT[6], # Orange
  Bothriocroton = CSIRO_CAT[3], # Gold
  Argas = CSIRO_CAT[8], # Fuschia
  Ornithodoros = CSIRO_CAT[9], # Plum
  Alectorobius = CSIRO_CAT[4], # Midday Blue
  Cryptocroton = CSIRO_CAT[7], # Light Mint
  Rhipicephalus = CSIRO_CAT[5], # Midnight Blue
  unidentified = CSIRO_OTHER
)

HOST_TYPE_SHAPE <- c(
  Mammal = 21,
  Reptile = 24,
  Bird = 22,
  Questing = 23,
  unknown = 25
)

# Ordered factor of tick genus with NA -> "unidentified", in palette order.
tick_genus_factor <- function(x) {
  x <- as.character(x)
  x <- dplyr::if_else(is.na(x) | x == "", "unidentified", x)
  factor(
    x,
    levels = intersect(names(TICK_GENUS_PAL), c(unique(x), "unidentified"))
  )
}
host_type_factor <- function(x) {
  x <- as.character(x)
  x <- dplyr::if_else(is.na(x) | x == "", "unknown", x)
  factor(x, levels = names(HOST_TYPE_SHAPE))
}

# ---- mitogenome tree helpers (steps 26, 26b) --------------------------------
# Published reference mitogenomes (step 20b) appear in the tree as tips whose
# id starts "REF_". One row per reference that made it into the analysis.
mito_ref_table <- function(params = load_params()) {
  p <- file.path(
    ROOT,
    params$mito$out_dir,
    "ref_genomes",
    "ref_genome_summary.tsv"
  )
  if (!file.exists(p)) {
    return(tibble::tibble(
      ref_id = character(),
      family = character(),
      genus = character(),
      species = character(),
      accession = character()
    ))
  }
  data.table::fread(p, na.strings = c("", "NA")) |>
    tibble::as_tibble() |>
    dplyr::filter(status == "ok") |>
    dplyr::select(ref_id, accession, family, genus, species)
}

# Root the mitogenome tree on the branch that separates the reference tips of
# `family` (default Argasidae) from all OTHER reference tips. Only references
# define that branch; samples may fall on either side of it (soft-tick samples
# belong inside the outgroup clade), so a mislabelled sample cannot move the
# root. Returns the tree with attribute "rooting" saying what was done, and
# attribute "outgroup_samples" listing the samples that fell in the outgroup
# clade; falls back to midpoint rooting when no single branch separates the two
# reference sets.
root_mito_tree <- function(tree, refs, family) {
  sid <- stringr::str_split_i(tree$tip.label, stringr::fixed("|"), 1)
  is_rf <- sid %in% refs$ref_id
  og <- tree$tip.label[is_rf & sid %in% refs$ref_id[refs$family %in% family]]
  other <- tree$tip.label[
    is_rf & !(sid %in% refs$ref_id[refs$family %in% family])
  ]
  out <- NULL
  if (length(og) >= 1 && length(other) >= 1) {
    out <- tryCatch(
      {
        t0 <- ape::root(tree, outgroup = other[1], resolve.root = TRUE)
        node <- if (length(og) == 1) {
          which(t0$tip.label == og)
        } else {
          ape::getMRCA(t0, og)
        }
        clade <- if (node <= ape::Ntip(t0)) {
          t0$tip.label[node]
        } else {
          ape::extract.clade(t0, node)$tip.label
        }
        if (any(other %in% clade)) {
          stop(
            "no branch separates the ",
            paste(family, collapse = "/"),
            " references from the other references"
          )
        }
        t <- ape::root(t0, node = node, resolve.root = TRUE)
        og_samples <- setdiff(clade, og)
        attr(t, "rooting") <- sprintf(
          "rooted on the %s clade (%d references + %d samples)",
          paste(family, collapse = "/"),
          length(og),
          length(og_samples)
        )
        attr(t, "outgroup_samples") <- stringr::str_remove(og_samples, "\\|.*")
        t
      },
      error = function(e) {
        message(
          "outgroup rooting failed (",
          conditionMessage(e),
          "); using midpoint"
        )
        NULL
      }
    )
  }
  if (is.null(out)) {
    out <- phangorn::midpoint(tree)
    attr(out, "rooting") <- "midpoint-rooted (no usable outgroup)"
  }
  out
}
