# Shared helpers for the R side of the tick microbiome analysis.
# source() this at the top of every Rscript in code/analysis/.

# ---- project root ----------------------------------------------------------
# `sys.frame(1)$ofile` is NULL under Rscript, so derive the root from the
# --file= argument instead, falling back to TICK_ROOT then getwd().
.script_path <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(f[1], mustWork = FALSE) else NA_character_
}

project_root <- function() {
  env <- Sys.getenv("TICK_ROOT", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = TRUE))
  sp <- .script_path()
  if (!is.na(sp)) {
    r <- normalizePath(file.path(dirname(sp), "..", ".."), mustWork = FALSE)
    if (dir.exists(file.path(r, "results"))) return(r)
  }
  d <- normalizePath(getwd(), mustWork = TRUE)
  while (!dir.exists(file.path(d, "results")) && dirname(d) != d) d <- dirname(d)
  if (!dir.exists(file.path(d, "results")))
    stop("cannot locate project root; set TICK_ROOT")
  d
}

ROOT <- project_root()

load_params <- function(root = ROOT) {
  yaml::read_yaml(file.path(root, "config", "analysis_params.yaml"))
}

# ---- I/O -------------------------------------------------------------------
read_matrix <- function(path) {
  # taxa x samples; first column is the taxon label
  dt <- data.table::fread(path, sep = "\t", header = TRUE)
  lab <- dt[[1]]
  m <- as.matrix(dt[, -1, with = FALSE])
  rownames(m) <- lab
  m
}

fig_path <- function(name, root = ROOT) {
  d <- file.path(root, "results", "figures")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.path(d, name)
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
  "Contaminant-flagged taxa are annotated, never removed.")

say <- function(...) cat(sprintf(...), "\n", sep = "")
