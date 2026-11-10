#!/usr/bin/env Rscript
# Step 7 -- assemble the TreeSummarizedExperiment.
#
#   conda activate tick_micro
#   Rscript code/analysis/07_build_tse.R
#
# ASSAYS
#   tpm_microbial   PRIMARY abundance. TPM renormalised within the microbial
#                   fraction, sums to 1e6 per sample. Host dominates these
#                   libraries, so without this renormalisation between-sample
#                   comparisons are driven by how much tick DNA each sample
#                   happened to yield rather than by microbial composition.
#   reads           EVIDENCE ONLY. Integer mapped-read support per taxon. Use
#                   to gate presence calls and to judge confidence. Do NOT use
#                   as an abundance and do NOT run count-based differential
#                   abundance on it -- these are assembly-derived contig read
#                   counts, not a designed count matrix.
#   n_contigs       EVIDENCE ONLY. How many contigs support the taxon.
#
# rowData = results/taxonomy/taxon_lineage.tsv  (lineage + scope + contaminant flag)
#
# Unidentifiable taxa and contaminant genera listed under taxon_removal in
# config/analysis_params.yaml are REMOVED here and TPM renormalised.
# colData = metadata/sample_metadata.tsv + per-sample QC from fraction_summary.tsv

# tidyverse is attached LAST so dplyr's verbs (filter, rename, desc, ...) are
# not masked by the Bioconductor packages.
suppressPackageStartupMessages({
  library(yaml)
  library(TreeSummarizedExperiment)
  library(mia)
  library(S4Vectors)
  library(tidyverse)
})

.here <- commandArgs(FALSE) |>
  str_subset("^--file=") |>
  first() |>
  str_remove("^--file=") |>
  dirname()
source(file.path(.here, "common.R"))
params <- load_params()

RANK <- params$scope$primary_rank            # "species"
say("Building TSE at rank: %s", RANK)

# ---- assays ----------------------------------------------------------------
mpath <- function(kind) file.path(ROOT, "results", "matrices",
                                  sprintf("%s_%s.tsv", kind, RANK))
tpm    <- read_matrix(mpath("tpm_microbial"))
reads  <- read_matrix(mpath("evidence_reads"))
ncont  <- read_matrix(mpath("evidence_ncontigs"))

# align evidence assays to the abundance assay
reads <- reads[rownames(tpm), colnames(tpm), drop = FALSE]
ncont <- ncont[rownames(tpm), colnames(tpm), drop = FALSE]
say("  assay dims: %d taxa x %d samples", nrow(tpm), ncol(tpm))

stopifnot(identical(dim(tpm), dim(reads)), identical(dim(tpm), dim(ncont)))

# ---- rowData ---------------------------------------------------------------
# Matrix rows are rank LABELS, not taxids -- several taxids collapse into one
# label at genus/family level, and unresolved bins have no taxid at all. So
# rowData is built by collapsing the lineage table to unique labels.
tax <- data.table::fread(file.path(ROOT, "results", "taxonomy", "taxon_lineage.tsv"),
                         na.strings = "") |>
  as_tibble()
rank_cols <- c("domain", "kingdom", "phylum", "class", "order", "family",
               "genus", "species")
idx <- match(RANK, rank_cols)

label_of <- function(df) {
  lab <- as.character(df[[RANK]])
  miss <- is.na(lab) | lab == ""
  for (up in rank_cols[seq_len(idx - 1)][rev(seq_len(idx - 1))]) {
    if (!any(miss)) break
    anc <- as.character(df[[up]])
    fill <- miss & !is.na(anc) & anc != ""
    lab[fill] <- paste0(anc[fill], " (unresolved)")
    miss <- miss & !fill
  }
  lab[miss] <- "unresolved"
  lab
}
tax <- tax |> mutate(label = label_of(pick(everything())))

# Restrict to microbial taxa BEFORE deduplicating by label.
#
# NCBI reuses names across kingdoms -- see results/taxonomy/homonyms.tsv. In
# this dataset "Bacillus" is both a bacterial genus and a stick-insect genus,
# and "Planococcus" is both a bacterium and a mealybug. The abundance matrix
# only ever contains microbial taxa, so deduplicating over the full lineage
# table could attach an Arthropoda lineage to a bacterial row. Filtering first
# makes that impossible.
tax <- tax |> filter(scope_class == "microbial")

# one row per label; keep the lineage down to RANK plus the annotations.
# Joining onto the matrix labels puts rowData in matrix row order, and a label
# with no lineage match keeps its name with NA lineage.
keep_ranks <- rank_cols[seq_len(idx)]
rd <- tax |>
  filter(label %in% rownames(tpm)) |>
  select(label, all_of(keep_ranks), superkingdom_group, scope_class,
         contaminant_flag, contaminant_confidence) |>
  distinct(label, .keep_all = TRUE)
rd <- tibble(label = rownames(tpm)) |>
  left_join(rd, by = "label")
rd <- as.data.frame(rd); rownames(rd) <- rownames(tpm)
say("  rowData: %d rows, %d contaminant-flagged",
    nrow(rd), sum(rd$contaminant_flag %in% 1, na.rm = TRUE))
stopifnot(!any(rd$scope_class %in% c("host", "bloodmeal")))
if (any(is.na(rd$scope_class))) {
  miss <- rownames(rd)[is.na(rd$scope_class)]
  say("  WARN: %d matrix label(s) had no lineage match: %s",
      length(miss), paste(head(miss, 5), collapse = ", "))
  say("        (expected for the 'unresolved' bin, which is a placeholder")
  say("         rather than a taxon and so has no lineage)")
}

# ---- taxon removal ---------------------------------------------------------
# See config taxon_removal. Unidentifiable bins and named contaminant genera are
# removed here, once, so every downstream step sees the same taxon set. TPM is
# renormalised to 1e6 over the remaining taxa.
rm_cfg <- params$taxon_removal
is_blank <- function(v) is.na(v) | v == ""
placeholder <- grepl(rm_cfg$placeholder_pattern, rd$species, ignore.case = TRUE)
drop_unres <- isTRUE(rm_cfg$drop_unresolved) & is_blank(rd$genus) &
  (is_blank(rd$species) | placeholder)
drop_excl <- rd$genus %in% unlist(rm_cfg$exclude_genera)
drop <- drop_unres | drop_excl
removed <- data.frame(
  label = rownames(rd)[drop],
  reason = ifelse(drop_excl[drop], "excluded_genus", "unresolved"),
  total_tpm = rowSums(tpm[drop, , drop = FALSE]),
  n_samples = rowSums(tpm[drop, , drop = FALSE] > 0)
)
lost <- colSums(tpm[drop, , drop = FALSE]) / pmax(colSums(tpm), 1)
say("  removed %d unresolved + %d excluded-genus taxa (%s)",
    sum(drop_unres & !drop_excl), sum(drop_excl),
    paste(unlist(rm_cfg$exclude_genera), collapse = ", "))
say("    median %.1f%% (max %.1f%%) of per-sample microbial TPM removed",
    100 * median(lost), 100 * max(lost))
tpm <- tpm[!drop, , drop = FALSE]
reads <- reads[!drop, , drop = FALSE]
ncont <- ncont[!drop, , drop = FALSE]
rd <- rd[!drop, , drop = FALSE]
cs0 <- colSums(tpm)
tpm <- sweep(tpm, 2, ifelse(cs0 > 0, cs0 / 1e6, 1), "/")

# ---- colData ---------------------------------------------------------------
meta_path <- file.path(ROOT, "metadata", "sample_metadata.tsv")
if (file.exists(meta_path)) {
  meta <- data.table::fread(meta_path, na.strings = c("", "NA")) |> as_tibble()
} else {
  say("  NOTE: metadata/sample_metadata.tsv absent -- falling back to the")
  say("        TEMPLATE. colData will carry the id-string guesses only.")
  meta <- data.table::fread(file.path(ROOT, "metadata", "sample_metadata_TEMPLATE.tsv"),
                            na.strings = c("", "NA")) |>
    as_tibble()
}
qc <- data.table::fread(file.path(ROOT, "results", "tables", "fraction_summary.tsv")) |>
  as_tibble()
# every QC sample, then only the matrix samples, in matrix column order
cd <- right_join(meta, qc, by = "sample_id") |>
  relocate(sample_id) |>
  filter(sample_id %in% colnames(tpm)) |>
  arrange(match(sample_id, colnames(tpm)))
stopifnot(identical(as.character(cd$sample_id), colnames(tpm)))
cd <- as.data.frame(cd); rownames(cd) <- cd$sample_id
say("  colData: %d samples x %d fields", nrow(cd), ncol(cd))

# ---- assemble --------------------------------------------------------------
tse <- TreeSummarizedExperiment(
  assays  = list(tpm_microbial = tpm, reads = reads, n_contigs = ncont),
  rowData = rd,
  colData = cd
)
metadata(tse) <- list(
  project      = params$project$name,
  built        = as.character(Sys.Date()),
  rank         = RANK,
  params       = params,
  primary_assay = "tpm_microbial",
  evidence_assays = c("reads", "n_contigs"),
  caveats      = CAVEAT
)

outdir <- file.path(ROOT, "results", "tse"); dir.create(outdir, showWarnings = FALSE)
saveRDS(tse, file.path(outdir, "tick_microbiome_tse.rds"))
say("  wrote results/tse/tick_microbiome_tse.rds")

# ---- agglomerated views ----------------------------------------------------
# mia::agglomerateByRank needs the lineage columns present in rowData, which
# they are (down to RANK). Higher ranks are derived here rather than in Python
# so there is a single source of truth.
for (r in intersect(params$scope$ranks_to_export, keep_ranks)) {
  if (r == RANK) next
  ok <- tryCatch({
    # At genus rank, a named species with no NCBI genus (common for tick
    # viruses) stands in as its own genus instead of being binned into
    # "<family> (unresolved)" -- or, with no family either, "unresolved".
    src <- tse
    if (r == "genus") {
      g <- rowData(src)$genus
      nog <- is_blank(g) & !is_blank(rowData(src)$species)
      g[nog] <- rowData(src)$species[nog]
      rowData(src)$genus <- g
      say("  genus rank: %d genus-less named species kept under their own name", sum(nog))
    }
    ag <- agglomerateByRank(src, rank = r, na.rm = FALSE)

    # mia names agglomerated rows by their FULL lineage path
    # ("Bacteria_Pseudomonadati_..._Acanthamoeba"), which correctly keeps
    # homonyms apart but does not match the short labels used by
    # results/matrices/. Relabel to the short convention so the two are
    # interchangeable, disambiguating only where a short label would collide.
    rdg <- as.data.frame(rowData(ag))
    short <- if (r %in% names(rdg)) as.character(rdg[[r]]) else rownames(ag)
    blank <- is.na(short) | short == ""
    if (any(blank)) {
      for (up in rank_cols[seq_len(match(r, rank_cols) - 1)][rev(seq_len(match(r, rank_cols) - 1))]) {
        if (!any(blank) || !(up %in% names(rdg))) next
        anc <- as.character(rdg[[up]]); f <- blank & !is.na(anc) & anc != ""
        short[f] <- paste0(anc[f], " (unresolved)"); blank <- blank & !f
      }
      short[blank] <- "unresolved"
    }
    dup <- duplicated(short) | duplicated(short, fromLast = TRUE)
    if (any(dup)) {
      say("  note: %d %s rows share a short label (distinct lineages); suffixed",
          sum(dup), r)
      short[dup] <- paste0(short[dup], " [", rownames(ag)[dup], "]")
    }
    rownames(ag) <- short

    saveRDS(ag, file.path(outdir, sprintf("tick_microbiome_tse_%s.rds", r)))
    say("  wrote results/tse/tick_microbiome_tse_%s.rds (%d taxa)", r, nrow(ag))

    # cross-check against the Python matrix: per-sample totals must agree even
    # where the row counts differ, because both renormalise to 1e6.
    pyf <- file.path(ROOT, "results", "matrices", sprintf("tpm_microbial_%s.tsv", r))
    if (file.exists(pyf)) {
      py <- read_matrix(pyf)
      shared <- intersect(colnames(ag), colnames(py))
      d <- max(abs(colSums(assay(ag, "tpm_microbial")[, shared, drop = FALSE]) -
                   colSums(py[, shared, drop = FALSE])))
      say("    cross-check vs results/matrices/: %d taxa (R) vs %d (python), "
          , nrow(ag), nrow(py))
      say("    max per-sample TPM difference = %.4f", d)
    }
    TRUE
  }, error = function(e) { say("  WARN: agglomeration to %s failed: %s", r, conditionMessage(e)); FALSE })
}

# ---- flat export so the data is usable outside R ---------------------------
flat <- file.path(outdir, "flat"); dir.create(flat, showWarnings = FALSE)
write.table(rd, file.path(flat, "rowData.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cd, file.path(flat, "colData.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(removed, file.path(flat, "removed_taxa.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("  wrote results/tse/flat/{rowData,colData,removed_taxa}.tsv")

# ---- verification ----------------------------------------------------------
cs <- colSums(assay(tse, "tpm_microbial"))
say("")
say("VERIFY  colSums(tpm_microbial): min=%.1f max=%.1f", min(cs), max(cs))
say("        (expect 1e6; 0 means a sample has no microbial signal)")
bad <- names(cs)[abs(cs - 1e6) > 1 & cs > 0]
if (length(bad)) say("  WARN: %d samples deviate from 1e6: %s",
                     length(bad), paste(head(bad, 5), collapse = ", "))
say("VERIFY  zero-signal samples: %d", sum(cs == 0))
say("")
print(tse)
