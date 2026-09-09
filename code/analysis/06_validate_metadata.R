#!/usr/bin/env Rscript
# Step 6 -- validate metadata/sample_metadata.tsv before it reaches the TSE.
#
# Run this after filling in the template. It fails loudly on the errors that
# would otherwise surface as silent NAs in colData or as samples quietly
# dropped from every figure.
#
#   Rscript code/analysis/06_validate_metadata.R
#
# Checks
#   * every matrix column has a metadata row, and vice versa
#   * sample_id unique
#   * controlled vocabularies for life_stage / state
#   * collection_date parses as a date
#   * numeric fields are numeric and in range
#   * tick_species disagreeing with the id-string guess -> REVIEW list, not an
#     error (the id string is only a guess; a mitogenome id overrides it)
#
# Output: results/audit/metadata_validation.txt

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

.here <- dirname(sub("^--file=", "",
                     commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
source(file.path(.here, "common.R"))
root <- ROOT

meta_path <- file.path(root, "metadata", "sample_metadata.tsv")
tmpl_path <- file.path(root, "metadata", "sample_metadata_TEMPLATE.tsv")

if (!file.exists(meta_path)) {
  stop(sprintf(paste0(
    "metadata/sample_metadata.tsv not found.\n",
    "  Fill in %s and save it without the _TEMPLATE suffix."), tmpl_path))
}

meta <- fread(meta_path, na.strings = c("", "NA"))
msg  <- character()
errs <- character()
note <- function(...) msg  <<- c(msg,  sprintf(...))
err  <- function(...) errs <<- c(errs, sprintf(...))

note("metadata rows: %d", nrow(meta))

# ---- required columns ------------------------------------------------------
required <- c("sample_id", "tick_species", "life_stage", "n_ticks_pooled",
              "host_species", "collection_site", "state", "collection_date")
missing_cols <- setdiff(required, names(meta))
if (length(missing_cols)) err("missing required columns: %s",
                              paste(missing_cols, collapse = ", "))

# ---- uniqueness ------------------------------------------------------------
dup <- meta$sample_id[duplicated(meta$sample_id)]
if (length(dup)) err("duplicate sample_id: %s", paste(unique(dup), collapse = ", "))

# ---- agreement with the abundance matrix -----------------------------------
mat_path <- file.path(root, "results", "matrices", "tpm_microbial_species.tsv")
if (file.exists(mat_path)) {
  hdr <- names(fread(mat_path, nrows = 0))
  mat_samples <- setdiff(hdr, c("label", "V1", ""))
  no_meta <- setdiff(mat_samples, meta$sample_id)
  no_data <- setdiff(meta[status == "complete"]$sample_id, mat_samples)
  note("samples in matrix: %d", length(mat_samples))
  if (length(no_meta)) err("in matrix but not in metadata (%d): %s",
                           length(no_meta), paste(head(no_meta, 10), collapse = ", "))
  if (length(no_data)) note("complete in audit but absent from matrix (%d): %s",
                            length(no_data), paste(head(no_data, 10), collapse = ", "))
} else {
  note("matrix not built yet -- skipping matrix/metadata cross-check")
}

# ---- controlled vocabularies ----------------------------------------------
ok_stage <- c("larva", "nymph", "female", "male", "mixed", "unknown")
bad <- setdiff(na.omit(unique(meta$life_stage)), ok_stage)
if (length(bad)) err("unrecognised life_stage: %s  (allowed: %s)",
                     paste(bad, collapse = ", "), paste(ok_stage, collapse = "/"))

ok_state <- c("QLD", "NSW", "VIC", "SA", "WA", "NT", "TAS", "ACT", "unknown")
bad <- setdiff(na.omit(unique(meta$state)), ok_state)
if (length(bad)) err("unrecognised state: %s", paste(bad, collapse = ", "))

# ---- types -----------------------------------------------------------------
if ("collection_date" %in% names(meta)) {
  d <- suppressWarnings(as.Date(as.character(meta$collection_date)))
  badd <- meta$sample_id[!is.na(meta$collection_date) & is.na(d)]
  if (length(badd)) err("collection_date does not parse (want YYYY-MM-DD): %s",
                        paste(head(badd, 10), collapse = ", "))
}
for (col in c("n_ticks_pooled", "dna_conc_ng_ul", "latitude", "longitude")) {
  if (!col %in% names(meta)) next
  v <- suppressWarnings(as.numeric(as.character(meta[[col]])))
  badn <- meta$sample_id[!is.na(meta[[col]]) & is.na(v)]
  if (length(badn)) err("%s is not numeric for: %s", col,
                        paste(head(badn, 10), collapse = ", "))
}
if ("latitude" %in% names(meta)) {
  v <- suppressWarnings(as.numeric(meta$latitude))
  if (any(!is.na(v) & (v < -45 | v > -9)))
    err("latitude outside the Australian range (-45..-9)")
}

# ---- completeness of the manual fields ------------------------------------
for (col in required[-1]) {
  if (!col %in% names(meta)) next
  n <- sum(is.na(meta[[col]]) | meta[[col]] == "")
  if (n > 0) note("%-18s blank for %d/%d samples", col, n, nrow(meta))
}

# ---- species review list (not an error) ------------------------------------
if (all(c("tick_species", "tick_species_guess") %in% names(meta))) {
  cmp <- meta[!is.na(tick_species) & !is.na(tick_species_guess) &
                tick_species != tick_species_guess]
  if (nrow(cmp)) {
    note("")
    note("REVIEW -- tick_species differs from the id-string guess for %d samples.", nrow(cmp))
    note("  This is expected wherever a real identification supersedes the guess.")
    for (i in seq_len(min(nrow(cmp), 20)))
      note("    %-34s supplied=%-30s guess=%s",
           cmp$sample_id[i], cmp$tick_species[i], cmp$tick_species_guess[i])
  }
}

# ---- report ----------------------------------------------------------------
out <- c("Metadata validation", strrep("=", 60), msg)
if (length(errs)) {
  out <- c(out, "", "ERRORS", strrep("-", 60), paste0("  ", errs))
} else {
  out <- c(out, "", "No errors.")
}
dir.create(file.path(root, "results", "audit"), showWarnings = FALSE, recursive = TRUE)
writeLines(out, file.path(root, "results", "audit", "metadata_validation.txt"))
cat(paste(out, collapse = "\n"), "\n")
if (length(errs)) quit(status = 1)
