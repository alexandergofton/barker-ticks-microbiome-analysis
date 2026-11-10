#!/usr/bin/env Rscript
# Step 8 -- descriptive composition, contamination diagnostics, cross-validation.
#
#   conda activate tick_micro
#   Rscript code/analysis/08_composition.R
#
# Everything here is DESCRIPTION, not inference. This study sequenced no
# extraction blanks and no no-template controls, so absolute composition cannot
# be validated. Reagent-contaminant-flagged taxa are annotated and plotted, not
# removed. Unidentifiable taxa and excluded genera (config taxon_removal) are
# already removed from the TSE at step 07.

# tidyverse is attached LAST so dplyr's verbs (filter, rename, desc, ...) are
# not masked by the Bioconductor packages.
suppressPackageStartupMessages({
  library(yaml)
  library(TreeSummarizedExperiment)
  library(mia)
  library(pheatmap)
  library(RColorBrewer)
  library(tidyverse)
})

.here <- commandArgs(FALSE) |>
  str_subset("^--file=") |>
  first() |>
  str_remove("^--file=") |>
  dirname()
source(file.path(.here, "common.R"))
params <- load_params()
# theme: theme_tick() from common.R (theme_classic + light-grey major gridlines)

tse <- readRDS(file.path(ROOT, "results", "tse", "tick_microbiome_tse.rds"))
say("loaded TSE: %d taxa x %d samples", nrow(tse), ncol(tse))

gpath <- file.path(ROOT, "results", "tse", "tick_microbiome_tse_genus.rds")
gtse <- if (file.exists(gpath)) readRDS(gpath) else tse
G <- assay(gtse, "tpm_microbial")
Greads <- assay(gtse, "reads")
rdG <- as.data.frame(rowData(gtse))
cd <- as.data.frame(colData(tse))

grp_var <- if (!all(is.na(cd$tick_species))) {
  "tick_species"
} else {
  "tick_species_guess"
}
say("grouping figures by: %s", grp_var)

# ---------------------------------------------------------------------------
# 1. Where does the sequence actually go? host vs microbial vs unassigned.
#    This is the first thing a reader of a whole-tick-homogenate study asks.
# ---------------------------------------------------------------------------
frac <- data.table::fread(file.path(
  ROOT,
  "results",
  "tables",
  "fraction_summary.tsv"
)) |>
  as_tibble()

# One fixed colour per category, shared with panel 1b below, so a category
# means the same colour whether or not no_lca is in the plot.
frac_cats <- c(
  "host",
  "bloodmeal",
  "microbial",
  "other_euk",
  "ambiguous",
  "no_lca",
  "filtered_out"
)
# CSIRO categorical slots in fixed order; no_lca (the unidentifiable bulk) is
# the neutral Mist so the eye goes to the identified fractions.
frac_pal <- setNames(c(csiro_cat(5), CSIRO_OTHER, CSIRO_MUTED), frac_cats)

# `fraction` is a factor in column order, so the legend and stacking follow
# frac_cats rather than the alphabet.
# Rows are grouped by fraction, which sets the order bars are drawn in.
fl <- frac |>
  select(
    sample_id,
    host = tpm_pct_host,
    bloodmeal = tpm_pct_bloodmeal,
    microbial = tpm_pct_microbial,
    other_euk = tpm_pct_other_eukaryote,
    ambiguous = tpm_pct_ambiguous,
    no_lca = tpm_pct_no_lca,
    filtered_out = tpm_pct_filtered_out
  ) |>
  pivot_longer(
    -sample_id,
    names_to = "fraction",
    values_to = "pct",
    names_transform = list(fraction = fct_inorder)
  ) |>
  arrange(fraction)
ord <- frac |> arrange(desc(tpm_pct_microbial)) |> pull(sample_id)
fl <- fl |> mutate(sample_id = factor(sample_id, levels = ord))

p <- ggplot(fl, aes(sample_id, pct, fill = fraction)) +
  geom_col(width = 1) +
  scale_fill_manual(values = frac_pal, name = "Taxon category") +
  labs(
    title = "Where the sequence goes: TPM fraction per sample",
    subtitle = "ordered by microbial share",
    x = NULL,
    y = "% of sample TPM",
    caption = strwrap(CAVEAT, 110) |> paste(collapse = "\n")
  ) +
  scale_y_continuous(
    #limits = c(0, 100),
    expand = c(0, 0)
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

save_plot(
  fig_path("01_sequence_fractions.png"),
  p,
  width = 11,
  height = 5,
  dpi = 150
)
say("  wrote 01_sequence_fractions.png")

# ---------------------------------------------------------------------------
# 1b. Same view, but with no-LCA contigs (no BLAST hit at all -- see
#     results/audit diagnostics) dropped and the remaining categories
#     renormalised to 100%. Panel 1 asks "where does ALL the sequence go";
#     this asks "where does IDENTIFIABLE sequence go". Same sample order and
#     colours as panel 1, so the two are directly comparable.
# ---------------------------------------------------------------------------
assigned_cats <- setdiff(frac_cats, "no_lca")
fbl <- frac |>
  select(
    sample_id,
    host = tpm_pct_host,
    bloodmeal = tpm_pct_bloodmeal,
    microbial = tpm_pct_microbial,
    other_euk = tpm_pct_other_eukaryote,
    ambiguous = tpm_pct_ambiguous,
    filtered_out = tpm_pct_filtered_out,
    no_lca = tpm_pct_no_lca
  ) |>
  mutate(
    denom = 100 - no_lca,
    across(all_of(assigned_cats), \(x) if_else(denom > 0, x / denom * 100, 0))
  ) |>
  select(sample_id, all_of(assigned_cats)) |>
  pivot_longer(
    -sample_id,
    names_to = "fraction",
    values_to = "pct",
    names_transform = list(fraction = fct_inorder)
  ) |>
  arrange(fraction) |>
  mutate(sample_id = factor(sample_id, levels = ord))

p <- ggplot(fbl, aes(sample_id, pct, fill = fraction)) +
  geom_col(width = 1) +
  scale_fill_manual(values = frac_pal, name = "Taxon category") +
  labs(
    title = "Where IDENTIFIABLE sequence goes: TPM fraction per sample",
    subtitle = "no-LCA contigs excluded, remaining fractions renormalised to 100%; same order as panel 1",
    x = NULL,
    y = "% of identifiable-contig TPM",
    caption = strwrap(CAVEAT, 110) |> paste(collapse = "\n")
  ) +
  scale_y_continuous(
    #limits = c(0, 100),
    expand = c(0, 0)
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank()
  )
save_plot(
  fig_path("01b_sequence_fractions_assigned.png"),
  p,
  width = 11,
  height = 5,
  dpi = 150
)
say("  wrote 01b_sequence_fractions_assigned.png")

# ---------------------------------------------------------------------------
# 1c. Same view again, but also dropping filtered_out (contigs that DID get an
#     LCA call but failed the evidence filter -- config/analysis_params.yaml
#     min_contig_len/min_mapped_reads/min_top_bitscore) on top of no_lca.
#     These are exactly the two categories that never reach tpm_microbial /
#     the TSE (see code/analysis/__build_matrices.py: filtered_out fails the
#     `keep` mask before scope is even attached; ambiguous fails the
#     is_microbial mask a few lines later). ambiguous is left in here since it
#     IS a scope-classification outcome, not a filtering one -- it means the
#     LCA taxid has no scope_class in the lineage table, not that the contig
#     lacked evidence. Same sample order and colours as panels 1/1b.
# ---------------------------------------------------------------------------
kept_cats <- setdiff(assigned_cats, "filtered_out")
fcl <- frac |>
  select(
    sample_id,
    host = tpm_pct_host,
    bloodmeal = tpm_pct_bloodmeal,
    microbial = tpm_pct_microbial,
    other_euk = tpm_pct_other_eukaryote,
    ambiguous = tpm_pct_ambiguous,
    no_lca = tpm_pct_no_lca,
    filtered_out = tpm_pct_filtered_out
  ) |>
  mutate(
    denom_c = 100 - no_lca - filtered_out,
    across(all_of(kept_cats), \(x) if_else(denom_c > 0, x / denom_c * 100, 0))
  ) |>
  select(sample_id, all_of(kept_cats)) |>
  pivot_longer(
    -sample_id,
    names_to = "fraction",
    values_to = "pct",
    names_transform = list(fraction = fct_inorder)
  ) |>
  arrange(fraction) |>
  mutate(sample_id = factor(sample_id, levels = ord))

p <- ggplot(fcl, aes(sample_id, pct, fill = fraction)) +
  geom_col(width = 1) +
  scale_fill_manual(values = frac_pal, name = "Taxon category") +
  labs(
    title = "Where RETAINED sequence goes: TPM fraction per sample",
    subtitle = "no-LCA and filtered_out contigs excluded, remaining fractions renormalised to 100%; same order as panel 1",
    x = NULL,
    y = "% of retained-contig TPM",
    caption = strwrap(CAVEAT, 110) |> paste(collapse = "\n")
  ) +
  scale_y_continuous(
    #limits = c(0, 100),
    expand = c(0, 0)
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank()
  )

save_plot(
  fig_path("01c_sequence_fractions_retained.png"),
  p,
  width = 11,
  height = 5,
  dpi = 150
)
say("  wrote 01c_sequence_fractions_retained.png")

# ---------------------------------------------------------------------------
# 2. Top genera, stacked, faceted by tick species.
# ---------------------------------------------------------------------------
# Nine named genera + "other": the most a categorical palette can carry and
# stay distinguishable under colour-vision deficiency (common.R).
topn <- length(CSIRO_CAT)
gm <- rowMeans(G)
top <- names(sort(gm, decreasing = TRUE))[seq_len(min(topn, length(gm)))]
genus_long <- G |>
  as_tibble(rownames = "taxon") |>
  pivot_longer(-taxon, names_to = "sample_id", values_to = "tpm") |>
  mutate(taxon2 = if_else(taxon %in% top, taxon, "other")) |>
  summarise(tpm = sum(tpm), .by = c(sample_id, taxon2)) |>
  left_join(
    as_tibble(cd) |> select(sample_id, grp = all_of(grp_var)),
    by = "sample_id"
  ) |>
  arrange(sample_id) |>
  mutate(grp = if_else(is.na(grp) | grp == "", "unassigned", grp))

pal <- c(csiro_cat(length(top)), CSIRO_OTHER)
names(pal) <- c(top, "other")
genus_long <- genus_long |>
  mutate(taxon2 = factor(taxon2, levels = c(top, "other")))
p <- ggplot(genus_long, aes(sample_id, tpm / 1e4, fill = taxon2)) +
  geom_col(width = 1) +
  facet_grid(~grp, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal, name = "genus") +
  labs(
    title = sprintf(
      "Top %d microbial genera (TPM, renormalised within the microbial fraction)",
      topn
    ),
    x = NULL,
    y = "% of microbial TPM",
    caption = strwrap(CAVEAT, 110) |> paste(collapse = "\n")
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(angle = 90, hjust = 0, size = 6),
    legend.text = element_text(size = 7)
  )
save_plot(
  fig_path("02_top_genera_by_species.png"),
  p,
  width = 14,
  height = 7,
  dpi = 150
)
say("  wrote 02_top_genera_by_species.png")

# ---------------------------------------------------------------------------
# 3. Prevalence x abundance, coloured by contaminant flag.
#    With no blanks this is the main contamination diagnostic available:
#    reagent taxa are high-prevalence, and their abundance is inversely
#    related to sample DNA input.
# ---------------------------------------------------------------------------
prev <- rowSums(Greads >= params$taxon_filter$min_reads_per_taxon) / ncol(G)
mean_ab <- rowMeans(G)
pv <- tibble(
  taxon = rownames(G),
  prevalence = prev,
  mean_tpm = mean_ab,
  contaminant = rdG$contaminant_flag %in% 1,
  confidence = rdG$contaminant_confidence
) |>
  mutate(
    label = if_else(
      prevalence > 0.5 | mean_tpm > quantile(mean_tpm, 0.999),
      taxon,
      ""
    )
  )
data.table::fwrite(
  pv |> arrange(desc(prevalence)),
  tab_path("genus_prevalence_abundance.tsv"),
  sep = "\t"
)

p <- ggplot(pv, aes(prevalence, mean_tpm + 1, colour = contaminant)) +
  geom_point(aes(shape = contaminant), alpha = 0.7, size = 1.6) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
  ggrepel::geom_text_repel(
    aes(label = label),
    size = 2.4,
    max.overlaps = 25,
    show.legend = FALSE
  ) +
  scale_y_log10() +
  scale_colour_manual(
    values = c(
      `FALSE` = unname(CSIRO["Steel"]),
      `TRUE` = unname(CSIRO["Orange"])
    ),
    name = "reagent-contaminant\nflagged"
  ) +
  labs(
    title = "Genus prevalence vs mean abundance",
    subtitle = "flagged taxa are annotated from metadata/reagent_contaminants.tsv; none are removed",
    x = "prevalence (fraction of samples above the read threshold)",
    y = "mean TPM + 1 (log10)",
    caption = strwrap(CAVEAT, 110) |> paste(collapse = "\n")
  )
save_plot(
  fig_path("03_prevalence_vs_abundance.png"),
  p,
  width = 9,
  height = 6.5,
  dpi = 150
)
say("  wrote 03_prevalence_vs_abundance.png")

# ---------------------------------------------------------------------------
# 4. decontam frequency method -- runs only if DNA concentrations were supplied.
#    It needs concentrations, not blanks, so it is the one decontam mode this
#    study can actually support.
# ---------------------------------------------------------------------------
if ("dna_conc_ng_ul" %in% names(cd) && sum(!is.na(cd$dna_conc_ng_ul)) >= 20) {
  if (requireNamespace("decontam", quietly = TRUE)) {
    conc <- as.numeric(cd$dna_conc_ng_ul)
    ok <- !is.na(conc) & conc > 0
    say(
      "running decontam frequency method on %d samples with DNA conc",
      sum(ok)
    )
    ct <- decontam::isContaminant(
      t(G[, ok, drop = FALSE]),
      conc = conc[ok],
      method = "frequency"
    )
    ct <- ct |> mutate(taxon = rownames(ct))
    data.table::fwrite(
      as_tibble(ct) |> arrange(p),
      tab_path("decontam_frequency_scores.tsv"),
      sep = "\t"
    )
    say(
      "  wrote decontam_frequency_scores.tsv (%d flagged at p<0.1)",
      sum(ct$contaminant, na.rm = TRUE)
    )
    say("  NOTE: scores are reported as an annotation. No taxa are removed.")
  } else {
    say("decontam not installed -- skipping frequency method")
  }
} else {
  say(
    "dna_conc_ng_ul not supplied for >=20 samples -- decontam frequency method skipped."
  )
  say("  Supplying DNA concentrations would enable the one contamination test")
  say("  this design can support without blanks.")
}

# ---------------------------------------------------------------------------
# 5. Genus x sample heatmap.
# ---------------------------------------------------------------------------
hm_top <- names(sort(gm, decreasing = TRUE))[seq_len(min(20, length(gm)))]
hm <- log10(G[hm_top, , drop = FALSE] + 1)
# Columns annotated by tick GENUS (the fixed TICK_GENUS_PAL), not species:
# ~30 species cannot be told apart by colour; the genus palette can.
ann_r <- data.frame(
  contaminant = ifelse(rdG[hm_top, "contaminant_flag"] %in% 1, "yes", "no"),
  row.names = hm_top
)
ann_c <- data.frame(
  tick_genus = as.character(tick_genus_factor(cd[colnames(G), "tick_genus"])),
  row.names = colnames(G)
)

# Group samples (columns) by tick genus, in TICK_GENUS_PAL order; within each
# genus, samples are ordered by hierarchical clustering (dendrogram not drawn).
tg_f <- tick_genus_factor(cd[colnames(hm), "tick_genus"])
col_ord <- unlist(
  lapply(split(colnames(hm), tg_f, drop = TRUE), function(s) {
    if (length(s) < 3) {
      return(s)
    }
    s[hclust(dist(t(hm[, s, drop = FALSE])))$order]
  }),
  use.names = FALSE
)
hm <- hm[, col_ord, drop = FALSE]
ann_c <- ann_c[col_ord, , drop = FALSE]
hm_gaps <- head(cumsum(table(droplevels(tg_f))), -1)

ann_cols <- list(
  contaminant = c(no = CSIRO_OTHER, yes = unname(CSIRO["Orange"])),
  tick_genus = TICK_GENUS_PAL[intersect(
    names(TICK_GENUS_PAL),
    unique(ann_c$tick_genus)
  )]
)
# Colour: ~80% of cells are 0 (not detected). With the white-anchored ramp
# those read as blank page, so 0 gets its own bin in CSIRO_OTHER (Mist, the
# "not detected" colour) and detected values (log10 TPM >= ~1) use the CSIRO
# blue ramp from its pale-blue step, so every cell is filled.
hm_max <- max(1.01, max(hm))
hm_breaks <- c(0, 0.5, seq(1, hm_max, length.out = 100))
hm_cols <- c(CSIRO_OTHER, grDevices::colorRampPalette(CSIRO_SEQ[-1])(99))
# 1500 x 900 px at 130 dpi
png_svg(
  fig_path("04_genus_heatmap.png"),
  width = 1500 / 130,
  height = 900 / 130,
  res = 130,
  draw = function() {
    pheatmap(
      hm,
      #annotation_row = ann_r,
      annotation_col = ann_c,
      annotation_colors = ann_cols,
      show_colnames = FALSE,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      gaps_col = hm_gaps,
      fontsize = 10,
      fontsize_row = 9,
      border_color = NA,
      #cellwidth = 25,
      #cellheight = 18,
      color = hm_cols,
      breaks = hm_breaks,
      #main = "log10(TPM+1), top genera",
      treeheight_col = 0,
      treeheight_row = 0
    )
  }
)
say("  wrote 04_genus_heatmap.png")

# ---------------------------------------------------------------------------
# 6. Cross-check against the independent read-level Bracken track.
#    Assembly+BLAST+LCA and read-level Kraken2/Bracken are methodologically
#    independent. Agreement is the strongest validation available in the
#    absence of controls; disagreement is itself a result worth reporting.
# ---------------------------------------------------------------------------
kdir <- file.path(ROOT, "..", "kraken2")
bfiles <- list.files(kdir, pattern = "\\.bracken\\.report$", full.names = TRUE)
if (length(bfiles)) {
  say("cross-checking against %d Bracken reports", length(bfiles))
  # bracken.report is a kraken-style report: pct, clade, direct, rank, taxid, name
  bg <- bfiles |>
    map(function(f) {
      d <- tryCatch(
        data.table::fread(f, header = FALSE, sep = "\t", quote = "") |>
          as_tibble(),
        error = function(e) NULL
      )
      if (is.null(d) || ncol(d) < 6) {
        return(NULL)
      }
      d <- d |>
        rename_with(
          \(x) c("pct", "clade", "direct", "rank", "taxid", "name"),
          .cols = 1:6
        ) |>
        filter(rank == "G")
      # strip the run suffix: SAMPLE_FLOWCELL_INDEX_LANE.bracken.report
      s <- basename(f) |>
        str_remove("\\.bracken\\.report$") |>
        str_remove("_[A-Z0-9]{9,}_[ACGT]+-[ACGT]+_L\\d{3}$")
      tibble(sample_id = s, genus = str_trim(d$name), bracken_pct = d$pct)
    }) |>
    list_rbind()

  shared <- intersect(unique(bg$sample_id), colnames(G))
  say("  %d samples shared with the LCA matrix", length(shared))
  if (length(shared) >= 5) {
    rho <- map_dbl(
      shared,
      function(s) {
        b <- bg |> filter(sample_id == s)
        gv <- G[, s]
        common <- intersect(b$genus, names(gv))
        if (length(common) < 5) {
          return(NA_real_)
        }
        suppressWarnings(cor(
          b$bracken_pct[match(common, b$genus)],
          gv[common],
          method = "spearman"
        ))
      }
    )
    res <- tibble(sample_id = shared, spearman_rho = rho)
    data.table::fwrite(
      res,
      tab_path("lca_vs_bracken_concordance.tsv"),
      sep = "\t"
    )
    say(
      "  median Spearman rho = %.3f (n=%d samples with >=5 shared genera)",
      median(rho, na.rm = TRUE),
      sum(!is.na(rho))
    )
    p <- ggplot(res |> filter(!is.na(spearman_rho)), aes(spearman_rho)) +
      geom_histogram(bins = 30, fill = CSIRO_CAT[1], colour = "white") +
      geom_vline(
        xintercept = median(rho, na.rm = TRUE),
        linetype = 2,
        colour = unname(CSIRO["Orange"])
      ) +
      labs(
        title = "Concordance: assembly LCA vs read-level Bracken",
        subtitle = "per-sample Spearman correlation across shared genera; dashed line = median",
        x = "Spearman rho",
        y = "samples"
      )
    save_plot(
      fig_path("05_lca_vs_bracken.png"),
      p,
      width = 7,
      height = 5,
      dpi = 150
    )
    say("  wrote 05_lca_vs_bracken.png")
  } else {
    say("  too few shared samples for a concordance check")
  }
} else {
  say("no Bracken reports found at %s -- cross-check skipped", kdir)
}

say("\ndone")
