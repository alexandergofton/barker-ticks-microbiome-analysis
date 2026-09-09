#!/usr/bin/env Rscript
# Step 8 -- descriptive composition, contamination diagnostics, cross-validation.
#
#   conda activate tick_micro
#   Rscript code/analysis/08_composition.R
#
# Everything here is DESCRIPTION, not inference. This study sequenced no
# extraction blanks and no no-template controls, so absolute composition cannot
# be validated. Contaminant-flagged taxa are annotated and plotted, never
# removed -- removal decisions stay with the analyst, per-analysis.

suppressPackageStartupMessages({
  library(data.table); library(yaml)
  library(TreeSummarizedExperiment); library(mia)
  library(ggplot2); library(pheatmap); library(RColorBrewer)
})

.here <- dirname(sub("^--file=", "",
                     commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
source(file.path(.here, "common.R"))
params <- load_params()
theme_set(theme_bw(base_size = 10))

tse <- readRDS(file.path(ROOT, "results", "tse", "tick_microbiome_tse.rds"))
say("loaded TSE: %d taxa x %d samples", nrow(tse), ncol(tse))

gpath <- file.path(ROOT, "results", "tse", "tick_microbiome_tse_genus.rds")
gtse  <- if (file.exists(gpath)) readRDS(gpath) else tse
G     <- assay(gtse, "tpm_microbial")
Greads <- assay(gtse, "reads")
rdG   <- as.data.frame(rowData(gtse))
cd    <- as.data.frame(colData(tse))

grp_var <- if (!all(is.na(cd$tick_species))) "tick_species" else "tick_species_guess"
say("grouping figures by: %s", grp_var)

# ---------------------------------------------------------------------------
# 1. Where does the sequence actually go? host vs microbial vs unassigned.
#    This is the first thing a reader of a whole-tick-homogenate study asks.
# ---------------------------------------------------------------------------
frac <- fread(file.path(ROOT, "results", "tables", "fraction_summary.tsv"))

# One fixed colour per category, shared with panel 1b below, so a category
# means the same colour whether or not no_lca is in the plot.
frac_cats <- c("host", "bloodmeal", "microbial", "other_euk", "ambiguous",
               "no_lca", "filtered_out")
frac_pal <- setNames(brewer.pal(length(frac_cats), "Set2"), frac_cats)

fl <- melt(frac[, .(sample_id, host = tpm_pct_host, bloodmeal = tpm_pct_bloodmeal,
                    microbial = tpm_pct_microbial, other_euk = tpm_pct_other_eukaryote,
                    ambiguous = tpm_pct_ambiguous, no_lca = tpm_pct_no_lca,
                    filtered_out = tpm_pct_filtered_out)],
           id.vars = "sample_id", variable.name = "fraction", value.name = "pct")
ord <- frac[order(-tpm_pct_microbial)]$sample_id
fl[, sample_id := factor(sample_id, levels = ord)]

p <- ggplot(fl, aes(sample_id, pct, fill = fraction)) +
  geom_col(width = 1) +
  scale_fill_manual(values = frac_pal, name = "fraction") +
  labs(title = "Where the sequence goes: TPM fraction per sample",
       subtitle = "ordered by microbial share", x = NULL, y = "% of sample TPM",
       caption = strwrap(CAVEAT, 110) |> paste(collapse = "\n")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(fig_path("01_sequence_fractions.png"), p, width = 11, height = 5, dpi = 150)
say("  wrote 01_sequence_fractions.png")

# ---------------------------------------------------------------------------
# 1b. Same view, but with no-LCA contigs (no BLAST hit at all -- see
#     results/audit diagnostics) dropped and the remaining categories
#     renormalised to 100%. Panel 1 asks "where does ALL the sequence go";
#     this asks "where does IDENTIFIABLE sequence go". Same sample order and
#     colours as panel 1, so the two are directly comparable.
# ---------------------------------------------------------------------------
assigned_cats <- setdiff(frac_cats, "no_lca")
fb <- frac[, .(sample_id, host = tpm_pct_host, bloodmeal = tpm_pct_bloodmeal,
               microbial = tpm_pct_microbial, other_euk = tpm_pct_other_eukaryote,
               ambiguous = tpm_pct_ambiguous, filtered_out = tpm_pct_filtered_out,
               no_lca = tpm_pct_no_lca)]
denom <- 100 - fb$no_lca
fb[, (assigned_cats) := lapply(.SD, function(x) fifelse(denom > 0, x / denom * 100, 0)),
   .SDcols = assigned_cats]

fbl <- melt(fb[, c("sample_id", assigned_cats), with = FALSE],
            id.vars = "sample_id", variable.name = "fraction", value.name = "pct")
fbl[, sample_id := factor(sample_id, levels = ord)]

p <- ggplot(fbl, aes(sample_id, pct, fill = fraction)) +
  geom_col(width = 1) +
  scale_fill_manual(values = frac_pal, name = "fraction") +
  labs(title = "Where IDENTIFIABLE sequence goes: TPM fraction per sample",
       subtitle = "no-LCA contigs excluded, remaining fractions renormalised to 100%; same order as panel 1",
       x = NULL, y = "% of identifiable-contig TPM",
       caption = strwrap(CAVEAT, 110) |> paste(collapse = "\n")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(fig_path("01b_sequence_fractions_assigned.png"), p, width = 11, height = 5, dpi = 150)
say("  wrote 01b_sequence_fractions_assigned.png")

# ---------------------------------------------------------------------------
# 1c. Same view again, but also dropping filtered_out (contigs that DID get an
#     LCA call but failed the evidence filter -- config/analysis_params.yaml
#     min_contig_len/min_mapped_reads/min_top_bitscore) on top of no_lca.
#     These are exactly the two categories that never reach tpm_microbial /
#     the TSE (see code/analysis/04_build_matrices.py: filtered_out fails the
#     `keep` mask before scope is even attached; ambiguous fails the
#     is_microbial mask a few lines later). ambiguous is left in here since it
#     IS a scope-classification outcome, not a filtering one -- it means the
#     LCA taxid has no scope_class in the lineage table, not that the contig
#     lacked evidence. Same sample order and colours as panels 1/1b.
# ---------------------------------------------------------------------------
kept_cats <- setdiff(assigned_cats, "filtered_out")
fc <- frac[, .(sample_id, host = tpm_pct_host, bloodmeal = tpm_pct_bloodmeal,
               microbial = tpm_pct_microbial, other_euk = tpm_pct_other_eukaryote,
               ambiguous = tpm_pct_ambiguous, no_lca = tpm_pct_no_lca,
               filtered_out = tpm_pct_filtered_out)]
denom_c <- 100 - fc$no_lca - fc$filtered_out
fc[, (kept_cats) := lapply(.SD, function(x) fifelse(denom_c > 0, x / denom_c * 100, 0)),
   .SDcols = kept_cats]

fcl <- melt(fc[, c("sample_id", kept_cats), with = FALSE],
            id.vars = "sample_id", variable.name = "fraction", value.name = "pct")
fcl[, sample_id := factor(sample_id, levels = ord)]

p <- ggplot(fcl, aes(sample_id, pct, fill = fraction)) +
  geom_col(width = 1) +
  scale_fill_manual(values = frac_pal, name = "fraction") +
  labs(title = "Where RETAINED sequence goes: TPM fraction per sample",
       subtitle = "no-LCA and filtered_out contigs excluded, remaining fractions renormalised to 100%; same order as panel 1",
       x = NULL, y = "% of retained-contig TPM",
       caption = strwrap(CAVEAT, 110) |> paste(collapse = "\n")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(fig_path("01c_sequence_fractions_retained.png"), p, width = 11, height = 5, dpi = 150)
say("  wrote 01c_sequence_fractions_retained.png")

# ---------------------------------------------------------------------------
# 2. Top genera, stacked, faceted by tick species.
# ---------------------------------------------------------------------------
topn <- 20
gm <- rowMeans(G)
top <- names(sort(gm, decreasing = TRUE))[seq_len(min(topn, length(gm)))]
dt <- as.data.table(as.table(G))
setnames(dt, c("taxon", "sample_id", "tpm"))
dt[, taxon2 := ifelse(taxon %in% top, as.character(taxon), "other")]
dt <- dt[, .(tpm = sum(tpm)), by = .(sample_id, taxon2)]
dt <- merge(dt, as.data.table(cd)[, .(sample_id, grp = get(grp_var))],
            by = "sample_id", all.x = TRUE)
dt[is.na(grp) | grp == "", grp := "unassigned"]

pal <- c(colorRampPalette(brewer.pal(12, "Paired"))(length(top)), "grey80")
names(pal) <- c(top, "other")
p <- ggplot(dt, aes(sample_id, tpm / 1e4, fill = taxon2)) +
  geom_col(width = 1) +
  facet_grid(~ grp, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal, name = "genus") +
  labs(title = sprintf("Top %d microbial genera (TPM, renormalised within the microbial fraction)", topn),
       x = NULL, y = "% of microbial TPM",
       caption = strwrap(CAVEAT, 110) |> paste(collapse = "\n")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text.x = element_text(angle = 90, hjust = 0, size = 6),
        legend.text = element_text(size = 7))
ggsave(fig_path("02_top_genera_by_species.png"), p, width = 14, height = 7, dpi = 150)
say("  wrote 02_top_genera_by_species.png")

# ---------------------------------------------------------------------------
# 3. Prevalence x abundance, coloured by contaminant flag.
#    With no blanks this is the main contamination diagnostic available:
#    reagent taxa are high-prevalence, and their abundance is inversely
#    related to sample DNA input.
# ---------------------------------------------------------------------------
prev <- rowSums(Greads >= params$taxon_filter$min_reads_per_taxon) / ncol(G)
mean_ab <- rowMeans(G)
pv <- data.table(taxon = rownames(G), prevalence = prev, mean_tpm = mean_ab,
                 contaminant = rdG$contaminant_flag %in% 1,
                 confidence = rdG$contaminant_confidence)
pv[, label := ifelse(prevalence > 0.5 | mean_tpm > quantile(mean_tpm, 0.999), taxon, "")]
fwrite(pv[order(-prevalence)], tab_path("genus_prevalence_abundance.tsv"), sep = "\t")

p <- ggplot(pv, aes(prevalence, mean_tpm + 1, colour = contaminant)) +
  geom_point(alpha = 0.6, size = 1.4) +
  ggrepel::geom_text_repel(aes(label = label), size = 2.4, max.overlaps = 25,
                           show.legend = FALSE) +
  scale_y_log10() +
  scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "firebrick"),
                      name = "reagent-contaminant\nflagged") +
  labs(title = "Genus prevalence vs mean abundance",
       subtitle = "flagged taxa are annotated from metadata/reagent_contaminants.tsv; none are removed",
       x = "prevalence (fraction of samples above the read threshold)",
       y = "mean TPM + 1 (log10)",
       caption = strwrap(CAVEAT, 110) |> paste(collapse = "\n"))
ggsave(fig_path("03_prevalence_vs_abundance.png"), p, width = 9, height = 6.5, dpi = 150)
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
    say("running decontam frequency method on %d samples with DNA conc", sum(ok))
    ct <- decontam::isContaminant(t(G[, ok, drop = FALSE]), conc = conc[ok],
                                  method = "frequency")
    ct$taxon <- rownames(ct)
    fwrite(as.data.table(ct)[order(p)], tab_path("decontam_frequency_scores.tsv"), sep = "\t")
    say("  wrote decontam_frequency_scores.tsv (%d flagged at p<0.1)",
        sum(ct$contaminant, na.rm = TRUE))
    say("  NOTE: scores are reported as an annotation. No taxa are removed.")
  } else say("decontam not installed -- skipping frequency method")
} else {
  say("dna_conc_ng_ul not supplied for >=20 samples -- decontam frequency method skipped.")
  say("  Supplying DNA concentrations would enable the one contamination test")
  say("  this design can support without blanks.")
}

# ---------------------------------------------------------------------------
# 5. Genus x sample heatmap.
# ---------------------------------------------------------------------------
hm <- log10(G[top, , drop = FALSE] + 1)
ann_r <- data.frame(phylum = rdG[top, "phylum"],
                    contaminant = ifelse(rdG[top, "contaminant_flag"] %in% 1, "yes", "no"),
                    row.names = top)
ann_c <- data.frame(group = as.character(cd[colnames(G), grp_var]),
                    row.names = colnames(G))
ann_c$group[is.na(ann_c$group) | ann_c$group == ""] <- "unassigned"
png(fig_path("04_genus_heatmap.png"), width = 1500, height = 900, res = 130)
pheatmap(hm, annotation_row = ann_r, annotation_col = ann_c,
         show_colnames = FALSE, fontsize_row = 7,
         color = colorRampPalette(brewer.pal(9, "YlGnBu"))(50),
         main = "log10(TPM+1), top genera")
dev.off()
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
  bg <- rbindlist(lapply(bfiles, function(f) {
    d <- tryCatch(fread(f, header = FALSE, sep = "\t", quote = ""), error = function(e) NULL)
    if (is.null(d) || ncol(d) < 6) return(NULL)
    setnames(d, 1:6, c("pct", "clade", "direct", "rank", "taxid", "name"))
    d <- d[rank == "G"]
    # strip the run suffix: SAMPLE_FLOWCELL_INDEX_LANE.bracken.report
    s <- sub("\\.bracken\\.report$", "", basename(f))
    s <- sub("_[A-Z0-9]{9,}_[ACGT]+-[ACGT]+_L\\d{3}$", "", s)
    data.table(sample_id = s, genus = trimws(d$name), bracken_pct = d$pct)
  }), fill = TRUE)

  shared <- intersect(unique(bg$sample_id), colnames(G))
  say("  %d samples shared with the LCA matrix", length(shared))
  if (length(shared) >= 5) {
    rho <- vapply(shared, function(s) {
      b <- bg[sample_id == s]
      gv <- G[, s]
      common <- intersect(b$genus, names(gv))
      if (length(common) < 5) return(NA_real_)
      suppressWarnings(cor(b[match(common, genus)]$bracken_pct,
                           gv[common], method = "spearman"))
    }, numeric(1))
    res <- data.table(sample_id = shared, spearman_rho = rho)
    fwrite(res, tab_path("lca_vs_bracken_concordance.tsv"), sep = "\t")
    say("  median Spearman rho = %.3f (n=%d samples with >=5 shared genera)",
        median(rho, na.rm = TRUE), sum(!is.na(rho)))
    p <- ggplot(res[!is.na(spearman_rho)], aes(spearman_rho)) +
      geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
      geom_vline(xintercept = median(rho, na.rm = TRUE), linetype = 2, colour = "firebrick") +
      labs(title = "Concordance: assembly LCA vs read-level Bracken",
           subtitle = "per-sample Spearman correlation across shared genera; dashed line = median",
           x = "Spearman rho", y = "samples")
    ggsave(fig_path("05_lca_vs_bracken.png"), p, width = 7, height = 5, dpi = 150)
    say("  wrote 05_lca_vs_bracken.png")
  } else say("  too few shared samples for a concordance check")
} else say("no Bracken reports found at %s -- cross-check skipped", kdir)

say("\ndone")
