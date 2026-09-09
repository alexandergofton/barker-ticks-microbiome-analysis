#!/usr/bin/env Rscript
# Step 9 -- alpha and beta diversity, and the batch check.
#
#   conda activate tick_micro
#   Rscript code/analysis/09_diversity.R
#
# CAVEATS BUILT INTO THE OUTPUT, not just the caption:
#   * TPM is compositional. Richness on assembly-derived TPM is bounded by how
#     deeply each library was sequenced and how well it assembled, so figure 06
#     plots richness AGAINST contigs-passing-filter to make that confound
#     visible rather than leaving it implicit.
#   * With no controls, the batch terms (flowcell, lane, seq_submission) are the
#     only handle on technical variation. A significant batch term in the
#     PERMANOVA is a finding to report, not a nuisance to absorb.

suppressPackageStartupMessages({
  library(data.table); library(yaml)
  library(TreeSummarizedExperiment); library(mia); library(vegan)
  library(ggplot2)
})

.here <- dirname(sub("^--file=", "",
                     commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
source(file.path(.here, "common.R"))
params <- load_params()
theme_set(theme_bw(base_size = 10))

gpath <- file.path(ROOT, "results", "tse", "tick_microbiome_tse_genus.rds")
tse <- readRDS(if (file.exists(gpath)) gpath else
               file.path(ROOT, "results", "tse", "tick_microbiome_tse.rds"))
G  <- assay(tse, "tpm_microbial")
cd <- as.data.frame(colData(tse))
say("diversity on %d taxa x %d samples", nrow(G), ncol(G))

keep <- colSums(G) > 0
if (any(!keep)) say("dropping %d samples with no microbial signal", sum(!keep))
G <- G[, keep, drop = FALSE]; cd <- cd[keep, , drop = FALSE]
X <- t(G)                                   # vegan wants samples x taxa

# ---------------------------------------------------------------------------
# Alpha diversity
# ---------------------------------------------------------------------------
alpha <- data.table(
  sample_id   = rownames(X),
  observed    = rowSums(X > 0),
  shannon     = vegan::diversity(X, "shannon"),
  invsimpson  = vegan::diversity(X, "invsimpson"),
  contigs_passing = cd$contigs_passing_filter,
  reads_mapped    = cd$reads_mapped_total,
  pct_microbial   = cd$tpm_pct_microbial
)
grp_var <- if (!all(is.na(cd$tick_species))) "tick_species" else "tick_species_guess"
alpha[, grp := as.character(cd[[grp_var]])]
alpha[is.na(grp) | grp == "", grp := "unassigned"]
fwrite(alpha, tab_path("alpha_diversity.tsv"), sep = "\t")
say("  wrote alpha_diversity.tsv")

# The confound, plotted rather than hidden.
rho <- cor(alpha$observed, alpha$contigs_passing, method = "spearman", use = "complete")
p <- ggplot(alpha, aes(contigs_passing, observed)) +
  geom_point(alpha = 0.7, colour = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, colour = "firebrick", linewidth = 0.5) +
  scale_x_log10() +
  labs(title = "Observed richness is bounded by assembly depth",
       subtitle = sprintf(paste("Spearman rho = %.2f. Richness differences between samples",
                                "partly reflect sequencing/assembly depth, not biology."), rho),
       x = "contigs passing the evidence filter (log10)", y = "observed genera",
       caption = strwrap(CAVEAT, 110) |> paste(collapse = "\n"))
ggsave(fig_path("06_richness_vs_depth.png"), p, width = 7.5, height = 5.5, dpi = 150)
say("  wrote 06_richness_vs_depth.png  (rho = %.2f)", rho)

al <- melt(alpha, id.vars = c("sample_id", "grp"),
           measure.vars = c("observed", "shannon", "invsimpson"))
p <- ggplot(al, aes(reorder(grp, value, median), value)) +
  geom_boxplot(outlier.size = 0.5, fill = "grey90") +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.5, colour = "steelblue") +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  coord_flip() +
  labs(title = "Alpha diversity by tick species", x = NULL, y = NULL,
       caption = strwrap(paste(CAVEAT, "Richness is depth-confounded -- see figure 06."),
                         110) |> paste(collapse = "\n"))
ggsave(fig_path("07_alpha_by_species.png"), p, width = 8, height = 10, dpi = 150)
say("  wrote 07_alpha_by_species.png")

# ---------------------------------------------------------------------------
# Beta diversity: Bray-Curtis and Aitchison (CLR)
# ---------------------------------------------------------------------------
pseudo <- params$diversity$clr_pseudocount
clr <- t(apply(X + pseudo, 1, function(v) log(v) - mean(log(v))))

dists <- list(
  bray      = vegan::vegdist(X, method = "bray"),
  aitchison = dist(clr)
)

for (nm in names(dists)) {
  d <- dists[[nm]]
  pc <- cmdscale(d, k = 2, eig = TRUE)
  ev <- round(100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]), 1)
  pd <- data.table(sample_id = rownames(X), Axis1 = pc$points[, 1], Axis2 = pc$points[, 2],
                   grp = alpha$grp, submission = cd$seq_submission, flowcell = cd$flowcell)
  fwrite(pd, tab_path(sprintf("ordination_%s.tsv", nm)), sep = "\t")

  p <- ggplot(pd, aes(Axis1, Axis2, colour = grp)) +
    geom_point(size = 2, alpha = 0.8) +
    labs(title = sprintf("PCoA -- %s distance", nm),
         x = sprintf("Axis 1 (%.1f%%)", ev[1]), y = sprintf("Axis 2 (%.1f%%)", ev[2]),
         colour = "tick species",
         caption = strwrap(CAVEAT, 110) |> paste(collapse = "\n")) +
    theme(legend.text = element_text(size = 6))
  ggsave(fig_path(sprintf("08_pcoa_%s_species.png", nm)), p, width = 9, height = 6.5, dpi = 150)

  # the same ordination coloured by batch -- if these separate, that is the story
  p <- ggplot(pd, aes(Axis1, Axis2, colour = submission)) +
    geom_point(size = 2, alpha = 0.8) +
    labs(title = sprintf("PCoA -- %s distance, coloured by sequencing submission", nm),
         subtitle = "separation here would indicate a technical, not biological, gradient",
         x = sprintf("Axis 1 (%.1f%%)", ev[1]), y = sprintf("Axis 2 (%.1f%%)", ev[2]))
  ggsave(fig_path(sprintf("09_pcoa_%s_batch.png", nm)), p, width = 8, height = 6, dpi = 150)
  say("  wrote 08/09_pcoa_%s_*.png", nm)
}

# ---------------------------------------------------------------------------
# PERMANOVA -- biology AND batch, in one table
# ---------------------------------------------------------------------------
terms <- intersect(params$diversity$permanova_terms, names(cd))
terms <- c(terms, intersect("seq_submission", names(cd)))
# Until metadata/sample_metadata.tsv is filled in, tick_species is empty and
# only the technical terms are testable. Fall back to the id-string guess so a
# preliminary biological signal is visible -- clearly labelled as provisional,
# because the guess is derived from the sample name, not an identification.
if ("tick_species" %in% names(cd) && all(is.na(cd$tick_species) | cd$tick_species == "") &&
    "tick_species_guess" %in% names(cd)) {
  say("NOTE: tick_species is empty -- using tick_species_guess (PROVISIONAL).")
  terms <- c(terms, "tick_species_guess")
}
# 0.5, not 0.7: the id-string species guess covers ~59% of samples, and a
# preliminary result on 59% with the n stated is more useful than no result at
# all. Coverage is reported per term so a thin one is never mistaken for a
# complete one.
MIN_COVERAGE <- 0.5
usable <- terms[vapply(terms, function(t) {
  v <- as.character(cd[[t]]); v <- v[!is.na(v) & v != ""]
  length(v) >= MIN_COVERAGE * nrow(cd) && length(unique(v)) > 1 &&
    length(unique(v)) < length(v)
}, logical(1))]
for (t in setdiff(terms, usable)) {
  v <- as.character(cd[[t]]); nb <- sum(!is.na(v) & v != "")
  say("  skipping term '%s': %d/%d samples populated, %d levels",
      t, nb, nrow(cd), length(unique(v[!is.na(v) & v != ""])))
}

if (length(usable)) {
  say("PERMANOVA terms usable: %s", paste(usable, collapse = ", "))
  res <- rbindlist(lapply(names(dists), function(nm) {
    rbindlist(lapply(usable, function(t) {
      v <- as.character(cd[[t]]); ok <- !is.na(v) & v != ""
      if (sum(ok) < 10 || length(unique(v[ok])) < 2) return(NULL)
      dsub <- as.dist(as.matrix(dists[[nm]])[ok, ok])
      a <- tryCatch(vegan::adonis2(dsub ~ v[ok],
                                   permutations = params$diversity$permanova_permutations),
                    error = function(e) NULL)
      if (is.null(a)) return(NULL)
      data.table(distance = nm, term = t, n = sum(ok),
                 coverage_pct = round(100 * sum(ok) / nrow(cd), 1),
                 levels = length(unique(v[ok])),
                 R2 = a$R2[1], F = a$F[1], p = a$`Pr(>F)`[1])
    }))
  }))
  res[, kind := ifelse(term %in% c("flowcell", "lane", "seq_submission"),
                       "TECHNICAL",
                ifelse(term == "tick_species_guess", "biological_PROVISIONAL",
                       "biological"))]
  setorder(res, distance, -R2)
  fwrite(res, tab_path("permanova.tsv"), sep = "\t")
  say("  wrote permanova.tsv")
  print(res)
  tech <- res[kind == "TECHNICAL" & p < 0.05]
  if (nrow(tech)) {
    say("")
    say("  ATTENTION: a technical term is significant (%s).",
        paste(unique(tech$term), collapse = ", "))
    say("  With no controls in this design, report this alongside the biological")
    say("  results rather than adjusting it away.")
  }
} else say("no usable PERMANOVA terms -- fill in metadata/sample_metadata.tsv")

say("\ndone")
