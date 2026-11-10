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

# tidyverse is attached LAST so dplyr's verbs (filter, rename, desc, ...) are
# not masked by the Bioconductor packages.
suppressPackageStartupMessages({
  library(yaml)
  library(TreeSummarizedExperiment)
  library(mia)
  library(vegan)
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

gpath <- file.path(ROOT, "results", "tse", "tick_microbiome_tse_genus.rds")
tse <- readRDS(
  if (file.exists(gpath)) {
    gpath
  } else {
    file.path(ROOT, "results", "tse", "tick_microbiome_tse.rds")
  }
)
G <- assay(tse, "tpm_microbial")
cd <- as.data.frame(colData(tse))
say("diversity on %d taxa x %d samples", nrow(G), ncol(G))

keep <- colSums(G) > 0
if (any(!keep)) {
  say("dropping %d samples with no microbial signal", sum(!keep))
}
G <- G[, keep, drop = FALSE]
cd <- cd[keep, , drop = FALSE]
X <- t(G) # vegan wants samples x taxa

# ---------------------------------------------------------------------------
# Alpha diversity
# ---------------------------------------------------------------------------
grp_var <- if (!all(is.na(cd$tick_species))) {
  "tick_species"
} else {
  "tick_species_guess"
}
alpha <- tibble(
  sample_id = rownames(X),
  observed = rowSums(X > 0),
  shannon = vegan::diversity(X, "shannon"),
  invsimpson = vegan::diversity(X, "invsimpson"),
  contigs_passing = cd$contigs_passing_filter,
  reads_mapped = cd$reads_mapped_total,
  pct_microbial = cd$tpm_pct_microbial
) |>
  mutate(
    grp = as.character(cd[[grp_var]]),
    grp = if_else(is.na(grp) | grp == "", "unassigned", grp),
    genus = tick_genus_factor(cd$tick_genus)
  )
data.table::fwrite(alpha, tab_path("alpha_diversity.tsv"), sep = "\t")
say("  wrote alpha_diversity.tsv")

# The confound, plotted rather than hidden.
rho <- cor(
  alpha$observed,
  alpha$contigs_passing,
  method = "spearman",
  use = "complete"
)
p <- ggplot(alpha, aes(contigs_passing, observed)) +
  geom_point(alpha = 0.7, colour = CSIRO_CAT[1], size = 3) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = unname(CSIRO["Orange"]),
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  scale_x_log10() +
  labs(
    title = "Observed richness is bounded by assembly depth",
    subtitle = sprintf(
      paste(
        "Spearman rho = %.2f. Richness differences between samples",
        "partly reflect sequencing/assembly depth, not biology."
      ),
      rho
    ),
    x = "contigs passing the evidence filter (log10)",
    y = "observed genera",
    caption = strwrap(CAVEAT, 110) |> paste(collapse = "\n")
  )

save_plot(
  fig_path("06_richness_vs_depth.png"),
  p,
  width = 7.5,
  height = 5.5,
  dpi = 150
)
say("  wrote 06_richness_vs_depth.png  (rho = %.2f)", rho)

# Species represented by a single sample give a box with no spread, so they
# are left out of figure 07 (they stay in alpha_diversity.tsv).
grp_n <- count(alpha, grp)
singletons <- grp_n$grp[grp_n$n < 2]
say(
  "  figure 07: dropping %d single-sample species: %s",
  length(singletons),
  paste(singletons, collapse = ", ")
)

# `variable` is a factor in column order, so the facets read observed,
# shannon, invsimpson rather than alphabetically.
al <- alpha |>
  filter(!grp %in% singletons) |>
  select(sample_id, grp, genus, observed, shannon, invsimpson) |>
  pivot_longer(
    c(observed, shannon, invsimpson),
    names_to = "variable",
    values_to = "value",
    names_transform = list(variable = fct_inorder)
  )
p <- ggplot(al, aes(reorder(grp, value, median), value)) +
  geom_boxplot(outlier.shape = NA, fill = CSIRO_OTHER, colour = CSIRO_MUTED) +
  geom_jitter(
    aes(fill = genus),
    width = 0.15,
    size = 1.3,
    shape = 21,
    stroke = 0.2,
    colour = MARK_OUTLINE,
    alpha = 0.85
  ) +
  scale_fill_manual(values = TICK_GENUS_PAL, name = "tick genus", drop = TRUE) +
  facet_wrap(~variable, scales = "free_x", nrow = 1) +
  coord_flip() +
  labs(
    title = "Alpha diversity by tick species",
    subtitle = sprintf(
      "Species with only one sample omitted (%d species).",
      length(singletons)
    ),
    x = NULL,
    y = NULL,
    caption = strwrap(
      paste(CAVEAT, "Richness is depth-confounded -- see figure 06."),
      110
    ) |>
      paste(collapse = "\n")
  )

save_plot(
  fig_path("07_alpha_by_species.png"),
  p,
  width = 8,
  height = 10,
  dpi = 150
)

say("  wrote 07_alpha_by_species.png")

# ---------------------------------------------------------------------------
# Beta diversity: Bray-Curtis and Aitchison (CLR)
# ---------------------------------------------------------------------------
pseudo <- params$diversity$clr_pseudocount
clr <- t(apply(X + pseudo, 1, function(v) log(v) - mean(log(v))))

dists <- list(
  bray = vegan::vegdist(X, method = "bray"),
  aitchison = dist(clr)
)

for (nm in names(dists)) {
  d <- dists[[nm]]
  pc <- cmdscale(d, k = 2, eig = TRUE)
  ev <- round(100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]), 1)
  pd <- tibble(
    sample_id = rownames(X),
    Axis1 = pc$points[, 1],
    Axis2 = pc$points[, 2],
    grp = alpha$grp,
    genus = alpha$genus,
    submission = cd$seq_submission,
    flowcell = cd$flowcell
  )
  data.table::fwrite(pd, tab_path(sprintf("ordination_%s.tsv", nm)), sep = "\t")

  # ~30 species cannot share one colour scale and stay distinguishable, so:
  # one panel per tick genus, every other sample in grey behind, species
  # coloured WITHIN the panel (CSIRO slots in order of size, past nine folded
  # to Mist) and named directly at their centroid -- no legend to decode.
  # slot = rank of the species within its genus by sample count (ties in
  # order of appearance)
  pd <- pd |>
    mutate(n_sp = n(), .by = grp) |>
    mutate(slot = match(grp, unique(grp[order(-n_sp)])), .by = genus) |>
    mutate(col = if_else(
      slot <= length(CSIRO_CAT),
      CSIRO_CAT[pmin(slot, length(CSIRO_CAT))],
      CSIRO_OTHER
    ))
  cen <- pd |>
    summarise(
      Axis1 = median(Axis1),
      Axis2 = median(Axis2),
      lab = str_remove(first(grp), "^\\S+ "),
      .by = c(genus, grp)
    )
  p <- ggplot(pd, aes(Axis1, Axis2)) +
    geom_point(data = select(pd, -genus), colour = "grey85", size = 1) +
    geom_point(
      aes(fill = col),
      shape = 21,
      size = 2.2,
      stroke = 0.25,
      colour = MARK_OUTLINE
    ) +
    ggrepel::geom_text_repel(
      data = cen,
      aes(label = lab),
      size = 2.6,
      fontface = "italic",
      min.segment.length = 0,
      max.overlaps = 30
    ) +
    scale_fill_identity() +
    facet_wrap(~genus) +
    labs(
      title = sprintf(
        "PCoA -- %s distance, by tick species (one panel per genus)",
        nm
      ),
      x = sprintf("Axis 1 (%.1f%%)", ev[1]),
      y = sprintf("Axis 2 (%.1f%%)", ev[2]),
      caption = strwrap(CAVEAT, 110) |> paste(collapse = "\n")
    ) +
    theme(strip.text = element_text(face = "italic"))
  save_plot(
    fig_path(sprintf("08_pcoa_%s_species.png", nm)),
    p,
    width = 11,
    height = 9,
    dpi = 150
  )

  # the same ordination coloured by batch -- if these separate, that is the story
  nsub <- length(unique(pd$submission))
  p <- ggplot(pd, aes(Axis1, Axis2, fill = submission, shape = submission)) +
    geom_point(size = 2.2, alpha = 0.9, stroke = 0.25, colour = MARK_OUTLINE) +
    scale_fill_manual(values = csiro_cat(nsub)) +
    scale_shape_manual(values = c(21, 24, 22, 23, 25)[seq_len(nsub)]) +
    labs(
      title = sprintf(
        "PCoA -- %s distance, coloured by sequencing submission",
        nm
      ),
      subtitle = "separation here would indicate a technical, not biological, gradient",
      x = sprintf("Axis 1 (%.1f%%)", ev[1]),
      y = sprintf("Axis 2 (%.1f%%)", ev[2])
    )
  save_plot(
    fig_path(sprintf("09_pcoa_%s_batch.png", nm)),
    p,
    width = 8,
    height = 6,
    dpi = 150
  )
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
if (
  "tick_species" %in%
    names(cd) &&
    all(is.na(cd$tick_species) | cd$tick_species == "") &&
    "tick_species_guess" %in% names(cd)
) {
  say("NOTE: tick_species is empty -- using tick_species_guess (PROVISIONAL).")
  terms <- c(terms, "tick_species_guess")
}
# 0.5, not 0.7: the id-string species guess covers ~59% of samples, and a
# preliminary result on 59% with the n stated is more useful than no result at
# all. Coverage is reported per term so a thin one is never mistaken for a
# complete one.
MIN_COVERAGE <- 0.5
usable <- terms[map_lgl(
  terms,
  function(t) {
    v <- as.character(cd[[t]])
    v <- v[!is.na(v) & v != ""]
    length(v) >= MIN_COVERAGE * nrow(cd) &&
      length(unique(v)) > 1 &&
      length(unique(v)) < length(v)
  }
)]
for (t in setdiff(terms, usable)) {
  v <- as.character(cd[[t]])
  nb <- sum(!is.na(v) & v != "")
  say(
    "  skipping term '%s': %d/%d samples populated, %d levels",
    t,
    nb,
    nrow(cd),
    length(unique(v[!is.na(v) & v != ""]))
  )
}

if (length(usable)) {
  say("PERMANOVA terms usable: %s", paste(usable, collapse = ", "))
  res <- map(names(dists), function(nm) {
    map(usable, function(t) {
      v <- as.character(cd[[t]])
      ok <- !is.na(v) & v != ""
      if (sum(ok) < 10 || length(unique(v[ok])) < 2) {
        return(NULL)
      }
      dsub <- as.dist(as.matrix(dists[[nm]])[ok, ok])
      a <- tryCatch(
        vegan::adonis2(
          dsub ~ v[ok],
          permutations = params$diversity$permanova_permutations
        ),
        error = function(e) NULL
      )
      if (is.null(a)) {
        return(NULL)
      }
      tibble(
        distance = nm,
        term = t,
        n = sum(ok),
        coverage_pct = round(100 * sum(ok) / nrow(cd), 1),
        levels = length(unique(v[ok])),
        R2 = a$R2[1],
        F = a$F[1],
        p = a$`Pr(>F)`[1]
      )
    }) |>
      list_rbind()
  }) |>
    list_rbind()
  res <- res |>
    mutate(kind = case_when(
      term %in% c("flowcell", "lane", "seq_submission") ~ "TECHNICAL",
      term == "tick_species_guess" ~ "biological_PROVISIONAL",
      .default = "biological"
    )) |>
    arrange(distance, desc(R2))
  data.table::fwrite(res, tab_path("permanova.tsv"), sep = "\t")
  say("  wrote permanova.tsv")
  print(res, n = Inf, width = Inf)
  tech <- res |> filter(kind == "TECHNICAL", p < 0.05)
  if (nrow(tech)) {
    say("")
    say(
      "  ATTENTION: a technical term is significant (%s).",
      paste(unique(tech$term), collapse = ", ")
    )
    say(
      "  With no controls in this design, report this alongside the biological"
    )
    say("  results rather than adjusting it away.")
  }
} else {
  say("no usable PERMANOVA terms -- fill in metadata/sample_metadata.tsv")
}

say("\ndone")
