#!/usr/bin/env Rscript
# Step 10b -- compare microbial composition and endosymbionts between tick
# genera and species, using the curated metadata.
#
#   conda activate tick_micro
#   Rscript code/analysis/10b_tick_group_comparison.R
#
# Needs steps 7 (TSE) and 10 (target presence calls) to have run.
#
# WHAT IS COMPARED
#   1. who was sampled: tick species x genus x host type          fig 10
#   2. mean composition per tick genus / species (heatmaps)       fig 11a-d
#   3. beta diversity: PCoA, PERMANOVA (genus, species-within-
#      genus, host type within genus, a marginal model), PERMDISP,
#      pairwise genus tests                                       fig 12
#   4. endosymbionts ("keystone" symbionts): prevalence per
#      genus/species, Fisher tests, abundance when present,
#      dominant microbial genus per tick, co-occurrence            fig 13-15
#   5. indicator genera per tick genus (indicspecies)
#   6. geography: Mantel test of composition vs distance
#
# DESIGN LIMITS that shape every test here
#   * Tick genus and host type are confounded (Amblyomma is mostly on reptiles,
#     Haemaphysalis/Ixodes on mammals). Host effects are therefore tested WITHIN
#     genus (permutations restricted to genus blocks), and a marginal model
#     reports each term's unique contribution.
#   * Groups below `tick_groups.min_samples_per_group` are shown in figures but
#     not tested.
#   * PERMANOVA is sensitive to dispersion as well as location; PERMDISP is
#     reported beside every single-factor test so the two can be read together.
#   * Rickettsia, Francisella and Coxiella detections are GENUS level. Each genus
#     holds endosymbionts and pathogens; see README_markers.md / README_trees.md
#     for the per-species resolution.
#   * No controls: composition is descriptive (see CAVEAT).
#
# Tables: results/tables/tick_groups/   Figures: results/figures/10..15

# tidyverse is attached LAST so dplyr's verbs (filter, rename, desc, ...) are
# not masked by the Bioconductor packages.
suppressPackageStartupMessages({
  library(yaml)
  library(TreeSummarizedExperiment); library(mia); library(vegan)
  library(patchwork)
  library(tidyverse)
})

.here <- commandArgs(FALSE) |>
  str_subset("^--file=") |>
  first() |>
  str_remove("^--file=") |>
  dirname()
source(file.path(.here, "common.R"))
params <- load_params()
tg <- params$tick_groups
MIN_N  <- tg$min_samples_per_group
NPERM  <- params$diversity$permanova_permutations
set.seed(tg$seed)
# theme: theme_tick() from common.R (theme_classic + light-grey major gridlines)

tdir <- file.path(ROOT, "results", "tables", "tick_groups")
dir.create(tdir, showWarnings = FALSE, recursive = TRUE)
tpath <- function(n) file.path(tdir, n)
cap <- function(...) paste(strwrap(paste(..., CAVEAT), 120), collapse = "\n")

# ---------------------------------------------------------------------------
# 0. Data
# ---------------------------------------------------------------------------
gpath <- file.path(ROOT, "results", "tse", "tick_microbiome_tse_genus.rds")
if (!file.exists(gpath)) stop("genus TSE missing -- run step 7")
tse <- readRDS(gpath)
G   <- assay(tse, "tpm_microbial")
Gr  <- assay(tse, "reads")
cd  <- as_tibble(as.data.frame(colData(tse)))
if (!"tick_genus" %in% names(cd) || all(is.na(cd$tick_genus)))
  stop("colData has no tick_genus -- run step 5b, then step 7")

keep <- colSums(G) > 0
if (any(!keep)) say("dropping %d samples with no microbial signal", sum(!keep))
G <- G[, keep, drop = FALSE]; Gr <- Gr[, keep, drop = FALSE]; cd <- filter(cd, keep)
stopifnot(identical(cd$sample_id, colnames(G)))
RA <- G / 1e4                                     # % of microbial TPM

cd <- cd |>
  mutate(genus   = tick_genus_factor(tick_genus),
         host    = host_type_factor(host_type),
         species = if_else(is.na(tick_species) | tick_species == "",
                           "unidentified", tick_species),
         sp_ok   = species_identified %in% "yes")
# group counts keep groups in order of first appearance (`.by`), as the
# GEN_TEST / SP_TEST order sets the order of every test below
gen_n <- cd |> summarise(N = n(), .by = genus)
sp_n  <- cd |> filter(sp_ok) |> summarise(N = n(), .by = species)
GEN_TEST <- gen_n |> filter(genus != "unidentified", N >= MIN_N) |> pull(genus) |> as.character()
SP_TEST  <- sp_n |> filter(N >= MIN_N) |> pull(species)
say("samples analysed: %d | genera tested (n>=%d): %s", nrow(cd), MIN_N,
    paste(GEN_TEST, collapse = ", "))
say("species tested (n>=%d, identified): %d", MIN_N, length(SP_TEST))

# species ordered within genus (palette order), then by n
sp_order <- cd |>
  summarise(N = n(), .by = c(genus, species)) |>
  arrange(genus, desc(N), species)
cd <- cd |> mutate(species_f = factor(species, levels = unique(sp_order$species)))
lab_n <- function(x, n) sprintf("%s (n=%d)", x, n)

# ---------------------------------------------------------------------------
# 1. Who was sampled
# ---------------------------------------------------------------------------
counts <- cd |> summarise(N = n(), .by = c(genus, species_f, host))
data.table::fwrite(counts |> arrange(genus, species_f, host), tpath("sample_counts.tsv"), sep = "\t")

p1 <- ggplot(cd |> summarise(N = n(), .by = c(genus, species_f)),
             aes(N, fct_rev(species_f), fill = genus)) +
  geom_col(width = 0.75, colour = MARK_OUTLINE, linewidth = 0.2) +
  geom_vline(xintercept = MIN_N - 0.5, linetype = 3, colour = CSIRO_MUTED) +
  scale_fill_manual(values = TICK_GENUS_PAL, guide = "none") +
  facet_grid(genus ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "samples analysed", y = NULL,
       subtitle = sprintf("dotted line: n = %d, the minimum for a group to be tested", MIN_N)) +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0))
gh <- cd |> summarise(N = n(), .by = c(genus, host))
p2 <- ggplot(gh, aes(host, fct_rev(genus), fill = N)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = N, colour = N > 0.55 * max(N)), size = 3.2) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
  scale_fill_gradientn(colours = CSIRO_SEQ, guide = "none") +
  labs(x = "host type", y = NULL,
       subtitle = "tick genus x host type: the two are confounded")
p <- (p1 | p2) + plot_layout(widths = c(1.4, 1)) +
  plot_annotation(title = "Samples analysed, by tick genus, species and host type",
                  caption = cap(""))
save_plot(fig_path("10_samples_by_tick_and_host.png"), p, width = 13, height = 9, dpi = 150)
say("  wrote 10_samples_by_tick_and_host.png")

# ---------------------------------------------------------------------------
# 2. Mean composition per group -- heatmaps, not 20-colour stacked bars
# ---------------------------------------------------------------------------
# (lookup tables are named *_count so they can't be confused with the new columns)
genus_count   <- cd |> summarise(N = n(), .by = genus) |> arrange(genus)
species_count <- cd |> summarise(N = n(), .by = species_f)

# As in figure 04: an absent genus (mean 0) is filled CSIRO_OTHER (Mist, the
# "not detected" colour) via na.value rather than white, and detected values
# use the CSIRO blue ramp from its pale-blue step, so every cell is filled.
heat <- function(d, x, title, topn, facet = FALSE) {
  p <- ggplot(d, aes(.data[[x]], taxon, fill = if_else(mean_pct > 0, mean_pct, NA_real_))) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = label, colour = mean_pct > 0.4 * max(mean_pct)), size = 2.4) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
    scale_fill_gradientn(colours = CSIRO_SEQ[-1], na.value = CSIRO_OTHER,
                         values = c(0, 0.1, 0.4, 1), name = "mean %\nmicrobial TPM") +
    labs(x = NULL, y = NULL, title = title,
         subtitle = sprintf("top %d microbial genera by overall mean; cells show mean %% (>= 1%%)", topn)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (facet) p <- p + facet_grid(~ genus, scales = "free_x", space = "free_x") +
    theme(strip.text.x = element_text(angle = 90, hjust = 0))
  p
}

# One genus heatmap + one species heatmap for the top `topn` taxa. The figure
# height grows with the number of rows so the 50-taxon versions stay legible.
composition_heatmaps <- function(topn, fig_genus, fig_species, write_tables = FALSE) {
  top <- names(sort(rowMeans(RA), decreasing = TRUE))[seq_len(min(topn, nrow(RA)))]
  long <- RA[top, , drop = FALSE] |>
    as_tibble(rownames = "taxon") |>
    pivot_longer(-taxon, names_to = "sample_id", values_to = "pct") |>
    inner_join(cd |> select(sample_id, genus, species_f), by = "sample_id") |>
    arrange(sample_id) |>
    mutate(taxon = factor(taxon, levels = rev(top)))

  grp_mean <- function(by) {
    long |>
      summarise(mean_pct = mean(pct), prevalence_pct = 100 * mean(pct > 0), n = n(),
                .by = all_of(c(by, "taxon"))) |>
      mutate(label = if_else(mean_pct >= 1, sprintf("%.0f", mean_pct), ""))
  }
  mg <- grp_mean("genus")
  ms <- grp_mean(c("genus", "species_f"))
  if (write_tables) {
    data.table::fwrite(mg, tpath("mean_composition_by_genus.tsv"), sep = "\t")
    data.table::fwrite(ms, tpath("mean_composition_by_species.tsv"), sep = "\t")
  }
  h <- 2.5 + 0.22 * length(top)
  mg <- mg |>
    mutate(gl = factor(lab_n(genus, genus_count$N[match(genus, genus_count$genus)]),
                       levels = lab_n(genus_count$genus, genus_count$N)))
  save_plot(fig_path(fig_genus),
         heat(mg, "gl", "Mean microbial composition by tick genus", topn) + labs(caption = cap("")),
         width = 9, height = max(7, h), dpi = 150)
  ms <- ms |>
    mutate(sl = factor(lab_n(species_f, species_count$N[match(species_f, species_count$species_f)]),
                       levels = lab_n(levels(cd$species_f),
                                      species_count$N[match(levels(cd$species_f), species_count$species_f)])))
  save_plot(fig_path(fig_species),
         heat(ms, "sl", "Mean microbial composition by tick species", topn, facet = TRUE) +
           labs(caption = cap("Groups with n < 3 are single ticks: read as individuals, not species means.")),
         width = 17, height = max(8, h + 1), dpi = 150)
  say("  wrote %s, %s", fig_genus, fig_species)
}
composition_heatmaps(tg$top_taxa, "11a_composition_by_tick_genus.png",
                     "11b_composition_by_tick_species.png", write_tables = TRUE)
composition_heatmaps(tg$top_taxa_extended, "11c_composition_by_tick_genus_top50.png",
                     "11d_composition_by_tick_species_top50.png")

# ---------------------------------------------------------------------------
# 3. Beta diversity
# ---------------------------------------------------------------------------
X <- t(G)
clr <- t(apply(X + params$diversity$clr_pseudocount, 1, function(v) log(v) - mean(log(v))))
DIST <- list(bray = vegdist(X, "bray"), aitchison = dist(clr))
sub_d <- function(d, ids) as.dist(as.matrix(d)[ids, ids])

disp_p <- function(d, g) {
  bd <- tryCatch(betadisper(d, g), error = function(e) NULL)
  if (is.null(bd)) return(NA_real_)
  permutest(bd, permutations = NPERM)$tab$`Pr(>F)`[1]
}

perm_rows <- list()
add_perm <- function(test, dname, a, term_rows, n, levels, disp = NA_real_, note = "") {
  for (t in term_rows) {
    # vegan >= 2.6-6 reports a single-term model as one row named "Model"
    # unless `by =` is given; per-term rows are named after the term.
    r <- if (t %in% rownames(a)) t else "Model"
    stopifnot(r %in% rownames(a))
    perm_rows[[length(perm_rows) + 1]] <<- tibble(
      test = test, distance = dname, term = t, n = n, levels = levels,
      R2 = a[r, "R2"], F = a[r, "F"], p = a[r, "Pr(>F)"], dispersion_p = disp, note = note)
  }
}

for (dn in names(DIST)) {
  D <- DIST[[dn]]

  # 3a. tick genus
  s <- cd |> filter(genus %in% GEN_TEST)
  a <- adonis2(sub_d(D, s$sample_id) ~ genus, data = droplevels(s), permutations = NPERM)
  add_perm("tick genus", dn, a, "genus", nrow(s), length(GEN_TEST),
           disp_p(sub_d(D, s$sample_id), droplevels(s$genus)))

  # 3a'. sensitivity: the same test without ticks whose mitogenome places them
  #      in a different genus from the curated one (step 26b). If those ticks
  #      are mixed or mislabelled, the genus effect should hold or strengthen
  #      without them.
  mvf <- file.path(ROOT, "results", "mito", "tables", "mito_vs_metadata.tsv")
  if (file.exists(mvf)) {
    disc <- data.table::fread(mvf) |>
      as_tibble() |>
      filter(verdict == "GENUS_DISCORDANT") |>
      pull(sample_id)
    s2 <- cd |> filter(genus %in% GEN_TEST, !sample_id %in% disc)
    a <- adonis2(sub_d(D, s2$sample_id) ~ genus, data = droplevels(s2), permutations = NPERM)
    add_perm("tick genus, excl. mito-discordant", dn, a, "genus", nrow(s2), length(GEN_TEST),
             disp_p(sub_d(D, s2$sample_id), droplevels(s2$genus)),
             sprintf("sensitivity: %d genus-discordant ticks (step 26b) removed",
                     sum(cd$sample_id %in% disc)))
  }

  # 3b. tick species, all tested species (marginal)
  s <- cd |> filter(species %in% SP_TEST)
  a <- adonis2(sub_d(D, s$sample_id) ~ species, data = s, permutations = NPERM)
  add_perm("tick species", dn, a, "species", nrow(s), length(SP_TEST),
           disp_p(sub_d(D, s$sample_id), factor(s$species)))

  # 3c. species WITHIN genus: sequential model, species after genus, and
  #     permutations restricted to within-genus blocks
  a <- adonis2(sub_d(D, s$sample_id) ~ genus + species, data = droplevels(s),
               by = "terms", permutations = how(blocks = droplevels(s$genus), nperm = NPERM))
  add_perm("tick species | genus", dn, a, "species", nrow(s), length(SP_TEST),
           note = "species after genus; permutations within genus")

  # 3d. species within each genus separately
  for (g in GEN_TEST) {
    sg <- s |> filter(genus == g)
    if (n_distinct(sg$species) < 2) next
    a <- adonis2(sub_d(D, sg$sample_id) ~ species, data = sg, permutations = NPERM)
    add_perm(sprintf("species within %s", g), dn, a, "species", nrow(sg),
             n_distinct(sg$species), disp_p(sub_d(D, sg$sample_id), factor(sg$species)),
             paste(sort(unique(sg$species)), collapse = "; "))
  }

  # 3e. host type -- marginal, and within genus (the confound)
  hs <- cd |>
    filter(host != "unknown") |>
    filter(n() >= MIN_N, .by = host)
  a <- adonis2(sub_d(D, hs$sample_id) ~ host, data = droplevels(hs), permutations = NPERM)
  add_perm("host type", dn, a, "host", nrow(hs), n_distinct(hs$host),
           disp_p(sub_d(D, hs$sample_id), droplevels(hs$host)),
           "NOT adjusted for tick genus")
  hg <- hs |> filter(genus %in% GEN_TEST)
  a <- adonis2(sub_d(D, hg$sample_id) ~ host, data = droplevels(hg),
               permutations = how(blocks = droplevels(hg$genus), nperm = NPERM))
  add_perm("host type | genus", dn, a, "host", nrow(hg), n_distinct(hg$host),
           note = "permutations within tick genus")

  # 3f. marginal model: each term's unique contribution given the others
  mm <- cd |>
    filter(genus %in% GEN_TEST, host != "unknown", !is.na(seq_submission)) |>
    droplevels()
  a <- adonis2(sub_d(D, mm$sample_id) ~ genus + host + seq_submission, data = mm,
               by = "margin", permutations = NPERM)
  add_perm("marginal model", dn, a, c("genus", "host", "seq_submission"), nrow(mm), NA,
           note = "~ genus + host + seq_submission, by = margin")

  # 3g. life stage within genus (adults vs immatures etc.)
  ls <- cd |>
    filter(!is.na(life_stage), life_stage != "", genus %in% GEN_TEST) |>
    filter(n() >= MIN_N, .by = life_stage) |>
    droplevels()
  a <- adonis2(sub_d(D, ls$sample_id) ~ life_stage, data = ls,
               permutations = how(blocks = ls$genus, nperm = NPERM))
  add_perm("life stage | genus", dn, a, "life_stage", nrow(ls), n_distinct(ls$life_stage),
           note = "permutations within tick genus")
}
perm <- list_rbind(perm_rows)
data.table::fwrite(perm, tpath("permanova.tsv"), sep = "\t")
say("  wrote tick_groups/permanova.tsv")
perm |>
  transmute(test, distance, term, n, R2 = round(R2, 3), p, dispersion_p) |>
  print(n = Inf, width = Inf)

# pairwise genus PERMANOVA, BH-adjusted within each distance
pw <- map(names(DIST), function(dn) {
  cmb <- combn(GEN_TEST, 2, simplify = FALSE)
  map(cmb, function(pr) {
    s <- droplevels(filter(cd, genus %in% pr))
    a <- adonis2(sub_d(DIST[[dn]], s$sample_id) ~ genus, data = s, permutations = NPERM)
    tibble(distance = dn, genus_a = pr[1], genus_b = pr[2], n = nrow(s),
           R2 = a$R2[1], F = a$F[1], p = a$`Pr(>F)`[1])
  }) |>
    list_rbind() |>
    mutate(p_adj_BH = p.adjust(p, "BH"))
}) |>
  list_rbind()
data.table::fwrite(pw, tpath("permanova_pairwise_genus.tsv"), sep = "\t")
say("  wrote tick_groups/permanova_pairwise_genus.tsv")

# PCoA, coloured by tick genus, ellipses + direct labels per genus
for (dn in names(DIST)) {
  pc <- cmdscale(DIST[[dn]], k = 2, eig = TRUE)
  ev <- round(100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]), 1)
  pd <- tibble(sample_id = rownames(X), A1 = pc$points[, 1], A2 = pc$points[, 2]) |>
    inner_join(cd |> select(sample_id, genus, host, species), by = "sample_id") |>
    arrange(sample_id)
  data.table::fwrite(pd, tpath(sprintf("pcoa_%s.tsv", dn)), sep = "\t")
  cen <- pd |>
    filter(genus %in% GEN_TEST) |>
    summarise(A1 = median(A1), A2 = median(A2), .by = genus)
  st <- perm |> filter(distance == dn, test == "tick genus")
  p <- ggplot(pd, aes(A1, A2)) +
    stat_ellipse(data = filter(pd, genus %in% GEN_TEST), aes(colour = genus),
                 type = "t", level = 0.68, linewidth = 0.6, show.legend = FALSE) +
    geom_point(aes(fill = genus), shape = 21, size = 2.6, colour = MARK_OUTLINE,
               stroke = 0.3, alpha = 0.9) +
    # label text in dark ink (several CSIRO hues are unreadable as text); the
    # label box is tinted with the genus colour (alpha 0.3) to tie it to its
    # ellipse. Fill is precomputed so the alpha applies to the box only.
    ggrepel::geom_label_repel(data = cen, aes(label = genus), colour = "black",
                              fill = scales::alpha(TICK_GENUS_PAL[as.character(cen$genus)], 0.3),
                              size = 3, fontface = "italic", label.size = 0.2,
                              show.legend = FALSE, seed = 1) +
    scale_fill_manual(values = TICK_GENUS_PAL, name = "tick genus",
                      guide = guide_legend(override.aes = list(size = 3))) +
    scale_colour_manual(values = TICK_GENUS_PAL) +
    labs(title = sprintf("PCoA, %s distance, genus-level microbial composition", dn),
         subtitle = sprintf(paste0("tick genus: R2 = %.2f, p = %.3f (PERMDISP p = %.3f)\n",
                                   "ellipses: 68%% for genera with n >= %d"),
                            st$R2, st$p, st$dispersion_p, MIN_N),
         x = sprintf("Axis 1 (%.1f%%)", ev[1]), y = sprintf("Axis 2 (%.1f%%)", ev[2]),
         caption = cap(""))
  save_plot(fig_path(sprintf("12_pcoa_%s_tick_genus.png", dn)), p, width = 10, height = 7.5, dpi = 150)
}
say("  wrote 12_pcoa_*_tick_genus.png")

# ---------------------------------------------------------------------------
# 4. Endosymbionts
# ---------------------------------------------------------------------------
tp <- data.table::fread(file.path(ROOT, "results", "targets", "target_prevalence.tsv")) |>
  as_tibble()
ENDO <- tp |> filter(category %in% tg$symbiont_categories) |> pull(target) |> unique()
tp <- tp |> filter(target %in% ENDO, sample_id %in% cd$sample_id)
# every sample x endosymbiont, absent unless step 10 made a presence call
es <- expand_grid(sample_id = cd$sample_id, target = ENDO) |>
  arrange(sample_id, target) |>
  left_join(tp |> transmute(sample_id, target, present, pct = tpm_microbial / 1e4, category),
            by = c("sample_id", "target")) |>
  mutate(pct         = if_else(is.na(present), 0, pct),
         present     = if_else(is.na(present), 0L, present),
         pct_present = if_else(present == 1, pct, NA_real_)) |>
  inner_join(cd |> select(sample_id, genus, species, species_f, sp_ok, host), by = "sample_id")
seen <- es |>
  summarise(n = sum(present), .by = target) |>
  filter(n > 0) |>
  arrange(desc(n))
ENDO <- seen$target
es <- es |> filter(target %in% ENDO)
both <- tp |> filter(category == "both") |> pull(target) |> unique()
endo_lab <- function(x) if_else(x %in% both, paste0(x, "*"), x)
es <- es |> mutate(target_f = factor(endo_lab(target), levels = endo_lab(ENDO)))
say("endosymbiont targets detected: %s", paste(ENDO, collapse = ", "))
BOTH_NOTE <- paste("* genus contains pathogens as well as endosymbionts: a genus-level",
                   "call does not say which (see the marker-gene trees).")

wilson <- function(k, n, z = 1.96) {
  p <- k / n; d <- 1 + z^2 / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d)
}
prev_tab <- function(by) {
  es |>
    summarise(n_present = sum(present), n = n(),
              median_pct_when_present = median(pct_present, na.rm = TRUE),
              .by = all_of(c(by, "target", "target_f"))) |>
    mutate(prevalence_pct = 100 * n_present / n,
           ci = map2(n_present, n, wilson),
           ci_low_pct  = 100 * map_dbl(ci, 1),
           ci_high_pct = 100 * map_dbl(ci, 2)) |>
    select(-ci)
}
pg <- prev_tab("genus"); ps <- prev_tab(c("genus", "species_f"))
data.table::fwrite(select(pg, -target_f), tpath("endosymbiont_prevalence_by_genus.tsv"), sep = "\t")
data.table::fwrite(select(ps, -target_f), tpath("endosymbiont_prevalence_by_species.tsv"), sep = "\t")

prev_plot <- function(d, y, title, facet = FALSE) {
  d <- d |> mutate(lab = sprintf("%d/%d", n_present, n))
  p <- ggplot(d, aes(target_f, .data[[y]], fill = prevalence_pct)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = lab, colour = prevalence_pct > 55), size = 2.5) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
    scale_fill_gradientn(colours = CSIRO_SEQ, limits = c(0, 100),
                        name = "prevalence %") +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL, title = title,
         subtitle = "cells: samples with a presence call / samples (step-10 read + contig thresholds)",
         caption = cap(BOTH_NOTE)) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, face = "italic"))
  if (facet) p <- p + facet_grid(genus ~ ., scales = "free_y", space = "free_y") +
    theme(strip.text.y = element_text(angle = 0))
  p
}
pg <- pg |> mutate(gl = factor(lab_n(genus, n), levels = unique(lab_n(genus, n)[order(genus)])))
save_plot(fig_path("13a_endosymbiont_prevalence_by_genus.png"),
       prev_plot(pg, "gl", "Endosymbiont prevalence by tick genus"),
       width = 10, height = 5.5, dpi = 150)
ps <- ps |> mutate(sl = factor(lab_n(species_f, n), levels = unique(lab_n(species_f, n)[order(species_f)])))
save_plot(fig_path("13b_endosymbiont_prevalence_by_species.png"),
       prev_plot(ps, "sl", "Endosymbiont prevalence by tick species", facet = TRUE),
       width = 11, height = 13, dpi = 150)
say("  wrote 13a/13b endosymbiont prevalence")

# Fisher's exact tests on presence/absence
fisher_rows <- list()
fish <- function(test, target, x, y, note = "") {
  tb <- table(x, y)
  if (nrow(tb) < 2 || ncol(tb) < 2) return(invisible())
  p <- if (all(dim(tb) == 2)) fisher.test(tb)$p.value else
    fisher.test(tb, simulate.p.value = TRUE, B = tg$fisher_B)$p.value
  fisher_rows[[length(fisher_rows) + 1]] <<- tibble(
    test = test, target = target, n = length(x), groups = nrow(tb),
    n_present = sum(y == 1), p = p, note = note)
}
for (t in ENDO) {
  e <- es |> filter(target == t)
  s <- e |> filter(genus %in% GEN_TEST); fish("tick genus", t, droplevels(s$genus), s$present)
  s <- e |> filter(species %in% SP_TEST); fish("tick species", t, s$species, s$present)
  for (g in GEN_TEST) {
    s <- e |> filter(genus == g, species %in% SP_TEST)
    if (n_distinct(s$species) >= 2) fish(sprintf("species within %s", g), t, s$species, s$present)
    s <- e |>
      filter(genus == g, host != "unknown") |>
      filter(n() >= MIN_N, .by = host)
    if (n_distinct(s$host) >= 2) fish(sprintf("host type within %s", g), t, droplevels(s$host),
                                      s$present, paste(sort(unique(as.character(s$host))), collapse = " vs "))
  }
}
fisher <- list_rbind(fisher_rows) |>
  mutate(p_adj_BH = p.adjust(p, "BH")) |>
  arrange(p)
data.table::fwrite(fisher, tpath("endosymbiont_fisher_tests.tsv"), sep = "\t")
say("  wrote tick_groups/endosymbiont_fisher_tests.tsv (%d tests, %d at BH q < 0.05)",
    nrow(fisher), sum(fisher$p_adj_BH < 0.05))

# abundance when present, by genus; Kruskal-Wallis on all samples (zeros kept)
MINP <- tg$min_detections_for_abundance
ab_t <- es |>
  summarise(n = sum(present), .by = target) |>
  filter(n >= MINP) |>
  pull(target)
kw <- map(ab_t, function(t) {
  s <- es |> filter(target == t, genus %in% GEN_TEST)
  k <- kruskal.test(pct ~ droplevels(genus), data = s)
  tibble(target = t, n = nrow(s), statistic = unname(k$statistic), p = k$p.value)
}) |>
  list_rbind() |>
  mutate(p_adj_BH = p.adjust(p, "BH"))
data.table::fwrite(kw, tpath("endosymbiont_abundance_kruskal.tsv"), sep = "\t")
pa <- es |> filter(present == 1, target %in% ab_t)
p <- ggplot(pa, aes(genus, pct_present)) +
  geom_boxplot(outlier.shape = NA, fill = CSIRO_OTHER, colour = CSIRO_MUTED, width = 0.6) +
  geom_jitter(aes(fill = genus, shape = host), width = 0.15, size = 2, colour = MARK_OUTLINE,
              stroke = 0.3) +
  scale_fill_manual(values = TICK_GENUS_PAL, guide = "none") +
  scale_shape_manual(values = HOST_TYPE_SHAPE, name = "host type", drop = FALSE,
                     guide = guide_legend(override.aes = list(fill = CSIRO_MUTED, size = 3))) +
  scale_y_log10(labels = function(x) paste0(x, "%")) +
  facet_wrap(~ target_f, scales = "free_x") +
  labs(title = "Endosymbiont relative abundance in the ticks that carry it",
       subtitle = sprintf("%% of microbial TPM, samples with a presence call only; targets with >= %d detections", MINP),
       x = NULL, y = "% of microbial TPM (log10)", caption = cap(BOTH_NOTE)) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, face = "italic"),
        strip.text = element_text(face = "italic"))
save_plot(fig_path("14_endosymbiont_abundance_by_genus.png"), p, width = 12, height = 8, dpi = 150)
say("  wrote 14_endosymbiont_abundance_by_genus.png")

# dominant microbial genus per tick -- which symbiont "owns" each tick lineage
dom <- tibble(sample_id = colnames(RA),
              dominant = rownames(RA)[apply(RA, 2, which.max)],
              dominant_pct = apply(RA, 2, max)) |>
  inner_join(cd |> select(sample_id, genus, species, species_f, host), by = "sample_id") |>
  arrange(sample_id)
data.table::fwrite(dom |> arrange(genus, species, desc(dominant_pct)),
                   tpath("dominant_genus_per_sample.tsv"), sep = "\t")
TOPD <- tg$top_dominant
dtop <- dom |>
  summarise(N = n(), .by = dominant) |>
  arrange(desc(N)) |>
  slice_head(n = TOPD) |>
  pull(dominant)
dom <- dom |>
  mutate(dom_f = factor(if_else(dominant %in% dtop, dominant, "other"), levels = c(dtop, "other")))
dc <- dom |>
  summarise(N = n(), .by = c(genus, species_f, dom_f)) |>
  mutate(tot = sum(N), .by = species_f)
p <- ggplot(dc, aes(dom_f, species_f, fill = N / tot)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = N, colour = N / tot > 0.55), size = 2.6) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
  scale_fill_gradientn(colours = CSIRO_SEQ[-1], limits = c(0, 1),
                      labels = scales::percent, name = "share of\nthat species'\nticks") +
  scale_y_discrete(limits = rev) +
  facet_grid(genus ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Dominant microbial genus in each tick",
       subtitle = "count of ticks whose single most abundant microbial genus is the column",
       x = NULL, y = NULL, caption = cap("")) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, face = "italic"),
        strip.text.y = element_text(angle = 0))
save_plot(fig_path("15_dominant_genus_by_tick_species.png"), p, width = 10, height = 13, dpi = 150)
say("  wrote 15_dominant_genus_by_tick_species.png")

# co-occurrence: overall Fisher, and Mantel-Haenszel stratified by tick genus so
# a pair that only co-occurs because both favour the same genus is not reported
# as an association
co_t <- es |>
  summarise(n = sum(present), .by = target) |>
  filter(n >= MINP) |>
  pull(target)
# one row per sample, one 0/1 column per target
W <- es |>
  filter(target %in% co_t) |>
  select(sample_id, genus, target, present) |>
  pivot_wider(names_from = target, values_from = present) |>
  arrange(sample_id, genus)
co <- map(combn(co_t, 2, simplify = FALSE), function(pr) {
  a <- W[[pr[1]]]; b <- W[[pr[2]]]
  ft <- fisher.test(table(factor(a, 0:1), factor(b, 0:1)))
  ok <- W$genus %in% GEN_TEST
  arr <- table(factor(a[ok], 0:1), factor(b[ok], 0:1), droplevels(W$genus[ok]))
  arr <- arr[, , apply(arr, 3, sum) >= 2, drop = FALSE]
  mh <- tryCatch(mantelhaen.test(arr, exact = TRUE), error = function(e) NULL)
  tibble(target_a = pr[1], target_b = pr[2], n_both = sum(a & b),
         n_a_only = sum(a & !b), n_b_only = sum(!a & b), n_neither = sum(!a & !b),
         odds_ratio = unname(ft$estimate), p = ft$p.value,
         odds_ratio_genus_stratified = if (is.null(mh)) NA_real_ else unname(mh$estimate),
         p_genus_stratified = if (is.null(mh)) NA_real_ else mh$p.value)
}) |>
  list_rbind() |>
  mutate(p_adj_BH = p.adjust(p, "BH"),
         p_genus_stratified_adj_BH = p.adjust(p_genus_stratified, "BH"))
data.table::fwrite(co |> arrange(p), tpath("endosymbiont_cooccurrence.tsv"), sep = "\t")
say("  wrote tick_groups/endosymbiont_cooccurrence.tsv")

# ---------------------------------------------------------------------------
# 5. Indicator genera per tick genus
# ---------------------------------------------------------------------------
if (requireNamespace("indicspecies", quietly = TRUE)) {
  s <- cd |> filter(genus %in% GEN_TEST)
  pres <- Gr[, s$sample_id] >= params$taxon_filter$min_reads_per_taxon
  cand <- rownames(RA)[rowMeans(pres) >= tg$indicator_min_prevalence]
  mp <- indicspecies::multipatt(t(RA[cand, s$sample_id]), droplevels(s$genus),
                                func = "IndVal.g", duleg = TRUE,
                                control = how(nperm = NPERM))
  A <- mp$A; B <- mp$B
  ind <- as_tibble(mp$sign, rownames = "taxon")
  grp_cols <- str_subset(names(ind), "^s\\.")
  ind <- ind |>
    # each taxon's indicator group = the s.<genus> column holding its 1
    mutate(tick_genus = str_remove(grp_cols, "^s\\.")[apply(pick(all_of(grp_cols)), 1, which.max)]) |>
    transmute(taxon, tick_genus, indval = stat, p = p.value) |>
    mutate(p_adj_BH = p.adjust(p, "BH"),
           specificity_A = map2_dbl(taxon, tick_genus, \(t, g) A[t, g]),
           fidelity_B    = map2_dbl(taxon, tick_genus, \(t, g) B[t, g]),
           flagged_contaminant = rowData(tse)[taxon, "contaminant_flag"] %in% 1) |>
    arrange(p_adj_BH, desc(indval))
  data.table::fwrite(ind, tpath("indicator_genera_by_tick_genus.tsv"), sep = "\t")
  say("  wrote tick_groups/indicator_genera_by_tick_genus.tsv (%d candidate genera, %d at BH q < 0.05)",
      length(cand), sum(ind$p_adj_BH < 0.05, na.rm = TRUE))
} else say("indicspecies not installed -- indicator analysis skipped")

# ---------------------------------------------------------------------------
# 6. Geography: does composition track distance?
# ---------------------------------------------------------------------------
hav <- function(lat, lon) {
  r <- pi / 180; n <- length(lat); m <- matrix(0, n, n)
  for (i in seq_len(n)) {
    dl <- (lat - lat[i]) * r; dn <- (lon - lon[i]) * r
    a <- sin(dl / 2)^2 + cos(lat[i] * r) * cos(lat * r) * sin(dn / 2)^2
    m[i, ] <- 2 * 6371 * asin(pmin(1, sqrt(a)))
  }
  m
}
cd <- cd |>
  mutate(lat = suppressWarnings(as.numeric(latitude)),
         lon = suppressWarnings(as.numeric(longitude)))
geo_sets <- c(list(`all samples` = filter(cd, !is.na(lat))$sample_id),
              map(set_names(GEN_TEST), \(g) filter(cd, genus == g, !is.na(lat))$sample_id),
              map(set_names(SP_TEST), \(s) filter(cd, species == s, !is.na(lat))$sample_id))
mant <- map(names(geo_sets), function(nm) {
  ids <- geo_sets[[nm]]
  if (length(ids) < tg$mantel_min_samples) return(NULL)
  pos <- cd |> slice(match(ids, sample_id))
  gd <- as.dist(hav(pos$lat, pos$lon))
  map(names(DIST), function(dn) {
    m <- mantel(sub_d(DIST[[dn]], ids), gd, method = "spearman", permutations = NPERM)
    tibble(set = nm, distance = dn, n = length(ids), mantel_r = m$statistic, p = m$signif)
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  mutate(p_adj_BH = p.adjust(p, "BH"), .by = distance)
data.table::fwrite(mant, tpath("mantel_composition_vs_geography.tsv"), sep = "\t")
say("  wrote tick_groups/mantel_composition_vs_geography.tsv")
say("  NOTE: across all samples, distance is confounded with tick genus and host;")
say("        the within-genus and within-species rows are the interpretable ones.")

say("\ndone")
