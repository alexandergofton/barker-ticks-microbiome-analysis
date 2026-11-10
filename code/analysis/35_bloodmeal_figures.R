#!/usr/bin/env Rscript
# Step 35 -- blood-meal figures.
#
#   Rscript code/analysis/35_bloodmeal_figures.R
#
# Reads results/bloodmeal/tables/bloodmeal_calls.tsv and bloodmeal_host_taxa.tsv
# (step 34) and draws:
#
#   bloodmeal_concordance.png   per metadata host class: how far the sequence
#                               call agrees with the recorded host (stacked
#                               counts; match depth in one blue ramp, mismatch
#                               in orange, untestable in grey)
#   bloodmeal_class_matrix.png  recorded host class x detected class (counts)
#   bloodmeal_per_sample.png    every sample: detected host families (point
#                               size = share of that sample's host reads), with
#                               the recorded host's family ringed
#
# All figures go to results/bloodmeal/figures/, png + svg.

suppressMessages({
  library(tidyverse)
})

.args <- commandArgs(trailingOnly = FALSE)
.this <- .args |> str_subset("^--file=") |> str_remove("^--file=")
source(file.path(dirname(normalizePath(.this)), "common.R"))
params <- load_params()
bm_dir <- file.path(ROOT, params$bloodmeal$out_dir)
fig_dir <- file.path(bm_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
read_t <- function(f) {
  as_tibble(data.table::fread(
    file.path(bm_dir, "tables", f),
    na.strings = c("", "NA")
  ))
}

calls <- read_t("bloodmeal_calls.tsv")
taxa <- read_t("bloodmeal_host_taxa.tsv")

CLASS_LEVELS <- c("Mammalia", "Aves", "Lepidosauria", "other", "none recorded")
class_lab <- function(x) {
  x <- if_else(is.na(x), "none recorded", x)
  factor(if_else(x %in% CLASS_LEVELS, x, "other"), levels = CLASS_LEVELS)
}

# ---- 1. concordance by recorded host class ---------------------------------
OUTCOMES <- c(
  "Match: species",
  "Match: genus",
  "Match: family",
  "Match: order / class",
  "Mismatch: species or genus",
  "Mismatch: family or above",
  "Lab DNA / lane cross-talk only",
  "No dominant host",
  "No blood meal detected"
)
outcome_lab <- function(x) {
  case_when(
    x == "MATCH_species" ~ OUTCOMES[1],
    x == "MATCH_genus" ~ OUTCOMES[2],
    x == "MATCH_family" ~ OUTCOMES[3],
    x %in% c("MATCH_order", "MATCH_class") ~ OUTCOMES[4],
    x %in% c("MISMATCH_species", "MISMATCH_genus") ~ OUTCOMES[5],
    str_starts(x, "MISMATCH") ~ OUTCOMES[6],
    x == "lab_or_crosstalk_only" ~ OUTCOMES[7],
    x == "no_dominant_host" ~ OUTCOMES[8],
    x == "no_bloodmeal_detected" ~ OUTCOMES[9],
    TRUE ~ NA_character_ # no host recorded / questing: nothing to compare
  ) |>
    factor(levels = OUTCOMES)
}
# match depth: one sequential hue (darkest = species); mismatch: Orange
# (darker = deeper conflict); lab/cross-talk only: Lavender; untestable: greys.
OUTCOME_PAL <- setNames(
  c(
    CSIRO[["MidnightBlue"]],
    CSIRO[["OceanBlue"]],
    CSIRO[["MiddayBlue"]],
    "#CCEEF5",
    "#F6B98A",
    CSIRO[["Orange"]],
    CSIRO[["Lavender"]],
    CSIRO[["Steel"]],
    CSIRO[["Mist"]]
  ),
  OUTCOMES
)

# Samples with no recorded host (or questing) have nothing to compare and are
# left out; they are in bloodmeal_calls.tsv and figure 3.
d1 <- calls |>
  mutate(
    meta_class = droplevels(class_lab(meta_class)),
    outcome = outcome_lab(comparison),
    evidence = factor(
      if_else(evidence == "weak", "weak (< min_reads)", "strong"),
      levels = c("strong", "weak (< min_reads)")
    )
  ) |>
  filter(!is.na(outcome)) |>
  count(meta_class, outcome, evidence)
n_cls <- d1 |> group_by(meta_class) |> summarise(n = sum(n))
p1 <- ggplot(
  d1,
  aes(y = fct_rev(meta_class), x = n, fill = outcome, alpha = evidence)
) +
  geom_col(colour = "white", linewidth = 0.5, width = 0.7) +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 2.8,
    colour = if_else(d1$outcome %in% OUTCOMES[c(1, 2, 6)], "white", "black")
  ) +
  geom_text(
    data = n_cls,
    aes(y = meta_class, x = n, label = paste0("n=", n)),
    inherit.aes = FALSE,
    hjust = -0.2,
    size = 3
  ) +
  scale_fill_manual(values = OUTCOME_PAL, drop = FALSE, name = NULL) +
  scale_alpha_manual(
    values = c(strong = 1, `weak (< min_reads)` = 0.45),
    name = "evidence",
    drop = FALSE
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  guides(
    fill = guide_legend(
      ncol = 3,
      order = 1,
      override.aes = list(colour = MARK_OUTLINE)
    ),
    alpha = guide_legend(ncol = 1, order = 2)
  ) +
  labs(
    x = "samples",
    y = "recorded host class",
    title = "Blood-meal host from sequence vs recorded host",
    subtitle = paste(
      "Deepest rank at which the primary sequence call agrees with the",
      "metadata host,\nor the first rank at which they differ",
      "(samples with no recorded host omitted)"
    )
  ) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())
save_plot(
  file.path(fig_dir, "bloodmeal_concordance.png"),
  p1,
  width = 9,
  height = 4.8,
  dpi = 200
)

# ---- 2. class matrix -------------------------------------------------------
d2 <- calls |>
  mutate(
    recorded = class_lab(meta_class),
    detected = case_when(
      !is.na(primary_class) ~ primary_class,
      tier_used == "none" ~ "not detected",
      TRUE ~ "no dominant host"
    )
  ) |>
  count(recorded, detected) |>
  mutate(
    detected = fct_relevel(
      fct_infreq(detected, n),
      "no dominant host",
      "not detected",
      after = Inf
    )
  )
p2 <- ggplot(d2, aes(x = detected, y = fct_rev(recorded), fill = n)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(
    aes(label = n, colour = n > max(n) / 2),
    size = 3.5,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_fill_gradientn(colours = CSIRO_SEQ[-1], name = "samples") +
  labs(
    x = "host class detected (primary call)",
    y = "recorded host class",
    title = "Recorded vs detected host class"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid = element_blank()
  )
save_plot(
  file.path(fig_dir, "bloodmeal_class_matrix.png"),
  p2,
  width = 7,
  height = 4,
  dpi = 200
)

# ---- 3. per sample ---------------------------------------------------------
used <- calls |> select(sample_id, tier_used)
d3 <- taxa |>
  semi_join(
    used |> filter(tier_used != "none"),
    by = c("sample_id", tier = "tier_used")
  ) |>
  # weak-evidence samples have no taxon over the floor: show what they have
  left_join(calls |> select(sample_id, evidence), by = "sample_id") |>
  filter(coalesce(passes_floor, "") == "yes" | evidence == "weak") |>
  mutate(
    family = coalesce(
      call_family,
      paste0(coalesce(call_class, "Vertebrata"), " (unresolved)")
    ),
    flag = case_when(
      coalesce(lab_contaminant, "") == "yes" ~ "lab",
      !is.na(cross_talk) ~ "xt",
      TRUE ~ "host"
    )
  ) |>
  group_by(sample_id, family, call_class) |>
  summarise(
    share = sum(share),
    # a family is only as clean as its least clean taxon
    flag = c("host", "xt", "lab")[max(match(flag, c("host", "xt", "lab")))],
    .groups = "drop"
  ) |>
  left_join(
    calls |>
      select(
        sample_id,
        meta_family,
        meta_class,
        host_species,
        tick_genus,
        comparison
      ),
    by = "sample_id"
  ) |>
  mutate(agrees = !is.na(meta_family) & family == meta_family)
# samples with nothing detected still get a row (their recorded host is ringed)
rows <- calls |>
  mutate(
    meta_class = class_lab(meta_class),
    label = paste0(
      sample_id,
      "  [",
      coalesce(host_species, "no host recorded"),
      "]"
    )
  ) |>
  arrange(meta_class, meta_family, sample_id)
fam_order <- d3 |>
  distinct(family, call_class) |>
  arrange(factor(call_class, levels = CLASS_LEVELS), family) |>
  pull(family) |>
  unique()
fam_levels <- unique(c(
  fam_order,
  sort(setdiff(na.omit(rows$meta_family), fam_order))
))
d3 <- d3 |>
  mutate(
    label = rows$label[match(sample_id, rows$sample_id)],
    meta_class = rows$meta_class[match(sample_id, rows$sample_id)]
  )
rings <- rows |> filter(!is.na(meta_family))
y_levels <- rev(rows$label)
p3 <- ggplot() +
  geom_point(
    data = rings,
    aes(
      x = factor(meta_family, levels = fam_levels),
      y = factor(label, levels = y_levels)
    ),
    shape = 0,
    size = 4.2,
    colour = CSIRO[["Steel"]],
    stroke = 0.7
  ) +
  geom_point(
    data = d3,
    aes(
      x = factor(family, levels = fam_levels),
      y = factor(label, levels = y_levels),
      size = share,
      fill = agrees,
      shape = flag
    ),
    colour = MARK_OUTLINE,
    stroke = 0.4
  ) +
  scale_size_area(
    max_size = 3.6,
    name = "share of host reads",
    breaks = c(.1, .5, 1)
  ) +
  scale_fill_manual(
    values = c(`TRUE` = CSIRO[["Mint"]], `FALSE` = CSIRO[["Orange"]]),
    labels = c(`TRUE` = "recorded host's family", `FALSE` = "other family"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(host = 21, lab = 24, xt = 25),
    labels = c(
      host = "host",
      lab = "likely lab DNA (human/chicken/cow/pig)",
      xt = "lane cross-talk?"
    ),
    name = NULL
  ) +
  scale_x_discrete(drop = FALSE) +
  facet_grid(meta_class ~ ., scales = "free_y", space = "free_y") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
  labs(
    x = "host family detected (evidence tier used for the call)",
    y = NULL,
    title = "Blood-meal hosts detected per sample",
    subtitle = paste(
      "Grey square = family of the recorded host; facets = recorded host class.",
      "\nSamples with no points: no vertebrate mito / 18S contig assembled"
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 7),
    axis.text.y = element_text(size = 5.5),
    legend.position = "top",
    legend.box = "vertical",
    strip.text.y = element_text(angle = 0)
  )
save_plot(
  file.path(fig_dir, "bloodmeal_per_sample.png"),
  p3,
  width = 10,
  height = max(7, 0.11 * nrow(rows) + 2.5),
  dpi = 200,
  limitsize = FALSE
)

say("figures in %s", fig_dir)
