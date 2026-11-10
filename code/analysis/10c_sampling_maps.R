#!/usr/bin/env Rscript
# Step 10c -- maps of where the ticks were collected.
#
#   conda activate tick_micro
#   Rscript code/analysis/10c_sampling_maps.R
#
# OUTPUTS
#   results/figures/16a_sampling_map.png              all samples, colour = tick genus,
#                                                     shape = host type, size = n at the site
#   results/figures/16b_sampling_map_by_genus.png     one panel per tick genus
#   results/figures/16c_endosymbiont_map.png          where each endosymbiont was detected
#   results/maps/sampling_map.html                    interactive (Leaflet): toggle genera,
#                                                     click a tick for its metadata and
#                                                     microbial summary
#   results/maps/sampling_sites.tsv                   the mapped table
#
# Every sample with coordinates is mapped, INCLUDING samples that never reached
# the abundance matrix (they were still collected); the popup and the
# `analysed` column say which. Samples outside the static extent (the Falklands
# I. uriae) are named in the caption; the interactive map shows them all.
#
# Base layers: Natural Earth countries (rnaturalearthdata, offline). State
# borders are downloaded from Natural Earth once and cached under
# results/maps/cache/; if that fails the maps are drawn without them.

# Runs in the tick_report env, which has these tidyverse packages but not the
# tidyverse meta-package (nor readr/forcats), so they are attached one by one.
suppressPackageStartupMessages({
  library(yaml); library(ggplot2)
  library(sf); library(rnaturalearth)
  library(dplyr); library(tidyr); library(purrr); library(stringr); library(tibble)
})

.here <- commandArgs(FALSE) |>
  str_subset("^--file=") |>
  first() |>
  str_remove("^--file=") |>
  dirname()
source(file.path(.here, "common.R"))
params <- load_params()
mp <- params$maps
# theme: theme_tick() from common.R (theme_classic + light-grey major gridlines)
sf_use_s2(FALSE)

mdir <- file.path(ROOT, "results", "maps"); dir.create(mdir, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Data: every sample with coordinates + what the analysis found in it
# ---------------------------------------------------------------------------
# (tables are read/written with data.table's fread/fwrite: see common.R)
meta <- data.table::fread(file.path(ROOT, "metadata", "sample_metadata.tsv"),
                          na.strings = c("", "NA")) |>
  as_tibble()
mat_ids <- names(data.table::fread(file.path(ROOT, "results", "matrices", "tpm_microbial_genus.tsv"),
                                   nrows = 0))[-1]
meta <- meta |>
  mutate(analysed = sample_id %in% mat_ids,
         genus    = tick_genus_factor(tick_genus),
         host     = host_type_factor(host_type),
         species  = if_else(is.na(tick_species), "unidentified", tick_species))

tp <- data.table::fread(file.path(ROOT, "results", "targets", "target_prevalence.tsv")) |>
  as_tibble()
cats <- params$tick_groups$symbiont_categories
endo <- tp |>
  filter(present == 1, category %in% cats) |>
  summarise(endosymbionts = paste(sprintf("%s (%.1f%%)", target, tpm_microbial / 1e4)[order(-tpm_microbial)],
                                  collapse = ", "),
            .by = sample_id)
path_ <- tp |>
  filter(present == 1, category %in% c("pathogen", "uncertain")) |>
  summarise(other_targets = paste(target[order(-tpm_microbial)], collapse = ", "),
            .by = sample_id)
meta <- meta |>
  left_join(endo, by = "sample_id") |>
  left_join(path_, by = "sample_id")
domf <- file.path(ROOT, "results", "tables", "tick_groups", "dominant_genus_per_sample.tsv")
if (file.exists(domf)) {
  dom <- data.table::fread(domf) |>
    as_tibble() |>
    select(sample_id, dominant, dominant_pct)
  meta <- meta |> left_join(dom, by = "sample_id")
} else meta <- meta |> mutate(dominant = NA_character_, dominant_pct = NA_real_)
meta <- meta |> relocate(sample_id) |> arrange(sample_id)

pts <- meta |> filter(!is.na(latitude), !is.na(longitude))
data.table::fwrite(
  pts |> select(sample_id, status, analysed, tick_family, tick_genus, tick_species,
                life_stage, n_ticks_pooled, host_type, host_species, collection_site,
                latitude, longitude, endosymbionts, other_targets, dominant, dominant_pct),
  file.path(mdir, "sampling_sites.tsv"), sep = "\t")
say("%d of %d samples have coordinates (%d analysed)", nrow(pts), nrow(meta), sum(pts$analysed))
no_xy <- meta |> filter(is.na(latitude))
say("not mapped (no coordinates): %d -- %s", nrow(no_xy), paste(no_xy$sample_id, collapse = ", "))

# ---------------------------------------------------------------------------
# Base layers
# ---------------------------------------------------------------------------
land <- ne_countries(scale = 50, returnclass = "sf")
cache <- file.path(mdir, "cache"); dir.create(cache, showWarnings = FALSE)
states <- tryCatch({
  f <- list.files(cache, pattern = "admin_1_states_provinces_lines\\.(gpkg|shp)$", full.names = TRUE)
  if (!length(f)) ne_download(scale = 50, type = "admin_1_states_provinces_lines",
                              category = "cultural", destdir = cache, load = FALSE)
  st_read(list.files(cache, pattern = "admin_1_states_provinces_lines\\.(gpkg|shp)$",
                     full.names = TRUE)[1], quiet = TRUE)
}, error = function(e) { say("  state borders unavailable (%s) -- drawn without", conditionMessage(e)); NULL })
if (!is.null(states) && "ADM0_NAME" %in% names(states))
  states <- states |> filter(ADM0_NAME %in% c("Australia", "Papua New Guinea"))

XL <- unlist(mp$xlim); YL <- unlist(mp$ylim)
inside <- pts |>
  filter(longitude >= XL[1], longitude <= XL[2], latitude >= YL[1], latitude <= YL[2])
outside <- pts |> filter(!sample_id %in% inside$sample_id)
off_note <- if (nrow(outside)) sprintf("Outside the map extent: %s.",
  paste(sprintf("%s (%s, %s; %.2f, %.2f)", outside$sample_id, outside$species,
                outside$collection_site, outside$latitude, outside$longitude), collapse = "; ")) else ""
nomap_note <- sprintf("%d samples have no coordinates and are not shown.", nrow(no_xy))

basemap <- function() {
  g <- list(geom_sf(data = land, fill = "#efeee9", colour = "#b9b8b2", linewidth = 0.25))
  if (!is.null(states)) g <- c(g, list(geom_sf(data = states, colour = "#c9c8c2",
                                               linewidth = 0.2, linetype = 2)))
  c(g, list(coord_sf(xlim = XL, ylim = YL, expand = FALSE),
            theme(panel.background = element_rect(fill = "#f4f8fc"),
                  panel.grid.major = element_line(colour = "grey85", linewidth = 0.3)),
            labs(x = NULL, y = NULL)))
}
has_spatial <- requireNamespace("ggspatial", quietly = TRUE)
decor <- function() if (has_spatial) list(
  ggspatial::annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.7),
  ggspatial::annotation_north_arrow(location = "tr", height = unit(0.9, "cm"),
                                    width = unit(0.7, "cm"),
                                    style = ggspatial::north_arrow_minimal())) else list()

# ---------------------------------------------------------------------------
# 16a -- all samples. Sites with several ticks of the same genus + host become
# one larger point rather than a hidden stack.
# ---------------------------------------------------------------------------
site <- inside |> summarise(n = n(), .by = c(latitude, longitude, genus, host))
p <- ggplot() + basemap() +
  geom_point(data = arrange(site, desc(n)), aes(longitude, latitude, fill = genus, shape = host, size = n),
             colour = MARK_OUTLINE, stroke = 0.3, alpha = 0.9) +
  scale_fill_manual(values = TICK_GENUS_PAL, name = "tick genus", drop = TRUE,
                    guide = guide_legend(override.aes = list(shape = 21, size = 3.5))) +
  scale_shape_manual(values = HOST_TYPE_SHAPE, name = "host type", drop = TRUE,
                     guide = guide_legend(override.aes = list(fill = CSIRO_MUTED, size = 3))) +
  scale_size_area(max_size = 6, breaks = c(1, 2, 4, 6), name = "samples at site",
                  guide = guide_legend(override.aes = list(shape = 21, fill = CSIRO_MUTED))) +
  decor() +
  labs(title = sprintf("Tick sampling locations (%d samples with coordinates)", nrow(pts)),
       subtitle = "colour = tick genus; shape = host type; size = samples sharing a site, genus and host type",
       caption = paste(strwrap(paste(off_note, nomap_note), 130), collapse = "\n"))
# Inset: south-east Queensland, where a third of the sites sit within ~150 km
# and overlap at continental scale. Same encodings; the inset's extent is drawn
# as a box on the main map.
ZX <- unlist(mp$zoom_xlim); ZY <- unlist(mp$zoom_ylim)
zsite <- site |>
  filter(longitude >= ZX[1], longitude <= ZX[2], latitude >= ZY[1], latitude <= ZY[2])
pz <- ggplot() +
  geom_sf(data = land, fill = "#efeee9", colour = "#b9b8b2", linewidth = 0.25) +
  { if (!is.null(states)) geom_sf(data = states, colour = "#c9c8c2", linewidth = 0.2, linetype = 2) } +
  geom_point(data = arrange(zsite, desc(n)), aes(longitude, latitude, fill = genus, shape = host, size = n),
             colour = MARK_OUTLINE, stroke = 0.3, alpha = 0.9, show.legend = FALSE) +
  scale_fill_manual(values = TICK_GENUS_PAL) + scale_shape_manual(values = HOST_TYPE_SHAPE) +
  scale_size_area(max_size = 6, limits = c(1, max(site$n))) +
  coord_sf(xlim = ZX, ylim = ZY, expand = FALSE) +
  labs(x = NULL, y = NULL, title = sprintf("SE Queensland (%d sites)", nrow(distinct(zsite, latitude, longitude)))) +
  theme(panel.background = element_rect(fill = "#f4f8fc"), axis.text = element_text(size = 6),
        plot.title = element_text(size = 8), plot.background = element_rect(colour = "grey40"))
p <- p + annotate("rect", xmin = ZX[1], xmax = ZX[2], ymin = ZY[1], ymax = ZY[2],
                  fill = NA, colour = "grey30", linewidth = 0.4) +
  patchwork::inset_element(pz, left = 0.02, bottom = 0.03, right = 0.36, top = 0.42)
save_plot(fig_path("16a_sampling_map.png"), p, width = 10, height = 9, dpi = 200)
say("  wrote 16a_sampling_map.png")

# ---------------------------------------------------------------------------
# 16b -- one panel per genus; every other site in grey for context
# ---------------------------------------------------------------------------
bg <- inside |> distinct(latitude, longitude)
gn <- inside |> summarise(N = n(), .by = genus)
inside <- inside |>
  mutate(glab = factor(sprintf("%s (n=%d)", genus, gn$N[match(genus, gn$genus)]),
                       levels = sprintf("%s (n=%d)", levels(droplevels(gn$genus)),
                                        gn$N[match(levels(droplevels(gn$genus)), gn$genus)])))
p <- ggplot() + basemap() +
  geom_point(data = bg, aes(longitude, latitude), colour = "#bdbcb6", size = 0.8) +
  geom_point(data = inside, aes(longitude, latitude, fill = genus, shape = host),
             colour = MARK_OUTLINE, stroke = 0.3, size = 2.4) +
  scale_fill_manual(values = TICK_GENUS_PAL, guide = "none") +
  scale_shape_manual(values = HOST_TYPE_SHAPE, name = "host type", drop = TRUE,
                     guide = guide_legend(override.aes = list(fill = CSIRO_MUTED, size = 3))) +
  facet_wrap(~ glab, ncol = 3) +
  theme(strip.text = element_text(face = "italic"), axis.text = element_text(size = 6),
        legend.position = "bottom") +
  labs(title = "Sampling locations by tick genus",
       subtitle = "grey: all other sampled sites",
       caption = paste(strwrap(paste(off_note, nomap_note), 150), collapse = "\n"))
save_plot(fig_path("16b_sampling_map_by_genus.png"), p, width = 12, height = 12, dpi = 200)
say("  wrote 16b_sampling_map_by_genus.png")

# ---------------------------------------------------------------------------
# 16c -- endosymbiont detections (analysed samples only)
# ---------------------------------------------------------------------------
an <- inside |> filter(analysed)
et <- tp |>
  filter(present == 1, category %in% cats, sample_id %in% an$sample_id) |>
  summarise(N = n(), .by = target) |>
  arrange(desc(N))
et <- et |>
  filter(N >= params$tick_groups$min_detections_for_abundance) |>
  slice_head(n = 6)
if (nrow(et)) {
  both <- tp |> filter(category == "both") |> pull(target) |> unique()
  em <- map(et$target, function(t) {
    pos <- tp |> filter(target == t, present == 1) |> pull(sample_id)
    an |> transmute(longitude, latitude, genus, target = t, detected = sample_id %in% pos)
  }) |>
    list_rbind()
  em <- em |>
    mutate(tlab = factor(sprintf("%s%s (%d/%d)", target, if_else(target %in% both, "*", ""),
                                 et$N[match(target, et$target)], nrow(an)),
                         levels = sprintf("%s%s (%d/%d)", et$target,
                                          if_else(et$target %in% both, "*", ""), et$N, nrow(an))))
  p <- ggplot() + basemap() +
    geom_point(data = filter(em, !detected), aes(longitude, latitude),
               shape = 21, fill = NA, colour = CSIRO_MUTED, size = 1.4, stroke = 0.4) +
    geom_point(data = filter(em, detected), aes(longitude, latitude, fill = genus),
               shape = 21, colour = MARK_OUTLINE, size = 2.6, stroke = 0.3) +
    scale_fill_manual(values = TICK_GENUS_PAL, name = "tick genus (detected)", drop = TRUE,
                      guide = guide_legend(override.aes = list(size = 3.5))) +
    facet_wrap(~ tlab, ncol = 3) +
    theme(strip.text = element_text(face = "italic"), axis.text = element_text(size = 6),
          legend.position = "bottom") +
    labs(title = "Where each endosymbiont was detected",
         subtitle = "filled: presence call (step 10 thresholds); open grey: analysed, not detected",
         caption = paste(strwrap(paste(
           "* genus contains pathogens as well as endosymbionts; a genus-level call does not say which.",
           "Mapped samples only; denominators count analysed samples inside the map extent."), 150),
           collapse = "\n"))
  save_plot(fig_path("16c_endosymbiont_map.png"), p, width = 12, height = 9, dpi = 200)
  say("  wrote 16c_endosymbiont_map.png")
}

# ---------------------------------------------------------------------------
# Interactive map
# ---------------------------------------------------------------------------
if (!requireNamespace("leaflet", quietly = TRUE) || !requireNamespace("htmlwidgets", quietly = TRUE)) {
  say("leaflet/htmlwidgets not installed -- interactive map skipped"); quit(status = 0)
}
library(leaflet)
esc <- function(x) htmltools::htmlEscape(if_else(is.na(x) | x == "", "NA", as.character(x)))

# fan out samples that share exact coordinates, so each can be clicked
pts <- pts |>
  mutate(k = row_number(), nk = n(), .by = c(latitude, longitude)) |>
  mutate(ang  = 2 * pi * (k - 1) / nk,
         plat = latitude + if_else(nk > 1, mp$spread_deg * sin(ang), 0),
         plon = longitude + if_else(nk > 1, mp$spread_deg * cos(ang), 0))

pts <- pts |> mutate(popup = sprintf(paste0(
  "<div style='font:12px/1.45 system-ui,sans-serif;min-width:240px'>",
  "<b style='font-size:13px'>%s</b><br><i>%s</i> &middot; %s &middot; %s tick(s)<br>",
  "<span style='color:#52514e'>host</span> %s: %s<br>",
  "<span style='color:#52514e'>site</span> %s<br>",
  "<span style='color:#52514e'>coords</span> %.4f, %.4f%s<hr style='margin:4px 0'>",
  "%s</div>"),
  esc(sample_id), esc(species), esc(life_stage), esc(n_ticks_pooled),
  esc(host_type), esc(host_species), esc(collection_site), latitude, longitude,
  if_else(nk > 1, sprintf(" <span style='color:#52514e'>(%d samples here)</span>", nk), ""),
  if_else(analysed,
          sprintf(paste0("<span style='color:#52514e'>dominant microbial genus</span> %s<br>",
                         "<span style='color:#52514e'>endosymbionts</span> %s<br>",
                         "<span style='color:#52514e'>other targets</span> %s"),
                  if_else(is.na(dominant), "NA", sprintf("<i>%s</i> (%.0f%%)", esc(dominant), dominant_pct)),
                  esc(if_else(is.na(endosymbionts), "none called", endosymbionts)),
                  esc(if_else(is.na(other_targets), "none called", other_targets))),
          sprintf("<span style='color:#a33'>not analysed</span> (status: %s)", esc(status)))))

m <- leaflet(options = leafletOptions(preferCanvas = FALSE)) |>
  addProviderTiles(providers$CartoDB.Positron, group = "Light") |>
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") |>
  addProviderTiles(providers$OpenTopoMap, group = "Topographic")
groups <- character()
for (g in levels(droplevels(pts$genus))) {
  d <- pts |> filter(genus == g)
  grp <- sprintf("%s (n=%d)", g, nrow(d)); groups <- c(groups, grp)
  m <- m |> addCircleMarkers(
    lng = d$plon, lat = d$plat, group = grp, radius = 6,
    fillColor = TICK_GENUS_PAL[[g]], fillOpacity = if_else(d$analysed, 0.9, 0.35),
    color = MARK_OUTLINE, weight = 1, opacity = 1,
    popup = d$popup, label = sprintf("%s - %s", d$sample_id, d$species))
}
m <- m |>
  addLegend("bottomright", colors = unname(TICK_GENUS_PAL[levels(droplevels(pts$genus))]),
            labels = groups, title = "Tick genus", opacity = 0.9) |>
  addLayersControl(baseGroups = c("Light", "Satellite", "Topographic"),
                   overlayGroups = groups, options = layersControlOptions(collapsed = FALSE)) |>
  addScaleBar("bottomleft") |>
  fitBounds(XL[1], YL[1], XL[2], YL[2]) |>
  addControl(sprintf(paste0("<div style='font:12px system-ui,sans-serif;max-width:260px'>",
                            "<b>Barker tick microbiome: sampling sites</b><br>",
                            "%d samples with coordinates; %d analysed. Faded points were ",
                            "collected but not analysed. Ticks sharing a site are fanned ",
                            "out slightly; the popup gives true coordinates. %s</div>"),
                     nrow(pts), sum(pts$analysed), esc(nomap_note)), position = "topright")

# the widget object itself, so the Quarto report can embed the same map
saveRDS(m, file.path(mdir, "sampling_map_widget.rds"))

# saveWidget(selfcontained = TRUE) needs pandoc; the env ships one next to R
envbin <- normalizePath(file.path(R.home(), "..", "..", "bin"), mustWork = FALSE)
Sys.setenv(PATH = paste(envbin, Sys.getenv("PATH"), sep = .Platform$path.sep))
out <- file.path(mdir, "sampling_map.html")
ok <- tryCatch({
  htmlwidgets::saveWidget(m, out, selfcontained = TRUE, title = "Tick sampling map"); TRUE
}, error = function(e) { say("  self-contained save failed (%s); writing with lib dir", conditionMessage(e)); FALSE })
if (!ok) htmlwidgets::saveWidget(m, out, selfcontained = FALSE, libdir = "lib",
                                 title = "Tick sampling map")
unlink(file.path(mdir, "sampling_map_files"), recursive = TRUE)
say("  wrote results/maps/sampling_map.html  (open in a browser; tiles need internet)")
say("\ndone")
