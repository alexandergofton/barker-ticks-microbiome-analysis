#!/usr/bin/env Rscript
# Step 26b -- cross-check the curated tick identifications against the
# mitogenome phylogeny.
#
#   Rscript code/analysis/26b_mito_vs_metadata.R
#
# The mitogenome tree (steps 20-25) is independent molecular evidence of what
# each tick is. The tip labels baked into it come from whatever species name
# existed when step 21 ran (often the id-string guess). This step ignores those
# and re-labels every tip from the CURATED metadata/sample_metadata.tsv, then
# asks two questions:
#
#   0. Per sample: what is the selected contig by its own nt-BLAST LCA? A
#      blood-meal (vertebrate) or bacterial contig gets verdict NOT_TICK_CONTIG.
#   1. Per sample: are its nearest relatives in the tree the same tick genus /
#      species as the metadata says? (k nearest tips by patristic distance)
#   2. Per curated species with >= 2 tips: is it monophyletic in the tree?
#
# Published reference mitogenomes (step 20b, tips "REF_<accession>") take part
# as neighbours of known identity -- they are never judged themselves -- and
# each sample also gets its single nearest reference (nearest_ref*), the
# closest thing this tree offers to a named-species identification.
#
# A disagreement is NOT automatically a metadata error. The mitogenome
# candidate is chosen by homology from a whole-tick assembly (step 21), so a
# disagreement can also be: the wrong contig picked (a contaminating or
# index-hopped tick's mitogenome), a NUMT, or a sparse sample on a long branch.
# Read the `n_genes` and `nn_dist` columns before deciding which.
#
# OUTPUTS
#   results/mito/tables/mito_vs_metadata.tsv          one row per tip
#   results/mito/tables/mito_species_monophyly.tsv    one row per curated species
#   results/mito/figures/mito_tree_curated.png        tips coloured by curated
#                                                     genus; discordant tips flagged

suppressMessages({ library(ape); library(phangorn); library(tidyverse) })

.args <- commandArgs(trailingOnly = FALSE)
.this <- .args |> str_subset("^--file=") |> str_remove("^--file=")
source(file.path(dirname(normalizePath(.this)), "common.R"))
params <- load_params()
mito_dir <- file.path(ROOT, params$mito$out_dir)
K <- 3   # nearest neighbours consulted

refs <- mito_ref_table(params)
tree <- read.tree(file.path(mito_dir, "tree", "mito.treefile"))
tree <- root_mito_tree(tree, refs, params$mito$tree$outgroup_family)
say("tree %s", attr(tree, "rooting"))
sid  <- str_split_i(tree$tip.label, fixed("|"), 1)
is_ref <- str_starts(sid, fixed("REF_"))

# (tables are read/written with data.table's fread/fwrite: see common.R)
meta <- data.table::fread(file.path(ROOT, "metadata", "sample_metadata.tsv"),
                          na.strings = c("", "NA")) |>
  as_tibble()
sel  <- data.table::fread(file.path(mito_dir, "mito_selection.tsv")) |> as_tibble()
refsum <- data.table::fread(file.path(mito_dir, "ref_genomes", "ref_genome_summary.tsv")) |>
  as_tibble()
# one row per tree tip, in tip order (all-NA rows for tips not in the metadata,
# which the references then fill in)
m <- meta[match(sid, meta$sample_id), ] |>
  select(sample_id, tick_genus, tick_species, species_identified) |>
  mutate(n_genes = sel$n_genes_found[match(sid, sel$sample_id)]) |>
  mutate(sample_id          = if_else(is_ref, sid, sample_id),
         tick_genus         = if_else(is_ref, refs$genus[match(sid, refs$ref_id)], tick_genus),
         tick_species       = if_else(is_ref, refs$species[match(sid, refs$ref_id)], tick_species),
         species_identified = if_else(is_ref, "yes", species_identified),
         n_genes            = if_else(is_ref, refsum$n_pcg[match(sid, refsum$ref_id)], n_genes)) |>
  mutate(genus   = if_else(is.na(tick_genus), "unidentified", tick_genus),
         species = if_else(is.na(tick_species), "unidentified", tick_species),
         ident   = species_identified %in% "yes")

# ---- 1. nearest neighbours ---------------------------------------------------
D <- cophenetic(tree)
rownames(D) <- colnames(D) <- sid
res <- map(which(!is_ref), function(i) {
  d <- D[i, -i]; o <- order(d)[seq_len(min(K, length(d)))]
  nb <- names(d)[o]; nbm <- m[match(nb, sid), ]
  dr <- d[str_starts(names(d), fixed("REF_"))]
  r1 <- names(dr)[which.min(dr)]
  r1sp <- if (length(r1)) m$species[match(r1, sid)] else NA_character_
  # only neighbours with an identified genus/species can vote
  g_vote <- nbm$genus[nbm$genus != "unidentified"]
  s_vote <- nbm$species[nbm$ident]
  tibble(
    sample_id = sid[i], tick_genus = m$genus[i], tick_species = m$species[i],
    n_genes = m$n_genes[i],
    nn1 = nb[1], nn1_species = nbm$species[1], nn1_dist = round(d[o[1]], 4),
    knn_genera = paste(nbm$genus, collapse = "; "),
    knn_species = paste(nbm$species, collapse = "; "),
    nearest_ref = if (length(r1)) str_remove(r1, "^REF_") else NA_character_,
    nearest_ref_species = r1sp,
    nearest_ref_dist = if (length(r1)) round(min(dr), 4) else NA_real_,
    ref_genus_agrees = if (m$genus[i] == "unidentified" || !length(r1)) NA
                       else m$genus[i] == str_remove(r1sp, " .*"),
    genus_agrees = if (m$genus[i] == "unidentified" || !length(g_vote)) NA
                   else m$genus[i] %in% g_vote,
    species_agrees = if (!m$ident[i] || !length(s_vote)) NA
                     else m$species[i] %in% s_vote)
}) |>
  list_rbind() |>
  mutate(verdict = case_when(
    is.na(genus_agrees)       ~ "no_call",
    !genus_agrees             ~ "GENUS_DISCORDANT",
    species_agrees %in% FALSE ~ "species_discordant",
    .default                  = "concordant"))
# ---- 1b. what IS the selected contig, per the whole-assembly nt BLAST? -------
# Step 21 picks the contig with the most mitochondrial genes by tblastn against
# TICK probes, but a vertebrate (blood-meal host) mitogenome also carries all
# 13 genes and is hit at ~40-60% protein identity -- comfortably past the
# E-value gate. The contig's own nt-BLAST LCA (lca/, classified by step 03)
# says what it really is. scope_class "host" = tick; "bloodmeal" = vertebrate
# host; "microbial" = a bacterial contig with spurious weak hits.
lin <- data.table::fread(file.path(ROOT, "results", "taxonomy", "taxon_lineage.tsv"),
                         select = c("taxid", "lca_name", "scope_class", "genus")) |>
  as_tibble()
ali <- tryCatch(data.table::fread(file.path(ROOT, "metadata", "sample_aliases.tsv"), header = FALSE,
                                  col.names = c("raw", "canonical")) |> as_tibble(),
                error = function(e) tibble(raw = character(), canonical = character()))
contig_lca <- function(s, contig) {
  if (is.na(contig) || contig == "") return(c(NA, NA, NA))
  for (nm in unique(c(s, filter(ali, canonical == s)$raw))) {
    f <- file.path(ROOT, "lca", paste0(nm, ".blast.lca"))
    if (!file.exists(f)) next
    hit <- suppressWarnings(system2("grep", c("-m1", "-F", shQuote(paste0(contig, "\t")), shQuote(f)),
                                    stdout = TRUE))
    hit <- str_subset(hit, paste0("^", contig, "\t"))
    if (!length(hit)) return(c("no_blast_hit", NA, NA))
    k <- match(as.integer(str_split_i(hit[1], "\t", 2)), lin$taxid)
    return(c(lin$lca_name[k], lin$scope_class[k], lin$genus[k]))
  }
  c("no_lca_file", NA, NA)
}
lc <- map(res$sample_id, function(s)
  contig_lca(s, sel$selected_contig[match(s, sel$sample_id)]))
res <- res |>
  mutate(contig_lca = map_chr(lc, 1), contig_class = map_chr(lc, 2), contig_genus = map_chr(lc, 3))

# The k-nearest-neighbour vote has a blind spot: two samples carrying the SAME
# foreign tick's mitogenome sit together and vote for each other's curated
# genus (e.g. two "Amblyomma" tubes that both hold an Argas mitogenome). So a
# sample is also GENUS_DISCORDANT when two independent lines agree against the
# curated genus: the contig's own nt-BLAST genus, and the genus of the nearest
# published reference in the tree.
res <- res |>
  mutate(ref_genus = str_remove(nearest_ref_species, " .*"),
         two_line_discord = contig_class %in% "host" & !is.na(contig_genus) & contig_genus != "" &
           tick_genus != "unidentified" & contig_genus != tick_genus & !is.na(ref_genus) &
           ref_genus == contig_genus,
         verdict = if_else(two_line_discord %in% TRUE, "GENUS_DISCORDANT", verdict),
         verdict = if_else(contig_class %in% c("bloodmeal", "microbial"), "NOT_TICK_CONTIG", verdict)) |>
  arrange(verdict, tick_genus, tick_species)
data.table::fwrite(res, file.path(mito_dir, "tables", "mito_vs_metadata.tsv"), sep = "\t")

# ---- 2. species monophyly ---------------------------------------------------
# Tested over samples AND references of that name: a curated species that is
# monophyletic only because no reference of a sister species was in the tree
# is a weaker result than one that holds with the references present. Only
# taxa with at least one sample are tested -- this is about the curated data.
sp_tab <- m |>
  filter(ident) |>
  summarise(N = n(), n_samp = sum(!str_starts(sample_id, fixed("REF_"))), .by = species) |>
  filter(N >= 2, n_samp >= 1)
mono <- map(sp_tab$species, function(s) {
  w <- m$species == s & m$ident
  tibble(tick_species = s, n_tips = sum(w), n_refs = sum(w & is_ref),
         monophyletic = is.monophyletic(tree, tree$tip.label[w]),
         tips = paste(sid[w], collapse = ", "))
}) |>
  list_rbind()
gen_mono <- map(setdiff(unique(m$genus[!is_ref]), "unidentified"), function(g) {
  w <- m$genus == g
  if (sum(w) < 2) return(NULL)
  tibble(tick_species = paste0(g, " (genus)"), n_tips = sum(w), n_refs = sum(w & is_ref),
         monophyletic = is.monophyletic(tree, tree$tip.label[w]), tips = "")
}) |>
  list_rbind()
mono <- bind_rows(gen_mono, mono)
data.table::fwrite(mono, file.path(mito_dir, "tables", "mito_species_monophyly.tsv"), sep = "\t")

# ---- figure -----------------------------------------------------------------
flag <- res$verdict[match(sid, res$sample_id)]
col  <- unname(TICK_GENUS_PAL[if_else(m$genus %in% names(TICK_GENUS_PAL), m$genus, "unidentified")])
col  <- if_else(is_ref, "black", col)
mark <- case_when(flag %in% "GENUS_DISCORDANT"   ~ " <<< GENUS",
                  flag %in% "NOT_TICK_CONTIG"    ~ " ### NOT TICK",
                  flag %in% "species_discordant" ~ " < sp.",
                  .default = "")
tree$tip.label <- if_else(is_ref, sprintf("REF  %s  %s", m$species, str_remove(sid, "^REF_")),
                          sprintf("%s  [%s]%s", sid, m$species, mark))
h <- max(8, 0.13 * length(sid))
png_svg(file.path(mito_dir, "figures", "mito_tree_curated.png"), width = 11, height = h,
        res = 150, draw = function() {
par(mar = c(3, 1, 3, 1))
# Label text stays dark (several CSIRO hues are too light to read as text);
# the genus colour is a filled square at each tip, outlined for contrast.
off <- 0.012 * max(ape::node.depth.edgelength(tree))
plot(tree, cex = 0.55, tip.color = if_else(is_ref, "black", "grey20"),
     font = case_when(is_ref ~ 4,
                      flag %in% c("GENUS_DISCORDANT", "species_discordant", "NOT_TICK_CONTIG") ~ 2,
                      .default = 1),
     label.offset = off, no.margin = FALSE)
tiplabels(pch = 22, bg = col, col = MARK_OUTLINE, cex = 0.75, lwd = 0.4, offset = off / 2)
add.scale.bar(cex = 0.6)
gl <- intersect(names(TICK_GENUS_PAL), m$genus)
legend("bottomleft", legend = c(gl, "REF = published mitogenome"),
       pch = 22, pt.bg = c(TICK_GENUS_PAL[gl], "black"), col = MARK_OUTLINE, pt.cex = 1.2,
       bty = "n", cex = 0.7, title = "curated tick genus")
title(main = sprintf("Mitogenome ML tree, tips labelled from the CURATED metadata (%s)", attr(tree, "rooting")),
      sub = sprintf("bold '<<< GENUS': none of the %d nearest tips share the curated genus; '< sp.': genus agrees, species does not; '### NOT TICK': selected contig is a blood-meal or bacterial sequence by nt BLAST", K),
      cex.main = 0.9, cex.sub = 0.7)
})

say("tips: %d", nrow(res))
res |> summarise(N = n(), .by = verdict) |> print(n = Inf)
say("discordant samples:")
res |>
  filter(verdict %in% c("GENUS_DISCORDANT", "species_discordant", "NOT_TICK_CONTIG")) |>
  select(sample_id, tick_species, n_genes, nn1, nn1_species, nn1_dist, nearest_ref_species) |>
  print(n = Inf, width = Inf)
say("wrote results/mito/tables/mito_vs_metadata.tsv, mito_species_monophyly.tsv,")
say("      results/mito/figures/mito_tree_curated.{png,svg}")
