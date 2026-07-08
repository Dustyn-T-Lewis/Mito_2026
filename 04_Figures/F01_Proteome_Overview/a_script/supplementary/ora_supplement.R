#!/usr/bin/env Rscript
# Supplementary ORA for F01: over-representation of the significant (Pi<0.05) proteins in
# each Venn contrast, split by direction, against the 5-DB lens (per-DB BH, EnrichmentMap
# dedup). Per contrast the top up- and down-regulated pathways sit back to back, named
# inside the bars with FDR at each tip; the overlap genes behind every shown pathway go to
# the supplementary workbook.

library(here)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
fns <- here::here("04_Figures", "functions")
source(file.path(fns, "F01-F03_data_paths_and_loaders.R"))
source(file.path(fns, "F01-F03_pathway_enrichment_dedup_ora.R"))
source(file.path(fns, "F01-F03_supplementary_workbook.R"))

BASE <- here::here("04_Figures", "F01_Proteome_Overview")
PNG <- file.path(BASE, "b_reports", "supplementary")
DAT <- file.path(BASE, "c_data")
dir.create(PNG, recursive = TRUE, showWarnings = FALSE)

# the three Venn contrasts, named via the shared map
CONTRASTS <- CONTRAST_DISPLAY_MAP[c("CTLvPHE", "CTLvMITO", "PHEvPHE_MITO")]
TOP_N <- 15L

de <- load_dep_long() |> filter(is_real_symbol(gene))
universe <- unique(de$gene)
collection <- build_harmonized_collection()
pw_by_db <- split(names(collection), classify_database(names(collection)))

# per-DB fora (each BH-corrected on its own) for one gene set
ora_set <- function(genes) {
  g <- intersect(genes, universe)
  bind_rows(lapply(names(pw_by_db), function(db) {
    dp <- collection[pw_by_db[[db]]]
    r <- as.data.frame(fgsea::fora(dp, genes = g, universe = universe, minSize = 10, maxSize = 500))
    r$leadingEdge <- vapply(r$overlapGenes, paste, character(1), collapse = ";")
    transmute(r, pathway, padj, n_overlap = overlap, set_size = size, database = db, leading_edge = leadingEdge)
  }))
}

ora <- bind_rows(lapply(names(CONTRASTS), function(ctr) {
  d <- filter(de, contrast == ctr, pi_score < H9C2_PI_THRESH, !is.na(logFC))
  bind_rows(
    up = ora_set(d$gene[d$logFC > 0]),
    down = ora_set(d$gene[d$logFC < 0]),
    .id = "direction"
  ) |>
    mutate(contrast = unname(CONTRASTS[ctr]))
})) |>
  filter(!pathway %in% MITO_DROP_SETS, !grepl(DISEASE_VIRAL_RE, pathway, ignore.case = TRUE))

# top distinct pathways one contrast-direction at a time, EnrichmentMap dedup then by FDR
top_pathways <- function(g) {
  deduplicate_enrichment(filter(g, padj < 0.10), collection) |>
    slice_min(padj, n = TOP_N, with_ties = FALSE) |>
    arrange(padj)
}

# every contrast in one figure: bars signed by direction (down left, up right), faceted
# by contrast, names inside the bars, FDR at each tip, coloured by database
shown <- bind_rows(lapply(unname(CONTRASTS), function(cn) {
  g <- filter(ora, contrast == cn)
  bind_rows(
    mutate(top_pathways(filter(g, direction == "down")), dir = -1),
    mutate(top_pathways(filter(g, direction == "up")), dir = 1)
  ) |>
    mutate(contrast = cn)
})) |>
  mutate(
    contrast = factor(contrast, levels = unname(CONTRASTS)),
    database = factor(database, levels = CANONICAL_DBS),
    lp = -log10(padj) * dir,
    name = clean_display_label(pathway),
    name = ifelse(nchar(name) > 32, paste0(substr(name, 1, 30), "…"), name),
    fdr = ifelse(padj < 0.001, formatC(padj, format = "e", digits = 0), formatC(padj, format = "f", digits = 3))
  ) |>
  arrange(contrast, lp) |>
  mutate(rid = factor(paste(contrast, database, pathway), levels = paste(contrast, database, pathway)))

missing <- setdiff(unname(CONTRASTS), as.character(unique(shown$contrast)))
sub <- "left = down-regulated, right = up-regulated; bar = -log10 FDR; name inside the bar; colour = database"
if (length(missing)) sub <- paste0(sub, " (no enrichment: ", paste(missing, collapse = ", "), ")")
xmax <- max(abs(shown$lp)) * 1.15
# name inside the bar (white) when it fits; outside (dark) when the bar is
# too short, with the FDR moving inside so the two never collide at the tip
shown <- shown |> mutate(fits = nchar(name) <= abs(lp) / xmax * 70)

fig <- ggplot(shown, aes(lp, rid, fill = database)) +
  geom_col(width = 0.78, colour = "black", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  shadowtext::geom_shadowtext(
    data = filter(shown, fits),
    aes(label = name, x = dir * 0.04, hjust = ifelse(dir > 0, 0, 1)),
    size = 2.4, fontface = "bold", colour = "white", bg.colour = "grey25", bg.r = 0.08
  ) +
  geom_text(
    data = filter(shown, !fits),
    aes(label = name, x = lp + dir * 0.04, hjust = ifelse(dir > 0, 0, 1)),
    size = 2.4, fontface = "bold", colour = "grey20"
  ) +
  geom_text(
    data = filter(shown, fits),
    aes(label = fdr, hjust = ifelse(dir > 0, -0.18, 1.18)),
    size = 1.9, fontface = "bold", colour = "grey20"
  ) +
  shadowtext::geom_shadowtext(
    data = filter(shown, !fits),
    aes(label = fdr, hjust = ifelse(dir > 0, 1.15, -0.15)),
    size = 1.9, fontface = "bold", colour = "white", bg.colour = "grey25", bg.r = 0.08
  ) +
  facet_grid(contrast ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = ORA_DB_COLORS, name = NULL, limits = CANONICAL_DBS) +
  scale_y_discrete(labels = NULL) +
  scale_x_continuous(limits = c(-xmax, xmax), expand = expansion(mult = 0.02)) +
  coord_cartesian(clip = "off") +
  labs(title = "ORA of the significant proteins per contrast", subtitle = sub, x = expression(-log[10] ~ FDR), y = NULL) +
  FIG_THEME +
  theme(
    plot.subtitle = element_text(size = FIG_AXIS_TEXT, colour = "grey35"),
    strip.text.y = element_text(angle = 0, face = "bold", size = FIG_AXIS_TEXT + 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(), legend.position = "bottom"
  )
ggsave(file.path(PNG, "SUPP_F01_ora.png"), fig, width = PANEL_MD, height = 130, units = "mm", dpi = 300)

# overlap behind every shown pathway -> supplementary workbook
overlap <- bind_rows(lapply(unname(CONTRASTS), function(cn) {
  g <- mutate(ora, database = factor(database, levels = CANONICAL_DBS)) |> filter(contrast == cn)
  bind_rows(lapply(c("up", "down"), function(dir) {
    mutate(top_pathways(filter(g, direction == dir)), contrast = cn, direction = dir)
  }))
})) |>
  transmute(contrast, direction, database,
    pathway = clean_display_label(pathway),
    FDR = signif(padj, 3), n_overlap, set_size,
    overlap_coef = round(n_overlap / set_size, 3), overlap_genes = leading_edge
  )

append_workbook(
  file.path(DAT, "F01_supplementary.xlsx"),
  sheet_specs = list(list(
    name = "ora_overlap", df = overlap,
    role = "Supplement: ORA of the significant proteins per contrast",
    contents = "shown pathways with FDR, overlap count, set size, overlap coefficient, and the overlapping genes"
  ))
)

mirror_to_box(file.path(PNG, "SUPP_F01_ora.png"), "02_Figures/F01_Proteome_Overview/supp")
mirror_to_box(file.path(DAT, "F01_supplementary.xlsx"), "03_Supplementary")

message(sprintf("F01 ORA supplement: 1 figure, %d shown pathways", nrow(overlap)))
