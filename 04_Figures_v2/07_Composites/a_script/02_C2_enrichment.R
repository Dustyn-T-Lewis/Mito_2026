#!/usr/bin/env Rscript
# C2 enrichment — the four contrast rings (Disease / Transplant / Rescue /
# Interaction) as a 2x2, drawn by enrichVolcano::volcano_ring_grid() from the
# reported DE table (pi_score for the volcano) and the fgsea cache (NES arcs,
# deduped across the 5-DB lens). Writes MAIN_C2_enrichment.png.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(enrichVolcano)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))

BASE <- here::here("04_Figures_v2", "07_Composites")
OUT <- file.path(BASE, "b_reports", "main", "png")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

dep_df <- load_combined_wide()
fgsea_all <- readr::read_csv(
  here::here("04_Figures", "shared", "fgsea_tstat_all_h9c2.csv"),
  show_col_types = FALSE
)
rat_gene_sets <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))
set_pool <- do.call(c, unname(rat_gene_sets[CANONICAL_DBS]))

# old contrast key -> role used as the ring title; pi_score drives volcano
# significance (suite-wide Π < 0.05), fgsea padj drives ring selection.
ring_map <- c(
  CTLvPHE = "Disease", CTLvMITO = "Transplant",
  PHEvPHE_MITO = "Rescue", Interaction = "Interaction"
)
RING_N <- 12

volc_long <- bind_rows(lapply(names(ring_map), function(ctr) {
  dep_df |>
    transmute(
      gene,
      logFC = .data[[paste0("logFC_", ctr)]],
      P.Value = .data[[paste0("P.Value_", ctr)]],
      padj = .data[[paste0("pi_score_", ctr)]],
      contrast = ring_map[[ctr]]
    ) |>
    filter(!is.na(logFC), !is.na(P.Value))
}))

# Per contrast: keep significant terms, collapse redundant biology across DBs
# (EnrichmentMap combined coefficient), then take the top RING_N by FDR.
enrich_long <- bind_rows(lapply(names(ring_map), function(ctr) {
  sig <- fgsea_all |>
    filter(contrast == ctr, database %in% CANONICAL_DBS, padj < 0.05) |>
    as.data.frame()
  deduplicate_enrichment(sig, set_pool, cutoff = 0.375) |>
    arrange(padj) |>
    slice_head(n = RING_N) |>
    mutate(contrast = ring_map[[ctr]])
}))

grid <- volcano_ring_grid(
  volc_dfs = volc_long,
  enrich_dfs = enrich_long,
  contrasts = unname(ring_map),
  subtitles = unname(CONTRAST_MATH_BRIEF[names(ring_map)]),
  ncol = 2,
  gene_col = "gene", logfc_col = "logFC", pval_col = "P.Value", padj_col = "padj",
  term_col = "pathway", nes_col = "NES", size_col = "size",
  genes_col = "leadingEdge", genes_sep = ";",
  p_threshold = 0.05, logfc_threshold = log2(1.5),
  volcano_radius = 4.2,
  label_size = 2.7, point_size = 2.2, point_alpha = 0.7,
  count_x_mult = 0.55, count_y_mult = 0.55,
  theme = volcano_ring_theme(base_size = 13)
)

# Pull the four rings together: trim per-panel margins and move the shared
# NES bar to the bottom as a horizontal strip so the panels fill the canvas
# instead of floating in white space.
fig <- grid$plot &
  ggplot2::theme(
    legend.position = "bottom",
    legend.key.width = ggplot2::unit(18, "mm"),
    legend.key.height = ggplot2::unit(3, "mm"),
    plot.margin = ggplot2::margin(1, 1, 1, 1, "mm")
  )
fig <- fig & ggplot2::guides(
  fill = ggplot2::guide_colorbar(direction = "horizontal", title.position = "top")
)

ggplot2::ggsave(
  file.path(OUT, "MAIN_C2_enrichment.png"), fig,
  width = 220, height = 230, units = "mm", dpi = 300, limitsize = FALSE
)
message("C2 enrichment composite built")
