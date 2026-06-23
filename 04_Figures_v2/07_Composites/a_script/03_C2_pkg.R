#!/usr/bin/env Rscript
# C2 (package variant) — the four enrichment rings (Disease / Transplant /
# Rescue / Interaction) drawn by enrichVolcano::volcano_ring_grid() instead of
# the in-suite builder, fed the same tidy DE + fgsea tables. Parallel to
# 02_C2_enrichment.R for comparison; writes MAIN_C2_pkg_enrichment.png.

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

enrich_long <- fgsea_all |>
  filter(contrast %in% names(ring_map), database %in% CANONICAL_DBS, padj < 0.05) |>
  mutate(contrast = ring_map[contrast]) |>
  group_by(contrast) |>
  arrange(padj, .by_group = TRUE) |>
  slice_head(n = RING_N) |>
  ungroup()

grid <- volcano_ring_grid(
  volc_dfs = volc_long,
  enrich_dfs = enrich_long,
  contrasts = unname(ring_map),
  ncol = 2,
  gene_col = "gene", logfc_col = "logFC", pval_col = "P.Value", padj_col = "padj",
  term_col = "pathway", nes_col = "NES", size_col = "size",
  genes_col = "leadingEdge", genes_sep = ";",
  p_threshold = 0.05, logfc_threshold = log2(1.5)
)

ggplot2::ggsave(
  file.path(OUT, "MAIN_C2_pkg_enrichment.png"), grid$plot,
  width = PANEL_MD, height = PANEL_MD, units = "mm", dpi = 300, limitsize = FALSE
)
message("C2 (package) composite built")
