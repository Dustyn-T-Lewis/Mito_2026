#!/usr/bin/env Rscript
# F02 enrich volcanoes — per-contrast volcano-in-ring panels + the 2x2 composite,
# drawn by enrichVolcano from the reported DE (pi_score volcano) and the fgsea
# cache (NES arcs, deduped across the 5-DB lens). Individual rings go to
# b_reports/panels; the four-ring composite to b_reports/main/png.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(enrichVolcano)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))

BASE <- here::here("04_Figures_v2", "F02_Enrich_Volcanoes")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
PANELS <- file.path(BASE, "b_reports", "panels")
DAT <- file.path(BASE, "c_data")
for (d in c(RPT_PNG, PANELS, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

dep_df <- load_combined_wide()
fgsea_all <- readr::read_csv(
  here::here("04_Figures", "shared", "fgsea_tstat_all_h9c2.csv"),
  show_col_types = FALSE
)
rat_gene_sets <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))
set_pool <- do.call(c, unname(rat_gene_sets[CANONICAL_DBS]))

contrasts <- tibble::tribble(
  ~ctr, ~role, ~tag,
  "CTLvPHE", "Disease", "disease",
  "CTLvMITO", "Transplant", "transplant",
  "PHEvPHE_MITO", "Rescue", "rescue",
  "Interaction", "Interaction", "interaction"
)
RING_N <- 12

# Per contrast: tidy volcano (pi_score = significance) + ring terms
# (significant -> EnrichmentMap dedup -> top RING_N by FDR).
prep <- lapply(seq_len(nrow(contrasts)), function(i) {
  ctr <- contrasts$ctr[i]
  pi_col <- paste0("pi_score_", ctr)
  n_dep <- sum(dep_df[[pi_col]] < H9C2_PI_THRESH, na.rm = TRUE)
  volc <- dep_df |>
    transmute(
      gene,
      logFC = .data[[paste0("logFC_", ctr)]],
      P.Value = .data[[paste0("P.Value_", ctr)]],
      padj = .data[[pi_col]]
    ) |>
    filter(!is.na(logFC), !is.na(P.Value))
  sig <- fgsea_all |>
    filter(
      contrast == ctr, database %in% CANONICAL_DBS, padj < 0.05,
      !pathway %in% MITO_DROP_SETS
    ) |>
    as.data.frame()
  enrich <- deduplicate_enrichment(sig, set_pool, cutoff = 0.375) |>
    arrange(padj) |>
    slice_head(n = RING_N)
  list(
    volc = volc, enrich = enrich,
    subtitle = sprintf(
      "%s | %d DEPs, %d pathways", CONTRAST_MATH_BRIEF[ctr], n_dep, nrow(enrich)
    )
  )
})
names(prep) <- contrasts$role

# Shared ring styling for the panels and the composite.
ring_args <- list(
  gene_col = "gene", logfc_col = "logFC", pval_col = "P.Value", padj_col = "padj",
  term_col = "pathway", nes_col = "NES", size_col = "size",
  genes_col = "leadingEdge", genes_sep = ";",
  p_threshold = H9C2_PI_THRESH, logfc_threshold = log2(1.5),
  ring_radius = 5.5, volcano_radius = 5.2, arc_height_range = c(0.1, 3.2),
  label_size = 2.7, point_size = 2.0, point_alpha = 0.7,
  count_x_mult = 0.55, count_y_mult = 0.55
)

# Individual ring panels.
for (i in seq_len(nrow(contrasts))) {
  p <- do.call(volcano_ring, c(
    list(prep[[i]]$volc, prep[[i]]$enrich,
      title = contrasts$role[i], subtitle = prep[[i]]$subtitle
    ),
    ring_args
  ))
  ggplot2::ggsave(
    file.path(PANELS, sprintf("MAIN_F02_%s_ring.png", contrasts$tag[i])), p,
    width = 150, height = 150, units = "mm", dpi = 300
  )
}

# Four-ring 2x2 composite, shared NES legend collected as a bottom strip.
grid <- do.call(volcano_ring_grid, c(
  list(
    volc_dfs = lapply(prep, `[[`, "volc"),
    enrich_dfs = lapply(prep, `[[`, "enrich"),
    contrasts = contrasts$role,
    subtitles = vapply(prep, `[[`, character(1), "subtitle"),
    ncol = 2, theme = volcano_ring_theme(base_size = 13)
  ),
  ring_args
))
fig <- grid$plot &
  ggplot2::theme(
    legend.position = "bottom",
    legend.key.width = ggplot2::unit(13, "mm"),
    legend.key.height = ggplot2::unit(2.5, "mm"),
    legend.title = ggplot2::element_text(size = 8),
    legend.text = ggplot2::element_text(size = 7),
    legend.box.spacing = ggplot2::unit(0, "mm"),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    plot.margin = ggplot2::margin(1, 1, 1, 1, "mm")
  )
fig <- fig & ggplot2::guides(
  fill = ggplot2::guide_colorbar(direction = "horizontal", title.position = "top")
)
ggplot2::ggsave(
  file.path(RPT_PNG, "MAIN_F02_enrich_volcanoes.png"), fig,
  width = 270, height = 285, units = "mm", dpi = 300, limitsize = FALSE
)

# Supplementary workbook: contrast key + every tested pathway per contrast.
contrast_map <- contrasts |>
  mutate(
    brief = vapply(ctr, contrast_brief, character(1)),
    definition = unname(CONTRAST_MATH_BRIEF[ctr])
  ) |>
  select(old_name = ctr, file_tag = tag, role, brief, definition) |>
  as.data.frame()
pathway_sheets <- lapply(seq_len(nrow(contrasts)), function(i) {
  ctr <- contrasts$ctr[i]
  full <- fgsea_all |>
    filter(
      contrast == ctr, database %in% CANONICAL_DBS, !pathway %in% MITO_DROP_SETS
    ) |>
    transmute(pathway, database, padj, pval, NES, size,
      shown = pathway %in% prep[[i]]$enrich$pathway
    ) |>
    arrange(padj) |>
    as.data.frame()
  list(
    name = contrasts$role[i], df = full,
    role = sprintf("All tested pathways for the %s contrast", contrasts$role[i]),
    contents = "pathway, database, padj, pval, NES, size, shown (drawn on the ring); sorted by padj"
  )
})
build_workbook(
  file.path(DAT, "F02_supplementary.xlsx"),
  figure_title = "F02: Enrichment volcano-in-ring panels (5-DB lens)",
  sheet_specs = c(
    list(list(
      name = "contrast_map", df = contrast_map,
      role = "Key linking file tags / roles / contrast algebra",
      contents = "old combined-results name, file tag, role, brief, contrast definition"
    )),
    pathway_sheets
  )
)

message(sprintf("F02: %d ring panels + composite -> b_reports", nrow(contrasts)))
