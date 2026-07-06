#!/usr/bin/env Rscript
# Supplementary GSEA figures for the enrich-volcano contrasts, drawn from the shared
# fgsea cache. Per contrast the top 30 up- and down-regulated pathways across the ring
# databases sit back to back, named inside their bars with FDR at each tip. The
# leading-edge genes behind every shown pathway append to the F02 workbook.

library(here)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
fns <- here::here("04_Figures", "functions")
source(file.path(fns, "02_data_paths_and_loaders.R"))
source(file.path(fns, "03_pathway_enrichment_dedup_ora.R"))
source(file.path(fns, "06_supplementary_workbook.R"))

BASE <- here::here("04_Figures", "F02_Enrich_Volcanoes")
PNG <- file.path(BASE, "b_reports", "supplementary")
DAT <- file.path(BASE, "c_data")
dir.create(PNG, recursive = TRUE, showWarnings = FALSE)

# the four primary 2x2 contrasts, named via the shared map (drops the Secondary contrast)
CONTRASTS <- CONTRAST_DISPLAY_MAP[c("CTLvPHE", "CTLvMITO", "PHEvPHE_MITO", "Interaction")]
TOP_N <- 30L

gsea <- readr::read_csv(
  here::here("04_Figures", "shared", "c_data", "fgsea_tstat_all_h9c2.csv"),
  show_col_types = FALSE
) |>
  filter(
    contrast %in% names(CONTRASTS), database %in% CANONICAL_DBS,
    !pathway %in% MITO_DROP_SETS, !grepl(DISEASE_VIRAL_RE, pathway, ignore.case = TRUE)
  ) |>
  mutate(
    contrast = factor(unname(CONTRASTS[contrast]), levels = CONTRASTS),
    database = factor(database, levels = CANONICAL_DBS)
  )

# top pathways one direction at a time, ordered so the strongest sits at the top
top_pathways <- function(g, direction) {
  g <- if (direction == "up") filter(g, NES > 0) else filter(g, NES < 0)
  g |>
    slice_max(if (direction == "up") NES else -NES, n = TOP_N, with_ties = FALSE) |>
    arrange(if (direction == "up") NES else desc(NES))
}

side_plot <- function(g, direction) {
  d <- mutate(g,
    rid = factor(paste(database, pathway), levels = paste(database, pathway)),
    name = clean_display_label(pathway),
    name = ifelse(nchar(name) > 32, paste0(substr(name, 1, 30), "…"), name),
    fdr = ifelse(padj < 0.001, formatC(padj, format = "e", digits = 0), formatC(padj, format = "f", digits = 3))
  )
  up <- direction == "up"
  ggplot(d, aes(NES, rid, fill = database)) +
    geom_col(width = 0.8, colour = "black", linewidth = 0.25) +
    geom_text(aes(label = name),
      x = if (up) 0.05 else -0.05, hjust = if (up) 0 else 1,
      size = 2.5, fontface = "bold", colour = "white"
    ) +
    geom_text(aes(label = fdr), hjust = if (up) -0.18 else 1.18, size = 1.95, fontface = "bold", colour = "grey20") +
    scale_fill_manual(values = DB_COLORS, name = NULL, drop = FALSE) +
    scale_y_discrete(labels = NULL) +
    scale_x_continuous(
      limits = if (up) c(0, max(d$NES)) else c(min(d$NES), 0),
      expand = expansion(mult = if (up) c(0.02, 0.16) else c(0.16, 0.02))
    ) +
    coord_cartesian(clip = "off") +
    labs(title = if (up) "Up-regulated" else "Down-regulated", x = "NES", y = NULL) +
    FIG_THEME +
    theme(
      plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), panel.grid.major.y = element_blank()
    )
}

for (cn in unname(CONTRASTS)) {
  g <- filter(gsea, contrast == cn)
  fig <- (side_plot(top_pathways(g, "down"), "down") | side_plot(top_pathways(g, "up"), "up")) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("GSEA: %s — top %d up / %d down by NES, no dedup", cn, TOP_N, TOP_N),
      theme = theme(plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE))
    ) &
    theme(legend.position = "bottom")
  ggsave(file.path(PNG, sprintf("SUPP_F02_gsea_%s.png", cn)), fig,
    width = PANEL_MD, height = 235, units = "mm", dpi = 300
  )
}

# leading-edge overlap behind every shown pathway
overlap <- bind_rows(lapply(unname(CONTRASTS), function(cn) {
  g <- filter(gsea, contrast == cn)
  bind_rows(up = top_pathways(g, "up"), down = top_pathways(g, "down"), .id = "direction") |>
    transmute(
      contrast = cn, direction, database,
      pathway = clean_display_label(pathway),
      NES = round(NES, 2), FDR = signif(padj, 3),
      set_size = size, leading_edge_n = lengths(strsplit(leadingEdge, ";")),
      overlap_coef = round(leading_edge_n / set_size, 3),
      leading_edge = leadingEdge
    )
}))

append_workbook(
  file.path(DAT, "F02_supplementary.xlsx"),
  sheet_specs = list(
    list(
      name = "gsea_overlap", df = overlap,
      role = "Supplement: per-contrast GSEA figures",
      contents = "shown pathways with NES, FDR, set size, and the leading-edge genes driving each enrichment"
    ),
    list(
      name = "gsea_full",
      df = gsea |>
        transmute(
          contrast = as.character(contrast), database, pathway,
          NES = round(NES, 3), pval = signif(pval, 3), padj = signif(padj, 3),
          size, leading_edge = leadingEdge
        ) |>
        arrange(contrast, database, padj),
      role = "Supplement: all fgsea results",
      contents = "every pathway tested per contrast and database (no dedup)"
    )
  )
)

message(sprintf("F02 GSEA supplement: %d figures, %d shown pathways", length(CONTRASTS), nrow(overlap)))
