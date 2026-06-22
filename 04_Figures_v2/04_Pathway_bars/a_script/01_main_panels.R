#!/usr/bin/env Rscript
# F04 MAIN — Panel-D-style pathway count summary across 5 DBs.
# Per contrast, Up bar and Down bar drawn from 0; mito subset overdrawn darker.
# Pool = CANONICAL_DBS (Hallmark + Reactome + KEGG + MitoCarta + GO Slim);
# padj<0.05, size>=10, MITO_DROP_SETS excluded; cross-DB Jaccard+Overlap dedup
# at 0.375 per contrast (EnrichmentMap, Merico 2010 / Reimand 2019).
# Mito flag: database == MitoCarta OR pathway matches MITO_PATHWAY_REGEX.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures_v2", "functions", "04_mitocarta_lens_lookup.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))
source(here::here("04_Figures_v2", "04_Pathway_bars", "a_script", "_build.R"))

BASE <- here::here("04_Figures_v2", "04_Pathway_bars")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
DAT <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

out <- build_pathway_bar_panel()
p <- out$plot
bar_df <- out$bar_df
sig_pw <- out$sig_pw
FILL_TOTAL <- out$FILL_TOTAL
FILL_MITO <- out$FILL_MITO

# Inline 4-key legend (matches Panel D).
key_df <- tibble::tibble(
  y    = 4:1,
  lab  = c("Up total", "Up mito", "Down total", "Down mito"),
  fill = c(FILL_TOTAL["Up"], FILL_MITO["Up"], FILL_TOTAL["Down"], FILL_MITO["Down"])
)
p_key <- ggplot2::ggplot(key_df) +
  ggplot2::geom_point(
    ggplot2::aes(0, y),
    shape = 22, size = 2, fill = key_df$fill,
    color = "grey30", stroke = 0.3
  ) +
  ggplot2::geom_text(
    ggplot2::aes(0.22, y, label = lab),
    hjust = 0, size = 1.6, color = "grey20"
  ) +
  ggplot2::scale_x_continuous(limits = c(-0.1, 2.6)) +
  ggplot2::scale_y_continuous(limits = c(0.5, 4.5)) +
  ggplot2::theme_void()
p <- p + patchwork::inset_element(p_key, left = 0.60, right = 0.99, top = 0.99, bottom = 0.62)

FIG_W <- 120
FIG_H <- 70
ggsave(file.path(RPT_PDF, "MAIN_F04_pathway_bars.pdf"), p,
  width = FIG_W, height = FIG_H, units = "mm", device = pdf_dev, limitsize = FALSE
)
ggsave(file.path(RPT_PNG, "MAIN_F04_pathway_bars.png"), p,
  width = FIG_W, height = FIG_H, units = "mm", dpi = 300, limitsize = FALSE
)

# ---- Workbook ----------------------------------------------------------------
CORE <- H9C2_CONTRAST_ORDER

counts_tab <- bar_df |>
  dplyr::transmute(
    contrast = contrast_brief(as.character(contrast)),
    direction, total, mito,
    mito_fraction = ifelse(total > 0, mito / total, NA_real_)
  ) |>
  dplyr::arrange(contrast, direction)

per_ctr_tabs <- lapply(CORE, function(ctr) {
  sig_pw |>
    dplyr::filter(contrast == ctr) |>
    dplyr::transmute(
      pathway, database,
      clean_label = clean_display_label(pathway),
      NES = round(NES, 3), padj = signif(padj, 3),
      size, direction, is_mito
    ) |>
    dplyr::arrange(padj)
})

sheet_specs <- c(
  list(list(
    name = "dep_pathway_counts",
    df = counts_tab,
    role = "F04 bar heights — total significant pathways and mito subset per contrast x direction",
    contents = paste0(
      "contrast, direction (Up/Down), total significant pathways, mito subset, mito_fraction. ",
      "Pool: Hallmark + Reactome + KEGG + MitoCarta + GO Slim (CANONICAL_DBS); ",
      "padj<0.05, size>=10, MITO_DROP_SETS excluded; cross-DB Jaccard+Overlap dedup at 0.375 ",
      "(EnrichmentMap, Merico 2010 / Reimand 2019). Mito flag: database==MitoCarta OR ",
      "pathway matches MITO_PATHWAY_REGEX (mitochondria-keyword regex in 04_mitocarta_lens_lookup.R)."
    )
  )),
  Map(
    function(ctr, tab) {
      list(
        name     = paste0(contrast_brief(ctr), "_sig_pathways"),
        df       = tab,
        role     = sprintf("Post-dedup significant pathways behind the %s bars", contrast_brief(ctr)),
        contents = "pathway, database, clean_label, NES, padj, size, direction (Up/Down), is_mito"
      )
    },
    CORE, per_ctr_tabs
  )
)

build_workbook(
  file.path(DAT, "F04_supplementary.xlsx"),
  figure_title = "F04 — Panel-D pathway count summary (5-DB pool, mito subset overlaid)",
  sheet_specs  = sheet_specs
)

# Remove the legacy shown_pathways.csv if a prior run left it behind — no consumer now.
old_csv <- file.path(DAT, "shown_pathways.csv")
if (file.exists(old_csv)) file.remove(old_csv)

message(sprintf("F04 done | %d sig pathways shown | output %s", nrow(sig_pw), RPT_PNG))
