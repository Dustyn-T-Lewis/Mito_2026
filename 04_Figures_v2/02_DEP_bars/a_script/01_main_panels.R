#!/usr/bin/env Rscript
# F02 MAIN — DEP bars + effect-size companion. Standalone (no composite).
# Reads existing pipeline outputs only; never re-runs 01-03.
#   Left  : horizontal DEP-count bars, one per contrast, x = % of proteome,
#           nested p/FDR/Π drawn as overlapping identity bars, contrast-colored.
#   Right : effect-size histograms faceted by contrast (contrast-colored),
#           median |log2FC| over all proteins annotated per facet.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))
source(here::here("04_Figures_v2", "02_DEP_bars", "a_script", "_build.R"))

BASE <- here::here("04_Figures_v2", "02_DEP_bars")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
DAT <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

# ---- Build panels ------------------------------------------------------------
pDEP <- build_dep_count_panel()
pHIST <- build_dep_effect_panel()

# ---- Compose -----------------------------------------------------------------
fig <- pDEP + pHIST + plot_layout(widths = c(1, 0.62))

FIG_W <- 178
FIG_H <- 120
ggsave(file.path(RPT_PDF, "MAIN_F02_dep_bars.pdf"), fig,
  width = FIG_W, height = FIG_H, units = "mm",
  device = pdf_dev, limitsize = FALSE
)
ggsave(file.path(RPT_PNG, "MAIN_F02_dep_bars.png"), fig,
  width = FIG_W, height = FIG_H, units = "mm", dpi = 300, limitsize = FALSE
)

# ---- Supplementary workbook --------------------------------------------------
dep_tabs <- lapply(CORE, function(c) {
  dep_results[[c]] |>
    transmute(uniprot_id, gene, logFC, P.Value, adj.P.Val, pi_score) |>
    arrange(pi_score)
})
sheet_specs <- c(
  list(list(
    name     = "dep_counts",
    df       = dep_count_data() |> select(contrast, threshold, n, pct),
    role     = "Panel A bar heights — DEP counts per contrast",
    contents = "Significant-protein counts and % of proteome at each independent threshold (p<0.05, FDR, Π<0.05) per contrast"
  )),
  Map(
    function(c, tab) {
      list(
        name = substr(contrast_brief(c), 1, 28),
        df = tab,
        role = "Underlying per-protein DE table feeding both panels",
        contents = sprintf(
          "Every protein tested in the %s contrast: uniprot_id, gene, logFC, P.Value, adj.P.Val, pi_score (sorted by Π)",
          contrast_brief(c)
        )
      )
    },
    CORE, dep_tabs
  )
)
build_workbook(
  file.path(DAT, "F02_supplementary.xlsx"),
  figure_title = "F02 — DEP counts by contrast + effect-size (|log2FC|) companion",
  sheet_specs  = sheet_specs
)

pi_summary <- dep_count_data() |>
  filter(threshold == THR_LEVELS[3]) |>
  select(contrast, n, pct)
message(sprintf("F02 MAIN (%dx%d mm) | proteome n=%d", FIG_W, FIG_H, n_total))
print(as.data.frame(pi_summary))
