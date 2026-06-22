#!/usr/bin/env Rscript
# F02 MAIN — DEP bars + effect-size companion. Standalone (no composite).
# Reads existing pipeline outputs only; never re-runs 01-03.
#   Left  : DEP counts as % of proteome, per contrast, at three independent
#           thresholds (p<0.05 / FDR<0.10 / Π<0.05), Up/Down dodged,
#           faint->solid by threshold.
#   Right : |log2FC| distribution per contrast — histogram + density faceted
#           by contrast, median |log2FC| (over all proteins) annotated per facet.

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
p_counts <- build_dep_count_panel()
p_eff <- build_dep_effect_panel()

# ---- Compose -----------------------------------------------------------------
fig <- p_counts + p_eff + plot_layout(widths = c(1, 0.55))

FIG_W <- 178
FIG_H <- 100
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
    df       = dep_count_data() |> select(contrast, threshold, direction, n, pct),
    role     = "Panel A bar heights — DEP counts per contrast",
    contents = "Up/Down significant-protein counts and % of proteome at each independent threshold (p<0.05, FDR, Π<0.05) per contrast"
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
  select(contrast, direction, n) |>
  pivot_wider(names_from = direction, values_from = n)
message(sprintf("F02 MAIN (%dx%d mm) | proteome n=%d", FIG_W, FIG_H, n_total))
print(as.data.frame(pi_summary))
