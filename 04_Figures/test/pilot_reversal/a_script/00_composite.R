#!/usr/bin/env Rscript
# Pilot: is the PHE proteome shift reversed by transplant? Three panels test the
# reversal at rising strictness -- per-protein return to baseline, whole-proteome
# distance, and the formal interaction -- then stitch into one composite.

fns <- here::here("04_Figures", "functions")
source(file.path(fns, "01_style_palettes_theme.R"))
source(file.path(fns, "02_data_paths_and_loaders.R"))
source(file.path(fns, "06_supplementary_workbook.R"))
source(file.path(fns, "08_composite_layout.R"))
pacman::p_load(dplyr, tibble, patchwork)

BASE <- here::here("04_Figures", "test", "pilot_reversal")
for (f in list.files(file.path(BASE, "a_script", "panels"), full.names = TRUE)) source(f)
RPT <- file.path(BASE, "b_reports")
DAT <- file.path(BASE, "c_data")
for (d in c(RPT, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

comb <- load_combined_wide()
dal <- readRDS(P05$imp_rds)

pa <- build_return_to_baseline(comb)
pb <- build_proteome_distance(as.matrix(dal$data), dal$metadata$group)
pc <- build_interaction_null(comb)

fig <- add_tag(pa$plot, "A") | add_tag(pb$plot, "B") | add_tag(pc$plot, "C")
ggplot2::ggsave(file.path(RPT, "PILOT_reversal.png"), fig,
  width = 250, height = 95, units = "mm", dpi = 300, limitsize = FALSE
)

build_workbook(
  file.path(DAT, "pilot_reversal.xlsx"),
  figure_title = "Pilot: is the PHE proteome shift reversed by transplant?",
  sheet_specs = list(
    list(
      name = "return_to_baseline", df = pa$table, role = "Panel A",
      contents = "disease-signature proteins: disease vs recovery log2FC, attenuation, returned flag"
    ),
    list(
      name = "proteome_distance", df = pb$table, role = "Panel B",
      contents = "centroid distances to Ctl, the ratio, and the label-permutation p"
    ),
    list(
      name = "interaction", df = pc$table, role = "Panel C",
      contents = "interaction contrast per protein: log2FC, p, FDR, hit flag"
    )
  )
)

message("pilot_reversal built")
