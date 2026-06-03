#!/usr/bin/env Rscript
# F04 orchestrator: renders main panels, bundles per-panel CSVs into the
# supplementary xlsx, optional Box copy. Mirrors F03 stitch.

source(here::here("05_Figures", "F04_rescue", "a_script", "01_main_panels.R"))
source(here::here("04_Figures", "shared", "figure_supplement_helpers.R"))

DAT <- here::here("05_Figures", "F04_rescue", "c_data")
f04_specs <- list(
  list(name = "panel_A_quadrant_proteins", path = file.path(DAT, "panel_A", "quadrant_proteins.csv")),
  list(name = "panel_A_labelled",          path = file.path(DAT, "panel_A", "labelled_proteins.csv")),
  list(name = "panel_B_pattern_class",     path = file.path(DAT, "panel_B_heatmap", "pattern_classification.csv")),
  list(name = "panel_B_sankey",            path = file.path(DAT, "panel_B_heatmap", "sankey_links.csv")),
  list(name = "panel_B_bar",               path = file.path(DAT, "panel_B_heatmap", "bar_data.csv")),
  list(name = "panel_D_nes_scatter",       path = file.path(DAT, "panel_D", "nes_scatter.csv")),
  list(name = "panel_E_rrho2_summary",     path = file.path(DAT, "panel_E", "rrho2_summary.csv")),
  list(name = "panel_E_rrho2_hotspot",     path = file.path(DAT, "panel_E", "rrho2_hotspot_genes.csv")),
  list(name = "panel_E_rrho2_ora_reversed",path = file.path(DAT, "panel_E", "rrho2_ora_concordant.csv")))
build_workbook(file.path(DAT, "F04_supplementary.xlsx"), sheet_specs = f04_specs)

box <- Sys.getenv("H9C2_BOX_DIR", "")
if (nzchar(box) && dir.exists(box)) {
  RPT <- here::here("05_Figures", "F04_rescue", "b_reports")
  for (sub in c("figures/pdf", "figures/png", "tables")) dir.create(file.path(box, "02_Figures", sub), recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(RPT, "main", "pdf", "MAIN_F04_composite.pdf"), file.path(box, "02_Figures", "figures", "pdf", "Figure_4_rescue.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "main", "png", "MAIN_F04_composite.png"), file.path(box, "02_Figures", "figures", "png", "Figure_4_rescue.png"), overwrite = TRUE)
  file.copy(file.path(DAT, "F04_supplementary.xlsx"), file.path(box, "02_Figures", "tables", "Table_F04.xlsx"), overwrite = TRUE)
  message("Copied F04 outputs to Box")
}
message("F04 complete")
