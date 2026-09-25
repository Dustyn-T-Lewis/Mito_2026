#!/usr/bin/env Rscript
# Copies the finished figures and workbooks into Supplementary_Material/ under the names
# the manuscript cites. Pipeline F01-F03 are manuscript Figures 4-6; Figures 1-3 come from
# the non-proteomic experiments and are not produced here. Run after every figure script.

pacman::p_load(here, purrr)

fig <- \(...) here::here("04_Figures", ...)
out <- here::here("Supplementary_Material")

MAIN <- c(
  Figure_4.png = fig("F01_Proteome_Overview", "b_reports", "MAIN_F01_proteome_overview.png"),
  Figure_5.png = fig("F02_Enrich_Volcanoes", "b_reports", "MAIN_F02_enrich_volcanoes.png"),
  Figure_6.png = fig("F03_WGCNA", "b_reports", "MAIN_F03_wgcna.png")
)
SUPP_FIGURES <- c(
  S1_Figure.pdf = fig("F01_Proteome_Overview", "b_reports", "S1_Figure.pdf"),
  S2_Figure.pdf = fig("F02_Enrich_Volcanoes", "b_reports", "S2_Figure.pdf"),
  S3_Figure.pdf = fig("F03_WGCNA", "b_reports", "S3_Figure.pdf"),
  S4_Figure.pdf = fig("F03_WGCNA", "b_reports", "S4_Figure.pdf")
)
SUPP_TABLES <- c(
  S1_Table.xlsx = fig("F01_Proteome_Overview", "c_data", "F01_supplementary.xlsx"),
  S2_Table.xlsx = fig("F02_Enrich_Volcanoes", "c_data", "F02_supplementary.xlsx"),
  S3_Table.xlsx = fig("F03_WGCNA", "c_data", "F03_supplementary.xlsx")
)

sets <- list(Main_Figures = MAIN, Supplementary_Figures = SUPP_FIGURES, Supplementary_Tables = SUPP_TABLES)
missing <- unlist(sets)[!file.exists(unlist(sets))]
if (length(missing)) stop("run the figure scripts first; missing:\n  ", paste(missing, collapse = "\n  "))

iwalk(sets, function(files, subdir) {
  dest <- file.path(out, subdir)
  unlink(dest, recursive = TRUE)
  dir.create(dest, recursive = TRUE)
  file.copy(files, file.path(dest, names(files)))
})
message("Supplementary_Material: ", paste(lengths(sets), names(sets), collapse = ", "))
