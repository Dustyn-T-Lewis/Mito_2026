#!/usr/bin/env Rscript
# Pilot: did transplant just add mitochondrial proteins, or reprogram the cell?
# Three panels -- mito vs non-mito, program breakdown, context dependence -- read
# the reported DE and the rat MitoCarta lens, then stitch into one composite.

fns <- here::here("04_Figures", "functions")
source(file.path(fns, "01_style_palettes_theme.R"))
source(file.path(fns, "02_data_paths_and_loaders.R"))
source(file.path(fns, "06_supplementary_workbook.R"))
source(file.path(fns, "08_composite_layout.R"))
pacman::p_load(dplyr, tidyr, tibble, purrr, patchwork)

BASE <- here::here("04_Figures", "test", "pilot_mito_program")
for (f in list.files(file.path(BASE, "a_script", "panels"), full.names = TRUE)) source(f)
RPT <- file.path(BASE, "b_reports")
DAT <- file.path(BASE, "c_data")
for (d in c(RPT, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

comb <- load_combined_wide()

mito <- readRDS(here::here("04_Figures", "shared", "mitocarta3_rat.rds"))
mito_genes <- toupper(mito$MITOCARTA_ALL)
cat_map <- mito[grepl("__", names(mito))] |>
  imap(function(genes, nm) {
    tibble(gene = toupper(genes), category = sub("__.*$", "", sub("^MITOCARTA_", "", nm)))
  }) |>
  bind_rows() |>
  distinct(gene, category)

pa <- build_mito_vs_nonmito(comb, mito_genes)
pb <- build_program_breakdown(comb, cat_map)
pc <- build_context_dependence(comb, mito_genes)

fig <- add_tag(pa$plot, "A") | add_tag(pb$plot, "B") | add_tag(pc$plot, "C")
ggplot2::ggsave(file.path(RPT, "PILOT_mito_program.png"), fig,
  width = 250, height = 100, units = "mm", dpi = 300, limitsize = FALSE
)

build_workbook(
  file.path(DAT, "pilot_mito_program.xlsx"),
  figure_title = "Pilot: added mitochondrial cargo or proteome reprogramming?",
  sheet_specs = list(
    list(
      name = "mito_vs_nonmito", df = pa$table, role = "Panel A",
      contents = "median log2FC and n for MitoCarta vs other proteins per contrast"
    ),
    list(
      name = "program_breakdown", df = pb$table, role = "Panel B",
      contents = "median log2FC per MitoPathway category (structural vs regulatory), Transplant and Rescue"
    ),
    list(
      name = "context_dependence", df = pc$table, role = "Panel C",
      contents = "per mito protein: transplant effect at baseline vs under stress"
    )
  )
)

message("pilot_mito_program built")
