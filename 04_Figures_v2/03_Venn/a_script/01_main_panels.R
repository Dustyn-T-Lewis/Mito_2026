#!/usr/bin/env Rscript
# F03 Venn — driver. Sources _build.R for all logic; this file only wires
# outputs to disk.

suppressPackageStartupMessages({
  library(ggplot2)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))
source(here::here("04_Figures_v2", "03_Venn", "a_script", "_build.R"))

BASE <- here::here("04_Figures_v2", "03_Venn")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
DAT <- file.path(BASE, "c_data")
for (d in c(RPT_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

out <- build_venn_panels()
venn <- out$venn
strip <- out$strip

have_ggplotify <- requireNamespace("ggplotify", quietly = TRUE)
have_patchwork <- requireNamespace("patchwork", quietly = TRUE)

if (have_ggplotify && have_patchwork) {
  combined <- patchwork::wrap_plots(venn, strip, widths = c(1.6, 1))
  ggsave(file.path(RPT_PNG, "MAIN_F03_venn.png"), combined,
    width = 140, height = 100, units = "mm", dpi = 300
  )
} else {
  # fallback: venn is still the raw grob from eulerr
  sub_txt <- sprintf(
    "Disease: %s  |  Transplant: %s  |  Rescue: %s",
    CONTRAST_MATH_BRIEF[["CTLvPHE"]],
    CONTRAST_MATH_BRIEF[["CTLvMITO"]],
    CONTRAST_MATH_BRIEF[["PHEvPHE_MITO"]]
  )
  png(file.path(RPT_PNG, "MAIN_F03_venn.png"),
    width = 140, height = 100,
    units = "mm", res = 300
  )
  grid::grid.newpage()
  grid::grid.text("DEP overlap (Π < 0.05)",
    y = 0.97,
    gp = grid::gpar(fontface = "bold", fontsize = 8)
  )
  grid::grid.text(sub_txt, y = 0.93, gp = grid::gpar(fontsize = 5, col = "grey30"))
  print(venn, newpage = FALSE)
  dev.off()
  ggsave(file.path(RPT_PNG, "MAIN_F03_venn_direction.png"), strip,
    width = 80, height = 70, units = "mm", dpi = 300
  )
}

build_workbook(
  file.path(DAT, "F03_supplementary.xlsx"),
  figure_title = "F03 — DEP overlap (Π < 0.05) across Disease / Transplant / Rescue",
  sheet_specs = list(
    list(
      name     = "membership",
      df       = out$membership,
      role     = "Per-protein set membership behind the Euler diagram",
      contents = "One row per significant protein (Π<0.05): TRUE/FALSE membership in Disease/Transplant/Rescue, the Up/Down direction in each set, and its region_key"
    ),
    list(
      name     = "region_counts",
      df       = out$region_counts,
      role     = "The 7 Euler region areas + the directional strip counts",
      contents = "Protein count per overlap region (region_key), sorted descending"
    ),
    list(
      name     = "euler_fit",
      df       = out$fit_stats,
      role     = "Euler fit quality — stress and diagError; high values mean the area-proportional layout is approximate",
      contents = "Two rows: stress and diagError from the eulerr fit object"
    )
  )
)

print(out$region_counts)
