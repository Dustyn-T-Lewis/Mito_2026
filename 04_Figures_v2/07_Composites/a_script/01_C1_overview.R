#!/usr/bin/env Rscript
# C1 overview — PCA + DEP bars + effect histogram + pathway bars + Venn. Who differs.
library(here)
source(here::here("04_Figures_v2", "functions", "08_composite_layout.R"))
source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures_v2", "functions", "04_mitocarta_lens_lookup.R"))
source(here::here("04_Figures_v2", "01_PCA", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "02_DEP_bars", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "03_Venn", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "04_Pathway_bars", "a_script", "_build.R"))

BASE <- here::here("04_Figures_v2", "07_Composites")

# A: PCA (fills its cell so the top row aligns top and bottom)
pca <- build_pca_panel()$plot

# B: DEP bars (clear baked tag)
dep <- build_dep_count_panel() + ggplot2::labs(tag = NULL)

# C: effect histogram (clear baked tag)
eff <- build_dep_effect_panel() + ggplot2::labs(tag = NULL)

# D: pathway bars
pw <- build_pathway_bar_panel()$plot

# E + F: Venn + direction strip
venn <- build_venn_panels()
venn$venn <- venn$venn + ggplot2::labs(subtitle = NULL)

# Layout: square 178×178 like YvO F01.
# Top row: PCA · Venn · direction strip. Bottom row: DEP counts · effect size · pathways.
design <- "
AABBCC
DDEEFF
"

fig <- add_tag(pca, "A") +
  add_tag(venn$venn, "B") +
  add_tag(venn$strip, "C") +
  add_tag(dep, "D") +
  add_tag(eff, "E") +
  add_tag(pw, "F") +
  patchwork::plot_layout(
    design  = design,
    widths  = c(1.05, 1.0, 1.0),
    heights = c(1.0, 1.0)
  ) +
  patchwork::plot_annotation(
    title = "H9c2 mito-transplant proteome — overview",
    subtitle = "2×2 factorial · n = 24 (6/group) · 4,806 proteins (DIA-MS)",
    caption = paste(
      "Disease = PHE − Ctl; Transplant = Mito − Ctl; Rescue = PHE_Mito − PHE.",
      "Significance Π < 0.05 (Xiao 2014). Up = red, Down = blue."
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 7),
      plot.subtitle = ggplot2::element_text(
        face = "italic", size = 5, color = "grey30"
      ),
      plot.caption = ggplot2::element_text(
        size = 5, color = "grey35", hjust = 0, lineheight = 1.1
      )
    )
  )

save_composite(fig, BASE, "MAIN_C1_overview", width_mm = PANEL_MD, height_mm = PANEL_MD)
message("C1 composite built")
