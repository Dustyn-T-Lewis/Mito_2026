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

# A: PCA — force square aspect ratio in composite cell
pca <- build_pca_panel()$plot +
  ggplot2::theme(aspect.ratio = 1)

# B: horizontal DEP bars — strip the baked standalone tag before re-tagging
dep <- build_dep_count_panel() + ggplot2::labs(tag = NULL)

# C: effect-size histogram — strip baked tag "B"
eff <- build_dep_effect_panel() + ggplot2::labs(tag = NULL)

# D: pathway count bars
pw <- build_pathway_bar_panel()$plot

# E + F: Venn diagram and direction strip
venn <- build_venn_panels()
venn$venn <- venn$venn + ggplot2::labs(subtitle = NULL)

# Layout: 2 rows x 3 cols
# Top:    A (PCA square)  |  B (DEP bars)  |  C (histogram)
# Bottom: D (pathway bars) | E (Venn)       | F (direction strip)
# Widths: PCA and Venn get equal share; DEP bars slightly wider than histogram;
# pathway bars match the PCA column.
# Heights: top row slightly taller to give PCA room to be square.
design <- "
AABBCC
DDEEFF
"

fig <- add_tag(pca, "A") +
  add_tag(dep, "B") +
  add_tag(eff, "C") +
  add_tag(pw, "D") +
  add_tag(venn$venn, "E") +
  add_tag(venn$strip, "F") +
  patchwork::plot_layout(
    design  = design,
    widths  = c(1.1, 1.1, 0.8),
    heights = c(1.15, 1.0)
  ) +
  composite_caption(paste(
    "Disease = PHE − Ctl; Transplant = Mito − Ctl; Rescue = PHE_Mito − PHE.",
    "Significance Π < 0.05 (Xiao 2014). Up = red, Down = blue."
  ))

save_composite(fig, BASE, "MAIN_C1_overview", width_mm = PANEL_MD, height_mm = 210)
message("C1 composite built")
