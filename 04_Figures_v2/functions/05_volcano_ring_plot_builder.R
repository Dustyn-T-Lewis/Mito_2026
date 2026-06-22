# Compatibility stub: volcano-ring composite plot utility split across siblings.
# Standard Cartesian ggplot with ggforce::geom_arc_bar(); NO coord_polar().
#   05a — clean_ring_label (MitoCarta + generic label cleaner)
#   05b — prepare_ring_data, build_tick_data
#   05c — build_volcano_layers
#   05d — build_ring_layers, build_label_layer, select/center utilities,
#         build_ring_180_split, build_ring_with_gaps
#   05e — make_volcano_ring (top-level composer), build_nes_legend_bar

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggforce)
  library(patchwork)
  library(scales)
})

if (!exists("FIG_THEME")) {
  source(here::here("04_Figures_v2", "functions", "01_style_palettes_theme.R"))
}

source(here::here("04_Figures_v2", "functions", "05a_ring_labels.R"))
source(here::here("04_Figures_v2", "functions", "05b_ring_data.R"))
source(here::here("04_Figures_v2", "functions", "05c_volcano_layer.R"))
source(here::here("04_Figures_v2", "functions", "05d_ring_layers.R"))
source(here::here("04_Figures_v2", "functions", "05e_composite.R"))
