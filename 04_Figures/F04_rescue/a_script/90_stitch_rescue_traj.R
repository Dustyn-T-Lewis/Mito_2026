# F04 Rescue — Companion composite: trajectory clusters + residual volcano
# Pairs the rescue trajectory (what moved back toward control) with the residual
# (what the transplant did NOT return to control). Reads the panel PNGs (rebuilds
# them first so it's always current). Output: PNG only (cairo-less env).
library(here)
source(here::here("04_Figures/shared/config.R"))         # sources style.R
suppressPackageStartupMessages({ library(patchwork); library(ggplot2); library(png); library(grid) })

PANELS <- here::here("04_Figures/F04_rescue/b_reports/main/png/panels")
RPT    <- here::here("04_Figures/F04_rescue/b_reports/main/png")
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

source(here::here("04_Figures/F04_rescue/a_script/panel_traj.R"))
source(here::here("04_Figures/F04_rescue/a_script/panel_resid.R"))

read_panel <- function(f) {
  p <- file.path(PANELS, f)
  if (!file.exists(p)) stop("missing panel PNG: ", p)
  ggplot() + annotation_custom(rasterGrob(readPNG(p), interpolate = TRUE)) +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
}

fig <- read_panel("MAIN_panel_traj_composite.png") /
       read_panel("MAIN_panel_resid_composite.png") +
  plot_layout(heights = c(150, 120)) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave(file.path(RPT, "MAIN_F04_rescue_trajectory_composite.png"), fig,
       width = 220, height = 280, units = "mm", dpi = 300)
message("Mito rescue trajectory+residual companion composite saved")
