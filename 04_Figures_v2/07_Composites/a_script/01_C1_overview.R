#!/usr/bin/env Rscript
# C1 overview — PCA + DEP bars + Venn. Who differs.
library(here)
source(here::here("04_Figures_v2", "functions", "08_composite_layout.R"))
source(here::here("04_Figures_v2", "01_PCA", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "02_DEP_bars", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "03_Venn", "a_script", "_build.R"))

BASE <- here::here("04_Figures_v2", "07_Composites")

pca <- build_pca_panel()$plot
# Strip the embedded tags the builders hard-code (tag = "A"/"B") before
# add_tag() re-applies the composite letters C, D. Both builders return
# patchworks (inset_element wraps the bar chart); tag lives on [[1]].
dep_raw <- build_dep_count_panel()
dep_raw[[1]] <- dep_raw[[1]] + labs(tag = NULL)
dep <- dep_raw
eff <- build_dep_effect_panel() + labs(tag = NULL)
venn <- build_venn_panels()
venn$venn <- venn$venn + labs(subtitle = NULL)

# Layout: PCA dominant top-left (3 rows x 2 cols), Venn top-right,
# DEP counts mid-right, direction strip bottom spanning full width.
# Heights: PCA/Venn row slightly taller; effect histograms compact; strip narrow.
design <- "
AABB
AACC
DDEE
"

fig <- add_tag(pca, "A") + add_tag(venn$venn, "B") +
  add_tag(dep, "C") + add_tag(eff, "D") + add_tag(venn$strip, "E") +
  plot_layout(
    design  = design,
    heights = c(1.05, 1.00, 0.75)
  ) +
  composite_caption(paste(
    "Disease = PHE − Ctl; Transplant = Mito − Ctl; Rescue = PHE_Mito − PHE.",
    "Significance Π < 0.05 (Xiao 2014). Up = red, Down = blue."
  ))

save_composite(fig, BASE, "MAIN_C1_overview", width_mm = PANEL_MD, height_mm = 195)
message("C1 composite built")
