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

# Panels — add_tag bakes the letter into each title for uniform spacing.
pca <- add_tag(build_pca_panel()$plot, "A")

# B: DEP counts + a stringency key — terms boxed tight and shaded, ordered
# dark -> light (Π / FDR / p) to read parallel with the bars, bottom-right.
dep_key_df <- tibble::tibble(
  x = c(1, 2.15, 3.2),
  lab = c("Π < 0.05", "FDR < 0.10", "p < 0.05"),
  a = c(1, 0.45, 0.18),
  txt = c("white", "grey15", "grey15")
)
dep_key <- ggplot2::ggplot(dep_key_df) +
  ggplot2::geom_label(
    ggplot2::aes(x, 0, label = lab, alpha = a, color = I(txt)),
    fill = "grey20", size = 1.8, fontface = "bold",
    label.padding = ggplot2::unit(0.5, "mm"), label.size = 0.25
  ) +
  ggplot2::scale_alpha_identity() +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.4)) +
  ggplot2::theme_void()
dep <- add_tag(build_dep_count_panel(), "B") +
  patchwork::inset_element(
    dep_key,
    left = 0.4, right = 1.0, top = 0.16, bottom = 0.02
  )

eff <- add_tag(build_dep_effect_panel(), "C")

# Strip shares column 2 with the DEP-counts panel, so patchwork aligns it to
# the same plot width and left edge as panel B.
venn <- build_venn_panels()
venn_gg <- add_tag(venn$venn, "D")
strip <- add_tag(venn$strip, "E")

pw <- add_tag(build_pathway_bar_panel()$plot, "F")

# Layout: square 178×178 like YvO F01.
# Top row: PCA · DEP counts · effect size. Bottom row: Venn · direction strip · pathways.
design <- "
AABBCC
DDEEFF
"

fig <- pca + dep + eff + venn_gg + strip + pw +
  patchwork::plot_layout(
    design  = design,
    widths  = c(1.3, 0.88, 0.82),
    heights = c(1.0, 1.0)
  )

save_composite(fig, BASE, "MAIN_C1_overview", width_mm = PANEL_MD, height_mm = 150)
message("C1 composite built")
