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

# B: DEP counts + a small stringency key (light -> dark = p / FDR / Π), top-right
dep_key_df <- tibble::tibble(
  y = 3:1,
  lab = c("p < 0.05", "FDR < 0.10", "Π < 0.05"),
  a = c(0.18, 0.45, 1)
)
dep_key <- ggplot2::ggplot(dep_key_df) +
  ggplot2::geom_tile(
    ggplot2::aes(0, y),
    width = 0.5, height = 0.6,
    fill = "grey20", alpha = dep_key_df$a, color = "black", linewidth = 0.2
  ) +
  ggplot2::geom_text(
    ggplot2::aes(0.38, y, label = lab),
    hjust = 0, size = 1.5, fontface = "bold", color = "grey15"
  ) +
  ggplot2::scale_x_continuous(limits = c(-0.35, 3.1)) +
  ggplot2::scale_y_continuous(limits = c(0.45, 3.55)) +
  ggplot2::theme_void()
dep <- add_tag(build_dep_count_panel(), "B") +
  patchwork::inset_element(
    dep_key,
    left = 0.67, right = 1.0, top = 1.0, bottom = 0.75
  )

eff <- add_tag(build_dep_effect_panel(), "C")

# Strip stays aligned with the DEP-counts panel above it (same column width).
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
    widths  = c(1.05, 1.0, 1.0),
    heights = c(1.0, 1.0)
  ) +
  patchwork::plot_annotation(
    caption = paste(
      "Disease = PHE − Ctl; Transplant = Mito − Ctl; Rescue = PHE_Mito − PHE.",
      "Significance Π < 0.05 (Xiao 2014). Up = red, Down = blue."
    ),
    theme = ggplot2::theme(
      plot.caption.position = "plot",
      plot.caption = ggplot2::element_text(
        size = 5, color = "grey35", hjust = 0, lineheight = 1.1
      )
    )
  )

save_composite(fig, BASE, "MAIN_C1_overview", width_mm = PANEL_MD, height_mm = 150)
message("C1 composite built")
