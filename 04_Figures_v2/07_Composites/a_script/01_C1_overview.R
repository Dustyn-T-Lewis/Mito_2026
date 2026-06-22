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
# free() the PCA left side so its narrow y-axis isn't stretched to match the
# wide contrast labels of the DEP-counts panel below it in the same column.
pca <- patchwork::free(add_tag(build_pca_panel()$plot, "A"), side = "l")

venn <- build_venn_panels()
venn_gg <- add_tag(venn$venn, "B")
strip <- add_tag(venn$strip, "C")

# D: DEP counts + a small stringency key (light -> dark = p / FDR / Π), top-right
dep_key_df <- tibble::tibble(
  y = 3:1,
  lab = c("p < 0.05", "FDR < 0.10", "Π < 0.05"),
  a = c(0.18, 0.45, 1)
)
dep_key <- ggplot2::ggplot(dep_key_df) +
  ggplot2::geom_tile(
    ggplot2::aes(0, y),
    width = 0.7, height = 0.78,
    fill = "grey20", alpha = dep_key_df$a, color = "black", linewidth = 0.2
  ) +
  ggplot2::geom_text(
    ggplot2::aes(0.5, y, label = lab),
    hjust = 0, size = 1.5, color = "grey15"
  ) +
  ggplot2::scale_x_continuous(limits = c(-0.5, 3.4)) +
  ggplot2::scale_y_continuous(limits = c(0.3, 3.7)) +
  ggplot2::theme_void()
dep <- add_tag(build_dep_count_panel(), "D") +
  patchwork::inset_element(
    dep_key,
    left = 0.6, right = 0.995, top = 0.99, bottom = 0.72
  )

eff <- add_tag(build_dep_effect_panel(), "E")
pw <- add_tag(build_pathway_bar_panel()$plot, "F")

# Layout: square 178×178 like YvO F01.
# Top row: PCA · Venn · direction strip. Bottom row: DEP counts · effect size · pathways.
design <- "
AABBCC
DDEEFF
"

fig <- pca + venn_gg + strip + dep + eff + pw +
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
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = ggplot2::element_text(
        face = "bold", size = 7, hjust = 0, margin = ggplot2::margin(l = 0, b = 1)
      ),
      plot.subtitle = ggplot2::element_text(
        face = "italic", size = 5, color = "grey30",
        hjust = 0, margin = ggplot2::margin(l = 0)
      ),
      plot.caption = ggplot2::element_text(
        size = 5, color = "grey35", hjust = 0, lineheight = 1.1
      )
    )
  )

save_composite(fig, BASE, "MAIN_C1_overview", width_mm = PANEL_MD, height_mm = PANEL_MD)
message("C1 composite built")
