#!/usr/bin/env Rscript
# F01 proteome overview — PCA + DEP bars + effect histogram + Venn + pathway bars.
# Saves each panel standalone under b_reports/panels, then the 6-panel composite.
library(here)
source(here::here("04_Figures_v2", "functions", "08_composite_layout.R"))
source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures_v2", "functions", "04_mitocarta_lens_lookup.R"))
source(here::here("04_Figures_v2", "01_PCA", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "02_DEP_bars", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "03_Venn", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "04_Pathway_bars", "a_script", "_build.R"))

BASE <- here::here("04_Figures_v2", "F01_Proteome_Overview")

# Base panels, built once — saved standalone below, then tagged for the composite.
pca_p <- build_pca_panel()$plot
dep_p <- build_dep_count_panel()
eff_p <- build_dep_effect_panel()
venn <- build_venn_panels()
pw_p <- build_pathway_bar_panel()$plot

# Standalone named panels, each with its own tight axes.
PANELS <- file.path(BASE, "b_reports", "panels")
dir.create(PANELS, recursive = TRUE, showWarnings = FALSE)
panel_specs <- list(
  list(name = "panel_a_pca", plot = pca_p, w = 95, h = 80),
  list(name = "panel_b_dep_counts", plot = dep_p, w = 95, h = 72),
  list(name = "panel_c_effect_size", plot = eff_p, w = 72, h = 88),
  list(name = "panel_d_venn", plot = venn$venn, w = 88, h = 80),
  list(name = "panel_e_direction", plot = venn$strip, w = 82, h = 76),
  list(name = "panel_f_pathway", plot = pw_p, w = 92, h = 80)
)
for (s in panel_specs) {
  ggplot2::ggsave(file.path(PANELS, paste0(s$name, ".png")), s$plot,
    width = s$w, height = s$h, units = "mm", dpi = 300
  )
}

# Composite — add_tag bakes the letter into each title for uniform spacing.
pca <- add_tag(pca_p, "A")

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
    fill = "grey40", size = 1.8, fontface = "bold",
    label.padding = ggplot2::unit(0.5, "mm"), label.size = 0.25
  ) +
  ggplot2::scale_alpha_identity() +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.4)) +
  ggplot2::theme_void()
dep <- add_tag(dep_p, "B") +
  patchwork::inset_element(
    dep_key,
    left = 0.32, right = 0.92, top = 0.13, bottom = 0.0
  )

eff <- add_tag(eff_p, "C")

# Strip shares column 2 with panel B, so patchwork aligns it to the same plot
# width and left edge automatically. Bottom-row panels get a 0 top margin to
# close the gap to the row above.
venn_gg <- add_tag(venn$venn, "D") +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 2, 1, 2))
# free(type = "label") lets E's y-title sit tight to its axis while the panel
# stays aligned to panel B's width above it.
strip <- patchwork::free(
  add_tag(venn$strip, "E") + ggplot2::theme(plot.margin = ggplot2::margin(0, 2, 1, 2)),
  type = "label", side = "l"
)

pw <- add_tag(pw_p, "F") +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 4, 1, 2))

# Layout: square 178×178 like YvO F01.
# Top row: PCA · DEP counts · effect size. Bottom row: Venn · direction strip · pathways.
design <- "
AABBCC
DDEEFF
"

fig <- pca + dep + eff + venn_gg + strip + pw +
  patchwork::plot_layout(
    design  = design,
    widths  = c(1.2, 1.02, 0.78),
    heights = c(1.0, 1.0)
  )

save_composite(fig, BASE, "MAIN_F01_proteome_overview", width_mm = PANEL_MD, height_mm = 150)
message("F01 proteome overview built")
