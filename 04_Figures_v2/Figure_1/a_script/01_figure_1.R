#!/usr/bin/env Rscript
# Figure 1 — overview. Builds each panel once, writes the standalone panels to
# b_reports/panels/, and the assembled composite to the b_reports root.
library(here)
source(here::here("04_Figures_v2", "functions", "08_composite_layout.R"))
source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures_v2", "functions", "04_mitocarta_lens_lookup.R"))
source(here::here("04_Figures_v2", "01_PCA", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "02_DEP_bars", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "03_Venn", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "04_Pathway_bars", "a_script", "_build.R"))

BASE <- here::here("04_Figures_v2", "Figure_1")
RPT <- file.path(BASE, "b_reports")
PANELS <- file.path(RPT, "panels")
for (d in c(RPT, PANELS, file.path(BASE, "c_data"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
pdf_dev <- get_pdf_device()

# ---- build each panel once --------------------------------------------------
pca_p <- build_pca_panel()$plot
dep_p <- build_dep_count_panel()
eff_p <- build_dep_effect_panel()
venn_l <- build_venn_panels()
venn_p <- venn_l$venn
strip_p <- venn_l$strip
pw_p <- build_pathway_bar_panel()$plot

# ---- standalone panel PNGs (each with its own tight axes) -------------------
panel_specs <- list(
  list(name = "panel_a_pca", plot = pca_p, w = 95, h = 80),
  list(name = "panel_b_dep_counts", plot = dep_p, w = 95, h = 72),
  list(name = "panel_c_effect_size", plot = eff_p, w = 72, h = 88),
  list(name = "panel_d_venn", plot = venn_p, w = 88, h = 80),
  list(name = "panel_e_direction", plot = strip_p, w = 82, h = 76),
  list(name = "panel_f_pathway", plot = pw_p, w = 92, h = 80)
)
for (s in panel_specs) {
  ggplot2::ggsave(
    file.path(PANELS, paste0(s$name, ".png")), s$plot,
    width = s$w, height = s$h, units = "mm", dpi = 300
  )
}
message(sprintf("Figure 1: wrote %d standalone panels", length(panel_specs)))

# ---- composite (letters baked into titles via add_tag) ----------------------
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

pca <- add_tag(pca_p, "A")
dep <- add_tag(dep_p, "B") +
  patchwork::inset_element(
    dep_key,
    left = 0.32, right = 0.92, top = 0.13, bottom = 0.0
  )
eff <- add_tag(eff_p, "C")
venn_gg <- add_tag(venn_p, "D") +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 2, 1, 2))
# free(type = "label") lets E's y-title sit tight to its axis while the panel
# stays aligned to panel B's width above it.
strip <- patchwork::free(
  add_tag(strip_p, "E") + ggplot2::theme(plot.margin = ggplot2::margin(0, 2, 1, 2)),
  type = "label", side = "l"
)
pw <- add_tag(pw_p, "F") +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 4, 1, 2))

# Top row: PCA · DEP counts · effect size. Bottom row: Venn · direction · pathways.
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

COMP_W <- PANEL_MD
COMP_H <- 150
ggplot2::ggsave(
  file.path(RPT, "MAIN_figure_1.png"), fig,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300, limitsize = FALSE
)
pdf_path <- file.path(RPT, "MAIN_figure_1.pdf")
pdf_ok <- tryCatch(
  {
    ggplot2::ggsave(pdf_path, fig,
      width = COMP_W, height = COMP_H, units = "mm",
      device = pdf_dev, limitsize = FALSE
    )
    file.exists(pdf_path) && file.info(pdf_path)$size > 1000
  },
  error = function(e) FALSE
)
if (!pdf_ok) {
  message("preferred pdf device failed — falling back to base grDevices::pdf()")
  grDevices::pdf(pdf_path, width = COMP_W / 25.4, height = COMP_H / 25.4)
  tryCatch(print(fig), finally = grDevices::dev.off())
}
message("Figure 1 composite built")
