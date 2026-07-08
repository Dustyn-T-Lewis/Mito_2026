#!/usr/bin/env Rscript
# F01 proteome overview — PCA + DEP bars + effect histogram + Venn + pathway bars.
# Builds each panel once: saves them standalone under b_reports/panels, assembles
# the 6-panel composite, and writes one supplementary workbook holding the data
# behind every panel (PCA scores, DE tables, overlap membership, pathway counts).

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
})

source(here::here("04_Figures", "functions", "08_composite_layout.R"))
source(here::here("04_Figures", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures", "functions", "04_mitocarta_lens_lookup.R"))
source(here::here("04_Figures", "functions", "06_supplementary_workbook.R"))
panels <- here::here("04_Figures", "F01_Proteome_Overview", "a_script", "panels")
for (f in list.files(panels, full.names = TRUE)) source(f)

BASE <- here::here("04_Figures", "F01_Proteome_Overview")
PANELS <- file.path(BASE, "b_reports", "panels")
DAT <- file.path(BASE, "c_data")
for (d in c(PANELS, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Base panels + their underlying data, built once.
pca_res <- build_pca_panel()
pca_p <- pca_res$plot
dep_p <- build_dep_count_panel()
eff_p <- build_dep_effect_panel()
venn <- build_venn_panels()
pw_out <- build_pathway_bar_panel()
pw_p <- pw_out$plot

# Standalone named panels, each with its own tight axes.
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

# B: DEP counts + a stringency key — square swatch then the bold term, stacked
# dark -> light (Π / FDR / p) top to bottom to read with the bars, bottom-right.
dep_key_df <- tibble::tibble(
  y = c(3, 2, 1),
  a = c(1, 0.45, 0.18),
  lab = c(
    paste0("Π < ", H9C2_PI_THRESH),
    "FDR < 0.10",
    "p < 0.05"
  )
)
dep_key <- ggplot2::ggplot(dep_key_df, ggplot2::aes(0, y)) +
  ggplot2::geom_point(
    ggplot2::aes(alpha = a),
    shape = 22, size = 3, fill = "grey40", color = "black", stroke = 0.4
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = 0.18, label = lab),
    hjust = 0, size = 1.6, fontface = "bold", color = "grey15"
  ) +
  ggplot2::scale_alpha_identity() +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = c(0.15, 1.4))) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = 1.1)) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_void()
dep <- add_tag(dep_p, "B") +
  patchwork::inset_element(
    dep_key,
    left = 0.7, right = 1.0, top = 0.24, bottom = 0.0
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

composite <- save_composite(fig, BASE, "MAIN_F01_proteome_overview", width_mm = PANEL_MD, height_mm = 150)
mirror_to_box(composite, "02_Figures/F01_Proteome_Overview")

# Supplementary workbook — the data behind each panel, one sheet group per panel.
CORE <- H9C2_CONTRAST_ORDER
dep_tabs <- lapply(CORE, function(ctr) {
  dep_results[[ctr]] |>
    transmute(uniprot_id, gene, logFC, P.Value, adj.P.Val, pi_score) |>
    arrange(pi_score)
})
pw_tabs <- lapply(CORE, function(ctr) {
  pw_out$sig_pw |>
    filter(contrast == ctr) |>
    transmute(pathway, database,
      clean_label = clean_display_label(pathway),
      NES = round(NES, 3), padj = signif(padj, 3), size, direction, is_mito
    ) |>
    arrange(padj)
})
pw_counts <- pw_out$bar_df |>
  filter(n > 0) |>
  transmute(
    contrast = contrast_brief(as.character(contrast)),
    direction, database = as.character(database), n
  ) |>
  arrange(contrast, direction, desc(n))

sheet_specs <- c(
  list(
    list(
      name = "pca_scores", df = pca_res$scores,
      role = "PCA scatter coordinates (panel A)",
      contents = "PC1/PC2 scores per sample (Col_ID) with Group"
    ),
    list(
      name = "permanova", df = pca_res$permanova,
      role = "Panel A stats annotation",
      contents = "adonis2 R2/p for overall Group + Disease/Transplant/Rescue pairwise, plus betadisper p"
    ),
    list(
      name = "dep_counts", df = dep_count_data() |> select(contrast, threshold, n, pct),
      role = "DEP count bar heights (panel B)",
      contents = "Counts and % of proteome at each threshold (p<0.05, FDR, Π<0.05) per contrast"
    )
  ),
  Map(function(ctr, tab) {
    list(
      name = paste0(contrast_brief(ctr), "_DE"), df = tab,
      role = sprintf("Per-protein DE table for the %s contrast", contrast_brief(ctr)),
      contents = "uniprot_id, gene, logFC, P.Value, adj.P.Val, pi_score (sorted by Π)"
    )
  }, CORE, dep_tabs),
  list(
    list(
      name = "venn_membership", df = venn$membership,
      role = "Per-protein set membership (panels D/E)",
      contents = "One row per Π<0.05 protein: membership + direction in Disease/Transplant/Rescue, region_key"
    ),
    list(
      name = "venn_regions", df = venn$region_counts,
      role = "Panel D Euler region areas + directional strip counts",
      contents = "Protein count per overlap region (region_key), sorted descending"
    ),
    list(
      name = "euler_fit", df = venn$fit_stats,
      role = "Panel D Euler fit quality",
      contents = "stress and diagError from the eulerr fit"
    ),
    list(
      name = "pathway_counts", df = pw_counts,
      role = "Pathway count bar heights (panel F)",
      contents = "contrast, direction, database, n (5-DB lens, padj<0.05, all sets, no dedup)"
    )
  ),
  Map(function(ctr, tab) {
    list(
      name = paste0(contrast_brief(ctr), "_pathways"), df = tab,
      role = sprintf("All significant pathways behind the %s bars", contrast_brief(ctr)),
      contents = "pathway, database, clean_label, NES, padj, size, direction, is_mito"
    )
  }, CORE, pw_tabs)
)

build_workbook(
  file.path(DAT, "F01_supplementary.xlsx"),
  figure_title = "F01: Proteome overview (PCA, DEP counts, effect size, overlap, pathway counts)",
  sheet_specs = sheet_specs
)

message("F01 proteome overview built")
