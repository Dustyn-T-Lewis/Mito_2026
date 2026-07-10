# Panel B: module preservation against the control reference. Zsummary vs module size,
# faceted by condition. A module that drops under PHE and climbs back under PHE+MitoTx is
# the rescue signature. Z=10 (strong) and Z=2 (weak) guides; size is shown because small
# modules score low by power alone (Langfelder 2011). Featured modules are labelled.
pacman::p_load(ggplot2, dplyr)

panel_preservation <- function(preservation, featured) {
  pres <- preservation |>
    mutate(is_featured = module %in% featured)
  ggplot(pres, aes(mod_size, Zsummary)) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    geom_hline(yintercept = 2, linetype = "dotted", color = "grey55", linewidth = 0.3) +
    geom_point(aes(fill = I(module), size = is_featured), shape = 21, color = "grey25", stroke = 0.3) +
    ggrepel::geom_text_repel(
      data = filter(pres, is_featured), aes(label = module),
      size = 1.7, color = "grey15", fontface = "bold", max.overlaps = Inf, segment.size = 0.12
    ) +
    facet_wrap(~condition, nrow = 1) +
    scale_size_manual(values = c(`FALSE` = 1.7, `TRUE` = 2.8), guide = "none") +
    scale_x_log10() +
    labs(
      title = "Module preservation vs control",
      x = "module size (proteins, log scale)", y = "Zsummary"
    ) +
    FIG_THEME +
    theme(
      strip.text = element_text(face = "bold", size = FIG_AXIS_TEXT + 1),
      panel.spacing = unit(3, "mm")
    )
}
