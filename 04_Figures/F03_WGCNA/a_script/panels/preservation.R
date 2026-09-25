# Supplementary: module preservation of the control-reference modules in each treatment
# group (WGCNA modulePreservation, from 00_build_wgcna.R). Zsummary < 2 is no evidence
# of preservation and > 10 is strong (Langfelder et al. 2011); the dashed lines mark both.
pacman::p_load(ggplot2, dplyr, ggrepel)

panel_preservation <- function(preservation) {
  d <- preservation |>
    mutate(condition = factor(H9C2_GROUP_LABELS[as.character(condition)],
      levels = H9C2_GROUP_LABELS[levels(condition)]
    ))
  ggplot(d, aes(mod_size, Zsummary)) +
    geom_hline(yintercept = c(2, 10), linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    geom_point(aes(fill = module), shape = 21, size = 2, colour = "grey20", stroke = 0.2) +
    ggrepel::geom_text_repel(aes(label = module),
      size = 1.6, colour = "grey20", segment.size = 0.2, max.overlaps = Inf, seed = 42
    ) +
    scale_fill_identity() +
    scale_x_log10() +
    facet_wrap(~condition, nrow = 1) +
    labs(
      title = "Module preservation in each treatment group",
      subtitle = "Control-reference modules; dashed lines at Zsummary 2 and 10",
      x = "module size (proteins, log scale)", y = "Zsummary"
    ) +
    FIG_THEME
}
