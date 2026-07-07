# Panel C: regulated, not bulk cargo. Bulk delivery adds the same mitochondrial
# proteins regardless of stress, so the MT x PHE interaction is null; a regulated
# response makes the mito programs behave differently under stress, so the
# interaction is significant. Each MitoPathway is tested on the interaction
# contrast with fry; programs right of the FDR line are context-dependent.

build_context_dependence <- function(cat_tests) {
  d <- cat_tests |>
    dplyr::filter(contrast == "Interaction") |>
    dplyr::mutate(sig = fry_fdr < 0.05)

  p <- ggplot2::ggplot(d, ggplot2::aes(-log10(fry_fdr), reorder(label, -log10(fry_fdr)), color = arm)) +
    ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey55", linewidth = 0.3) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_color_manual(
      values = c(structural = "#8073AC", regulatory = "#009988"), name = NULL
    ) +
    ggplot2::labs(
      title = "Stress-dependence (MT × PHE)",
      subtitle = "interaction per program (fry); right of the line = context-dependent, not additive",
      x = "−log10 fry FDR (interaction)", y = NULL
    ) +
    FIG_THEME +
    ggplot2::theme(legend.position = "bottom")

  list(plot = p, table = dplyr::arrange(d, fry_fdr))
}
