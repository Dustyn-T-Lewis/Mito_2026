# Panel A: set-level reversal. Each set is tested with fry (self-contained,
# valid at n=6) on Disease (Phe - Ctl) and Recovery (Phe+Mito - Ctl); the segment
# runs from the Disease enrichment to the Recovery enrichment, so a set pulled
# back toward zero has reversed. Disease and Recovery share Ctl, not PHE. The
# disease sets are defined in-sample (Disease is a circular anchor there); the
# cardiac-stress Hallmark pathways are a priori, so both contrasts are honest.

build_signature_reversal <- function(set_rev) {
  wide <- set_rev |>
    dplyr::select(label, set_type, contrast, signed) |>
    tidyr::pivot_wider(names_from = contrast, values_from = signed)

  p <- ggplot2::ggplot(set_rev, ggplot2::aes(signed, label)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    ggplot2::geom_segment(
      data = wide, ggplot2::aes(x = Disease, xend = Recovery, y = label, yend = label),
      color = "grey70", linewidth = 0.4, inherit.aes = FALSE
    ) +
    ggplot2::geom_point(ggplot2::aes(color = contrast), size = 1.4) +
    ggplot2::facet_grid(set_type ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_color_manual(values = c(Disease = "#D6604D", Recovery = "#4DAF4A"), name = NULL) +
    ggplot2::labs(
      title = "Signature reversal (fry)",
      subtitle = "Disease → Recovery; toward zero = returned to control",
      x = "signed −log10 fry FDR", y = NULL
    ) +
    FIG_THEME +
    ggplot2::theme(legend.position = "bottom")

  list(plot = p, table = dplyr::arrange(set_rev, set_type, label, contrast))
}
