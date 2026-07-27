# F01 panel C: effect-size companion (contrast-colored). Builds from the shared
# dep_data() object (dep.R) the composite passes in, so nothing loads at source time.

build_dep_effect_panel <- function(dep) {
  dep_results <- dep$results
  CORE <- dep$core
  CTR_LAB <- dep$ctr_lab
  ctr_levels <- dep$ctr_levels

  # All logFC (no window filter) for the median stat
  lfc_long_all <- dplyr::bind_rows(lapply(CORE, \(c)
  tibble::tibble(contrast = CTR_LAB[c], logFC = dep_results[[c]]$logFC))) |>
    dplyr::filter(!is.na(logFC)) |>
    dplyr::mutate(contrast = factor(contrast, levels = ctr_levels))

  # Median over ALL proteins; label clarifies scope
  lfc_stats <- lfc_long_all |>
    dplyr::summarise(med_abs = median(abs(logFC), na.rm = TRUE), .by = contrast) |>
    dplyr::mutate(lab = sprintf("median |log2FC| %.2f", med_abs))

  # Histogram restricted to ±0.8 window
  lfc_plot <- dplyr::filter(lfc_long_all, abs(logFC) <= 0.8)

  hbw <- 1.6 / 40

  ggplot2::ggplot(lfc_plot, ggplot2::aes(logFC)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, color = "grey55") +
    ggplot2::geom_histogram(
      ggplot2::aes(fill = contrast),
      breaks = seq(-0.8, 0.8, by = hbw),
      color = "white",
      linewidth = 0.1,
      alpha = 0.85
    ) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(count) * hbw),
      color = "grey20", linewidth = 0.4
    ) +
    ggplot2::geom_label(
      data = lfc_stats,
      ggplot2::aes(x = -0.78, y = Inf, label = lab),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1.0,
      size = scale_text(BASE_STAT, 50) + 0.3,
      fontface = "bold", color = "grey15",
      fill = scales::alpha("white", 0.8), linewidth = 0,
      label.padding = ggplot2::unit(0.3, "mm")
    ) +
    ggplot2::facet_wrap(~contrast, ncol = 1, scales = "free_y", strip.position = "top") +
    ggplot2::scale_fill_manual(
      values = stats::setNames(unname(CONTRAST_COLORS[CORE]), unname(CTR_LAB)),
      guide  = "none"
    ) +
    ggplot2::scale_x_continuous(breaks = c(-0.8, 0, 0.8)) +
    ggplot2::scale_y_continuous(breaks = NULL) +
    ggplot2::labs(
      title = "Effect size",
      subtitle = "log2FC per contrast (median over all proteins)",
      x = expression(bold(log[2] ~ FC)),
      y = NULL
    ) +
    FIG_THEME +
    ggplot2::theme(
      strip.text          = ggplot2::element_blank(),
      strip.background    = ggplot2::element_blank(),
      axis.text.y         = ggplot2::element_blank(),
      axis.ticks.y        = ggplot2::element_blank(),
      axis.text.x         = ggplot2::element_text(size = FIG_AXIS_TEXT),
      panel.grid          = ggplot2::element_blank(),
      panel.spacing.y     = ggplot2::unit(1.5, "pt"),
      plot.margin         = ggplot2::margin(5, 2, 1, 0)
    )
}
