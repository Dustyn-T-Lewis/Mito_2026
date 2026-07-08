# Panel A: the pooled transplant main effect per protein. Colour marks MitoCarta
# membership; the Fisher test asks whether the pi-score movers are enriched for
# delivered cargo over host proteins, which is what separates "added protein"
# from "changed the host".

build_mito_main_volcano <- function(main, pi_thresh = H9C2_PI_THRESH) {
  n_hit <- sum(main$hit, na.rm = TRUE)
  n_mito_hit <- sum(main$hit & main$lens == "MitoCarta", na.rm = TRUE)
  ft <- fisher.test(table(
    factor(main$hit, c(FALSE, TRUE)),
    factor(main$lens, c("other", "MitoCarta"))
  ))

  p <- ggplot2::ggplot(main, ggplot2::aes(logFC, -log10(P.Value))) +
    ggplot2::geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(color = lens, alpha = hit, size = hit)) +
    ggplot2::scale_color_manual(
      values = c(MitoCarta = "#009988", other = "grey70"), name = NULL
    ) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.25), guide = "none") +
    ggplot2::scale_size_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.4), guide = "none") +
    ggplot2::labs(
      title = "Transplant main effect",
      subtitle = sprintf(
        "%d movers at Π < %.2f, %d MitoCarta; cargo enrichment OR %.1f, %s",
        n_hit, pi_thresh, n_mito_hit, ft$estimate, fmt_p(ft$p.value)
      ),
      x = "Main-effect log2FC  (Mito, Phe+Mito) − (Ctl, Phe)",
      y = "−log10 p"
    ) +
    FIG_THEME +
    ggplot2::theme(legend.position = "bottom")

  list(plot = p, table = dplyr::arrange(main, pi_score))
}
