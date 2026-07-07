# Panel C: is the mitochondrial gain the same with and without stress? Plot the
# transplant effect at baseline (Mito - Ctl) against the transplant effect under
# stress (Phe+Mito - Phe) for mito proteins. Pure cargo delivery adds the same
# mass either way (slope 1); a context-sensitive response bends off the diagonal.

build_context_dependence <- function(comb, mito_genes) {
  d <- comb |>
    dplyr::filter(toupper(gene) %in% mito_genes) |>
    dplyr::transmute(gene, baseline = logFC_CTLvMITO, stressed = logFC_PHEvPHE_MITO)

  slope <- unname(stats::coef(stats::lm(stressed ~ baseline, d))[2])
  r <- stats::cor(d$baseline, d$stressed, use = "complete.obs")

  p <- ggplot2::ggplot(d, ggplot2::aes(baseline, stressed)) +
    ggplot2::geom_abline(slope = 1, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = 0, color = "grey85", linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = 0, color = "grey85", linewidth = 0.2) +
    ggplot2::geom_point(color = "#009988", size = 0.6, alpha = 0.6) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#B2182B", linewidth = 0.4) +
    ggplot2::labs(
      title = "Context dependence",
      subtitle = sprintf("mito gain: slope %.2f, r %.2f (slope 1 = stress-independent delivery)", slope, r),
      x = "Transplant at baseline  (Mito − Ctl)",
      y = "Transplant under stress  (Phe+Mito − Phe)"
    ) +
    FIG_THEME

  list(plot = p, table = dplyr::arrange(d, dplyr::desc(abs(baseline - stressed))))
}
