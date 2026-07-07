# Panel B: which mitochondrial programs move, and how coordinately. Each
# MitoPathway is tested with fry (self-contained rotation test, valid at n=6);
# the signed statistic is the fry direction times -log10 FDR, so a point far from
# zero is a program that moves strongly and coherently. Structural machinery
# (OXPHOS, mtDNA/translation) is split from the regulatory programs; a regulated
# response moves the regulatory arm, not only the structural cargo.

build_program_breakdown <- function(cat_tests) {
  d <- cat_tests |>
    dplyr::filter(contrast %in% c("Transplant", "Rescue")) |>
    dplyr::mutate(
      signed = ifelse(fry_dir == "Up", 1, -1) * -log10(fry_fdr),
      arm = factor(arm, c("regulatory", "structural")),
      sig = fry_fdr < 0.05
    )

  p <- ggplot2::ggplot(d, ggplot2::aes(signed, reorder(label, signed), color = contrast, shape = sig)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::facet_grid(arm ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_color_manual(values = c(Transplant = "#4393C3", Rescue = "#4DAF4A"), name = NULL) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    ggplot2::labs(
      title = "Which programs move (fry)",
      subtitle = "self-contained test; filled = FDR < 0.05, right = up",
      x = "signed −log10 fry FDR", y = NULL
    ) +
    FIG_THEME +
    ggplot2::theme(legend.position = "bottom")

  list(plot = p, table = dplyr::arrange(cat_tests, contrast, set))
}
