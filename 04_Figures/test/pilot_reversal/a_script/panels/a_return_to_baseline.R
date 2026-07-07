# Panel A: do proteins moved by PHE return toward the control baseline after
# transplant? Disease sets the signature (Phe - Ctl); Recovery (Phe+Mito - Ctl)
# measures how much of it remains. The two contrasts share Ctl, not PHE, so this
# avoids the shared-PHE coupling that biases a direct Disease-vs-Rescue scatter.

build_return_to_baseline <- function(comb) {
  d <- comb |>
    dplyr::filter(sig_pi_CTLvPHE != 0) |>
    dplyr::transmute(
      uniprot_id, gene,
      disease = logFC_CTLvPHE,
      recovery = logFC_CTLvPHE_MITO,
      attenuation = abs(disease) - abs(recovery),
      returned = attenuation > 0
    )
  ret_frac <- mean(d$returned)
  residual <- median(abs(d$recovery)) / median(abs(d$disease))

  p <- ggplot2::ggplot(d, ggplot2::aes(disease, recovery, color = returned)) +
    ggplot2::geom_abline(slope = 1, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = 0, color = H9C2_PAL_GROUP[["Ctl"]], linewidth = 0.3) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6) +
    ggplot2::scale_color_manual(
      values = c(`TRUE` = "#4DAF4A", `FALSE` = "grey60"),
      labels = c(`TRUE` = "toward baseline", `FALSE` = "no return"), name = NULL
    ) +
    ggplot2::labs(
      title = "Return to baseline",
      subtitle = sprintf(
        "%.0f%% of disease proteins move back; residual %.2f (0 = full reversal)",
        100 * ret_frac, residual
      ),
      x = "Disease  (Phe − Ctl) log2FC",
      y = "Recovery  (Phe+Mito − Ctl) log2FC"
    ) +
    FIG_THEME +
    ggplot2::theme(legend.position = "bottom")

  list(plot = p, table = dplyr::arrange(d, dplyr::desc(attenuation)))
}
