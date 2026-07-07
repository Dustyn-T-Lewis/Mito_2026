# Panel A: do proteins moved by PHE return toward the control baseline after
# transplant? Disease sets the signature (Phe - Ctl, Pi < 0.05); Recovery
# (Phe+Mito - Ctl) measures what remains. Disease and Recovery share Ctl, not
# PHE, so this avoids the shared-PHE coupling that biases a Disease-vs-Rescue
# scatter. Each disease protein is classed by its Recovery Pi-score: normalized
# (no longer significant vs Ctl), persistent (still off in the same direction),
# or overshoot (crossed past baseline).

build_return_to_baseline <- function(comb, pi_thresh = H9C2_PI_THRESH) {
  d <- comb |>
    dplyr::filter(sig_pi_CTLvPHE != 0) |>
    dplyr::transmute(
      uniprot_id, gene,
      disease = logFC_CTLvPHE,
      recovery = logFC_CTLvPHE_MITO,
      recovery_pi = pi_score_CTLvPHE_MITO,
      class = dplyr::case_when(
        pi_score_CTLvPHE_MITO >= pi_thresh ~ "normalized",
        sign(logFC_CTLvPHE_MITO) == sign(logFC_CTLvPHE) ~ "persistent",
        TRUE ~ "overshoot"
      )
    )
  d$class <- factor(d$class, c("normalized", "persistent", "overshoot"))
  norm_frac <- mean(d$class == "normalized")
  residual <- median(abs(d$recovery)) / median(abs(d$disease))

  p <- ggplot2::ggplot(d, ggplot2::aes(disease, recovery, color = class)) +
    ggplot2::geom_abline(slope = 1, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = 0, color = H9C2_PAL_GROUP[["Ctl"]], linewidth = 0.3) +
    ggplot2::geom_point(size = 0.8, alpha = 0.8) +
    ggplot2::scale_color_manual(
      values = c(normalized = "#4DAF4A", persistent = "#E08214", overshoot = "#8073AC"),
      drop = FALSE, name = NULL
    ) +
    ggplot2::labs(
      title = "Return to baseline",
      subtitle = sprintf(
        "%.0f%% of disease proteins normalized (Π ≥ %.2f vs Ctl); residual %.2f",
        100 * norm_frac, pi_thresh, residual
      ),
      x = "Disease  (Phe − Ctl) log2FC",
      y = "Recovery  (Phe+Mito − Ctl) log2FC"
    ) +
    FIG_THEME +
    ggplot2::theme(legend.position = "bottom")

  list(plot = p, table = dplyr::arrange(d, recovery_pi))
}
