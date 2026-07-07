# Panel C: the honest cap. The formal 2x2 interaction is the strict per-protein
# test of transplant x stress. If it carries few FDR hits, the reversal is a
# coordinated shift toward baseline, not a significant protein-by-protein rescue.

build_interaction_null <- function(comb, fdr = H9C2_FDR_EXPLOR) {
  d <- comb |>
    dplyr::transmute(
      uniprot_id, gene,
      logFC = logFC_Interaction,
      p = P.Value_Interaction,
      fdr = adj.P.Val_Interaction,
      hit = adj.P.Val_Interaction < fdr
    )
  n_hit <- sum(d$hit, na.rm = TRUE)

  p <- ggplot2::ggplot(d, ggplot2::aes(logFC, -log10(p), color = hit)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.5) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#D6604D", `FALSE` = "grey70"), guide = "none") +
    ggplot2::labs(
      title = "Interaction is null",
      subtitle = sprintf("%d proteins at FDR < %.2f: a coordinated shift, not per-protein rescue", n_hit, fdr),
      x = "Interaction log2FC  (Phe+Mito − Phe) − (Mito − Ctl)", y = "−log10 p"
    ) +
    FIG_THEME

  list(plot = p, table = dplyr::arrange(d, p))
}
