# Panel A: does transplant move only mitochondrial proteins (cargo delivery) or
# also the wider proteome (remodeling)? Compare log2FC of MitoCarta vs other
# proteins across Disease, Transplant, and Rescue.

build_mito_vs_nonmito <- function(comb, mito_genes) {
  d <- comb |>
    dplyr::transmute(
      gene,
      lens = ifelse(toupper(gene) %in% mito_genes, "MitoCarta", "other"),
      Disease = logFC_CTLvPHE,
      Transplant = logFC_CTLvMITO,
      Rescue = logFC_PHEvPHE_MITO
    ) |>
    tidyr::pivot_longer(c(Disease, Transplant, Rescue), names_to = "contrast", values_to = "logFC") |>
    dplyr::mutate(contrast = factor(contrast, c("Disease", "Transplant", "Rescue")))

  p <- ggplot2::ggplot(d, ggplot2::aes(contrast, logFC, fill = lens)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3) +
    ggplot2::geom_boxplot(outlier.size = 0.15, linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = c(MitoCarta = "#009988", other = "grey80"), name = NULL) +
    ggplot2::labs(
      title = "Mito vs non-mito",
      subtitle = "cargo delivery moves only mito; remodeling also shifts the rest",
      x = NULL, y = "log2FC"
    ) +
    FIG_THEME +
    ggplot2::theme(legend.position = "bottom")

  tab <- d |>
    dplyr::group_by(contrast, lens) |>
    dplyr::summarise(median_logFC = median(logFC, na.rm = TRUE), n = dplyr::n(), .groups = "drop")
  list(plot = p, table = tab)
}
