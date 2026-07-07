# Panel A: is the mitochondrial set shifted relative to the rest of the proteome?
# Boxplots show the reported log2FC of MitoCarta vs other proteins; the annotation
# is the competitive camera test (Wu & Smyth 2012) on the whole MitoCarta set,
# which asks the same question with inter-protein correlation controlled.

build_mito_vs_nonmito <- function(comb, mito_genes, all_tests) {
  keep <- c("Transplant", "Rescue")
  d <- comb |>
    dplyr::transmute(
      gene,
      lens = ifelse(toupper(gene) %in% mito_genes, "MitoCarta", "other"),
      Transplant = logFC_CTLvMITO,
      Rescue = logFC_PHEvPHE_MITO
    ) |>
    tidyr::pivot_longer(dplyr::all_of(keep), names_to = "contrast", values_to = "logFC") |>
    dplyr::mutate(contrast = factor(contrast, keep))

  cam <- all_tests |>
    dplyr::filter(contrast %in% keep) |>
    dplyr::mutate(
      contrast = factor(contrast, keep),
      lab = sprintf(
        "mito %s\ncamera FDR %s", tolower(camera_dir),
        ifelse(camera_fdr < 0.001, "<0.001", sprintf("%.3f", camera_fdr))
      )
    )

  p <- ggplot2::ggplot(d, ggplot2::aes(lens, logFC, fill = lens)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3) +
    ggplot2::geom_boxplot(outlier.size = 0.15, linewidth = 0.25) +
    ggplot2::geom_text(
      data = cam, ggplot2::aes(x = 1.5, y = Inf, label = lab),
      inherit.aes = FALSE, vjust = 1.1, size = 1.6, color = "grey20", lineheight = 0.9
    ) +
    ggplot2::facet_wrap(~contrast) +
    ggplot2::scale_fill_manual(values = c(MitoCarta = "#009988", other = "grey80"), guide = "none") +
    ggplot2::labs(
      title = "Mito set vs the proteome",
      subtitle = "MitoCarta shifted relative to the rest (competitive camera)",
      x = NULL, y = "log2FC"
    ) +
    FIG_THEME

  tab <- d |>
    dplyr::group_by(contrast, lens) |>
    dplyr::summarise(median_logFC = median(logFC, na.rm = TRUE), n = dplyr::n(), .groups = "drop") |>
    dplyr::left_join(dplyr::select(cam, contrast, camera_dir, camera_fdr), by = "contrast")
  list(plot = p, table = tab)
}
