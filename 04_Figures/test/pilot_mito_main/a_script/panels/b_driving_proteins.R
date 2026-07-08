# Panel B: which proteins drive the pooled effect. Rank the pi-score movers by
# moderated |t|; the lollipop shows signed main-effect log2FC and whether each
# driver is delivered MitoCarta cargo or a host protein. A cargo-dominated bar is
# the added-protein signature; host proteins among the drivers are the transplant
# reaching past its own cargo.

build_driving_proteins <- function(main, top_n = 25L) {
  d <- main |>
    dplyr::filter(hit) |>
    dplyr::slice_max(abs(t), n = top_n) |>
    dplyr::mutate(
      gene = dplyr::if_else(is.na(gene) | gene == "", uniprot_id, gene),
      gene = reorder(gene, logFC)
    )

  p <- ggplot2::ggplot(d, ggplot2::aes(logFC, gene, color = lens)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = logFC, yend = gene), linewidth = 0.3) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_color_manual(
      values = c(MitoCarta = "#009988", other = "grey60"), name = NULL
    ) +
    ggplot2::labs(
      title = sprintf("Top %d drivers", nrow(d)),
      subtitle = "pi-score movers ranked by |moderated t|; right = up with transplant",
      x = "main-effect log2FC", y = NULL
    ) +
    FIG_THEME +
    ggplot2::theme(legend.position = "bottom")

  tab <- main |>
    dplyr::filter(hit) |>
    dplyr::arrange(dplyr::desc(abs(t)))
  list(plot = p, table = tab)
}
