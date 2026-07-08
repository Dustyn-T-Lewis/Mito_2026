# Panel C: which disease proteins return, with a key. Each disease-signature
# protein is plotted as its disease effect (Phe - Ctl) against its residual after
# transplant (Phe+Mito - Ctl); the diagonal is no rescue, the y = 0 line is full
# return to baseline. Colour marks MitoCarta membership and shape marks the
# disease direction, so the reversal can be read against mito identity.

build_return_key <- function(prot) {
  n_mito <- sum(prot$lens == "MitoCarta")

  p <- ggplot2::ggplot(prot, ggplot2::aes(disease, recovery)) +
    ggplot2::geom_abline(slope = 1, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = 0, color = H9C2_PAL_GROUP[["Ctl"]], linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(color = lens, shape = direction), size = 1, alpha = 0.85) +
    ggplot2::scale_color_manual(values = c(MitoCarta = "#009988", other = "grey55"), name = NULL) +
    ggplot2::scale_shape_manual(values = c(`disease-up` = 17, `disease-down` = 16), name = NULL) +
    ggplot2::labs(
      title = "Which proteins return",
      subtitle = sprintf(
        "%d of %d disease proteins are MitoCarta; diagonal = no rescue, y = 0 = baseline",
        n_mito, nrow(prot)
      ),
      x = "Disease  (Phe − Ctl) log2FC", y = "Recovery  (Phe+Mito − Ctl) log2FC"
    ) +
    FIG_THEME +
    ggplot2::theme(legend.position = "bottom")

  list(plot = p, table = dplyr::arrange(prot, recovery))
}
