# Panel B: proteome geometry. Sample scores and group centroids in PC space, with
# arrows from the PHE and Phe+Mito centroids to the control centroid. A shorter
# arrow means closer to baseline; the annotated ratio and permutation p are the
# full-space distances, not the 2D projection.

build_pca_trajectory <- function(scores, centroids, var_exp, dist_stat) {
  ctl <- centroids[centroids$group == "Ctl", ]
  arrows <- centroids[centroids$group %in% c("PHE", "PHE_Mito"), ]
  arrows$xend <- ctl$PC1
  arrows$yend <- ctl$PC2

  p <- ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2, color = group)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.55) +
    ggplot2::geom_segment(
      data = arrows, ggplot2::aes(xend = xend, yend = yend),
      arrow = grid::arrow(length = grid::unit(1.5, "mm")), linewidth = 0.4, color = "grey35"
    ) +
    ggplot2::geom_point(data = centroids, size = 2.6, shape = 21, fill = "white", stroke = 0.7) +
    ggplot2::scale_color_manual(values = GROUP_COLORS, labels = GROUP_LABELS, name = NULL) +
    ggplot2::labs(
      title = "Proteome trajectory to control",
      subtitle = sprintf(
        "d(Phe+Mito,Ctl)/d(Phe,Ctl) = %.2f, %s (full space)",
        dist_stat$ratio, fmt_p(dist_stat$p)
      ),
      x = sprintf("PC1 (%.1f%%)", var_exp[1]), y = sprintf("PC2 (%.1f%%)", var_exp[2])
    ) +
    FIG_THEME +
    ggplot2::theme(legend.position = "bottom")

  tab <- tibble::tibble(
    statistic = c("d(Phe,Ctl)", "d(Phe+Mito,Ctl)", "ratio", "perm_p"),
    value = c(dist_stat$d_phe, dist_stat$d_res, dist_stat$ratio, dist_stat$p)
  )
  list(plot = p, table = tab)
}
