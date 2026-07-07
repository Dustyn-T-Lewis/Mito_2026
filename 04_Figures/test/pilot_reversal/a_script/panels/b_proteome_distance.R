# Panel B: is the transplanted-stressed proteome closer to control than the
# stressed one? Compare centroid distances to Ctl in full protein space. The
# permutation shuffles the PHE / Phe+Mito labels among the 12 stressed samples;
# under the null (transplant does not move cells toward Ctl) the ratio sits at 1.

build_proteome_distance <- function(mat, groups, n_perm = 2000L, seed = 42L) {
  centroid <- function(cols) rowMeans(mat[, cols, drop = FALSE])
  ctl <- centroid(groups == "Ctl")
  dist_to_ctl <- function(v) sqrt(sum((v - ctl)^2))

  d_phe <- dist_to_ctl(centroid(groups == "PHE"))
  d_rescue <- dist_to_ctl(centroid(groups == "PHE_Mito"))
  obs <- d_rescue / d_phe

  stressed <- which(groups %in% c("PHE", "PHE_Mito"))
  n_rescue <- sum(groups == "PHE_Mito")
  set.seed(seed)
  null <- vapply(seq_len(n_perm), function(i) {
    rescue <- sample(stressed, n_rescue)
    dist_to_ctl(centroid(rescue)) / dist_to_ctl(centroid(setdiff(stressed, rescue)))
  }, numeric(1))
  p_val <- (1 + sum(null <= obs)) / (n_perm + 1)

  p <- ggplot2::ggplot(data.frame(ratio = null), ggplot2::aes(ratio)) +
    ggplot2::geom_histogram(bins = 40, fill = "grey78", color = NA) +
    ggplot2::geom_vline(xintercept = obs, color = "#D6604D", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    ggplot2::labs(
      title = "Proteome distance to control",
      subtitle = sprintf(
        "d(Phe+Mito,Ctl) / d(Phe,Ctl) = %.2f, %s (label permutation)",
        obs, fmt_p(p_val)
      ),
      x = "distance ratio under permuted stress labels", y = "permutations"
    ) +
    FIG_THEME

  tab <- tibble::tibble(
    statistic = c("d(Phe,Ctl)", "d(Phe+Mito,Ctl)", "ratio", "perm_p", "n_perm"),
    value = c(d_phe, d_rescue, obs, p_val, n_perm)
  )
  list(plot = p, table = tab)
}
