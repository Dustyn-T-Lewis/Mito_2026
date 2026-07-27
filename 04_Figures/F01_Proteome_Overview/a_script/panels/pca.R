# PCA scores + PERMANOVA stats for F01 panel A.
pacman::p_load(dplyr, tidyr, tibble, readr, ggplot2, vegan, patchwork)
source(here::here("04_Figures", "functions", "shared_data_loaders.R"))
source(here::here("04_Figures", "functions", "f01_pca_stats.R"))

build_pca_panel <- function() {
  dal_imp <- readRDS(P05$imp_rds)
  imp_mat <- as.matrix(dal_imp$data)
  imp_meta <- as_tibble(dal_imp$metadata)
  imp_meta$Group <- factor(imp_meta$Group, levels = H9C2_GROUP_LEVELS)
  imp_mat <- imp_mat[, imp_meta$Col_ID]

  pca <- prcomp(t(imp_mat), center = TRUE, scale. = TRUE)
  var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)
  pca_df <- as.data.frame(pca$x[, 1:2]) |>
    mutate(Col_ID = rownames(pca$x)) |>
    left_join(dplyr::select(imp_meta, Col_ID, Group), by = "Col_ID") |>
    mutate(Group = factor(Group, levels = H9C2_GROUP_LEVELS))

  pairs_to_test <- list(
    Disease = c("Ctl", "PHE"),
    Transplant = c("Ctl", "Mito"),
    Rescue = c("PHE", "PHE_Mito")
  )
  stats <- pca_group_stats(imp_mat, imp_meta$Group, pairs_to_test)
  perm_R2 <- stats$overall_R2
  perm_p <- stats$overall_p
  pair_res <- stats$pairwise
  bd_p <- stats$betadisper_p
  permanova_out <- stats$table
  bd_txt <- if (is.na(bd_p)) "betadisper p = NA" else paste0("betadisper ", fmt_p(bd_p))

  GRP_LAB <- GROUP_LABELS
  GRP_SHP <- c(Ctl = 16, Mito = 17, PHE = 15, PHE_Mito = 18)
  fmt_perm <- function(role, r2, q) {
    sprintf("%s R²=%.2f, %s", role, r2, sub("^p", "q", fmt_p(q)))
  }
  stat_lines <- paste(c(
    "pairwise PERMANOVA (BH q):",
    fmt_perm("Disease", pair_res$R2[1], pair_res$q[1]),
    fmt_perm("Transplant", pair_res$R2[2], pair_res$q[2]),
    fmt_perm("Rescue", pair_res$R2[3], pair_res$q[3])
  ), collapse = "\n")

  p <- ggplot(pca_df, aes(PC1, PC2, color = Group, shape = Group)) +
    stat_ellipse(aes(fill = Group),
      geom = "polygon", alpha = 0.08,
      level = 0.80, show.legend = FALSE
    ) +
    stat_ellipse(
      level = 0.80, linewidth = 0.3, linetype = "dashed",
      show.legend = FALSE
    ) +
    geom_point(size = 1.8, alpha = 0.9) +
    annotate("label",
      x = -Inf, y = Inf, label = stat_lines, hjust = 0, vjust = 1,
      size = 1.6, color = "grey15", lineheight = 0.95,
      fill = alpha("white", 0.7),
      label.padding = unit(0.6, "mm")
    ) +
    scale_color_manual(
      values = GROUP_COLORS, labels = GRP_LAB, name = NULL,
      guide = guide_legend(nrow = 1, override.aes = list(size = 2))
    ) +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    scale_shape_manual(
      values = GRP_SHP, labels = GRP_LAB, name = NULL,
      guide = guide_legend(nrow = 1)
    ) +
    scale_x_continuous(expand = expansion(mult = 0.04)) +
    scale_y_continuous(expand = expansion(mult = 0.04)) +
    labs(
      title = "Sample PCA",
      subtitle = sprintf(
        "PERMANOVA Group R²=%.2f, %s  |  %s",
        perm_R2, fmt_p(perm_p), bd_txt
      ),
      x = sprintf("PC1 (%.1f%%)", var_pct[1]),
      y = sprintf("PC2 (%.1f%%)", var_pct[2])
    ) +
    FIG_THEME +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.5, 0.01),
      legend.justification = c(0.5, 0),
      legend.direction = "horizontal",
      legend.background = element_rect(
        fill = alpha("white", 0.75), color = "grey60", linewidth = 0.3
      ),
      legend.key = element_blank(), legend.key.size = unit(2.6, "mm"),
      legend.margin = margin(1, 1, 1, 1),
      plot.margin = margin(4, 3, 2, 2)
    )

  list(plot = p, scores = pca_df, permanova = permanova_out)
}
