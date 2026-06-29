# PCA scores + PERMANOVA stats for F01 panel A.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(vegan)
  library(patchwork)
})
source(here::here("04_Figures", "functions", "02_data_paths_and_loaders.R"))

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

  set.seed(42)
  dist_mat <- dist(scale(t(imp_mat)))
  set.seed(42)
  perm <- adonis2(dist_mat ~ Group, data = imp_meta, permutations = 9999, by = "terms")
  perm_R2 <- perm["Group", "R2"]
  perm_p <- perm["Group", "Pr(>F)"]

  pairs_to_test <- list(
    `Disease` = c("Ctl", "PHE"),
    `Transplant` = c("Ctl", "Mito"),
    `Rescue` = c("PHE", "PHE_Mito")
  )
  pair_res <- bind_rows(lapply(names(pairs_to_test), function(nm) {
    pr <- pairs_to_test[[nm]]
    keep <- imp_meta$Group %in% pr
    sm <- imp_mat[, keep]
    smeta <- imp_meta[keep, ]
    smeta$Group <- droplevels(smeta$Group)
    set.seed(42)
    r <- adonis2(dist(scale(t(sm))) ~ Group,
      data = smeta,
      permutations = 9999, by = "terms"
    )
    tibble(role = nm, R2 = r$R2[1], p = r$`Pr(>F)`[1])
  }))

  set.seed(42)
  bd_p <- permutest(betadisper(dist_mat, imp_meta$Group),
    permutations = 999
  )$tab$`Pr(>F)`[1]
  if (!is.na(bd_p) && bd_p < 0.05) {
    warning("Heterogeneous group dispersions (betadisper p < 0.05)")
  }

  permanova_out <- bind_rows(
    tibble(role = "Group (overall)", R2 = perm_R2, p = perm_p),
    pair_res,
    tibble(role = "dispersion (betadisper)", R2 = NA_real_, p = bd_p)
  )

  GRP_LAB <- GROUP_LABELS
  GRP_SHP <- c(Ctl = 16, Mito = 17, PHE = 15, PHE_Mito = 18)
  fmt_perm <- function(role, r2, p) sprintf("%s R²=%.2f, %s", role, r2, fmt_p(p))
  stat_lines <- paste(c(
    fmt_perm("Disease", pair_res$R2[1], pair_res$p[1]),
    fmt_perm("Transplant", pair_res$R2[2], pair_res$p[2]),
    fmt_perm("Rescue", pair_res$R2[3], pair_res$p[3])
  ), collapse = "\n")

  n <- nrow(imp_meta)
  np <- format(nrow(imp_mat), big.mark = ",")
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
      fill = alpha("white", 0.7), label.size = 0.3,
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
        "PERMANOVA Group R²=%.2f, %s  |  n=%d, %s proteins (imputed)",
        perm_R2, fmt_p(perm_p), n, np
      ),
      x = sprintf("PC1 (%.1f%%)", var_pct[1]),
      y = sprintf("PC2 (%.1f%%)", var_pct[2])
    ) +
    FIG_THEME +
    theme(
      legend.position = c(0.5, 0.01),
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
