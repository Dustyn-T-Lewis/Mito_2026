#!/usr/bin/env Rscript
# F01 PCA — standalone sample PCA + PERMANOVA on the imputed matrix.
# Reads existing pipeline outputs only; never re-runs 01-03.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(readr); library(ggplot2); library(vegan)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))

BASE    <- here::here("04_Figures_v2", "01_PCA")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
DAT     <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

dal_imp  <- readRDS(P05$imp_rds)
imp_mat  <- as.matrix(dal_imp$data)
imp_meta <- as_tibble(dal_imp$metadata)
imp_meta$Group <- factor(imp_meta$Group, levels = H9C2_GROUP_LEVELS)
imp_mat  <- imp_mat[, imp_meta$Col_ID]

pca <- prcomp(t(imp_mat), center = TRUE, scale. = TRUE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)
pca_df <- as.data.frame(pca$x[, 1:2]) |>
  mutate(Col_ID = rownames(pca$x)) |>
  left_join(dplyr::select(imp_meta, Col_ID, Group), by = "Col_ID") |>
  mutate(Group = factor(Group, levels = H9C2_GROUP_LEVELS))

set.seed(42); dist_mat <- dist(scale(t(imp_mat)))
set.seed(42)
perm <- adonis2(dist_mat ~ Group, data = imp_meta, permutations = 9999, by = "terms")
perm_R2 <- perm["Group", "R2"]; perm_p <- perm["Group", "Pr(>F)"]

pairs_to_test <- list(`Disease (Ctl|PHE)` = c("Ctl", "PHE"),
                      `Transplant (Ctl|Mito)` = c("Ctl", "Mito"),
                      `Rescue (PHE|PHE+Mito)` = c("PHE", "PHE_Mito"))
pair_res <- bind_rows(lapply(names(pairs_to_test), function(nm) {
  pr <- pairs_to_test[[nm]]; keep <- imp_meta$Group %in% pr
  sm <- imp_mat[, keep]; smeta <- imp_meta[keep, ]; smeta$Group <- droplevels(smeta$Group)
  set.seed(42)
  r <- adonis2(dist(scale(t(sm))) ~ Group, data = smeta, permutations = 9999, by = "terms")
  tibble(role = nm, R2 = r$R2[1], p = r$`Pr(>F)`[1])
}))

# Dispersion homogeneity check — PERMANOVA assumes equal multivariate spread.
set.seed(42)
bd_p <- permutest(betadisper(dist_mat, imp_meta$Group), permutations = 999)$tab$`Pr(>F)`[1]
if (!is.na(bd_p) && bd_p < 0.05)
  warning("Heterogeneous group dispersions (betadisper p < 0.05) — interpret PERMANOVA cautiously")

permanova_out <- bind_rows(
  tibble(role = "Group (overall)", R2 = perm_R2, p = perm_p),
  pair_res,
  tibble(role = "dispersion (betadisper)", R2 = NA_real_, p = bd_p))

GRP_SHP <- c(Ctl = 16, Mito = 17, PHE = 15, PHE_Mito = 18)
GRP_LAB <- c(Ctl = "Ctl", Mito = "Mito", PHE = "PHE", PHE_Mito = "PHE+Mito")
fmt_perm <- function(role, r2, p) sprintf("%s: R²=%.3f, %s", role, r2, fmt_p(p))
perm_label <- paste(c("PERMANOVA",
  fmt_perm("Group", perm_R2, perm_p),
  fmt_perm("Disease (Ctl|PHE)", pair_res$R2[1], pair_res$p[1]),
  fmt_perm("Transplant (Ctl|Mito)", pair_res$R2[2], pair_res$p[2]),
  fmt_perm("Rescue (PHE|PHE+Mito)", pair_res$R2[3], pair_res$p[3]),
  sprintf("dispersion (betadisper): %s", fmt_p(bd_p))), collapse = "\n")

xr <- range(pca_df$PC1); yr <- range(pca_df$PC2); ytop <- yr[2] + 0.42 * diff(yr)
pA <- ggplot(pca_df, aes(PC1, PC2, color = Group, shape = Group)) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.08, level = 0.80, show.legend = FALSE) +
  stat_ellipse(level = 0.80, linewidth = 0.3, linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 1.8, alpha = 0.9) +
  annotate("label", x = xr[1] - 0.02 * diff(xr), y = ytop + 0.03 * diff(yr), label = perm_label,
           hjust = 0, vjust = 1, size = 1.7, color = "grey20", fontface = "bold", lineheight = 0.95,
           fill = alpha("white", 0.85), label.size = 0, label.padding = unit(0.15, "lines")) +
  scale_color_manual(values = H9C2_PAL_GROUP, labels = GRP_LAB, name = NULL,
                     guide = guide_legend(ncol = 1, override.aes = list(size = 2))) +
  scale_fill_manual(values = H9C2_PAL_GROUP, guide = "none") +
  scale_shape_manual(values = GRP_SHP, labels = GRP_LAB, name = NULL, guide = guide_legend(ncol = 1)) +
  coord_cartesian(ylim = c(yr[1] - 0.02 * diff(yr), ytop + 0.03 * diff(yr))) +
  labs(title = "Sample PCA",
       subtitle = sprintf("n = %d, %s proteins (imputed)", nrow(imp_meta), format(nrow(imp_mat), big.mark = ",")),
       x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC2 (%.1f%%)", var_pct[2])) +
  FIG_THEME +
  theme(legend.position = c(0.015, 0.58), legend.justification = c(0, 0.5),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.key = element_blank(), legend.key.size = unit(2.8, "mm"),
        legend.text = element_text(size = FIG_LEGEND_TEXT + 0.5), legend.spacing.y = unit(0.3, "mm"),
        plot.margin = margin(5, 2, 1, 1))

ggsave(file.path(RPT_PDF, "MAIN_F01_pca.pdf"), pA, width = 120, height = 110, units = "mm", device = pdf_dev)
ggsave(file.path(RPT_PNG, "MAIN_F01_pca.png"), pA, width = 120, height = 110, units = "mm", dpi = 300)

build_workbook(
  file.path(DAT, "F01_supplementary.xlsx"),
  figure_title = "F01 — Sample PCA + PERMANOVA on the imputed protein matrix",
  sheet_specs = list(
    list(name = "pca_scores", df = pca_df,
         role     = "Panel coordinates — the PCA scatter points",
         contents = "PC1/PC2 scores per sample (Col_ID) with Group; % variance shown in the axis titles"),
    list(name = "permanova",  df = permanova_out,
         role     = "Stats annotation block on the figure",
         contents = "adonis2 R² and p for overall Group + Disease/Transplant/Rescue pairwise, plus betadisper dispersion p")))

message(sprintf("F01 PCA | PC1=%.1f%% PC2=%.1f%% | PERMANOVA Group R²=%.3f %s | betadisper %s",
                var_pct[1], var_pct[2], perm_R2, fmt_p(perm_p), fmt_p(bd_p)))
