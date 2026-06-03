#!/usr/bin/env Rscript
# Stage 01 diagnostics: builds 04_diagnostics.pdf from
# 00_report_intermediates.rds.

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)

source(here::here("01_normalization", "a_script", "_setup.R"))

int <- readRDS(file.path(DAT, "00_report_intermediates.rds"))
list2env(int, envir = environment())

pal_group <- H9C2_PAL_GROUP
theme_qc  <- theme_minimal(base_size = 12)

# Page 1: filtering & missingness

p_filter <- ggplot(filter_bar_data, aes(step, n, fill = status)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c(Retained = "#2166AC", Removed = "#B2182B")) +
  labs(x = NULL, y = "Proteins", fill = NULL, title = "Protein retention") +
  theme_qc + theme(axis.text.x = element_text(angle = 25, hjust = 1))

p_miss <- ggplot(miss_bar_data,
                 aes(reorder(Col_ID, -n * (status == "Detected")),
                     n, fill = status)) +
  geom_col(aes(alpha = is_outlier), width = 0.8) +
  scale_fill_manual(values = c(Detected = "#2166AC", Missing = "#D6604D")) +
  scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.4), guide = "none") +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Proteins", fill = NULL,
       title = "Per-sample detection (outliers faded)") +
  theme_qc + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
                   strip.text = element_text(face = "bold"))

# Page 2: outlier diagnostics (4-method consensus)

p_out_miss <- ggplot(outlier_diag, aes(reorder(Sample_ID, pct_missing),
                                       pct_missing, color = Group)) +
  geom_point(size = 3) +
  geom_hline(yintercept = miss_thresh, linetype = "dashed",
             color = "red", alpha = 0.5) +
  geom_text_repel(data = \(d) filter(d, miss_flag),
                  aes(label = Sample_ID), size = 2.5, show.legend = FALSE) +
  scale_color_manual(values = pal_group) +
  labs(x = "Sample", y = "Missingness (%)",
       title = "A: Sample missingness",
       subtitle = sprintf("IQR fence: %.1f%% | %d flagged",
                          miss_thresh, sum(outlier_diag$miss_flag))) +
  theme_qc + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))

pca_outlier_df <- pca_pre$scores |>
  left_join(outlier_diag |> select(Col_ID, pca_flag), by = join_by(Col_ID))

p_out_pca <- ggplot(pca_outlier_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text_repel(data = \(d) filter(d, pca_flag),
                  aes(label = Sample_ID), size = 2.5, show.legend = FALSE) +
  scale_color_manual(values = pal_group) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_pre$var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", pca_pre$var_exp[2]),
       title = "B: PCA Mahalanobis",
       subtitle = sprintf("p < %.2f | %d flagged",
                          mahal_p, sum(pca_outlier_df$pca_flag))) +
  coord_fixed() + theme_qc

p_out_mad <- ggplot(outlier_diag,
                    aes(reorder(Sample_ID, sample_median), sample_median,
                        color = Group)) +
  geom_point(size = 2.5) +
  geom_text_repel(data = \(d) filter(d, mad_flag),
                  aes(label = Sample_ID), size = 2.5, show.legend = FALSE) +
  geom_hline(yintercept = global_med) +
  geom_hline(yintercept = global_med + c(-1, 1) * mad_k * mad_val,
             linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = pal_group) +
  labs(x = "Sample", y = "Median log2 intensity",
       title = "C: MAD median intensity",
       subtitle = sprintf("%dx MAD | %d flagged",
                          mad_k, sum(outlier_diag$mad_flag))) +
  theme_qc + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))

p_out_cor <- ggplot(outlier_diag,
                    aes(reorder(Sample_ID, median_cor), median_cor,
                        color = Group)) +
  geom_point(size = 2.5) +
  geom_text_repel(data = \(d) filter(d, cor_flag),
                  aes(label = Sample_ID), size = 2.5, show.legend = FALSE) +
  geom_hline(yintercept = median(outlier_diag$median_cor)) +
  geom_hline(yintercept = median(outlier_diag$median_cor) -
               mad_k * mad(outlier_diag$median_cor),
             linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = pal_group) +
  labs(x = "Sample", y = "Median pairwise correlation",
       title = "D: Inter-sample correlation",
       subtitle = sprintf("%dx MAD | %d flagged",
                          mad_k, sum(outlier_diag$cor_flag))) +
  theme_qc + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))

# Page 3: post-normalization PCA

p_pca_post <- ggplot(pca_post$scores, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(aes(group = Group), type = "norm", level = 0.68, linewidth = 0.7) +
  scale_color_manual(values = pal_group) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_post$var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", pca_post$var_exp[2]),
       title = "Post-normalization PCA") +
  theme_qc + theme(legend.position = "bottom")

# Page 4: variability

p_var <- ggplot(samp_var, aes(reorder(Sample_ID, iqr), iqr, color = Group)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pal_group) +
  labs(x = "Sample", y = "IQR (log2 intensity)",
       title = "Per-sample variability",
       subtitle = sprintf("%d samples", nrow(samp_var))) +
  theme_qc + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
                   legend.position = "bottom")

eta2_df <- data.frame(eta2 = eta2_vals[!is.na(eta2_vals)])
p_eta2 <- ggplot(eta2_df, aes(eta2)) +
  geom_histogram(bins = 50, fill = "#2166AC", color = "white", alpha = 0.8) +
  geom_vline(xintercept = median(eta2_df$eta2), linetype = "dashed",
             color = "red") +
  annotate("text", x = median(eta2_df$eta2) + 0.02, y = Inf,
           vjust = 2, size = 3.5, color = "red",
           label = sprintf("median = %.2f", median(eta2_df$eta2))) +
  labs(x = expression(eta^2), y = "Proteins",
       title = "Variance partition by group") +
  theme_qc

# Assemble PDF

pdf(file.path(RPT, "04_diagnostics.pdf"), width = 20, height = 10)

print(
  p_filter / p_miss + plot_layout(heights = c(1, 1.2)) +
    plot_annotation(
      title = "Protein Filtering & Detection",
      subtitle = sprintf("%d raw -> %d retained | %d samples",
                         n_raw, dal_nrow, dal_ncol),
      theme = theme(plot.title = element_text(size = 18, face = "bold"),
                    plot.subtitle = element_text(size = 14))))

shared_legend <- get_legend(
  p_out_miss + theme(legend.position = "bottom",
                     legend.justification = "center"))
print(
  ((p_out_miss + theme(legend.position = "none")) |
   (p_out_pca  + theme(legend.position = "none"))) /
  ((p_out_mad  + theme(legend.position = "none")) |
   (p_out_cor  + theme(legend.position = "none"))) /
  wrap_elements(shared_legend) +
    plot_layout(heights = c(1, 1, 0.08)) +
    plot_annotation(
      title = "Outlier Diagnostics (4-method consensus)",
      subtitle = sprintf(">=%d/4 agreement | %d removed", outlier_k, n_outliers),
      theme = theme(plot.title = element_text(size = 18, face = "bold"),
                    plot.subtitle = element_text(size = 13))))

print(p_pca_post + plot_annotation(
  title = "Post-Normalization QC",
  theme = theme(plot.title = element_text(size = 18, face = "bold"))))

print(
  (p_var | p_eta2) + plot_annotation(
    title = "Variability Summary",
    theme = theme(plot.title = element_text(size = 18, face = "bold"))))

dev.off()
message("Saved: ", file.path(RPT, "04_diagnostics.pdf"))

# Optional Box delivery (skipped if $H9C2_BOX_SUPP is unset).

box <- Sys.getenv("H9C2_BOX_SUPP", "")
if (nzchar(box) && dir.exists(box)) {
  box_tbl <- file.path(box, "tables")
  dir.create(box_tbl, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(DAT, "01_normalization.xlsx"),
            file.path(box_tbl, "S01_Table_Normalization.xlsx"), overwrite = TRUE)
  message("Copied to Box: tables/S01_Table_Normalization.xlsx")
}
