#!/usr/bin/env Rscript
# Stage 02 reports: builds 01_missingness_report.pdf + 02_imputation_report.pdf
# from 00_report_intermediates.rds. Benchmark page is optional (only if the
# benchmark ranking CSV exists).

library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)
library(readr)
library(scales)

source(here::here("02_Imputation", "a_script", "_setup.R"))

rpt <- readRDS(file.path(DAT, "00_report_intermediates.rds"))
list2env(rpt, envir = environment())

mc <- miss_class  # alias for brevity

THM <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank())

# Report 1: Missingness (1 page)

p_miss_hist <- mc |>
  filter(classification != "Complete") |>
  ggplot(aes(pct_miss, fill = classification)) +
  geom_histogram(binwidth = 5, boundary = 0, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = PAL_CLASS, name = NULL) +
  scale_x_continuous(labels = label_percent(scale = 1)) +
  labs(title = "A. Per-protein missingness", x = "% missing", y = "Count") +
  THM + theme(legend.position = "top")

class_counts <- mc |>
  count(classification) |>
  mutate(classification = factor(classification,
                                 levels = c("Complete", "MAR", "MNAR")))
p_class_bar <- ggplot(class_counts, aes(1, n, fill = classification)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = paste0(classification, "\n", n)),
            position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = PAL_CLASS) + coord_flip() +
  labs(title = "B. Classification", x = NULL, y = "Proteins") +
  THM + theme(legend.position = "none", axis.text.y = element_blank())

p_int_vs_miss <- mc |>
  filter(classification != "Complete") |>
  ggplot(aes(mean_intensity, pct_miss, color = classification)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = PAL_CLASS, name = NULL) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  labs(title = "C. Intensity vs missingness",
       x = "Mean log2 intensity", y = "% missing") +
  THM + theme(legend.position = "top")

sample_miss <- tibble(
  sample = colnames(mat), n_miss = colSums(is.na(mat))) |>
  left_join(meta, by = join_by(sample == Col_ID)) |>
  mutate(Group = factor(Group, levels = names(PAL_GROUP)))

p_samp_miss <- ggplot(sample_miss, aes(reorder(Sample_ID, n_miss), n_miss,
                                       fill = Group)) +
  geom_col() + scale_fill_manual(values = PAL_GROUP, name = NULL) + coord_flip() +
  labs(title = "D. Per-sample missingness", x = NULL, y = "Missing proteins") +
  THM + theme(legend.position = "top", axis.text.y = element_text(size = 6))

pdf(file.path(RPT, "01_missingness_report.pdf"), width = 16, height = 10)
print(
  (p_miss_hist | p_class_bar) / (p_int_vs_miss | p_samp_miss) +
    plot_annotation(
      title = "Missingness Report",
      subtitle = sprintf("%s proteins x %d samples | %.1f%% missing | %s",
                         comma(nrow(mat)), ncol(mat), pct_miss,
                         classification_method),
      theme = theme(plot.title = element_text(face = "bold", size = 14))))
dev.off()
message("Saved: ", file.path(RPT, "01_missingness_report.pdf"))

# Report 2: Imputation quality
# Page 1 (benchmark) is OPTIONAL — only if a benchmark ranking CSV is present.
# Page 2 (quality) is always produced.

pages <- list()

if (file.exists(BENCH_RANKING)) {
  bm <- read_csv(BENCH_RANKING, show_col_types = FALSE)
  bench_classes <- sort(unique(bm$class))
  PAL_BENCH <- setNames(
    c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628",
      "#F781BF", "#999999")[seq_along(bench_classes)], bench_classes)

  top10 <- bm |> slice_min(rank, n = 10)
  p_bench_rank <- top10 |>
    mutate(method = factor(method, levels = rev(method)),
           is_mf = method == "missForest") |>
    ggplot(aes(composite, method, fill = class, alpha = is_mf)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.3f", composite)), hjust = -0.1, size = 3) +
    scale_fill_manual(values = PAL_BENCH, name = "Method class") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.6), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "A. Composite ranking (top 10)", x = "Composite score", y = NULL) +
    THM

  p_bench_scatter <- bm |>
    mutate(is_mf = method == "missForest") |>
    ggplot(aes(nrmse_mcar, fc_rho, color = class)) +
    geom_point(aes(size = is_mf), alpha = 0.7) +
    geom_text_repel(aes(label = method), size = 2.8, max.overlaps = 15) +
    scale_color_manual(values = PAL_BENCH, name = "Method class") +
    scale_size_manual(values = c(`TRUE` = 3.5, `FALSE` = 2), guide = "none") +
    labs(title = "B. Reconstruction vs fold-change fidelity",
         x = "NRMSE (MCAR)", y = expression("FC " * rho)) +
    THM

  pages$benchmark <- (p_bench_rank | p_bench_scatter) +
    plot_annotation(
      title = "Imputation Benchmark",
      subtitle = sprintf("#1 %s %.3f | #2 %s %.3f | OOB = %.3f",
                         bm$method[1], bm$composite[1],
                         bm$method[2], bm$composite[2], oob_error),
      theme = theme(plot.title = element_text(face = "bold", size = 14)))
} else {
  message("No benchmark ranking found — benchmark page skipped (optional).")
}

# Page 2: quality (always)
obs_vals <- as.numeric(mat[!was_na])
imp_vals <- as.numeric(mat_imp[was_na])
dens_df <- bind_rows(
  tibble(value = obs_vals, type = "Observed"),
  tibble(value = imp_vals, type = "Imputed"))

p_dens <- ggplot(dens_df, aes(value, fill = type, color = type)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c(Observed = "#377EB8", Imputed = "#E41A1C"), name = NULL) +
  scale_color_manual(values = c(Observed = "#377EB8", Imputed = "#E41A1C"), name = NULL) +
  labs(title = "A. Observed vs imputed distributions",
       x = "log2 intensity", y = "Density") +
  THM + theme(legend.position = "top")

p_mnar_audit <- ggplot(mnar_audit, aes(pre_mean, post_mean,
                                       color = imputation_reliable)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c(`TRUE` = "#4DAF4A", `FALSE` = "#E41A1C"),
                     labels = c(`TRUE` = "Reliable", `FALSE` = "Unreliable"),
                     name = NULL) +
  labs(title = "B. MNAR audit", x = "Pre-imputation mean",
       y = "Post-imputation mean") +
  THM + theme(legend.position = "top")

d_med <- median(mnar_audit$effect_d, na.rm = TRUE)
p_effect_d <- ggplot(mnar_audit, aes(effect_d)) +
  geom_histogram(binwidth = 0.1, fill = "#984EA3", color = "white", alpha = 0.7) +
  geom_vline(xintercept = d_med, linetype = "dashed") +
  annotate("text", x = d_med, y = Inf, vjust = 2, hjust = -0.1, size = 3.2,
           label = sprintf("median = %.2f", d_med)) +
  labs(title = "C. Imputation shift (MNAR)", x = "Cohen's d", y = "Count") +
  THM

pages$quality <- p_dens / (p_mnar_audit | p_effect_d) +
  plot_annotation(
    title = "Imputation Quality",
    subtitle = sprintf("missForest | %s values imputed | %d unreliable (>50%% missing)",
                       comma(total_miss_vals),
                       sum(!mnar_audit$imputation_reliable)),
    theme = theme(plot.title = element_text(face = "bold", size = 14)))

pdf(file.path(RPT, "02_imputation_report.pdf"), width = 14, height = 10)
if (!is.null(pages$benchmark)) print(pages$benchmark)
print(pages$quality)
dev.off()
message("Saved: ", file.path(RPT, "02_imputation_report.pdf"))

# Optional Box delivery

box <- Sys.getenv("H9C2_BOX_SUPP", "")
if (nzchar(box) && dir.exists(box)) {
  box_tbl <- file.path(box, "tables")
  dir.create(box_tbl, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(DAT, "02_imputation.xlsx"),
            file.path(box_tbl, "S02_Table_Imputation.xlsx"), overwrite = TRUE)
  message("Copied to Box: tables/S02_Table_Imputation.xlsx")
}
