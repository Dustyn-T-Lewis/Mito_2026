#!/usr/bin/env Rscript
# Stage 03 reports: per-contrast volcano + top-25 table PDFs, overview bar
# chart, and outlier-sensitivity refit (only when Stage 01 removed outliers).

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(gridExtra)
library(limma)

source(here::here("03_DEP", "a_script", "_setup.R"))

SUM <- file.path(RPT, "03_contrast_summaries")
dir.create(SUM, recursive = TRUE, showWarnings = FALSE)

theme_dep <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom")
pal_dir <- H9C2_PAL_DIR

# Load

dal <- readRDS(file.path(DAT, "01_limma_DAList.rds"))
contrast_names <- names(dal$results)

results_list <- lapply(contrast_names, \(cn) {
  as.data.frame(read_excel(XLSX, sheet = cn))
})
names(results_list) <- contrast_names

da_summary <- as.data.frame(read_excel(XLSX, sheet = "DA_summary"))

# Per-contrast volcano + top-25 table

for (cname in contrast_names) {
  res <- results_list[[cname]] |>
    mutate(
      nlog10_pval = -log10(pmax(P.Value, 1e-300)),
      nlog10_adj  = -log10(pmax(adj.P.Val, 1e-300)),
      nlog10_pi   = -log10(pmax(pi_score, 1e-300)),
      dir_pi = factor(case_when(
        sig_pi ==  1 ~ "Up", sig_pi == -1 ~ "Down", TRUE ~ "NS"),
        levels = c("Up", "Down", "NS")))

  make_vol <- function(df, ycol, ylab, top_n, thresh) {
    top <- slice_min(df, .data[[ycol]] * -1, n = top_n, with_ties = FALSE)
    ggplot(df, aes(logFC, .data[[ycol]], color = dir_pi)) +
      geom_point(alpha = 0.35, size = 1) +
      geom_hline(yintercept = thresh, linetype = "dashed", color = "grey40",
                 linewidth = 0.4) +
      geom_text_repel(data = top, aes(label = gene), size = 2.5,
                      max.overlaps = 15, show.legend = FALSE, seed = 42) +
      scale_color_manual(values = pal_dir, drop = FALSE) +
      labs(x = expression(log[2]~FC), y = ylab) +
      theme_dep + theme(legend.position = "none")
  }

  p1 <- make_vol(res, "nlog10_pval", expression(-log[10](P)), 10,
                 -log10(0.01)) + ggtitle("Nominal P < 0.01")
  p2 <- make_vol(res, "nlog10_adj", expression(-log[10](FDR)), 10,
                 -log10(H9C2_FDR_EXPLOR)) +
    ggtitle(sprintf("FDR < %.2f", H9C2_FDR_EXPLOR))
  p3 <- make_vol(res, "nlog10_pi", expression(-log[10](Pi)), 10,
                 -log10(H9C2_PI_THRESH)) +
    ggtitle(sprintf("Pi < %.2f", H9C2_PI_THRESH)) +
    theme(legend.position = "right") + labs(color = NULL)

  tbl_data <- res |>
    slice_min(pi_score, n = 25, with_ties = FALSE) |>
    transmute(Gene = gene, logFC = sprintf("%.3f", logFC),
              P = formatC(P.Value, format = "e", digits = 2),
              FDR = formatC(adj.P.Val, format = "e", digits = 2),
              Pi = formatC(pi_score, format = "e", digits = 2))

  p_tbl <- tableGrob(tbl_data, rows = NULL,
    theme = ttheme_minimal(base_size = 8,
      core    = list(fg_params = list(hjust = 0, x = 0.02)),
      colhead = list(fg_params = list(hjust = 0, x = 0.02, fontface = "bold"))))

  n_fdr <- sum(res$adj.P.Val < H9C2_FDR_EXPLOR, na.rm = TRUE)
  n_pi  <- sum(res$sig_pi != 0)
  cdir  <- file.path(SUM, cname)
  dir.create(cdir, showWarnings = FALSE)

  # Clean biological label for this contrast (no aging/training wording).
  clabel <- H9C2_CONTRAST_LABELS[[cname]] %||% cname

  pdf(file.path(cdir, "summary.pdf"), width = 16, height = 14)
  print((p1 | p2 | p3) + plot_annotation(
    title = sprintf("%s  (%s)", cname, clabel),
    subtitle = sprintf("FDR<%.2f: %d | Pi<%.2f: %d | %d proteins",
                       H9C2_FDR_EXPLOR, n_fdr, H9C2_PI_THRESH, n_pi, nrow(res)),
    theme = theme(plot.title = element_text(face = "bold", size = 14))))
  grid::grid.newpage()
  grid::grid.draw(arrangeGrob(p_tbl,
    top = grid::textGrob(sprintf("%s — Top 25 by Pi", cname),
                         gp = grid::gpar(fontface = "bold", fontsize = 14))))
  dev.off()
  message(sprintf("  %s: FDR<%.2f=%d, Pi<%.2f=%d",
                  cname, H9C2_FDR_EXPLOR, n_fdr, H9C2_PI_THRESH, n_pi))
}

# Overview bar chart

sc <- list_rbind(lapply(contrast_names, \(cname) {
  res <- results_list[[cname]]
  bind_rows(
    tibble(contrast = cname, criterion = sprintf("FDR < %.2f", H9C2_FDR_EXPLOR),
           up = sum(res$adj.P.Val < H9C2_FDR_EXPLOR & res$logFC > 0, na.rm = TRUE),
           down = sum(res$adj.P.Val < H9C2_FDR_EXPLOR & res$logFC < 0, na.rm = TRUE)),
    tibble(contrast = cname, criterion = sprintf("Pi < %.2f", H9C2_PI_THRESH),
           up = sum(res$sig_pi == 1), down = sum(res$sig_pi == -1)))
})) |>
  pivot_longer(c(up, down), names_to = "direction", values_to = "count") |>
  mutate(signed = if_else(direction == "down", -count, count),
         direction = factor(stringr::str_to_title(direction),
                            levels = c("Up", "Down")),
         contrast = factor(contrast, levels = contrast_names),
         criterion = factor(criterion))

p_bar <- ggplot(sc, aes(contrast, signed, fill = direction)) +
  geom_col(position = "identity", width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(data = filter(sc, signed > 0),
            aes(label = count), vjust = -0.3, size = 3) +
  geom_text(data = filter(sc, signed < 0),
            aes(label = count), vjust = 1.3, size = 3) +
  facet_wrap(~criterion) +
  scale_fill_manual(values = c(Up = "#D6604D", Down = "#4393C3")) +
  labs(title = "H9c2 Differential Abundance Overview",
       x = NULL, y = "DA proteins (Up / Down)", fill = NULL) +
  theme_dep + theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Outlier sensitivity (only if Stage 01 removed outliers)

int <- readRDS(file.path(NORM_DAT, "00_report_intermediates.rds"))
outlier_ids <- int$outlier_ids

run_limma_sens <- function(mat, meta) {
  meta$Group <- factor(meta$Group, levels = H9C2_GROUP_LEVELS)
  design <- model.matrix(~ 0 + Group, data = meta)
  colnames(design) <- gsub("^Group", "", colnames(design))
  fit <- lmFit(mat, design)
  cv  <- h9c2_parse_contrasts()
  cm  <- makeContrasts(contrasts = cv, levels = design)
  colnames(cm) <- names(cv)
  eBayes(contrasts.fit(fit, cm))
}

if (length(outlier_ids) > 0) {
  # Pre-outlier matrix is already cycloess-normalized in Stage 01 (m12).
  full_fit <- run_limma_sens(int$data_pre_outlier_norm, int$meta_pre_outlier)

  red_mat <- as.matrix(readRDS(NORM_DALIST)$data)
  red_fit <- run_limma_sens(red_mat, readRDS(NORM_DALIST)$metadata)

  sens_compare <- list_rbind(lapply(contrast_names, \(cname) {
    full_res <- topTable(full_fit, coef = cname, number = Inf, sort.by = "none")
    red_res  <- topTable(red_fit,  coef = cname, number = Inf, sort.by = "none")
    shared   <- intersect(rownames(full_res), rownames(red_res))
    tibble(Contrast = cname,
           FDR_full    = sum(full_res$adj.P.Val < H9C2_FDR_EXPLOR, na.rm = TRUE),
           FDR_reduced = sum(red_res$adj.P.Val < H9C2_FDR_EXPLOR, na.rm = TRUE),
           Pi_full     = sum(full_res$P.Value^abs(full_res$logFC) < H9C2_PI_THRESH, na.rm = TRUE),
           Pi_reduced  = sum(red_res$P.Value^abs(red_res$logFC) < H9C2_PI_THRESH, na.rm = TRUE),
           Pearson_r   = cor(full_res[shared, "logFC"], red_res[shared, "logFC"],
                             use = "complete.obs"),
           Spearman_rho = cor(full_res[shared, "logFC"], red_res[shared, "logFC"],
                              use = "complete.obs", method = "spearman"))
  }))
} else {
  message(sprintf("Stage 01 flagged %d consensus outliers — sensitivity refit N/A.",
                  length(outlier_ids)))
  sens_compare <- tibble(note = "no consensus outliers; full = reduced dataset")
}

# Add sensitivity sheet to xlsx
wb <- loadWorkbook(XLSX)
if ("outlier_sensitivity" %in% names(wb)) removeWorksheet(wb, "outlier_sensitivity")
addWorksheet(wb, "outlier_sensitivity")
writeData(wb, "outlier_sensitivity", sens_compare,
  headerStyle = createStyle(textDecoration = "bold", fgFill = "#DCE6F1"))
freezePane(wb, "outlier_sensitivity", firstRow = TRUE)
setColWidths(wb, "outlier_sensitivity", cols = seq_len(ncol(sens_compare)),
             widths = "auto")
saveWorkbook(wb, XLSX, overwrite = TRUE)

# Assemble overview PDF

if ("Pearson_r" %in% names(sens_compare)) {
  sens_display <- sens_compare |>
    mutate(across(c(Pearson_r, Spearman_rho), \(x) sprintf("%.3f", x)))
} else {
  sens_display <- sens_compare
}
p_sens <- tableGrob(sens_display, rows = NULL,
  theme = ttheme_minimal(base_size = 10,
    core    = list(fg_params = list(hjust = 0.5)),
    colhead = list(fg_params = list(fontface = "bold", hjust = 0.5))))

pdf(file.path(RPT, "02_dep_overview.pdf"), width = 14, height = 10)
print(p_bar + plot_annotation(
  title = "H9c2 Differential Abundance Overview",
  theme = theme(plot.title = element_text(face = "bold", size = 16))))
print(wrap_elements(p_sens) + plot_annotation(
  title = "Outlier-Removal Sensitivity",
  subtitle = if (length(outlier_ids) > 0)
    sprintf("Removed: %s", paste(outlier_ids, collapse = ", "))
  else "No outliers removed in Stage 01",
  theme = theme(plot.title = element_text(face = "bold", size = 14),
                plot.subtitle = element_text(size = 11))))
dev.off()
message("Saved: ", file.path(RPT, "02_dep_overview.pdf"))
