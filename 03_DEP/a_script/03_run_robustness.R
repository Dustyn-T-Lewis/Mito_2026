#!/usr/bin/env Rscript
# Stage 03 robustness: mito-effect magnitude (paired-by-protein |logFC| of
# CTLvMITO vs PHEvPHE_MITO), bootstrap CIs, two-sample power analysis, and
# imputation-sensitivity refit. Appends sheets to 03_DEP_results.xlsx.

library(dplyr)
library(tibble)
library(purrr)
library(readxl)
library(openxlsx)
library(boot)
library(pwr)
library(proteoDA)

set.seed(42)
source(here::here("03_DEP", "a_script", "_setup.R"))

# Load

dal <- readRDS(file.path(DAT, "01_limma_DAList.rds"))
contrast_names <- names(dal$results)
meta <- as.data.frame(dal$metadata)

results_list <- lapply(contrast_names, \(cn) {
  as.data.frame(read_excel(XLSX, sheet = cn))
})
names(results_list) <- contrast_names

# 1. Mito-effect magnitude (CTLvMITO vs PHEvPHE_MITO |logFC|)
# YvO compared training-response magnitude in young vs old muscle. The H9c2
# analogue asks whether mitochondrial transplantation moves the proteome more
# in healthy cells (CTLvMITO) or in PHE-stressed cells (PHEvPHE_MITO).

mag_df <- tibble(
  uniprot_id = results_list[["CTLvMITO"]]$uniprot_id,
  abs_lfc_mito_alone = abs(results_list[["CTLvMITO"]]$logFC),
  abs_lfc_mito_under_phe = abs(results_list[["PHEvPHE_MITO"]]$logFC)
) |> filter(!is.na(abs_lfc_mito_alone) & !is.na(abs_lfc_mito_under_phe))

ks_res <- ks.test(mag_df$abs_lfc_mito_alone, mag_df$abs_lfc_mito_under_phe)
fk_res <- fligner.test(
  abs_lfc ~ contrast,
  data = data.frame(
    abs_lfc  = c(mag_df$abs_lfc_mito_alone, mag_df$abs_lfc_mito_under_phe),
    contrast = rep(c("Mito_alone", "Mito_under_PHE"), each = nrow(mag_df))))
wx_res <- wilcox.test(mag_df$abs_lfc_mito_alone,
                      mag_df$abs_lfc_mito_under_phe, paired = TRUE)
diffs  <- mag_df$abs_lfc_mito_alone - mag_df$abs_lfc_mito_under_phe
cliff  <- (sum(diffs > 0) - sum(diffs < 0)) / length(diffs)
cliff_mag <- case_when(
  abs(cliff) < 0.147 ~ "negligible", abs(cliff) < 0.33 ~ "small",
  abs(cliff) < 0.474 ~ "medium",     TRUE ~ "large")

mito_magnitude <- tibble(
  test = c("KS", "Fligner-Killeen", "Wilcoxon signed-rank", "Cliff's delta"),
  statistic = c(ks_res$statistic, fk_res$statistic, wx_res$statistic, cliff),
  p_value = c(ks_res$p.value, fk_res$p.value, wx_res$p.value, NA),
  interpretation = c(
    ifelse(ks_res$p.value < 0.05, "Distributions differ", "No difference"),
    ifelse(fk_res$p.value < 0.05, "Variance differs", "No difference"),
    ifelse(wx_res$p.value < 0.05, "Paired shift significant", "No shift"),
    sprintf("%s (d=%.3f; >0 = larger mito effect in healthy/Ctl cells)",
            stringr::str_to_title(cliff_mag), cliff)))

message(sprintf("Mito-effect magnitude: KS p=%.2g, Cliff d=%.3f (%s)",
                ks_res$p.value, cliff, cliff_mag))

# 2. Bootstrap CI (median |logFC|, BCa, 10k reps)
# NOTE: boot.ci(type = "bca")'s default influence estimator regresses the
# bootstrap statistics on the R x n frequency matrix (here ~1e4 x ~4.8e3) — a
# pathologically slow least-squares solve. Supply jackknife influence values
# explicitly instead: O(n) median computations, BCa interval unchanged.

median_stat <- function(d, i) median(d[i])

boot_df <- list_rbind(lapply(contrast_names, \(cname) {
  vals <- abs(results_list[[cname]]$logFC)
  vals <- vals[!is.na(vals)]
  b  <- boot(vals, median_stat, R = 10000)
  L  <- empinf(data = vals, statistic = median_stat, stype = "i", type = "jack")
  ci <- tryCatch(boot.ci(b, type = "bca", L = L),
                 error = \(e) boot.ci(b, type = "perc"))
  ci_lo <- if (!is.null(ci$bca)) ci$bca[4] else ci$percent[4]
  ci_hi <- if (!is.null(ci$bca)) ci$bca[5] else ci$percent[5]
  tibble(contrast = cname, median_absLFC = median(vals),
         ci_lower = ci_lo, ci_upper = ci_hi,
         boot_se = sd(b$t), n_proteins = length(vals))
}))

# 3. Power analysis (min detectable logFC at 80% power)
# All H9c2 contrasts are two-sample comparisons of independent groups (n = 6).
# YvO's paired branch / within-subject correlation are removed. The Interaction
# is a difference-of-differences across 4 groups — its effective contrast SD is
# larger (factor 2 vs sqrt(2)); its power figure is therefore approximate.

fit       <- dal$eBayes_fit
sigma_res <- sqrt(mean(fit$sigma^2, na.rm = TRUE))
n_pergrp  <- min(table(meta$group))

power_df <- list_rbind(lapply(contrast_names, \(cname) {
  is_interaction <- cname == "Interaction"
  eff_sig <- if (is_interaction) sigma_res * 2 else sigma_res * sqrt(2)
  pw <- pwr.t.test(n = n_pergrp, d = NULL, sig.level = H9C2_PVAL_THRESH,
                   power = 0.80, type = "two.sample")
  tibble(contrast = cname, n_per_group = n_pergrp,
         design = if (is_interaction) "difference-of-differences (approx.)"
                  else "two-sample",
         effective_sigma = round(eff_sig, 4),
         min_detectable_d = round(pw$d, 4),
         min_detectable_logFC = round(pw$d * eff_sig, 4),
         power = 0.80, alpha = H9C2_PVAL_THRESH)
}))

# 4. Imputation sensitivity
# Refit the same contrasts on the missForest-imputed matrix and correlate
# t-statistics against the non-imputed fit.

sens_df <- tibble(contrast = character(), spearman_rho = numeric(),
                  p_value = numeric(), n_proteins = integer())

if (file.exists(IMP_DALIST)) {
  dal_imp_raw <- readRDS(IMP_DALIST)
  imp_mat <- as.matrix(dal_imp_raw$data)   # rownames = uniprot_id

  # Align imputed-matrix rows to the DEP annotation row order (Stage 02 sorted
  # by uniprot_id for missForest determinism).
  ann_dep <- as.data.frame(dal$annotation)
  ord <- match(ann_dep$uniprot_id, rownames(imp_mat))
  if (any(is.na(ord)))
    stop("Imputed DAList is missing proteins from the DEP annotation.")
  imp_mat <- imp_mat[ord, , drop = FALSE]
  stopifnot(identical(rownames(imp_mat), ann_dep$uniprot_id))

  shared_samps <- intersect(colnames(dal$data), colnames(imp_mat))

  dal_imp <- DAList(
    data       = imp_mat[, shared_samps],
    annotation = ann_dep,
    metadata   = meta[meta$sample_id %in% shared_samps, ],
    tags       = list(norm_method = "cycloess_imputed"))
  dal_imp <- add_design(dal_imp, H9C2_DESIGN_FORMULA)
  colnames(dal_imp$design$design_matrix) <- gsub("^group", "",
    colnames(dal_imp$design$design_matrix))
  dal_imp <- add_contrasts(dal_imp, contrasts_vector = H9C2_CONTRASTS)
  dal_imp <- fit_limma_model(dal_imp)
  dal_imp <- extract_DA_results(dal_imp, pval_thresh = H9C2_PVAL_THRESH,
                                lfc_thresh = 0, adj_method = "BH")

  comb <- readr::read_csv(file.path(DAT, "03_combined_results.csv"),
                          show_col_types = FALSE)

  sens_df <- list_rbind(lapply(contrast_names, \(cname) {
    t_col <- paste0("t_", cname)
    if (!(t_col %in% names(comb)) || !(cname %in% names(dal_imp$results)))
      return(NULL)
    imp_t <- dal_imp$results[[cname]] |>
      rownames_to_column("uniprot_id") |>
      select(uniprot_id, t_imp = t)
    merged <- inner_join(
      comb |> select(uniprot_id, t_nonimp = all_of(t_col)),
      imp_t, by = join_by(uniprot_id)) |>
      filter(!is.na(t_nonimp) & !is.na(t_imp))
    sp <- suppressWarnings(
      cor.test(merged$t_nonimp, merged$t_imp, method = "spearman"))
    tibble(contrast = cname, spearman_rho = round(sp$estimate, 4),
           p_value = sp$p.value, n_proteins = nrow(merged))
  }))
  message("Imputation sensitivity: ", paste(sprintf("%s rho=%.3f",
    sens_df$contrast, sens_df$spearman_rho), collapse = " | "))
} else {
  message("Imputed DAList not found — skipping sensitivity")
}

# 4b. Reinjection + robust-eBayes sensitivity
# Reinjection (r-suffixed re-runs) is imbalanced across groups (Ctl 1, Mito 0,
# PHE 2, PHE_Mito 4) and is NOT in the primary `~ 0 + group` model. Refit with
# reinjected as an additive covariate and correlate t-stats against the primary
# fit — high concordance means the imbalance does not drive the calls. Separately
# refit the primary design with robust eBayes: CTLvMITO and PHEvPHE_MITO call
# >20% of proteins DA, stretching the standard eBayes "most unchanged" prior.

prim_mat <- as.matrix(dal$data)
grp      <- factor(meta$group, levels = H9C2_GROUP_LEVELS)
cv       <- h9c2_parse_contrasts()

fit_contrasts <- function(design, robust = FALSE) {
  colnames(design) <- gsub("^grp", "", colnames(design))
  cm <- limma::makeContrasts(contrasts = cv, levels = design)
  colnames(cm) <- names(cv)
  limma::eBayes(limma::contrasts.fit(limma::lmFit(prim_mat, design), cm),
                robust = robust)
}

comb_t <- readr::read_csv(file.path(DAT, "03_combined_results.csv"),
                          show_col_types = FALSE)

norm_meta <- as.data.frame(readRDS(NORM_DALIST)$metadata)
reinj_vec <- as.logical(norm_meta$Reinjected[match(meta$sample_id, norm_meta$Col_ID)])

reinjection_sensitivity <- tibble(contrast = character(), spearman_rho = numeric(),
                                  p_value = numeric(), n_proteins = integer())
if (length(unique(reinj_vec)) == 2) {
  fit_reinj <- fit_contrasts(model.matrix(~ 0 + grp + reinj_vec))
  reinjection_sensitivity <- list_rbind(lapply(contrast_names, \(cn) {
    tc <- paste0("t_", cn)
    if (!(tc %in% names(comb_t))) return(NULL)
    rt <- limma::topTable(fit_reinj, coef = cn, number = Inf, sort.by = "none")
    m  <- inner_join(comb_t |> select(uniprot_id, t_prim = all_of(tc)),
                     rt |> rownames_to_column("uniprot_id") |> select(uniprot_id, t_reinj = t),
                     by = join_by(uniprot_id)) |>
      filter(!is.na(t_prim) & !is.na(t_reinj))
    sp <- suppressWarnings(cor.test(m$t_prim, m$t_reinj, method = "spearman"))
    tibble(contrast = cn, spearman_rho = round(sp$estimate, 4),
           p_value = sp$p.value, n_proteins = nrow(m))
  }))
  message("Reinjection sensitivity: ", paste(sprintf("%s rho=%.3f",
    reinjection_sensitivity$contrast, reinjection_sensitivity$spearman_rho), collapse = " | "))
} else {
  message("Reinjected has <2 levels — skipping reinjection sensitivity")
}

fit_def <- fit_contrasts(model.matrix(~ 0 + grp), robust = FALSE)
fit_rob <- fit_contrasts(model.matrix(~ 0 + grp), robust = TRUE)
ebayes_robust_sensitivity <- list_rbind(lapply(contrast_names, \(cn) {
  d <- limma::topTable(fit_def, coef = cn, number = Inf, sort.by = "none")
  r <- limma::topTable(fit_rob, coef = cn, number = Inf, sort.by = "none")
  tibble(contrast = cn,
         FDR_default = sum(d$adj.P.Val < H9C2_FDR_EXPLOR, na.rm = TRUE),
         FDR_robust  = sum(r$adj.P.Val < H9C2_FDR_EXPLOR, na.rm = TRUE),
         Pi_default  = sum(d$P.Value^abs(d$logFC) < H9C2_PI_THRESH, na.rm = TRUE),
         Pi_robust   = sum(r$P.Value^abs(r$logFC) < H9C2_PI_THRESH, na.rm = TRUE))
}))
message("Robust-eBayes sensitivity computed for ", length(contrast_names), " contrasts")

# 5. Add robustness sheets to xlsx

write_sheet <- function(wb, name, data) {
  addWorksheet(wb, name)
  writeData(wb, name, data,
    headerStyle = createStyle(textDecoration = "bold", fgFill = "#DCE6F1"))
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
}

wb <- loadWorkbook(XLSX)
robustness_sheets <- c("mito_effect_magnitude", "bootstrap_ci",
                       "power_analysis", "imputation_sensitivity",
                       "reinjection_sensitivity", "ebayes_robust_sensitivity")
for (s in intersect(robustness_sheets, names(wb))) removeWorksheet(wb, s)
write_sheet(wb, "mito_effect_magnitude", mito_magnitude)
write_sheet(wb, "bootstrap_ci",          boot_df)
write_sheet(wb, "power_analysis",        power_df)
if (nrow(sens_df) > 0) {
  write_sheet(wb, "imputation_sensitivity", sens_df)
}
if (nrow(reinjection_sensitivity) > 0) {
  write_sheet(wb, "reinjection_sensitivity", reinjection_sensitivity)
}
write_sheet(wb, "ebayes_robust_sensitivity", ebayes_robust_sensitivity)
saveWorkbook(wb, XLSX, overwrite = TRUE)

# Optional Box delivery

box <- Sys.getenv("H9C2_BOX_SUPP", "")
if (nzchar(box) && dir.exists(box)) {
  box_tbl <- file.path(box, "tables")
  dir.create(box_tbl, recursive = TRUE, showWarnings = FALSE)
  file.copy(XLSX, file.path(box_tbl, "S03_Table_DEP.xlsx"), overwrite = TRUE)
  message("Copied to Box: tables/S03_Table_DEP.xlsx")
}

message("Done: robustness analyses added to ", basename(XLSX))
