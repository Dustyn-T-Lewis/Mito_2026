#!/usr/bin/env Rscript
# Stage 03: limma differential expression on cycloess-normalized non-imputed
# matrix (limma handles NAs per-protein). Independent 4-group design
# (~ 0 + group), 5 contrasts, Pi-score.

library(dplyr)
library(tibble)
library(purrr)
library(proteoDA)
library(openxlsx)

set.seed(42)
source(here::here("03_DEP", "a_script", "_setup.R"))

# 1. Load from Stage 01
# Numeric matrix from CSV for cross-machine float reproducibility.

df  <- readr::read_csv(NORM_CSV, show_col_types = FALSE)
mat <- as.matrix(df[, -(1:4)])
rownames(mat) <- df$uniprot_id    # stable unique key (cf. Stage 02)

dal_norm <- readRDS(NORM_DALIST)
ann  <- as.data.frame(dal_norm$annotation)

# 4-group design blocked by Replicate (paired passage/plate/day). No age,
# timepoint, or subject columns (YvO had all three).
meta <- tibble(
  sample_id = dal_norm$metadata$Col_ID,
  group     = factor(dal_norm$metadata$Group, levels = H9C2_GROUP_LEVELS),
  Replicate = factor(dal_norm$metadata$Replicate))

stopifnot(setequal(colnames(mat), meta$sample_id))
mat <- mat[, meta$sample_id]
message(sprintf("Loaded: %d proteins x %d samples | %.1f%% missing",
                nrow(mat), ncol(mat), 100 * mean(is.na(mat))))

# 2. Build DAList

meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta$sample_id

dal <- DAList(data = mat, annotation = ann, metadata = meta_df,
              tags = list(norm_method = "cycloess"))

# 3. Design + contrasts + fit
# Replicate is a real block (paired passage/plate/day). proteoDA picks up
# the `(1|Replicate)` term and routes through duplicateCorrelation
# (Smyth, Michaud & Scott 2005).

dal <- add_design(dal, H9C2_DESIGN_FORMULA)
colnames(dal$design$design_matrix) <- gsub("^group", "",
                                           colnames(dal$design$design_matrix))
stopifnot("design columns must be the 4 H9c2 groups" =
            setequal(colnames(dal$design$design_matrix), H9C2_GROUP_LEVELS))

dal <- add_contrasts(dal, contrasts_vector = H9C2_CONTRASTS)

dal <- fit_limma_model(dal)

# NOTE: proteoDA/limma will warn that >20% of proteins are differentially
# abundant for CTLvMITO and PHEvPHE_MITO. Empirical Bayes assumes most proteins
# are unchanged; that assumption is stretched here. This is treated as real
# biology — mitochondrial transplantation adds organelles, so a genome-wide
# shift is expected — not as a bug. The imputation-sensitivity check in
# 03_run_robustness.R (Spearman rho ~0.99) confirms the fit is stable. If a
# more conservative variance model is wanted, re-fit with robust/trended eBayes.
dal <- extract_DA_results(dal, pval_thresh = H9C2_PVAL_THRESH,
                          lfc_thresh = 0, adj_method = "BH")

saveRDS(dal, file.path(DAT, "01_limma_DAList.rds"))

# 4. proteoDA reports

tryCatch(
  write_limma_plots(dal, grouping_column = "group", output_dir = PDA,
                    table_columns = c("uniprot_id", "gene", "protein"),
                    title_column = "gene", overwrite = TRUE),
  error = \(e) message("write_limma_plots: ", conditionMessage(e)))

# 5. Extract results + Pi-score

contrast_names <- names(dal$results)
ann_df <- as.data.frame(dal$annotation)

results_list <- lapply(contrast_names, \(cname) {
  dal$results[[cname]] |>
    rownames_to_column("uniprot_id") |>
    left_join(ann_df |> select(uniprot_id, gene, protein, description),
              by = join_by(uniprot_id)) |>
    mutate(pi_score = P.Value ^ abs(logFC),  # Pi-score: Xiao 2014 PMID 22321699
           sig_pi = case_when(
             pi_score < H9C2_PI_THRESH & logFC > 0 ~  1L,
             pi_score < H9C2_PI_THRESH & logFC < 0 ~ -1L,
             TRUE ~ 0L),
           contrast = cname) |>
    select(-any_of(c("sig.PVal", "sig.FDR")))
})
names(results_list) <- contrast_names

# 6. Combined results (wide format)

data_df <- as.data.frame(dal$data)
base_df <- bind_cols(
  ann_df |> select(any_of(c("uniprot_id", "protein", "gene", "description"))),
  data_df)

for (cname in contrast_names) {
  stat_cols <- results_list[[cname]] |>
    select(uniprot_id, logFC, CI.L, CI.R, average_intensity,
           t, B, P.Value, adj.P.Val, pi_score, sig_pi)
  names(stat_cols)[-1] <- paste0(names(stat_cols)[-1], "_", cname)
  base_df <- left_join(base_df, stat_cols, by = join_by(uniprot_id))
}

readr::write_csv(base_df, file.path(DAT, "03_combined_results.csv"))

# 7. DA summary

da_summary <- list_rbind(lapply(contrast_names, \(cname) {
  res <- results_list[[cname]]
  bind_rows(
    tibble(contrast = cname, type = "up",
           sig.PVal = sum(res$P.Value < H9C2_PVAL_THRESH & res$logFC > 0, na.rm = TRUE),
           sig.FDR  = sum(res$adj.P.Val < H9C2_FDR_EXPLOR & res$logFC > 0, na.rm = TRUE),
           sig.Pi   = sum(res$sig_pi == 1, na.rm = TRUE),
           sig.FDR.05 = sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm = TRUE)),
    tibble(contrast = cname, type = "down",
           sig.PVal = sum(res$P.Value < H9C2_PVAL_THRESH & res$logFC < 0, na.rm = TRUE),
           sig.FDR  = sum(res$adj.P.Val < H9C2_FDR_EXPLOR & res$logFC < 0, na.rm = TRUE),
           sig.Pi   = sum(res$sig_pi == -1, na.rm = TRUE),
           sig.FDR.05 = sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm = TRUE)),
    tibble(contrast = cname, type = "nonsig",
           sig.PVal = sum(res$P.Value >= H9C2_PVAL_THRESH, na.rm = TRUE),
           sig.FDR  = sum(res$adj.P.Val >= H9C2_FDR_EXPLOR, na.rm = TRUE),
           sig.Pi   = sum(res$sig_pi == 0, na.rm = TRUE),
           sig.FDR.05 = sum(!(res$adj.P.Val < 0.05), na.rm = TRUE)))
}))

# 8. Build xlsx

wb <- createWorkbook()
write_h9c2_sheet(wb, "combined_results", base_df)
write_h9c2_sheet(wb, "DA_summary", da_summary)
for (cname in contrast_names) {
  res <- results_list[[cname]] |> arrange(pi_score)
  write_h9c2_sheet(wb, cname, res)
}
saveWorkbook(wb, XLSX, overwrite = TRUE)

message(sprintf("Done: %d contrasts | %d proteins -> %s/",
                length(contrast_names), nrow(dal$data), DAT))
