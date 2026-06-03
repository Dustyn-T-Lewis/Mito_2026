#!/usr/bin/env Rscript
# Stage 01: DIA-NN load, dedup, missingness filter, 4-method outlier
# consensus, cycloess normalization.

library(proteoDA)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)

set.seed(42)
source(here::here("01_normalization", "a_script", "_setup.R"))

# Helpers

run_pca <- function(mat, metadata, log_transform = TRUE) {
  for (j in seq_len(ncol(mat)))
    mat[is.na(mat[, j]), j] <- median(mat[, j], na.rm = TRUE)
  if (log_transform) mat <- log2(mat)
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  ve  <- round(summary(pca)$importance[2, 1:3] * 100, 1)
  list(pca = pca,
       scores = as.data.frame(pca$x[, 1:3]) |>
         mutate(Col_ID = rownames(pca$x)) |>
         left_join(metadata, by = join_by(Col_ID)),
       var_exp = ve)
}

write_sheet <- function(wb, name, data) {
  addWorksheet(wb, name)
  writeData(wb, name, data,
    headerStyle = createStyle(textDecoration = "bold", fgFill = "#DCE6F1"))
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
}

# 1. Load DIA-NN protein summary

raw <- read_excel(RAW_XLSX)

# Map DIA-NN annotation columns to a stable schema.
annotation <- raw[, ANNOT_COLS_RAW] |>
  rename(uniprot_id    = `Protein.Group`,
         protein       = `Protein.Names`,
         gene          = `Genes`,
         description   = `First.Protein.Description`,
         n_seq         = `N.Sequences`,
         n_proteotypic = `N.Proteotypic.Sequences`)

# Backfill missing gene symbols (23 proteins) with the accession.
na_gene <- is.na(annotation$gene) | annotation$gene == ""
annotation$gene[na_gene] <- annotation$uniprot_id[na_gene]
message(sprintf("Gene symbol backfilled from accession for %d proteins", sum(na_gene)))

intensity <- raw[, setdiff(names(raw), ANNOT_COLS_RAW)]
names(intensity) <- h9c2_strip_raw(names(intensity))

metadata <- load_h9c2_metadata(META_CSV)
assert_h9c2_group_sizes(metadata, min_n = 3L)
stopifnot("Sample mismatch: intensity columns vs metadata Col_ID" =
            setequal(colnames(intensity), metadata$Col_ID))
intensity <- intensity[, metadata$Col_ID]

n_raw <- nrow(annotation)
filter_log <- tibble(step = "Raw input", n_before = NA_integer_,
                     n_after = n_raw, n_removed = NA_integer_)
message(sprintf("Raw: %d proteins x %d samples", n_raw, ncol(intensity)))

# 2. Contaminant removal (keratins + FBS)
# Patterns and gene list live in 00_input/h9c2_design.R.

n_before <- nrow(annotation)
is_keratin <- grepl(H9C2_KERATIN_GENE_PATTERN, annotation$gene, ignore.case = TRUE) |
              grepl(H9C2_KERATIN_DESC_PATTERN, annotation$description, ignore.case = TRUE)
is_serum   <- toupper(annotation$gene) %in% toupper(H9C2_SERUM_CONTAMINANTS)
keep_contam <- !(is_keratin | is_serum)

removed_contam <- annotation[!keep_contam, c("uniprot_id", "gene", "description")] |>
  mutate(removal_step = if_else(is_keratin[!keep_contam],
                                "Keratin contaminant", "FBS serum contaminant"))

intensity  <- intensity[keep_contam, ]
annotation <- annotation[keep_contam, ]

filter_log <- bind_rows(filter_log, tibble(
  step = "Contaminant removal (keratin + FBS)", n_before = n_before,
  n_after = nrow(annotation), n_removed = n_before - nrow(annotation)))
message(sprintf("Contaminants removed: %d keratin + %d FBS serum",
                sum(is_keratin), sum(is_serum & !is_keratin)))

# 3. Deduplicate by accession
# Guard only — DIA-NN Protein.Group is unique, so this normally no-ops.

if (any(duplicated(annotation$uniprot_id))) {
  n_before <- nrow(annotation)
  annotation$row_mean <- rowMeans(data.matrix(intensity), na.rm = TRUE)
  keep_idx <- annotation |>
    mutate(row_idx = row_number()) |>
    slice_max(row_mean, n = 1, with_ties = FALSE, by = uniprot_id) |>
    pull(row_idx)
  annotation <- annotation[keep_idx, ]
  intensity  <- intensity[keep_idx, ]
  annotation$row_mean <- NULL
  filter_log <- bind_rows(filter_log, tibble(
    step = "Deduplication", n_before = n_before,
    n_after = nrow(annotation), n_removed = n_before - nrow(annotation)))
}

# 4. Assemble DAList + missingness filter

int_mat <- as.data.frame(data.matrix(intensity))
rownames(int_mat) <- annotation$uniprot_id
annot_df <- as.data.frame(annotation); rownames(annot_df) <- annotation$uniprot_id
meta_df  <- as.data.frame(metadata);   rownames(meta_df)  <- metadata$Col_ID

dal <- DAList(data = int_mat, annotation = annot_df, metadata = meta_df) |>
  zero_to_missing()

n_before <- nrow(dal$data)
dal <- filter_proteins_by_group(dal, min_reps = H9C2_MIN_REPS, min_groups = 1,
                                grouping_column = "Group")

filter_log <- bind_rows(filter_log, tibble(
  step = sprintf("Missingness (>=%d in >=1 group)", H9C2_MIN_REPS),
  n_before = n_before, n_after = nrow(dal$data),
  n_removed = n_before - nrow(dal$data))) |>
  mutate(pct_of_raw = round(n_after / n_raw * 100, 1))

filtered_proteins <- bind_rows(
  removed_contam,
  as_tibble(annot_df) |>
    filter(!uniprot_id %in% rownames(dal$data)) |>
    select(uniprot_id, gene, description) |>
    mutate(removal_step = "Missingness")) |>
  distinct(uniprot_id, .keep_all = TRUE)

# 5. Outlier detection (4-method consensus, >= OUTLIER_K / 4)
# Methods: missingness, PCA-Mahalanobis, MAD intensity, inter-sample correlation.

pct_missing <- colMeans(is.na(dal$data)) * 100

miss_info <- dal$metadata |>
  select(Col_ID, Sample_ID, Group) |>
  mutate(pct_missing = pct_missing[Col_ID])

miss_thresh <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
miss_info$miss_flag <- miss_info$pct_missing > miss_thresh

complete_mat <- dal$data[rowSums(is.na(dal$data)) == 0, ]
pca_pre <- run_pca(complete_mat, dal$metadata, log_transform = TRUE)
pc3     <- pca_pre$pca$x[, 1:3]
mahal   <- mahalanobis(pc3, colMeans(pc3), cov(pc3))
pca_flags <- tibble(Col_ID = colnames(dal$data), mahal_dist = mahal,
                    pca_flag = mahal > qchisq(1 - H9C2_MAHAL_P, df = 3))

samp_med   <- apply(log2(dal$data), 2, median, na.rm = TRUE)
global_med <- median(samp_med)
mad_val    <- mad(samp_med)
mad_flags  <- tibble(Col_ID = names(samp_med), sample_median = samp_med,
                     mad_flag = abs(samp_med - global_med) > H9C2_MAD_K * mad_val)

cor_mat <- cor(log2(dal$data), use = "pairwise.complete.obs")
med_cor <- apply(cor_mat, 2, \(x) median(x[x < 1], na.rm = TRUE))
cor_flags <- tibble(Col_ID = names(med_cor), median_cor = med_cor,
                    cor_flag = med_cor < median(med_cor) - H9C2_MAD_K * mad(med_cor))

outlier_diag <- miss_info |>
  left_join(pca_flags, by = join_by(Col_ID)) |>
  left_join(mad_flags, by = join_by(Col_ID)) |>
  left_join(cor_flags, by = join_by(Col_ID)) |>
  mutate(n_flags = miss_flag + pca_flag + mad_flag + cor_flag,
         consensus_outlier = n_flags >= H9C2_OUTLIER_K)

outlier_ids <- outlier_diag |> filter(consensus_outlier) |> pull(Col_ID)
n_outliers <- length(outlier_ids)
message(sprintf("Outliers: %d flagged (>=%d/4)", n_outliers, H9C2_OUTLIER_K))

data_pre_outlier <- dal$data
meta_pre_outlier <- dal$metadata

if (n_outliers > 0) {
  dal <- filter_samples(dal, !(Col_ID %in% outlier_ids))
  message(sprintf("Removed: %s (%d remain)",
                  paste(outlier_ids, collapse = ", "), ncol(dal$data)))
}

# 6. Normalize (cycloess via proteoDA)

write_norm_report(dal, grouping_column = "Group",
                  output_dir = RPT, filename = "01_norm_comparison.pdf",
                  overwrite = TRUE)
write_qc_report(dal, color_column = "Group",
                output_dir = RPT, filename = "02_qc_pre.pdf", overwrite = TRUE)

dal <- normalize_data(dal, norm_method = "cycloess")

write_qc_report(dal, color_column = "Group",
                output_dir = RPT, filename = "03_qc_post.pdf", overwrite = TRUE)

message(sprintf("Normalized: %d proteins x %d samples",
                nrow(dal$data), ncol(dal$data)))

# 7. Build xlsx

norm_df <- bind_cols(
  as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
  as_tibble(dal$data))

meta_out <- as.data.frame(metadata) |>
  mutate(QC_Status = if_else(Col_ID %in% outlier_ids, "Excluded", "Retained"))

sample_cols <- dal$metadata$Col_ID
prot_miss <- tibble(
  gene     = norm_df$gene,
  n_miss   = rowSums(is.na(norm_df[, sample_cols])),
  pct_miss = round(100 * n_miss / length(sample_cols), 2))

samp_miss <- tibble(
  Col_ID   = sample_cols,
  n_miss   = colSums(is.na(dal$data)),
  pct_miss = round(100 * n_miss / nrow(dal$data), 2)) |>
  left_join(dal$metadata |> select(Col_ID, Sample_ID, Group),
            by = join_by(Col_ID))

wb <- createWorkbook()
write_sheet(wb, "sample_metadata",     meta_out)
write_sheet(wb, "normalized_matrix",   norm_df)
write_sheet(wb, "filter_log",          filter_log)
write_sheet(wb, "outlier_diagnostics", outlier_diag)
write_sheet(wb, "filtered_proteins",   filtered_proteins)
write_sheet(wb, "protein_missingness", prot_miss)
write_sheet(wb, "sample_missingness",  samp_miss)
saveWorkbook(wb, file.path(DAT, "01_normalization.xlsx"), overwrite = TRUE)

# CSV handoff for stages 02-03. Cols 1-4 must stay
# uniprot_id/protein/gene/description — downstream reads the matrix as df[, -(1:4)].
readr::write_csv(norm_df, file.path(DAT, "02_normalized.csv"))

# 8. Save R objects

saveRDS(dal, file.path(DAT, "03_DAList_normalized.rds"))

# Diagnostic intermediates for the report PDF and F00.
filter_bar_data <- filter_log |>
  filter(!is.na(n_removed)) |>
  mutate(step = factor(step, levels = step)) |>
  pivot_longer(c(n_after, n_removed), names_to = "status", values_to = "n") |>
  mutate(status = recode(status, n_after = "Retained", n_removed = "Removed"))

miss_bar_data <- meta_pre_outlier |>
  select(Col_ID, Group) |>
  mutate(detected = colSums(!is.na(data_pre_outlier[, Col_ID])),
         missing  = nrow(data_pre_outlier) - detected,
         is_outlier = Col_ID %in% outlier_ids) |>
  pivot_longer(c(detected, missing), names_to = "status", values_to = "n") |>
  mutate(status = str_to_title(status))

# Per-sample variability (IQR by Group).
samp_var <- dal$metadata |>
  mutate(iqr = apply(dal$data[, Col_ID], 2, IQR, na.rm = TRUE)) |>
  select(Col_ID, Sample_ID, Group, iqr)

grp_vec <- dal$metadata$Group[match(colnames(dal$data), dal$metadata$Col_ID)]
eta2_vals <- apply(dal$data, 1, function(x) {
  ok <- !is.na(x)
  if (sum(ok) < 4) return(NA_real_)
  xk <- x[ok]; gk <- grp_vec[ok]
  ss_b <- sum(tapply(xk, gk, length) * (tapply(xk, gk, mean) - mean(xk))^2)
  ss_t <- sum((xk - mean(xk))^2)
  if (ss_t > 0) ss_b / ss_t else NA_real_
})

pca_post <- run_pca(dal$data, dal$metadata, log_transform = FALSE)

saveRDS(list(
  filter_log        = filter_log,
  filter_bar_data   = filter_bar_data,
  miss_bar_data     = miss_bar_data,
  n_raw             = n_raw,
  n_outliers        = n_outliers,
  outlier_diag      = outlier_diag,
  outlier_ids       = outlier_ids,
  miss_thresh       = miss_thresh,
  pca_pre           = pca_pre,
  pca_post          = pca_post,
  global_med        = global_med,
  mad_val           = mad_val,
  samp_var          = samp_var,
  eta2_vals         = eta2_vals,
  filtered_proteins = filtered_proteins,
  data_pre_outlier  = data_pre_outlier,
  meta_pre_outlier  = meta_pre_outlier,
  dal_nrow          = nrow(dal$data),
  dal_ncol          = ncol(dal$data),
  mahal_p           = H9C2_MAHAL_P,
  mad_k             = H9C2_MAD_K,
  outlier_k         = H9C2_OUTLIER_K,
  pal_group         = H9C2_PAL_GROUP
), file.path(DAT, "00_report_intermediates.rds"))

message(sprintf("Done: %d proteins x %d samples -> %s/",
                nrow(dal$data), ncol(dal$data), DAT))
