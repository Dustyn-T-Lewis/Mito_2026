#!/usr/bin/env Rscript
# =============================================================================
# 01_run_filtering.R  --  Mito Stage 01: Filtering (proteoDA-native)
#
# H9c2 rat cardiac cell line, 4 groups x 6 replicates (Ctl/Mito/PHE/PHE_Mito).
# load -> contaminant removal (keratin regex + reference-derived FBS serum panel)
# -> dedup -> build DAList -> native filter_proteins_by_group (missingness) ->
# 4-method outlier consensus -> filter_samples. Un-normalized output for Stage 02.
#
# Contaminants are reference-derived (R/contaminant_filter.R), not hand-typed:
# keratins by regex; FBS serum = bovine Hao cell-culture contaminants that are
# HPA "Secreted to blood" (excludes abundant-cellular bovine proteins = real biology).
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(readxl); library(readr)
  library(dplyr); library(tidyr); library(stringr)
})
set.seed(42)
source(here("00_input", "h9c2_design.R"))      # h9c2_strip_raw, load_h9c2_metadata, H9C2_* constants
source(here("R", "contaminant_filter.R"))

ANNOT_COLS_RAW <- c("Protein.Group", "Protein.Names", "Genes",
                    "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences")
cfg <- list(
  raw = here("00_input", "H9c2_raw.xlsx"), meta = here("00_input", "H9c2_meta.csv"),
  hao = here("00_input", "hao_cellculture_contaminants.fasta"),
  hpa_secreted = here("00_input", "hpa_secreted_to_blood.txt"),
  data_dir = here("01_Filtering", "c_data"),
  min_reps = H9C2_MIN_REPS, outlier_k = H9C2_OUTLIER_K, mad_k = H9C2_MAD_K, mahal_p = H9C2_MAHAL_P
)
dir.create(cfg$data_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load + metadata ------------------------------------------------------
raw <- read_excel(cfg$raw)
annotation <- raw[, ANNOT_COLS_RAW] |>
  rename(uniprot_id = `Protein.Group`, protein = `Protein.Names`, gene = `Genes`,
         description = `First.Protein.Description`, n_seq = `N.Sequences`,
         n_proteotypic = `N.Proteotypic.Sequences`)
na_gene <- is.na(annotation$gene) | annotation$gene == ""
annotation$gene[na_gene] <- annotation$uniprot_id[na_gene]
intensity <- raw[, setdiff(names(raw), ANNOT_COLS_RAW)]; names(intensity) <- h9c2_strip_raw(names(intensity))
metadata <- load_h9c2_metadata(cfg$meta); metadata$group <- metadata$Group   # lowercase for formula
stopifnot(setequal(colnames(intensity), metadata$Col_ID)); intensity <- intensity[, metadata$Col_ID]
n_raw <- nrow(annotation)

# --- 2. Contaminant removal (keratin regex + FBS serum panel) ----------------
serum <- serum_contaminant_genes(cfg$hao, cfg$hpa_secreted)
flag  <- flag_contaminants(annotation$gene, annotation$description, serum)
flog  <- tibble(step = "Raw input", n_after = n_raw, n_removed = NA_integer_)
flog  <- bind_rows(flog, tibble(step = "Contaminant removal (keratin + FBS serum)",
                                n_after = sum(!flag$is_contaminant), n_removed = sum(flag$is_contaminant)))
write_csv(annotation[flag$is_contaminant, c("uniprot_id", "gene", "description")] |>
            mutate(reason = flag$reason[flag$is_contaminant]),
          file.path(cfg$data_dir, "04_contaminants_removed.csv"))
annotation <- annotation[!flag$is_contaminant, ]; intensity <- intensity[!flag$is_contaminant, ]
cat(sprintf("Contaminants removed: %d keratin + %d FBS serum\n",
            sum(flag$reason == "keratin", na.rm = TRUE), sum(grepl("serum", flag$reason), na.rm = TRUE)))

if (any(duplicated(annotation$uniprot_id))) {            # guard (DIA-NN groups are unique)
  rm_mean <- rowMeans(data.matrix(intensity), na.rm = TRUE)
  keep_idx <- tibble(i = seq_along(rm_mean), id = annotation$uniprot_id, m = rm_mean) |>
    group_by(id) |> slice_max(m, n = 1, with_ties = FALSE) |> pull(i)
  annotation <- annotation[keep_idx, ]; intensity <- intensity[keep_idx, ]
}

# --- 3. Build DAList + native missingness filter -----------------------------
int_mat <- as.data.frame(data.matrix(intensity)); rownames(int_mat) <- annotation$uniprot_id
annot_df <- as.data.frame(annotation); rownames(annot_df) <- annotation$uniprot_id
meta_df  <- as.data.frame(metadata);   rownames(meta_df)  <- metadata$Col_ID
dal <- zero_to_missing(DAList(data = int_mat, annotation = annot_df, metadata = meta_df))
n0 <- nrow(dal$data)
dal <- filter_proteins_by_group(dal, min_reps = cfg$min_reps, min_groups = 1, grouping_column = "group")
flog <- bind_rows(flog, tibble(step = sprintf("Missingness (>=%d in >=1 group)", cfg$min_reps),
                               n_after = nrow(dal$data), n_removed = n0 - nrow(dal$data)))
flog <- flog |> mutate(pct_of_raw = round(n_after / n_raw * 100, 1))

# --- 4. Outlier consensus (4-method, >=K/4) -> native filter_samples ---------
pct_missing <- colMeans(is.na(dal$data)) * 100
miss_flag <- pct_missing > quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
sm <- apply(log2(dal$data), 2, median, na.rm = TRUE); mad_flag <- abs(sm - median(sm)) > cfg$mad_k * mad(sm)
pc3 <- prcomp(t(log2(dal$data[rowSums(is.na(dal$data)) == 0, ])), center = TRUE, scale. = TRUE)$x[, 1:3]
pca_flag <- mahalanobis(pc3, colMeans(pc3), cov(pc3)) > qchisq(1 - cfg$mahal_p, df = 3)
cmx <- cor(log2(dal$data), use = "pairwise.complete.obs"); mc <- apply(cmx, 2, function(x) median(x[x < 1], na.rm = TRUE))
cor_flag <- mc < median(mc) - cfg$mad_k * mad(mc)
outlier_diag <- tibble(Col_ID = colnames(dal$data), miss_flag = miss_flag[Col_ID], mad_flag = mad_flag[Col_ID],
                       pca_flag = pca_flag[Col_ID], cor_flag = cor_flag[Col_ID]) |>
  mutate(n_flags = miss_flag + mad_flag + pca_flag + cor_flag, consensus_outlier = n_flags >= cfg$outlier_k)
outlier_ids <- outlier_diag$Col_ID[outlier_diag$consensus_outlier]
cat(sprintf("Outliers (>=%d/4): %s\n", cfg$outlier_k, if (length(outlier_ids)) paste(outlier_ids, collapse = ", ") else "none"))
if (length(outlier_ids)) dal <- filter_samples(dal, !(Col_ID %in% outlier_ids))

# --- 5. Export ---------------------------------------------------------------
saveRDS(dal, file.path(cfg$data_dir, "01_DAList_filtered.rds"))
write_csv(bind_cols(as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
                    as_tibble(dal$data)), file.path(cfg$data_dir, "02_filtered_matrix.csv"))
write_csv(flog, file.path(cfg$data_dir, "03_filter_log.csv"))
saveRDS(list(filter_log = flog, outlier_diag = outlier_diag, serum_panel = serum, n_raw = n_raw),
        file.path(cfg$data_dir, "00_filter_intermediates.rds"))
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
cat(sprintf("Done: %d proteins x %d samples -> %s/\n", nrow(dal$data), ncol(dal$data), cfg$data_dir))
print(as.data.frame(flog))
