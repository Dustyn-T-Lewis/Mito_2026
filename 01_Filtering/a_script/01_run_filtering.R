#!/usr/bin/env Rscript

pacman::p_load(proteoDA, here, readxl, readr, dplyr, tidyr, stringr, openxlsx, ggplot2)

group_levels <- c("Ctl", "Mito", "PHE", "PHE_Mito") # factor order used everywhere downstream
annot_cols_raw <- c(
  "Protein.Group", "Protein.Names", "Genes",
  "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences"
)
data_dir <- here("01_Filtering", "c_data") # outputs; cleared here, written at the end
report_dir <- here("01_Filtering", "b_reports")

# stage owns its outputs, so wipe stale runs
clear_dir <- function(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(d, full.names = TRUE), recursive = TRUE)
}
clear_dir(data_dir)
clear_dir(report_dir)

raw <- read_excel(here("00_input", "H9c2_raw.xlsx"))
annotation <- raw[, annot_cols_raw] |>
  rename(
    uniprot_id = `Protein.Group`, protein = `Protein.Names`, gene = `Genes`,
    description = `First.Protein.Description`, n_seq = `N.Sequences`,
    n_proteotypic = `N.Proteotypic.Sequences`
  )
na_gene <- is.na(annotation$gene) | annotation$gene == "" # fall back to accession so no row is label-less
annotation$gene[na_gene] <- annotation$uniprot_id[na_gene]
intensity <- raw[, setdiff(names(raw), annot_cols_raw)]
names(intensity) <- sub("\\.raw$", "", names(intensity))

metadata <- as.data.frame(read_csv(here("00_input", "H9c2_meta.csv"), show_col_types = FALSE))
stopifnot(setequal(unique(metadata$Group), group_levels), !anyDuplicated(metadata$Col_ID))
metadata$Group <- factor(metadata$Group, levels = group_levels)
metadata$group <- metadata$Group # lowercase alias the model formula reads
rownames(metadata) <- metadata$Col_ID
stopifnot(setequal(colnames(intensity), metadata$Col_ID))
intensity <- intensity[, metadata$Col_ID]
n_raw <- nrow(annotation)

# FBS serum panel is reference-derived, not hand-typed: bovine entries in the Hao
# cell-culture FASTA (Frankenfield 2022) intersected with HPA "Secreted to blood".
# The search DB is rat-only, so no bovine protein is ever identified as such: this
# drops rat groups whose gene symbol matches a bovine secreted entry, an orthology
# proxy for FBS peptides landing on conserved rat sequence. Deliberately over-broad.
# Of the 19 it takes, Alb/Tf/Ltf/C3/A2m are the plasma proteins it is aimed at;
# Ppia/Lgals1/Gsn/Fn1 are cellular and go with them.
hao_hdr <- grep("^>", readLines(here("00_input", "hao_cellculture_contaminants.fasta"), warn = FALSE), value = TRUE)
bovine <- toupper(na.omit(str_match(hao_hdr, "GN=([^ ]+)")[, 2][
  grepl("Bos taurus", str_match(hao_hdr, "OS=(.*?) OX=")[, 2])
]))
secreted <- toupper(readLines(here("00_input", "hpa_secreted_to_blood.txt"), warn = FALSE))
serum <- sort(intersect(unique(bovine), secreted))

is_keratin <- grepl("^Krt[0-9]|^Krtap", annotation$gene, ignore.case = TRUE) |
  grepl("^Keratin", annotation$description, ignore.case = TRUE)
gene_tokens <- strsplit(toupper(annotation$gene), ";", fixed = TRUE) # DIA-NN joins grouped genes with ";"
is_serum <- vapply(gene_tokens, function(g) any(g %in% serum), logical(1)) & !is_keratin
reason <- ifelse(is_keratin, "keratin", ifelse(is_serum, "FBS serum (bovine plasma)", NA_character_))
is_contaminant <- is_keratin | is_serum

flog <- tibble(step = "Raw input", n_after = n_raw, n_removed = NA_integer_)
flog <- bind_rows(flog, tibble(
  step = "Contaminant removal (keratin + FBS serum)",
  n_after = sum(!is_contaminant), n_removed = sum(is_contaminant)
))
removed_contam <- annotation[is_contaminant, c("uniprot_id", "gene", "description")] |>
  mutate(reason = reason[is_contaminant])
annotation <- annotation[!is_contaminant, ]
intensity <- intensity[!is_contaminant, ]
cat(sprintf("Contaminants removed: %d keratin + %d FBS serum\n", sum(is_keratin), sum(is_serum)))

if (any(duplicated(annotation$uniprot_id))) { # guard only; DIA-NN protein groups are unique
  rm_mean <- rowMeans(data.matrix(intensity), na.rm = TRUE)
  keep_idx <- tibble(i = seq_along(rm_mean), id = annotation$uniprot_id, m = rm_mean) |>
    group_by(id) |>
    slice_max(m, n = 1, with_ties = FALSE) |>
    pull(i)
  annotation <- annotation[keep_idx, ]
  intensity <- intensity[keep_idx, ]
}

int_mat <- as.data.frame(data.matrix(intensity))
rownames(int_mat) <- annotation$uniprot_id
annot_df <- as.data.frame(annotation)
rownames(annot_df) <- annotation$uniprot_id
meta_df <- as.data.frame(metadata)
rownames(meta_df) <- metadata$Col_ID
dal <- zero_to_missing(DAList(data = int_mat, annotation = annot_df, metadata = meta_df)) # DIA-NN 0 = non-detection

n0 <- nrow(dal$data)
# keep a protein quantified in >= 4 of 6 reps in at least one group (detectable somewhere)
dal <- filter_proteins_by_group(dal, min_reps = 4, min_groups = 1, grouping_column = "group")
flog <- bind_rows(flog, tibble(
  step = "Missingness (>=4 in >=1 group)",
  n_after = nrow(dal$data), n_removed = n0 - nrow(dal$data)
))
flog <- flog |> mutate(pct_of_raw = round(n_after / n_raw * 100, 1))

# Four independent QC heuristics; a sample is dropped only on >= 3/4 agreement.
pct_missing <- colMeans(is.na(dal$data)) * 100
miss_flag <- pct_missing > quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing) # Tukey upper fence
sm <- apply(log2(dal$data), 2, median, na.rm = TRUE)
mad_flag <- abs(sm - median(sm)) > 3 * mad(sm) # Hampel 3xMAD on median intensity
pc3 <- prcomp(t(log2(dal$data[rowSums(is.na(dal$data)) == 0, ])), center = TRUE, scale. = TRUE)$x[, 1:3]
pca_flag <- mahalanobis(pc3, colMeans(pc3), cov(pc3)) > qchisq(0.99, df = 3) # chi-sq tail on PC1-3 position
cmx <- cor(log2(dal$data), use = "pairwise.complete.obs")
mc <- apply(cmx, 2, function(x) median(x[x < 1], na.rm = TRUE)) # x<1 drops the self-correlation
cor_flag <- mc < median(mc) - 3 * mad(mc)
outlier_diag <- tibble(
  Col_ID = colnames(dal$data), miss_flag = miss_flag[Col_ID], mad_flag = mad_flag[Col_ID],
  pca_flag = pca_flag[Col_ID], cor_flag = cor_flag[Col_ID]
) |>
  mutate(n_flags = miss_flag + mad_flag + pca_flag + cor_flag, consensus_outlier = n_flags >= 3)
outlier_ids <- outlier_diag$Col_ID[outlier_diag$consensus_outlier]
cat(sprintf("Outliers (>=3/4): %s\n", if (length(outlier_ids)) paste(outlier_ids, collapse = ", ") else "none"))
if (length(outlier_ids)) dal <- filter_samples(dal, !(Col_ID %in% outlier_ids))

saveRDS(dal, file.path(data_dir, "DAList_filtered.rds")) # un-normalized handoff to Stage 02
write.xlsx(
  list(
    filter_log = flog,
    contaminants_removed = removed_contam,
    outlier_diagnostics = outlier_diag,
    serum_panel = tibble(gene = serum)
  ),
  file.path(data_dir, "filtering_report.xlsx"),
  overwrite = TRUE
)

ggsave(file.path(report_dir, "contaminants.pdf"),
  count(removed_contam, reason) |> ggplot(aes(reorder(reason, n), n, fill = reason)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = n), hjust = -0.2, size = 4) +
    scale_fill_manual(values = c("keratin" = "#4393C3", "FBS serum (bovine plasma)" = "#D6604D")) +
    coord_flip(clip = "off") +
    labs(title = "Contaminants filtered out", x = NULL, y = "proteins") +
    theme_minimal(base_size = 11),
  width = 7, height = 3.2
)
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
cat(sprintf("Done: %d proteins x %d samples -> %s/\n", nrow(dal$data), ncol(dal$data), data_dir))
print(as.data.frame(flog))
