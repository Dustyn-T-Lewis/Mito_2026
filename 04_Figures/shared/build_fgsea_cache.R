#!/usr/bin/env Rscript
# Build the H9c2 fGSEA cache (long-format) — one row per pathway-contrast-database
# Matches the YvO cache schema (fgsea_tstat_all_v2.csv).
# Inputs : 03_DEP/c_data/03_combined_results.csv, shared/rat_gene_sets.rds
# Output : 04_Figures/shared/fgsea_tstat_all_h9c2.csv

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr); library(fgsea)
})

source(here::here("00_input", "h9c2_design.R"))

set.seed(42)

OUT  <- here::here("04_Figures", "shared", "fgsea_tstat_all_h9c2.csv")
SETS <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))
DEP  <- read_csv(here::here("03_DEP", "c_data", "03_combined_results.csv"),
                 show_col_types = FALSE)

# Run fgsea on every contrast x database, ranking by the moderated t-statistic.
all_contrasts <- c(H9C2_CORE_CONTRASTS, "MITOvPHE_MITO")

run_one <- function(contrast, db_name, sets) {
  t_col <- paste0("t_", contrast)
  stats <- DEP[[t_col]]
  names(stats) <- DEP$gene
  stats <- stats[!is.na(stats) & !is.na(names(stats)) & names(stats) != ""]
  # Average duplicate gene symbols (DIA-NN occasionally has duplicates)
  if (anyDuplicated(names(stats))) {
    stats <- tapply(stats, names(stats), mean)
  }
  stats <- sort(stats)
  res <- fgsea(pathways = sets, stats = stats,
               minSize  = 10, maxSize = 500, eps = 0)
  if (nrow(res) == 0) return(NULL)
  res |>
    mutate(database = db_name, contrast = contrast,
           leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";")) |>
    select(pathway, pval, padj, log2err, ES, NES, size, leadingEdge,
           database, contrast)
}

cache <- list()
for (ctr in all_contrasts) {
  for (db in names(SETS)) {
    key <- paste(ctr, db, sep = "::")
    cache[[key]] <- run_one(ctr, db, SETS[[db]])
    n_sig <- if (!is.null(cache[[key]])) sum(cache[[key]]$padj < 0.05, na.rm = TRUE) else 0
    n_tot <- if (!is.null(cache[[key]])) nrow(cache[[key]]) else 0
    cat(sprintf("  %-15s | %-10s : %4d pathways, %3d with FDR<0.05\n",
                ctr, db, n_tot, n_sig))
  }
}

out <- bind_rows(cache)
write_csv(out, OUT)
cat(sprintf("\nSaved: %s\n  rows: %d, pathways: %d, contrasts: %d, databases: %d\n",
            OUT, nrow(out), length(unique(out$pathway)),
            length(unique(out$contrast)), length(unique(out$database))))
