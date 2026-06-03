#!/usr/bin/env Rscript
# Additive rat GO Slim fGSEA (for F02 Hallmark + GO Slim ring exploration).
# Does NOT touch the frozen main cache (fgsea_tstat_all_h9c2.csv). Rat-specific:
# GO Slim gene sets are built via org.Rn.eg.db (Rattus norvegicus). Ranking and
# fgsea params mirror build_fgsea_cache.R so results are comparable.
#   Inputs : 03_DEP/c_data/03_combined_results.csv  (moderated t per contrast)
#   Outputs: 04_Figures/shared/fgsea_goslim_h9c2.csv  (cache-schema rows)
#            04_Figures/shared/goslim_rat_gene_sets.rds (sets, for ring dedup)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(fgsea)
})

source(here::here("00_input", "h9c2_design.R"))
source(here::here("04_Figures", "shared", "pathway_utils.R"))

set.seed(42)

OUT_CSV <- here::here("04_Figures", "shared", "fgsea_goslim_h9c2.csv")
OUT_RDS <- here::here("04_Figures", "shared", "goslim_rat_gene_sets.rds")
DEP <- read_csv(here::here("03_DEP", "c_data", "03_combined_results.csv"),
                show_col_types = FALSE)

# Rat (org.Rn.eg.db) GO Slim Generic BP gene sets.
goslim <- build_goslim_gene_sets(species = "Rattus norvegicus",
                                 min_size = 10, max_size = 500)
saveRDS(goslim, OUT_RDS)

all_contrasts <- c(H9C2_CORE_CONTRASTS, "MITOvPHE_MITO")

run_one <- function(contrast) {
  t_col <- paste0("t_", contrast)
  stats <- DEP[[t_col]]
  names(stats) <- DEP$gene
  stats <- stats[!is.na(stats) & !is.na(names(stats)) & names(stats) != ""]
  if (anyDuplicated(names(stats))) stats <- tapply(stats, names(stats), mean)
  stats <- sort(stats)
  res <- fgsea(pathways = goslim, stats = stats, minSize = 10, maxSize = 500, eps = 0)
  if (nrow(res) == 0) return(NULL)
  res |>
    mutate(database = "GO Slim", contrast = contrast,
           leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";")) |>
    select(pathway, pval, padj, log2err, ES, NES, size, leadingEdge, database, contrast)
}

out <- bind_rows(lapply(all_contrasts, run_one))
write_csv(out, OUT_CSV)
for (ctr in all_contrasts) {
  sub <- out[out$contrast == ctr, ]
  cat(sprintf("  %-15s | GO Slim : %3d sets, %2d with FDR<0.05\n",
              ctr, nrow(sub), sum(sub$padj < 0.05, na.rm = TRUE)))
}
cat(sprintf("\nSaved %s (%d rows) + %s (%d sets)\n",
            OUT_CSV, nrow(out), OUT_RDS, length(goslim)))
