#!/usr/bin/env Rscript
# Build the H9c2 fGSEA cache (long: one row per pathway x contrast x database).
# Ranks each contrast by the limma moderated t. Databases come from rat_gene_sets.rds
# (Hallmark, Reactome, KEGG, GO Slim, MitoCarta).
#   Input : 03_DEP/a_non_imputed/c_data/combined_results_pi.csv, shared/c_data/rat_gene_sets.rds
#   Output: 04_Figures/shared/c_data/fgsea_tstat_all_h9c2.csv

pacman::p_load(readr, dplyr, tidyr)
source(here::here("04_Figures", "functions", "shared_gene_set_helpers.R"))

CONTRASTS <- c("CTLvPHE", "CTLvMITO", "PHEvPHE_MITO", "CTLvPHE_MITO", "Interaction", "MITOvPHE_MITO")
OUT <- here::here("04_Figures", "shared", "c_data", "fgsea_tstat_all_h9c2.csv")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
SETS <- readRDS(here::here("04_Figures", "shared", "c_data", "rat_gene_sets.rds"))
DEP <- read_csv(
  here::here("03_DEP", "a_non_imputed", "c_data", "combined_results_pi.csv"),
  show_col_types = FALSE
) |>
  pivot_wider(
    id_cols = c(uniprot_id, gene), names_from = contrast,
    values_from = t, names_glue = "t_{contrast}"
  )
stopifnot(all(paste0("t_", CONTRASTS) %in% names(DEP)))

out <- bind_rows(lapply(names(SETS), \(db) run_fgsea_cache(DEP, SETS[[db]], db, CONTRASTS)))
write_csv(out, OUT)
cat(sprintf(
  "Saved %s: %d rows, %d pathways, %d databases\n",
  OUT, nrow(out), length(unique(out$pathway)), length(unique(out$database))
))
