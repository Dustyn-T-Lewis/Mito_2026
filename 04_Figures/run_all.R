#!/usr/bin/env Rscript
# Build every figure in dependency order. Each main script creates its figure's
# single workbook; the supplement scripts append their sheets to it, so run them
# after their main. Assumes shared/c_data/rat_gene_sets.rds already exists.

root <- rprojroot::find_root(rprojroot::has_dir("04_Figures"))

scripts <- c(
  "shared/a_script/02_build_fgsea_cache.R",
  "F01_Proteome_Overview/a_script/01_proteome_overview.R",
  "F01_Proteome_Overview/a_script/supplementary/ora_supplement.R",
  "F02_Enrich_Volcanoes/a_script/01_enrich_volcanoes.R",
  "F02_Enrich_Volcanoes/a_script/supplementary/gsea_supplement.R",
  "F03_Clustering/a_script/01_clustering.R"
)

for (s in scripts) {
  message("\n----- ", s, " -----")
  if (system2("Rscript", shQuote(file.path(root, "04_Figures", s))) != 0L) {
    stop("failed: ", s, call. = FALSE)
  }
}
message("\nall figures built")
