#!/usr/bin/env Rscript
# Build a gene -> subcellular compartment lookup for the F03/F04 protein scatters
# (ring encoding: Mitochondrial / Nuclear / Cytosolic / Other). Run ONCE; writes
# 04_Figures/shared/protein_localization_rat.csv.
#
# Sources (data-derived, not curated by hand):
#   Mitochondrial = MitoCarta (rat_gene_sets.rds) UNION GO:0005739 (mitochondrion)
#   Nuclear       = GO:0005634 (nucleus)
#   Cytosolic     = GO:0005829 (cytosol)
# GOALL keytype returns genes annotated to a term OR any of its descendants.
# GO "nucleus" is heavily over-annotated, so a naive priority labels most of the
# proteome nuclear and washes out the signal. Assignment is therefore:
#   Mitochondrial  = MitoCarta UNION GO:0005739 (curated + GO; most specific)
#   Nuclear        = nucleus AND NOT cytosol AND NOT mito (compartment-specific)
#   Cytosolic      = cytosol AND NOT nucleus AND NOT mito (compartment-specific)
#   Other          = multi-localised (nucleus+cytosol) or unannotated
# i.e. only confidently single-compartment proteins get a nuclear/cytosolic ring.

setwd(rprojroot::find_rstudio_root_file())
suppressPackageStartupMessages({
  library(org.Rn.eg.db)
  library(AnnotationDbi)
  library(GO.db)
})

GO_IDS <- c(Mitochondrial = "GO:0005739", Nuclear = "GO:0005634", Cytosolic = "GO:0005829")
# self-verify the GO term names so a wrong ID can't pass silently
for (nm in names(GO_IDS))
  message(sprintf("  %-14s %s = %s", nm, GO_IDS[nm], AnnotationDbi::Term(GOTERM[[GO_IDS[nm]]])))

go_symbols <- function(go_id) {
  s <- suppressMessages(AnnotationDbi::select(
    org.Rn.eg.db, keys = go_id, keytype = "GOALL", columns = "SYMBOL"))
  unique(s$SYMBOL[!is.na(s$SYMBOL)])
}
mito_go <- go_symbols(GO_IDS["Mitochondrial"])
nuc     <- go_symbols(GO_IDS["Nuclear"])
cyto    <- go_symbols(GO_IDS["Cytosolic"])

mito_carta <- unique(unlist(readRDS("04_Figures/shared/rat_gene_sets.rds")$MitoCarta,
                            use.names = FALSE))
mito <- union(mito_go, mito_carta)

all_genes <- unique(c(mito, nuc, cyto))
loc <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
loc$localization <- dplyr::case_when(
  loc$gene %in% mito                              ~ "Mitochondrial",
  loc$gene %in% nuc  & !loc$gene %in% cyto        ~ "Nuclear",
  loc$gene %in% cyto & !loc$gene %in% nuc         ~ "Cytosolic",
  TRUE                                            ~ "Other")

readr::write_csv(loc, "04_Figures/shared/protein_localization_rat.csv")
message(sprintf("\nWrote %d genes: %s",
                nrow(loc),
                paste(sprintf("%s=%d", names(table(loc$localization)),
                              as.integer(table(loc$localization))), collapse = " | ")))

# Distribution within the actual DEP dataset (what the scatter will show).
dep <- readr::read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
dep_loc <- merge(data.frame(gene = unique(dep$gene)), loc, all.x = TRUE)
dep_loc$localization[is.na(dep_loc$localization)] <- "Other"
message(sprintf("Within DEP dataset (%d genes): %s", nrow(dep_loc),
                paste(sprintf("%s=%d", names(table(dep_loc$localization)),
                              as.integer(table(dep_loc$localization))), collapse = " | ")))
