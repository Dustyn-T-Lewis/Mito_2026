#!/usr/bin/env Rscript
# Build the unified rat gene-set collection every enrichment panel consumes.
#   Output: 04_Figures/shared/c_data/rat_gene_sets.rds
#
# Five databases, each deduplicated within itself by Jaccard >= 0.5 (larger set
# kept) so the tested universe carries no near-duplicate pairs. Overlap-coefficient
# (containment) removal is deliberately NOT applied here: at collection scope it
# collapses the curated general-superset-of-specific hierarchies (every MitoCarta
# sub-set sits inside MITOCARTA_ALL), erasing the specific biology. Containment
# redundancy is handled at runtime in the dedup engine, where both sets being
# significant makes it a genuine duplicate finding (Merico 2010 PMID:21085593;
# Reimand 2019 PMID:30664679).
#
#   Hallmark / KEGG (legacy) / Reactome  -- msigdbr, Rattus norvegicus
#   GO Slim                              -- GO Slim Generic BP via org.Rn.eg.db
#   MitoCarta                            -- MitoCarta3.0 rat, tracked source rds
#
# PROVENANCE: MSigDB sets are pinned to the installed msigdbr (recorded in the rds
# attributes). The MitoCarta3.0 rat block is read from the tracked source
# mitocarta3_rat.rds (the Broad Rat.MitoCarta3.0.xls is not in this repo).
# See PROVENANCE.md.

suppressPackageStartupMessages(library(msigdbr))
source(here::here("04_Figures", "shared", "a_script", "00_gene_set_helpers.R"))

SP <- "Rattus norvegicus"
OUT <- here::here("04_Figures", "shared", "c_data", "rat_gene_sets.rds")

# Within-database near-duplicate collapse, Jaccard only. Largest set first so the
# kept entry is the representative; a set is dropped when its Jaccard to an
# already-kept set reaches the cutoff.
dedup_within <- function(sets, jaccard_cutoff = 0.5) {
  sets <- sets[order(lengths(sets), decreasing = TRUE)]
  kept <- vector("list", length(sets))
  keep <- logical(length(sets))
  n_kept <- 0L
  for (i in seq_along(sets)) {
    g <- sets[[i]]
    redundant <- FALSE
    for (j in seq_len(n_kept)) {
      inter <- length(intersect(g, kept[[j]]))
      if (inter == 0L) next
      if (inter / (length(g) + length(kept[[j]]) - inter) >= jaccard_cutoff) {
        redundant <- TRUE
        break
      }
    }
    if (!redundant) {
      keep[i] <- TRUE
      n_kept <- n_kept + 1L
      kept[[n_kept]] <- g
    }
  }
  sets[keep]
}

msig_block <- function(collection, sub = NULL) {
  d <- msigdbr(species = SP, collection = collection, subcollection = sub)
  lapply(split(d$gene_symbol, d$gs_name), unique)
}

hallmark <- msig_block("H")
kegg <- msig_block("C2", "CP:KEGG_LEGACY") # LEGACY has rat orthologs; MEDICUS is human-only
reactome <- msig_block("C2", "CP:REACTOME")
goslim <- build_goslim_gene_sets(species = SP, min_size = 10, max_size = 500)
mitocarta <- readRDS(here::here("04_Figures", "shared", "mitocarta3_rat.rds"))

raw_n <- c(
  Hallmark = length(hallmark), Reactome = length(reactome),
  KEGG = length(kegg), `GO Slim` = length(goslim),
  MitoCarta = length(mitocarta)
)

sets <- list(
  Hallmark = dedup_within(hallmark),
  Reactome = dedup_within(reactome),
  KEGG = dedup_within(kegg),
  `GO Slim` = dedup_within(goslim),
  MitoCarta = dedup_within(mitocarta)
)

attr(sets, "provenance") <- list(
  msigdbr_version = as.character(packageVersion("msigdbr")),
  species = SP, built = as.character(Sys.Date()),
  dedup = "within-DB Jaccard>=0.5, larger kept",
  mitocarta = "MitoCarta3.0 rat from tracked mitocarta3_rat.rds"
)

saveRDS(sets, OUT)

kept_n <- vapply(sets, length, integer(1))
for (db in names(sets)) {
  cat(sprintf(
    "%-10s %4d -> %4d sets (dedup removed %d)\n",
    db, raw_n[[db]], kept_n[[db]], raw_n[[db]] - kept_n[[db]]
  ))
}
cat(sprintf("Saved %s (msigdbr %s)\n", OUT, packageVersion("msigdbr")))
