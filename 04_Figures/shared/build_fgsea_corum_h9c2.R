#!/usr/bin/env Rscript
# Additive rat CORUM protein-complex fGSEA (for F02 CORUM ring composite).
# Does NOT touch the frozen main cache. Source: CORUM allComplexes (Human +
# Mouse complexes), subunit gene symbols mapped to RAT via babelgene
# (consistent with fetch_rat_gene_sets.R: Human→rat direct; Mouse→human→rat).
# Rat-native CORUM complexes are skipped to avoid symbol-case ambiguity; Human
# (n≈5638) + Mouse (n≈2245) dominate coverage anyway.
#   Inputs : /Users/dtl0018/Documents/corum_allComplexes.txt
#            03_DEP/c_data/03_combined_results.csv (moderated t per contrast)
#   Outputs: 04_Figures/shared/fgsea_corum_h9c2.csv  (cache-schema rows)
#            04_Figures/shared/corum_rat_gene_sets.rds (sets, for ring dedup)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(fgsea); library(babelgene)
})

source(here::here("00_input", "h9c2_design.R"))
set.seed(42)

CORUM_TXT <- "/Users/dtl0018/Documents/corum_allComplexes.txt"
OUT_CSV   <- here::here("04_Figures", "shared", "fgsea_corum_h9c2.csv")
OUT_RDS   <- here::here("04_Figures", "shared", "corum_rat_gene_sets.rds")
DEP <- read_csv(here::here("03_DEP", "c_data", "03_combined_results.csv"), show_col_types = FALSE)

corum <- read.delim(CORUM_TXT, sep = "\t", quote = "", stringsAsFactors = FALSE) |>
  filter(organism %in% c("Human", "Mouse"),
         !is.na(subunits_gene_name), nzchar(subunits_gene_name))

# Ortholog maps to rat
split_genes <- function(x) unique(unlist(strsplit(x, ";")))
human_genes <- split_genes(corum$subunits_gene_name[corum$organism == "Human"])
mouse_genes <- split_genes(corum$subunits_gene_name[corum$organism == "Mouse"])

# Human → rat (input genes ARE human → human = TRUE; `symbol` is rat).
h2r <- babelgene::orthologs(genes = human_genes, species = "rat", human = TRUE) |>
  distinct(human_symbol, .keep_all = TRUE)
h2r_map <- setNames(h2r$symbol, h2r$human_symbol)

# Mouse → human (human = FALSE; `symbol` is mouse) → human → rat.
m2h <- babelgene::orthologs(genes = mouse_genes, species = "mouse", human = FALSE) |>
  distinct(symbol, .keep_all = TRUE)
mh2r <- babelgene::orthologs(genes = unique(m2h$human_symbol), species = "rat", human = TRUE) |>
  distinct(human_symbol, .keep_all = TRUE)
m2r_map <- setNames(mh2r$symbol[match(m2h$human_symbol, mh2r$human_symbol)], m2h$symbol)
m2r_map <- m2r_map[!is.na(m2r_map)]

to_rat <- function(genes, org) {
  m <- if (org == "Human") h2r_map else m2r_map
  unique(unname(m[genes[genes %in% names(m)]]))
}

# Build complex gene sets (rat)
MIN_SUBUNITS <- 3
sets <- list()
for (i in seq_len(nrow(corum))) {
  g <- to_rat(split_genes(corum$subunits_gene_name[i]), corum$organism[i])
  if (length(g) < MIN_SUBUNITS) next
  nm <- str_squish(corum$complex_name[i])
  if (!nzchar(nm)) next
  key <- nm
  k <- 1L; while (key %in% names(sets)) { k <- k + 1L; key <- paste0(nm, " (", k, ")") }
  sets[[key]] <- g
}
# Drop complexes with identical rat membership (ortholog collapse of paralogous
# human/mouse complexes); keep the first (largest source name).
memb_key <- vapply(sets, function(g) paste(sort(g), collapse = "|"), character(1))
sets <- sets[!duplicated(memb_key)]
saveRDS(sets, OUT_RDS)
cat(sprintf("CORUM: %d complexes -> %d rat sets (>=%d subunits, dedup membership)\n",
            nrow(corum), length(sets), MIN_SUBUNITS))

# fGSEA per core contrast (t-stat ranks; small minSize for complexes)
all_contrasts <- c(H9C2_CORE_CONTRASTS, "MITOvPHE_MITO")
run_one <- function(contrast) {
  stats <- DEP[[paste0("t_", contrast)]]; names(stats) <- DEP$gene
  stats <- stats[!is.na(stats) & !is.na(names(stats)) & names(stats) != ""]
  if (anyDuplicated(names(stats))) stats <- tapply(stats, names(stats), mean)
  stats <- sort(stats)
  res <- fgsea(pathways = sets, stats = stats, minSize = 3, maxSize = 500, eps = 0)
  if (nrow(res) == 0) return(NULL)
  res |>
    mutate(database = "CORUM", contrast = contrast,
           leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";")) |>
    select(pathway, pval, padj, log2err, ES, NES, size, leadingEdge, database, contrast)
}
out <- bind_rows(lapply(all_contrasts, run_one))
write_csv(out, OUT_CSV)
for (ctr in all_contrasts) {
  sub <- out[out$contrast == ctr, ]
  cat(sprintf("  %-15s | CORUM : %4d sets, %3d with FDR<0.05\n",
              ctr, nrow(sub), sum(sub$padj < 0.05, na.rm = TRUE)))
}
cat(sprintf("\nSaved %s (%d rows) + %s (%d sets)\n", OUT_CSV, nrow(out), OUT_RDS, length(sets)))
