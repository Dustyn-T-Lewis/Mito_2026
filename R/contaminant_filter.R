# =============================================================================
# contaminant_filter.R  --  exogenous contaminant flags for H9c2 (rat cardiac cell line)
#
# H9c2 is cultured in fetal bovine serum (FBS); the exogenous contaminants are
# keratins (skin/handling) and FBS bovine serum proteins. Both are derived from
# references, not hand-typed:
#   * Keratins  -> regex (^Krt..., ^Keratin description)
#   * FBS serum -> BOVINE entries in the Hao cell-culture contaminant FASTA
#                  (Frankenfield 2022) that are ALSO HPA "Secreted to blood"
#                  plasma proteins. Intersecting with the secretome EXCLUDES the
#                  abundant-cellular bovine proteins (GAPDH, actin, tubulin, HSPs)
#                  that are genuine H9c2 biology, keeping only true plasma proteins.
#
# Generic MS contaminant libraries (cRAP / Hao universal) are NOT subtracted
# wholesale: they list abundant cellular proteins (myoglobin, actin, GAPDH) that
# are real cardiac signal. Tissue-RNA references over-remove and miss keratins.
# =============================================================================

suppressPackageStartupMessages({ library(readr); library(stringr) })

KERATIN_GENE_PATTERN <- "^Krt[0-9]|^Krtap"
KERATIN_DESC_PATTERN <- "^Keratin"

# FBS serum panel = bovine Hao cell-culture contaminants intersected with HPA
# "Secreted to blood". Returns UPPER-CASE gene symbols.
serum_contaminant_genes <- function(hao_fasta, hpa_secreted_file) {
  h  <- grep("^>", readLines(hao_fasta, warn = FALSE), value = TRUE)
  os <- str_match(h, "OS=(.*?) OX=")[, 2]
  gn <- str_match(h, "GN=([^ ]+)")[, 2]
  bovine   <- toupper(na.omit(gn[grepl("Bos taurus", os)]))
  secreted <- toupper(readLines(hpa_secreted_file, warn = FALSE))
  sort(intersect(unique(bovine), secreted))
}

# Per-protein contaminant flag + reason. genes/descriptions are character vectors.
flag_contaminants <- function(genes, descriptions, serum_genes) {
  is_keratin <- grepl(KERATIN_GENE_PATTERN, genes, ignore.case = TRUE) |
                grepl(KERATIN_DESC_PATTERN, descriptions, ignore.case = TRUE)
  is_serum   <- toupper(genes) %in% serum_genes & !is_keratin
  reason <- ifelse(is_keratin, "keratin",
            ifelse(is_serum, "FBS serum (bovine plasma)", NA_character_))
  data.frame(reason = reason, is_contaminant = is_keratin | is_serum,
             stringsAsFactors = FALSE)
}
