#!/usr/bin/env Rscript
# F06 analysis — mitochondrial complex stoichiometry from LFQ abundance data.
# Two validated, parsimonious approaches (no fractionation needed):
#   (1) AlteredPQR  — subunit-pair quantity-ratio remodeling per complex
#       (Buljan et al. 2023, Nat Methods, doi:10.1038/s41592-023-02011-w).
#   (2) Mitonuclear OXPHOS balance — mtDNA- vs nuclear-encoded subunit ratio
#       (Houtkooper et al. 2013, Nature, PMID 23698443).
# Input: imputed (complete) rat LFQ log2 abundance matrix; MitoCarta OXPHOS
# subunit sets (rat). Reads existing pipeline outputs only; additive.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr); library(stringr)
  library(limma); library(AlteredPQR)
})

OUT <- here::here("05_Figures", "F06_complex_mito", "c_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

da  <- readRDS(here::here("02_Imputation", "c_data", "01_DAList_imputed.rds"))
gs  <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))$MitoCarta
comb <- read_csv(here::here("03_DEP", "c_data", "03_combined_results.csv"), show_col_types = FALSE)

# Gene-level log2 abundance matrix (average duplicate gene symbols)
gene <- da$annotation$gene[match(rownames(da$data), da$annotation$uniprot_id)]
keep <- !is.na(gene) & nzchar(gene)
mat  <- limma::avereps(da$data[keep, ], ID = gene[keep])   # genes × 24 samples
grp  <- da$metadata$Group
stopifnot(identical(colnames(mat), da$metadata$Col_ID))

# OXPHOS complex subunit sets (rat) → labelled complexes
COMPLEX_SETS <- c(
  "Complex I"   = "MITOCARTA_OXPHOS__COMPLEX_I__CI_SUBUNITS",
  "Complex III" = "MITOCARTA_OXPHOS__COMPLEX_III__CIII_SUBUNITS",
  "Complex IV"  = "MITOCARTA_OXPHOS__COMPLEX_IV__CIV_SUBUNITS",
  "Complex V"   = "MITOCARTA_OXPHOS__COMPLEX_V__CV_SUBUNITS",
  "Mitoribosome"= "MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA__TRANSLATION__MITOCHONDRIAL_RIBOSOME")
complex_members <- lapply(COMPLEX_SETS, function(s) intersect(gs[[s]], rownames(mat)))
for (cx in names(complex_members))
  message(sprintf("  %s: %d subunits measured", cx, length(complex_members[[cx]])))

# int_pairs: all within-complex subunit pairs (gene symbols)
# AlteredPQR builds/splits pair names on "-", so sanitize hyphenated rat mtDNA
# symbols (Mt-co2 -> Mt.co2) for the AlteredPQR matrix + pairs only; map back later.
san <- function(x) gsub("-", ".", x)
mat_apqr <- mat; rownames(mat_apqr) <- san(rownames(mat))
orig_of  <- setNames(rownames(mat), san(rownames(mat)))

int_pairs <- bind_rows(lapply(names(complex_members), function(cx) {
  g <- san(complex_members[[cx]])
  if (length(g) < 2) return(NULL)
  cb <- t(combn(g, 2))
  tibble(ProtA = cb[, 1], ProtB = cb[, 2], complex = cx)
}))
pair_to_complex <- setNames(int_pairs$complex,
                            paste(int_pairs$ProtA, int_pairs$ProtB, sep = "-"))
message(sprintf("  %d within-complex subunit pairs across %d complexes",
                nrow(int_pairs), length(complex_members)))

# AlteredPQR per contrast (ref vs test groups; defaults per Buljan 2023)
# AlteredPQR_RB() reads globals: quant_data_all, cols_with_reference_data, int_pairs.
CONTRASTS <- list(
  list(name = "Transplant_Mito",  ref = "Ctl", test = "Mito"),
  list(name = "Disease_Phe",      ref = "Ctl", test = "PHE"),
  list(name = "Rescue_Mito+Phe",  ref = "PHE", test = "PHE_Mito"))

# AlteredPQR reads bare globals; assign directly into the global env (the
# package ships locked datasets of the same names, so <<- would error).
assign("int_pairs", as.data.frame(int_pairs[, c("ProtA", "ProtB")]), envir = globalenv())

run_apqr <- function(cc) {
  ref_cols  <- which(grp == cc$ref)
  test_cols <- which(grp == cc$test)
  qd <- mat_apqr[, c(ref_cols, test_cols)]
  assign("quant_data_all", as.data.frame(qd), envir = globalenv())
  assign("cols_with_reference_data", seq_along(ref_cols), envir = globalenv())
  rp <- AlteredPQR_RB(modif_z_score_threshold = 3.5,
                      fraction_of_samples_threshold = 0.1)
  if (is.null(rp) || nrow(rp) == 0) return(NULL)
  rev_key <- function(pp) sapply(strsplit(pp, "-"), function(x) paste(rev(x), collapse = "-"))
  rp |>
    mutate(contrast = cc$name,
           complex  = pair_to_complex[as.character(Protein_pair)],
           complex  = ifelse(is.na(complex),
                             pair_to_complex[rev_key(as.character(Protein_pair))], complex),
           pair_readable = sapply(strsplit(as.character(Protein_pair), "-"),
                                  function(p) paste(orig_of[p], collapse = " / ")))
}
apqr <- bind_rows(lapply(CONTRASTS, run_apqr))
if (nrow(apqr) == 0) message("  AlteredPQR: no altered pairs at z>3.5 / 10% samples")
write_csv(apqr, file.path(OUT, "alteredpqr_results.csv"))
apqr_counts <- apqr |>
  count(contrast, complex, Change, name = "n_pairs")
write_csv(apqr_counts, file.path(OUT, "alteredpqr_counts.csv"))
message(sprintf("  AlteredPQR: %d altered subunit pairs total", nrow(apqr)))

# Mitonuclear OXPHOS balance (Houtkooper 2013)
oxphos_all <- intersect(gs[["MITOCARTA_OXPHOS__OXPHOS_SUBUNITS"]], rownames(mat))
mtdna_canonical <- c("Mt-nd1","Mt-nd2","Mt-nd3","Mt-nd4","Mt-nd4l","Mt-nd5","Mt-nd6",
                     "Mt-co1","Mt-co2","Mt-co3","Mt-cyb","Mt-atp6","Mt-atp8")
mtdna_genes <- intersect(mtdna_canonical, rownames(mat))
nuclear_oxphos <- setdiff(oxphos_all, mtdna_genes)
message(sprintf("  Mitonuclear: %d mtDNA-encoded + %d nuclear-encoded OXPHOS subunits",
                length(mtdna_genes), length(nuclear_oxphos)))

balance <- tibble(
  Col_ID  = colnames(mat),
  Group   = grp,
  mtDNA_mean    = colMeans(mat[mtdna_genes, , drop = FALSE]),
  nuclear_mean  = colMeans(mat[nuclear_oxphos, , drop = FALSE]),
  mitonuclear_log2ratio = mtDNA_mean - nuclear_mean)   # log2(mtDNA / nuclear)
write_csv(balance, file.path(OUT, "mitonuclear_balance.csv"))
# Inferential stats for this metric live in 02_lmm_dynamics.R (design-aware LMM
# + emmeans contrasts, BH-FDR) — the single canonical framework reported.

# Per-complex subunit logFC (coordinated-stoichiometry view)
contrast_old <- c(Transplant_Mito = "CTLvMITO", Disease_Phe = "CTLvPHE",
                  `Rescue_Mito+Phe` = "PHEvPHE_MITO")
subunit_lfc <- bind_rows(lapply(names(complex_members), function(cx) {
  g <- complex_members[[cx]]
  bind_rows(lapply(names(contrast_old), function(cn) {
    lfc_col <- paste0("logFC_", contrast_old[cn]); pi_col <- paste0("pi_score_", contrast_old[cn])
    sub <- comb |> filter(gene %in% g) |>
      transmute(gene, complex = cx, contrast = cn,
                logFC = .data[[lfc_col]], pi_score = .data[[pi_col]])
    sub
  }))
}))
write_csv(subunit_lfc, file.path(OUT, "complex_subunit_logfc.csv"))

saveRDS(list(complex_members = complex_members, mtdna_genes = mtdna_genes,
             nuclear_oxphos = nuclear_oxphos),
        file.path(OUT, "analysis_meta.rds"))
message("F06 analysis done")
