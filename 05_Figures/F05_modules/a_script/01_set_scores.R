#!/usr/bin/env Rscript
# F05 (1/2) — knowledge-driven per-sample set scores tested in the 2x2 LMM.
# Replaces WGCNA's data-derived modules as the headline grouping: sets are defined
# EXTERNALLY (CORUM complexes + MitoCarta hierarchy + Reactome), so there is no
# double-dipping (modules are not learned from the same 24 samples they are tested on).
# GSVA (Hanzelmann 2013) is the primary score; singscore (Foroutan 2018) is a
# rank-based cross-check. Each set score is tested with
# lmer(score ~ PHE*Mito + (1|Replicate)) + emmeans + BH-FDR.
# Reads frozen pipeline outputs only; additive; seeded.
#
# pkgs: GSVA 2.6.2, singscore 1.32.0

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr); library(stringr)
  library(limma); library(lme4); library(lmerTest); library(emmeans)
  library(GSVA); library(singscore)
})
emm_options(lmerTest.limit = 50, pbkrtest.limit = 50)
set.seed(42)

OUT <- here::here("05_Figures", "F05_modules", "c_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

da <- readRDS(here::here("02_Imputation", "c_data", "01_DAList_imputed.rds"))
gs <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))
corum <- readRDS(here::here("04_Figures", "shared", "corum_rat_gene_sets.rds"))

# Gene-level log2 matrix (average duplicate symbols), mirrors F06
gene <- da$annotation$gene[match(rownames(da$data), da$annotation$uniprot_id)]
keep <- !is.na(gene) & nzchar(gene)
mat  <- limma::avereps(da$data[keep, ], ID = gene[keep])
md   <- da$metadata
stopifnot(identical(colnames(mat), md$Col_ID))

# Harmonized external set collection
# Source-tagged; MitoCarta drops the all-mito set and bare compartments (IMM/IMS/
# MATRIX/OMM) per stage convention; display = most-specific term.
MIN_MEASURED <- 5L; MAX_MEASURED <- 300L
mito_drop <- c("MITOCARTA_ALL", "MITOCARTA_IMM", "MITOCARTA_IMS",
               "MITOCARTA_MATRIX", "MITOCARTA_OMM")
mito_sets <- gs$MitoCarta[setdiff(names(gs$MitoCarta), mito_drop)]

clean_set_label <- function(x) {
  x <- str_replace_all(x, "_", " ") |> str_to_title()
  x <- str_replace_all(x, c("Oxphos" = "OXPHOS", "Ci " = "CI ", "Ciii " = "CIII ",
                            "Civ " = "CIV ", "Cv " = "CV ", "Dna" = "DNA",
                            "Rna" = "RNA", "Tca" = "TCA", "Micos" = "MICOS"))
  str_squish(x)
}
last_seg <- function(x) {
  seg <- str_split(str_remove(x, "^MITOCARTA_"), "__")
  vapply(seg, function(s) clean_set_label(tail(s, 1)), character(1))
}
clean_reactome <- function(x) {
  x <- str_remove(x, "^REACTOME_") |> str_replace_all("_", " ") |> str_to_title()
  str_squish(x)
}

raw_sets <- c(
  setNames(mito_sets, paste0("MITO::", names(mito_sets))),
  setNames(gs$Reactome, paste0("REAC::", names(gs$Reactome))),
  setNames(corum, paste0("CORUM::", names(corum))))

source_of  <- function(id) sub("::.*$", "", id)
display_of <- function(id) {
  body <- sub("^[^:]+::", "", id)
  src  <- source_of(id)
  if (src == "MITO")  return(last_seg(body))
  if (src == "REAC")  return(clean_reactome(body))
  str_squish(body)                       # CORUM names are already readable
}

# Restrict every set to measured genes; keep sets in the size window
genes_in <- rownames(mat)
sets_measured <- lapply(raw_sets, function(g) intersect(unique(g), genes_in))
nmemb <- lengths(sets_measured)
keep_set <- nmemb >= MIN_MEASURED & nmemb <= MAX_MEASURED
sets_use <- sets_measured[keep_set]
# Drop duplicate gene-content sets (keep first) to reduce redundant testing.
sig <- vapply(sets_use, function(g) paste(sort(g), collapse = ""), character(1))
sets_use <- sets_use[!duplicated(sig)]

set_meta <- tibble(set_id = names(sets_use),
                   source = source_of(names(sets_use)),
                   display = vapply(names(sets_use), display_of, character(1)),
                   n_members = lengths(sets_use))
message(sprintf("  harmonized sets: %d (CORUM %d / MitoCarta %d / Reactome %d) | %d genes",
                nrow(set_meta), sum(set_meta$source == "CORUM"),
                sum(set_meta$source == "MITO"), sum(set_meta$source == "REAC"),
                length(genes_in)))

# Per-sample scores: GSVA (primary) + singscore (robustness)
gp <- gsvaParam(mat, sets_use, minSize = MIN_MEASURED, maxSize = MAX_MEASURED,
                kcdf = "Gaussian")
gsva_mat <- gsva(gp, verbose = FALSE)
gsva_mat <- gsva_mat[intersect(names(sets_use), rownames(gsva_mat)), , drop = FALSE]

# singscore: each set treated as one directional signature (knownDirection = TRUE,
# default) so its TotalScore is a per-sample set-abundance score directly comparable
# to GSVA. knownDirection = FALSE would instead score bidirectional rank dispersion,
# a different quantity that does not track GSVA's directional contrasts.
rank_data <- rankGenes(mat)
ss_mat <- t(vapply(rownames(gsva_mat), function(id)
  simpleScore(rank_data, upSet = sets_use[[id]])$TotalScore,
  numeric(ncol(mat))))
colnames(ss_mat) <- colnames(mat)

write_csv(as_tibble(gsva_mat, rownames = "set_id"), file.path(OUT, "set_scores_gsva.csv"))
write_csv(as_tibble(ss_mat,   rownames = "set_id"), file.path(OUT, "set_scores_singscore.csv"))

# 2x2 factorial LMM with Replicate block (mirrors F06 02_lmm_dynamics.R)
design <- tibble(
  Col_ID    = md$Col_ID,
  PHE       = factor(as.integer(md$Group %in% c("PHE", "PHE_Mito"))),
  Mito      = factor(as.integer(md$Group %in% c("Mito", "PHE_Mito"))),
  Replicate = factor(md$Replicate))
# emmeans cell order ~ Mito + PHE: M0P0, M1P0, M0P1, M1P1
CONTRAST_DEFS <- list(
  Transplant_Mito   = c(-1, 1, 0, 0),
  Disease_Phe       = c(-1, 0, 1, 0),
  `Rescue_Mito+Phe` = c(0, 0, -1, 1),
  Interaction       = c(1, -1, -1, 1))

lmm_contrasts <- function(y) {
  d <- design |> mutate(y = y)
  m <- suppressMessages(suppressWarnings(
    lmerTest::lmer(y ~ PHE * Mito + (1 | Replicate), data = d,
                   control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4)))))
  vc <- as.data.frame(lme4::VarCorr(m))
  icc <- vc$vcov[vc$grp == "Replicate"] / (vc$vcov[vc$grp == "Replicate"] + vc$vcov[vc$grp == "Residual"])
  emm <- emmeans(m, ~ Mito + PHE)
  ct  <- summary(contrast(emm, method = CONTRAST_DEFS), infer = c(TRUE, TRUE))
  tibble(contrast = ct$contrast, estimate = ct$estimate, se = ct$SE,
         p_value = ct$p.value, replicate_icc = icc, singular = lme4::isSingular(m))
}

run_lmm_over <- function(score_mat) {
  bind_rows(lapply(rownames(score_mat), function(id)
    lmm_contrasts(score_mat[id, ]) |> mutate(set_id = id, .before = 1))) |>
    left_join(set_meta, by = "set_id") |>
    group_by(contrast) |> mutate(p_fdr = p.adjust(p_value, "BH")) |> ungroup()
}

gsva_lmm <- run_lmm_over(gsva_mat)
ss_lmm   <- run_lmm_over(ss_mat)
write_csv(gsva_lmm, file.path(OUT, "set_lmm_gsva.csv"))
write_csv(ss_lmm,   file.path(OUT, "set_lmm_singscore.csv"))

# Method-robustness: per-contrast Spearman of GSVA vs singscore effect sizes
concord <- inner_join(
  gsva_lmm |> select(set_id, contrast, est_gsva = estimate),
  ss_lmm   |> select(set_id, contrast, est_ss = estimate),
  by = c("set_id", "contrast")) |>
  group_by(contrast) |>
  summarise(n = n(), spearman = cor(est_gsva, est_ss, method = "spearman"),
            .groups = "drop")
write_csv(concord, file.path(OUT, "set_method_concordance.csv"))

saveRDS(list(set_meta = set_meta, sets_use = sets_use,
             min_measured = MIN_MEASURED, max_measured = MAX_MEASURED),
        file.path(OUT, "harmonized_sets.rds"))

message(sprintf("  GSVA-vs-singscore effect-size Spearman per contrast: %s",
                paste(sprintf("%s %.2f", concord$contrast, concord$spearman), collapse = " | ")))
message(sprintf("  FDR<0.05 sets per contrast (GSVA): %s",
                paste(sprintf("%s %d", unique(gsva_lmm$contrast),
                              tapply(gsva_lmm$p_fdr < 0.05, gsva_lmm$contrast, sum)[unique(gsva_lmm$contrast)]),
                      collapse = " | ")))
message("F05 set-score analysis done")
