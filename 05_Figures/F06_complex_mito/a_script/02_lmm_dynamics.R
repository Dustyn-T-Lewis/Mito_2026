#!/usr/bin/env Rscript
# F06 (2/3) — design-aware stats + mitochondrial dynamics + content.
# Refactors the F06 stats to respect the 2×2 factorial (PHE × Mito) and the
# paired-passage block (Replicate 1-6 across all four groups):
#   • Mitonuclear balance: lmer(y ~ PHE*Mito + (1|Replicate)) + emmeans contrasts
#     (Transplant/Disease/Rescue/Interaction), BH-FDR. Reports Replicate variance
#     / singularity so the block assumption is data-checked, not assumed.
#   • Per-pair complex stoichiometry: blocked lm(log-ratio ~ PHE*Mito + Replicate)
#     per within-complex subunit pair; 4 effects; BH-FDR. Model-based, design-aware
#     alternative to AlteredPQR (Buljan 2023) — appropriate at n=6/group.
#   • Mito content (mass) per sample + mitochondrial-dynamics axis module scores.
# Additive; reads existing pipeline outputs only.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr); library(stringr)
  library(limma); library(lme4); library(lmerTest); library(emmeans)
})
emm_options(lmerTest.limit = 50, pbkrtest.limit = 50)

OUT <- here::here("05_Figures", "F06_complex_mito", "c_data")
da  <- readRDS(here::here("02_Imputation", "c_data", "01_DAList_imputed.rds"))
gs  <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))$MitoCarta
meta_an <- readRDS(file.path(OUT, "analysis_meta.rds"))

# Gene-level abundance + factorial/block design
gene <- da$annotation$gene[match(rownames(da$data), da$annotation$uniprot_id)]
keep <- !is.na(gene) & nzchar(gene)
mat  <- limma::avereps(da$data[keep, ], ID = gene[keep])
md <- da$metadata
design <- tibble(
  Col_ID    = md$Col_ID,
  PHE       = as.integer(md$Group %in% c("PHE", "PHE_Mito")),
  Mito      = as.integer(md$Group %in% c("Mito", "PHE_Mito")),
  Replicate = factor(md$Replicate),
  Group     = factor(md$Group, levels = c("Ctl", "Mito", "PHE", "PHE_Mito")))
stopifnot(identical(design$Col_ID, colnames(mat)))

# emmeans cell order for ~ Mito + PHE (Mito fastest): M0P0, M1P0, M0P1, M1P1
CONTRAST_DEFS <- list(
  Transplant_Mito = c(-1, 1, 0, 0),    # Mito − Ctl   (Mito effect at PHE 0)
  Disease_Phe     = c(-1, 0, 1, 0),    # PHE − Ctl    (PHE effect at Mito 0)
  `Rescue_Mito+Phe` = c(0, 0, -1, 1),  # PHE_Mito − PHE (Mito effect at PHE 1)
  Interaction     = c(1, -1, -1, 1))   # (PHE_Mito−PHE) − (Mito−Ctl)

# LMM with Replicate block; emmeans factorial contrasts; report singularity
lmm_contrasts <- function(y) {
  d <- design |> mutate(y = y, Mito = factor(Mito), PHE = factor(PHE))
  m <- suppressMessages(suppressWarnings(
    lmerTest::lmer(y ~ PHE * Mito + (1 | Replicate), data = d,
                   control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4)))))
  vc <- as.data.frame(lme4::VarCorr(m))
  rep_var <- vc$vcov[vc$grp == "Replicate"]; resid_var <- vc$vcov[vc$grp == "Residual"]
  icc <- rep_var / (rep_var + resid_var)
  emm <- emmeans(m, ~ Mito + PHE)
  ct  <- summary(contrast(emm, method = CONTRAST_DEFS), infer = c(TRUE, TRUE))
  tibble(contrast = ct$contrast, estimate = ct$estimate, se = ct$SE,
         p_value = ct$p.value, replicate_icc = icc, singular = lme4::isSingular(m))
}

# (a) Mitonuclear balance
oxphos_all <- intersect(gs[["MITOCARTA_OXPHOS__OXPHOS_SUBUNITS"]], rownames(mat))
mtdna <- intersect(c("Mt-nd1","Mt-nd2","Mt-nd3","Mt-nd4","Mt-nd4l","Mt-nd5","Mt-nd6",
                     "Mt-co1","Mt-co2","Mt-co3","Mt-cyb","Mt-atp6","Mt-atp8"), rownames(mat))
nuc_ox <- setdiff(oxphos_all, mtdna)
mn_ratio <- colMeans(mat[mtdna, , drop = FALSE]) - colMeans(mat[nuc_ox, , drop = FALSE])
mn_lmm <- lmm_contrasts(mn_ratio) |> mutate(metric = "mitonuclear_balance", .before = 1)
mn_lmm$p_fdr <- p.adjust(mn_lmm$p_value, "BH")
write_csv(mn_lmm, file.path(OUT, "mitonuclear_lmm.csv"))
message(sprintf("  Mitonuclear LMM: Replicate ICC = %.3f%s",
                mn_lmm$replicate_icc[1], if (mn_lmm$singular[1]) " (singular → block uninformative)" else ""))

# (b) Mito content (mass) — mean of all MitoCarta proteins
mito_all <- intersect(unique(unlist(gs)), rownames(mat))
mito_content <- colMeans(mat[mito_all, , drop = FALSE])
content_df <- design |> mutate(mito_content = mito_content)
write_csv(content_df, file.path(OUT, "mito_content.csv"))
mc_lmm <- lmm_contrasts(mito_content) |> mutate(metric = "mito_content", .before = 1)
mc_lmm$p_fdr <- p.adjust(mc_lmm$p_value, "BH")
write_csv(mc_lmm, file.path(OUT, "mito_content_lmm.csv"))

# (c) Mitochondrial-dynamics axes — module z-scores + LMM contrasts
DYN_SETS <- c(
  Fission     = "MITOCARTA_MITOCHONDRIAL_DYNAMICS_AND_SURVEILLANCE__FISSION",
  Fusion      = "MITOCARTA_MITOCHONDRIAL_DYNAMICS_AND_SURVEILLANCE__FUSION",
  Mitophagy   = "MITOCARTA_MITOCHONDRIAL_DYNAMICS_AND_SURVEILLANCE__MITOPHAGY",
  `Cristae/MICOS` = "MITOCARTA_MITOCHONDRIAL_DYNAMICS_AND_SURVEILLANCE__CRISTAE_FORMATION__MICOS_COMPLEX",
  `Proteostasis (UPRmt)` = "MITOCARTA_PROTEIN_IMPORT_SORTING_AND_HOMEOSTASIS__PROTEIN_HOMEOSTASIS")
z <- t(scale(t(mat)))   # gene-wise z across samples
dyn_scores <- lapply(DYN_SETS, function(s) {
  g <- intersect(gs[[s]], rownames(z)); if (length(g) < 3) return(NULL)
  colMeans(z[g, , drop = FALSE])
})
dyn_scores <- dyn_scores[!sapply(dyn_scores, is.null)]
dyn_long <- bind_rows(lapply(names(dyn_scores), function(ax)
  design |> mutate(axis = ax, score = dyn_scores[[ax]]))) |>
  select(Col_ID, Group, axis, score)
write_csv(dyn_long, file.path(OUT, "dynamics_scores.csv"))
dyn_lmm <- bind_rows(lapply(names(dyn_scores), function(ax)
  lmm_contrasts(dyn_scores[[ax]]) |> mutate(axis = ax, .before = 1)))
dyn_lmm <- dyn_lmm |> group_by(contrast) |> mutate(p_fdr = p.adjust(p_value, "BH")) |> ungroup()
write_csv(dyn_lmm, file.path(OUT, "dynamics_lmm.csv"))
for (ax in names(dyn_scores))
  message(sprintf("  Dynamics axis %-20s: %d members", ax, length(intersect(gs[[DYN_SETS[ax]]], rownames(z)))))

# (d) Per-pair complex stoichiometry — blocked lm, design-aware, BH-FDR
# Blocked RCBD lm (Replicate fixed block) ≈ random-block inference for a balanced
# design, but fast/robust across thousands of pairs. log-ratio = logA − logB.
complex_members <- meta_an$complex_members
pairs <- bind_rows(lapply(names(complex_members), function(cx) {
  g <- complex_members[[cx]]; if (length(g) < 2) return(NULL)
  cb <- t(combn(g, 2)); tibble(A = cb[, 1], B = cb[, 2], complex = cx)
}))
PHEn <- design$PHE; Miton <- design$Mito; Repl <- design$Replicate
pair_effect <- function(i) {
  lr <- mat[pairs$A[i], ] - mat[pairs$B[i], ]
  fit <- lm(lr ~ PHEn * Miton + Repl)
  b <- coef(fit); V <- vcov(fit); rdf <- df.residual(fit)
  est <- c(Transplant_Mito = unname(b["Miton"]),
           Disease_Phe     = unname(b["PHEn"]),
           Interaction     = unname(b["PHEn:Miton"]),
           `Rescue_Mito+Phe`= unname(b["Miton"] + b["PHEn:Miton"]))
  se  <- c(sqrt(V["Miton","Miton"]), sqrt(V["PHEn","PHEn"]),
           sqrt(V["PHEn:Miton","PHEn:Miton"]),
           sqrt(V["Miton","Miton"] + V["PHEn:Miton","PHEn:Miton"] + 2*V["Miton","PHEn:Miton"]))
  p <- 2 * pt(abs(est / se), rdf, lower.tail = FALSE)
  tibble(complex = pairs$complex[i], pair = paste(pairs$A[i], pairs$B[i], sep = "/"),
         contrast = names(est), estimate = est, p_value = p)
}
message(sprintf("  Fitting blocked lm for %d subunit pairs ...", nrow(pairs)))
pair_res <- bind_rows(lapply(seq_len(nrow(pairs)), pair_effect))
pair_res <- pair_res |> group_by(contrast) |> mutate(p_fdr = p.adjust(p_value, "BH")) |> ungroup()
write_csv(pair_res, file.path(OUT, "pair_stoich_lmm.csv"))
pair_counts <- pair_res |> filter(p_fdr < 0.05) |>
  count(complex, contrast, name = "n_altered")
total_pairs <- sapply(complex_members, function(g) choose(length(g), 2))
pair_counts <- pair_counts |> mutate(total = total_pairs[complex], pct = 100 * n_altered / total)
write_csv(pair_counts, file.path(OUT, "pair_stoich_counts.csv"))
message(sprintf("  Stoichiometry LMM: %d pair-effects FDR<0.05",
                sum(pair_res$p_fdr < 0.05, na.rm = TRUE)))

# (e) Coordinated set shift — limma::camera (correlation-aware self-contained test)
# Replaces a per-subunit t.test(logFC vs 0), which treats co-regulated subunits as
# independent. camera accounts for inter-subunit correlation. Design matches the DEP
# (~ 0 + group, no block) so set p-values are consistent with 03_combined_results.csv.
g_fac <- design$Group
des   <- model.matrix(~ 0 + g_fac); colnames(des) <- levels(g_fac)
cont  <- limma::makeContrasts(
  Transplant_Mito   = Mito - Ctl,
  Disease_Phe       = PHE - Ctl,
  `Rescue_Mito+Phe` = PHE_Mito - PHE,
  levels = des)
camera_sets <- c(
  lapply(complex_members, function(g) which(rownames(mat) %in% g)),
  list(`OXPHOS_all` = which(rownames(mat) %in% unlist(complex_members[OXPHOS_C <- c(
    "Complex I", "Complex III", "Complex IV", "Complex V")]))),
  lapply(DYN_SETS, function(s) which(rownames(mat) %in% gs[[s]])))
camera_sets <- camera_sets[lengths(camera_sets) >= 3]
camera_res <- bind_rows(lapply(colnames(cont), function(cn) {
  cr <- limma::camera(mat, camera_sets, design = des, contrast = cont[, cn])
  tibble(set = rownames(cr), contrast = cn, n_genes = cr$NGenes,
         direction = cr$Direction, p_value = cr$PValue, fdr = cr$FDR)
}))
oxphos_lab <- c("Complex I", "Complex III", "Complex IV", "Complex V", "Mitoribosome", "OXPHOS_all")
camera_res |> filter(set %in% oxphos_lab) |>
  write_csv(file.path(OUT, "complex_camera.csv"))
camera_res |> filter(set %in% names(DYN_SETS)) |>
  write_csv(file.path(OUT, "dynamics_camera.csv"))
message(sprintf("  camera: OXPHOS_all p = %s across contrasts",
                paste(sprintf("%.3g", camera_res$p_value[camera_res$set == "OXPHOS_all"]), collapse = " / ")))

# (f) MitoCarta-pathway camera — broad hierarchy for the F06 dot plot
# Same correlation-aware competitive test, over the MitoCarta hierarchy (drop the
# all-mito set + bare compartments; keep most-specific terms with >=3 measured).
mito_drop <- c("MITOCARTA_ALL", "MITOCARTA_IMM", "MITOCARTA_IMS", "MITOCARTA_MATRIX", "MITOCARTA_OMM")
mito_path_sets <- gs[setdiff(names(gs), mito_drop)] |>
  lapply(function(g) which(rownames(mat) %in% g))
mito_path_sets <- mito_path_sets[lengths(mito_path_sets) >= 3]
mito_path_camera <- bind_rows(lapply(colnames(cont), function(cn) {
  cr <- limma::camera(mat, mito_path_sets, design = des, contrast = cont[, cn])
  tibble(set = rownames(cr), contrast = cn, n_genes = cr$NGenes,
         direction = cr$Direction, p_value = cr$PValue, fdr = cr$FDR)
}))
write_csv(mito_path_camera, file.path(OUT, "mito_pathway_camera.csv"))
message(sprintf("  MitoCarta-pathway camera: %d sets x %d contrasts; %d at FDR<0.05",
                length(mito_path_sets), ncol(cont), sum(mito_path_camera$fdr < 0.05)))

message("F06 LMM + dynamics analysis done")
