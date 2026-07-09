#!/usr/bin/env Rscript
# Module analysis step 01: module-level inference on the WGCNA modules, on the same
# design as the protein DE (~ 0 + group, paired Replicate block via duplicateCorrelation).
#   - eigengene limma per contrast: the module's shift under Disease/Transplant/Rescue/
#     Interaction, with the same random Replicate block the protein DE uses.
#   - moderated omnibus F per module: any group effect at all (3-df, groups vs Ctl); the
#     principled responsiveness gate, FDR across modules.
#   - fry (self-contained) per contrast: do the module's member proteins respond, and which
#     way, using every protein and the residual covariance rather than the eigengene PC1.
#   - camera (competitive) per contrast: shifted relative to the rest of the proteome.
#     Reported for completeness only; competitive tests are anti-conservative on large
#     co-expression modules (the "background" is itself structured), so fry is the gate.
# At n=6/group these are exploratory; confirmatory signal is the protein DE in a_/b_.
pacman::p_load(here, limma, dplyr, tibble)

GROUP_LEVELS <- c("Ctl", "Mito", "PHE", "PHE_Mito")
CONTRASTS <- c(
  Disease = "PHE - Ctl",
  Transplant = "Mito - Ctl",
  Rescue = "PHE_Mito - PHE",
  Interaction = "(PHE_Mito - PHE) - (Mito - Ctl)"
)

dat <- here("03_DEP", "c_modules", "c_data")
net_cache <- file.path(dat, "wgcna_network.rds")
if (!file.exists(net_cache)) stop("run 00_build_wgcna.R first; missing ", net_cache)
w <- readRDS(net_cache)

dal <- readRDS(here("02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"))
expr <- as.matrix(dal$data)
meta <- as.data.frame(dal$metadata)
grp <- factor(meta$Group[match(colnames(expr), meta$Col_ID)], levels = GROUP_LEVELS)
rep_block <- meta$Replicate[match(colnames(expr), meta$Col_ID)]
design <- stats::model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)
cm <- limma::makeContrasts(contrasts = unname(CONTRASTS), levels = design)
colnames(cm) <- names(CONTRASTS)

non_grey <- setdiff(sub("^ME", "", colnames(w$MEs)), "grey")
me <- t(as.matrix(w$MEs[, paste0("ME", non_grey), drop = FALSE]))
rownames(me) <- non_grey

# eigengene limma with the paired Replicate block, matching the protein DE
eig_rho <- limma::duplicateCorrelation(me, design, block = rep_block)$consensus
eig_fit <- limma::lmFit(me, design, block = rep_block, correlation = eig_rho)
eig_fit <- limma::eBayes(limma::contrasts.fit(eig_fit, cm))
mod_stats <- bind_rows(lapply(names(CONTRASTS), function(cn) {
  tt <- limma::topTable(eig_fit, coef = cn, number = Inf, sort.by = "none")
  tibble(module = non_grey, contrast = cn, logFC = tt$logFC, p = tt$P.Value, fdr = tt$adj.P.Val)
}))

# omnibus moderated F: any group effect on the eigengene (3-df basis = each group vs Ctl)
omni_cm <- limma::makeContrasts(PHE - Ctl, Mito - Ctl, PHE_Mito - Ctl, levels = design)
omni_fit <- limma::eBayes(limma::contrasts.fit(
  limma::lmFit(me, design, block = rep_block, correlation = eig_rho), omni_cm
))
omnibus <- tibble(
  module = non_grey, F = omni_fit$F, p = omni_fit$F.p.value,
  fdr = p.adjust(omni_fit$F.p.value, "BH")
)

# set tests on the member proteins. fry takes the same block; camera does not accept one,
# but at this rho the block is immaterial to it.
mc <- w$module_colors[rownames(expr)]
idx <- split(seq_len(nrow(expr)), mc)[non_grey]
set_rho <- limma::duplicateCorrelation(expr, design, block = rep_block)$consensus
settests <- bind_rows(lapply(names(CONTRASTS), function(cn) {
  fr <- limma::fry(expr, idx, design, contrast = cm[, cn], block = rep_block, correlation = set_rho)
  ca <- limma::camera(expr, idx, design, contrast = cm[, cn])
  tibble(
    module = rownames(fr), contrast = cn, direction = fr$Direction,
    fry_p = fr$PValue, fry_fdr = fr$FDR,
    camera_p = ca[rownames(fr), "PValue"], camera_fdr = ca[rownames(fr), "FDR"]
  )
}))

mod_size <- count(tibble(module = mc), module) |> filter(module %in% non_grey)

res <- list(
  mod_stats = mod_stats, omnibus = omnibus, settests = settests,
  mod_size = mod_size, eig_rho = eig_rho, set_rho = set_rho, contrasts = CONTRASTS
)
saveRDS(res, file.path(dat, "module_stats.rds"))
readr::write_csv(mod_stats, file.path(dat, "module_eigengene_stats.csv"))
readr::write_csv(omnibus, file.path(dat, "module_omnibus_F.csv"))
readr::write_csv(settests, file.path(dat, "module_set_tests.csv"))

cat(sprintf(
  "module stats: %d modules x %d contrasts | eigengene rho %.3f, set rho %.3f\n",
  length(non_grey), length(CONTRASTS), eig_rho, set_rho
))
print(omnibus |> arrange(p), n = Inf)
