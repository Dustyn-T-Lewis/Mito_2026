#!/usr/bin/env Rscript
# Build the signed WGCNA co-abundance network the module statistics and the F03 render
# consume, into this figure's c_data (regenerated every run; the .rds is just a cache).
# Best-practice settings (WGCNA package docs + Horvath-lab FAQ), tuned for proteomics:
#   - signed network + signed TOM: direction-coherent modules (Horvath recommendation).
#   - biweight midcorrelation (bicor, maxPOutliers = 0.05): robust to outlier samples,
#     which matters at n=6/group and on an imputed (correlation-inflating) matrix.
#   - soft power by the scale-free criterion (lowest power with R^2 >= 0.85, lands on 8);
#     if none is reached, the signed-network sample-size fallback (n in [20,30) -> 16).
#   - minModuleSize 20, mergeCutHeight 0.25, deepSplit 2: 15 modules. Power 8 vs 16 and
#     merge 0.25 vs 0.15 leave the solution unchanged or grey-heavier (sensitivity check).
#   - single block (maxBlockSize > n proteins): global, deterministic module detection.
# Preservation (Zsummary) of the full-cohort modules is tested in each condition against
# the control reference (Langfelder 2011): the n=6-robust rescue statistic, also cached.
pacman::p_load(WGCNA, dplyr, tibble, here)
source(here::here("04_Figures", "functions", "shared_data_loaders.R"))
WGCNA::disableWGCNAThreads() # single-threaded for a reproducible result
set.seed(42)
RSQ_CUT <- 0.85
PRESERVATION_PERMS <- 200L
out_dir <- here::here("04_Figures", "F03_WGCNA", "c_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(P05$imp_rds)
expr <- as.matrix(dal$data) # proteins x samples, log2, no gaps
datExpr <- t(expr) # WGCNA wants samples x genes
stopifnot(WGCNA::goodSamplesGenes(datExpr, verbose = 0)$allOK)

# soft power: scale-free fit with bicor; fall back to the signed-network table if unmet
sft <- WGCNA::pickSoftThreshold(
  datExpr,
  networkType = "signed", corFnc = WGCNA::bicor, corOptions = list(maxPOutliers = 0.05),
  RsquaredCut = RSQ_CUT, powerVector = c(1:10, seq(12, 20, 2)), verbose = 0
)
# fall back to the WGCNA signed-network value for n in [20,30) if no power reaches the cut
power <- if (!is.na(sft$powerEstimate)) sft$powerEstimate else 16L
power_basis <- if (!is.na(sft$powerEstimate)) sprintf("scale-free R2>=%.2f", RSQ_CUT) else "signed sample-size fallback"

net <- WGCNA::blockwiseModules(
  datExpr,
  power = power, networkType = "signed", TOMType = "signed",
  corType = "bicor", maxPOutliers = 0.05,
  minModuleSize = 20, mergeCutHeight = 0.25, deepSplit = 2,
  numericLabels = FALSE, maxBlockSize = 5000,
  saveTOMs = FALSE, verbose = 0
)

module_colors <- net$colors # one colour per protein, in datExpr column order
mes <- net$MEs # samples x modules, colnames ME<colour>
kme <- WGCNA::bicor(datExpr, mes, maxPOutliers = 0.05) # proteins x modules
colnames(kme) <- sub("^ME", "kME_", colnames(mes)) # match the kME_<colour> contract

ann <- dal$annotation[rownames(expr), c("uniprot_id", "gene", "description"), drop = FALSE]
non_grey <- setdiff(unique(module_colors), "grey")

hubs <- tibble(uniprot_id = rownames(expr), module = module_colors, gene = ann$gene, description = ann$description) |>
  filter(module %in% non_grey) |>
  mutate(kME = vapply(seq_len(n()), function(i) kme[uniprot_id[i], paste0("kME_", module[i])], numeric(1))) |>
  group_by(module) |>
  slice_max(kME, n = 5, with_ties = FALSE) |>
  ungroup() |>
  select(module, uniprot_id, kME, gene, description)

# Module preservation: do the full-cohort modules, evaluated in the control samples,
# reproduce within each treatment's samples? Low under PHE and recovered under
# PHE+MitoTx is the rescue signal. The permutation null is built within the reference
# network, so this is the most defensible cross-condition statistic at n=6.
meta <- as.data.frame(dal$metadata)
groups <- c("Ctl", "PHE", "Mito", "PHE_Mito")
multi_expr <- lapply(groups, function(g) {
  list(data = t(expr[, meta$Col_ID[meta$Group == g], drop = FALSE]))
})
names(multi_expr) <- groups
mp <- WGCNA::modulePreservation(
  multi_expr, list(Ctl = module_colors),
  referenceNetworks = 1, networkType = "signed",
  nPermutations = PRESERVATION_PERMS, randomSeed = 42, maxModuleSize = 5000,
  savePermutedStatistics = FALSE, verbose = 0
)
ref_z <- mp$preservation$Z[["ref.Ctl"]]
ref_rank <- mp$preservation$observed[["ref.Ctl"]]
test_slots <- grep("^inColumnsAlsoPresentIn\\.", names(ref_z), value = TRUE)
preservation <- bind_rows(lapply(test_slots, function(s) {
  z <- ref_z[[s]]
  if (!is.data.frame(z)) {
    return(NULL) # reference-vs-itself slot is NA
  }
  condition <- sub("^inColumnsAlsoPresentIn\\.", "", s)
  tibble(
    module = rownames(z), mod_size = z$moduleSize,
    Zsummary = z$Zsummary.pres, median_rank = ref_rank[[s]]$medianRank.pres,
    condition = condition
  )
})) |>
  filter(!module %in% c("gold", "grey"), condition != "Ctl") |>
  mutate(condition = factor(condition, levels = c("PHE", "Mito", "PHE_Mito")))

w <- list(
  net = net,
  module_colors = module_colors,
  MEs = mes,
  kME = kme,
  hubs = hubs,
  preservation = preservation,
  chosen_power = power,
  sft_df = sft$fitIndices,
  sample_order = colnames(expr),
  gene_order = rownames(expr),
  ann = ann
)
saveRDS(w, file.path(out_dir, "wgcna_network.rds"))

cat(sprintf(
  "WGCNA rebuilt (bicor, signed): %d proteins x %d samples | power %d (%s) | %d modules (+grey)\n",
  nrow(expr), ncol(expr), power, power_basis, length(non_grey)
))
print(sort(table(module_colors), decreasing = TRUE))
