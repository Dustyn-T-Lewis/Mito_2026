#!/usr/bin/env Rscript
# WGCNA network build (signed Pearson; n = 24 is exploratory) for F05_modules.
# Run before 05_Figures/F05_modules/a_script/91_wgcna_supp.R. ~30 s.
# WGCNA is the EXPLORATORY de-novo cross-check (supplement); the F05 headline is
# the knowledge-driven set-score modules (01_set_scores.R + 02_main_panels.R).

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(WGCNA); library(readr)
})

source(here::here("00_input", "h9c2_design.R"))

BASE <- here::here("05_Figures", "F05_modules")
DAT  <- file.path(BASE, "c_data")
RPT  <- file.path(BASE, "b_reports")
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(RPT, "wgcna_run.log")
log_msg <- function(...) {
  msg <- sprintf(...)
  message(msg); cat(msg, "\n", file = LOG, append = TRUE)
}

cat("WGCNA build log\n",
    sprintf("Started: %s\n\n", Sys.time()), file = LOG)

set.seed(42)
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 4)

# 1. Load imputed matrix + metadata

dal <- readRDS(here::here("02_Imputation", "c_data", "01_DAList_imputed.rds"))
mat <- as.matrix(dal$data)                # rows = proteins, cols = samples
meta <- as.data.frame(dal$metadata)

# WGCNA expects datExpr: rows = samples, cols = genes/proteins
datExpr <- t(mat)
log_msg("Input: %d samples x %d proteins", nrow(datExpr), ncol(datExpr))

stopifnot(setequal(rownames(datExpr), meta$Col_ID))
datExpr <- datExpr[meta$Col_ID, , drop = FALSE]

# Drop near-constant proteins (WGCNA needs variance per protein)
gsg <- goodSamplesGenes(datExpr, verbose = 0)
log_msg("goodSamplesGenes: %d/%d proteins kept, %d/%d samples kept",
        sum(gsg$goodGenes), length(gsg$goodGenes),
        sum(gsg$goodSamples), length(gsg$goodSamples))
if (!all(gsg$allOK))
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# 2. Soft-threshold scan (signed Pearson)

powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = "signed", corFnc = "cor",
                         verbose = 0)
sft_df <- sft$fitIndices |>
  rename(power = Power, R2 = SFT.R.sq, slope = slope,
         mean_k = mean.k., median_k = median.k., max_k = max.k.)
write_csv(sft_df, file.path(DAT, "wgcna_softthresh.csv"))

# Smallest power achieving R^2 >= 0.8 with negative slope; cap at 14 for n < 30
# (signed-network rule of thumb for small samples).
ok_pwrs <- with(sft_df, power[R2 >= 0.80 & slope < 0])
chosen_power <- if (length(ok_pwrs) > 0) min(ok_pwrs) else 14L
if (chosen_power > 14L) chosen_power <- 14L
log_msg("Soft-threshold scan: chose power = %d (R^2 = %.2f, mean k = %.1f)",
        chosen_power,
        sft_df$R2[sft_df$power == chosen_power],
        sft_df$mean_k[sft_df$power == chosen_power])

# 3. Block-wise modules (signed Pearson, minModuleSize = 30)

net <- blockwiseModules(
  datExpr,
  power           = chosen_power,
  networkType     = "signed",
  corType         = "pearson",
  TOMType         = "signed",
  minModuleSize   = 30,
  reassignThreshold = 0,
  mergeCutHeight  = 0.25,
  numericLabels   = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs        = FALSE,
  verbose         = 0,
  maxBlockSize    = max(5000, ncol(datExpr) + 1))

module_colors <- labels2colors(net$colors)
log_msg("Module detection: %d modules (incl. grey), %d non-grey",
        length(unique(net$colors)),
        length(unique(net$colors[net$colors != 0])))
mod_tab <- table(module_colors)
log_msg("Module size summary:")
for (nm in names(sort(mod_tab, decreasing = TRUE)))
  cat(sprintf("  %-15s %5d\n", nm, mod_tab[nm]), file = LOG, append = TRUE)

# 4. Module eigengenes

MEs <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
MEs <- orderMEs(MEs)
log_msg("Eigengenes: %d modules x %d samples", ncol(MEs), nrow(MEs))

# Module membership (kME): correlation of each protein with each module's
# eigengene. Used for hub-protein selection (top |kME| within each module).
kME <- as.data.frame(cor(datExpr, MEs, use = "pairwise.complete.obs"))
colnames(kME) <- paste0("kME_", sub("^ME", "", colnames(MEs)))

# Hub proteins per module: top 5 by |kME| within their own module
hubs <- bind_rows(lapply(setdiff(unique(module_colors), "grey"), \(col) {
  prots <- colnames(datExpr)[module_colors == col]
  if (length(prots) == 0) return(NULL)
  k_self <- kME[prots, paste0("kME_", col), drop = TRUE]
  ord <- order(-abs(k_self))
  tibble(module = col, uniprot_id = prots[ord][seq_len(min(5, length(ord)))],
         kME = k_self[ord][seq_len(min(5, length(ord)))])
}))

ann <- as.data.frame(dal$annotation) |> select(uniprot_id, gene, description)
hubs <- hubs |> left_join(ann, by = "uniprot_id")

# 5. Save network

saveRDS(list(
  net = net, module_colors = module_colors, MEs = MEs, kME = kME, hubs = hubs,
  chosen_power = chosen_power, sft_df = sft_df,
  sample_order = rownames(datExpr),
  gene_order   = colnames(datExpr),
  ann = ann),
  file.path(DAT, "wgcna_network.rds"))

log_msg("\nSaved: %s\nDone: %s",
        file.path(DAT, "wgcna_network.rds"), Sys.time())
message("WGCNA build done — see ", LOG)
