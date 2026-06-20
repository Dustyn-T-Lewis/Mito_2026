#!/usr/bin/env Rscript
# F06 MAIN — multi-pilot cluster framework. Each pilot is one figure with a
# per-cluster row layout (trajectory left | Hallmark ORA right, color-coded by
# cluster). Pilots 1-3 (this task): fuzzy c-means on group means, gated on
# p<0.05 / Pi<0.05 / FDR<0.10 in >=1 core contrast. Fixed c = 6, seed 42.
# Cluster-selection Dmin sweeps still emitted to b_reports/supp for record.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(readr); library(ggplot2); library(patchwork); library(scales)
  library(e1071); library(limma)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))
source(here::here("04_Figures_v2", "functions", "07_cluster_row_layout.R"))

BASE     <- here::here("04_Figures_v2", "06_Cluster")
MAIN_PDF <- file.path(BASE, "b_reports", "main", "pdf")
MAIN_PNG <- file.path(BASE, "b_reports", "main", "png")
SUPP_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
SUPP_PNG <- file.path(BASE, "b_reports", "supp", "png")
DAT      <- file.path(BASE, "c_data")
for (d in c(MAIN_PDF, MAIN_PNG, SUPP_PDF, SUPP_PNG, DAT))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

FIG_W      <- PANEL_MD       # 178 mm
FIXED_C    <- 6L
SEED       <- 42L
CORE_MEMB  <- 0.5
C_RANGE    <- 2:12
SEEDS_DMIN <- 1:5

# ---------------------------------------------------------------------------
# Inputs: matrices + tables
# ---------------------------------------------------------------------------
dal     <- readRDS(P05$imp_rds)
prot    <- dal$data
ann     <- dal$annotation
meta    <- dal$metadata
gene_v  <- ann$gene[match(rownames(prot), ann$uniprot_id)]
keep_g  <- !is.na(gene_v) & nzchar(gene_v)
prot    <- prot[keep_g, , drop = FALSE]
gene_v  <- gene_v[keep_g]
gene_mat <- limma::avereps(prot, ID = gene_v)
ALL_GENES <- rownames(gene_mat)

grp <- meta$Group[match(colnames(gene_mat), meta$Col_ID)]
grp <- factor(grp, levels = H9C2_GROUP_LEVELS)
group_mat <- vapply(H9C2_GROUP_LEVELS, function(g)
  rowMeans(gene_mat[, grp == g, drop = FALSE]), numeric(nrow(gene_mat)))
rownames(group_mat) <- ALL_GENES

comb_long <- read_csv(P05$comb, show_col_types = FALSE)
CORE      <- H9C2_CORE_CONTRASTS

# ---------------------------------------------------------------------------
# c-means engine (carried over, simplified)
# ---------------------------------------------------------------------------
standardise_genes <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.finite(rowSums(z)), , drop = FALSE]
}
mestimate_fuzzifier <- function(z) {
  N <- nrow(z); D <- ncol(z)
  1 + (1418 / N + 22.05) * D^(-2) +
    (12.33 / N + 0.243) * D^(-0.0406 * log(N) - 0.1134)
}
mfuzz_cmeans <- function(z, c, m, seed = SEED) {
  set.seed(seed)
  e1071::cmeans(z, centers = c, m = m, method = "cmeans", iter.max = 200)
}
mean_dmin <- function(z, c, m, seeds = SEEDS_DMIN) {
  mean(vapply(seeds, function(s) min(dist(mfuzz_cmeans(z, c, m, seed = s)$centers)),
              numeric(1)))
}

# ---------------------------------------------------------------------------
# Per-pilot runner (c-means flavor)
# ---------------------------------------------------------------------------
run_cmeans_pilot <- function(key, gene_set, gate_label) {
  message(sprintf("\n=== %s (n_genes = %d) ===", key, length(gene_set)))
  mat <- group_mat[intersect(gene_set, rownames(group_mat)), , drop = FALSE]
  if (nrow(mat) < FIXED_C * 2)
    stop(sprintf("pilot %s: too few genes (%d) for c = %d", key, nrow(mat), FIXED_C))
  z <- standardise_genes(mat)
  m <- mestimate_fuzzifier(z)

  # selection sweep (diagnostic only; not used to pick c)
  dmin_tbl <- tibble(c = C_RANGE,
                     mean_Dmin = vapply(C_RANGE, function(c) mean_dmin(z, c, m),
                                        numeric(1)))
  fit  <- mfuzz_cmeans(z, FIXED_C, m, seed = SEED)
  hard <- fit$cluster
  max_mem <- apply(fit$membership, 1, max)
  memb <- tibble(gene = rownames(z), cluster = as.integer(hard),
                 membership = as.numeric(max_mem),
                 core = max_mem > CORE_MEMB)

  pal <- cluster_palette(FIXED_C)
  rows <- lapply(sort(unique(memb$cluster)), function(cl) {
    g_in_cl <- memb$gene[memb$cluster == cl]
    z_cl   <- z[g_in_cl, , drop = FALSE]
    ora    <- run_hallmark_ora(genes = memb$gene[memb$cluster == cl & memb$core],
                               universe = ALL_GENES)
    color  <- pal[as.character(cl)]
    hdr    <- sprintf("Cluster %d  |  n = %d (core %d)  |  Hallmark ORA",
                      cl, length(g_in_cl), sum(memb$cluster == cl & memb$core))
    build_cluster_row(
      traj_plot = build_trajectory_panel(z_cl, cluster = cl,
                                         x_levels = H9C2_GROUP_LEVELS,
                                         x_lab = "condition (group mean)",
                                         color = color, kind = "line"),
      ora_plot  = build_ora_bar_panel(ora, color = color, max_n = 6),
      header_text = hdr, color = color)
  })

  fig <- stack_cluster_rows(rows,
    title    = sprintf("F06 %s — c = %d fuzzy c-means", key, FIXED_C),
    subtitle = sprintf("gate: %s in >=1 of {%s}; m = %.3f; rows = clusters; right = Hallmark top-6",
                       gate_label, paste(CORE, collapse = " / "), m))
  h_mm <- 32 + 32 * FIXED_C
  ggsave(file.path(MAIN_PDF, sprintf("MAIN_F06_%s.pdf", key)), fig,
         width = FIG_W, height = h_mm, units = "mm", device = pdf_dev, limitsize = FALSE)
  ggsave(file.path(MAIN_PNG, sprintf("MAIN_F06_%s.png", key)), fig,
         width = FIG_W, height = h_mm, units = "mm", dpi = 300, limitsize = FALSE)

  # supp: Dmin-vs-c diagnostic
  sup <- ggplot(dmin_tbl, aes(c, mean_Dmin)) +
    geom_line(color = "grey50", linewidth = 0.4) +
    geom_point(size = 1.4, color = "grey30") +
    geom_vline(xintercept = FIXED_C, color = "#D6604D",
               linetype = "dashed", linewidth = 0.4) +
    scale_x_continuous(breaks = C_RANGE) +
    labs(title = sprintf("F06 %s cluster-selection diagnostic", key),
         subtitle = sprintf("fixed c = %d (no auto-pick); m = %.3f", FIXED_C, m),
         x = "number of clusters (c)",
         y = "mean min centroid distance (Dmin)") + FIG_THEME
  ggsave(file.path(SUPP_PDF, sprintf("MAIN_F06_%s_selection.pdf", key)), sup,
         width = 120, height = 80, units = "mm", device = pdf_dev)
  ggsave(file.path(SUPP_PNG, sprintf("MAIN_F06_%s_selection.png", key)), sup,
         width = 120, height = 80, units = "mm", dpi = 300)

  # ORA across all clusters (for the workbook)
  ora_all <- bind_rows(lapply(sort(unique(memb$cluster)), function(cl) {
    g <- memb$gene[memb$cluster == cl & memb$core]
    o <- run_hallmark_ora(g, universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) return(NULL)
    mutate(o, cluster = cl)
  }))

  list(key = key, m = m, fixed_c = FIXED_C, n_genes = nrow(z),
       memb = memb, ora = ora_all, dmin_tbl = dmin_tbl)
}

# ---------------------------------------------------------------------------
# PILOTS 1-3 (c-means × p / Π / FDR)
# ---------------------------------------------------------------------------
PILOTS_CMEANS <- list(
  list(key = "pilot_p",   col = "P.Value",   threshold = 0.05,
       gate_label = "p < 0.05"),
  list(key = "pilot_pi",  col = "pi_score",  threshold = 0.05,
       gate_label = "Π < 0.05"),
  list(key = "pilot_fdr", col = "adj.P.Val", threshold = 0.10,
       gate_label = "FDR < 0.10")
)

results <- list()
for (p in PILOTS_CMEANS) {
  genes <- filter_sig_in_any_contrast(comb_long, col = p$col,
                                      threshold = p$threshold, contrasts = CORE)
  results[[p$key]] <- run_cmeans_pilot(key = p$key, gene_set = genes,
                                       gate_label = p$gate_label)
}

# ---------------------------------------------------------------------------
# Single supplementary workbook (Pilots 1-3 only at this stage)
# ---------------------------------------------------------------------------
overview <- tibble(
  Pilot = vapply(results, `[[`, character(1), "key"),
  Method = "fuzzy c-means on group means",
  Gate = vapply(PILOTS_CMEANS, `[[`, character(1), "gate_label"),
  N_genes = vapply(results, `[[`, integer(1), "n_genes"),
  Fuzzifier_m = vapply(results, function(r) round(r$m, 3), numeric(1)),
  Cluster_c = vapply(results, `[[`, integer(1), "fixed_c"))

sheet_specs <- list(list(
  name = "Overview", df = overview,
  role = "F06 pilots — run-level summary",
  contents = "pilot key, method, gate, gene count, fuzzifier m, cluster count c"))

for (r in results) {
  sheet_specs <- c(sheet_specs, list(
    list(name = paste0(r$key, "_membership"), df = r$memb,
         role = sprintf("Soft-cluster assignment for %s", r$key),
         contents = "gene, hard cluster, max membership, core flag (>0.5)"),
    list(name = paste0(r$key, "_ora"),
         df = if (is.null(r$ora)) tibble(cluster = integer(), pathway = character(),
                                         padj = numeric(), overlap = integer(),
                                         size = integer(), odds_ratio = numeric())
              else r$ora,
         role = sprintf("Hallmark ORA per cluster for %s", r$key),
         contents = "cluster, pathway, padj, overlap, size, odds_ratio")))
}

build_workbook(file.path(DAT, "F06_supplementary.xlsx"),
               figure_title = "F06 — multi-pilot cluster framework",
               sheet_specs = sheet_specs)

message(sprintf("F06 c-means pilots done: %s", paste(names(results), collapse = ", ")))
