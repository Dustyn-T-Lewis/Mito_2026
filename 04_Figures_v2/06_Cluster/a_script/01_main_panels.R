#!/usr/bin/env Rscript
# F06 MAIN — multi-pilot cluster framework. Six pilots, each emitting one figure
# with a per-cluster row layout (trajectory left | Hallmark ORA middle |
# MitoCarta ORA right). Pilots 1-3 are c-means at p / Π / FDR gates; pilot 4 is
# WGCNA modules; pilot 5 is k-means on per-protein logFC vectors; pilot 6 is
# RRHO2 Disease<->Rescue. The c-means runner lives here; pilots 4-6 live in
# sibling files sourced after the c-means block.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(e1071)
  library(limma)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))
source(here::here("04_Figures_v2", "functions", "07_cluster_row_layout.R"))

BASE <- here::here("04_Figures_v2", "06_Cluster")
MAIN_PDF <- file.path(BASE, "b_reports", "main", "pdf")
MAIN_PNG <- file.path(BASE, "b_reports", "main", "png")
SUPP_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
SUPP_PNG <- file.path(BASE, "b_reports", "supp", "png")
DAT <- file.path(BASE, "c_data")
for (d in c(MAIN_PDF, MAIN_PNG, SUPP_PDF, SUPP_PNG, DAT)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
pdf_dev <- get_pdf_device()

# Wipe stale per-pilot outputs from previous designs so the directory only
# contains the current six pilots' figures.
PILOT_KEYS <- c("pilot_p", "pilot_pi", "pilot_fdr", "pilot_wgcna", "pilot_logfc", "pilot_rrho2", "pilot_rrho2_heatmap")
for (d in c(MAIN_PDF, MAIN_PNG, SUPP_PDF, SUPP_PNG)) {
  existing <- list.files(d, pattern = "^MAIN_F06_.*\\.(pdf|png)$", full.names = TRUE)
  valid_re <- paste0("MAIN_F06_(", paste(PILOT_KEYS, collapse = "|"), ")(_selection|_me_traits)?\\.(pdf|png)$")
  stale <- existing[!grepl(valid_re, basename(existing))]
  if (length(stale)) {
    message(sprintf("Removing %d stale F06 output(s) from %s", length(stale), d))
    file.remove(stale)
  }
}

FIG_W <- PANEL_MD # 178 mm
DEFAULT_C <- 6L # fallback if no flattening detected in Dmin/gap
SEED <- 42L
CORE_MEMB <- 0.5
C_RANGE <- 2:12
SEEDS_DMIN <- 1:5
GAP_B <- 50L # bootstraps for clusGap (Tibshirani 2001)

# Inputs: matrices + tables
dal <- readRDS(P05$imp_rds)
prot <- dal$data
ann <- dal$annotation
meta <- dal$metadata
gene_v <- ann$gene[match(rownames(prot), ann$uniprot_id)]
keep_g <- !is.na(gene_v) & nzchar(gene_v)
prot <- prot[keep_g, , drop = FALSE]
gene_v <- gene_v[keep_g]
gene_mat <- limma::avereps(prot, ID = gene_v)
ALL_GENES <- rownames(gene_mat)

grp <- meta$Group[match(colnames(gene_mat), meta$Col_ID)]
grp <- factor(grp, levels = H9C2_GROUP_LEVELS)
group_mat <- vapply(H9C2_GROUP_LEVELS, function(g) {
  rowMeans(gene_mat[, grp == g, drop = FALSE])
}, numeric(nrow(gene_mat)))
rownames(group_mat) <- ALL_GENES

comb_long <- read_csv(P05$comb, show_col_types = FALSE)
CORE <- H9C2_CORE_CONTRASTS

# c-means engine — standardise_genes() is defined in 07_cluster_row_layout.R.
mestimate_fuzzifier <- function(z) {
  N <- nrow(z)
  D <- ncol(z)
  1 + (1418 / N + 22.05) * D^(-2) +
    (12.33 / N + 0.243) * D^(-0.0406 * log(N) - 0.1134)
}
mfuzz_cmeans <- function(z, c, m, seed = SEED) {
  set.seed(seed)
  e1071::cmeans(z, centers = c, m = m, method = "cmeans", iter.max = 200)
}
mean_dmin <- function(z, c, m, seeds = SEEDS_DMIN) {
  mean(vapply(
    seeds, function(s) min(dist(mfuzz_cmeans(z, c, m, seed = s)$centers)),
    numeric(1)
  ))
}

# Per-pilot runner (c-means flavor)
run_cmeans_pilot <- function(key, gene_set, gate_label) {
  message(sprintf("%s: n_genes = %d", key, length(gene_set)))
  mat <- group_mat[intersect(gene_set, rownames(group_mat)), , drop = FALSE]
  if (nrow(mat) < 6L) {
    stop(sprintf("pilot %s: too few genes (%d) for clustering", key, nrow(mat)))
  }
  z <- standardise_genes(mat)
  m <- mestimate_fuzzifier(z)

  # Per-pilot c selection via Dmin elbow (Schwämmle & Jensen 2010 PMID 20880957)
  # capped at sqrt(N/2) per Mardia, Kent & Bibby 1979 small-N rule.
  dmin_tbl <- tibble(
    c = C_RANGE,
    mean_Dmin = vapply(
      C_RANGE, function(c) mean_dmin(z, c, m),
      numeric(1)
    )
  )
  pick <- pick_c_dmin(dmin_tbl,
    n_genes = nrow(z),
    drop_frac = 0.10, default_c = DEFAULT_C
  )
  pilot_c <- pick$c
  message(sprintf("  %s", pick$basis))

  fit <- mfuzz_cmeans(z, pilot_c, m, seed = SEED)
  hard <- fit$cluster
  max_mem <- apply(fit$membership, 1, max)
  memb <- tibble(
    gene = rownames(z), cluster = as.integer(hard),
    membership = as.numeric(max_mem),
    core = max_mem > CORE_MEMB
  )

  pal <- cluster_palette(pilot_c)
  rows <- lapply(sort(unique(memb$cluster)), function(cl) {
    g_in_cl <- memb$gene[memb$cluster == cl]
    z_cl <- z[g_in_cl, , drop = FALSE]
    core_g <- memb$gene[memb$cluster == cl & memb$core]
    hall <- run_hallmark_ora(genes = core_g, universe = ALL_GENES)
    mito <- run_mitocarta_ora(genes = core_g, universe = ALL_GENES)
    color <- pal[as.character(cl)]
    hdr <- sprintf(
      "Cluster %d  |  n = %d (core %d)  |  Hallmark + MitoCarta",
      cl, length(g_in_cl), sum(memb$cluster == cl & memb$core)
    )
    build_cluster_row(
      traj_plot = build_trajectory_panel(z_cl,
        x_levels = H9C2_GROUP_LEVELS,
        x_lab = "condition (group mean)",
        color = color, kind = "line"
      ),
      ora_plot = build_ora_bar_panel(hall,
        color = color, max_n = 6,
        db_name = "Hallmark"
      ),
      ora_plot2 = build_ora_bar_panel(mito,
        color = color, max_n = 6,
        db_name = "MitoCarta"
      ),
      header_text = hdr, color = color
    )
  })

  fig <- stack_cluster_rows(rows,
    title = sprintf("F06 %s — c = %d fuzzy c-means (%s)", key, pilot_c, pick$basis),
    subtitle = sprintf(
      "gate: %s in >=1 of {%s}; m = %.3f; rows = clusters; middle = Hallmark top-6, right = MitoCarta top-6",
      gate_label, paste(CORE, collapse = " / "), m
    )
  )
  h_mm <- 32 + 32 * pilot_c
  ggsave(file.path(MAIN_PDF, sprintf("MAIN_F06_%s.pdf", key)), fig,
    width = FIG_W, height = h_mm, units = "mm", device = pdf_dev, limitsize = FALSE
  )
  ggsave(file.path(MAIN_PNG, sprintf("MAIN_F06_%s.png", key)), fig,
    width = FIG_W, height = h_mm, units = "mm", dpi = 300, limitsize = FALSE
  )

  # supp: Dmin-vs-c diagnostic with the chosen c highlighted
  sup <- ggplot(dmin_tbl, aes(c, mean_Dmin)) +
    geom_line(color = "grey50", linewidth = 0.4) +
    geom_point(size = 1.4, color = "grey30") +
    geom_vline(
      xintercept = pilot_c, color = "#D6604D",
      linetype = "dashed", linewidth = 0.4
    ) +
    geom_point(
      data = dmin_tbl[dmin_tbl$c == pilot_c, ],
      color = "#D6604D", size = 2.8
    ) +
    scale_x_continuous(breaks = C_RANGE) +
    labs(
      title = sprintf("F06 %s cluster-selection diagnostic", key),
      subtitle = sprintf("%s; m = %.3f", pick$basis, m),
      x = "number of clusters (c)",
      y = "mean min centroid distance (Dmin)"
    ) +
    FIG_THEME
  ggsave(file.path(SUPP_PDF, sprintf("MAIN_F06_%s_selection.pdf", key)), sup,
    width = 120, height = 80, units = "mm", device = pdf_dev
  )
  ggsave(file.path(SUPP_PNG, sprintf("MAIN_F06_%s_selection.png", key)), sup,
    width = 120, height = 80, units = "mm", dpi = 300
  )

  # ORA across all clusters (for the workbook) — Hallmark + MitoCarta separately
  ora_all <- bind_rows(lapply(sort(unique(memb$cluster)), function(cl) {
    g <- memb$gene[memb$cluster == cl & memb$core]
    o <- run_hallmark_ora(g, universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) {
      return(NULL)
    }
    mutate(o, cluster = cl)
  }))
  ora_mito_all <- bind_rows(lapply(sort(unique(memb$cluster)), function(cl) {
    g <- memb$gene[memb$cluster == cl & memb$core]
    o <- run_mitocarta_ora(g, universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) {
      return(NULL)
    }
    mutate(o, cluster = cl)
  }))

  list(
    key = key, m = m, pilot_c = pilot_c, n_genes = nrow(z),
    selection_basis = pick$basis,
    sheets = list(
      membership   = memb,
      ora          = ora_all,
      ora_mito     = ora_mito_all,
      dmin_tbl     = dmin_tbl
    )
  )
}

# PILOTS 1-3 (c-means × p / Π / FDR)
PILOTS_CMEANS <- list(
  list(
    key = "pilot_p", col = "P.Value", threshold = 0.05,
    gate_label = "p < 0.05"
  ),
  list(
    key = "pilot_pi", col = "pi_score", threshold = 0.05,
    gate_label = "Π < 0.05"
  ),
  list(
    key = "pilot_fdr", col = "adj.P.Val", threshold = 0.10,
    gate_label = "FDR < 0.10"
  )
)

results <- list()
for (p in PILOTS_CMEANS) {
  genes <- filter_sig_in_any_contrast(comb_long,
    col = p$col,
    threshold = p$threshold, contrasts = CORE
  )
  results[[p$key]] <- run_cmeans_pilot(
    key = p$key, gene_set = genes,
    gate_label = p$gate_label
  )
}

# Single supplementary workbook (Pilots 1-3 only at this stage)
overview <- tibble(
  Pilot = vapply(results, `[[`, character(1), "key"),
  Method = "fuzzy c-means on group means",
  Gate = vapply(PILOTS_CMEANS, `[[`, character(1), "gate_label"),
  N_genes = vapply(results, `[[`, integer(1), "n_genes"),
  Fuzzifier_m = vapply(results, function(r) round(r$m, 3), numeric(1)),
  Cluster_c = vapply(results, `[[`, integer(1), "pilot_c"),
  Selection = vapply(results, `[[`, character(1), "selection_basis")
)

sheet_specs <- list(list(
  name = "Pilot_summary", df = overview,
  role = "F06 pilots — run-level summary",
  contents = "pilot key, method, gate, gene count, fuzzifier m, cluster count c"
))

for (r in results) {
  if (is.null(r$sheets)) next
  sh <- r$sheets
  empty_ora <- tibble(
    cluster = integer(), pathway = character(),
    padj = numeric(), overlap = integer(),
    size = integer(), odds_ratio = numeric()
  )
  sheet_specs <- c(sheet_specs, list(
    list(
      name = paste0(r$key, "_membership"), df = sh$membership,
      role = sprintf("Soft-cluster assignment for %s", r$key),
      contents = "gene, hard cluster, max membership, core flag (>0.5)"
    ),
    list(
      name = paste0(r$key, "_ora"),
      df = if (is.null(sh$ora)) empty_ora else sh$ora,
      role = sprintf("Hallmark ORA per cluster for %s", r$key),
      contents = "cluster, pathway, padj, overlap, size, odds_ratio"
    ),
    list(
      name = paste0(r$key, "_ora_mito"),
      df = if (is.null(sh$ora_mito)) empty_ora else sh$ora_mito,
      role = sprintf("MitoCarta ORA per cluster for %s", r$key),
      contents = "cluster, pathway, padj, overlap, size, odds_ratio"
    )
  ))
}

# Pilots 4-6 live in sibling files and continue to grow `results` + `sheet_specs`.
source(here::here("04_Figures_v2", "06_Cluster", "a_script", "02_pilot_wgcna.R"))
source(here::here("04_Figures_v2", "06_Cluster", "a_script", "03_pilot_logfc.R"))
source(here::here("04_Figures_v2", "06_Cluster", "a_script", "04_pilot_rrho2.R"))

build_workbook(file.path(DAT, "F06_supplementary.xlsx"),
  figure_title = "F06 — multi-pilot cluster framework",
  sheet_specs = sheet_specs
)

message(sprintf("F06 done: %s", paste(names(results), collapse = ", ")))
