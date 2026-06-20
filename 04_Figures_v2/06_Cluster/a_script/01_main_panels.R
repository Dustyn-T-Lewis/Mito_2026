#!/usr/bin/env Rscript
# F06 MAIN — multi-pilot cluster framework. Each pilot is one figure with a
# per-cluster row layout (trajectory left | Hallmark ORA right, color-coded by
# cluster). Pilots 1-3 (this task): fuzzy c-means on group means, gated on
# p<0.05 / Pi<0.05 / FDR<0.10 in >=1 core contrast. Fixed c = 6, seed 42.
# Cluster-selection Dmin sweeps still emitted to b_reports/supp for record.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(readr); library(ggplot2); library(patchwork)
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

# Wipe stale per-pilot outputs from previous designs so the directory only
# contains the current six pilots' figures.
PILOT_KEYS <- c("pilot_p","pilot_pi","pilot_fdr","pilot_wgcna","pilot_logfc","pilot_rrho2","pilot_rrho2_heatmap")
for (d in c(MAIN_PDF, MAIN_PNG, SUPP_PDF, SUPP_PNG)) {
  existing <- list.files(d, pattern = "^MAIN_F06_.*\\.(pdf|png)$", full.names = TRUE)
  valid_re <- paste0("MAIN_F06_(", paste(PILOT_KEYS, collapse = "|"), ")(_selection|_me_traits)?\\.(pdf|png)$")
  stale <- existing[!grepl(valid_re, basename(existing))]
  if (length(stale) > 0) {
    message(sprintf("Removing %d stale F06 output(s) from %s", length(stale), d))
    file.remove(stale)
  }
}

FIG_W      <- PANEL_MD       # 178 mm
DEFAULT_C  <- 6L             # fallback if no flattening detected in Dmin/gap
SEED       <- 42L
CORE_MEMB  <- 0.5
C_RANGE    <- 2:12
SEEDS_DMIN <- 1:5
GAP_B      <- 50L            # bootstraps for clusGap (Tibshirani 2001)

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
# standardise_genes() is defined in 07_cluster_row_layout.R (generic utility).
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
  if (nrow(mat) < 6L)
    stop(sprintf("pilot %s: too few genes (%d) for clustering", key, nrow(mat)))
  z <- standardise_genes(mat)
  m <- mestimate_fuzzifier(z)

  # Per-pilot c selection via Dmin elbow (Schwämmle & Jensen 2010, PMID 20880957)
  # capped at sqrt(N/2) per Mardia 1979 small-N rule.
  dmin_tbl <- tibble(c = C_RANGE,
                     mean_Dmin = vapply(C_RANGE, function(c) mean_dmin(z, c, m),
                                        numeric(1)))
  pick <- pick_c_dmin(dmin_tbl, n_genes = nrow(z),
                      drop_frac = 0.10, default_c = DEFAULT_C)
  pilot_c <- pick$c
  message(sprintf("  %s", pick$basis))

  fit  <- mfuzz_cmeans(z, pilot_c, m, seed = SEED)
  hard <- fit$cluster
  max_mem <- apply(fit$membership, 1, max)
  memb <- tibble(gene = rownames(z), cluster = as.integer(hard),
                 membership = as.numeric(max_mem),
                 core = max_mem > CORE_MEMB)

  pal <- cluster_palette(pilot_c)
  rows <- lapply(sort(unique(memb$cluster)), function(cl) {
    g_in_cl <- memb$gene[memb$cluster == cl]
    z_cl   <- z[g_in_cl, , drop = FALSE]
    core_g <- memb$gene[memb$cluster == cl & memb$core]
    hall   <- run_hallmark_ora(genes = core_g, universe = ALL_GENES)
    mito   <- run_mitocarta_ora(genes = core_g, universe = ALL_GENES)
    color  <- pal[as.character(cl)]
    hdr    <- sprintf("Cluster %d  |  n = %d (core %d)  |  Hallmark + MitoCarta",
                      cl, length(g_in_cl), sum(memb$cluster == cl & memb$core))
    build_cluster_row(
      traj_plot = build_trajectory_panel(z_cl,
                                         x_levels = H9C2_GROUP_LEVELS,
                                         x_lab = "condition (group mean)",
                                         color = color, kind = "line"),
      ora_plot  = build_ora_bar_panel(hall, color = color, max_n = 6,
                                      db_name = "Hallmark"),
      ora_plot2 = build_ora_bar_panel(mito, color = color, max_n = 6,
                                      db_name = "MitoCarta"),
      header_text = hdr, color = color)
  })

  fig <- stack_cluster_rows(rows,
    title    = sprintf("F06 %s — c = %d fuzzy c-means (%s)", key, pilot_c, pick$basis),
    subtitle = sprintf("gate: %s in >=1 of {%s}; m = %.3f; rows = clusters; middle = Hallmark top-6, right = MitoCarta top-6",
                       gate_label, paste(CORE, collapse = " / "), m))
  h_mm <- 32 + 32 * pilot_c
  ggsave(file.path(MAIN_PDF, sprintf("MAIN_F06_%s.pdf", key)), fig,
         width = FIG_W, height = h_mm, units = "mm", device = pdf_dev, limitsize = FALSE)
  ggsave(file.path(MAIN_PNG, sprintf("MAIN_F06_%s.png", key)), fig,
         width = FIG_W, height = h_mm, units = "mm", dpi = 300, limitsize = FALSE)

  # supp: Dmin-vs-c diagnostic with the chosen c highlighted
  sup <- ggplot(dmin_tbl, aes(c, mean_Dmin)) +
    geom_line(color = "grey50", linewidth = 0.4) +
    geom_point(size = 1.4, color = "grey30") +
    geom_vline(xintercept = pilot_c, color = "#D6604D",
               linetype = "dashed", linewidth = 0.4) +
    geom_point(data = dmin_tbl[dmin_tbl$c == pilot_c, ],
               color = "#D6604D", size = 2.8) +
    scale_x_continuous(breaks = C_RANGE) +
    labs(title = sprintf("F06 %s cluster-selection diagnostic", key),
         subtitle = sprintf("%s; m = %.3f", pick$basis, m),
         x = "number of clusters (c)",
         y = "mean min centroid distance (Dmin)") + FIG_THEME
  ggsave(file.path(SUPP_PDF, sprintf("MAIN_F06_%s_selection.pdf", key)), sup,
         width = 120, height = 80, units = "mm", device = pdf_dev)
  ggsave(file.path(SUPP_PNG, sprintf("MAIN_F06_%s_selection.png", key)), sup,
         width = 120, height = 80, units = "mm", dpi = 300)

  # ORA across all clusters (for the workbook) — Hallmark + MitoCarta separately
  ora_all <- bind_rows(lapply(sort(unique(memb$cluster)), function(cl) {
    g <- memb$gene[memb$cluster == cl & memb$core]
    o <- run_hallmark_ora(g, universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) return(NULL)
    mutate(o, cluster = cl)
  }))
  ora_mito_all <- bind_rows(lapply(sort(unique(memb$cluster)), function(cl) {
    g <- memb$gene[memb$cluster == cl & memb$core]
    o <- run_mitocarta_ora(g, universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) return(NULL)
    mutate(o, cluster = cl)
  }))

  list(key = key, m = m, pilot_c = pilot_c, n_genes = nrow(z),
       selection_basis = pick$basis,
       sheets = list(
         membership   = memb,
         ora          = ora_all,
         ora_mito     = ora_mito_all,
         dmin_tbl     = dmin_tbl))
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
  Cluster_c = vapply(results, `[[`, integer(1), "pilot_c"),
  Selection = vapply(results, `[[`, character(1), "selection_basis"))

sheet_specs <- list(list(
  name = "Overview", df = overview,
  role = "F06 pilots — run-level summary",
  contents = "pilot key, method, gate, gene count, fuzzifier m, cluster count c"))

for (r in results) {
  if (is.null(r$sheets)) next   # skip pilots with different structure (4-6)
  sh <- r$sheets
  empty_ora <- tibble(cluster = integer(), pathway = character(),
                      padj = numeric(), overlap = integer(),
                      size = integer(), odds_ratio = numeric())
  sheet_specs <- c(sheet_specs, list(
    list(name = paste0(r$key, "_membership"), df = sh$membership,
         role = sprintf("Soft-cluster assignment for %s", r$key),
         contents = "gene, hard cluster, max membership, core flag (>0.5)"),
    list(name = paste0(r$key, "_ora"),
         df = if (is.null(sh$ora)) empty_ora else sh$ora,
         role = sprintf("Hallmark ORA per cluster for %s", r$key),
         contents = "cluster, pathway, padj, overlap, size, odds_ratio"),
    list(name = paste0(r$key, "_ora_mito"),
         df = if (is.null(sh$ora_mito)) empty_ora else sh$ora_mito,
         role = sprintf("MitoCarta ORA per cluster for %s", r$key),
         contents = "cluster, pathway, padj, overlap, size, odds_ratio")))
}

# ---------------------------------------------------------------------------
# PILOT 4 — WGCNA modules
# ---------------------------------------------------------------------------
WGCNA_RDS <- here::here("04_Figures", "F05_modules", "c_data", "wgcna_network.rds")
if (file.exists(WGCNA_RDS)) {
  message("\n=== pilot_wgcna ===")
  w <- load_wgcna_modules(WGCNA_RDS)
  mods <- w$modules |> filter(.data$module != "grey")
  MEs  <- w$MEs

  # Build the per-sample (Group) indicator vectors over the 4 contrasts the
  # eigengenes are correlated against. Group ordering must match MEs rownames.
  me_meta <- tibble(Col_ID = rownames(MEs)) |>
    left_join(as_tibble(meta) |> select(Col_ID, Group), by = "Col_ID")
  contrast_pairs <- list(
    Disease    = c("Ctl", "PHE"),
    Transplant = c("Ctl", "Mito"),
    Rescue     = c("PHE", "PHE_Mito"))
  me_corr <- compute_me_contrast_correlations(MEs, me_meta, contrast_pairs)

  # Interaction contrast: (PHE_Mito - Mito) - (PHE - Ctl) = (-1, +1, +1, -1)
  # aligned to (Ctl, Mito, PHE, PHE_Mito). Cannot use compute_me_contrast_correlations
  # (binary indicator only), so compute per-module cor.test directly.
  int_order <- H9C2_GROUP_LEVELS   # c("Ctl","Mito","PHE","PHE_Mito")
  int_vec_map <- c(Ctl = -1, Mito = 1, PHE = 1, PHE_Mito = -1)
  int_indic <- int_vec_map[me_meta$Group]
  int_keep  <- !is.na(int_indic)
  if (sum(int_keep) >= 4) {
    int_rows <- lapply(colnames(MEs), function(mod) {
      ct <- suppressWarnings(cor.test(MEs[int_keep, mod], int_indic[int_keep],
                                      method = "pearson"))
      tibble(module = mod, contrast = "Interaction",
             r = unname(ct$estimate), p = ct$p.value)
    })
    me_corr <- bind_rows(me_corr, bind_rows(int_rows))
  }

  # r_floor 0.25 (vs helper default 0.40): with n = 24 samples and 11 modules,
  # the strict 0.40 floor collapses every module into "Other"; 0.25 surfaces
  # the full reversal cohort — turquoise (D=+0.34, R=-0.61), brown (D=-0.28, R=+0.74),
  # black (D=+0.27, R=-0.40), blue (D=-0.30, R=+0.49) — plus red as Concordant down
  # (D=-0.29, R=-0.51). At n=24 these correlations are noisier than they'd be in a
  # larger cohort, so treat the sign-pattern column as a ranking aid rather than an
  # inferential claim.
  signs <- me_corr |>
    pivot_wider(id_cols = module, names_from = contrast, values_from = r,
                names_prefix = "r_") |>
    mutate(sign_pattern = classify_module_sign_pattern(r_Disease, r_Rescue,
                                                       r_floor = 0.25))
  # Row order: Reversal first, then Concordant up/down, then Other, then by |r_Rescue|
  mod_order <- signs |>
    arrange(.data$sign_pattern, desc(abs(.data$r_Rescue))) |>
    pull(.data$module)

  pal <- setNames(viridis::turbo(length(mod_order) + 2)[seq_along(mod_order) + 1],
                  mod_order)

  rows <- lapply(mod_order, function(mod) {
    g_in <- mods$gene[mods$module == mod]
    g_in <- intersect(g_in, rownames(gene_mat))
    if (length(g_in) < 5) return(NULL)
    # Trajectory uses group means (4 condition centroids), not individual sample
    # profiles. Deliberate deviation from spec wording: n=24 sample-level lines
    # add noise that obscures module shape, and group means keep visual consistency
    # with the c-means pilots (Pilots 1-3) that also use group means.
    z_cl <- standardise_genes(group_mat[g_in, , drop = FALSE])
    color <- pal[mod]
    sp   <- signs$sign_pattern[match(mod, signs$module)]
    rD   <- signs$r_Disease[match(mod, signs$module)]
    rR   <- signs$r_Rescue [match(mod, signs$module)]
    rT   <- signs$r_Transplant[match(mod, signs$module)]
    hdr  <- sprintf("Module %s  |  n = %d  |  %s  |  r(D)=%.2f r(R)=%.2f r(T)=%.2f",
                    mod, length(g_in), sp, rD, rR, rT)
    hall <- run_hallmark_ora (g_in, universe = ALL_GENES)
    mito <- run_mitocarta_ora(g_in, universe = ALL_GENES)
    build_cluster_row(
      traj_plot = build_trajectory_panel(z_cl,
                                         x_levels = H9C2_GROUP_LEVELS,
                                         x_lab = "condition (group mean)",
                                         color = color, kind = "line"),
      ora_plot  = build_ora_bar_panel(hall, color = color, max_n = 6,
                                      db_name = "Hallmark"),
      ora_plot2 = build_ora_bar_panel(mito, color = color, max_n = 6,
                                      db_name = "MitoCarta"),
      header_text = hdr, color = color)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]

  fig <- stack_cluster_rows(rows,
    title    = "F06 pilot_wgcna — modules from F05 WGCNA artifact",
    subtitle = "rows = modules; ordered by Disease<->Rescue sign pattern (reversal first); middle = Hallmark top-6, right = MitoCarta top-6")
  h_mm <- 32 + 32 * length(rows)
  ggsave(file.path(MAIN_PDF, "MAIN_F06_pilot_wgcna.pdf"), fig,
         width = FIG_W, height = h_mm, units = "mm", device = pdf_dev, limitsize = FALSE)
  ggsave(file.path(MAIN_PNG, "MAIN_F06_pilot_wgcna.png"), fig,
         width = FIG_W, height = h_mm, units = "mm", dpi = 300, limitsize = FALSE)

  # supp ME-trait heatmap
  me_hm <- ggplot(me_corr, aes(.data$contrast, .data$module, fill = .data$r)) +
    geom_tile(color = "white", linewidth = 0.1) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 1.6, color = "grey15") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
    labs(title = "F06 pilot_wgcna — module eigengene × contrast indicator",
         x = NULL, y = NULL) + FIG_THEME +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = FIG_AXIS_TEXT))
  ggsave(file.path(SUPP_PDF, "MAIN_F06_pilot_wgcna_me_traits.pdf"), me_hm,
         width = 100, height = 6 + 4 * length(unique(me_corr$module)),
         units = "mm", device = pdf_dev, limitsize = FALSE)
  ggsave(file.path(SUPP_PNG, "MAIN_F06_pilot_wgcna_me_traits.png"), me_hm,
         width = 100, height = 6 + 4 * length(unique(me_corr$module)),
         units = "mm", dpi = 300, limitsize = FALSE)

  # ORA across all modules for the workbook — Hallmark + MitoCarta separately
  ora_w <- bind_rows(lapply(mod_order, function(mod) {
    g <- mods$gene[mods$module == mod]
    o <- run_hallmark_ora(g, universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) return(NULL)
    mutate(o, module = mod)
  }))
  ora_mito_w <- bind_rows(lapply(mod_order, function(mod) {
    g <- mods$gene[mods$module == mod]
    o <- run_mitocarta_ora(g, universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) return(NULL)
    mutate(o, module = mod)
  }))

  results$pilot_wgcna <- list(
    key = "pilot_wgcna",
    sheets = list(
      pilot_wgcna_membership = mods,
      pilot_wgcna_me_traits  = signs,
      pilot_wgcna_ora        = if (is.null(ora_w))      tibble() else ora_w,
      pilot_wgcna_ora_mito   = if (is.null(ora_mito_w)) tibble() else ora_mito_w))
} else {
  warning("WGCNA artifact not found at ", WGCNA_RDS, " — skipping pilot_wgcna")
}

if (!is.null(results$pilot_wgcna)) {
  wg <- results$pilot_wgcna$sheets
  sheet_specs <- c(sheet_specs,
    list(list(name = "pilot_wgcna_membership", df = wg$pilot_wgcna_membership,
              role = "Gene -> WGCNA module assignment (grey excluded)",
              contents = "gene, module (color label from F05 build)")),
    list(list(name = "pilot_wgcna_me_traits", df = wg$pilot_wgcna_me_traits,
              role = "Module eigengene Pearson r vs contrast indicator vectors",
              contents = "module, r_Disease, r_Transplant, r_Rescue, r_Interaction, sign_pattern")),
    list(list(name = "pilot_wgcna_ora", df = wg$pilot_wgcna_ora,
              role = "Hallmark ORA per module",
              contents = "module, pathway, padj, overlap, size, odds_ratio")),
    list(list(name = "pilot_wgcna_ora_mito", df = wg$pilot_wgcna_ora_mito,
              role = "MitoCarta ORA per module",
              contents = "module, pathway, padj, overlap, size, odds_ratio")))
  # Add a row to Overview
  overview2 <- tibble(
    Pilot = "pilot_wgcna", Method = "WGCNA modules (existing artifact)",
    Gate = "all genes (network gate)",
    N_genes = nrow(wg$pilot_wgcna_membership),
    Fuzzifier_m = NA_real_,
    Cluster_c = length(unique(wg$pilot_wgcna_membership$module)))
  sheet_specs[[1]]$df <- bind_rows(sheet_specs[[1]]$df, overview2)
}

# ---------------------------------------------------------------------------
# PILOT 5 — k-means on per-protein 4-D logFC vector
# ---------------------------------------------------------------------------
message("\n=== pilot_logfc ===")
comb_wide <- load_combined_wide()
lf_cols  <- paste0("logFC_", CORE)
stopifnot(all(lf_cols %in% names(comb_wide)))
lf <- comb_wide |>
  select(gene, all_of(lf_cols)) |>
  filter(!is.na(.data$gene), nzchar(.data$gene),
         if_all(all_of(lf_cols), ~ !is.na(.x)))
lf_mat <- as.matrix(lf[, lf_cols])
rownames(lf_mat) <- lf$gene
if (nrow(lf_mat) < 12L)
  stop(sprintf("pilot_logfc: too few genes (%d) for clustering", nrow(lf_mat)))
# k selection: gap statistic (Tibshirani 2001 doi:10.1111/1467-9868.00293)
# via cluster::clusGap with firstSEmax — smallest k whose gap is within 1 SE
# of the maximum, the published anti-overfit choice for k-means on omics.
gap_pick <- pick_c_gap(lf_mat, k_range = 2:10, B = GAP_B, seed = SEED,
                       default_c = DEFAULT_C)
pilot_c_lf <- gap_pick$c
message(sprintf("pilot_logfc: %s", gap_pick$basis))
set.seed(SEED)
km <- kmeans(lf_mat, centers = pilot_c_lf, nstart = 50, iter.max = 100)

# Quadrant-style label from (logFC_Disease, logFC_Rescue) centroid signs.
disease_idx <- which(lf_cols == "logFC_CTLvPHE")
rescue_idx  <- which(lf_cols == "logFC_PHEvPHE_MITO")
quadrant_label <- function(dis, res) {
  if (abs(dis) < 0.1 && abs(res) < 0.1) return("Neutral")
  if (dis > 0 && res < 0) return("Reversed Down")     # disease up, rescue brings down
  if (dis < 0 && res > 0) return("Reversed Up")
  if (dis > 0 && res > 0) return("Concordant Up")
  if (dis < 0 && res < 0) return("Concordant Down")
  "Other"
}
centroids <- as_tibble(km$centers, rownames = "cluster") |>
  mutate(label = mapply(quadrant_label, km$centers[, disease_idx],
                                         km$centers[, rescue_idx]),
         n = km$size)

pal <- cluster_palette(pilot_c_lf)
rows <- lapply(seq_len(pilot_c_lf), function(cl) {
  g_in <- rownames(lf_mat)[km$cluster == cl]
  color <- pal[as.character(cl)]
  # trajectory = mean-logFC bar per contrast (kind = "barlogfc")
  z_cl  <- lf_mat[g_in, , drop = FALSE]   # already in logFC space
  colnames(z_cl) <- gsub("^logFC_", "", colnames(z_cl))
  univ  <- unique(rownames(lf_mat))
  hall  <- run_hallmark_ora (g_in, universe = univ)
  mito  <- run_mitocarta_ora(g_in, universe = univ)
  cent  <- centroids[centroids$cluster == as.character(cl), ]
  hdr   <- sprintf("Cluster %d  |  n = %d  |  %s",
                   cl, cent$n, cent$label)
  build_cluster_row(
    traj_plot = build_trajectory_panel(z_cl,
                                       x_levels = colnames(z_cl),
                                       x_lab = "contrast (centroid logFC)",
                                       color = color, kind = "barlogfc"),
    ora_plot  = build_ora_bar_panel(hall, color = color, max_n = 6,
                                    db_name = "Hallmark"),
    ora_plot2 = build_ora_bar_panel(mito, color = color, max_n = 6,
                                    db_name = "MitoCarta"),
    header_text = hdr, color = color)
})

fig <- stack_cluster_rows(rows,
  title    = sprintf("F06 pilot_logfc — k-means on per-protein 4-D logFC vector (c = %d; %s)",
                     pilot_c_lf, gap_pick$basis),
  subtitle = "rows = clusters; cluster labels derived from (Disease, Rescue) sign quadrant; middle = Hallmark top-6, right = MitoCarta top-6")
h_mm <- 32 + 32 * pilot_c_lf
ggsave(file.path(MAIN_PDF, "MAIN_F06_pilot_logfc.pdf"), fig,
       width = FIG_W, height = h_mm, units = "mm", device = pdf_dev, limitsize = FALSE)
ggsave(file.path(MAIN_PNG, "MAIN_F06_pilot_logfc.png"), fig,
       width = FIG_W, height = h_mm, units = "mm", dpi = 300, limitsize = FALSE)

# supp: gap-statistic diagnostic — gap(k) ± 1 SE bars with chosen c highlighted
gap_tbl <- gap_pick$gap_tbl
sup <- ggplot(gap_tbl, aes(.data$k, .data$gap)) +
  geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim),
                color = "grey55", width = 0.2, linewidth = 0.3) +
  geom_line(color = "grey50", linewidth = 0.4) +
  geom_point(size = 1.4, color = "grey30") +
  geom_vline(xintercept = pilot_c_lf, color = "#D6604D",
             linetype = "dashed", linewidth = 0.4) +
  geom_point(data = gap_tbl[gap_tbl$k == pilot_c_lf, ],
             color = "#D6604D", size = 2.8) +
  scale_x_continuous(breaks = 2:10) +
  labs(title = "F06 pilot_logfc cluster-selection diagnostic",
       subtitle = sprintf("gap statistic (Tibshirani 2001); %s", gap_pick$basis),
       x = "number of clusters (k)", y = "gap(k) ± 1 SE") + FIG_THEME
ggsave(file.path(SUPP_PDF, "MAIN_F06_pilot_logfc_selection.pdf"), sup,
       width = 120, height = 80, units = "mm", device = pdf_dev)
ggsave(file.path(SUPP_PNG, "MAIN_F06_pilot_logfc_selection.png"), sup,
       width = 120, height = 80, units = "mm", dpi = 300)

memb_lf <- tibble(gene = rownames(lf_mat), cluster = as.integer(km$cluster))
univ_lf <- unique(rownames(lf_mat))
ora_lf <- bind_rows(lapply(seq_len(pilot_c_lf), function(cl) {
  g <- memb_lf$gene[memb_lf$cluster == cl]
  o <- run_hallmark_ora(g, universe = univ_lf)
  if (is.null(o) || nrow(o) == 0) return(NULL)
  mutate(o, cluster = cl)
}))
ora_lf_mito <- bind_rows(lapply(seq_len(pilot_c_lf), function(cl) {
  g <- memb_lf$gene[memb_lf$cluster == cl]
  o <- run_mitocarta_ora(g, universe = univ_lf)
  if (is.null(o) || nrow(o) == 0) return(NULL)
  mutate(o, cluster = cl)
}))

results$pilot_logfc <- list(
  key = "pilot_logfc",
  sheets = list(
    pilot_logfc_membership = memb_lf,
    pilot_logfc_centroids  = centroids,
    pilot_logfc_ora        = if (is.null(ora_lf))      tibble() else ora_lf,
    pilot_logfc_ora_mito   = if (is.null(ora_lf_mito)) tibble() else ora_lf_mito))

if (!is.null(results$pilot_logfc)) {
  lg <- results$pilot_logfc$sheets
  sheet_specs <- c(sheet_specs,
    list(list(name = "pilot_logfc_membership", df = lg$pilot_logfc_membership,
              role = "Gene -> k-means cluster on per-protein 4-D logFC vector",
              contents = "gene, cluster (1..c)")),
    list(list(name = "pilot_logfc_centroids", df = lg$pilot_logfc_centroids,
              role = "k-means centroids in logFC space + quadrant label",
              contents = "cluster, logFC_<contrast>... mean centroid, quadrant label, n")),
    list(list(name = "pilot_logfc_ora", df = lg$pilot_logfc_ora,
              role = "Hallmark ORA per logFC cluster",
              contents = "cluster, pathway, padj, overlap, size, odds_ratio")),
    list(list(name = "pilot_logfc_ora_mito", df = lg$pilot_logfc_ora_mito,
              role = "MitoCarta ORA per logFC cluster",
              contents = "cluster, pathway, padj, overlap, size, odds_ratio")))
  sheet_specs[[1]]$df <- bind_rows(
    sheet_specs[[1]]$df,
    tibble(Pilot = "pilot_logfc", Method = "k-means on per-protein 4-D logFC vector",
           Gate = "all genes with non-NA logFC across core contrasts",
           N_genes = nrow(lg$pilot_logfc_membership),
           Fuzzifier_m = NA_real_, Cluster_c = pilot_c_lf))
}

# ---------------------------------------------------------------------------
# PILOT 6 — RRHO2 Disease <-> Rescue threshold-free concordance map
# ---------------------------------------------------------------------------
if (requireNamespace("RRHO2", quietly = TRUE)) {
  message("\n=== pilot_rrho2 ===")
  cw <- load_combined_wide()
  # Rank by signed limma moderated t (Smyth 2004, doi:10.2202/1544-6115.1027) —
  # variance-stabilized at small n and the same statistic fgsea consumes upstream,
  # so RRHO2 and fgsea agree on which gene "leads" each contrast. Replaces the
  # original Cahill 2018 (PMID 29442049) signed -log10(P) ranking; Cahill notes
  # the rank, not the metric, drives the hypergeometric tail.
  rank_vec <- function(tstat) {
    tstat[!is.finite(tstat)] <- 0
    tstat
  }
  d_rank <- rank_vec(cw$t_CTLvPHE)
  r_rank <- rank_vec(cw$t_PHEvPHE_MITO)
  keep <- !is.na(d_rank) & !is.na(r_rank) & !is.na(cw$gene) & nzchar(cw$gene)
  rrho_df <- slice_max(
    tibble(gene = cw$gene[keep],
           score_d = d_rank[keep],
           score_r = r_rank[keep],
           abs_d   = abs(d_rank[keep])),
    abs_d, n = 1, by = gene, with_ties = FALSE)
  list1 <- data.frame(gene = rrho_df$gene, score = rrho_df$score_d)
  list2 <- data.frame(gene = rrho_df$gene, score = rrho_df$score_r)
  set.seed(SEED)
  rr <- RRHO2::RRHO2_initialize(list1, list2,
                                labels = c("Disease (signed mod-t)", "Rescue (signed mod-t)"),
                                log10.ind = TRUE,
                                stepsize = ceiling(sqrt(nrow(list1))),
                                boundary = 0.025)  # per spec; trims 2.5% from each tail

  # Heatmap (saved separately; RRHO2_heatmap uses base graphics so needs its own device)
  # Axis convention: Disease (x) and Rescue (y) both run from low rank (most UP,
  # top-left) to high rank (most DOWN, bottom-right).
  # UU = top-left (both up), DD = bottom-right (both down),
  # UD = top-right (Disease up / Rescue down), DU = bottom-left (Disease down / Rescue up).
  # RRHO2_heatmap uses fields::image.plot internally, which leaves the device
  # active on the color-bar panel, so plain text() after the call writes onto the
  # legend strip, not the heatmap. Reserving outer margins before the call and
  # using mtext(outer=TRUE) puts the four corner labels on the device's outer
  # frame instead, independent of whatever sub-panel is active.
  .add_rrho2_quadrant_labels <- function() {
    mtext("UU", side = 3, line = 0.4, adj = 0.04, outer = TRUE,
          font = 2, cex = 1.3, col = "grey15")
    mtext("UD  (reversed)", side = 3, line = 0.4, adj = 0.78, outer = TRUE,
          font = 2, cex = 1.3, col = "#B2182B")
    mtext("DU  (reversed)", side = 1, line = 0.4, adj = 0.04, outer = TRUE,
          font = 2, cex = 1.3, col = "#D6604D")
    mtext("DD", side = 1, line = 0.4, adj = 0.78, outer = TRUE,
          font = 2, cex = 1.3, col = "grey15")
  }
  hm_pdf <- file.path(MAIN_PDF, "MAIN_F06_pilot_rrho2_heatmap.pdf")
  hm_png <- file.path(MAIN_PNG, "MAIN_F06_pilot_rrho2_heatmap.png")
  pdf(hm_pdf, width = 7, height = 7)
  par(oma = c(2, 1, 2, 1))   # bottom, left, top, right outer-margin lines
  RRHO2::RRHO2_heatmap(rr)
  .add_rrho2_quadrant_labels()
  dev.off()
  png(hm_png, width = 1800, height = 1800, res = 300)
  par(oma = c(2, 1, 2, 1))
  RRHO2::RRHO2_heatmap(rr)
  .add_rrho2_quadrant_labels()
  dev.off()

  # Quadrant gene lists — extract consistent top-fraction overlap per quadrant.
  # rr$genelist_* is a sub-list with gene_list1_*, gene_list2_*, gene_list_overlap_*.
  # When the peak-overlap intersection is very small (common with sparse concordance),
  # fall back to all genes that rank in the top RRHO2_PCT of BOTH lists in the
  # direction implied by each quadrant. This yields biologically meaningful gene
  # sets for trajectory + ORA even when the strict peak-cell overlap is tiny.
  RRHO2_PCT <- 0.20   # top 20% of ranked list defines each quadrant boundary
  n_top <- ceiling(RRHO2_PCT * nrow(rrho_df))
  d_asc  <- order(rrho_df$score_d)            # ascending (most negative first)
  d_desc <- order(rrho_df$score_d, decreasing = TRUE)
  r_asc  <- order(rrho_df$score_r)
  r_desc <- order(rrho_df$score_r, decreasing = TRUE)
  top_d_up   <- rrho_df$gene[d_desc[seq_len(n_top)]]
  top_d_down <- rrho_df$gene[d_asc [seq_len(n_top)]]
  top_r_up   <- rrho_df$gene[r_desc[seq_len(n_top)]]
  top_r_down <- rrho_df$gene[r_asc [seq_len(n_top)]]
  # RRHO2 peak-overlap sets (use if available and large enough, else use pct sets)
  uu_peak <- rr$genelist_uu$gene_list_overlap_uu
  dd_peak <- rr$genelist_dd$gene_list_overlap_dd
  ud_peak <- rr$genelist_ud$gene_list_overlap_ud
  du_peak <- rr$genelist_du$gene_list_overlap_du
  # pick_genes returns a named list with $genes (character) and $source (label)
  # so downstream can record which genes came from RRHO2 peak overlap vs the
  # top-20% percentile fallback.
  pick_genes <- function(peak, pct_fallback) {
    if (length(peak) >= 5) {
      list(genes = as.character(peak), source = "rrho2_peak")
    } else {
      cands <- if (is.null(pct_fallback)) character(0) else as.character(pct_fallback)
      list(genes = cands, source = "pct_fallback")
    }
  }
  quad_res <- list(
    UU = pick_genes(uu_peak, intersect(top_d_up,   top_r_up)),
    DD = pick_genes(dd_peak, intersect(top_d_down, top_r_down)),
    UD = pick_genes(ud_peak, intersect(top_d_up,   top_r_down)),  # reversed
    DU = pick_genes(du_peak, intersect(top_d_down, top_r_up)))    # reversed
  quad_lists  <- lapply(quad_res, `[[`, "genes")
  quad_source <- vapply(quad_res, `[[`, character(1), "source")
  message(sprintf("RRHO2 quadrant sizes: UU=%d DD=%d UD=%d DU=%d",
    length(quad_lists$UU), length(quad_lists$DD),
    length(quad_lists$UD), length(quad_lists$DU)))
  pal_q <- c(UU = "#2E7D32", DD = "#1565C0", UD = "#B2182B", DU = "#D6604D")
  quad_role <- c(UU = "Concordant Up", DD = "Concordant Down",
                 UD = "Reversed (Disease Up / Rescue Down)",
                 DU = "Reversed (Disease Down / Rescue Up)")
  rows <- lapply(names(quad_lists), function(q) {
    g_in <- quad_lists[[q]]
    if (length(g_in) < 5) return(NULL)
    g_in <- intersect(g_in, rownames(group_mat))
    if (length(g_in) < 5) return(NULL)
    z_cl <- standardise_genes(group_mat[g_in, , drop = FALSE])
    color <- pal_q[q]
    hall  <- run_hallmark_ora (g_in, universe = ALL_GENES)
    mito  <- run_mitocarta_ora(g_in, universe = ALL_GENES)
    hdr   <- sprintf("Quadrant %s  |  n = %d  |  %s", q, length(g_in), quad_role[q])
    build_cluster_row(
      traj_plot = build_trajectory_panel(z_cl,
                                         x_levels = H9C2_GROUP_LEVELS,
                                         x_lab = "condition (group mean)",
                                         color = color, kind = "line"),
      ora_plot  = build_ora_bar_panel(hall, color = color, max_n = 6,
                                      db_name = "Hallmark"),
      ora_plot2 = build_ora_bar_panel(mito, color = color, max_n = 6,
                                      db_name = "MitoCarta"),
      header_text = hdr, color = color)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]

  fig <- stack_cluster_rows(rows,
    title    = "F06 pilot_rrho2 — Disease<->Rescue quadrants (per-quadrant ORA)",
    subtitle = "heatmap saved separately; rows = significant RRHO2 quadrants; middle = Hallmark top-6, right = MitoCarta top-6")
  h_mm <- 32 + 32 * length(rows)
  ggsave(file.path(MAIN_PDF, "MAIN_F06_pilot_rrho2.pdf"), fig,
         width = FIG_W, height = h_mm, units = "mm", device = pdf_dev, limitsize = FALSE)
  ggsave(file.path(MAIN_PNG, "MAIN_F06_pilot_rrho2.png"), fig,
         width = FIG_W, height = h_mm, units = "mm", dpi = 300, limitsize = FALSE)

  # Workbook payloads — Hallmark + MitoCarta separately
  genelist_long <- bind_rows(lapply(names(quad_lists), function(q)
    tibble(quadrant = q, role = quad_role[q], gene = quad_lists[[q]],
           source = quad_source[q])))
  ora_rr <- bind_rows(lapply(names(quad_lists), function(q) {
    o <- run_hallmark_ora(quad_lists[[q]], universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) return(NULL)
    mutate(o, quadrant = q, role = quad_role[q])
  }))
  ora_rr_mito <- bind_rows(lapply(names(quad_lists), function(q) {
    o <- run_mitocarta_ora(quad_lists[[q]], universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) return(NULL)
    mutate(o, quadrant = q, role = quad_role[q])
  }))

  results$pilot_rrho2 <- list(
    key = "pilot_rrho2",
    sheets = list(
      pilot_rrho2_genelists = genelist_long,
      pilot_rrho2_ora       = if (is.null(ora_rr))      tibble() else ora_rr,
      pilot_rrho2_ora_mito  = if (is.null(ora_rr_mito)) tibble() else ora_rr_mito))
} else {
  warning("RRHO2 package missing — skipping pilot_rrho2")
}

if (!is.null(results$pilot_rrho2)) {
  rr <- results$pilot_rrho2$sheets
  sheet_specs <- c(sheet_specs,
    list(list(name = "pilot_rrho2_genelists", df = rr$pilot_rrho2_genelists,
              role = "Per-quadrant gene lists from RRHO2 (UU/DD/UD/DU)",
              contents = "quadrant (UU=concordant up, DD=concordant down, UD/DU=reversed), role, gene, source (rrho2_peak=from gene_list_overlap_*, pct_fallback=top-20% rank intersection)")),
    list(list(name = "pilot_rrho2_ora", df = rr$pilot_rrho2_ora,
              role = "Hallmark ORA per RRHO2 quadrant",
              contents = "quadrant, role, pathway, padj, overlap, size, odds_ratio")),
    list(list(name = "pilot_rrho2_ora_mito", df = rr$pilot_rrho2_ora_mito,
              role = "MitoCarta ORA per RRHO2 quadrant",
              contents = "quadrant, role, pathway, padj, overlap, size, odds_ratio")))
  sheet_specs[[1]]$df <- bind_rows(
    sheet_specs[[1]]$df,
    tibble(Pilot = "pilot_rrho2",
           Method = "RRHO2 threshold-free Disease<->Rescue map",
           Gate = "all genes with non-NA P.Value+logFC in Disease and Rescue",
           N_genes = nrow(rrho_df),   # input universe, not long-format output rows
           Fuzzifier_m = NA_real_, Cluster_c = 4L))
}

build_workbook(file.path(DAT, "F06_supplementary.xlsx"),
               figure_title = "F06 — multi-pilot cluster framework",
               sheet_specs = sheet_specs)

message(sprintf("F06 c-means pilots done: %s", paste(names(results), collapse = ", ")))
