#!/usr/bin/env Rscript
# F02 — All-contrast volcano composite (uniform 2×2, YvO-F03 parity).
#   A Transplant_Mito (Mito−Ctl)   B Rescue_Mito+Phe (PHE_Mito−PHE)   [primary, top]
#   C Disease_Phe (PHE−Ctl)        D Interaction_Mito                 [secondary, bottom]
# Volcano-in-ring per contrast: central volcano + gene ticks + NES-coloured
# enrichment arcs (Up right / Down left, symmetric). Reads existing pipeline
# outputs + frozen rat fGSEA cache + additive rat caches; never re-runs 01-03.
#
# THREE composites are emitted, one per enrichment lens (specialized DB alone):
#   _goslim : Hallmark + GO Slim   (coarse functional themes; rat org.Rn.eg.db)
#   _corum  : CORUM complexes      (protein complexes; Human+Mouse→rat babelgene)
#   _mito   : MitoCarta            (mitochondrial pathways; Mouse→rat babelgene)
# All gene sets are RAT-specific (msigdbr rat orthologs / org.Rn.eg.db / babelgene).
# Additive caches built by: build_fgsea_goslim_h9c2.R, build_fgsea_corum_h9c2.R.
#
# Biological framing verified against literature (see c_data/F02_citations.csv):
#   PHE -> H9c2 hypertrophy: Hahn 2014 (PMID 24794531); Jeong 2009 (PMID 19299911)
#   Mito transplant cardiac rescue: Masuzawa/Pacak 2013 (PMID 23355340);
#     Doulamis 2024 (PMID 39732955). DBs: MitoCarta3.0 Rath 2021 (33174596),
#     Hallmark Liberzon 2015 (26771021), Reactome Gillespie 2022 (34788843).

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr); library(readxl)
  library(stringr)
  library(ggplot2); library(ggrepel); library(ggforce); library(patchwork); library(cowplot)
})

source(here::here("05_Figures", "shared", "config.R"))
source(here::here("04_Figures", "shared", "volcano_ring.R"))
source(here::here("04_Figures", "shared", "pathway_utils.R"))
source(here::here("04_Figures", "shared", "figure_supplement_helpers.R"))
source(here::here("04_Figures", "shared", "mitocarta_utils.R"))

# MitoCarta is the canonical Figure-2 lens (OXPHOS/mitoribosome is the headline
# of a mito-transplant rescue study); GO-Slim and CORUM go to the supplement.
CANONICAL_LENS <- "mito"

BASE    <- fig05_base("F02_volcano")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
PNL_PNG <- file.path(RPT_PNG, "panels")
SUP_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
SUP_PNG <- file.path(BASE, "b_reports", "supp", "png")
SUP_PNL <- file.path(SUP_PNG, "panels")
DAT     <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, PNL_PNG, SUP_PDF, SUP_PNG, SUP_PNL, DAT))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

lens_pdf_dir <- function(cfg) if (cfg$suffix == CANONICAL_LENS) RPT_PDF else SUP_PDF
lens_png_dir <- function(cfg) if (cfg$suffix == CANONICAL_LENS) RPT_PNG else SUP_PNG
lens_pnl_dir <- function(cfg) if (cfg$suffix == CANONICAL_LENS) PNL_PNG else SUP_PNL

dep_df    <- read_csv(P05$comb, show_col_types = FALSE)
fgsea_all <- read_csv(here::here("04_Figures", "shared", "fgsea_tstat_all_h9c2.csv"), show_col_types = FALSE)
rat_gene_sets <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))

# Additive rat caches (GO Slim, CORUM) — bound in so they're selectable as DBs.
load_extra <- function(csv) if (file.exists(csv)) read_csv(csv, show_col_types = FALSE) else NULL
goslim_rds <- here::here("04_Figures", "shared", "goslim_rat_gene_sets.rds")
corum_rds  <- here::here("04_Figures", "shared", "corum_rat_gene_sets.rds")
goslim_sets <- if (file.exists(goslim_rds)) readRDS(goslim_rds) else list()
corum_sets  <- if (file.exists(corum_rds))  readRDS(corum_rds)  else list()
fgsea_all <- bind_rows(
  fgsea_all,
  load_extra(here::here("04_Figures", "shared", "fgsea_goslim_h9c2.csv")),
  load_extra(here::here("04_Figures", "shared", "fgsea_corum_h9c2.csv")))

RING_N_EACH  <- 6           # ≤6 Up + ≤6 Down per panel → mirror-symmetric ring
RING_JACCARD <- 0.5
RING_PADJ    <- 0.05

# 0-row ring template (correct columns). Passed as ring_data_override for panels
# with NO FDR-significant terms, so make_volcano_ring renders the volcano with an
# EMPTY ring instead of falling back to its own (unfiltered) term selection.
EMPTY_RING <- {
  tt <- fgsea_all |>
    filter(contrast == "PHEvPHE_MITO", database == "Hallmark", padj < 0.05, size >= 15) |>
    arrange(padj) |> head(4)
  build_ring_180_split(tt, "PHEvPHE_MITO", fgsea_all, databases = "Hallmark")[0, ]
}
# Bare MitoCarta sub-compartment sets are localizations, not pathways — exclude.
MITO_COMPARTMENT_SETS <- c("MITOCARTA_IMM", "MITOCARTA_IMS",
                           "MITOCARTA_MATRIX", "MITOCARTA_OMM")

# Three enrichment-lens composites (specialized DB alone). set_pool feeds dedup.
COMPOSITE_CONFIGS <- list(
  goslim = list(suffix = "goslim", dbs = c("Hallmark", "GO Slim"), min_size = 15,
                set_pool = c(rat_gene_sets$Hallmark, goslim_sets),
                lens = "Hallmark + GO Slim (rat)"),
  corum  = list(suffix = "corum", dbs = c("CORUM"), min_size = 3,
                set_pool = corum_sets,
                lens = "CORUM protein complexes (rat)"),
  mito   = list(suffix = "mito", dbs = c("MitoCarta"), min_size = 10,
                set_pool = rat_gene_sets$MitoCarta,
                lens = "MitoCarta mitochondrial pathways (rat)"))

# Tidy a few over-long Reactome/CORUM/GO labels (keys = engine-cleaned text).
LABEL_SHORTEN <- c(
  "Cargo Recognition For Clathrin Mediated Endocytosis" = "Clathrin\nEndocytosis",
  "Assembly Of Collagen Fibrils & Other Multimeric Structures" = "Collagen Fibril\nAssembly",
  "Collagen Chain Trimerization"               = "Collagen\nTrimerization",
  "Processing Of Capped Intron Containing Pre mRNA" = "Pre-mRNA\nProcessing",
  "Respiratory Chain Complex I (Holoenzyme), Mitochondrial" = "Respiratory\nComplex I",
  "Respiratory Chain Complex I, Mitochondrial" = "Respiratory\nComplex I",
  "Mitochondrial Ribosome, Large Subunit"      = "Mitoribosome\n(Large)",
  "Mitochondrial Ribosome, Small Subunit"      = "Mitoribosome\n(Small)")
shorten_label <- function(x) {
  key <- gsub("\n", " ", x)
  out <- unname(LABEL_SHORTEN[key])
  ifelse(is.na(out), x, out)
}

pick_symmetric <- function(pool, n_each) {
  if (nrow(pool) == 0) return(pool[FALSE, ])
  pool <- arrange(pool, padj)
  bind_rows(filter(pool, NES > 0) |> slice_head(n = n_each),
            filter(pool, NES < 0) |> slice_head(n = n_each))
}

contrast_stats <- function(ctr, dbs) {
  pi_col <- paste0("pi_score_", ctr)
  n_dep  <- if (pi_col %in% names(dep_df)) sum(dep_df[[pi_col]] < H9C2_PI_THRESH, na.rm = TRUE) else 0
  rows   <- fgsea_all |> filter(contrast == ctr, database %in% dbs)
  sprintf("%d DEPs (Π < 0.05)  |  %d / %d pathways (FDR < 0.05)",
          n_dep, sum(rows$padj < RING_PADJ, na.rm = TRUE), sum(!is.na(rows$padj)))
}

build_ring_panel <- function(ctr, tag, cfg) {
  pi_col <- paste0("pi_score_", ctr)
  n_dep  <- if (pi_col %in% names(dep_df)) sum(dep_df[[pi_col]] < H9C2_PI_THRESH, na.rm = TRUE) else 0

  sig_pool <- fgsea_all |>
    filter(contrast == ctr, database %in% cfg$dbs, !is.na(padj), padj < RING_PADJ,
           size >= cfg$min_size, !pathway %in% MITO_COMPARTMENT_SETS) |>
    arrange(padj)
  n_pre <- nrow(sig_pool)
  if (n_pre > 1) {
    sig_pool <- deduplicate_enrichment(
      as.data.frame(sig_pool), pathways = cfg$set_pool,
      jaccard_cutoff = RING_JACCARD, cross_db = TRUE) |> as_tibble()
  }
  top_terms <- pick_symmetric(sig_pool, RING_N_EACH)
  n_path <- nrow(top_terms); n_up <- sum(top_terms$NES > 0); n_dn <- sum(top_terms$NES < 0)

  ring_data <- if (n_path == 0) EMPTY_RING else
    build_ring_180_split(top_terms, ctr, fgsea_all, databases = cfg$dbs)
  if (!is.null(ring_data) && nrow(ring_data) > 0 && "clean_label" %in% names(ring_data))
    ring_data$clean_label <- shorten_label(ring_data$clean_label)
  adaptive_gap <- if (!is.null(ring_data) && nrow(ring_data) > 0)
    0.7 + 0.3 * (max(ring_data$arc_r1_var, na.rm = TRUE) - 4.8) / 1.6 else 0.7

  p <- make_volcano_ring(
    de_df = dep_df, go_df = fgsea_all, contrast = ctr,
    title = NULL, contrast_title = contrast_brief(ctr),
    contrast_subtitle = sprintf("%s | %d DEPs, %d pathways",
                                CONTRAST_MATH_BRIEF[ctr], n_dep, n_path),
    databases = cfg$dbs, ring_data_override = ring_data,
    label_size = 2.7, label_gap = adaptive_gap, title_size = 5, subtitle_size = 3.5,
    point_size = 0.5, point_alpha = 0.55,
    count_label_size = scale_text(BASE_COUNT, 89) + 0.4,
    count_y_mult = 0.75, count_x_mult = 0.85,
    bg_color = unname(CONTRAST_COLORS[ctr]), bg_alpha = 0.20,
    show_legend = FALSE) +
    labs(tag = tag)

  if (!is.null(ring_data) && nrow(ring_data) > 0)
    write_csv(top_terms |> select(pathway, database, padj, NES, size, any_of(c("ES", "log2err"))),
              file.path(DAT, sprintf("panel_%s_ring_terms_%s.csv", tag, cfg$suffix)))

  ggsave(file.path(lens_pnl_dir(cfg), sprintf("MAIN_panel_%s_%s_%s.png", tag, tolower(ctr), cfg$suffix)),
         p, width = 89, height = 89, units = "mm", dpi = 300)
  message(sprintf("  [%s] panel %s (%s): %d pre-dedup -> %d shown (%d up, %d down)",
                  cfg$suffix, tag, contrast_brief(ctr), n_pre, n_path, n_up, n_dn))
  strip_for_composite(p)
}

make_composite <- function(cfg) {
  message(sprintf("=== F02 composite: %s ===", cfg$lens))
  pA <- build_ring_panel("CTLvMITO",     "A", cfg)
  pB <- build_ring_panel("PHEvPHE_MITO", "B", cfg)
  pC <- build_ring_panel("CTLvPHE",      "C", cfg)
  pD <- build_ring_panel("Interaction",  "D", cfg)

  pA <- pA + theme(plot.margin = margin(4, -9, 0, 9, "mm"))
  pB <- pB + theme(plot.margin = margin(4,  0, 0, 0, "mm"))
  pC <- pC + theme(plot.margin = margin(0, -9, 0, 9, "mm"))
  pD <- pD + theme(plot.margin = margin(0,  0, 0, 0, "mm"))
  composite <- ((pA | pB) / (pC | pD)) + plot_layout(heights = c(1, 1))

  nes_legend <- build_nes_legend_bar(text_size = 5, title_size = 5,
                                     bar_margin = margin(0, 0, 0, 0, "mm"))

  COMP_W <- 178; COMP_H <- 180
  txt <- composite_text_sizes(COMP_H)
  TAG_SZ <- txt$tag + 4; TTL_SZ <- txt$title + 2; SUB_SZ <- txt$subtitle + 2
  X_L <- 0.070; X_R <- 0.510; X_TTL <- 0.040
  Y_TOP <- 0.960; Y_BOT <- 0.505; SUB_OFF <- 0.022

  titles <- c(contrast_brief("CTLvMITO"), contrast_brief("PHEvPHE_MITO"),
              contrast_brief("CTLvPHE"),  contrast_brief("Interaction"))
  subs   <- c(contrast_stats("CTLvMITO", cfg$dbs), contrast_stats("PHEvPHE_MITO", cfg$dbs),
              contrast_stats("CTLvPHE", cfg$dbs),  contrast_stats("Interaction", cfg$dbs))
  tags <- LETTERS[1:4]
  xs <- c(X_L, X_R, X_L, X_R); ys <- c(Y_TOP, Y_TOP, Y_BOT, Y_BOT)

  composite <- ggdraw(composite)
  for (i in 1:4) {
    composite <- composite +
      draw_label(tags[i], x = xs[i], y = ys[i] + 0.002, size = TAG_SZ,
                 fontface = "bold", hjust = 0, vjust = 1) +
      draw_label(titles[i], x = xs[i] + X_TTL, y = ys[i], size = TTL_SZ,
                 fontface = "bold", hjust = 0, vjust = 1) +
      draw_label(subs[i], x = xs[i] + X_TTL, y = ys[i] - SUB_OFF, size = SUB_SZ,
                 fontface = "bold.italic", colour = "grey40", hjust = 0, vjust = 1)
  }
  composite <- composite +
    draw_label(paste0(cfg$lens, " | N=24 (n=6/group); Interaction underpowered"),
               x = 0.5, y = 0.012, size = SUB_SZ, fontface = "italic",
               colour = "grey45", hjust = 0.5, vjust = 0) +
    draw_plot(nes_legend, x = 0.35, y = 0.030, width = 0.30, height = 0.028)

  ggsave(file.path(lens_pdf_dir(cfg), sprintf("MAIN_F02_composite_%s.pdf", cfg$suffix)), composite,
         width = COMP_W, height = COMP_H, units = "mm", device = pdf_dev, limitsize = FALSE)
  ggsave(file.path(lens_png_dir(cfg), sprintf("MAIN_F02_composite_%s.png", cfg$suffix)), composite,
         width = COMP_W, height = COMP_H, units = "mm", dpi = 300, limitsize = FALSE)
  message(sprintf("  saved MAIN_F02_composite_%s.{pdf,png} (%s)", cfg$suffix,
                  if (cfg$suffix == CANONICAL_LENS) "MAIN" else "supp"))
}

for (cfg in COMPOSITE_CONFIGS) make_composite(cfg)

# Verified-citation log + supplementary workbook
citations <- tribble(
  ~claim, ~pmid, ~ref,
  "PHE induces hypertrophy in H9c2 cardiomyoblasts", "24794531", "Hahn et al. 2014, Cell Signal",
  "Catecholamine/PHE-induced cardiac hypertrophy in H9c2", "19299911", "Jeong et al. 2009, Exp Mol Med",
  "Autologous mitochondrial transplant protects heart from I/R", "23355340", "Masuzawa/Pacak et al. 2013, AJP-Heart",
  "Mito transplant normalizes I/R proteomic shift (rescue)", "39732955", "Doulamis et al. 2024, Sci Rep",
  "Pi-score significance (P.Value^|logFC|)", "22321699", "Xiao et al. 2014, Bioinformatics",
  "MitoCarta 3.0 mitochondrial inventory", "33174596", "Rath et al. 2021, Nucleic Acids Res",
  "MSigDB Hallmark gene-set collection", "26771021", "Liberzon et al. 2015, Cell Syst",
  "Reactome pathway knowledgebase", "34788843", "Gillespie et al. 2022, Nucleic Acids Res")
write_csv(citations, file.path(DAT, "F02_citations.csv"))

sheets <- lapply(c("CTLvMITO", "PHEvPHE_MITO", "CTLvPHE", "Interaction"), \(c) {
  r <- as.data.frame(read_excel(P05$dep_xlsx, sheet = c))
  r |> select(uniprot_id, gene, description, logFC, P.Value, adj.P.Val, pi_score, sig_pi)
})
ring_csvs <- list.files(DAT, pattern = "^panel_._ring_terms_.*\\.csv$", full.names = TRUE)
ring_sheets <- lapply(ring_csvs, \(f) list(name = substr(sub(".csv$", "", basename(f)), 1, 31),
                                           df = as.data.frame(read_csv(f, show_col_types = FALSE))))
build_workbook(
  file.path(DAT, "F02_supplementary.xlsx"),
  sheet_specs = c(
    list(list(name = "contrast_map", df = fig05_contrast_table()),
         list(name = "citations", df = as.data.frame(citations))),
    lapply(seq_along(sheets), \(i) list(
      name = sprintf("DEP_%s", c("Transplant", "Rescue", "Disease", "Interaction")[i]),
      df = sheets[[i]])),
    ring_sheets))

message("F02: mito lens -> b_reports/main; goslim + corum -> b_reports/supp")
