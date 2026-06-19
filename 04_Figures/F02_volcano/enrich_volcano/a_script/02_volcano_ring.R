#!/usr/bin/env Rscript
# Mito volcano rings — clean YvO-F03 engine + enrichVolcano dedup. One clear
# parameter set, no per-DB cap.
#
#   DBs        Hallmark + Reactome + MitoCarta + GO Slim  (curated backbone +
#              detailed + mito-specific + broad GO themes)
#   rank       moderated t (limma)            ranking already in the fgsea cache
#   sig DEP    Xiao pi-score < 0.05
#   sig path   BH FDR < 0.05, gene-set size >= 10
#   dedup      drop a pathway whose leading edge is >50% contained in a better-
#              ranked one (overlap coefficient >= 0.5), within each database then
#              across — the same rule as F01/F03-F06 (Merico 2010 EnrichmentMap)
#   ring       top 12 pathways by FDR per panel (any up/down mix)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2)
  library(patchwork); library(cowplot); library(enrichVolcano)
})
source(here::here("04_Figures", "shared", "config.R"))
source(here::here("04_Figures", "shared", "volcano_ring.R"))
source(here::here("04_Figures", "shared", "figure_supplement_helpers.R"))

OUT  <- here::here("04_Figures", "F02_volcano", "enrich_volcano", "b_reports")
DBS  <- c("Hallmark", "Reactome", "MitoCarta", "GO Slim")
DROP <- c("MITOCARTA_ALL", "MITOCARTA_IMM", "MITOCARTA_IMS", "MITOCARTA_MATRIX", "MITOCARTA_OMM")  # localizations, not pathways
PANELS <- list(  # contrast, brief title, tag — YvO-F03 panel order
  list(c = "CTLvMITO",     t = "Transplant (Mito − Ctl)",  tag = "A"),
  list(c = "PHEvPHE_MITO", t = "Rescue (PHE_Mito − PHE)",   tag = "B"),
  list(c = "CTLvPHE",      t = "Disease (PHE − Ctl)",       tag = "C"),
  list(c = "Interaction",  t = "Interaction (Mito × PHE)",  tag = "D"))

# Tidy the few long Reactome / MitoCarta leaf names (keys = engine-cleaned text).
SHORTEN <- c("Cargo Recognition For Clathrin Mediated Endocytosis"        = "Clathrin\nEndocytosis",
             "Assembly Of Collagen Fibrils & Other Multimeric Structures" = "Collagen Fibril\nAssembly",
             "Binding And Uptake Of Ligands By Scavenger Receptors"       = "Scavenger\nReceptors",
             "Collagen Chain Trimerization"                               = "Collagen\nTrimerization",
             "Small Molecule Transport SLC25A family"                     = "SLC25A\nCarriers",
             "ROS and Glutathione Metabolism"                             = "ROS /\nGlutathione",
             "Branched Chain AA Metabolism"                               = "Branched-Chain\nAA Metab")
shorten_label <- function(x) { k <- gsub("\n", " ", x); o <- unname(SHORTEN[k]); ifelse(is.na(o), x, o) }

fgsea_all <- read_csv(here::here("04_Figures", "shared", "fgsea_tstat_all_h9c2.csv"),
                      show_col_types = FALSE) |>
  filter(database %in% DBS, !pathway %in% DROP)
dep_df <- load_combined_wide()

dedup_kept_terms <- function(ctr) ev_collapse(
  fgsea_all |> filter(contrast == ctr) |>
    mutate(leading_edge = leadingEdge, mode = "fgsea", direction = ifelse(NES > 0, "up", "down")),
  method = "jaccard", similarity = "overlap", cutoff = 0.5,
  scope = c("within_db", "cross_db"), sig_threshold = 0.05) |> filter(dedup_kept)

build_panel <- function(p) {
  pool    <- fgsea_all |> filter(contrast == p$c, database %in% DBS, !is.na(padj), size >= 10)
  n_total <- nrow(pool)                  # pathways tested in these DBs
  n_sig   <- sum(pool$padj < 0.05)       # significant (FDR < 0.05)
  top  <- dedup_kept_terms(p$c) |>       # ring = top 12 of the non-redundant significant set
    filter(database %in% DBS, padj < 0.05, size >= 10) |> arrange(padj) |> slice_head(n = 12)
  ring <- build_ring_180_split(top, p$c, fgsea_all, databases = DBS)
  if (nrow(ring) > 0 && "clean_label" %in% names(ring)) ring$clean_label <- shorten_label(ring$clean_label)
  gap  <- if (nrow(ring) > 0) 0.7 + 0.3 * (max(ring$arc_r1_var, na.rm = TRUE) - 4.8) / 1.6 else 0.7
  n_dep <- sum(dep_df[[paste0("pi_score_", p$c)]] < H9C2_PI_THRESH, na.rm = TRUE)
  panel <- make_volcano_ring(
    de_df = dep_df, go_df = fgsea_all, contrast = p$c,
    title = NULL, contrast_title = p$t,
    contrast_subtitle = sprintf("%s | %d DEPs, %d of %d pathways",
                                CONTRAST_MATH_BRIEF[p$c], n_dep, n_sig, n_total),
    databases = DBS, ring_data_override = ring,
    label_size = 2.7, label_gap = gap, point_size = 0.5, point_alpha = 0.55,
    count_label_size = scale_text(BASE_COUNT, 89) + 0.4,
    count_y_mult = 0.75, count_x_mult = 0.85,
    bg_color = unname(CONTRAST_COLORS[p$c]), bg_alpha = 0.20, show_legend = FALSE) +
    labs(tag = p$tag)
  list(plot = strip_for_composite(panel), n_sig = n_sig, n_total = n_total)
}

panels <- lapply(PANELS, build_panel)
pA <- panels[[1]]$plot + theme(plot.margin = margin(4, -9, 0, 9, "mm"))
pB <- panels[[2]]$plot + theme(plot.margin = margin(4,  0, 0, 0, "mm"))
pC <- panels[[3]]$plot + theme(plot.margin = margin(0, -9, 0, 9, "mm"))
pD <- panels[[4]]$plot + theme(plot.margin = margin(0,  0, 0, 0, "mm"))

COMP_W <- 178; COMP_H <- 180
txt <- composite_text_sizes(COMP_H)
xs <- c(0.070, 0.510, 0.070, 0.510); ys <- c(0.960, 0.960, 0.505, 0.505)
composite <- ggdraw((pA | pB) / (pC | pD))
for (i in 1:4) {
  composite <- composite +
    draw_label(PANELS[[i]]$tag, x = xs[i], y = ys[i] + 0.002, size = txt$tag + 4,
               fontface = "bold", hjust = 0, vjust = 1) +
    draw_label(PANELS[[i]]$t, x = xs[i] + 0.040, y = ys[i], size = txt$title + 2,
               fontface = "bold", hjust = 0, vjust = 1) +
    draw_label(sprintf("%d DEPs (Π<0.05) | %d of %d pathways (FDR<0.05)",
                       sum(dep_df[[paste0("pi_score_", PANELS[[i]]$c)]] < H9C2_PI_THRESH, na.rm = TRUE),
                       panels[[i]]$n_sig, panels[[i]]$n_total),
               x = xs[i] + 0.040, y = ys[i] - 0.022, size = txt$subtitle + 2,
               fontface = "bold.italic", colour = "grey40", hjust = 0, vjust = 1)
}
composite <- composite +
  draw_plot(build_nes_legend_bar(text_size = 5, title_size = 5), 0.35, 0.025, 0.30, 0.028)

ggsave(file.path(OUT, "MAIN_F02_volcano_ring.pdf"), composite, width = COMP_W, height = COMP_H,
       units = "mm", device = get_pdf_device(), limitsize = FALSE)
ggsave(file.path(OUT, "MAIN_F02_volcano_ring.png"), composite, width = COMP_W, height = COMP_H,
       units = "mm", dpi = 300, limitsize = FALSE)
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
message("Mito volcano rings (Hallmark + Reactome + MitoCarta + GO Slim, overlap dedup) -> ", OUT)
