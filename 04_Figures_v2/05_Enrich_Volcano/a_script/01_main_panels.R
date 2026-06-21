#!/usr/bin/env Rscript
# F05 — Standalone per-contrast enrichment volcano-in-ring panels.
# Same per-panel pipeline as F02's build_ring_panel(), but each of the 4
# contrasts is emitted as its OWN square ring (not a 2x2 composite). A shared
# NES colour-bar legend is emitted once. Reads existing pipeline outputs +
# frozen rat fGSEA cache; never re-runs 01-03.
#
# Curated 5-DB lens (Hallmark + Reactome + KEGG + MitoCarta + GO Slim, rat
# orthologs; per Reimand 2019 PMID 30664679). Each ring shows the 12 most
# significant pathways (lowest padj) across all DBs, after EnrichmentMap
# combined-coefficient redundancy collapse (0.5*overlap + 0.5*jaccard >= 0.375;
# Merico 2010 PMID 21085593) in 03_pathway_enrichment_dedup_ora.R. Full
# per-contrast fGSEA tables (every tested pathway + FDR) go to the workbook.
#
# Filename contrast tags use the OLD combined-results names, lowercased:
#   ctlvphe       Disease     (PHE - Ctl)
#   ctlvmito      Transplant  (Mito - Ctl)
#   phevphe_mito  Rescue      (PHE_Mito - PHE)
#   interaction   Interaction (orthogonal 2x2)

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr)
  library(stringr); library(ggplot2); library(ggforce)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "05_volcano_ring_plot_builder.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))
source(here::here("04_Figures_v2", "functions", "04_mitocarta_lens_lookup.R"))

BASE    <- here::here("04_Figures_v2", "05_Enrich_Volcano")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
DAT     <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

dep_df        <- load_combined_wide()
fgsea_all     <- read_csv(here::here("04_Figures", "shared", "fgsea_tstat_all_h9c2.csv"), show_col_types = FALSE)
rat_gene_sets <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))

# Canonical 5-DB lens (Hallmark + Reactome + KEGG + MitoCarta + GO Slim), shared
# with the pathway-bar and cluster figures via CANONICAL_DBS / MITO_DROP_SETS.
POOL_DBS  <- CANONICAL_DBS
SET_POOL  <- do.call(c, unname(rat_gene_sets[POOL_DBS]))
RING_N    <- 12      # most significant pathways drawn per ring (split Up/Down by NES)
RING_PADJ <- 0.05
SIM_CUT   <- 0.375   # EnrichmentMap combined-coefficient redundancy cutoff
RING_MIN_SZ <- 10

# 4 contrasts, OLD name -> filename tag + role brief.
CONTRASTS <- tribble(
  ~ctr,            ~tag,           ~role,
  "CTLvPHE",       "ctlvphe",      "Disease",
  "CTLvMITO",      "ctlvmito",     "Transplant",
  "PHEvPHE_MITO",  "phevphe_mito", "Rescue",
  "Interaction",   "interaction",  "Interaction")

# 0-row ring template (correct columns) for panels with no sig terms, so
# make_volcano_ring draws an empty ring instead of its own fallback selection.
EMPTY_RING <- build_ring_180_split(
  head(arrange(filter(fgsea_all, contrast == "PHEvPHE_MITO", database == "Hallmark"), padj), 2),
  "PHEvPHE_MITO", fgsea_all, databases = "Hallmark")[0, ]

# Tidy a few over-long Reactome/GO labels (keys = engine-cleaned text).
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

# The n most significant pathways (lowest padj); the ring splits them Up/Down by NES.
pick_top <- function(pool, n) slice_head(arrange(pool, padj), n = n)

# Build + save one standalone ring; return its shown terms (with contrast tag).
build_ring <- function(ctr, tag, role) {
  pi_col <- paste0("pi_score_", ctr)
  n_dep  <- if (pi_col %in% names(dep_df)) sum(dep_df[[pi_col]] < H9C2_PI_THRESH, na.rm = TRUE) else 0

  sig_pool <- fgsea_all |>
    filter(contrast == ctr, database %in% POOL_DBS, !is.na(padj), padj < RING_PADJ,
           size >= RING_MIN_SZ, !pathway %in% MITO_DROP_SETS) |>
    arrange(padj)
  n_pre <- nrow(sig_pool)
  if (n_pre > 1)
    sig_pool <- deduplicate_enrichment(as.data.frame(sig_pool), pathways = SET_POOL,
                                       cutoff = SIM_CUT, cross_db = TRUE) |> as_tibble()
  top_terms <- pick_top(sig_pool, RING_N)
  n_path <- nrow(top_terms); n_up <- sum(top_terms$NES > 0); n_dn <- sum(top_terms$NES < 0)

  # Full tested-pathway table for this contrast (every pooled-DB pathway + FDR);
  # shown = drawn on the ring. Filter padj < 0.05 in the workbook for the sig set.
  full_tbl <- fgsea_all |>
    filter(contrast == ctr, database %in% POOL_DBS, !pathway %in% MITO_DROP_SETS) |>
    transmute(pathway, database, padj, pval, NES, size,
              shown = pathway %in% top_terms$pathway) |>
    arrange(padj)

  ring_data <- if (n_path == 0) EMPTY_RING else
    build_ring_180_split(top_terms, ctr, fgsea_all, databases = POOL_DBS)
  if (!is.null(ring_data) && nrow(ring_data) > 0 && "clean_label" %in% names(ring_data))
    ring_data$clean_label <- shorten_label(ring_data$clean_label)
  adaptive_gap <- if (!is.null(ring_data) && nrow(ring_data) > 0)
    0.7 + 0.3 * (max(ring_data$arc_r1_var, na.rm = TRUE) - 4.8) / 1.6 else 0.7

  p <- make_volcano_ring(
    de_df = dep_df, go_df = fgsea_all, contrast = ctr,
    title = NULL, contrast_title = contrast_brief(ctr),
    contrast_subtitle = sprintf("%s | %d DEPs, %d pathways",
                                CONTRAST_MATH_BRIEF[ctr], n_dep, n_path),
    databases = POOL_DBS, ring_data_override = ring_data,
    label_size = 2.7, label_gap = adaptive_gap, title_size = 7, subtitle_size = 4.5,
    point_size = 0.5, point_alpha = 0.55,
    count_label_size = scale_text(BASE_COUNT, 89) + 0.4,
    count_y_mult = 0.75, count_x_mult = 0.85,
    bg_color = unname(CONTRAST_COLORS[ctr]), bg_alpha = 0.20,
    show_legend = FALSE)

  ggsave(file.path(RPT_PDF, sprintf("MAIN_F05_%s_ring.pdf", tag)), p,
         width = 110, height = 110, units = "mm", device = pdf_dev)
  ggsave(file.path(RPT_PNG, sprintf("MAIN_F05_%s_ring.png", tag)), p,
         width = 110, height = 110, units = "mm", dpi = 300)

  shown <- if (n_path > 0)
    top_terms |>
      select(pathway, database, padj, NES, size, any_of(c("ES", "log2err"))) |>
      arrange(desc(NES)) else
    top_terms[FALSE, , drop = FALSE]

  message(sprintf("  [%s/%s] %s: %d pre-dedup -> %d shown (%d up, %d down)",
                  tag, role, contrast_brief(ctr), n_pre, n_path, n_up, n_dn))
  list(tag = tag, ctr = ctr, role = role, n_up = n_up, n_dn = n_dn,
       terms = shown, full = full_tbl)
}

panels <- Map(build_ring, CONTRASTS$ctr, CONTRASTS$tag, CONTRASTS$role)

# Shared NES legend strip (~80 x 20 mm).
nes_legend <- build_nes_legend_bar(text_size = 9, title_size = 10,
                                   bar_margin = margin(2, 6, 2, 6, "mm"))
ggsave(file.path(RPT_PDF, "MAIN_F05_nes_legend.pdf"), nes_legend,
       width = 80, height = 20, units = "mm", device = pdf_dev)
ggsave(file.path(RPT_PNG, "MAIN_F05_nes_legend.png"), nes_legend,
       width = 80, height = 20, units = "mm", dpi = 300)

# shown_pathways.csv — union of all displayed ring terms across the 4 contrasts.
# F06 reads this to avoid re-showing the same biology.
shown_pathways <- bind_rows(lapply(panels, function(pn) {
  if (nrow(pn$terms) == 0) return(NULL)
  pn$terms |> transmute(pathway, database, contrast = pn$ctr, role = pn$role, NES, padj)
}))

# Supplementary workbook: contrast_map + one ring-term sheet per contrast.
contrast_map <- CONTRASTS |>
  mutate(brief = vapply(ctr, contrast_brief, character(1)),
         math  = unname(CONTRAST_MATH_BRIEF[ctr])) |>
  select(old_name = ctr, file_tag = tag, role, brief, definition = math) |>
  as.data.frame()
pathway_sheets <- lapply(panels, function(pn) list(
  name = pn$role, df = as.data.frame(pn$full),
  role     = sprintf("All tested pathways for the %s (%s) contrast", pn$role, contrast_brief(pn$ctr)),
  contents = "pathway, database, padj (BH FDR), pval, NES, size, shown (TRUE = drawn on the ring); sorted by padj — filter padj<0.05 for the significant set"))

build_workbook(
  file.path(DAT, "F05_supplementary.xlsx"),
  figure_title = "F05 — Per-contrast enrichment volcano-in-ring panels (4-DB pooled lens)",
  sheet_specs = c(list(list(name = "contrast_map", df = contrast_map,
    role     = "Key linking file tags / brief names / contrast algebra",
    contents = "old combined-results name, figure file_tag, role, brief display name, and the contrast definition (math)")), pathway_sheets))

# Downstream deliverable: F06 (cluster figure) reads shown_pathways.csv to avoid
# re-showing the same biology — keep it on disk as a loose CSV.
write_csv(shown_pathways, file.path(DAT, "shown_pathways.csv"))

message(sprintf("F05: %d standalone rings + NES legend -> b_reports/main; %d shown pathways logged",
                length(panels), nrow(shown_pathways)))
