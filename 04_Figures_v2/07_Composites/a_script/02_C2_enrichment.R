#!/usr/bin/env Rscript
# C2 enrichment — three enrichment rings (Disease / Transplant / Rescue) + NES legend.
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggforce)
  library(patchwork)
})

source(here::here("04_Figures_v2", "functions", "08_composite_layout.R"))
source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "05_volcano_ring_plot_builder.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures_v2", "functions", "04_mitocarta_lens_lookup.R"))

# Set up globals that build_ring() depends on. These mirror 05_Enrich_Volcano/a_script/01_main_panels.R.
# build_ring() has ggsave side-effects to RPT_PDF / RPT_PNG; point them at the
# existing F05 output dirs so the side-effect writes are harmless (re-renders the
# already-existing standalone ring files).
F05_BASE <- here::here("04_Figures_v2", "05_Enrich_Volcano")
RPT_PDF <- file.path(F05_BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(F05_BASE, "b_reports", "main", "png")
for (d in c(RPT_PDF, RPT_PNG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

dep_df <- load_combined_wide()
fgsea_all <- read_csv(
  here::here("04_Figures", "shared", "fgsea_tstat_all_h9c2.csv"),
  show_col_types = FALSE
)
rat_gene_sets <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))

POOL_DBS <- CANONICAL_DBS
SET_POOL <- do.call(c, unname(rat_gene_sets[POOL_DBS]))
RING_N <- 12
RING_PADJ <- 0.05
SIM_CUT <- 0.375
RING_MIN_SZ <- 10

EMPTY_RING <- build_ring_180_split(
  head(arrange(filter(fgsea_all, contrast == "PHEvPHE_MITO", database == "Hallmark"), padj), 2),
  "PHEvPHE_MITO", fgsea_all,
  databases = "Hallmark"
)[0, ]

source(here::here("04_Figures_v2", "05_Enrich_Volcano", "a_script", "_build.R"))

BASE <- here::here("04_Figures_v2", "07_Composites")

disease <- build_ring("CTLvPHE", "ctlvphe", "Disease")$plot
transplant <- build_ring("CTLvMITO", "ctlvmito", "Transplant")$plot
rescue <- build_ring("PHEvPHE_MITO", "phevphe_mito", "Rescue")$plot
nes <- build_nes_legend()

# Layout: three equal rings side by side (A / B / C) filling the full width.
# NES legend inset into the Disease ring (A) at the lower-right corner — that
# corner is consistently clear of ring-segment labels (the left-half labels run
# left of centre; the right-half labels run right of centre, leaving a gap near
# the 5 o'clock position of each ring).  Placing the legend in ring A rather
# than C avoids the dense Rescue label cluster on the right side of panel C.
disease_with_nes <- add_tag(disease, "A") + patchwork::inset_element(
  nes,
  left = 0.62, right = 0.98, top = 0.38, bottom = 0.01,
  align_to = "panel"
)

fig <- disease_with_nes + add_tag(transplant, "B") + add_tag(rescue, "C") +
  plot_layout(ncol = 3) +
  composite_caption(paste(
    "5-DB lens (Hallmark/Reactome/KEGG/MitoCarta/GO Slim), EnrichmentMap dedup.",
    "Ring = top pathways by FDR; centre = protein volcano. NES colour bar shared."
  ))

save_composite(fig, BASE, "MAIN_C2_enrichment", width_mm = PANEL_MD, height_mm = 178)
message("C2 composite built")
