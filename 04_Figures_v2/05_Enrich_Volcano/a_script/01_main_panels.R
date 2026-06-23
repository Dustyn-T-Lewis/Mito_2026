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
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggforce)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "05_volcano_ring_plot_builder.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))
source(here::here("04_Figures_v2", "functions", "04_mitocarta_lens_lookup.R"))

BASE <- here::here("04_Figures_v2", "05_Enrich_Volcano")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
DAT <- file.path(BASE, "c_data")
for (d in c(RPT_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

dep_df <- load_combined_wide()
fgsea_all <- read_csv(here::here("04_Figures", "shared", "fgsea_tstat_all_h9c2.csv"), show_col_types = FALSE)
rat_gene_sets <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))

# Canonical 5-DB lens (Hallmark + Reactome + KEGG + MitoCarta + GO Slim), shared
# with the pathway-bar and cluster figures via CANONICAL_DBS / MITO_DROP_SETS.
POOL_DBS <- CANONICAL_DBS
SET_POOL <- do.call(c, unname(rat_gene_sets[POOL_DBS]))
RING_N <- 12 # most significant pathways drawn per ring (split Up/Down by NES)
RING_PADJ <- 0.05
SIM_CUT <- 0.375 # EnrichmentMap combined-coefficient redundancy cutoff
RING_MIN_SZ <- 10

# 4 contrasts, OLD name -> filename tag + role brief.
CONTRASTS <- tribble(
  ~ctr,            ~tag,           ~role,
  "CTLvPHE",       "ctlvphe",      "Disease",
  "CTLvMITO",      "ctlvmito",     "Transplant",
  "PHEvPHE_MITO",  "phevphe_mito", "Rescue",
  "Interaction",   "interaction",  "Interaction"
)

# 0-row ring template (correct columns) for panels with no sig terms, so
# make_volcano_ring draws an empty ring instead of its own fallback selection.
EMPTY_RING <- build_ring_180_split(
  head(arrange(filter(fgsea_all, contrast == "PHEvPHE_MITO", database == "Hallmark"), padj), 2),
  "PHEvPHE_MITO", fgsea_all,
  databases = "Hallmark"
)[0, ]

source(here::here("04_Figures_v2", "05_Enrich_Volcano", "a_script", "_build.R"))

panels <- Map(build_ring, CONTRASTS$ctr, CONTRASTS$tag, CONTRASTS$role)

# Shared NES legend strip (~80 x 20 mm).
nes_legend <- build_nes_legend()
ggsave(file.path(RPT_PNG, "MAIN_F05_nes_legend.png"), nes_legend,
  width = 80, height = 20, units = "mm", dpi = 300
)

# shown_pathways.csv — union of all displayed ring terms across the 4 contrasts.
# F06 reads this to avoid re-showing the same biology.
shown_pathways <- bind_rows(lapply(panels, function(pn) {
  if (nrow(pn$terms) == 0) {
    return(NULL)
  }
  pn$terms |> transmute(pathway, database, contrast = pn$ctr, role = pn$role, NES, padj)
}))

# Supplementary workbook: contrast_map + one ring-term sheet per contrast.
contrast_map <- CONTRASTS |>
  mutate(
    brief = vapply(ctr, contrast_brief, character(1)),
    math = unname(CONTRAST_MATH_BRIEF[ctr])
  ) |>
  select(old_name = ctr, file_tag = tag, role, brief, definition = math) |>
  as.data.frame()
pathway_sheets <- lapply(panels, function(pn) {
  list(
    name = pn$role, df = as.data.frame(pn$full),
    role = sprintf("All tested pathways for the %s (%s) contrast", pn$role, contrast_brief(pn$ctr)),
    contents = "pathway, database, padj (BH FDR), pval, NES, size, shown (TRUE = drawn on the ring); sorted by padj — filter padj<0.05 for the significant set"
  )
})

build_workbook(
  file.path(DAT, "F05_supplementary.xlsx"),
  figure_title = "F05 — Per-contrast enrichment volcano-in-ring panels (4-DB pooled lens)",
  sheet_specs = c(list(list(
    name = "contrast_map", df = contrast_map,
    role = "Key linking file tags / brief names / contrast algebra",
    contents = "old combined-results name, figure file_tag, role, brief display name, and the contrast definition (math)"
  )), pathway_sheets)
)

# Downstream deliverable: F06 (cluster figure) reads shown_pathways.csv to avoid
# re-showing the same biology — keep it on disk as a loose CSV.
write_csv(shown_pathways, file.path(DAT, "shown_pathways.csv"))

message(sprintf(
  "F05: %d standalone rings + NES legend -> b_reports/main; %d shown pathways logged",
  length(panels), nrow(shown_pathways)
))
