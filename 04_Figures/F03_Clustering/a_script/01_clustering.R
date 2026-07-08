#!/usr/bin/env Rscript
# F03 clustering figure. The WGCNA network and module statistics are built upstream in
# 03_DEP/c_modules (00_build_wgcna.R + 01_module_stats.R); this script reads them and
# renders. Each module is a card: its 2x2 contrast effects (Disease/Transplant/Rescue/
# Interaction eigengene shift in SD units, FDR and nominal significance) beside its top
# ORA pathways. Design-responsive modules (nominal eigengene p <= 0.10 on Disease/
# Transplant/Rescue) lead the main figure; the rest, the full top-5 ORA, construction,
# preservation, and hub PPI go to the supplements. ORA enrichment (annotation, not the
# confirmatory test) is computed here because it needs the figures' pathway functions.
# At n=6/group only a few eigengenes clear FDR; WGCNA is the descriptive scaffold and
# confirmatory significance is carried at the protein level (F02).

library(here)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
fns <- here::here("04_Figures", "functions")
source(file.path(fns, "F01-F03_style_palettes_theme.R"))
source(file.path(fns, "F01-F03_data_paths_and_loaders.R"))
source(file.path(fns, "F01-F03_pathway_enrichment_dedup_ora.R"))
source(file.path(fns, "F01-F03_supplementary_workbook.R"))
source(file.path(fns, "F01-F03_composite_layout.R"))
panels <- here::here("04_Figures", "F03_Clustering", "a_script", "panels")
for (f in list.files(panels, full.names = TRUE)) source(f)

BASE <- here::here("04_Figures", "F03_Clustering")
MOD <- here::here("03_DEP", "c_modules", "c_data")
MAIN_PNG <- file.path(BASE, "b_reports", "main", "png")
SUPP_PNG <- file.path(BASE, "b_reports", "supp", "png")
DAT <- file.path(BASE, "c_data")
for (d in c(MAIN_PNG, SUPP_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
P_NOMINAL <- 0.10 # responsive = nominal eigengene p on Disease/Transplant/Rescue
PRIMARY_CONTRASTS <- c("CTLvPHE", "CTLvMITO", "PHEvPHE_MITO")

cache <- file.path(MOD, "wgcna_network.rds")
stats_cache <- file.path(MOD, "module_stats.rds")
if (!file.exists(cache) || !file.exists(stats_cache)) {
  stop("run 03_DEP/c_modules 00_build_wgcna.R then 01_module_stats.R first")
}
w <- readRDS(cache)
ms <- readRDS(stats_cache)
mod_stats <- ms$mod_stats
mod_size <- ms$mod_size
omnibus <- ms$omnibus
mod_order <- arrange(mod_size, desc(n))$module
meta <- as_tibble(readRDS(P05$imp_rds)$metadata)
expr <- as.matrix(readRDS(P05$imp_rds)$data)
non_grey <- setdiff(sub("^ME", "", colnames(w$MEs)), "grey")

# eigengene matrix kept for the group-mean trajectories (the tests live in c_modules)
me <- as.matrix(w$MEs[, paste0("ME", non_grey), drop = FALSE])
colnames(me) <- non_grey

# Design-responsive modules: nominal eigengene p<=0.10 on Disease, Transplant, or Rescue.
# This one set leads the main figure and seeds the hub-PPI/preservation supplements; the
# unresponsive remainder is the "other modules" supplement.
responsive_modules <- mod_stats |>
  filter(contrast %in% c("Disease", "Transplant", "Rescue")) |>
  group_by(module) |>
  summarise(min_p = min(p), .groups = "drop") |>
  filter(min_p <= P_NOMINAL) |>
  arrange(min_p) |>
  pull(module)
if (!length(responsive_modules)) {
  responsive_modules <- head(mod_order, 3)
  message("no module reached eigengene p<=", P_NOMINAL, "; featuring the 3 largest instead")
}
other_modules <- setdiff(non_grey, responsive_modules)

# group-mean eigengene (+/- SE) and the rescue index per module
group_eig <- as_tibble(me, rownames = "Col_ID") |>
  pivot_longer(-Col_ID, names_to = "module", values_to = "eig") |>
  mutate(Group = factor(meta$Group[match(Col_ID, meta$Col_ID)], levels = H9C2_GROUP_LEVELS)) |>
  group_by(module, Group) |>
  summarise(mean_eig = mean(eig), se = sd(eig) / sqrt(n()), .groups = "drop")

# per-protein and eigengene group-mean trajectories (standardized) for the cluster plot
grp_of <- function(ids) factor(meta$Group[match(ids, meta$Col_ID)], levels = H9C2_GROUP_LEVELS)
prot_traj <- as_tibble(t(scale(t(expr))), rownames = "uniprot_id") |>
  pivot_longer(-uniprot_id, names_to = "Col_ID", values_to = "z") |>
  mutate(
    module = w$module_colors[match(uniprot_id, names(w$module_colors))],
    Group = grp_of(Col_ID)
  ) |>
  filter(module %in% non_grey) |>
  group_by(module, uniprot_id, Group) |>
  summarise(zmean = mean(z), .groups = "drop")
eig_traj <- as_tibble(scale(me), rownames = "Col_ID") |>
  pivot_longer(-Col_ID, names_to = "module", values_to = "z") |>
  mutate(Group = grp_of(Col_ID)) |>
  group_by(module, Group) |>
  summarise(zmean = mean(z), .groups = "drop")
fdr_by <- function(cn) {
  setNames(mod_stats$fdr[mod_stats$contrast == cn], mod_stats$module[mod_stats$contrast == cn])
}
disease_fdr <- fdr_by("Disease")
rescue_fdr <- fdr_by("Rescue")

# per-module ORA: per-DB fora, BH within DB, keep every FDR<0.10 hit (full table), then
# EnrichmentMap dedup for the displayed top terms
pw_collection <- c(build_harmonized_collection(), load_go_sets())
pw_by_db <- split(names(pw_collection), classify_database(names(pw_collection)))
mod_genes <- tibble(gene = w$ann$gene, module = w$module_colors) |>
  filter(is_real_symbol(gene), module %in% non_grey)
universe <- unique(w$ann$gene[is_real_symbol(w$ann$gene)])
ora_raw <- bind_rows(lapply(mod_order, function(m) {
  g <- intersect(mod_genes$gene[mod_genes$module == m], universe)
  res <- run_fora_by_db(g, universe, pw_collection, pw_by_db)
  res <- res[!is.na(res$padj), ]
  transmute(arrange(res, padj), module = m, database, pathway, p = pval, padj, n_overlap = overlap, n_set = size)
}))
ora_full <- filter(ora_raw, padj < 0.10)
ora_top5 <- ora_raw |>
  group_by(module) |>
  slice_min(padj, n = 5, with_ties = FALSE) |>
  ungroup()
ora_dedup <- bind_rows(lapply(unique(ora_full$module), function(m) {
  deduplicate_enrichment(ora_full[ora_full$module == m, ], pw_collection) |>
    transmute(module, database, pathway, padj)
}))

# per-protein peak-contrast abundance change for the hub node fill
ab <- readr::read_csv(P05$comb, show_col_types = FALSE) |>
  filter(contrast %in% PRIMARY_CONTRASTS, !is.na(adj.P.Val)) |>
  group_by(uniprot_id) |>
  slice_min(adj.P.Val, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(uniprot_id, ab_logfc = logFC)

module_figure <- function(mods, title) {
  panel_module_card(
    filter(prot_traj, module %in% mods), filter(eig_traj, module %in% mods),
    filter(mod_size, module %in% mods), filter(mod_stats, module %in% mods), ora_top5, title
  )
}

# main figure: each responsive module is a card -- its contrast cells stacked over its
# eigengene trajectory -- beside the top ORA pathway
fig <- module_figure(responsive_modules, sprintf(
  "Design-responsive modules  (min eigengene p ≤ %.2f over Disease/Transplant/Rescue)", P_NOMINAL
))
ggsave(file.path(MAIN_PNG, "MAIN_F03_clustering.png"), fig,
  width = 290, height = 220, units = "mm", dpi = 300, limitsize = FALSE
)

# supplements: unresponsive modules, construction, preservation, hub map
s_other <- module_figure(other_modules, "Modules without an eigengene response")
ggsave(file.path(SUPP_PNG, "SUPP_F03_modules_other.png"), s_other,
  width = 290, height = 205, units = "mm", dpi = 300, limitsize = FALSE
)
s_construction <- (panel_scale_free(w$sft_df, w$chosen_power) / panel_dendro(w$net) +
  plot_layout(heights = c(0.8, 1))) +
  plot_annotation(
    title = "WGCNA network construction",
    theme = theme(plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE, hjust = 0))
  )
ggsave(file.path(SUPP_PNG, "SUPP_F03_construction.png"), s_construction,
  width = PANEL_MD, height = 150, units = "mm", dpi = 300
)
s_preservation <- panel_preservation(w$preservation, responsive_modules)
ggsave(file.path(SUPP_PNG, "SUPP_F03_preservation.png"), s_preservation,
  width = PANEL_MD, height = 120, units = "mm", dpi = 300, limitsize = FALSE
)
s_hub_map <- panel_hub_map(responsive_modules, w, ab)
ggsave(file.path(SUPP_PNG, "SUPP_F03_hub_map.png"), s_hub_map,
  width = PANEL_MD, height = 135, units = "mm", dpi = 300, limitsize = FALSE
)

# supplementary workbook: one sheet per figure component (set tests stay upstream in c_modules)
build_workbook(
  file.path(DAT, "F03_supplementary.xlsx"),
  figure_title = "F03: WGCNA module clustering and the mitochondrial-transplant rescue",
  sheet_specs = list(
    list(
      name = "module_summary",
      df = arrange(mod_size, desc(n)) |>
        left_join(representative_terms(ora_dedup), by = "module") |>
        left_join(transmute(omnibus, module, omnibus_fdr = signif(fdr, 3)), by = "module") |>
        transmute(module,
          n_proteins = n, responsive = module %in% responsive_modules, omnibus_fdr,
          representative_term = term, term_database = as.character(database)
        ),
      role = "Panels A / overview module landscape",
      contents = "module colour, protein count, responsive flag, omnibus-F FDR (any group effect), representative ORA term and its database"
    ),
    list(
      name = "eigengene_limma",
      df = transmute(mod_stats, module, contrast, delta_eig_sd = round(logFC, 3), p = signif(p, 3), fdr = signif(fdr, 3)),
      role = "Panels A / C eigengene response",
      contents = "module x contrast eigengene shift in SD units, p, FDR (limma on eigengenes, paired Replicate block; from 03_DEP/c_modules)"
    ),
    list(
      name = "preservation",
      df = transmute(arrange(w$preservation, condition, desc(Zsummary)),
        module, condition,
        mod_size,
        Zsummary = round(Zsummary, 2), median_rank = round(median_rank, 1)
      ),
      role = "Panel B module preservation",
      contents = "Zsummary and medianRank of each module in each condition vs the control reference"
    ),
    list(
      name = "trajectory",
      df = group_eig |>
        transmute(module,
          group = as.character(Group), eigengene_mean = round(mean_eig, 4),
          se = round(se, 4),
          rescue_fdr = signif(rescue_fdr[module], 3),
          disease_fdr = signif(disease_fdr[module], 3)
        ),
      role = "Panel C rescue trajectories",
      contents = "per-group mean eigengene, SE, Rescue and Disease FDR per module"
    ),
    list(
      name = "module_ora",
      df = ora_full |>
        mutate(shown = paste(module, pathway) %in% paste(ora_dedup$module, ora_dedup$pathway)) |>
        transmute(module, database, pathway, padj = signif(padj, 3), n_overlap, n_set, shown) |>
        arrange(module, padj),
      role = "Supplement S2 pathway landscape",
      contents = "per-module over-representation hits; shown = survived EnrichmentMap dedup"
    ),
    list(
      name = "top_hubs",
      df = bind_rows(lapply(non_grey, function(m) {
        module_hub_nodes(m, w, n_hub = 10L) |>
          left_join(ab, by = "uniprot_id") |>
          transmute(module = m, gene = label, kME = round(kME, 3), peak_log2FC = round(ab_logfc, 3))
      })),
      role = "Panel D hub networks",
      contents = "top-10 proteins per module by kME with peak-contrast log2FC"
    )
  )
)

mirror_to_box(file.path(MAIN_PNG, "MAIN_F03_clustering.png"), "02_Figures/F03_Clustering")
mirror_to_box(list.files(SUPP_PNG, "^SUPP_F03_.*\\.png$", full.names = TRUE), "02_Figures/F03_Clustering/supp")
mirror_to_box(file.path(DAT, "F03_supplementary.xlsx"), "03_Supplementary")

message(sprintf("F03 built | responsive modules: %s", paste(responsive_modules, collapse = ", ")))
print(mod_stats |>
  filter(module %in% responsive_modules, contrast %in% c("Disease", "Transplant", "Rescue")) |>
  arrange(module, contrast))
