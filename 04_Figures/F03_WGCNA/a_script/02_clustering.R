#!/usr/bin/env Rscript
# F03 WGCNA figure. The network and module statistics are built by 00_build_wgcna.R and
# 01_module_stats.R in this figure's a_script; this script reads their cache and renders.
# The main figure is one aligned card per featured module: its protein count, member-response
# heatmap, group-mean eigengene trajectory, and top ORA pathways. Featured = modules with a
# nominal omnibus moderated-F group effect (p < 0.05); the rest go to the washed-out supplement.
# ORA enrichment (annotation, not the confirmatory test) is computed here because it needs the
# figures' pathway functions. At n=6/group only a few modules clear the gate; WGCNA is the
# descriptive scaffold and confirmatory significance is carried at the protein level (F02).

pacman::p_load(here, dplyr, tidyr, tibble, ggplot2, patchwork)
fns <- here::here("04_Figures", "functions")
source(file.path(fns, "shared_theme_palettes.R"))
source(file.path(fns, "shared_data_loaders.R"))
source(file.path(fns, "shared_enrichment_ora.R"))
source(file.path(fns, "shared_workbook.R"))
panels <- here::here("04_Figures", "F03_WGCNA", "a_script", "panels")
for (f in list.files(panels, full.names = TRUE)) source(f)

BASE <- here::here("04_Figures", "F03_WGCNA")
MOD <- file.path(BASE, "c_data")
RPT <- file.path(BASE, "b_reports")
SUPP_PNG <- file.path(RPT, "supp")
DAT <- MOD
for (d in c(RPT, SUPP_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

P_NOMINAL <- 0.05 # main figure = modules with a nominal omnibus-F group effect
FDR_OMNIBUS <- 0.05 # the subset also clearing FDR is flagged in the workbook
# Design contrasts eligible for the trajectory's headline fry annotation. Interaction and
# DiseaseAfter are excluded: the headline is the module's response to disease or the
# transplant, not a difference-of-differences.
GATE_CONTRASTS <- c("Disease", "Transplant", "Rescue", "Recovery")
HEATMAP_CONTRASTS <- c("Disease", "Transplant", "Rescue", "Recovery", "Interaction")
PRIMARY_CONTRASTS <- c("CTLvPHE", "CTLvMITO", "PHEvPHE_MITO")

cache <- file.path(MOD, "wgcna_network.rds")
stats_cache <- file.path(MOD, "module_stats.rds")
if (!file.exists(cache) || !file.exists(stats_cache)) {
  stop("run 00_build_wgcna.R then 01_module_stats.R (this figure's a_script) first")
}
w <- readRDS(cache)
ms <- readRDS(stats_cache)
mod_stats <- ms$mod_stats
mod_size <- ms$mod_size
omnibus <- ms$omnibus
meta <- as_tibble(readRDS(P05$imp_rds)$metadata)
non_grey <- setdiff(sub("^ME", "", colnames(w$MEs)), "grey")
grey_n <- sum(w$module_colors == "grey")

me <- as.matrix(w$MEs[, paste0("ME", non_grey), drop = FALSE])
colnames(me) <- non_grey

# Main figure: modules with a nominal omnibus group effect (any group difference, F-test
# p < 0.05). The FDR-significant subset (turquoise, blue, brown, greenyellow, tan) is the
# confident core; the two nominal-only modules (red, pink) carry strong ORA and near-miss
# responsiveness, so they earn a card but no starred post-hoc contrast. Every module without
# a nominal effect moves to the washed-out supplement.
size_order <- function(m) m[order(match(m, arrange(mod_size, desc(n))$module))]
responsive_modules <- omnibus |>
  filter(fdr < FDR_OMNIBUS) |>
  pull(module)
featured_modules <- omnibus |>
  filter(p < P_NOMINAL) |>
  pull(module) |>
  size_order()
if (!length(featured_modules)) {
  featured_modules <- head(arrange(mod_size, desc(n))$module, 5)
  message("no module reached nominal significance; featuring the 5 largest instead")
}
other_modules <- setdiff(non_grey, featured_modules)

# Headline member-level gate annotated on each trajectory: the design contrast with the
# smallest fry FDR. fry is the self-contained rotation test (uses every member protein and
# the residual covariance), so it, not the eigengene PC1, decides whether the members move.
fry_head <- ms$settests |>
  filter(contrast %in% GATE_CONTRASTS) |>
  group_by(module) |>
  slice_min(fry_fdr, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(module, contrast, fdr = fry_fdr)

# group-mean eigengene (+/- SE) per module
group_eig <- as_tibble(me, rownames = "Col_ID") |>
  pivot_longer(-Col_ID, names_to = "module", values_to = "eig") |>
  mutate(Group = factor(meta$Group[match(Col_ID, meta$Col_ID)], levels = H9C2_GROUP_LEVELS)) |>
  group_by(module, Group) |>
  summarise(mean_eig = mean(eig), se = sd(eig) / sqrt(n()), .groups = "drop")

fdr_by <- function(cn) {
  setNames(mod_stats$fdr[mod_stats$contrast == cn], mod_stats$module[mod_stats$contrast == cn])
}
disease_fdr <- fdr_by("Disease")
rescue_fdr <- fdr_by("Rescue")

# per-module ORA: per-DB fora, BH within DB, keep every FDR<0.10 hit (full table), then
# EnrichmentMap dedup for the displayed representative terms
pw_collection <- c(build_harmonized_collection(), load_go_sets())
pw_by_db <- split(names(pw_collection), classify_database(names(pw_collection)))
mod_genes <- tibble(gene = w$ann$gene, module = w$module_colors) |>
  filter(is_real_symbol(gene), module %in% non_grey)
universe <- unique(w$ann$gene[is_real_symbol(w$ann$gene)])
mod_order <- arrange(mod_size, desc(n))$module
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

# Each module scored as a gene set by fgsea on each contrast's moderated-t ranking. NES is
# competitive and anti-conservative on modules this large, so it is shown as descriptive
# direction and magnitude only (fry is the gate); the heatmap uses it to make that failure
# mode visible beside the member-level fry.
source(file.path(fns, "shared_gene_set_helpers.R"))
NES_CONTRASTS <- c(
  Disease = "CTLvPHE", Transplant = "CTLvMITO",
  Rescue = "PHEvPHE_MITO", Recovery = "CTLvPHE_MITO", Interaction = "Interaction"
)
dep_wide <- readr::read_csv(P05$comb, show_col_types = FALSE) |>
  filter(is_real_symbol(gene), contrast %in% NES_CONTRASTS) |>
  select(gene, contrast, t) |>
  tidyr::pivot_wider(names_from = contrast, values_from = t, names_prefix = "t_", values_fn = mean)
module_nes_full <- run_fgsea_cache(
  dep_wide, lapply(split(mod_genes$gene, mod_genes$module), unique),
  "WGCNA", unname(NES_CONTRASTS),
  min_size = 10, max_size = 5000, eps = 1e-10
) |>
  rename(module = pathway) |>
  mutate(contrast = names(NES_CONTRASTS)[match(contrast, NES_CONTRASTS)])
module_nes <- select(module_nes_full, module, contrast, NES, padj)

module_card <- function(mods, emphasize = mods) {
  panel_module_card(
    filter(group_eig, module %in% mods), filter(mod_size, module %in% mods),
    ora_top5, ms$settests, module_nes, HEATMAP_CONTRASTS,
    mod_stats = mod_stats, emphasize = emphasize
  )
}

# Main figure: one aligned row per nominally responsive module — its protein-count bar,
# member-response heatmap, eigengene trajectory (with the post-hoc brackets), and top ORA
# pathways.
fig <- module_card(featured_modules) +
  plot_annotation(
    title = "WGCNA co-expression modules and their eigengene response across the design",
    subtitle = sprintf(
      "Rows: modules with a nominal omnibus group effect (F-test p < %.2f). Trajectory brackets = post-hoc contrasts; ✱ = FDR < %.2f. Other modules are supplementary.",
      P_NOMINAL, FDR_OMNIBUS
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, colour = "grey10"),
      plot.subtitle = element_text(face = "italic", size = 9, colour = "grey40")
    )
  )
ggsave(file.path(RPT, "MAIN_F03_wgcna.png"), fig,
  width = 300, height = 48 + 35 * length(featured_modules), units = "mm", dpi = 300, bg = "white", limitsize = FALSE
)

# Supplements: the remaining (non-nominal) modules in the same count/trajectory/ORA
# format, washed out; the full member-level fry grid; and the network construction.
s_other <- module_card(other_modules, emphasize = character(0)) +
  plot_annotation(
    title = "Modules without a nominal omnibus-F group effect (p ≥ 0.05)",
    subtitle = "Same layout as the main figure, washed out to mark them as not responsive: protein count, member-response heatmap, eigengene trajectory, top ORA pathways.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, colour = "grey15"),
      plot.subtitle = element_text(face = "italic", size = 9, colour = "grey40")
    )
  )
ggsave(file.path(SUPP_PNG, "SUPP_F03_modules_other.png"), s_other,
  width = 300, height = 48 + 35 * length(other_modules), units = "mm", dpi = 300, bg = "white", limitsize = FALSE
)
s_construction <- (panel_scale_free(w$sft_df, w$chosen_power) / panel_dendro(w$net) +
  plot_layout(heights = c(0.8, 1))) +
  plot_annotation(
    title = "WGCNA network construction",
    theme = theme(plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE, hjust = 0))
  )
ggsave(file.path(SUPP_PNG, "SUPP_F03_construction.png"), s_construction,
  width = PANEL_MD, height = 150, units = "mm", dpi = 300, bg = "white"
)

# Supplementary workbook: one sheet per figure component. The full module x contrast fry
# grid (S3) surfaces the member-level gate, including the Recovery and Disease-after-
# transplant contrasts the cards do not draw.
build_workbook(
  file.path(DAT, "F03_supplementary.xlsx"),
  figure_title = "F03: WGCNA module clustering and the mitochondrial-transplant rescue",
  sheet_specs = list(
    list(
      name = "module_summary",
      df = arrange(mod_size, desc(n)) |>
        left_join(representative_terms(ora_dedup), by = "module") |>
        left_join(transmute(omnibus, module, omnibus_F = round(.data$F, 2), omnibus_p = signif(p, 3), omnibus_fdr = signif(fdr, 3)), by = "module") |>
        left_join(transmute(fry_head, module, fry_contrast = contrast, fry_fdr = signif(fdr, 3)), by = "module") |>
        transmute(module,
          n_proteins = n,
          in_main_figure = module %in% featured_modules,
          omnibus_fdr_sig = module %in% responsive_modules,
          omnibus_F, omnibus_p, omnibus_fdr, fry_contrast, fry_fdr,
          representative_term = term, term_database = as.character(database)
        ),
      role = "Main figure module rows",
      contents = "module colour, protein count, in-main-figure (nominal omnibus-F p<0.05) and FDR-significant flags, omnibus moderated-F statistic with p and FDR, headline fry contrast and FDR, representative ORA term and its database"
    ),
    list(
      name = "eigengene_limma",
      df = transmute(mod_stats, module, contrast, delta_eig_sd = round(logFC, 3), p = signif(p, 3), fdr = signif(fdr, 3)),
      role = "Eigengene contrasts (main-figure post-hoc brackets)",
      contents = "module x contrast eigengene shift in SD units, nominal p and BH FDR (limma on eigengenes, paired Replicate block; from 01_module_stats.R); the Disease/Transplant/Rescue/Recovery rows are the post-hoc brackets drawn on the main figure"
    ),
    list(
      name = "fry_grid",
      df = ms$settests |>
        transmute(module, contrast,
          direction = as.character(direction),
          fry_p = signif(fry_p, 3), fry_fdr = signif(fry_fdr, 3),
          camera_p = signif(camera_p, 3), camera_fdr = signif(camera_fdr, 3)
        ) |>
        arrange(module, contrast),
      role = "Member-level set tests (the featured gate)",
      contents = "module x contrast self-contained fry and competitive camera results, all six contrasts (Disease, Transplant, Rescue, Recovery, Disease-after, Interaction); fry is the gate, camera is anti-conservative on large modules and reported for completeness"
    ),
    list(
      name = "module_fgsea",
      df = module_nes_full |>
        transmute(module, contrast, NES = round(NES, 3), padj = signif(padj, 3)) |>
        arrange(module, contrast),
      role = "Competitive enrichment (main-figure NES row)",
      contents = "module x contrast fgsea NES and adjusted p over each contrast's moderated-t ranking; competitive and anti-conservative on modules this large (padj collapses toward zero), so read as descriptive direction and magnitude only"
    ),
    list(
      name = "preservation",
      df = transmute(arrange(w$preservation, condition, desc(Zsummary)),
        module, condition,
        mod_size,
        Zsummary = round(Zsummary, 2), median_rank = round(median_rank, 1)
      ),
      role = "Module preservation",
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
      role = "Eigengene trajectories",
      contents = "per-group mean eigengene, SE, Rescue and Disease FDR per module"
    ),
    list(
      name = "module_ora",
      df = ora_full |>
        mutate(shown = paste(module, pathway) %in% paste(ora_dedup$module, ora_dedup$pathway)) |>
        transmute(module, database, pathway, padj = signif(padj, 3), n_overlap, n_set, shown) |>
        arrange(module, padj),
      role = "Pathway landscape",
      contents = "per-module over-representation hits; shown = survived EnrichmentMap dedup"
    ),
    list(
      name = "top_hubs",
      df = bind_rows(lapply(non_grey, function(m) {
        module_hub_nodes(m, w, n_hub = 10L) |>
          left_join(ab, by = "uniprot_id") |>
          transmute(module = m, gene = label, kME = round(kME, 3), peak_log2FC = round(ab_logfc, 3))
      })),
      role = "Hub networks",
      contents = "top-10 proteins per module by kME with peak-contrast log2FC"
    )
  )
)

message(sprintf(
  "F03 built | main (omnibus p<%.2f): %s | other: %d modules",
  P_NOMINAL, paste(featured_modules, collapse = ", "), length(other_modules)
))

# Optional local Box mirror (author's machine only; no-ops elsewhere).
mirror_to_box(file.path(RPT, "MAIN_F03_wgcna.png"), "02_Figures/F03_WGCNA", guard = TRUE)
mirror_to_box(list.files(SUPP_PNG, "^SUPP_F03_.*\\.png$", full.names = TRUE), "02_Figures/F03_WGCNA/supp", guard = TRUE)
mirror_to_box(file.path(DAT, "F03_supplementary.xlsx"), "03_Supplementary", guard = TRUE)
