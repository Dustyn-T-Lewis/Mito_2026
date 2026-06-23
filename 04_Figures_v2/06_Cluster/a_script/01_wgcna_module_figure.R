#!/usr/bin/env Rscript
# F06 WGCNA module figure — one row-aligned block, module = row:
#   [counts] [hub proteins] [eigengene heatmap, 3 contrasts] [trajectory] [ORA]
# Rows shaded by module colour so clusters read across the panels; modules
# ordered by protein-set size.
#
# Modules: signed WGCNA network (04_Figures/F05_modules; power 7, signed TOM,
# minModuleSize 30, mergeCutHeight 0.25, deepSplit 2). A module eigengene is the
# module's 1st PC — one value per sample — run through the same limma model as
# the proteins (~0 + Group) so heatmap cells are tested log2FC (+ FDR).
# Trajectory = eigengene group means (Ctl -> Mito -> PHE -> PHE_Mito).
# ORA = 8 DBs (GO BP/CC/MF, Hallmark, KEGG, Reactome, MitoCarta, GO Slim),
# per-DB BH at FDR<0.10; best hit per DB, top 4 distinct DBs / module.

library(here)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(limma)
})
source(here::here("04_Figures_v2", "functions", "01_style_palettes_theme.R"))
source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures_v2", "functions", "03a_dedup_engine.R"))
source(here::here("04_Figures_v2", "functions", "03b_enrichment_runners.R"))

BASE <- here::here("04_Figures_v2", "06_Cluster")
MAIN_PDF <- file.path(BASE, "b_reports", "main", "pdf")
MAIN_PNG <- file.path(BASE, "b_reports", "main", "png")
DAT <- file.path(BASE, "c_data")
for (d in c(MAIN_PDF, MAIN_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

DB_COLORS <- c(
  `GO:BP` = "#1B9E77", `GO:CC` = "#66C2A5", `GO:MF` = "#A6D854",
  Hallmark = "#FC8D62", KEGG = "#E78AC3", Reactome = "#8DA0CB",
  MitoCarta = "#FFD92F", `GO Slim` = "#B3B3B3"
)
CONTRASTS <- c(
  Disease    = "PHE - Ctl",
  Transplant = "Mito - Ctl",
  Rescue     = "PHE_Mito - PHE"
)
GROUP_SHORT <- c(Ctl = "Ctl", Mito = "Mito", PHE = "PHE", PHE_Mito = "PHE+Mito")
fmt_fdr <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "<.001",
      ifelse(p < 0.01, sub("^0", "", sprintf("%.3f", p)), sub("^0", "", sprintf("%.2f", p)))
    )
  )
}

# ---- network + eigengene-limma (3 core contrasts) ---------------------------
w <- readRDS(here::here("04_Figures", "F05_modules", "c_data", "wgcna_network.rds"))
meta <- as_tibble(readRDS(P05$imp_rds)$metadata)
non_grey <- setdiff(sub("^ME", "", colnames(w$MEs)), "grey")

me <- as.matrix(w$MEs[, paste0("ME", non_grey), drop = FALSE])
colnames(me) <- non_grey
samp_grp <- factor(meta$Group[match(rownames(me), meta$Col_ID)], levels = H9C2_GROUP_LEVELS)

design <- stats::model.matrix(~ 0 + samp_grp)
colnames(design) <- levels(samp_grp)
cm <- limma::makeContrasts(contrasts = unname(CONTRASTS), levels = design)
colnames(cm) <- names(CONTRASTS)
fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(t(me), design), cm))

mod_stats <- bind_rows(lapply(names(CONTRASTS), function(cn) {
  tt <- limma::topTable(fit, coef = cn, number = Inf, sort.by = "none")
  tibble(module = non_grey, contrast = cn, logFC = tt$logFC, p = tt$P.Value, fdr = tt$adj.P.Val)
}))
# order modules by protein-set size (largest at the top)
mod_size <- count(tibble(module = w$module_colors), module)
mod_order <- mod_size |>
  filter(module %in% non_grey) |>
  arrange(desc(n)) |>
  pull(module)
mod_levels <- rev(mod_order)
n_mod <- length(non_grey)
idx_of <- function(m) as.integer(factor(m, levels = mod_levels))

# ---- consolidated ORA — maximal DB coverage --------------------------------
# 5 canonical DBs (Hallmark/Reactome/KEGG/MitoCarta/GO Slim) + full GO BP/CC/MF.
go_cache <- file.path(DAT, "rat_go_bpccmf_sets.rds")
if (file.exists(go_cache)) {
  go_sets <- readRDS(go_cache)
} else {
  go_sets <- do.call(c, lapply(c("GO:BP", "GO:CC", "GO:MF"), function(sub) {
    m5 <- msigdbr::msigdbr(species = "Rattus norvegicus", collection = "C5", subcollection = sub)
    split(m5$gene_symbol, m5$gs_name)
  }))
  saveRDS(go_sets, go_cache)
}
pw_collection <- c(build_harmonized_collection(), go_sets)

mod_genes <- tibble(gene = w$ann$gene, module = w$module_colors) |>
  filter(!is.na(gene), nzchar(gene), module %in% non_grey)
universe <- unique(w$ann$gene[!is.na(w$ann$gene) & nzchar(w$ann$gene)])
pw_by_db <- split(names(pw_collection), classify_database(names(pw_collection)))

# Per-module ORA: per-DB fora (each DB BH-corrected on its own), then the best
# hit from each DB and the top 3 distinct DBs — so curated pathway DBs (KEGG,
# Reactome, Hallmark, MitoCarta) surface alongside GO rather than being deduped
# into their GO equivalents.
ora_all <- bind_rows(lapply(mod_order, function(m) {
  g <- intersect(mod_genes$gene[mod_genes$module == m], universe)
  res <- bind_rows(lapply(names(pw_by_db), function(db) {
    dp <- pw_collection[pw_by_db[[db]]]
    if (length(dp) < 2) {
      return(NULL)
    }
    r <- as.data.frame(fgsea::fora(pathways = dp, genes = g, universe = universe, minSize = 10, maxSize = 500))
    r$database <- db
    r
  }))
  sig <- res[!is.na(res$padj) & res$padj < 0.10, ]
  if (!nrow(sig)) {
    return(NULL)
  }
  sig |>
    group_by(database) |>
    slice_min(padj, n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(padj) |>
    head(4) |>
    transmute(module = m, database, pathway, padj)
}))

# ---- eigengene group-mean trajectory ----------------------------------------
traj <- as.data.frame(me) |>
  rownames_to_column("Col_ID") |>
  pivot_longer(-Col_ID, names_to = "module", values_to = "v") |>
  mutate(Group = factor(meta$Group[match(Col_ID, meta$Col_ID)], levels = H9C2_GROUP_LEVELS)) |>
  group_by(module, Group) |>
  summarise(m = mean(v), .groups = "drop") |>
  mutate(gx = as.integer(Group), mod_idx = idx_of(module))
traj <- mutate(traj, y = mod_idx + m * (0.4 / max(abs(traj$m))))

# ---- shared row geometry ----------------------------------------------------
row_bg <- tibble(y = seq_len(n_mod), fill = mod_levels)
module_rows <- function(alpha = 0.18) {
  ggplot2::geom_rect(
    data = row_bg,
    ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5, fill = I(fill)),
    alpha = alpha, inherit.aes = FALSE
  )
}
Y_SCALE <- ggplot2::scale_y_continuous(limits = c(0.5, n_mod + 0.5), expand = c(0, 0))
# thin separators between module rows (added to every panel for a clean grid)
ROW_LINES <- ggplot2::geom_hline(yintercept = seq_len(n_mod + 1) - 0.5, color = "grey70", linewidth = 0.2)

# ---- protein-count bars — grow right-to-left, sit against the heatmap --------
cnt_df <- tibble(module = mod_levels, y = seq_len(n_mod)) |>
  left_join(count(tibble(module = w$module_colors), module), by = "module")
p_counts <- ggplot(cnt_df, aes(n, y)) +
  ROW_LINES +
  geom_col(aes(fill = I(module)), width = 0.78, color = "grey40", linewidth = 0.2, orientation = "y") +
  geom_text(aes(label = n), hjust = 1.2, size = 1.6, fontface = "bold", color = "grey20") +
  scale_x_reverse(expand = expansion(mult = c(0.22, 0))) +
  scale_y_continuous(breaks = seq_len(n_mod), labels = mod_levels, limits = c(0.5, n_mod + 0.5), expand = c(0, 0)) +
  labs(title = "Proteins (n)", x = NULL, y = NULL) +
  FIG_THEME +
  theme(
    axis.text.y = element_text(face = "bold", size = FIG_AXIS_TEXT, color = "grey15"),
    axis.text.x = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), panel.border = element_blank(),
    plot.margin = margin(2, 0, 1, 2)
  )

# ---- A. eigengene-limma heatmap (log2FC over FDR; key horizontal below) ------
heat <- mod_stats |>
  mutate(
    contrast = factor(contrast, levels = names(CONTRASTS)),
    y = idx_of(module), x = as.integer(contrast)
  )
ax <- stats::quantile(abs(heat$logFC), 0.98, na.rm = TRUE)
heat$text_col <- ifelse(abs(heat$logFC) >= 0.55 * ax, "white", "grey15")
p_heat <- ggplot(heat, aes(x, y)) +
  geom_tile(aes(fill = logFC), color = "grey85", linewidth = 0.2) +
  ROW_LINES +
  geom_tile(data = filter(heat, fdr < 0.05), fill = NA, color = "black", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", logFC), color = I(text_col)), size = 1.7, fontface = "bold", vjust = -0.35) +
  geom_text(aes(label = fmt_fdr(fdr), color = I(text_col)), size = 1.4, vjust = 1.45) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-ax, ax), oob = scales::squish, name = "eigengene log2FC",
    guide = guide_colorbar(barwidth = 7, barheight = 0.45, title.position = "top", title.hjust = 0.5)
  ) +
  scale_x_continuous(breaks = seq_along(CONTRASTS), labels = names(CONTRASTS), expand = c(0, 0)) +
  Y_SCALE +
  labs(title = "Eigengene response", x = NULL, y = NULL) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = FIG_AXIS_TEXT, face = "bold"),
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), panel.border = element_blank(),
    legend.position = "bottom", legend.title = element_text(size = FIG_LEGEND_TITLE - 0.5),
    legend.margin = margin(t = -2), plot.margin = margin(2, 1, 1, 1)
  )

# ---- B. group trajectory ----------------------------------------------------
p_traj <- ggplot(traj, aes(gx, y, group = module)) +
  module_rows() +
  ROW_LINES +
  geom_segment(data = row_bg, aes(x = 0.7, xend = 4.3, y = y, yend = y), color = "grey85", linewidth = 0.15, inherit.aes = FALSE) +
  geom_line(aes(color = I(module)), linewidth = 0.5) +
  geom_point(aes(fill = I(module)), shape = 21, size = 1, color = "grey30", stroke = 0.2) +
  scale_x_continuous(breaks = 1:4, labels = unname(GROUP_SHORT), limits = c(0.7, 4.3), expand = c(0, 0)) +
  Y_SCALE +
  labs(title = "Trajectory", x = "Group", y = "Eigengene (per-module scaled)") +
  FIG_THEME +
  theme(
    axis.title.x = element_text(size = FIG_AXIS_TEXT, face = "bold", margin = margin(t = 1)),
    axis.title.y = element_text(size = FIG_AXIS_TEXT - 0.5, face = "bold", angle = 90),
    axis.text.x = element_text(angle = 35, hjust = 1, size = FIG_AXIS_TEXT - 0.5),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid = element_blank(), plot.margin = margin(2, 2, 1, 1)
  )

# ---- C. consolidated ORA — bars scaled within each module, top hit first -----
ora_d <- ora_all |>
  mutate(mod_idx = idx_of(module), neglp = -log10(padj)) |>
  group_by(mod_idx) |>
  arrange(padj, .by_group = TRUE) |>
  mutate(
    rank = row_number(), y = mod_idx + (2.5 - rank) * 0.2,
    bar = 0.26 * neglp / max(neglp), # per-module relative length (top hit = full bar)
    lab = clean_pathway_name(pathway),
    lab = ifelse(nchar(lab) > 40, paste0(substr(lab, 1, 38), "…"), lab)
  ) |>
  ungroup()
empty_df <- tibble(module = setdiff(mod_levels, ora_all$module)) |> mutate(y = idx_of(module))
p_ora <- ggplot(ora_d) +
  module_rows() +
  ROW_LINES +
  geom_col(aes(x = bar, y = y, fill = database), orientation = "y", width = 0.16, color = "grey25", linewidth = 0.12) +
  geom_text(aes(0.3, y, label = lab), hjust = 0, vjust = 0.5, size = 1.3, color = "grey10") +
  geom_text(data = empty_df, aes(0.02, y, label = "no enrichment"), hjust = 0, size = 1.4, fontface = "italic", color = "grey55") +
  scale_fill_manual(values = DB_COLORS, name = NULL, guide = guide_legend(nrow = 1)) +
  Y_SCALE +
  scale_x_continuous(limits = c(0, 1.55), expand = c(0, 0)) +
  labs(title = "Pathway enrichment (top hit first)", x = NULL, y = NULL) +
  FIG_THEME +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), legend.position = "bottom",
    legend.key.size = unit(2.2, "mm"), legend.margin = margin(t = -2), plot.margin = margin(2, 1, 1, 1)
  )

# ---- D. hub proteins — top-5 by |kME|, coloured by DE significance -----------
# Each hub's own differential abundance: min FDR across the 3 core contrasts.
hub_sig <- readr::read_csv(P05$comb, show_col_types = FALSE) |>
  filter(contrast %in% c("CTLvPHE", "CTLvMITO", "PHEvPHE_MITO")) |>
  group_by(gene) |>
  summarise(min_fdr = min(adj.P.Val, na.rm = TRUE), .groups = "drop") |>
  mutate(neglfdr = -log10(pmax(min_fdr, 1e-6)))
hub_layout <- w$hubs |>
  filter(module %in% non_grey) |>
  group_by(module) |>
  arrange(desc(abs(kME))) |>
  slice_head(n = 5) |>
  mutate(rank = row_number()) |>
  ungroup() |>
  left_join(hub_sig, by = "gene") |>
  mutate(
    mod_idx = idx_of(module),
    ny = mod_idx + (3 - rank) * 0.18, # rank-1 (top hub) at the top of the band
    face = ifelse(rank == 1, "bold", "plain"),
    neglfdr = ifelse(is.na(neglfdr), 0, neglfdr)
  )
p_hub <- ggplot(hub_layout) +
  module_rows() +
  ROW_LINES +
  geom_point(aes(0.05, ny, fill = neglfdr), shape = 21, size = 1.3, color = "grey40", stroke = 0.2) +
  geom_text(aes(0.12, ny, label = gene, fontface = face), hjust = 0, vjust = 0.5, size = 1.4, color = "grey15") +
  scale_fill_gradient(
    low = "grey90", high = "#B2182B", name = "hub −log10 FDR",
    guide = guide_colorbar(barwidth = 5, barheight = 0.4, title.position = "top", title.hjust = 0.5)
  ) +
  Y_SCALE +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Hub proteins (top-5)", x = NULL, y = NULL) +
  FIG_THEME +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), panel.border = element_blank(),
    legend.position = "bottom", legend.title = element_text(size = FIG_LEGEND_TITLE - 1),
    legend.margin = margin(t = -2), plot.margin = margin(2, 1, 1, 1)
  )

# ---- assemble ---------------------------------------------------------------
key_txt <- paste0(
  "Contrasts:  Disease = PHE−Ctl   |   Transplant = Mito−Ctl   |   Rescue = PHE_Mito−PHE.   ",
  "Cells: eigengene-limma log2FC over FDR; black box = FDR<0.05.\n",
  "ORA: 8 rat DBs (GO BP/CC/MF, Hallmark, KEGG, Reactome, MitoCarta, GO Slim; msigdbr Rattus norvegicus / org.Rn.eg.db), per-DB BH, FDR<0.10; best hit / DB, top 4 / module; bars scaled within module.  Hubs: top-5 |kME|, coloured by own DE significance (min FDR over the 3 contrasts)."
)
fig <- p_counts + p_heat + p_traj + p_ora + p_hub +
  plot_layout(widths = c(0.34, 0.7, 0.5, 1.1, 1.0)) +
  plot_annotation(
    caption = key_txt,
    theme = theme(plot.caption = element_text(size = 4.4, color = "grey35", hjust = 0, lineheight = 1.2))
  )

FIG_W <- PANEL_MD
FIG_H <- 160
ggsave(file.path(MAIN_PDF, "MAIN_F06_wgcna_modules.pdf"), fig,
  width = FIG_W, height = FIG_H, units = "mm", device = pdf_dev, limitsize = FALSE
)
ggsave(file.path(MAIN_PNG, "MAIN_F06_wgcna_modules.png"), fig,
  width = FIG_W, height = FIG_H, units = "mm", dpi = 300, limitsize = FALSE
)
message("F06 WGCNA module figure built")
