#!/usr/bin/env Rscript
# F06 WGCNA module figure — one row-aligned block, module = row:
#   [ cluster colour + name ] [ eigengene-limma heatmap, 5 contrasts ] [ GO BP / CC / MF ]
# Rows are faintly shaded by module colour so the clusters read across the panels.
#
# Modules: signed WGCNA network (04_Figures/F05_modules; power 7, signed TOM,
# minModuleSize 30, mergeCutHeight 0.25, deepSplit 2). A module eigengene is the
# module's 1st PC — one value per sample — run through the same limma model as
# the proteins (~0 + Group, the 5 study contrasts) so cells are tested log2FC.
# GO: clusterProfiler enrichGO + simplify (org.Rn.eg.db, rat), top 3 per ontology.

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
source(here::here("04_Figures_v2", "functions", "09_go_enrichment.R"))

BASE <- here::here("04_Figures_v2", "06_Cluster")
MAIN_PDF <- file.path(BASE, "b_reports", "main", "pdf")
MAIN_PNG <- file.path(BASE, "b_reports", "main", "png")
DAT <- file.path(BASE, "c_data")
for (d in c(MAIN_PDF, MAIN_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

ONT_COLORS <- c(BP = "#1B9E77", CC = "#7570B3", MF = "#D95F02")
CONTRASTS <- c(
  Disease     = "PHE - Ctl",
  Transplant  = "Mito - Ctl",
  Rescue      = "PHE_Mito - PHE",
  Interaction = "(PHE_Mito - Mito) - (PHE - Ctl)",
  Secondary   = "PHE_Mito - Mito"
)

# ---- network + eigengene-limma (tested log2FC + FDR per module x contrast) ---
w <- readRDS(here::here("04_Figures", "F05_modules", "c_data", "wgcna_network.rds"))
meta <- as_tibble(readRDS(P05$imp_rds)$metadata)
non_grey <- setdiff(sub("^ME", "", colnames(w$MEs)), "grey")

me_mat <- t(as.matrix(w$MEs[, paste0("ME", non_grey), drop = FALSE]))
rownames(me_mat) <- non_grey
samp_grp <- factor(meta$Group[match(rownames(w$MEs), meta$Col_ID)], levels = H9C2_GROUP_LEVELS)

design <- stats::model.matrix(~ 0 + samp_grp)
colnames(design) <- levels(samp_grp)
cm <- limma::makeContrasts(contrasts = unname(CONTRASTS), levels = design)
colnames(cm) <- names(CONTRASTS)
fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(me_mat, design), cm))

mod_stats <- bind_rows(lapply(names(CONTRASTS), function(cn) {
  tt <- limma::topTable(fit, coef = cn, number = Inf, sort.by = "none")
  tibble(module = rownames(me_mat), contrast = cn, logFC = tt$logFC, p = tt$P.Value, fdr = tt$adj.P.Val)
}))

mod_order <- mod_stats |>
  group_by(module) |>
  summarise(score = sum(abs(logFC) * -log10(p)), .groups = "drop") |>
  arrange(desc(score)) |>
  pull(module)
mod_levels <- rev(mod_order) # rev so the top-scoring module sits at the top
n_mod <- length(non_grey)

# ---- per-module GO enrichment (cached; simplify() is slow) ------------------
GO_CACHE <- file.path(DAT, "wgcna_go_enrichment.rds")
if (file.exists(GO_CACHE)) {
  go_all <- readRDS(GO_CACHE)
} else {
  mod_genes <- tibble(gene = w$ann$gene, module = w$module_colors) |>
    filter(!is.na(gene), nzchar(gene), module %in% non_grey)
  universe <- unique(w$ann$gene[!is.na(w$ann$gene) & nzchar(w$ann$gene)])
  go_all <- bind_rows(lapply(mod_order, function(m) {
    run_go_module(mod_genes$gene[mod_genes$module == m], universe, top_n = 3) |> mutate(module = m)
  }))
  saveRDS(go_all, GO_CACHE)
}

# ---- shared row geometry: each module is one integer y; faint colour band ----
row_bg <- tibble(y = seq_len(n_mod), fill = mod_levels)
module_rows <- function(alpha = 0.13) {
  ggplot2::geom_rect(
    data = row_bg,
    ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5, fill = I(fill)),
    alpha = alpha, inherit.aes = FALSE
  )
}
Y_SCALE <- ggplot2::scale_y_continuous(limits = c(0.5, n_mod + 0.5), expand = c(0, 0))

# ---- cluster strip: module colour swatch + name -----------------------------
clust_df <- tibble(
  y = seq_len(n_mod), module = mod_levels,
  txt = ifelse(vapply(mod_levels, is_light_color, logical(1)), "grey10", "white")
)
p_clust <- ggplot(clust_df) +
  geom_tile(aes(1, y, fill = I(module)), width = 0.92, height = 0.92, color = "grey40", linewidth = 0.25) +
  geom_text(aes(1, y, label = module, color = I(txt)), size = 1.6, fontface = "bold") +
  Y_SCALE +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  labs(title = "Module", x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE, hjust = 0.5), plot.margin = margin(2, 1, 1, 2))

# ---- A. eigengene-limma contrast heatmap ------------------------------------
heat <- mod_stats |>
  mutate(
    contrast = factor(contrast, levels = names(CONTRASTS)),
    y = as.integer(factor(module, levels = mod_levels)), x = as.integer(contrast),
    stars = sig_stars(fdr)
  )
ax <- stats::quantile(abs(heat$logFC), 0.98, na.rm = TRUE)

p_heat <- ggplot(heat, aes(x, y)) +
  module_rows(alpha = 0.10) +
  geom_tile(aes(fill = logFC), color = "grey85", linewidth = 0.2) +
  geom_tile(data = filter(heat, fdr < 0.05), fill = NA, color = "black", linewidth = 0.5) +
  geom_text(aes(label = stars), size = 2.3, fontface = "bold", color = "grey10", vjust = 0.72) +
  geom_vline(xintercept = 3.5, color = "grey55", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-ax, ax), oob = scales::squish, name = "eigengene\nlog2FC",
    guide = guide_colorbar(barwidth = 0.5, barheight = 3.2)
  ) +
  scale_x_continuous(breaks = seq_along(CONTRASTS), labels = names(CONTRASTS), position = "top", expand = c(0, 0)) +
  Y_SCALE +
  labs(title = "Eigengene response", x = NULL, y = NULL) +
  FIG_THEME +
  theme(
    axis.text.x.top = element_text(angle = 35, hjust = 0, size = FIG_AXIS_TEXT, face = "bold"),
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), panel.border = element_blank(),
    legend.position = "left", plot.margin = margin(2, 1, 1, 1)
  )

# ---- B. GO bar panels (faint module rows + aligned bars) --------------------
build_go_panel <- function(ont, title) {
  d <- go_all |>
    filter(ontology == ont, module %in% non_grey) |>
    mutate(mod_idx = as.integer(factor(module, levels = mod_levels))) |>
    group_by(mod_idx) |>
    mutate(
      rank = row_number(), y = mod_idx + (2 - rank) * 0.27, neglp = -log10(padj),
      lab = ifelse(nchar(term) > 40, paste0(substr(term, 1, 38), "…"), term)
    ) |>
    ungroup()
  x_max <- if (nrow(d)) max(d$neglp) * 1.04 else 1
  ggplot(d) +
    module_rows() +
    geom_col(aes(x = neglp, y = y), orientation = "y", fill = ONT_COLORS[[ont]], width = 0.22, color = "grey30", linewidth = 0.1) +
    geom_text(aes(0.04, y, label = lab), hjust = 0, vjust = 0.5, size = 1.5, color = "grey10") +
    Y_SCALE +
    scale_x_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0.02))) +
    labs(title = title, x = expression(-log[10] ~ FDR), y = NULL) +
    FIG_THEME +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(), plot.margin = margin(2, 2, 1, 1)
    )
}

# ---- assemble (no scatter -> the block is tall) -----------------------------
key_txt <- paste0(
  "Contrasts:  Disease = PHE−Ctl   |   Transplant = Mito−Ctl   |   Rescue = PHE_Mito−PHE",
  "   |   Interaction = (PHE_Mito−Mito)−(PHE−Ctl)   |   Secondary = PHE_Mito−Mito\n",
  "Stars: * FDR<0.05  ** <0.01  *** <0.001 (eigengene limma).  ",
  "GO: clusterProfiler enrichGO + simplify (org.Rn.eg.db, rat); top 3 / ontology / module."
)

fig <- p_clust + p_heat +
  build_go_panel("BP", "GO:BP") + build_go_panel("CC", "GO:CC") + build_go_panel("MF", "GO:MF") +
  plot_layout(widths = c(0.32, 1.15, 1, 1, 1)) +
  plot_annotation(
    caption = key_txt,
    theme = theme(plot.caption = element_text(size = 4.4, color = "grey35", hjust = 0, lineheight = 1.2))
  )

FIG_W <- PANEL_MD
FIG_H <- 135
ggsave(file.path(MAIN_PDF, "MAIN_F06_wgcna_modules.pdf"), fig,
  width = FIG_W, height = FIG_H, units = "mm", device = pdf_dev, limitsize = FALSE
)
ggsave(file.path(MAIN_PNG, "MAIN_F06_wgcna_modules.png"), fig,
  width = FIG_W, height = FIG_H, units = "mm", dpi = 300, limitsize = FALSE
)
message("F06 WGCNA module figure built")
