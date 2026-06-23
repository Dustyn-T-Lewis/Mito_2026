#!/usr/bin/env Rscript
# F06 WGCNA module figure.
#   A  module eigengene response to the 5 contrasts (limma log2FC + FDR stars)
#   B  GO:BP / CC / MF over-representation per module (enrichGO + simplify)
#   C  Disease <-> Rescue reversal scatter (which modules transplant corrects)
#
# Modules: signed WGCNA network (04_Figures/F05_modules; power 7, signed TOM,
# minModuleSize 30, mergeCutHeight 0.25, deepSplit 2). A module eigengene is the
# module's 1st PC — one value per sample. Each eigengene is run through the same
# limma model as the proteins (~0 + Group, the 5 study contrasts) so the cell
# values are tested log2 fold-changes, not descriptive correlations.

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
HEAT_LO <- "#2166AC"
HEAT_HI <- "#B2182B"

# Contrast names + algebra (used for columns, the bottom key, and the scatter).
CONTRASTS <- c(
  Disease     = "PHE - Ctl",
  Transplant  = "Mito - Ctl",
  Rescue      = "PHE_Mito - PHE",
  Interaction = "(PHE_Mito - Mito) - (PHE - Ctl)",
  Secondary   = "PHE_Mito - Mito"
)

# ---- network + sample metadata ----------------------------------------------
w <- readRDS(here::here("04_Figures", "F05_modules", "c_data", "wgcna_network.rds"))
meta <- as_tibble(readRDS(P05$imp_rds)$metadata)
non_grey <- setdiff(sub("^ME", "", colnames(w$MEs)), "grey")

# ---- eigengene-limma: one tested log2FC + FDR per module x contrast ----------
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

# module order: total signal across contrasts (most responsive at the top)
mod_order <- mod_stats |>
  group_by(module) |>
  summarise(score = sum(abs(logFC) * -log10(p)), .groups = "drop") |>
  arrange(desc(score)) |>
  pull(module)
mod_levels <- rev(mod_order) # rev so the top-scoring module sits at the top

# ---- per-module GO enrichment (cached; simplify() is slow) ------------------
GO_CACHE <- file.path(DAT, "wgcna_go_enrichment.rds")
if (file.exists(GO_CACHE)) {
  go_all <- readRDS(GO_CACHE)
  message("GO enrichment loaded from cache")
} else {
  message("running GO enrichment (enrichGO + simplify per module) ...")
  mod_genes <- tibble(gene = w$ann$gene, module = w$module_colors) |>
    filter(!is.na(gene), nzchar(gene), module %in% non_grey)
  universe <- unique(w$ann$gene[!is.na(w$ann$gene) & nzchar(w$ann$gene)])
  go_all <- bind_rows(lapply(mod_order, function(m) {
    run_go_module(mod_genes$gene[mod_genes$module == m], universe, top_n = 3) |>
      mutate(module = m)
  }))
  saveRDS(go_all, GO_CACHE)
}

n_mod <- length(non_grey)

# ---- A. eigengene-limma contrast heatmap ------------------------------------
heat <- mod_stats |>
  mutate(
    module = factor(module, levels = mod_levels),
    contrast = factor(contrast, levels = names(CONTRASTS)),
    y = as.integer(module), x = as.integer(contrast),
    stars = sig_stars(fdr)
  )
ax <- stats::quantile(abs(heat$logFC), 0.98, na.rm = TRUE)

p_heat <- ggplot(heat, aes(x, y)) +
  geom_tile(aes(fill = logFC), color = "grey85", linewidth = 0.2) +
  geom_tile(data = filter(heat, fdr < 0.05), fill = NA, color = "black", linewidth = 0.5) +
  geom_text(aes(label = stars), size = 2.3, fontface = "bold", color = "grey10", vjust = 0.72) +
  geom_vline(xintercept = 3.5, color = "grey55", linewidth = 0.3) +
  scale_fill_gradient2(
    low = HEAT_LO, mid = "white", high = HEAT_HI, midpoint = 0,
    limits = c(-ax, ax), oob = scales::squish, name = "eigengene\nlog2FC",
    guide = guide_colorbar(barwidth = 0.5, barheight = 3.2)
  ) +
  scale_x_continuous(
    breaks = seq_along(CONTRASTS), labels = names(CONTRASTS),
    position = "top", expand = c(0, 0)
  ) +
  scale_y_continuous(breaks = seq_len(n_mod), labels = mod_levels, expand = c(0, 0)) +
  labs(title = "Module eigengene response", x = NULL, y = NULL) +
  FIG_THEME +
  theme(
    axis.text.x.top = element_text(angle = 35, hjust = 0, size = FIG_AXIS_TEXT, face = "bold"),
    axis.text.y = element_text(size = FIG_AXIS_TEXT, face = "bold", color = "grey15"),
    axis.ticks = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
    legend.position = "left", plot.margin = margin(2, 1, 1, 2)
  )

# ---- B. GO bar panels (shared module y -> rows align with the heatmap) ------
build_go_panel <- function(ont, title) {
  d <- go_all |>
    filter(ontology == ont, module %in% non_grey) |>
    mutate(module = factor(module, levels = mod_levels), mod_idx = as.integer(module)) |>
    group_by(module) |>
    mutate(
      rank = row_number(), y = mod_idx + (2 - rank) * 0.28, neglp = -log10(padj),
      lab = ifelse(nchar(term) > 40, paste0(substr(term, 1, 38), "…"), term)
    ) |>
    ungroup()
  x_max <- if (nrow(d)) max(d$neglp) * 1.04 else 1
  ggplot(d) +
    geom_col(aes(x = neglp, y = y),
      orientation = "y", fill = ONT_COLORS[[ont]],
      width = 0.24, color = "grey30", linewidth = 0.1
    ) +
    geom_text(aes(0.04, y, label = lab), hjust = 0, vjust = 0.5, size = 1.4, color = "grey10") +
    scale_y_continuous(limits = c(0.5, n_mod + 0.5), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0.02))) +
    labs(title = title, x = expression(-log[10] ~ FDR), y = NULL) +
    FIG_THEME +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(), plot.margin = margin(2, 2, 1, 1)
    )
}
p_bp <- build_go_panel("BP", "GO:BP")
p_cc <- build_go_panel("CC", "GO:CC")
p_mf <- build_go_panel("MF", "GO:MF")

# ---- C. Disease <-> Rescue reversal scatter ---------------------------------
rev_w <- mod_stats |>
  filter(contrast %in% c("Disease", "Rescue")) |>
  select(module, contrast, logFC, fdr) |>
  pivot_wider(names_from = contrast, values_from = c(logFC, fdr)) |>
  mutate(
    reversed = sign(logFC_Disease) != sign(logFC_Rescue) &
      pmin(fdr_Disease, fdr_Rescue) < 0.05
  )
lim <- max(abs(c(rev_w$logFC_Disease, rev_w$logFC_Rescue)), na.rm = TRUE) * 1.12
mod_cols <- stats::setNames(non_grey, non_grey)

p_rev <- ggplot(rev_w, aes(logFC_Disease, logFC_Rescue)) +
  annotate("rect", xmin = 0, xmax = lim, ymin = -lim, ymax = 0, fill = "#4DAF4A", alpha = 0.07) +
  annotate("rect", xmin = -lim, xmax = 0, ymin = 0, ymax = lim, fill = "#4DAF4A", alpha = 0.07) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
  geom_point(aes(fill = I(mod_cols[as.character(module)])),
    shape = 21, size = 2.4, color = "grey25", stroke = 0.3
  ) +
  geom_text(aes(label = module),
    size = 1.7, fontface = "bold", color = "grey15", vjust = -0.9
  ) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(
    title = "Disease ↔ Rescue reversal",
    subtitle = "green quadrants = transplant reverses the disease shift",
    x = "Disease log2FC (PHE − Ctl)",
    y = "Rescue log2FC (PHE_Mito − PHE)"
  ) +
  FIG_THEME +
  theme(plot.margin = margin(2, 2, 1, 2))

# ---- assemble ---------------------------------------------------------------
key_txt <- paste0(
  "Contrasts:  Disease = PHE−Ctl   |   Transplant = Mito−Ctl   |   Rescue = PHE_Mito−PHE",
  "   |   Interaction = (PHE_Mito−Mito)−(PHE−Ctl)   |   Secondary = PHE_Mito−Mito\n",
  "Stars: * FDR<0.05  ** <0.01  *** <0.001 (eigengene limma).  ",
  "GO: clusterProfiler enrichGO + simplify (org.Rn.eg.db, rat); top 3 / ontology / module."
)

top_row <- p_heat + p_bp + p_cc + p_mf + plot_layout(widths = c(1.15, 1, 1, 1))
fig <- (top_row / p_rev) +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(
    caption = key_txt,
    theme = theme(plot.caption = element_text(size = 4.4, color = "grey35", hjust = 0, lineheight = 1.2))
  )

FIG_W <- PANEL_MD
FIG_H <- 165
ggsave(file.path(MAIN_PDF, "MAIN_F06_wgcna_modules.pdf"), fig,
  width = FIG_W, height = FIG_H, units = "mm", device = pdf_dev, limitsize = FALSE
)
ggsave(file.path(MAIN_PNG, "MAIN_F06_wgcna_modules.png"), fig,
  width = FIG_W, height = FIG_H, units = "mm", dpi = 300, limitsize = FALSE
)
message("F06 WGCNA module figure built")
