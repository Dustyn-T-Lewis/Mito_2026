#!/usr/bin/env Rscript
# F05 SUPPLEMENT — WGCNA composite (exploratory de-novo cross-check).
# The F05 headline is the knowledge-driven set-score modules (01_set_scores.R +
# 02_main_panels.R); this WGCNA view is supplementary because its modules are
# data-derived from the same 24 samples then tested on them (double-dipping) and
# N=24 is too small for module preservation. Kept as the unsupervised counterpart.
# Layout (YvO-F06 A + B, ported into 05_Figures):
# A: module × 5-contrast eigengene-limma heatmap
#    rows: modules (ordered by total signal score across the 5 contrasts)
#    cols: CTLvPHE · CTLvMITO · PHEvPHE_MITO · Interaction · MITOvPHE_MITO
#    cell: eigengene-limma logFC; * FDR<0.05  ** FDR<0.01  *** FDR<0.001
#    left annotation: module color block + protein count
#    right annotation: top-3 hub genes (|kME|)
# B: pathway-NES scatter diptych (reversal + concordance) — points coloured by module
#    (a) Disease vs Intervention  (b) Intervention vs Rescue  (c) Disease vs Rescue
# Soft-threshold + eigengene boxplots moved to supp.
# Run 05_Figures/shared/build_wgcna_network.R first to build the signed Pearson
# network (writes c_data/wgcna_network.rds).

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(ggrepel); library(patchwork); library(cowplot); library(grid)
  library(WGCNA); library(limma); library(openxlsx); library(scales)
  library(readr); library(png); library(ComplexHeatmap); library(circlize)
})

source(here::here("04_Figures", "shared", "style.R"))
source(here::here("04_Figures", "shared", "figure_supplement_helpers.R"))

BASE    <- here::here("05_Figures", "F05_modules")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png")
PNL_PNG <- file.path(RPT_PNG, "panels")
PNL_PDF <- file.path(RPT_PDF, "panels")
DAT     <- file.path(BASE, "c_data")
for (d in c(PNL_PNG, PNL_PDF, RPT_PDF, DAT))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

w        <- readRDS(file.path(DAT, "wgcna_network.rds"))
meta     <- read.csv(here::here("00_input", "H9c2_meta.csv"))
meta$Group <- factor(meta$Group, levels = H9C2_GROUP_LEVELS)
meta     <- meta[match(w$sample_order, meta$Col_ID), ]
fgsea_all <- read_csv(here::here("04_Figures", "shared",
                                 "fgsea_tstat_all_h9c2.csv"),
                      show_col_types = FALSE)
sets <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))

# Eigengene-limma against all 5 contrasts

ME_mat <- t(as.matrix(w$MEs))
colnames(ME_mat) <- meta$Col_ID
design <- model.matrix(~ 0 + meta$Group)
colnames(design) <- gsub("^meta\\$Group", "", colnames(design))

cv <- h9c2_parse_contrasts()
cm <- makeContrasts(contrasts = cv, levels = design)
colnames(cm) <- names(cv)
fit  <- lmFit(ME_mat, design)
fit2 <- eBayes(contrasts.fit(fit, cm))

mod_results <- bind_rows(lapply(names(cv), \(cname) {
  tt <- topTable(fit2, coef = cname, number = Inf, sort.by = "none")
  tibble(module    = sub("^ME", "", rownames(ME_mat)),
         contrast  = cname,
         logFC     = tt$logFC,
         P.Value   = tt$P.Value,
         adj.P.Val = tt$adj.P.Val)
}))

# Module ordering: total signal score across the 5 contrasts (descending)
mod_levels <- mod_results |>
  filter(module != "grey") |>
  group_by(module) |>
  summarise(score = sum(abs(logFC) * -log10(P.Value))) |>
  arrange(desc(score)) |>
  pull(module)

mod_sizes <- as.integer(table(w$module_colors)[mod_levels])
names(mod_sizes) <- mod_levels

hub_top3 <- w$hubs |>
  filter(module %in% mod_levels) |>
  group_by(module) |>
  slice_head(n = 3) |>
  summarise(hubs = paste(gene, collapse = ", "), .groups = "drop") |>
  mutate(module = factor(module, levels = mod_levels)) |>
  arrange(module)

# Panel A — Module × contrast tile heatmap (YvO _panel_A_module_heatmap.R port)
# Custom ggplot layout (not ComplexHeatmap): bracket -> count-bars | tiles | hubs -> legend.

library(ggnewscale)

contrasts_5 <- c("CTLvPHE", "CTLvMITO", "PHEvPHE_MITO", "Interaction", "MITOvPHE_MITO")
contrast_col_labels <- c(
  CTLvPHE       = "Disease\n(PHE − Ctl)",
  CTLvMITO      = "Intervention\n(Mito − Ctl)",
  PHEvPHE_MITO  = "Rescue\n(PHE_Mito − PHE)",
  Interaction   = "Interaction\n(Rescue − Intervention)",
  MITOvPHE_MITO = "Secondary\n(PHE_Mito − Mito)")

logfc_mat <- mod_results |>
  filter(module %in% mod_levels) |>
  select(module, contrast, logFC) |>
  pivot_wider(names_from = contrast, values_from = logFC) |>
  arrange(match(module, mod_levels)) |>
  column_to_rownames("module") |> as.matrix()
logfc_mat <- logfc_mat[, contrasts_5]

fdr_mat <- mod_results |>
  filter(module %in% mod_levels) |>
  select(module, contrast, adj.P.Val) |>
  pivot_wider(names_from = contrast, values_from = adj.P.Val) |>
  arrange(match(module, mod_levels)) |>
  column_to_rownames("module") |> as.matrix()
fdr_mat <- fdr_mat[, contrasts_5]

star_mat <- sig_stars(as.vector(fdr_mat))
dim(star_mat) <- dim(fdr_mat); dimnames(star_mat) <- dimnames(fdr_mat)

# Hub-derived biology theme per module (mirrors YvO's mod_bio_labels.csv).
mod_bio_labels <- c(
  yellow    = "Proteostasis /\nGlycolysis",
  turquoise = "ER Stress /\nCytoskeleton",
  green     = "RNA Processing /\nRibosome",
  red       = "Chromatin /\nNuclear",
  blue      = "Mitochondrial /\nTrafficking",
  black     = "ECM / Smooth\nMuscle",
  brown     = "Stress Response /\nLipid")

# Module display IDs (M1 = largest module by signal score).
mod_ids <- setNames(paste0("M", seq_along(mod_levels)), mod_levels)

# Order rows bottom-to-top by protein count (matches YvO `mod_order`).
mod_order <- names(sort(mod_sizes))

hub_text <- setNames(hub_top3$hubs, as.character(hub_top3$module))[mod_order]
hub_text[is.na(hub_text)] <- ""

# Column positions (gap between core contrasts and secondary)
xpos      <- c(1, 2, 3, 4, 5.5)
trait_x   <- setNames(xpos, contrasts_5)
xmin_all  <- min(xpos) - 0.55
xmax_all  <- max(xpos) + 0.55

# heat_df: one row per (module, contrast) cell
heat_df <- expand.grid(module = mod_order, contrast = contrasts_5,
                       stringsAsFactors = FALSE) |>
  mutate(logFC = as.vector(logfc_mat[mod_order, ]),
         padj  = as.vector(fdr_mat[mod_order, ]),
         stars = sig_stars(padj),
         label = sprintf("%.2f", logFC),
         xpos  = trait_x[contrast],
         module = factor(module, levels = mod_order))

# Sizing — matches YvO _panel_A_module_heatmap.R sizing block.
PA_W <- 400; PA_H <- 280
txt_cell  <- scale_text(BASE_GENE, PA_W) * 0.95 * 1.25 + 0.7
txt_count <- scale_text(BASE_COUNT, PA_W) * 0.85 * 1.25 + 2.0
txt_brack <- scale_text(BASE_GENE, PA_W) * 0.95 * 1.25

ax_max <- max(abs(logfc_mat), na.rm = TRUE)

# Brackets row (top section label)
brackets <- tibble(
  label = sprintf("Module-Level Contrasts (eigengene limma, BH)\n%d modules | n = %d wells",
                  length(mod_order), ncol(ME_mat)),
  start = xpos[1], end = xpos[4], mid = (xpos[1] + xpos[4]) / 2)
brackets2 <- tibble(label = "Secondary", start = xpos[5], end = xpos[5],
                    mid = xpos[5])

p_brackets <- ggplot() +
  geom_segment(data = brackets,
               aes(x = start - 0.4, xend = end + 0.4, y = 0.20, yend = 0.20),
               linewidth = 0.4, color = "grey30") +
  geom_segment(data = brackets,
               aes(x = start - 0.4, xend = start - 0.4, y = 0.20, yend = 0.10),
               linewidth = 0.4, color = "grey30") +
  geom_segment(data = brackets,
               aes(x = end + 0.4, xend = end + 0.4, y = 0.20, yend = 0.10),
               linewidth = 0.4, color = "grey30") +
  geom_text(data = brackets, aes(x = mid, y = 0.60, label = label),
            size = txt_brack, fontface = "bold", color = "grey25", lineheight = 0.85) +
  geom_segment(data = brackets2,
               aes(x = start - 0.4, xend = end + 0.4, y = 0.20, yend = 0.20),
               linewidth = 0.4, color = "grey30") +
  geom_segment(data = brackets2,
               aes(x = start - 0.4, xend = start - 0.4, y = 0.20, yend = 0.10),
               linewidth = 0.4, color = "grey30") +
  geom_segment(data = brackets2,
               aes(x = end + 0.4, xend = end + 0.4, y = 0.20, yend = 0.10),
               linewidth = 0.4, color = "grey30") +
  geom_text(data = brackets2, aes(x = mid, y = 0.60, label = label),
            size = txt_brack, fontface = "bold", color = "grey25", lineheight = 0.85) +
  scale_x_continuous(limits = c(xmin_all, xmax_all), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(2, 2, -12, 0))

# Count bars (left panel: protein count per module, bio_label on bar)
mod_counts <- tibble(
  module     = factor(mod_order, levels = mod_order),
  mod_color  = mod_order,
  n_proteins = mod_sizes[mod_order],
  bio_label  = stringr::str_wrap(
    gsub("\n", " ", mod_bio_labels[mod_order]), width = 14),
  module_id  = mod_ids[mod_order],
  bar_text_col = ifelse(vapply(mod_order, is_light_color, logical(1)),
                        "black", "white"))

p_counts <- ggplot(mod_counts, aes(x = sqrt(n_proteins), y = module)) +
  geom_col(fill = mod_counts$mod_color, color = "black",
           linewidth = 0.3, width = 0.65) +
  geom_text(aes(label = bio_label, x = sqrt(n_proteins) * 0.95),
            size = txt_count * 0.78, fontface = "bold",
            color = mod_counts$bar_text_col, hjust = 0, lineheight = 0.85) +
  geom_text(aes(label = sprintf("%s (n=%d)", module_id, n_proteins),
                x = sqrt(n_proteins) / 2, y = module),
            nudge_y = 0.42, size = txt_count * 0.7, fontface = "bold",
            color = "black", hjust = 0.5) +
  scale_x_reverse(expand = c(0, 0),
                  limits = c(sqrt(max(mod_counts$n_proteins)) * 1.10, 0),
                  breaks = sqrt(c(100, 300, 600, 1000)),
                  labels = c(100, 300, 600, 1000)) +
  scale_y_discrete(labels = NULL) +
  labs(y = NULL, x = "Protein Counts") +
  FIG_THEME +
  theme(axis.text.y       = element_blank(),
        axis.ticks.y      = element_blank(),
        axis.text.x       = element_text(size = txt_cell * 1.7 + 1, face = "bold"),
        axis.ticks.x      = element_line(color = "black", linewidth = 0.6),
        axis.ticks.length.x = unit(1.5, "mm"),
        axis.title.x      = element_text(size = txt_cell * 2.3 + 2, face = "bold",
                                          margin = margin(t = -9, unit = "mm")),
        panel.grid.major.y = element_blank(),
        panel.grid.minor  = element_blank(),
        panel.border      = element_blank(),
        axis.line.x       = element_line(color = "black", linewidth = 0.3),
        legend.position   = "none",
        plot.margin       = margin(2, 0, -6, 2, "mm"))

# Heatmap (centre): tile + cell value + stars + section divider
heat_df <- heat_df |>
  mutate(has_sig  = !is.na(padj) & padj < 0.05,
         is_large = !is.na(logFC) & abs(logFC) >= 0.6 * ax_max,
         text_col = ifelse(is_large, "white", "grey10"))

p_heat <- ggplot(heat_df, aes(x = xpos, y = module)) +
  geom_tile(aes(fill = logFC), color = "black", linewidth = 0.3) +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-ax_max, ax_max),
                       oob = scales::squish, name = "Eigengene logFC",
                       na.value = "grey90", guide = "none") +
  geom_tile(data = heat_df |> filter(has_sig),
            color = "black", linewidth = 1.0, fill = NA) +
  geom_text(aes(label = label),
            size = txt_cell - 0.5, fontface = "bold",
            color = heat_df$text_col,
            nudge_y = -0.10) +
  geom_text(data = heat_df |> filter(stars != ""),
            aes(label = stars), size = (txt_cell - 0.5) * 1.15,
            fontface = "bold", color = "grey10",
            nudge_y = 0.30) +
  geom_vline(xintercept = (xpos[4] + xpos[5]) / 2, color = "grey55",
             linewidth = 0.35) +
  scale_x_continuous(breaks = xpos, labels = contrast_col_labels[contrasts_5],
                     limits = c(xmin_all, xmax_all), expand = c(0, 0)) +
  scale_y_discrete(labels = NULL) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.x       = element_text(angle = 0, hjust = 0.5,
                                          size = txt_cell * 2.0 + 2,
                                          face = "bold", lineheight = 0.85),
        axis.text.y       = element_blank(),
        axis.ticks        = element_blank(),
        panel.grid        = element_blank(),
        panel.border      = element_blank(),
        legend.position   = "none",
        plot.margin       = margin(0, 2, 0, -3))

# Hub gene column (right)
hub_df <- tibble(module    = factor(mod_order, levels = mod_order),
                 hubs      = hub_text)
p_hubs <- ggplot(hub_df, aes(y = module)) +
  geom_text(aes(x = 0, label = hubs), hjust = 0, size = txt_count * 0.85,
            fontface = "italic", color = "grey20", lineheight = 0.9) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_discrete(labels = NULL) +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0, "mm"))

# Legend strip (bottom): eigengene-logFC color bar
df_leg <- data.frame(x = c(0, 0), y = c(0, 0), v = c(-ax_max, ax_max))
p_leg <- ggplot(df_leg, aes(x, y, fill = v)) +
  geom_tile(alpha = 0) +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-ax_max, ax_max),
                       oob = scales::squish, name = "Eigengene logFC",
                       guide = guide_colorbar(barwidth  = unit(52, "mm"),
                                              barheight = unit(4.5, "mm"),
                                              title.position = "left",
                                              title.vjust    = 0.8)) +
  theme_void() +
  theme(legend.position    = "top",
        legend.title       = element_text(size = txt_cell * 2.2,
                                          face = "bold.italic"),
        legend.text        = element_text(size = txt_cell * 1.8),
        legend.margin      = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0, "mm"),
        plot.margin        = margin(-48, 0, 0, 0, unit = "mm"))

# Assemble Panel A (patchwork area)
design <- c(
  area(1,  5,  1, 13),                   # brackets across heat columns
  area(2,  1, 10,  4),                   # counts (4 cols wide)
  area(2,  5, 10, 11),                   # heatmap (7 cols wide)
  area(2, 12, 10, 13),                   # hub-gene column
  area(11, 5, 11, 13))                   # legend strip
fig_A <- wrap_elements(p_brackets) + p_counts + p_heat + p_hubs + p_leg +
  plot_layout(design = design,
              heights = c(0.45, rep(1, 9), 0.30)) +
  plot_annotation(theme = theme(plot.margin = margin(12, 2, 2, 2)))

ht_png <- file.path(PNL_PNG, "MAIN_panel_A_module_heatmap.png")
ggsave(ht_png, fig_A, width = PA_W, height = PA_H, units = "mm",
       dpi = 300, limitsize = FALSE)
ggsave(file.path(PNL_PDF, "MAIN_panel_A_module_heatmap.pdf"),
       fig_A, width = PA_W, height = PA_H, units = "mm",
       device = get_pdf_device(), limitsize = FALSE)

# Cell-value table for supp xlsx (replaces ComplexHeatmap's mod_results sheet).
write_csv(heat_df |>
            select(module, contrast, logFC, padj, stars),
          file.path(DAT, "01_panel_A_heatmap_data.csv"))

pA <- wrap_elements(full = rasterGrob(readPNG(ht_png), interpolate = TRUE)) +
  theme(plot.margin = margin(2, 4, 2, 4, "mm"))

# Panel B — Module-NES scatter diptych (reversal + concordance)
# YvO F06 pattern: treat each WGCNA module as an fGSEA "pathway" (its members),
# rank proteins by t-stat per contrast, score each module via fgseaMultilevel.
# Result: ONE point per module per scatter, fill = module color, size = n proteins.

library(fgsea)

mod_genes_map <- split(w$ann$gene[match(w$gene_order, w$ann$uniprot_id)],
                       w$module_colors) |>
  lapply(\(g) unique(g[!is.na(g) & g != ""]))
mod_genes_map <- mod_genes_map[setdiff(names(mod_genes_map), "grey")]

dep <- read_csv(here::here("03_DEP", "c_data", "03_combined_results.csv"),
                show_col_types = FALSE)

build_ranks <- function(t_col) {
  vals <- dep[[t_col]]; names(vals) <- dep$gene
  vals <- vals[!is.na(vals) & names(vals) != "" & !is.na(names(vals))]
  sort(vals[!duplicated(names(vals))], decreasing = TRUE)
}
run_module_fgsea <- function(t_col) {
  set.seed(42)
  res <- fgseaMultilevel(pathways = mod_genes_map, stats = build_ranks(t_col),
                         minSize = 15, maxSize = 5000, nPermSimple = 10000, eps = 0)
  as_tibble(res) |>
    transmute(module = pathway, !!paste0("NES_",  t_col) := NES,
              !!paste0("padj_", t_col) := padj,
              !!paste0("size_", t_col) := size)
}

t_cols <- c("t_CTLvPHE", "t_CTLvMITO", "t_PHEvPHE_MITO")
mod_nes <- Reduce(\(a, b) left_join(a, b, by = "module"),
                  lapply(t_cols, run_module_fgsea)) |>
  mutate(n_proteins = mod_sizes[module])

# YvO-style M1/M2/... module IDs + biology theme (matches mod_bio_labels.csv
# pattern). Display label = "Mn: bio_label" wrapped for in-plot legibility.
mod_ids <- setNames(paste0("M", seq_along(mod_levels)), mod_levels)
mod_nes <- mod_nes |>
  mutate(bio_label = mod_bio_labels[module],
         module_id = mod_ids[module],
         display   = stringr::str_wrap(
           paste0(module_id, ": ", gsub("\n", " ", bio_label)),
           width = 16))

build_module_scatter <- function(df, c_x, c_y, x_lab, y_lab, quad_labels) {
  nes_x <- paste0("NES_t_",  c_x); nes_y <- paste0("NES_t_",  c_y)
  x_vals <- df[[nes_x]]; y_vals <- df[[nes_y]]
  nes_lim <- max(abs(c(x_vals, y_vals)), na.rm = TRUE) * 1.35

  sp <- cor.test(x_vals, y_vals, method = "spearman")
  subtitle_txt <- sprintf(
    "n = %d modules | rho = %.2f%s",
    nrow(df), sp$estimate,
    if (sp$p.value < 0.001) ", p < 0.001"
    else sprintf(", p = %.3f", sp$p.value))

  q_tr <- sum(x_vals > 0 & y_vals > 0, na.rm = TRUE)
  q_bl <- sum(x_vals < 0 & y_vals < 0, na.rm = TRUE)
  q_tl <- sum(x_vals < 0 & y_vals > 0, na.rm = TRUE)
  q_br <- sum(x_vals > 0 & y_vals < 0, na.rm = TRUE)

  df$label_col <- ifelse(vapply(df$module, is_light_color, logical(1)),
                         "black", "white")
  df$label_col[df$module %in% c("green", "yellow")] <- "black"

  ggplot(df, aes(x = .data[[nes_x]], y = .data[[nes_y]])) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
             fill = quad_labels$fill[1], alpha = 0.20) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
             fill = quad_labels$fill[2], alpha = 0.20) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
             fill = quad_labels$fill[3], alpha = 0.20) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
             fill = quad_labels$fill[4], alpha = 0.20) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_abline(slope = quad_labels$slope, intercept = 0,
                linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_point(aes(size = n_proteins), fill = df$module, color = "black",
               shape = 21, alpha = 0.85, stroke = 0.5) +
    geom_label_repel(aes(label = display),
                     fill = scales::alpha(df$module, 0.85),
                     color = df$label_col,
                     size = 5.0, fontface = "bold", lineheight = 0.7,
                     max.overlaps = Inf, segment.size = 0.3,
                     segment.color = "grey40", min.segment.length = 0,
                     box.padding = 1.4, point.padding = 1.0,
                     force = 100, force_pull = 0.04,
                     label.padding = unit(1.5, "pt"),
                     label.r = unit(1.5, "pt"),
                     label.size = 0, seed = 7,
                     show.legend = FALSE) +
    annotate("label", x = Inf, y = Inf,
             label = sprintf("%s  n=%d", quad_labels$label[1], q_tr),
             hjust = 1, vjust = 1, size = 5, fontface = "bold",
             color = quad_labels$text_col[1], fill = alpha("white", 0.92),
             label.padding = unit(2, "pt")) +
    annotate("label", x = -Inf, y = -Inf,
             label = sprintf("%s  n=%d", quad_labels$label[2], q_bl),
             hjust = 0, vjust = 0, size = 5, fontface = "bold",
             color = quad_labels$text_col[2], fill = alpha("white", 0.92),
             label.padding = unit(2, "pt")) +
    annotate("label", x = Inf, y = -Inf,
             label = sprintf("%s  n=%d", quad_labels$label[3], q_br),
             hjust = 1, vjust = 0, size = 5, fontface = "bold",
             color = quad_labels$text_col[3], fill = alpha("white", 0.92),
             label.padding = unit(2, "pt")) +
    annotate("label", x = -Inf, y = Inf,
             label = sprintf("%s  n=%d", quad_labels$label[4], q_tl),
             hjust = 0, vjust = 1, size = 5, fontface = "bold",
             color = quad_labels$text_col[4], fill = alpha("white", 0.92),
             label.padding = unit(2, "pt")) +
    scale_size_continuous(range = c(3, 9), name = "Proteins",
                          breaks = c(100, 300, 600, 1000)) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    coord_fixed(ratio = 1,
                xlim = c(-nes_lim, nes_lim),
                ylim = c(-nes_lim, nes_lim),
                clip = "off") +
    labs(title = NULL, subtitle = subtitle_txt, x = x_lab, y = y_lab) +
    FIG_THEME +
    theme(axis.text     = element_text(size = 11, face = "bold", color = "grey30"),
          axis.title.x  = element_text(size = 12, face = "bold"),
          axis.title.y  = element_text(size = 12, face = "bold",
                                       margin = margin(r = -4)),
          plot.subtitle = element_text(size = 9, color = "grey30",
                                       face = "bold.italic"),
          legend.position = "none",
          plot.margin   = margin(2, 4, 2, -8))
}

# Top scatter: Disease vs Rescue (reversal — mirrors YvO Aging vs Tr.(O)).
# Quadrant fills: warm (red) for exacerbation diagonal, green for reversal
# off-diagonal. Text colour echoes YvO's red/blue quad-label scheme.
quad_rev <- list(
  label    = c("Exacerbated", "Exacerbated", "Reversed", "Reversed"),
  fill     = c("#D6604D", "#D6604D", "#2E7D32", "#2E7D32"),
  text_col = c("#D6604D", "#D6604D", "#4393C3", "#4393C3"),
  slope    = -1)
pB_top <- build_module_scatter(
  mod_nes, "CTLvPHE", "PHEvPHE_MITO",
  "NES (Disease: PHE - Ctl)",
  "NES (Rescue: PHE_Mito - PHE)",
  quad_rev)

# Bottom scatter: Intervention vs Rescue (concordance — mirrors YvO Tr.(Y) vs
# Tr.(O)). Diagonal = concordant (warm fill), off-diagonal = discordant (cool).
quad_conc <- list(
  label    = c("Concordant Up", "Concordant Down", "Discordant", "Discordant"),
  fill     = c("#D6604D", "#D6604D", "#4393C3", "#4393C3"),
  text_col = c("#D6604D", "#D6604D", "#4393C3", "#4393C3"),
  slope    = 1)
pB_bot <- build_module_scatter(
  mod_nes, "CTLvMITO", "PHEvPHE_MITO",
  "NES (Intervention: Mito - Ctl)",
  "NES (Rescue: PHE_Mito - PHE)",
  quad_conc)

ggsave(file.path(PNL_PNG, "MAIN_panel_B_top_reversal.png"), pB_top,
       width = 220, height = 135, units = "mm", dpi = 300)
ggsave(file.path(PNL_PNG, "MAIN_panel_B_bot_concordance.png"), pB_bot,
       width = 220, height = 135, units = "mm", dpi = 300)

# Shared legend: protein-count size scale (matches YvO F06 panel B legend).
pB_leg_src <- pB_top +
  scale_size_continuous(range = c(3, 9), name = "Proteins per module",
                        breaks = c(100, 300, 600, 1000),
                        guide = guide_legend(nrow = 1,
                          override.aes = list(alpha = 0.7,
                                              fill = "grey60",
                                              stroke = 0))) +
  theme(legend.position   = "bottom",
        legend.title      = element_text(size = 9, face = "bold"),
        legend.text       = element_text(size = 8),
        legend.key        = element_rect(fill = NA, color = NA),
        legend.key.size   = unit(4, "mm"),
        legend.background = element_rect(fill = NA, color = NA))
legend_grob <- cowplot::get_plot_component(pB_leg_src, "guide-box-bottom",
                                            return_all = FALSE)
pB_leg <- wrap_elements(full = legend_grob)

pB_diptych <- (pB_top / pB_bot) + plot_layout(heights = c(1, 1))

# Save panel B (2-scatter diptych) + legend as separate PNGs for rasterGrob
# mount into the final composite.
pB_scatters_png <- file.path(PNL_PNG, "MAIN_panel_B_scatters.png")
pB_legend_png   <- file.path(PNL_PNG, "MAIN_panel_B_legend.png")
ggsave(pB_scatters_png, pB_diptych,
       width = 220, height = 270, units = "mm", dpi = 300)
ggsave(pB_legend_png, wrap_elements(full = legend_grob),
       width = 90, height = 16, units = "mm", dpi = 300)

# Composite (470 × 300 mm landscape, YvO F06 rasterGrob layout)

pA_grob     <- rasterGrob(readPNG(ht_png),         interpolate = TRUE)
pB_grob     <- rasterGrob(readPNG(pB_scatters_png), interpolate = TRUE)
pB_leg_grob <- rasterGrob(readPNG(pB_legend_png),   interpolate = TRUE)

layout_cfg <- list(
  comp_w = 470, comp_h = 300,
  tag_sz      = 14,
  title_sz    = 12,
  subtitle_sz = 9,
  grid_top    = 0.80,
  grid_bot    = 0.22,
  b_x         = 0.43,
  b_w         = 0.55,
  b_leg_x = 0.70, b_leg_y = 0.195, b_leg_w = 0.24, b_leg_h = 0.04,
  a_grob_x = 0.01, a_grob_y = 0.18, a_grob_w = 0.60, a_grob_h = 0.62,
  a_let_x = 15,  a_let_y = 248,
  a_ttl_x = 58,  a_ttl_y = 248,
  a_sub_x = 58,  a_sub_y = 243,
  b_let_x = 285, b_let_y = 248,
  b_ttl_x = 294, b_ttl_y = 248,
  b_sub_x = 294, b_sub_y = 243,
  crop_l = 0.01, crop_r = 0.80,
  crop_b = 0.17, crop_t = 0.85)

COMP_W      <- layout_cfg$comp_w
COMP_H      <- layout_cfg$comp_h
TAG_SZ      <- layout_cfg$tag_sz
TITLE_SZ    <- layout_cfg$title_sz
SUBTITLE_SZ <- layout_cfg$subtitle_sz
GRID_TOP    <- layout_cfg$grid_top
GRID_BOT    <- layout_cfg$grid_bot
B_X         <- layout_cfg$b_x
B_W         <- layout_cfg$b_w
B_GRID_H    <- GRID_TOP - GRID_BOT

mm2x <- function(mm) mm / COMP_W
mm2y <- function(mm) mm / COMP_H

CROP_L <- layout_cfg$crop_l; CROP_R <- layout_cfg$crop_r
CROP_B <- layout_cfg$crop_b; CROP_T <- layout_cfg$crop_t
SAVE_W <- COMP_W * (CROP_R - CROP_L)
SAVE_H <- COMP_H * (CROP_T - CROP_B)

# Honesty: modules are data-derived from the same 24-sample matrix then tested
# for differential eigengenes (descriptive, not a preservation/stability claim).
# N=24 is too small for modulePreservation/Zsummary, so none is asserted.
a_title    <- "Module × contrast (eigengene limma) — exploratory, N = 24"
a_subtitle <- sprintf("EXPLORATORY — %d modules, data-derived and tested on the same 24 samples (no preservation test feasible) | rows by total signal | * FDR<0.05  ** <0.01  *** <0.001",
                      length(mod_levels))
b_title    <- "Module-Level NES Scatters"
b_subtitle <- "fGSEA on module-member t-stat ranks | reversal (top) + concordance (bottom) | point size = proteins per module"

composite <- ggdraw(xlim = c(CROP_L, CROP_R), ylim = c(CROP_B, CROP_T)) +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  draw_grob(pB_grob, x = B_X, y = GRID_BOT, width = B_W, height = B_GRID_H,
            hjust = 0, vjust = 0) +
  draw_grob(pB_leg_grob, x = layout_cfg$b_leg_x, y = layout_cfg$b_leg_y,
            width = layout_cfg$b_leg_w, height = layout_cfg$b_leg_h,
            hjust = 0.5, vjust = 0.5) +
  draw_grob(pA_grob, x = layout_cfg$a_grob_x, y = layout_cfg$a_grob_y,
            width = layout_cfg$a_grob_w, height = layout_cfg$a_grob_h,
            hjust = 0, vjust = 0) +
  draw_label("A", x = mm2x(layout_cfg$a_let_x), y = mm2y(layout_cfg$a_let_y),
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(a_title,
             x = mm2x(layout_cfg$a_ttl_x), y = mm2y(layout_cfg$a_ttl_y),
             size = TITLE_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(a_subtitle,
             x = mm2x(layout_cfg$a_sub_x), y = mm2y(layout_cfg$a_sub_y),
             size = SUBTITLE_SZ, fontface = "bold.italic", colour = "grey40",
             hjust = 0, vjust = 1) +
  draw_label("B", x = mm2x(layout_cfg$b_let_x), y = mm2y(layout_cfg$b_let_y),
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(b_title,
             x = mm2x(layout_cfg$b_ttl_x), y = mm2y(layout_cfg$b_ttl_y),
             size = TITLE_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(b_subtitle,
             x = mm2x(layout_cfg$b_sub_x), y = mm2y(layout_cfg$b_sub_y),
             size = SUBTITLE_SZ, fontface = "bold.italic", colour = "grey40",
             hjust = 0, vjust = 1)

pdf_device <- get_pdf_device()
ggsave(file.path(RPT_PDF, "SUPP_F05_wgcna.pdf"), composite,
       width = SAVE_W, height = SAVE_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "SUPP_F05_wgcna.png"), composite,
       width = SAVE_W, height = SAVE_H, units = "mm", dpi = 300)
message("F05 WGCNA supplement (470 × 300, rasterGrob layout) saved")

# Supplementary xlsx

mod_assign <- w$ann |>
  mutate(module = w$module_colors[match(uniprot_id, w$gene_order)])

build_workbook(
  file.path(DAT, "F05_supplementary.xlsx"),
  sheet_specs = list(
    list(name = "panel_A_eigengene_limma",
         df = as.data.frame(mod_results)),
    list(name = "panel_A_module_size_hubs",
         df = as.data.frame(hub_top3 |>
           left_join(tibble(module = names(mod_sizes),
                            n_proteins = as.integer(mod_sizes)),
                     by = "module"))),
    list(name = "panel_B_module_NES",  df = as.data.frame(mod_nes)),
    list(name = "module_assignments",  df = as.data.frame(mod_assign)),
    list(name = "hubs",                df = as.data.frame(w$hubs))))

message(sprintf("F05 2-panel composite done | %d non-grey modules | %d hub records",
                length(mod_levels), nrow(w$hubs)))
