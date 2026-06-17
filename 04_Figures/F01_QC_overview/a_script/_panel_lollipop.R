# F01 — leading-edge lollipop, faceted across the factorial contrast set.
# Top recurrent leading-edge proteins (across significant pathways) shown as a
# horizontal lollipop of logFC, faceted by 4 contrasts, coloured by GO cellular
# compartment (org.Rn.eg.db GO-CC). Shared gene order across facets so a
# protein's behaviour can be compared across contrasts.
#
# Contrast set (user-specified factorial lens):
#   Mito - Ctl                    (Transplant_Mito)
#   PHE  - Ctl                    (Disease_Phe)
#   PHE_Mito - Mito               (PHE effect in transplanted cells)
#   (PHE_Mito - Mito) - (PHE - Ctl)  == Interaction_Mito
# Sourced from 01_main_panels.R. Inherits dep_df/comb, fgsea_all, DAT, PNL_PNG,
#   FIG_* sizes, contrast_brief, CONTRAST_MATH_BRIEF.

# GO Slim functional annotation (BP slim consolidated to 15 muscle-relevant
# categories; far more interpretable for leading-edge proteins than raw GO-CC).
source(here::here("04_Figures", "shared", "go_slim_categories.R"))
LOLLI_CORE <- c("CTLvMITO", "CTLvPHE", "MITOvPHE_MITO", "Interaction")
LOLLI_LABS <- setNames(sprintf("%s\n(%s)", contrast_brief(LOLLI_CORE),
                               c("Mito−Ctl", "PHE−Ctl", "PHE_Mito−Mito", "(PHE_Mito−Mito)−(PHE−Ctl)")),
                       LOLLI_CORE)

# Top leading-edge proteins PER contrast (faceted by which contrast they drive),
# ranked by recurrence within that contrast's significant pathways.
# DB set matches F02 + panel F (Hallmark + Reactome + MitoCarta); GO:BP + KEGG
# excluded for one consistent enrichment universe across the figure.
N_PER <- 7
le_pool <- fgsea_all |>
  filter(contrast %in% LOLLI_CORE, !is.na(padj), padj < 0.05,
         database %in% c("Hallmark", "Reactome", "MitoCarta")) |>
  tidyr::separate_rows(leadingEdge, sep = ";") |>
  filter(leadingEdge != "", !is.na(leadingEdge))
top_per <- bind_rows(lapply(LOLLI_CORE, function(ctr) {
  le_pool |> filter(contrast == ctr) |> count(gene = leadingEdge, sort = TRUE) |>
    filter(gene %in% dep_df$gene) |> slice_head(n = N_PER) |> mutate(contrast = ctr)
}))
genes_all <- unique(top_per$gene)
slim_tbl <- assign_go_slim_consolidated(genes_all, unique(dep_df$gene), min_cat_size = 1)
gene_comp <- setNames(as.character(slim_tbl$consolidated), slim_tbl$gene)

# logFC of each gene in its facet's contrast.
lf_long <- dep_df |> distinct(gene, .keep_all = TRUE) |>
  dplyr::select(gene, all_of(paste0("logFC_", LOLLI_CORE))) |>
  tidyr::pivot_longer(-gene, names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = sub("^logFC_", "", contrast))
lolli_long <- top_per |> left_join(lf_long, by = c("gene", "contrast")) |>
  mutate(compartment = factor(gene_comp[gene], levels = CONSOLIDATED_PATHWAY_ORDER),
         contrast = factor(contrast, levels = LOLLI_CORE))

# reorder_within: order genes by logFC inside each facet (manual, no tidytext).
lolli_long <- lolli_long |> arrange(contrast, logFC) |>
  mutate(gene_o = factor(paste(gene, contrast, sep = "___"),
                         levels = unique(paste(gene, contrast, sep = "___"))))
strip_key <- function(x) sub("___.*", "", x)
write_csv(lolli_long |> dplyr::select(gene, contrast, logFC, compartment),
          file.path(DAT, "panel_lollipop_leadingedge.csv"))
N_LOLLI <- length(genes_all)

PL_W <- 88; PL_H <- 120
# Per-facet background wash in each contrast's accent colour (replaces titles).
bg_lol <- tibble(contrast = factor(LOLLI_CORE, levels = LOLLI_CORE),
                 fill = unname(CONTRAST_COLORS[LOLLI_CORE]))
pLOL <- ggplot(lolli_long, aes(logFC, gene_o)) +
  geom_rect(data = bg_lol, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = bg_lol$fill, alpha = 0.16, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey55") +
  geom_segment(aes(x = 0, xend = logFC, yend = gene_o, color = compartment), linewidth = 0.45) +
  geom_point(aes(color = compartment), size = 1.5) +
  facet_wrap(~ contrast, ncol = 1, scales = "free_y") +
  scale_color_manual(values = CONSOLIDATED_COLORS, name = "GO Slim category", drop = TRUE) +
  scale_y_discrete(labels = strip_key) +
  scale_x_continuous(breaks = scales::pretty_breaks(4)) +
  labs(title = "Leading-edge proteins",
       subtitle = sprintf("top %d per contrast | x = log₂FC | colour = GO Slim function", N_PER),
       x = expression(bold(log[2]~FC)), y = NULL, tag = "B") +
  FIG_THEME +
  theme(strip.text = element_blank(), strip.background = element_blank(),
        axis.text.y = element_text(size = FIG_AXIS_TEXT - 1, face = "italic"),
        panel.spacing.y = unit(2, "pt"),
        legend.position = "bottom", legend.key.size = unit(2.4, "mm"),
        legend.title = element_text(size = FIG_LEGEND_TITLE - 0.5, hjust = 0.5),
        legend.text = element_text(size = FIG_LEGEND_TEXT),
        plot.margin = margin(5, 4, 1, 2)) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 2))

ggsave(file.path(PNL_PNG, "MAIN_panel_B_lollipop.png"), pLOL, width = PL_W, height = PL_H, units = "mm", dpi = 300)
message(sprintf("F01 leading-edge lollipop saved | %d proteins, compartments: %s",
                N_LOLLI, paste(names(table(gene_comp))[table(gene_comp) > 0], collapse = ", ")))
