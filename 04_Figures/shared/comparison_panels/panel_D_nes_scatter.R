# Shared NES Scatter Panel — Panel D in F04 (concordance) and F05 (reversal)

source(here::here("04_Figures", "shared", "style.R"))
source(here::here("04_Figures", "shared", "print_scale_apply_380mm.R"))
source(here::here("04_Figures", "shared", "pathway_utils.R"))
source(here::here("04_Figures", "shared", "mitocarta_utils.R"))  # MITO_PATHWAY_REGEX

library(tidyverse)
library(ggrepel)

PG_W <- cfg$panel_w %||% 146
RPT_PNG <- cfg$rpt_png
RPT_PDF <- cfg$rpt_pdf
DAT     <- cfg$dat
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_D"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

fgsea_cache <- cfg$fgsea_cache %||% here::here("04_Figures", "shared", "fgsea_tstat_all_h9c2.csv")
stopifnot("fGSEA cache missing" = file.exists(fgsea_cache))
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

databases_keep <- cfg$databases %||% c("Hallmark", "GO Slim")

# GO Slim lives in an additive cache (not the main fGSEA cache). Bind it in when
# requested so it is selectable like Hallmark/Reactome/MitoCarta.
goslim_csv <- here::here("04_Figures", "shared", "fgsea_goslim_h9c2.csv")
if ("GO Slim" %in% databases_keep && file.exists(goslim_csv))
  fgsea_all <- bind_rows(fgsea_all, read_csv(goslim_csv, show_col_types = FALSE))

cx <- cfg$contrast_x
cy <- cfg$contrast_y
nes_x <- paste0("NES_", cx)
nes_y <- paste0("NES_", cy)
padj_x <- paste0("padj_", cx)
padj_y <- paste0("padj_", cy)
size_x <- paste0("size_", cx)
size_y <- paste0("size_", cy)

fgsea_hg <- fgsea_all |>
  filter(database %in% databases_keep,
         contrast %in% c(cx, cy),
         # drop bare MitoCarta compartment sets (localizations, not pathways)
         !pathway %in% c("MITOCARTA_IMM", "MITOCARTA_IMS",
                         "MITOCARTA_MATRIX", "MITOCARTA_OMM"),
         # drop pathogen/disease sets (irrelevant; leak in as high-NES noise)
         !grepl(DISEASE_VIRAL_RE, pathway, ignore.case = TRUE))

# Optional mito-lens regex filter: keep MitoCarta entries unconditionally,
# keep other DBs only when pathway name matches cfg$pathway_regex.
if (!is.null(cfg$pathway_regex) && nzchar(cfg$pathway_regex)) {
  fgsea_hg <- fgsea_hg |>
    filter(database == "MitoCarta" |
             grepl(cfg$pathway_regex, pathway, perl = TRUE))
}

fgsea_wide <- fgsea_hg |>
  dplyr::select(pathway, contrast, NES, padj, size, database) |>
  pivot_wider(id_cols = c(pathway, database), names_from = contrast,
              values_from = c(NES, padj, size)) |>
  filter(!is.na(.data[[nes_x]]), !is.na(.data[[nes_y]])) |>
  mutate(set_size = coalesce(.data[[size_x]], .data[[size_y]]))

fgsea_wide <- fgsea_wide |>
  mutate(
    sig_1 = !is.na(.data[[padj_x]]) & .data[[padj_x]] < 0.05,
    sig_2 = !is.na(.data[[padj_y]]) & .data[[padj_y]] < 0.05,
    significance = case_when(
      sig_1 & sig_2 ~ cfg$quadrant_defs$sig_both_label,
      sig_1         ~ cfg$quadrant_defs$sig_x_label,
      sig_2         ~ cfg$quadrant_defs$sig_y_label,
      TRUE          ~ "NS"
    ) |> factor(levels = names(cfg$sig_colors)),
    pathway_label = clean_pathway_name(pathway),
    db_shape = ifelse(database == "Hallmark", 24, 21),
    # mitochondrial pathway? -> outline ring (matches the protein-scatter mito ring)
    is_mito_pw = grepl("^MITOCARTA_", pathway) |
                 grepl(MITO_PATHWAY_REGEX, pathway, perl = TRUE)
  )

fgsea_sig <- fgsea_wide |> filter(significance != "NS")

db_summary <- table(fgsea_wide$database)
message(sprintf("  %d total pathways (%s) | %d significant",
                nrow(fgsea_wide),
                paste(sprintf("%s: %d", names(db_summary), as.integer(db_summary)),
                      collapse = ", "),
                nrow(fgsea_sig)))

nes_cor_all <- cor.test(fgsea_wide[[nes_x]], fgsea_wide[[nes_y]], method = "spearman")
nes_ci_all  <- fisher_z_ci(nes_cor_all$estimate, nrow(fgsea_wide))
nes_cor_sig <- if (nrow(fgsea_sig) >= 3) {
  cor.test(fgsea_sig[[nes_x]], fgsea_sig[[nes_y]], method = "spearman")
} else NULL

# wider margin: pulls the extreme points inward so the corner quadrant labels sit
# in empty space (no point overlap)
nes_lim <- max(abs(c(fgsea_wide[[nes_x]], fgsea_wide[[nes_y]]))) * 1.35

qd <- cfg$quadrant_defs
n_q1 <- sum(fgsea_sig[[nes_x]] > 0 & fgsea_sig[[nes_y]] > 0)
n_q2 <- sum(fgsea_sig[[nes_x]] < 0 & fgsea_sig[[nes_y]] > 0)
n_q3 <- sum(fgsea_sig[[nes_x]] < 0 & fgsea_sig[[nes_y]] < 0)
n_q4 <- sum(fgsea_sig[[nes_x]] > 0 & fgsea_sig[[nes_y]] < 0)
# totals (all tested pathways) per quadrant -> annotation reads "sig / total"
nt_q1 <- sum(fgsea_wide[[nes_x]] > 0 & fgsea_wide[[nes_y]] > 0)
nt_q2 <- sum(fgsea_wide[[nes_x]] < 0 & fgsea_wide[[nes_y]] > 0)
nt_q3 <- sum(fgsea_wide[[nes_x]] < 0 & fgsea_wide[[nes_y]] < 0)
nt_q4 <- sum(fgsea_wide[[nes_x]] > 0 & fgsea_wide[[nes_y]] < 0)
n_metric <- qd$metric_count_fn(n_q1, n_q2, n_q3, n_q4)
n_total_sig <- nrow(fgsea_sig)
metric_frac <- if (n_total_sig > 0) n_metric / n_total_sig else 0

message(sprintf("  NES Spearman (all): rho = %.3f [%.3f, %.3f]",
                nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2]))

txt_pw   <- scale_text(BASE_PATHWAY, PG_W) * 0.80   # smaller — readable in the compact panel
txt_quad <- scale_text(BASE_QUADRANT, PG_W) * 1.15

# Label selection — UNBIASED & SYSTEMATIC (no hand-picked whitelist):
#   1. Collapse gene-set-redundant pathways across databases by Jaccard overlap
#      (single threshold), keeping the representative with the largest combined
#      effect. "OXPHOS / OXPHOS subunits / Complex V" merge because their genes
#      overlap, not because they were named.
#   2. Label the top-k survivors PER POPULATED QUADRANT (sign of NES_x, NES_y) by
#      combined magnitude |NES_x| + |NES_y| (ties -> lower p_adj). Quadrant-balanced
#      so the dominant corner can't monopolise the budget and the weaker arms
#      (discordant / down) still get labelled.
# Selection is a deterministic function of the data + (jaccard, k); both reported.
N_PER_QUAD <- cfg$n_label_per_quad %||% 4L
.gs_main <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))
.gs_goslim_path <- here::here("04_Figures", "shared", "goslim_rat_gene_sets.rds")
.gs_goslim <- if (file.exists(.gs_goslim_path)) readRDS(.gs_goslim_path) else list()
.pw_dedup <- c(.gs_main$Hallmark, .gs_main$Reactome, .gs_main$KEGG,
               .gs_main[["GO:BP"]], .gs_main$MitoCarta, .gs_goslim)

label_pool <- fgsea_sig |>
  mutate(padj = pmin(.data[[padj_x]], .data[[padj_y]], na.rm = TRUE),
         nes_mag = abs(.data[[nes_x]]) + abs(.data[[nes_y]]))
label_pool <- deduplicate_enrichment(as.data.frame(label_pool),
                                     pathways = .pw_dedup,
                                     jaccard_cutoff = cfg$label_jaccard %||% 0.4,
                                     cross_db = TRUE) |> as_tibble()

label_pw <- label_pool |>
  mutate(.quad = paste0(ifelse(.data[[nes_x]] >= 0, "xpos", "xneg"), "_",
                        ifelse(.data[[nes_y]] >= 0, "ypos", "yneg"))) |>
  group_by(.quad) |>
  arrange(desc(nes_mag), padj, .by_group = TRUE) |>
  slice_head(n = N_PER_QUAD) |>
  ungroup() |>
  select(-nes_mag, -padj, -.quad) |>
  mutate(
    label_fill     = cfg$sig_label_fill[as.character(significance)],
    label_text_col = cfg$sig_label_text[as.character(significance)],
    # shared cleaner: MitoCarta hierarchy -> last segment; DB overrides applied
    pathway_label  = clean_display_label(pathway, extra = cfg$display_overrides),
    pathway_label  = str_wrap(str_trunc(pathway_label, 26), width = 16)
  )
message(sprintf("  NES labels: %d shown (<= %d/quadrant, Jaccard %.2f collapse)",
                nrow(label_pw), N_PER_QUAD, cfg$label_jaccard %||% 0.4))

ns_df  <- fgsea_wide |> filter(significance == "NS")
sig_df <- fgsea_wide |> filter(significance != "NS") |>
  mutate(draw_order = factor(significance, levels = cfg$sig_draw_order)) |>
  arrange(draw_order)

rho_sig_str <- if (!is.null(nes_cor_sig)) sprintf(", \u03c1(sig) = %.2f", nes_cor_sig$estimate) else ""
db_label <- paste(databases_keep, collapse = " + ")
n_mito_pw <- sum(fgsea_wide$is_mito_pw)
subtitle_str <- sprintf(
  "%s | %d pathways (%d sig.) | \u03c1 = %.2f [%.2f, %.2f], %s%s\n%.0f%% %s | %s | orange ring = mitochondrial pathway (n=%d)",
  db_label, nrow(fgsea_wide), n_total_sig,
  nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
  ifelse(nes_cor_all$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", nes_cor_all$p.value)),
  rho_sig_str, metric_frac * 100,
  cfg$subtitle_metric, cfg$subtitle_interpretation, n_mito_pw
)

pD <- ggplot(mapping = aes(x = .data[[nes_x]], y = .data[[nes_y]])) +
  annotate("rect", xmin = qd$bg_red_1[1], xmax = qd$bg_red_1[2],
           ymin = qd$bg_red_1[3], ymax = qd$bg_red_1[4],
           fill = cfg$bg_concordant %||% "#D6604D", alpha = 0.20, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = qd$bg_red_2[1], xmax = qd$bg_red_2[2],
           ymin = qd$bg_red_2[3], ymax = qd$bg_red_2[4],
           fill = cfg$bg_concordant %||% "#D6604D", alpha = 0.20, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = qd$bg_blue_1[1], xmax = qd$bg_blue_1[2],
           ymin = qd$bg_blue_1[3], ymax = qd$bg_blue_1[4],
           fill = cfg$bg_discordant %||% "#4393C3", alpha = 0.20, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = qd$bg_blue_2[1], xmax = qd$bg_blue_2[2],
           ymin = qd$bg_blue_2[3], ymax = qd$bg_blue_2[4],
           fill = cfg$bg_discordant %||% "#4393C3", alpha = 0.20, color = "grey70", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = cfg$ref_slope, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_point(data = ns_df, aes(shape = database),
             size = 1.5, fill = "grey70",
             color = ifelse(ns_df$is_mito_pw, unname(LOC_COLORS["Mitochondrial"]), "grey55"),
             alpha = 0.40, stroke = ifelse(ns_df$is_mito_pw, 0.5, 0.2)) +
  # mitochondrial pathways get an orange ring (matches the protein-scatter mito ring)
  geom_point(data = sig_df, aes(fill = significance, size = set_size, shape = database),
             color = ifelse(sig_df$is_mito_pw, unname(LOC_COLORS["Mitochondrial"]), "grey55"),
             alpha = 0.80, stroke = ifelse(sig_df$is_mito_pw, 1.1, 0.35)) +
  scale_fill_manual(values = cfg$sig_colors, name = "Significance") +
  scale_shape_manual(values = cfg$db_shapes %||% c("Hallmark" = 24, "GO Slim" = 21),
                     name = "Database") +
  scale_size_continuous(range = c(2.5, 8), name = "Set size",
                        breaks = c(20, 50, 100, 200)) +
  geom_label_repel(data = label_pw, aes(label = pathway_label),
                   fill = label_pw$label_fill, color = label_pw$label_text_col,
                   size = txt_pw, fontface = "bold", lineheight = 0.85,
                   max.overlaps = Inf,
                   segment.size = 0.4, segment.color = "grey20",
                   min.segment.length = 0, show.legend = FALSE,
                   # stronger repel + larger point padding so labels clear the
                   # points; xlim/ylim keep them out of the corner bands reserved
                   # for the quadrant labels
                   box.padding = 1.9, point.padding = 1.0,
                   force = 90, force_pull = 0.08,
                   xlim = c(-nes_lim, nes_lim) * 0.86,
                   ylim = c(-nes_lim, nes_lim) * 0.86,
                   label.padding = unit(1, "pt"),
                   label.r = unit(0.5, "pt"),
                   label.size = cfg$label_border_size %||% 0.10, seed = 42) +
  annotate("label", x = nes_lim, y = nes_lim,
           label = sprintf("%s\n%d / %d sig", qd$label_tr, n_q1, nt_q1),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold", lineheight = 0.9,
           color = qd$color_tr, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -nes_lim, y = nes_lim,
           label = sprintf("%s\n%d / %d sig", qd$label_tl, n_q2, nt_q2),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold", lineheight = 0.9,
           color = qd$color_tl, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -nes_lim, y = -nes_lim,
           label = sprintf("%s\n%d / %d sig", qd$label_bl, n_q3, nt_q3),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold", lineheight = 0.9,
           color = qd$color_bl, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = nes_lim, y = -nes_lim,
           label = sprintf("%s\n%d / %d sig", qd$label_br, n_q4, nt_q4),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold", lineheight = 0.9,
           color = qd$color_br, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  coord_fixed(ratio = 1, xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
  labs(title = cfg$title,
       subtitle = subtitle_str,
       x = sprintf("NES (%s)", cfg$axis_x_label),
       y = sprintf("NES (%s)", cfg$axis_y_label)) +
  FIG_THEME +
  theme(
    axis.text         = element_text(size = FIG_AXIS_TEXT, face = "bold", color = "grey30"),
    axis.title        = element_text(size = FIG_AXIS_TEXT, face = "bold"),
    legend.position   = "bottom",
    legend.title      = element_text(size = FIG_LEGEND_TITLE, face = "bold", color = "grey25"),
    legend.text       = element_text(size = FIG_LEGEND_TEXT, color = "grey20"),
    legend.key.size   = unit(2 * PRINT_SCALE, "mm"),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box        = "horizontal",
    legend.box.just   = "center",
    legend.spacing.x  = unit(3 * PRINT_SCALE, "mm"),
    legend.box.margin = margin(t = -2),
    plot.margin       = margin(0, 0, 0, 0)
  ) +
  guides(fill  = "none",
         shape = guide_legend(nrow = 1, order = 1,
                               keyheight = unit(4 * PRINT_SCALE, "mm"),
                               keywidth  = unit(4 * PRINT_SCALE, "mm"),
                               override.aes = list(size = 3 * PRINT_SCALE, fill = "grey50")),
         size  = guide_legend(nrow = 1, order = 2,
                               keyheight = unit(4 * PRINT_SCALE, "mm"),
                               keywidth  = unit(4 * PRINT_SCALE, "mm")))

ggsave(file.path(RPT_PNG, "MAIN_panel_D_nes_scatter.png"), pD,
       width = PG_W, height = PG_W, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_D_nes_scatter.pdf"), pD,
       width = PG_W, height = PG_W, units = "mm", device = pdf_device)

export_df <- fgsea_wide |>
  transmute(
    pathway, pathway_label, database,
    !!paste0("NES_", cx) := round(.data[[nes_x]], 3),
    !!paste0("NES_", cy) := round(.data[[nes_y]], 3),
    !!paste0("padj_", cx) := signif(.data[[padj_x]], 4),
    !!paste0("padj_", cy) := signif(.data[[padj_y]], 4),
    significance = as.character(significance),
    set_size
  ) |>
  arrange(significance, desc(abs(.data[[nes_x]]) + abs(.data[[nes_y]])))
write_csv(export_df, file.path(DAT, "panel_D", "nes_scatter.csv"))

pD_legend_grob <- cowplot::get_plot_component(pD, "guide-box-bottom", return_all = FALSE)
if (!is.null(pD_legend_grob)) {
  pD_legend_plot <- cowplot::ggdraw(pD_legend_grob)
  ggsave(file.path(RPT_PNG, "MAIN_panel_D_legend.png"), pD_legend_plot,
         width = 120, height = 14, units = "mm", dpi = 300)
}

pD_title    <- cfg$title
pD_subtitle <- subtitle_str
pD_legend   <- NULL
# Strip titles but KEEP legend (legend provides shape/size key for composite)
pD          <- pD + labs(title = NULL, subtitle = NULL, tag = NULL)

pw_conc_frac <- metric_frac  # F04 stitcher
pw_rev_frac  <- metric_frac  # F05 stitcher

cat(sprintf("%s Panel D done\n", cfg$fig_id))
