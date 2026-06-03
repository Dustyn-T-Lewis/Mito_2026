#!/usr/bin/env Rscript
# F04 main composite: Disease_Phe (CTLvPHE) vs Rescue_Mito+Phe (PHEvPHE_MITO) REVERSAL.
# Question: does mito transplant reverse PHE-driven remodeling? Negative ρ = rescue.
# 4-panel 2-row layout: A protein quadrant (NEW), B pattern heatmap (top);
# D fGSEA NES scatter, E RRHO2 (bottom). Panel B sources AnnotationDbi via
# go_slim_categories.R, so it runs LAST to avoid select() S4 masking.
# Ported from 04_Figures/F05 (already Disease vs Rescue reversal); the ORA
# composite panel A is replaced by the top-5-per-quadrant protein scatter, fry C dropped.
# Mirrors YvO F05 (Aging Reversal) visual/format. Cite: Doulamis 2024 (PMID 39732955).

library(dplyr); library(tidyr); library(tibble); library(stringr); library(readr)
library(ggplot2); library(patchwork); library(cowplot)
source(here::here("04_Figures", "shared", "style.R"))

BASE    <- "05_Figures/F04_rescue"
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
PNL_PNG <- file.path(RPT_PNG, "panels")
PNL_PDF <- file.path(RPT_PDF, "panels")
DAT     <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, PNL_PNG, PNL_PDF, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()
message("=== F04 Composite (reversal): sourcing panels ===")

source(here::here("05_Figures", "F04_rescue", "a_script", "_panel_A_quadrant.R"))

# Panel D: fGSEA NES scatter (reversal; ref_slope = -1)
cfg <- list(
  fig_id = "F04", contrast_x = "CTLvPHE", contrast_y = "PHEvPHE_MITO",
  fgsea_cache = "04_Figures/shared/fgsea_tstat_all_h9c2.csv",
  databases = c("Hallmark", "GO Slim"), db_shapes = c("Hallmark" = 24, "GO Slim" = 21),
  title = "Pathway-Level Reversal (fGSEA)",
  axis_x_label = "Disease (PHE - Ctl)", axis_y_label = "Rescue (PHE+Mito - PHE)",
  subtitle_metric = "reversed",
  subtitle_interpretation = "negative ρ = mito-transplant opposes PHE remodeling",
  ref_slope = -1, panel_w = 146, label_border_size = 0.25,
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  sig_colors = SIG_COLORS_F3, sig_label_fill = SIG_LABEL_FILL_F3, sig_label_text = SIG_LABEL_TEXT_F3,
  sig_draw_order = c("Sig Rescue only", "Sig Disease only", "Sig Both"),
  quadrant_defs = list(
    sig_both_label = "Sig Both", sig_x_label = "Sig Disease only", sig_y_label = "Sig Rescue only",
    bg_blue_1 = c(0, Inf, -Inf, 0), bg_blue_2 = c(-Inf, 0, 0, Inf),
    bg_red_1 = c(0, Inf, 0, Inf), bg_red_2 = c(-Inf, 0, -Inf, 0),
    label_tr = "Exacerbated", color_tr = "#D6604D", label_tl = "Reversed", color_tl = "#4393C3",
    label_bl = "Exacerbated", color_bl = "#D6604D", label_br = "Reversed", color_br = "#4393C3",
    metric_count_fn = function(q1, q2, q3, q4) q2 + q4),
  display_overrides = c(
    "Oxidative Phosphorylation" = "OXPHOS",
    "Aerobic Respiration And Respiratory Electron Transport" = "Aerobic Resp. + ETC",
    "Citric Acid Cycle Tca Cycle" = "TCA Cycle", "Mitochondrial Translation" = "Mito Translation",
    "Mitochondrial Protein Degradation" = "Mito Protein Deg.", "Mitochondrial Calcium Ion Transport" = "Mito Ca²⁺ Transport",
    "Pink1 Prkn Mediated Mitophagy" = "PINK1/PRKN Mitophagy", "Fatty Acid Beta Oxidation" = "FA β-Oxidation",
    "Branched Chain Amino Acid Metabolism" = "BCAA Metabolism", "Cardiolipin Biosynthesis" = "Cardiolipin Biosynth.",
    "Mitochondrial Organization" = "Mito Org.", "Mitochondrial Transport" = "Mito Transport",
    "Mitochondrial Protein Import" = "Mito Protein Import", "Extracellular Matrix Organization" = "ECM Org.",
    "Protein Folding" = "Protein Folding", "Cytoplasmic Translation" = "Cytoplasmic Transl."))
source(here::here("04_Figures", "shared", "comparison_panels", "panel_D_nes_scatter.R"))
n_pw_D <- nrow(fgsea_wide); rho_D <- as.numeric(nes_cor_all$estimate); pw_rev_D <- pw_rev_frac

# Panel E: RRHO2 (reversal; warm off-diagonal)
cfg <- list(
  fig_id = "F04", t_col_1 = "t_CTLvPHE", t_col_2 = "t_PHEvPHE_MITO",
  rrho_labels = c("Disease", "Rescue"),
  title = "Threshold-Free Reversal (RRHO2)",
  subtitle_fmt = "Stratified hypergeometric | %d shared genes | warm off-diagonal = mito reverses PHE | No MTC (Cahill et al. 2018)",
  axis_label_1 = expression("Disease rank"~(Up %->% Down)), axis_label_2 = expression("Rescue rank"~(Up %->% Down)),
  quadrant_labels = list(UU = "Exacerbated Up", DD = "Exacerbated Down", UD = "Reversed (Dis↑ Re↓)", DU = "Reversed (Dis↓ Re↑)"),
  hotspot_export_names = list(UU = "Exacerbated_Up", DD = "Exacerbated_Down", UD = "Reversed_DisUp_ResDown", DU = "Reversed_DisDown_ResUp"),
  ora_min_size = 10,
  ora_quadrant_names = list(UU = "Exacerbated Up", DD = "Exacerbated Down", UD = "Reversed (Disease Up / Rescue Down)", DU = "Reversed (Disease Down / Rescue Up)"),
  ora_grouped = list(file_1_quads = c("ora_UD", "ora_DU"), file_2_quads = c("ora_UU", "ora_DD"),
                     note_if_empty_2 = "No pathways enriched in exacerbation quadrants (padj<0.05)"),
  ora_colors = ORA_QUAD_COLORS_F3,
  summary_quadrant_names = list(UU = "Exacerbated_Up", UU_slug = "exacerbated_up", DD = "Exacerbated_Down", DD_slug = "exacerbated_down",
                                UD = "Reversed_DisUp_ResDown", UD_slug = "reversed_dis_up", DU = "Reversed_DisDown_ResUp", DU_slug = "reversed_dis_down"),
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT, supp = NULL)
source(here::here("04_Figures", "shared", "comparison_panels", "panel_E_rrho2.R"))
n_shared_E <- n_shared; max_rev_E <- max(max_UD, max_DU); n_rev_E <- if (max_UD >= max_DU) n_UD else n_DU

# Panel B: pattern heatmap (LAST)
ROW_H <- 0.078
cfg <- list(
  fig_id = "F04", contrast_x = "CTLvPHE", contrast_y = "PHEvPHE_MITO",
  title = "Disease Reversal Patterns", col_headers = c("Dis.", "Re."), sort_col = "logFC_CTLvPHE",
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  classify_fn = function(dep_df) {
    dep_df |>
      dplyr::filter(!is.na(logFC_CTLvPHE), !is.na(logFC_PHEvPHE_MITO)) |>
      dplyr::filter(pi_score_CTLvPHE < H9C2_PI_THRESH | pi_score_PHEvPHE_MITO < H9C2_PI_THRESH) |>
      dplyr::mutate(
        quadrant = dplyr::case_when(
          logFC_CTLvPHE > 0 & logFC_PHEvPHE_MITO < 0 ~ "Reversed Up",
          logFC_CTLvPHE < 0 & logFC_PHEvPHE_MITO > 0 ~ "Reversed Down",
          TRUE ~ "Non-reversed"),
        sig_cat = dplyr::case_when(
          pi_score_CTLvPHE < H9C2_PI_THRESH & pi_score_PHEvPHE_MITO < H9C2_PI_THRESH ~ "Both",
          pi_score_CTLvPHE < H9C2_PI_THRESH ~ "Dis.",
          pi_score_PHEvPHE_MITO < H9C2_PI_THRESH ~ "Re.", TRUE ~ "NS")) },
  QUAD_ORDER = c("Reversed Up", "Reversed Down", "Non-reversed"),
  QUAD_COLORS = c("Reversed Up" = "#B2182B", "Reversed Down" = "#2166AC", "Non-reversed" = "#1B7837"),
  QUAD_BG = c("Reversed Up" = "#F4D9D2", "Reversed Down" = "#D5DEEF", "Non-reversed" = "#C8E0CD", "Tied" = "#EEEEEE"),
  ENDPOINT_COLORS = c("Reversed Up" = "#67001F", "Reversed Down" = "#053061", "Non-reversed" = "#00441B"),
  SIG_COLORS = c("Both" = "#2E7D32", "Dis." = "#E05A4E", "Re." = "#5DA5DA", "NS" = "grey70"),
  display_labels = c("Carbohydrate & Energy Metabolism" = "Carb. & Energy Metab.", "Amino Acid & Cofactor Metabolism" = "AA & Cofactor\nMetab."),
  col_header_colors = c(CONTRAST_COLORS["CTLvPHE"], CONTRAST_COLORS["PHEvPHE_MITO"]),
  bg_extend_right = 4.5, bar_scale = 0.20, bar_ref_width = 32,
  key_y_base = ROW_H * 15.5, key_dy = ROW_H * 3.8, key_x_sig = NULL,
  protein_count_x_mult = 15, count_tick_y_label = ROW_H * 2.6,
  count_tick_filter = function(df) dplyr::filter(df, val != 15),
  sig_cats = c("Res.", "Dis.", "Both"), sig_cat_labels = c("Sig Rescue", "Sig Disease", "Sig Both"),
  inset_legend = FALSE)
# Override panel B keys to F04 convention (Dis./Re./Both); SIG_COLORS at line 104
# already uses these keys, so no recolouring needed.
cfg$sig_cats <- c("Dis.", "Re.", "Both"); cfg$sig_cat_labels <- c("Sig Dis.", "Sig Re.", "Sig Both")
source(here::here("04_Figures", "shared", "comparison_panels", "panel_B_pattern_heatmap.R"))

RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")

nudge_idx3 <- 0.30
quad_legend <- ggplot(inset_quad_df) +
  geom_rect(aes(xmin = (as.integer(quadrant) - 1) * 3.5 + (as.integer(quadrant) == 3) * nudge_idx3,
                xmax = (as.integer(quadrant) - 1) * 3.5 + 0.7 + (as.integer(quadrant) == 3) * nudge_idx3, ymin = -0.35, ymax = 0.35),
            fill = inset_quad_df$bg_color, color = "black", linewidth = 0.5) +
  geom_rect(aes(xmin = (as.integer(quadrant) - 1) * 3.5 + 0.10 + (as.integer(quadrant) == 3) * nudge_idx3,
                xmax = (as.integer(quadrant) - 1) * 3.5 + 0.60 + (as.integer(quadrant) == 3) * nudge_idx3, ymin = -0.15, ymax = 0.15),
            fill = inset_quad_df$bar_color, color = "black", linewidth = 0.3) +
  geom_text(aes(x = (as.integer(quadrant) - 1) * 3.5 + 0.85 + (as.integer(quadrant) == 3) * nudge_idx3, y = 0, label = as.character(quadrant)),
            hjust = 0, size = 3.5, fontface = "bold", color = "grey15") +
  coord_cartesian(xlim = c(0, 11.5), ylim = c(-0.7, 0.7), clip = "off") + theme_void() +
  theme(plot.background = element_blank(), panel.background = element_blank(), plot.margin = margin(0, 0, 0, 0, "mm"))

# Composite layout: A full top row; B | C | D across the bottom
COMP_W <- 440; COMP_H <- 330
PRINT_SCALE_C <- 380 / 178
TAG_SZ <- round(10 * PRINT_SCALE_C * 0.85); TTL_SZ <- TAG_SZ; SUB_SZ <- round(7 * PRINT_SCALE_C * 0.85)

ttl_A <- "Protein Reversal (Quadrant ORA)"
sub_A <- sprintf("N = %d | %d DEPs (Π) | ρ = %.2f [%.2f, %.2f] | %.0f%% reversed | bold outline = mitochondrial protein", n_total_A, n_sig_A, r_spear_A, rho_lo_A, rho_hi_A, pct_rev_A)
ttl_B <- "Protein-to-Pathway"; sub_B <- sprintf("%d proteins | %d pathways", n_total, n_pw)
ttl_C <- "Pathway Reversal"; sub_C <- sprintf("Hallmark + GO-Slim | ρ = %.2f | %.0f%% reversed | orange ring = mito", rho_D, pw_rev_D * 100)
ttl_D <- "RRHO2 Reversal"; sub_D <- sprintf("%d genes | %s", n_shared_E, dens_reversed_cmp)

a_row   <- "AAAAAAAAAAAAA"
bcd_row <- "BBBBBCCCCDDDD"
gap_row <- "#############"
layout  <- paste(c(gap_row, rep(a_row, 8), gap_row, gap_row, rep(bcd_row, 8)), collapse = "\n")
heights <- c(0.55, rep(1.0, 8), 0.30, 1.15, rep(1.0, 8))

pD      <- pD + theme(plot.margin = margin(1, 2, 1, 2, "mm"))     # NES -> panel C
pE_heat <- pE_heat + theme(plot.margin = margin(1, 1, 1, 1, "mm"), axis.title = element_text(face = "bold", size = 8))
pB <- pB + coord_cartesian(xlim = c(-0.25, X_BAR_RIGHT + 0.3), ylim = c(BAR_YMAX + ROW_H * 6.5, -ROW_H * 0.05), expand = FALSE) +
           theme(plot.margin = margin(1, 1, 16, -1, "mm"))   # raised bottom -> aligns with C/D

fig <- wrap_elements(full = pA) + wrap_elements(full = pB) + wrap_elements(full = pD) + wrap_elements(full = pE_heat) +
       plot_layout(design = layout, widths = rep(1, 13), heights = heights)

X_TTL <- 0.020; TAG_DY <- -0.002; SUB_OFFSET <- 0.019
Y_TOP <- 0.985; Y_BOT <- 0.500
X_A <- 0.004; X_B <- 0.004; X_C <- 0.392; X_D <- 0.700

composite_final <- ggdraw(fig) +
  draw_label("A", x = X_A, y = Y_TOP - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(ttl_A, x = X_A + X_TTL, y = Y_TOP, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(sub_A, x = X_A + X_TTL, y = Y_TOP - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("B", x = X_B, y = Y_BOT - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(ttl_B, x = X_B + X_TTL, y = Y_BOT, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(sub_B, x = X_B + X_TTL, y = Y_BOT - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("C", x = X_C, y = Y_BOT - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(ttl_C, x = X_C + X_TTL, y = Y_BOT, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(sub_C, x = X_C + X_TTL, y = Y_BOT - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("D", x = X_D, y = Y_BOT - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(ttl_D, x = X_D + X_TTL, y = Y_BOT, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(sub_D, x = X_D + X_TTL, y = Y_BOT - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_plot(quad_legend, x = 0.03, y = 0.004, width = 0.34, height = 0.032) +
  draw_label("N=24 (n=6/group); Interaction underpowered",
             x = 0.99, y = 0.004, size = SUB_SZ, fontface = "italic",
             colour = "grey45", hjust = 1, vjust = 0)

ggsave(file.path(RPT_PDF, "MAIN_F04_composite.pdf"), composite_final, width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "MAIN_F04_composite.png"), composite_final, width = COMP_W, height = COMP_H, units = "mm", dpi = 300)
message("F04 composite (A top; B|C|D bottom) saved")
