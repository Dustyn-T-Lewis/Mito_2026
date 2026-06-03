#!/usr/bin/env Rscript
# F03 main composite: Transplant_Mito (CTLvMITO) vs Rescue_Mito+Phe (PHEvPHE_MITO).
# Question: does mito transplant engage the same biology in healthy vs diseased cells?
# 4-panel 2-row layout: A protein quadrant (NEW), B pattern heatmap (top);
# D fGSEA NES scatter, E RRHO2 (bottom). Panel B sources AnnotationDbi via
# go_slim_categories.R, so it runs LAST to avoid select() S4 masking.
#
# Ported almost verbatim from 04_Figures/F04/a_script/01_main_panels.R (already
# CTLvMITO vs PHEvPHE_MITO); the ORA composite panel A is replaced by the
# top-5-per-quadrant protein scatter, and the fry panel C is dropped.

setwd(rprojroot::find_rstudio_root_file())

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(ggplot2)
library(patchwork)
library(cowplot)

source("04_Figures/shared/style.R")

BASE    <- "05_Figures/F03_concordance"
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
PNL_PNG <- file.path(RPT_PNG, "panels")
PNL_PDF <- file.path(RPT_PDF, "panels")
DAT     <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, PNL_PNG, PNL_PDF, DAT))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

message("=== F03 Composite: sourcing panels ===")

# Panel A: top-5-per-quadrant protein concordance scatter (NEW)
source("05_Figures/F03_concordance/a_script/_panel_A_quadrant.R")

# Panel D: fGSEA NES scatter (config wrapper -> shared engine)
cfg <- list(
  fig_id      = "F03",
  contrast_x  = "CTLvMITO",
  contrast_y  = "PHEvPHE_MITO",
  fgsea_cache = "04_Figures/shared/fgsea_tstat_all_h9c2.csv",
  databases   = c("Hallmark", "GO Slim"),
  db_shapes   = c("Hallmark" = 24, "GO Slim" = 21),
  title       = "Pathway-Level Concordance (fGSEA)",
  axis_x_label = "Transplant_Mito (Mito - Ctl)",
  axis_y_label = "Rescue_Mito+Phe (PHE_Mito - PHE)",
  subtitle_metric = "concordant",
  subtitle_interpretation = "positive ρ = mito-transplant engages the same pathways in healthy and stressed cells",
  ref_slope   = 1,
  panel_w     = 146,
  label_border_size = 0.25,
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  sig_colors     = SIG_COLORS_F2,
  sig_label_fill = SIG_LABEL_FILL_F2,
  sig_label_text = SIG_LABEL_TEXT_F2,
  sig_draw_order = c("Sig Rescue only", "Sig Intervention only", "Sig Both", "Interaction"),
  quadrant_defs = list(
    sig_both_label = "Sig Both",
    sig_x_label    = "Sig Intervention only",
    sig_y_label    = "Sig Rescue only",
    bg_red_1  = c(0, Inf, 0, Inf),
    bg_red_2  = c(-Inf, 0, -Inf, 0),
    bg_blue_1 = c(0, Inf, -Inf, 0),
    bg_blue_2 = c(-Inf, 0, 0, Inf),
    label_tr = "Concordant Up",   color_tr = "#D6604D",
    label_tl = "Discordant",      color_tl = "#4393C3",
    label_bl = "Concordant Down", color_bl = "#D6604D",
    label_br = "Discordant",      color_br = "#4393C3",
    metric_count_fn = function(q1, q2, q3, q4) q1 + q3
  ),
  display_overrides = c(
    "Oxidative Phosphorylation"             = "OXPHOS",
    "Aerobic Respiration And Respiratory Electron Transport" = "Aerobic Resp. + ETC",
    "Citric Acid Cycle Tca Cycle"           = "TCA Cycle",
    "Mitochondrial Translation"             = "Mito Translation",
    "Mitochondrial Protein Degradation"     = "Mito Protein Deg.",
    "Mitochondrial Calcium Ion Transport"   = "Mito Ca²⁺ Transport",
    "Pink1 Prkn Mediated Mitophagy"         = "PINK1/PRKN Mitophagy",
    "Fatty Acid Beta Oxidation"             = "FA β-Oxidation",
    "Branched Chain Amino Acid Metabolism"  = "BCAA Metabolism",
    "Cardiolipin Biosynthesis"              = "Cardiolipin Biosynth.",
    "Ribosome Biogenesis"                   = "Ribo Bio",
    "Mitochondrial Organization"            = "Mito Org.",
    "Mitochondrial Transport"               = "Mito Transport",
    "Mitochondrial Protein Import"          = "Mito Protein Import",
    "Extracellular Matrix Organization"     = "ECM Org.",
    "Protein Folding"                       = "Protein Folding",
    "Amino Acid Metabolism"                 = "AA Metabolism",
    "Fatty Acid Metabolism"                 = "FA Metabolism",
    "Cytoplasmic Translation"               = "Cytoplasmic Transl."
  )
)
source("04_Figures/shared/comparison_panels/panel_D_nes_scatter.R")
n_pw_D     <- nrow(fgsea_wide)
n_sig_pw_D <- n_total_sig
rho_D      <- as.numeric(nes_cor_all$estimate)
pw_conc_D  <- pw_conc_frac

# Panel E: RRHO2 (config wrapper -> shared engine)
cfg <- list(
  fig_id     = "F03",
  t_col_1    = "t_CTLvMITO",
  t_col_2    = "t_PHEvPHE_MITO",
  rrho_labels = c("Transplant", "Rescue"),
  title       = "Threshold-Free Concordance (RRHO2)",
  subtitle_fmt = "Stratified hypergeometric | %d shared genes | warm corners = concordant gene regulation | No MTC (Cahill et al. 2018)",
  axis_label_1 = expression("Transplant rank"~(Up %->% Down)),
  axis_label_2 = expression("Rescue rank"~(Up %->% Down)),
  quadrant_labels = list(
    UU = "Concordant Up", DD = "Concordant Down",
    UD = "Discordant (Tr↑ Re↓)", DU = "Discordant (Tr↓ Re↑)"
  ),
  hotspot_export_names = list(UU = "UU", DD = "DD", UD = "UD", DU = "DU"),
  ora_min_size = 15,
  ora_quadrant_names = list(
    UU = "Concordant Up", DD = "Concordant Down",
    UD = "Discordant (Tr Up / Re Down)", DU = "Discordant (Tr Down / Re Up)"
  ),
  ora_grouped = list(
    file_1_quads = c("ora_UU", "ora_DD"),
    file_2_quads = c("ora_UD", "ora_DU")
  ),
  ora_colors = ORA_QUAD_COLORS_F2,
  summary_quadrant_names = list(
    UU = "Concordant_Up",       UU_slug = "concordant_up",
    DD = "Concordant_Down",     DD_slug = "concordant_down",
    UD = "Discordant_TrUp_ReDown",  UD_slug = "discordant_tr_up",
    DU = "Discordant_TrDown_ReUp",  DU_slug = "discordant_tr_down"
  ),
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  supp = NULL
)
source("04_Figures/shared/comparison_panels/panel_E_rrho2.R")
n_shared_E <- n_shared
n_UU_E     <- n_UU

# Panel B: pattern heatmap (LAST: loads AnnotationDbi, clobbers select)
ROW_H <- 0.078
cfg <- list(
  fig_id     = "F03",
  contrast_x = "CTLvMITO",
  contrast_y = "PHEvPHE_MITO",
  title      = "Transplant vs Rescue Response Patterns",
  col_headers = c("Tr.", "Re."),
  sort_col   = "logFC_CTLvMITO",
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  classify_fn = function(dep_df) {
    dep_df |>
      dplyr::filter(!is.na(logFC_CTLvMITO), !is.na(logFC_PHEvPHE_MITO)) |>
      dplyr::filter(pi_score_CTLvMITO < H9C2_PI_THRESH |
                    pi_score_PHEvPHE_MITO < H9C2_PI_THRESH |
                    pi_score_Interaction < H9C2_PI_THRESH) |>
      dplyr::mutate(
        quadrant = dplyr::case_when(
          logFC_CTLvMITO > 0 & logFC_PHEvPHE_MITO > 0 ~ "Concordant Up",
          logFC_CTLvMITO < 0 & logFC_PHEvPHE_MITO < 0 ~ "Concordant Down",
          TRUE ~ "Discordant"
        ),
        sig_cat = dplyr::case_when(
          pi_score_CTLvMITO < H9C2_PI_THRESH & pi_score_PHEvPHE_MITO < H9C2_PI_THRESH ~ "Both",
          pi_score_CTLvMITO < H9C2_PI_THRESH ~ "Tr.",
          pi_score_PHEvPHE_MITO < H9C2_PI_THRESH ~ "Re.",
          pi_score_Interaction < H9C2_PI_THRESH ~ "Inter.",
          TRUE ~ "NS"
        )
      )
  },
  QUAD_ORDER      = c("Concordant Up", "Concordant Down", "Discordant"),
  QUAD_COLORS     = c("Concordant Up" = "#B2182B", "Concordant Down" = "#2166AC",
                      "Discordant" = "#1B7837"),
  QUAD_BG         = c("Concordant Up" = "#F4D9D2", "Concordant Down" = "#D5DEEF",
                      "Discordant" = "#C8E0CD", "Tied" = "#EEEEEE"),
  ENDPOINT_COLORS = c("Concordant Up" = "#67001F", "Concordant Down" = "#053061",
                      "Discordant" = "#00441B"),
  SIG_COLORS      = c("Both" = "#2E7D32", "Tr." = "#E05A4E",
                      "Re." = "#5DA5DA", "Inter." = "#7B5EA7", "NS" = "grey70"),
  display_labels = c(
    "Carbohydrate & Energy Metabolism" = "Carb. & Energy Metab.",
    "Amino Acid & Cofactor Metabolism" = "AA & Cofactor\nMetab."
  ),
  col_header_colors = c(
    CONTRAST_COLORS["CTLvMITO"],
    CONTRAST_COLORS["PHEvPHE_MITO"]
  ),
  x_sig                = 0.27,
  tile_w               = 0.80,
  bar_scale            = 0.46,
  bar_ref_width        = 14,
  key_y_base           = ROW_H * 15.5,
  key_dy               = ROW_H * 3.8,
  key_x_sig            = NULL,
  protein_count_x_mult = 7.5,
  count_tick_y_label   = ROW_H * 2.6,
  count_ticks_max      = 16,
  count_tick_filter    = function(df) dplyr::filter(df, !val %in% c(15, 25)),
  sig_cats       = c("Tr.", "Re.", "Both", "Inter."),
  sig_cat_labels = c("Sig Tr.", "Sig Re.", "Sig Both", "Interaction"),
  inset_legend   = FALSE
)
source("04_Figures/shared/comparison_panels/panel_B_pattern_heatmap.R")

# Restore RPT paths (shared engines clobber RPT_PDF/RPT_PNG to panels subdir)
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")

# Quadrant legend (horizontal) — overlaid below panel B in the composite.
nudge_idx3 <- 0.30
quad_legend <- ggplot(inset_quad_df) +
  geom_rect(aes(xmin = (as.integer(quadrant) - 1) * 3.5 +
                       (as.integer(quadrant) == 3) * nudge_idx3,
                xmax = (as.integer(quadrant) - 1) * 3.5 + 0.7 +
                       (as.integer(quadrant) == 3) * nudge_idx3,
                ymin = -0.35, ymax = 0.35),
            fill = inset_quad_df$bg_color, color = "black", linewidth = 0.5) +
  geom_rect(aes(xmin = (as.integer(quadrant) - 1) * 3.5 + 0.10 +
                       (as.integer(quadrant) == 3) * nudge_idx3,
                xmax = (as.integer(quadrant) - 1) * 3.5 + 0.60 +
                       (as.integer(quadrant) == 3) * nudge_idx3,
                ymin = -0.15, ymax = 0.15),
            fill = inset_quad_df$bar_color, color = "black", linewidth = 0.3) +
  geom_text(aes(x = (as.integer(quadrant) - 1) * 3.5 + 0.85 +
                    (as.integer(quadrant) == 3) * nudge_idx3,
                y = 0, label = as.character(quadrant)),
            hjust = 0, size = 3.5, fontface = "bold", color = "grey15") +
  coord_cartesian(xlim = c(0, 10.5), ylim = c(-0.7, 0.7), clip = "off") +
  theme_void() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm"))

# Composite layout: A full top row; B | C | D across the bottom
# A = quadrant ORA composite (top, full width). Bottom: B protein-to-pathway
# (wider) | C pathway concordance (NES) | D RRHO2. C and D are equal-sized squares.
COMP_W <- 440
COMP_H <- 330
PRINT_SCALE_C <- 380 / 178
TAG_SZ <- round(10 * PRINT_SCALE_C * 0.85)
TTL_SZ <- TAG_SZ
SUB_SZ <- round(7 * PRINT_SCALE_C * 0.85)

ttl_A <- "Protein Concordance (Quadrant ORA)"
sub_A <- sprintf("N = %d | %d DEPs (Π) | ρ = %.2f [%.2f, %.2f] | %.0f%% concordant | bold outline = mitochondrial protein",
                 n_total_A, n_sig_A, r_spear_A, rho_lo_A, rho_hi_A, pct_conc_A)
ttl_B <- "Protein-to-Pathway"
sub_B <- sprintf("%d proteins | %d pathways", n_total, n_pw)
ttl_C <- "Pathway Concordance"
sub_C <- sprintf("Hallmark + GO-Slim | ρ = %.2f | %.0f%% concordant | orange ring = mito", rho_D, pw_conc_D * 100)
ttl_D <- "RRHO2 Concordance"
sub_D <- sprintf("%d genes | %s", n_shared_E, dens_concordant_cmp)

# 13-col grid: A spans full top; bottom row B(1-5) | C(6-9) | D(10-13).
a_row   <- "AAAAAAAAAAAAA"
bcd_row <- "BBBBBCCCCDDDD"
gap_row <- "#############"
layout  <- paste(c(gap_row, rep(a_row, 8), gap_row, gap_row, rep(bcd_row, 8)),
                 collapse = "\n")
heights <- c(0.55, rep(1.0, 8), 0.30, 1.15, rep(1.0, 8))

pD      <- pD + theme(plot.margin = margin(1, 2, 1, 2, "mm"))      # NES -> panel C
pE_heat <- pE_heat + theme(plot.margin = margin(1, 1, 1, 1, "mm"),
                           axis.title = element_text(face = "bold", size = 8))
pB <- pB + coord_cartesian(xlim = c(-0.25, X_BAR_RIGHT + 0.3),
                           ylim = c(BAR_YMAX + ROW_H * 6.5, -ROW_H * 0.05),
                           expand = FALSE) +
           theme(plot.margin = margin(1, 1, 16, -1, "mm"))   # raised bottom -> aligns with C/D

fig <- wrap_elements(full = pA) +      # region A (top, full width)
       wrap_elements(full = pB) +      # region B (bottom-left)
       wrap_elements(full = pD) +      # region C (bottom-mid, NES)
       wrap_elements(full = pE_heat) + # region D (bottom-right, RRHO2)
       plot_layout(design = layout, widths = rep(1, 13), heights = heights)

X_TTL <- 0.020; TAG_DY <- -0.002; SUB_OFFSET <- 0.019
Y_TOP <- 0.985; Y_BOT <- 0.500
X_A <- 0.004; X_B <- 0.004; X_C <- 0.392; X_D <- 0.700

composite_final <- ggdraw(fig) +
  draw_label("A",   x = X_A,         y = Y_TOP - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_A, x = X_A + X_TTL, y = Y_TOP,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_A, x = X_A + X_TTL, y = Y_TOP - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("B",   x = X_B,         y = Y_BOT - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_B, x = X_B + X_TTL, y = Y_BOT,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_B, x = X_B + X_TTL, y = Y_BOT - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("C",   x = X_C,         y = Y_BOT - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_C, x = X_C + X_TTL, y = Y_BOT,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_C, x = X_C + X_TTL, y = Y_BOT - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("D",   x = X_D,         y = Y_BOT - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_D, x = X_D + X_TTL, y = Y_BOT,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_D, x = X_D + X_TTL, y = Y_BOT - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_plot(quad_legend, x = 0.03, y = 0.004, width = 0.34, height = 0.032) +
  draw_label("N=24 (n=6/group); Interaction underpowered",
             x = 0.99, y = 0.004, size = SUB_SZ, fontface = "italic",
             colour = "grey45", hjust = 1, vjust = 0)

ggsave(file.path(RPT_PDF, "MAIN_F03_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "MAIN_F03_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F03 composite (A top; B|C|D bottom) saved")
