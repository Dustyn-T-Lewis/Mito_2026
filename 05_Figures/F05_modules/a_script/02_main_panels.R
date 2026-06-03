#!/usr/bin/env Rscript
# F05 (2/2) — headline composite for the knowledge-driven set-score modules.
# A: set x contrast heatmap — fill = GSVA emmeans contrast estimate, BH-FDR stars,
#    left strip = source (CORUM / MitoCarta / Reactome) + measured-member count.
# B: method robustness — GSVA vs singscore contrast estimates for the shown sets,
#    per-contrast Spearman (the headline does not depend on the scoring method).
# No circularity: sets are external; scores tested in the 2x2 LMM (see 01_set_scores.R).
# Reads c_data/* written by 01_set_scores.R.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(patchwork); library(cowplot)
})
source(here::here("05_Figures", "shared", "config.R"))
set.seed(42)

BASE    <- here::here("05_Figures", "F05_modules")
DAT     <- file.path(BASE, "c_data")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
for (d in c(RPT_PDF, RPT_PNG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

CORE   <- c("Transplant_Mito", "Disease_Phe", "Rescue_Mito+Phe", "Interaction")
COL_LAB <- c(Transplant_Mito = "Transplant\n(Mito − Ctl)",
             Disease_Phe     = "Disease\n(PHE − Ctl)",
             `Rescue_Mito+Phe` = "Rescue\n(PHE_Mito − PHE)",
             Interaction     = "Interaction\n(Rescue − Transplant)")
CONTRAST_PAL <- c(Transplant_Mito = unname(CONTRAST_COLORS["CTLvMITO"]),
                  Disease_Phe     = unname(CONTRAST_COLORS["CTLvPHE"]),
                  `Rescue_Mito+Phe` = unname(CONTRAST_COLORS["PHEvPHE_MITO"]),
                  Interaction     = unname(CONTRAST_COLORS["Interaction"]))
SRC_LAB <- c(CORUM = "CORUM", MITO = "MitoCarta", REAC = "Reactome")
SRC_PAL <- c(CORUM = "#8856A7", MITO = "#D55E00", REAC = "#1565C0")

gsva_lmm <- read_csv(file.path(DAT, "set_lmm_gsva.csv"), show_col_types = FALSE) |>
  filter(contrast %in% CORE) |>
  mutate(contrast = factor(contrast, levels = CORE))
ss_lmm <- read_csv(file.path(DAT, "set_lmm_singscore.csv"), show_col_types = FALSE) |>
  filter(contrast %in% CORE) |> mutate(contrast = factor(contrast, levels = CORE))
concord <- read_csv(file.path(DAT, "set_method_concordance.csv"), show_col_types = FALSE)

# Select displayed sets: FDR<0.05 in >=1 core contrast, top by peak signal, capped
# for legibility; ordered by source then signal.
N_SHOW <- 30L
sig_sets <- gsva_lmm |> group_by(set_id) |>
  summarise(source = first(source), display = first(display),
            n_members = first(n_members),
            any_sig = any(p_fdr < 0.05),
            peak = max(abs(estimate) * -log10(p_value)), .groups = "drop") |>
  filter(any_sig) |> arrange(desc(peak)) |> slice_head(n = N_SHOW)

row_ord <- sig_sets |> arrange(source, peak) |>
  mutate(row_lab = ifelse(nchar(display) > 42, paste0(substr(display, 1, 40), "…"), display),
         row_lab = make.unique(row_lab, sep = " "))
ROW_LV <- row_ord$row_lab
row_ord <- row_ord |> mutate(row_lab = factor(row_lab, levels = ROW_LV))

heat <- gsva_lmm |> filter(set_id %in% row_ord$set_id) |>
  left_join(row_ord |> select(set_id, row_lab), by = "set_id") |>
  mutate(row_lab = factor(as.character(row_lab), levels = ROW_LV),
         stars = sig_stars(p_fdr))
ax_max <- max(abs(heat$estimate), na.rm = TRUE)

pA <- ggplot(heat, aes(contrast, row_lab, fill = estimate)) +
  geom_tile(color = "grey88", linewidth = 0.25) +
  geom_tile(data = filter(heat, p_fdr < 0.05),
            color = "black", linewidth = 0.5, fill = NA) +
  geom_text(aes(label = stars), size = 2.0, fontface = "bold",
            color = "grey10", vjust = 0.75) +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-ax_max, ax_max),
                       oob = scales::squish, name = "GSVA contrast\nestimate") +
  scale_x_discrete(labels = COL_LAB[CORE], position = "top") +
  labs(x = NULL, y = NULL,
       title = "Knowledge-driven set-score modules — exploratory, N = 24",
       subtitle = "GSVA score per sample → lmer(score ~ PHE×Mito + (1|Replicate)) → emmeans contrast + BH-FDR | external sets (no circularity) | * FDR<0.05  ** <0.01  *** <0.001") +
  FIG_THEME +
  theme(axis.text.x.top = element_text(face = "bold", size = 5, lineheight = 0.85),
        axis.text.y     = element_text(size = 4.2, lineheight = 0.8),
        panel.grid      = element_blank(), panel.border = element_blank(),
        legend.position = "right", legend.key.width = unit(2.5, "mm"),
        legend.key.height = unit(5, "mm"))

src_strip <- ggplot(row_ord, aes(x = 1, y = row_lab)) +
  geom_tile(aes(fill = source), width = 1) +
  geom_text(aes(label = n_members), size = 2.0, color = "white", fontface = "bold") +
  scale_fill_manual(values = SRC_PAL, labels = SRC_LAB, name = "Set source") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, subtitle = "src | n") +
  FIG_THEME +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        plot.subtitle = element_text(size = 3.5, hjust = 0.5),
        legend.position = "right", legend.key.size = unit(2.5, "mm"))

# Panel B — GSVA vs singscore effect-size concordance for the displayed sets.
both <- inner_join(
  gsva_lmm |> filter(set_id %in% row_ord$set_id) |> select(set_id, contrast, gsva = estimate),
  ss_lmm   |> filter(set_id %in% row_ord$set_id) |> select(set_id, contrast, ss = estimate),
  by = c("set_id", "contrast"))
rho_lab <- concord |> filter(contrast %in% CORE) |>
  mutate(contrast = factor(contrast, levels = CORE),
         lab = sprintf("%s: ρ = %.2f", contrast, spearman))

pB <- ggplot(both, aes(gsva, ss, color = contrast)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.2) +
  geom_point(size = 1.0, alpha = 0.8) +
  scale_color_manual(values = CONTRAST_PAL, name = NULL) +
  labs(x = "GSVA emmeans estimate", y = "singscore emmeans estimate",
       title = "Method robustness (GSVA vs singscore)",
       subtitle = sprintf("Spearman ρ = %.2f–%.2f across contrasts (displayed sets)",
                          min(concord$spearman), max(concord$spearman))) +
  FIG_THEME +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.7), color = NA),
        legend.key.size = unit(2.5, "mm"))

heat_row <- (src_strip | pA) + plot_layout(widths = c(0.06, 1))

# Honesty + interpretation block for the right column (no over-claiming at N=24).
CAPTION <- paste(
  "Sets are defined externally (CORUM + MitoCarta + Reactome), scored per sample,",
  "then tested — modules are NOT learned from these 24 samples, so there is no",
  "double-dipping (contrast with the WGCNA supplement).",
  "",
  "Coordinated positive shifts concentrate in Transplant and Rescue; Disease and",
  "the orthogonal Interaction show no set-level signal.",
  "",
  "Caveats: N=24 (n=6/group), single timepoint; MaxLFQ is cross-sample relative",
  "(no absolute stoichiometry); transplanted mitochondria are skeletal-muscle, so",
  "mito-set increases may reflect donor admixture, not only recipient remodeling.",
  sep = "\n")

COMP_W <- 200; COMP_H <- 168
composite <- ggdraw() +
  draw_plot(heat_row, x = 0,    y = 0,    width = 0.66, height = 1) +
  draw_plot(pB,       x = 0.685, y = 0.40, width = 0.31, height = 0.42) +
  draw_label("A", x = 0.005, y = 0.992, size = 11, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = 0.685, y = 0.86, size = 11, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(CAPTION, x = 0.69, y = 0.34, size = 5.4, hjust = 0, vjust = 1,
             colour = "grey25", lineheight = 0.95)

ggsave(file.path(RPT_PDF, "MAIN_F05_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_F05_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300, limitsize = FALSE)

write_csv(heat |> select(set_id, source, display, contrast, estimate, se, p_value, p_fdr,
                         replicate_icc, singular, n_members),
          file.path(DAT, "01_panel_A_setscore_heatmap_data.csv"))
message(sprintf("F05 set-score main composite done | %d sets shown | %s",
                nrow(row_ord),
                paste(sprintf("%s ρ=%.2f", concord$contrast, concord$spearman), collapse = " ")))
