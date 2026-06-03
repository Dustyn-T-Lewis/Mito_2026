#!/usr/bin/env Rscript
# F06 (3/4) figure — mitochondrial complex & mito-pathway remodeling (design-aware).
# Field-standard, parsimonious encodings for LFQ expression proteomics:
#   A  subunit-resolved OXPHOS heatmap (headline) — rows = ETC subunits grouped by
#      complex (CI/CIII/CIV/CV), cols = contrasts, fill = logFC, BH-FDR stars.
#      (García-Poyatos 2020, EMBO Rep, doi:10.15252/embr.202050287 — per-complex
#      grouped-subunit abundance is the closest reproducible analogue to complexome
#      profiling, which needs BN-PAGE fractionation we do NOT have.)
#   B  mitochondrial content (mass) — mean log2 of all MitoCarta proteins / sample.
#   C  mitonuclear OXPHOS balance — log2(mtDNA / nuclear), mass-independent.
#   D  mito-pathway dot plot — limma::camera (correlation-aware; Wu & Smyth 2012,
#      PMID 22638577) over the MitoCarta hierarchy; size = −log10 FDR, colour = dir.
# B/C significance = design-aware LMM emmeans-contrast FDRs (02_lmm_dynamics.R) — the
# single canonical framework. A/D significance = BH-FDR. Caveats: N=24 (n=6/group),
# single timepoint, MaxLFQ relative-quant (no absolute stoichiometry), donor mito are
# skeletal-muscle (admixture confound). Houtkooper 2013 (PMID 23698443); MitoCarta 3.0
# Rath 2021 (PMID 33174596).

setwd(rprojroot::find_rstudio_root_file())
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(patchwork); library(cowplot)
})
source("04_Figures/shared/style.R")
set.seed(42)

DAT     <- "05_Figures/F06_complex_mito/c_data"
RPT_PDF <- "05_Figures/F06_complex_mito/b_reports/main/pdf"
RPT_PNG <- "05_Figures/F06_complex_mito/b_reports/main/png"
for (d in c(RPT_PDF, RPT_PNG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

GROUP_LV    <- c("Ctl", "Mito", "PHE", "PHE_Mito")
ETC_C       <- c("Complex I", "Complex III", "Complex IV", "Complex V")
CONTRAST_LV <- c("Transplant_Mito", "Disease_Phe", "Rescue_Mito+Phe")
CONTRAST_AB <- c(Transplant_Mito = "Tx", Disease_Phe = "Dz", `Rescue_Mito+Phe` = "Rx")
contrast_old <- c(Transplant_Mito = "CTLvMITO", Disease_Phe = "CTLvPHE", `Rescue_Mito+Phe` = "PHEvPHE_MITO")
COL_LAB <- c(Transplant_Mito = "Transplant\n(Mito−Ctl)", Disease_Phe = "Disease\n(PHE−Ctl)",
             `Rescue_Mito+Phe` = "Rescue\n(PHE_Mito−PHE)")
mark <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*",
                    ifelse(p < 0.10, ".", ""))))

comb    <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
meta    <- readRDS(file.path(DAT, "analysis_meta.rds"))
content <- read_csv(file.path(DAT, "mito_content.csv"), show_col_types = FALSE) |> mutate(Group = factor(Group, levels = GROUP_LV))
balance <- read_csv(file.path(DAT, "mitonuclear_balance.csv"), show_col_types = FALSE) |> mutate(Group = factor(Group, levels = GROUP_LV))
mc_lmm  <- read_csv(file.path(DAT, "mito_content_lmm.csv"), show_col_types = FALSE)
mn_lmm  <- read_csv(file.path(DAT, "mitonuclear_lmm.csv"), show_col_types = FALSE)
mp_cam  <- read_csv(file.path(DAT, "mito_pathway_camera.csv"), show_col_types = FALSE)

# B / C box panels (LMM-sourced significance)
LMM_BRACKET <- list(Transplant_Mito = c("Ctl", "Mito"), Disease_Phe = c("Ctl", "PHE"),
                    `Rescue_Mito+Phe` = c("PHE", "PHE_Mito"))
br1 <- function(g1, g2, y, lab, tip) list(
  annotate("segment", x = g1, xend = g2, y = y, yend = y, linewidth = 0.3),
  annotate("segment", x = g1, xend = g1, y = y, yend = y - tip, linewidth = 0.3),
  annotate("segment", x = g2, xend = g2, y = y, yend = y - tip, linewidth = 0.3),
  annotate("text", x = (g1 + g2) / 2, y = y + tip * 0.4, label = lab,
           size = scale_text(BASE_STAT, 85), fontface = "bold", vjust = 0))
lmm_brackets <- function(lmm, y0, step, tip) {
  out <- list(); y <- y0
  for (cn in names(LMM_BRACKET)) {
    row <- lmm[lmm$contrast == cn, ]; if (!nrow(row)) next
    x1 <- match(LMM_BRACKET[[cn]][1], GROUP_LV); x2 <- match(LMM_BRACKET[[cn]][2], GROUP_LV)
    out <- c(out, br1(min(x1, x2), max(x1, x2), y, mark(row$p_fdr), tip)); y <- y + step }
  out }
lmm_corner <- function(lmm, x, y) {
  it <- lmm[lmm$contrast == "Interaction", ]
  annotate("label", x = x, y = y, hjust = 1, vjust = 1,
           label = sprintf("LMM y ~ PHE×Mito + (1|Rep)\ninteraction FDR %s | Rep ICC %.2f%s",
                           fmt_p(it$p_fdr), it$replicate_icc[1],
                           if (isTRUE(it$singular[1])) " (sing.)" else ""),
           size = scale_text(BASE_STAT, 85) - 0.2, lineheight = 0.95,
           fill = scales::alpha("white", 0.85), label.padding = unit(1.2, "pt")) }
box_panel <- function(df, yv, lmm, ylab, title, sub) {
  r <- range(df[[yv]]); s <- diff(r)
  ggplot(df, aes(Group, .data[[yv]])) +
    geom_boxplot(aes(fill = Group), width = 0.6, outlier.shape = NA, alpha = 0.55, linewidth = 0.3, colour = "grey25") +
    geom_jitter(aes(fill = Group), width = 0.14, size = 1.1, shape = 21, stroke = 0.25, colour = "grey25", alpha = 0.9) +
    scale_fill_manual(values = H9C2_PAL_GROUP, guide = "none") +
    lmm_brackets(lmm, r[2] + s * 0.06, s * 0.11, s * 0.02) + lmm_corner(lmm, 4.6, r[1] - s * 0.02) +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.42))) +
    labs(title = title, subtitle = sub, x = NULL, y = ylab) +
    FIG_THEME + theme(axis.text.x = element_text(face = "bold", size = FIG_AXIS_TEXT)) }

pB <- box_panel(content, "mito_content", mc_lmm, "mean log₂ (mito proteome)",
                "Mitochondrial content (mass)",
                "mean log₂ of all MitoCarta proteins | transplant delivers mitochondria")
pC <- box_panel(balance, "mitonuclear_log2ratio", mn_lmm, expression(log[2]*"(mtDNA / nuclear)"),
                "Mitonuclear OXPHOS balance",
                sprintf("mass-independent | %d mtDNA / %d nuclear subunits | donor (muscle) admixture caveat",
                        length(meta$mtdna_genes), length(meta$nuclear_oxphos)))

# A subunit-resolved ETC heatmap (headline)
sub_long <- bind_rows(lapply(ETC_C, function(cx) {
  g <- meta$complex_members[[cx]]
  bind_rows(lapply(names(contrast_old), function(cn) comb |> filter(gene %in% g) |>
    transmute(complex = cx, gene, contrast = cn,
              logFC = .data[[paste0("logFC_", contrast_old[cn])]],
              fdr   = .data[[paste0("adj.P.Val_", contrast_old[cn])]])))})) |>
  mutate(complex = factor(complex, levels = ETC_C),
         contrast = factor(contrast, levels = CONTRAST_LV), stars = sig_stars(fdr))
ord <- sub_long |> filter(contrast == "Transplant_Mito") |> arrange(complex, logFC) |>
  pull(gene) |> unique()
sub_long <- sub_long |> mutate(gene = factor(gene, levels = ord))
ax_a <- max(abs(sub_long$logFC), na.rm = TRUE)

pA <- ggplot(sub_long, aes(contrast, gene, fill = logFC)) +
  geom_tile(color = "grey90", linewidth = 0.2) +
  geom_text(aes(label = stars), size = 1.7, fontface = "bold", colour = "grey10", vjust = 0.75) +
  facet_grid(complex ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D", midpoint = 0,
                       limits = c(-ax_a, ax_a), oob = scales::squish, name = "subunit log₂FC") +
  scale_x_discrete(labels = COL_LAB[CONTRAST_LV], position = "top") +
  labs(x = NULL, y = NULL, title = "Subunit-resolved OXPHOS abundance",
       subtitle = "per-subunit log₂FC grouped by complex | * FDR<0.05 ** <0.01 *** <0.001 | abundance, not subunit-ratio") +
  FIG_THEME +
  theme(axis.text.x.top = element_text(face = "bold", size = 4.5, lineheight = 0.8),
        axis.text.y     = element_text(size = 3.6, face = "italic"),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 4.2),
        strip.placement = "outside", panel.grid = element_blank(), panel.border = element_blank(),
        panel.spacing.y = unit(0.6, "mm"),
        legend.position = "bottom", legend.key.height = unit(2, "mm"), legend.key.width = unit(6, "mm"))

# D mito-pathway dot plot (camera over MitoCarta hierarchy)
mp_label <- function(x) {
  seg <- str_split(str_remove(x, "^MITOCARTA_"), "__")
  lab <- vapply(seg, function(s) tail(s, 1), character(1))
  lab <- str_replace_all(lab, "_", " ") |> str_to_title()
  str_replace_all(lab, c("Oxphos" = "OXPHOS", "Dna" = "DNA", "Rna" = "RNA", "Tca" = "TCA",
                         "Micos" = "MICOS", "Ros " = "ROS ")) |> str_squish()
}
mp <- mp_cam |> mutate(contrast = factor(contrast, levels = CONTRAST_LV),
                       label = mp_label(set), nlfdr = -log10(fdr),
                       dir = factor(direction, levels = c("Up", "Down")))
top_sets <- mp |> group_by(set) |> summarise(best = min(fdr), .groups = "drop") |>
  filter(best < 0.05) |> arrange(best) |> slice_head(n = 18) |> pull(set)
lab_order <- mp |> filter(set %in% top_sets) |> group_by(label) |>
  summarise(b = min(fdr), .groups = "drop") |> arrange(desc(b)) |> pull(label)
mp_top <- mp |> filter(set %in% top_sets) |> mutate(label = factor(label, levels = lab_order))

pD <- ggplot(mp_top, aes(contrast, label)) +
  geom_point(aes(size = nlfdr, fill = dir, colour = fdr < 0.05), shape = 21, stroke = 0.5) +
  scale_size_area(max_size = 5, name = expression(-log[10]~FDR), breaks = c(1, 2, 4)) +
  scale_fill_manual(values = c(Up = "#D6604D", Down = "#4393C3"), name = "Direction") +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey75"), guide = "none") +
  scale_x_discrete(labels = CONTRAST_AB[CONTRAST_LV], position = "top") +
  labs(x = NULL, y = NULL, title = "MitoCarta-pathway shifts (camera)",
       subtitle = "top pathways (FDR<0.05 in ≥1 contrast) | size = −log10 FDR | black ring FDR<0.05") +
  FIG_THEME +
  theme(axis.text.x.top = element_text(face = "bold", size = 5),
        axis.text.y     = element_text(size = 4.2),
        panel.grid.major = element_line(color = "grey93", linewidth = 0.2),
        legend.position = "right", legend.box = "vertical",
        legend.key.size = unit(2.5, "mm"), legend.spacing.y = unit(0.5, "mm"))

# Composite: A tall on left (headline) | B / C / D stacked on right
COMP_W <- 188; COMP_H <- 220
composite <- ggdraw() +
  draw_plot(pA, x = 0,     y = 0,    width = 0.46, height = 1) +
  draw_plot(pB, x = 0.475, y = 0.70, width = 0.525, height = 0.29) +
  draw_plot(pC, x = 0.475, y = 0.41, width = 0.525, height = 0.29) +
  draw_plot(pD, x = 0.475, y = 0,    width = 0.525, height = 0.40) +
  draw_label("A", x = 0.005, y = 0.995, size = 11, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = 0.480, y = 0.995, size = 11, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C", x = 0.480, y = 0.705, size = 11, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("D", x = 0.480, y = 0.415, size = 11, fontface = "bold", hjust = 0, vjust = 1)

ggsave(file.path(RPT_PDF, "MAIN_F06_complex_mito.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_F06_complex_mito.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300, limitsize = FALSE)
message("F06 figure: MAIN (subunit heatmap + content + balance + mito-pathway dots) done")
