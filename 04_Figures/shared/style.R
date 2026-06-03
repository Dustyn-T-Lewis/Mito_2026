# style.R — palettes, themes, sizing helpers for figure scripts.
# Single source of truth for palettes + thresholds: 00_input/h9c2_design.R.

library(ggplot2)
library(scales)
library(grid)

# Single source of truth for groups, contrasts, palettes, thresholds.
source(here::here("00_input", "h9c2_design.R"))

# suppress stray Rplots.pdf from implicit device opens
options(device = function(...) grDevices::pdf(file = nullfile(), ...))

# Palettes (from 00_input/h9c2_design.R)
GROUP_COLORS    <- H9C2_PAL_GROUP            # Ctl / Mito / PHE / PHE_Mito
GROUP_FILL      <- H9C2_PAL_GROUP
PCA_COLORS      <- H9C2_PAL_GROUP
DIR_COLORS      <- H9C2_PAL_DIR              # Up / Down / NS
CONTRAST_COLORS <- H9C2_PAL_CONTRAST          # clinical-intuition palette
CONTRAST_ROLES  <- H9C2_CONTRAST_ROLES
PAL_CLASS <- c(Complete = "#4DAF4A", MAR = "#377EB8", MNAR = "#E41A1C")

DB_COLORS <- c(Hallmark = "#AA336A", KEGG = "#E65100", Reactome = "#1565C0",
               "GO:BP" = "#00796B", "GO Slim" = "#8D6E63",
               BioCarta = "#795548", PID = "#455A64")

# Subcellular-compartment ring palette (Okabe-Ito; colourblind-safe). Used as the
# point OUTLINE on the F03/F04 protein scatters (fill = significance class) and
# the mitochondrial outline on the NES scatter. Lookup: protein_localization_rat.csv.
LOC_COLORS <- c(Mitochondrial = "#D55E00", Nuclear = "#332288",
                Cytosolic = "#117733", Other = "grey70")

# Contrast ordering / short labels
CONTRAST_ORDER <- c(H9C2_CORE_CONTRASTS, "MITOvPHE_MITO")
CTR_SHORT <- c(
  CTLvPHE       = "PHE",
  CTLvMITO      = "Mito",
  PHEvPHE_MITO  = "PHE+Mito",
  Interaction   = "Inter.",
  MITOvPHE_MITO = "PHE|Mito"
)
CTR_FACET <- CTR_SHORT

# Sizing (J Physiol double-column spec — design-agnostic, kept from YvO)
PANEL_MD      <- 178
BASE_PATHWAY  <- 2.8
BASE_GENE     <- 2.5
BASE_STAT     <- 2.5
BASE_QUADRANT <- 2.8
BASE_COUNT    <- 2.5
BASE_TAG      <- 8

scale_text <- function(base_size, panel_width_mm, ref_width = PANEL_MD) {
  base_size * sqrt(panel_width_mm / ref_width)
}

strip_for_composite <- function(p) {
  p + labs(title = NULL, subtitle = NULL, tag = NULL) +
    theme(legend.position = "none")
}

is_light_color <- function(color_name) {
  rgb_val <- col2rgb(color_name)
  (0.299 * rgb_val[1] + 0.587 * rgb_val[2] + 0.114 * rgb_val[3]) / 255 > 0.6
}

# Text hierarchy + theme
FIG_TITLE_SIZE    <- 7
FIG_SUBTITLE_SIZE <- 4
FIG_STRIP_SIZE    <- 5
FIG_AXIS_TEXT     <- 5
FIG_LEGEND_TITLE  <- 5
FIG_LEGEND_TEXT   <- 4

composite_text_sizes <- function(comp_h_mm) {
  list(
    title    = pmax(6, pmin(8, round(5 + comp_h_mm / 80))),
    subtitle = pmax(4, pmin(6, round(3 + comp_h_mm / 100))),
    tag      = 8
  )
}

FIG_THEME <- theme_bw(base_size = 6, base_family = "Helvetica") +
  theme(
    plot.title       = element_text(face = "bold", size = FIG_TITLE_SIZE,
                                    margin = margin(b = 1)),
    plot.subtitle    = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                    color = "grey30", margin = margin(t = 0, b = 2)),
    plot.tag         = element_text(face = "bold", size = BASE_TAG),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = FIG_STRIP_SIZE),
    axis.title.x     = element_text(face = "bold", size = 5, margin = margin(t = 0)),
    axis.title.y     = element_text(face = "bold", size = 5, margin = margin(r = -1)),
    axis.text        = element_text(size = FIG_AXIS_TEXT, color = "grey15"),
    legend.title     = element_text(face = "bold", size = FIG_LEGEND_TITLE,
                                    color = "grey20"),
    legend.text      = element_text(size = FIG_LEGEND_TEXT, color = "grey15"),
    legend.key.size  = unit(2.5, "mm"),
    panel.grid.minor = element_blank()
  )

# F04/F05 concordance helpers (ported from YvO style.R)
# F04 = Intervention (CTLvMITO) vs Rescue (PHEvPHE_MITO) — concordance
# F05 = Disease (CTLvPHE) vs Rescue (PHEvPHE_MITO) — reversal
# F2-style sig classification: 4 categories + NS, used in concordance/reversal
# scatter plots and ORA quadrant analyses.

SIG_COLORS_F2 <- c(
  "Interaction"                  = "#7B5EA7",
  "Sig Both"                     = "#2E7D32",
  "Sig Intervention only"        = "#E05A4E",
  "Sig Rescue only"              = "#5DA5DA",
  "NS"                           = "grey70"
)

SIG_LABEL_FILL_F2 <- c(
  "Interaction"                  = scales::alpha("#7B5EA7", 0.75),
  "Sig Both"                     = scales::alpha("#2E7D32", 0.75),
  "Sig Intervention only"        = scales::alpha("#E05A4E", 0.75),
  "Sig Rescue only"              = scales::alpha("#5DA5DA", 0.75),
  "NS"                           = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_F2 <- setNames(rep("white", 5), names(SIG_LABEL_FILL_F2))

# F05-style sig classification (reversal — Disease vs Rescue)
SIG_COLORS_F3 <- c(
  "Sig Both"             = "#2E7D32",
  "Sig Disease only"     = "#E05A4E",
  "Sig Rescue only"      = "#5DA5DA",
  "NS"                   = "grey70"
)
SIG_LABEL_FILL_F3 <- c(
  "Sig Both"             = scales::alpha("#2E7D32", 0.75),
  "Sig Disease only"     = scales::alpha("#E05A4E", 0.75),
  "Sig Rescue only"      = scales::alpha("#5DA5DA", 0.75),
  "NS"                   = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_F3 <- setNames(rep("white", 4), names(SIG_LABEL_FILL_F3))

# Quadrant colors for F04/F05 ORA panels (4 sig quadrants in concordance scatters)
ORA_QUAD_COLORS_F2 <- c(
  "Concordant Up"                = "#E57373",
  "Concordant Down"              = "#64B5F6",
  "Discordant (Int Up / Res Dn)" = "#FFB74D",
  "Discordant (Int Dn / Res Up)" = "#81C784"
)
ORA_QUAD_COLORS_F3 <- c(
  "Reversed (Disease Up / Rescue Down)" = "#E57373",
  "Reversed (Disease Down / Rescue Up)" = "#64B5F6",
  "Exacerbated Up"                      = "#FFB74D",
  "Exacerbated Down"                    = "#81C784"
)

# F2 classify: F04 concordance. pi_Int = Intervention (CTLvMITO),
# pi_Res = Rescue (PHEvPHE_MITO), pi_int = Interaction (H9C2 Interaction contrast).
classify_proteins_f2 <- function(pi_Int, pi_Res, pi_int, threshold = H9C2_PI_THRESH) {
  dplyr::case_when(
    !is.na(pi_int) & pi_int < threshold                              ~ "Interaction",
    !is.na(pi_Int) & !is.na(pi_Res) &
      pi_Int < threshold & pi_Res < threshold                        ~ "Sig Both",
    !is.na(pi_Int) & pi_Int < threshold                              ~ "Sig Intervention only",
    !is.na(pi_Res) & pi_Res < threshold                              ~ "Sig Rescue only",
    TRUE                                                             ~ "NS"
  ) |>
    factor(levels = c("Interaction", "Sig Both",
                      "Sig Intervention only", "Sig Rescue only", "NS"))
}

# F3 classify: F05 reversal. pi_Dis = Disease (CTLvPHE), pi_Res = Rescue.
classify_proteins_f3 <- function(pi_Dis, pi_Res, threshold = H9C2_PI_THRESH) {
  dplyr::case_when(
    !is.na(pi_Dis) & !is.na(pi_Res) &
      pi_Dis < threshold & pi_Res < threshold                        ~ "Sig Both",
    !is.na(pi_Dis) & pi_Dis < threshold                              ~ "Sig Disease only",
    !is.na(pi_Res) & pi_Res < threshold                              ~ "Sig Rescue only",
    TRUE                                                             ~ "NS"
  ) |>
    factor(levels = c("Sig Both", "Sig Disease only",
                      "Sig Rescue only", "NS"))
}

# Stats / formatting helpers (design-agnostic, kept from YvO)
fmt_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01)  return(sprintf("p = %.3f", p))
  sprintf("p = %.2f", p)
}

fmt_p_plot <- function(p, threshold = 0.05) {
  label <- fmt_p(p)
  if (p < threshold) paste0('bold("', label, '")') else paste0('"', label, '"')
}

# Bonett & Wright 2000 — Fisher z CI for r (k = number of covariates)
fisher_z_ci <- function(r, n, k = 0, level = 0.95) {
  n_eff <- n - k
  if (n_eff < 4 || is.na(r)) return(c(lo = NA_real_, hi = NA_real_))
  z <- atanh(r); se <- 1 / sqrt(n_eff - 3)
  crit <- qnorm(1 - (1 - level) / 2)
  c(lo = tanh(z - crit * se), hi = tanh(z + crit * se))
}

sig_stars <- function(padj) {
  dplyr::case_when(
    padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ "")
}

boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}

# Pathway-name cleaner (used by enrichment figures F02/F04/F05 when built).
.DB_PREFIXES <- c("^HALLMARK_", "^GOSLIM_", "^GOBP_", "^GOCC_", "^GOMF_",
                  "^REACTOME_", "^KEGG_MEDICUS_", "^KEGG_")
.SCI_CAPS <- c(
  "Mtorc1" = "mTORC1", "Myc " = "MYC ", "E2f " = "E2F ", "Dna " = "DNA ",
  "Rna " = "RNA ", "Tnfa " = "TNFa ", "Uv " = "UV ", "G2m " = "G2M ",
  "Il6 " = "IL6 ", "Il2 " = "IL2 ", "Kras " = "KRAS ", "P53 " = "p53 ",
  "Tgf " = "TGF ", "Nfkb" = "NF-kB", "Atp " = "ATP ", "Nadh " = "NADH ",
  "Oxidative Phosphorylation"                  = "OXPHOS",
  "External Encapsulating Structure Or.*"      = "Extracellular Matrix Organization",
  "Enzyme Linked Receptor Protein Signaling.*" = "Receptor Protein Signaling")

clean_pathway_name <- function(name) {
  out <- name
  for (pfx in .DB_PREFIXES) out <- stringr::str_remove(out, pfx)
  out <- stringr::str_replace_all(out, "_", " ")
  out <- stringr::str_to_title(out)
  for (i in seq_along(.SCI_CAPS))
    out <- stringr::str_replace(out, names(.SCI_CAPS)[i], .SCI_CAPS[i])
  out
}

make_sigmoid_ribbon <- function(x0, x1, y0_top, y0_bot, y1_top, y1_bot,
                                n_pts = 50, ribbon_id) {
  t <- seq(0, 1, length.out = n_pts)
  blend <- (1 - cos(pi * t)) / 2
  tibble::tibble(
    x = c(x0 + (x1 - x0) * t, rev(x0 + (x1 - x0) * t)),
    y = c(y0_top + (y1_top - y0_top) * blend,
          rev(y0_bot + (y1_bot - y0_bot) * blend)),
    ribbon_id = ribbon_id)
}

# YvO-style composite finalisation
# Match the polish of A_YvO_2026/04_Figures/F03/a_script/01_main_panels.R:
# strip per-panel labels from each ggplot, combine via patchwork, wrap in
# ggdraw, then draw_label each panel's tag/title/subtitle and a composite-
# level header at explicit normalized coordinates.
#
# Usage:
#   composite <- (sl(pA) | sl(pB)) / (sl(pC) | sl(pD))
#   out <- finalize_composite_yvo(
#     composite, comp_h_mm = 180,
#     panels = list(
#       list(tag = "A", title = "PHE vs Control",          x = 0.07, y = 0.92),
#       list(tag = "B", title = "Mito vs Control",          x = 0.51, y = 0.92),
#       ...),
#     header_title    = "H9c2 volcano composites",
#     header_subtitle = "Pi-score = P.Value^|logFC|")
#
# strip_for_composite() above strips tag/title/subtitle from a single panel.
finalize_composite_yvo <- function(composite, comp_h_mm,
                                    panels,
                                    header_title = NULL,
                                    header_subtitle = NULL,
                                    header_y_title = 0.99,
                                    header_y_sub   = 0.965) {
  txt <- composite_text_sizes(comp_h_mm)
  TAG_SZ <- txt$tag      + 4
  TTL_SZ <- txt$title    + 2
  SUB_SZ <- txt$subtitle + 2
  X_OFF  <- 0.040
  Y_SUB  <- 0.020

  out <- cowplot::ggdraw(composite)
  for (p in panels) {
    if (!is.null(p$tag) && nzchar(p$tag)) {
      out <- out + cowplot::draw_label(
        p$tag, x = p$x, y = p$y + 0.002,
        size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1)
    }
    if (!is.null(p$title) && nzchar(p$title)) {
      out <- out + cowplot::draw_label(
        p$title, x = p$x + X_OFF, y = p$y,
        size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1)
    }
    if (!is.null(p$subtitle) && nzchar(p$subtitle)) {
      out <- out + cowplot::draw_label(
        p$subtitle, x = p$x + X_OFF, y = p$y - Y_SUB,
        size = SUB_SZ, fontface = "bold.italic", colour = "grey40",
        hjust = 0, vjust = 1)
    }
  }
  if (!is.null(header_title) && nzchar(header_title)) {
    out <- out + cowplot::draw_label(
      header_title, x = 0.02, y = header_y_title,
      size = txt$title + 4, fontface = "bold", hjust = 0, vjust = 1)
  }
  if (!is.null(header_subtitle) && nzchar(header_subtitle)) {
    out <- out + cowplot::draw_label(
      header_subtitle, x = 0.02, y = header_y_sub,
      size = txt$subtitle + 1, fontface = "italic",
      colour = "grey30", hjust = 0, vjust = 1)
  }
  out
}

get_pdf_device <- function() {
  # cairo_pdf > quartz > base pdf
  tryCatch(
    { cairo_pdf(tempfile()); dev.off(); cairo_pdf },
    error = function(e) tryCatch(
      {
        fp <- tempfile(fileext = ".pdf")
        quartz(type = "pdf", file = fp); dev.off()
        function(filename, width, height, ...)
          quartz(file = filename, type = "pdf", width = width, height = height)
      },
      error = function(e) "pdf"))
}
