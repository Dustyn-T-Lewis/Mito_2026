# style.R — palettes, themes, sizing helpers for figure scripts.
# Groups, contrasts, palettes, and thresholds are defined inline here (the figure
# stage's source of truth); the 00-03 pipeline scripts carry their own copies.

library(ggplot2)
library(scales)
library(grid)

# --- groups / contrasts / palettes / thresholds (match the 00-03 scripts) ----
H9C2_GROUP_LEVELS   <- c("Ctl", "Mito", "PHE", "PHE_Mito")
H9C2_CORE_CONTRASTS <- c("CTLvPHE", "CTLvMITO", "PHEvPHE_MITO", "Interaction")
H9C2_CONTRAST_ROLES <- c(CTLvPHE = "Disease", CTLvMITO = "Intervention",
                         PHEvPHE_MITO = "Rescue", Interaction = "Interaction",
                         MITOvPHE_MITO = "Secondary")
H9C2_PAL_GROUP <- c(Ctl = "#4393C3", Mito = "#009E73", PHE = "#D6604D", PHE_Mito = "#984EA3")
H9C2_PAL_DIR   <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")
H9C2_PAL_CONTRAST <- c(CTLvPHE = "#D6604D", CTLvMITO = "#4393C3", PHEvPHE_MITO = "#4DAF4A",
                       Interaction = "#7B5EA7", MITOvPHE_MITO = "#FF8C00")
H9C2_PI_THRESH  <- 0.05
H9C2_FDR_EXPLOR <- 0.10

# suppress stray Rplots.pdf from implicit device opens
options(device = function(...) grDevices::pdf(file = nullfile(), ...))

# Palettes (aliases used by figure scripts)
GROUP_COLORS    <- H9C2_PAL_GROUP            # Ctl / Mito / PHE / PHE_Mito
DIR_COLORS      <- H9C2_PAL_DIR              # Up / Down / NS
CONTRAST_COLORS <- H9C2_PAL_CONTRAST
CONTRAST_ROLES  <- H9C2_CONTRAST_ROLES

# Subcellular-compartment ring palette (Okabe-Ito; colourblind-safe). Used as the
# point OUTLINE on the F03/F04 protein scatters (fill = significance class) and
# the mitochondrial outline on the NES scatter. Lookup: protein_localization_rat.csv.
LOC_COLORS <- c(Mitochondrial = "#D55E00", Nuclear = "#332288",
                Cytosolic = "#117733", Other = "grey70")

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

# Fisher z CI for a correlation (k = number of covariates). Spearman rho has an
# inflated sampling variance vs Pearson: Var(z) = (1 + r^2/2)/(n-3) (Bonett &
# Wright 2000), so callers reporting a Spearman rho must pass method="spearman".
fisher_z_ci <- function(r, n, k = 0, level = 0.95,
                        method = c("pearson", "spearman")) {
  method <- match.arg(method)
  n_eff <- n - k
  if (n_eff < 4 || is.na(r)) return(c(lo = NA_real_, hi = NA_real_))
  z <- atanh(r)
  var_z <- if (method == "spearman") (1 + r^2 / 2) / (n_eff - 3) else 1 / (n_eff - 3)
  crit <- qnorm(1 - (1 - level) / 2)
  c(lo = tanh(z - crit * sqrt(var_z)), hi = tanh(z + crit * sqrt(var_z)))
}

sig_stars <- function(padj) {
  dplyr::case_when(
    padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ "")
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
