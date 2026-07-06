# Palettes, theme, and sizing helpers. Groups, contrasts, and thresholds live
# here as the figure stage's source of truth; the 00-03 scripts carry their own.

library(ggplot2)
library(scales)
library(grid)

# Groups, contrasts, palettes, thresholds (match the 00-03 scripts).
H9C2_GROUP_LEVELS <- c("Ctl", "Mito", "PHE", "PHE_Mito")
# Display annotations for the four conditions (legends, axis strips, contrast math).
# Internal tokens above drive data joins and must not change; only these strings show.
H9C2_GROUP_LABELS <- c(Ctl = "CTL", Mito = "MitoTx", PHE = "Phe-only", PHE_Mito = "Phe+MitoTx")
H9C2_PAL_GROUP <- c(Ctl = "#3B7DB5", Mito = "#009E73", PHE = "#E08214", PHE_Mito = "#8073AC")
H9C2_PAL_DIR <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")
H9C2_PAL_DIR_MITO <- c(Up = "#B2182B", Down = "#2166AC") # darker mito subset
H9C2_CONTRAST_ORDER <- c("CTLvPHE", "CTLvMITO", "PHEvPHE_MITO", "Interaction")
H9C2_PAL_CONTRAST <- c(
  CTLvPHE = "#D6604D", CTLvMITO = "#4393C3", PHEvPHE_MITO = "#4DAF4A",
  Interaction = "#7B5EA7", MITOvPHE_MITO = "#FF8C00"
)
H9C2_PI_THRESH <- 0.05
H9C2_FDR_EXPLOR <- 0.10

# Enrichment-database palette, shared by every figure that colours pathways by source.
# Colourblind-safe (Paul Tol), all dark enough to carry white bold labels on the bars.
DB_COLORS <- c(
  Hallmark = "#CC3311", Reactome = "#0077BB", KEGG = "#AA3377",
  MitoCarta = "#009988", `GO Slim` = "#555555",
  `GO:BP` = "#117733", `GO:CC` = "#44AA99", `GO:MF` = "#882255"
)

# Lighter qualitative palette for the named ORA pathway-bar panels (F01's ORA
# supplement and Figure 1F's stacked counts, F03's per-module top-5) — DB_COLORS
# above stays as the F02 GSEA supplement's palette, which wasn't part of this restyle.
# Significance text inside these bars uses shadowtext so it stays legible at this
# lighter saturation.
ORA_DB_COLORS <- c(
  Hallmark = "#F7C08A", Reactome = "#8FC9BC", KEGG = "#B3A5D9",
  MitoCarta = "#9CCDEE", `GO Slim` = "#BFBFBF",
  `GO:BP` = "#A3CBB0", `GO:CC` = "#E8B8D0", `GO:MF` = "#E0A985"
)

# suppress stray Rplots.pdf from implicit device opens
options(device = function(...) grDevices::pdf(file = nullfile(), ...))

# Palettes (aliases used by figure scripts)
GROUP_COLORS <- H9C2_PAL_GROUP # Ctl / Mito / PHE / PHE_Mito
GROUP_LABELS <- H9C2_GROUP_LABELS
DIR_COLORS <- H9C2_PAL_DIR # Up / Down / NS
DIR_COLORS_MITO <- H9C2_PAL_DIR_MITO
CONTRAST_COLORS <- H9C2_PAL_CONTRAST

# Sizing (J Physiol double-column spec)
PANEL_MD <- 178
BASE_STAT <- 2.5
BASE_TAG <- 8

scale_text <- function(base_size, panel_width_mm, ref_width = PANEL_MD) {
  base_size * sqrt(panel_width_mm / ref_width)
}

# Text hierarchy + theme
FIG_TITLE_SIZE <- 7
FIG_SUBTITLE_SIZE <- 4
FIG_STRIP_SIZE <- 5
FIG_AXIS_TEXT <- 5
FIG_LEGEND_TITLE <- 5
FIG_LEGEND_TEXT <- 4

FIG_THEME <- theme_bw(base_size = 6, base_family = "Helvetica") +
  theme(
    plot.title = element_text(
      face = "bold", size = FIG_TITLE_SIZE,
      margin = margin(b = 1)
    ),
    plot.subtitle = element_text(
      face = "bold.italic", size = FIG_SUBTITLE_SIZE,
      color = "grey30", margin = margin(t = 0, b = 2)
    ),
    plot.tag = element_text(face = "bold", size = BASE_TAG),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = FIG_STRIP_SIZE),
    axis.title.x = element_text(face = "bold", size = 5, margin = margin(t = 0)),
    axis.title.y = element_text(face = "bold", size = 5, margin = margin(r = -1)),
    axis.text = element_text(size = FIG_AXIS_TEXT, color = "grey15"),
    legend.title = element_text(
      face = "bold", size = FIG_LEGEND_TITLE,
      color = "grey20"
    ),
    legend.text = element_text(size = FIG_LEGEND_TEXT, color = "grey15"),
    legend.key.size = unit(2.5, "mm"),
    panel.grid.minor = element_blank()
  )

fmt_p <- function(p) {
  if (p < 0.001) {
    return("p < 0.001")
  }
  if (p < 0.01) {
    return(sprintf("p = %.3f", p))
  }
  sprintf("p = %.2f", p)
}
