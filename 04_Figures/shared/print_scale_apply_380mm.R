# Print-scale override for 380mm-wide composites.
# Used only by F04 / F05 380mm-wide main composites.
#
# WARNING: this script MUTATES style.R's exported globals
# (FIG_TITLE_SIZE, FIG_SUBTITLE_SIZE, FIG_STRIP_SIZE, FIG_AXIS_TEXT,
#  FIG_LEGEND_TITLE, FIG_LEGEND_TEXT, FIG_THEME, BASE_*) by multiplying
# them by PRINT_SCALE (≈ 2.13). Any code sourced AFTER this script
# sees the upscaled values, not style.R's defaults.
#
# Source order for 380mm composites: style.R first, this second,
# then panel scripts.

PRINT_SCALE <- 380 / 178   # ≈ 2.13

FIG_TITLE_SIZE    <- round(FIG_TITLE_SIZE * PRINT_SCALE, 1)
FIG_SUBTITLE_SIZE <- round(FIG_SUBTITLE_SIZE * PRINT_SCALE, 1)
FIG_STRIP_SIZE    <- round(FIG_STRIP_SIZE * PRINT_SCALE, 1)
FIG_AXIS_TEXT     <- round(FIG_AXIS_TEXT * PRINT_SCALE, 1)
FIG_LEGEND_TITLE  <- round(FIG_LEGEND_TITLE * PRINT_SCALE, 1)
FIG_LEGEND_TEXT   <- round(FIG_LEGEND_TEXT * PRINT_SCALE, 1)

# drop 2pt to avoid overfilling
BASE_PATHWAY  <- round(BASE_PATHWAY * PRINT_SCALE, 1) - 2
BASE_GENE     <- round(BASE_GENE * PRINT_SCALE, 1) - 2
BASE_STAT     <- round(BASE_STAT * PRINT_SCALE, 1) - 2
BASE_QUADRANT <- round(BASE_QUADRANT * PRINT_SCALE, 1) - 2
BASE_COUNT    <- round(BASE_COUNT * PRINT_SCALE, 1) - 2

FIG_THEME <- theme_bw(base_size = round(6 * PRINT_SCALE, 1), base_family = "Helvetica") +
  theme(
    plot.title         = element_text(face = "bold", size = FIG_TITLE_SIZE,
                                      margin = margin(b = 1)),
    plot.subtitle      = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                      color = "grey30", margin = margin(t = 0, b = 2)),
    plot.tag           = element_text(face = "bold", size = FIG_TITLE_SIZE),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = FIG_STRIP_SIZE),
    axis.title.x       = element_text(face = "bold", size = FIG_AXIS_TEXT,
                                      margin = margin(t = 0)),
    axis.title.y       = element_text(face = "bold", size = FIG_AXIS_TEXT,
                                      margin = margin(r = -1)),
    axis.text          = element_text(size = FIG_AXIS_TEXT, color = "grey15"),
    legend.title       = element_text(face = "bold", size = FIG_LEGEND_TITLE,
                                      color = "grey20"),
    legend.text        = element_text(size = FIG_LEGEND_TEXT, color = "grey15"),
    legend.key.size    = unit(2.5 * PRINT_SCALE, "mm"),
    panel.grid.minor   = element_blank()
  )
