# Shared Pattern Heatmap Panel — Panel B in F04 (concordance) and F05 (reversal)
# Requires: cfg list — see F04/F05 wrappers for full specification.
#
# The cfg must provide: fig_id, contrast_x, contrast_y, title, col_headers,
# classify_fn (function(df) returning df with quadrant + sig_cat columns),
# QUAD_ORDER, QUAD_COLORS, QUAD_BG, ENDPOINT_COLORS, SIG_COLORS,
# display_labels, bar_scale, bar_ref_width, key_y_base, key_dy, key_x_sig,
# sig_cats, sig_cat_labels, count_tick_filter, protein_count_x_mult,
# count_tick_y_label, rpt_png, rpt_pdf, dat, sort_col

source(here::here("04_Figures", "shared", "style.R"))
source(here::here("04_Figures", "shared", "print_scale_apply_380mm.R"))
source(here::here("04_Figures", "shared", "go_slim_categories.R"))

library(tidyverse)
library(patchwork)

RPT_PNG <- cfg$rpt_png
RPT_PDF <- cfg$rpt_pdf
DAT     <- cfg$dat
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_B_heatmap"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# 1. LOAD & CLASSIFY
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

# cfg$classify_fn does the figure-specific filtering and classification
sig_df <- cfg$classify_fn(dep_df)

QUAD_ORDER      <- cfg$QUAD_ORDER
QUAD_COLORS     <- cfg$QUAD_COLORS
QUAD_BG         <- cfg$QUAD_BG
ENDPOINT_COLORS <- cfg$ENDPOINT_COLORS
SIG_COLORS      <- cfg$SIG_COLORS

go_result <- assign_go_slim_consolidated(sig_df$gene, dep_df$gene)
sig_df <- sig_df |>
  left_join(go_result |> select(gene, consolidated), by = "gene") |>
  mutate(pathway = ifelse(is.na(consolidated), "Other", as.character(consolidated)))

sig_df <- sig_df |>
  mutate(quadrant = factor(quadrant, levels = QUAD_ORDER)) |>
  arrange(quadrant, pathway, desc(.data[[cfg$sort_col]]))

n_total <- nrow(sig_df)
message(sprintf("  %d significant proteins across %d quadrants", n_total,
                n_distinct(sig_df$quadrant)))

# 2. Y-COORDINATE LAYOUT
ROW_H <- 0.078

quad_counts <- sig_df |> count(quadrant, .drop = FALSE) |>
  mutate(quadrant = factor(quadrant, levels = QUAD_ORDER)) |> arrange(quadrant)

y_pos <- numeric(n_total)
quad_starts <- numeric(nrow(quad_counts))
quad_ends   <- numeric(nrow(quad_counts))
idx <- 1; current_y <- 0

for (q in seq_len(nrow(quad_counts))) {
  nq <- quad_counts$n[q]
  quad_starts[q] <- current_y
  if (nq > 0) for (p in seq_len(nq)) {
    y_pos[idx] <- current_y + (p - 0.5) * ROW_H
    idx <- idx + 1
  }
  quad_ends[q] <- current_y + nq * ROW_H
  current_y <- current_y + nq * ROW_H
}
total_h <- current_y
sig_df$y <- y_pos
names(quad_starts) <- QUAD_ORDER
names(quad_ends)   <- QUAD_ORDER

BAR_FRAC <- 1.0
BAR_YMIN <- 0
BAR_YMAX <- total_h * BAR_FRAC

# 3. PATHWAY LAYOUT — rows GROUPED BY DOMINANT QUADRANT for interpretability
# (Concordant Up block, then Concordant Down, then Discordant), each sorted by
# protein count. The per-row background (QUAD_BG[dom_quad]) then forms contiguous
# quadrant-coloured bands instead of being interleaved.
cx <- cfg$contrast_x
cy <- cfg$contrast_y
lfc_x <- paste0("logFC_", cx)
lfc_y <- paste0("logFC_", cy)

pw_counts <- sig_df |>
  filter(pathway != "Other") |>
  count(pathway, name = "n_prot") |>
  filter(n_prot >= 2)
n_pw <- nrow(pw_counts)

# Dominant quadrant per pathway (most proteins; ties broken by summed |logFC|).
dom_quad <- sig_df |>
  filter(pathway %in% pw_counts$pathway) |>
  group_by(pathway, quadrant) |>
  summarise(n = n(),
            lfc_sum = sum(abs(.data[[lfc_x]]) + abs(.data[[lfc_y]])),
            .groups = "drop") |>
  group_by(pathway) |>
  arrange(pathway, desc(n), desc(lfc_sum)) |>
  slice_head(n = 1) |>
  ungroup() |>
  mutate(dom_quad = as.character(quadrant)) |>
  select(pathway, dom_quad)

pw_counts <- pw_counts |>
  left_join(dom_quad, by = "pathway") |>
  mutate(dom_quad = factor(dom_quad, levels = QUAD_ORDER)) |>
  arrange(dom_quad, desc(n_prot)) |>
  mutate(dom_quad = as.character(dom_quad))

row_height <- (BAR_YMAX - BAR_YMIN) / n_pw
pw_counts$y_center <- BAR_YMIN + row_height * (seq_len(n_pw) - 0.5)
pw_counts$y_top    <- BAR_YMIN + row_height * (seq_len(n_pw) - 1)
pw_counts$y_bot    <- BAR_YMIN + row_height * seq_len(n_pw)
BAR_H <- row_height * 0.78

# 4. X-COORDINATE LAYOUT
STRIP_W <- 0.10; TILE_W <- cfg$tile_w %||% 0.70

X_SIG   <- cfg$x_sig %||% 0.8
X_COL1  <- X_SIG + STRIP_W/2 + TILE_W/2 + 0.01
X_COL2  <- X_COL1 + TILE_W + 0.01
X_QUAD  <- X_COL2 + TILE_W/2 + STRIP_W/2 + 0.01
HEAT_LEFT  <- X_SIG - STRIP_W/2
HEAT_RIGHT <- X_QUAD + STRIP_W/2

X_SANK_L  <- HEAT_RIGHT + 0.08
X_SANK_R  <- 3.2
X_BAR_L   <- 3.3
BAR_SCALE <- cfg$bar_scale

count_max <- max(pw_counts$n_prot)
X_BAR_MAX <- max(X_BAR_L + cfg$bar_ref_width * BAR_SCALE, X_BAR_L + count_max * BAR_SCALE)

PW_OUT <- 178; PH_OUT <- 130

# 5. STACKED BAR DATA
bar_data <- sig_df |>
  filter(pathway %in% pw_counts$pathway) |>
  count(pathway, quadrant, name = "n_seg") |>
  left_join(pw_counts |> select(pathway, y_center, n_prot), by = "pathway") |>
  group_by(pathway) |>
  arrange(pathway, desc(n_seg)) |>
  mutate(
    cum_n = cumsum(n_seg) - n_seg,
    xmin = X_BAR_L + cum_n * BAR_SCALE,
    xmax = X_BAR_L + (cum_n + n_seg) * BAR_SCALE,
    ymin = y_center - BAR_H / 2,
    ymax = y_center + BAR_H / 2
  ) |> ungroup()

# Right extent: bars + room for the (wrapped) pathway labels; the shaded band
# runs all the way to this edge so labels sit on the quadrant colour.
X_BAR_RIGHT <- X_BAR_MAX + (cfg$bg_extend_right %||% 4.5)

bg_stripes <- pw_counts |>
  transmute(
    xmin = X_BAR_L - 0.05, xmax = X_BAR_RIGHT,
    ymin = y_top, ymax = y_bot,
    fill = QUAD_BG[dom_quad]
  )

pw_labels <- pw_counts |>
  mutate(display_pathway = ifelse(pathway %in% names(cfg$display_labels),
                                   cfg$display_labels[pathway], pathway),
         display_pathway = gsub("\n", " ", display_pathway),
         display_pathway = stringr::str_wrap(display_pathway, width = 15)) |>
  transmute(x = X_BAR_L + n_prot * BAR_SCALE + 0.12, y = y_center,
            label = display_pathway)

count_ticks_max <- cfg$count_ticks_max %||% count_max
count_ticks <- tibble(
  val = pretty(c(0, count_ticks_max), n = 4),
  x = X_BAR_L + val * BAR_SCALE,
  y_tick_top = BAR_YMAX, y_tick_bot = BAR_YMAX + ROW_H * 1.6,
  y_label = BAR_YMAX + cfg$count_tick_y_label
) |> filter(val >= 0, val <= count_ticks_max) |>
  (\(df) cfg$count_tick_filter(df))()

# 6. SANKEY
flow_df <- sig_df |>
  filter(pathway %in% pw_counts$pathway) |>
  count(quadrant, pathway, name = "n_flow") |>
  filter(n_flow > 0)

source_bands <- flow_df |>
  group_by(quadrant) |>
  mutate(
    total_q = sum(n_flow), frac = n_flow / total_q,
    q_start = quad_starts[as.character(quadrant)],
    q_end   = quad_ends[as.character(quadrant)],
    q_height = q_end - q_start
  ) |>
  arrange(quadrant, match(pathway, pw_counts$pathway)) |>
  mutate(
    cum_frac = cumsum(frac) - frac,
    src_top = q_start + cum_frac * q_height,
    src_bot = q_start + (cum_frac + frac) * q_height
  ) |> ungroup()

target_bands <- bar_data |>
  group_by(pathway) |>
  arrange(pathway, desc(n_seg)) |>
  mutate(
    frac = n_seg / sum(n_seg),
    cum_frac = cumsum(frac) - frac,
    tgt_top = ymin + cum_frac * (ymax - ymin),
    tgt_bot = ymin + (cum_frac + frac) * (ymax - ymin)
  ) |> ungroup() |>
  select(pathway, quadrant, tgt_top, tgt_bot)

ribbon_df <- source_bands |>
  select(quadrant, pathway, n_flow, src_top, src_bot) |>
  left_join(target_bands, by = c("quadrant", "pathway"))

all_ribbons <- pmap_dfr(ribbon_df, function(quadrant, pathway, n_flow,
                                              src_top, src_bot, tgt_top, tgt_bot) {
  rid <- paste(quadrant, pathway, sep = "___")
  df <- make_sigmoid_ribbon(X_SANK_L, X_SANK_R, src_top, src_bot, tgt_top, tgt_bot,
                            n_pts = 60, ribbon_id = rid)
  df$quadrant <- quadrant; df$pathway <- pathway; df
})

endpoint_bars <- bar_data |>
  transmute(xmin = X_SANK_R - 0.04, xmax = X_SANK_R + 0.04,
            ymin, ymax, quadrant = as.character(quadrant))

# 7. HEATMAP
fc_max <- max(abs(c(sig_df[[lfc_x]], sig_df[[lfc_y]])), na.rm = TRUE)

lfc_to_color <- function(v, fc_max) {
  v <- pmax(-fc_max, pmin(fc_max, v))
  ifelse(v >= 0,
         scales::seq_gradient_pal("#FFFFFF", "#B2182B")(v / fc_max),
         scales::seq_gradient_pal("#2166AC", "#FFFFFF")((v + fc_max) / fc_max))
}

heat_tiles <- bind_rows(
  sig_df |> transmute(x = X_COL1, y, w = TILE_W, h = ROW_H,
                        fill = lfc_to_color(.data[[lfc_x]], fc_max)),
  sig_df |> transmute(x = X_COL2, y, w = TILE_W, h = ROW_H,
                        fill = lfc_to_color(.data[[lfc_y]], fc_max))
)
sig_tiles <- sig_df |>
  transmute(x = X_SIG, y, w = STRIP_W, h = ROW_H, fill = SIG_COLORS[sig_cat])
quad_tiles <- sig_df |>
  transmute(x = X_QUAD, y, w = STRIP_W, h = ROW_H,
            fill = QUAD_COLORS[as.character(quadrant)])

divider_ys <- quad_ends[1:(length(QUAD_ORDER) - 1)]
divider_ys <- divider_ys[divider_ys > 0 & divider_ys < total_h]

col_headers <- tibble(x = c(X_COL1, X_COL2), y = total_h + ROW_H * 2.2,
                      label = cfg$col_headers,
                      color = unname(cfg$col_header_colors))

# 8. LEGENDS
n_g <- 50
HEAT_MID  <- (HEAT_LEFT + HEAT_RIGHT) / 2
GRAD_HALFW <- (HEAT_RIGHT - HEAT_LEFT) * 0.30
GRAD_L    <- HEAT_MID - GRAD_HALFW
GRAD_R    <- HEAT_MID + GRAD_HALFW
grad_xs <- seq(GRAD_L, GRAD_R, length.out = n_g)
grad_h_legend <- tibble(
  xmin = grad_xs,
  xmax = lead(grad_xs, default = max(grad_xs) + diff(grad_xs)[1]),
  fv = seq(-fc_max, fc_max, length.out = n_g),
  fill = lfc_to_color(seq(-fc_max, fc_max, length.out = n_g), fc_max)
)
GRAD_Y <- total_h + ROW_H * 2.9

FONT_UNI <- 3.0 * PRINT_SCALE - 2
FONT_BAR <- 2.0 * PRINT_SCALE - 1
FONT_PW  <- FONT_UNI * 0.58 + 1.0   # smaller; labels are wrapped to 2 lines

KEY_SQ_SZ   <- 3.0 * PRINT_SCALE
KEY_TXT_SZ  <- 2.3 * PRINT_SCALE
KEY_HDR_SZ  <- 3.0 * PRINT_SCALE
KEY_Y_BASE  <- BAR_YMAX + cfg$key_y_base
KEY_DY      <- cfg$key_dy

bar_mid     <- (X_BAR_L + X_BAR_MAX) / 2
KEY_X_SIG   <- cfg$key_x_sig %||% (bar_mid - 1.43)
KEY_X_QUAD  <- bar_mid + 0.6

sig_key_df <- tibble(
  x     = KEY_X_SIG,
  y     = KEY_Y_BASE + KEY_DY * seq_along(cfg$sig_cats),
  label = cfg$sig_cat_labels,
  fill  = SIG_COLORS[cfg$sig_cats]
)

quad_key_df <- tibble(
  x     = KEY_X_QUAD,
  y     = KEY_Y_BASE + KEY_DY * seq_along(QUAD_ORDER),
  label = QUAD_ORDER,
  fill  = QUAD_COLORS[QUAD_ORDER]
)

# 9. RENDER
pB <- ggplot() +
  geom_rect(data = bg_stripes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = bg_stripes$fill, color = "grey70", linewidth = 0.2) +
  geom_rect(data = heat_tiles,
            aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2),
            fill = heat_tiles$fill, color = NA) +
  geom_rect(data = sig_tiles,
            aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2),
            fill = sig_tiles$fill, color = NA) +
  geom_rect(data = quad_tiles,
            aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2),
            fill = quad_tiles$fill, color = NA) +
  geom_segment(data = tibble(y = divider_ys),
               aes(x = X_SIG - STRIP_W/2, xend = X_QUAD + STRIP_W/2,
                   y = y, yend = y),
               color = "grey30", linewidth = 0.4) +
  geom_text(data = col_headers, aes(x = x, y = y, label = label),
            size = FONT_UNI - 0.5, fontface = "bold", color = "grey20") +
  geom_polygon(data = all_ribbons, aes(x = x, y = y, group = ribbon_id),
               fill = QUAD_COLORS[all_ribbons$quadrant], alpha = 0.40, color = NA) +
  geom_rect(data = endpoint_bars,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = ENDPOINT_COLORS[endpoint_bars$quadrant], color = NA) +
  geom_rect(data = bar_data,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = QUAD_COLORS[as.character(bar_data$quadrant)],
            color = "black", linewidth = 0.3) +
  geom_text(data = bar_data,
            aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = n_seg),
            size = FONT_BAR, fontface = "bold", color = "white") +
  geom_text(data = pw_labels, aes(x = x, y = y, label = label),
            size = FONT_PW, hjust = 0, fontface = "bold", color = "grey15",
            lineheight = 0.8) +
  annotate("segment", x = X_BAR_L, xend = X_BAR_MAX,
           y = BAR_YMAX, yend = BAR_YMAX, color = "grey20", linewidth = 0.5) +
  geom_segment(data = count_ticks,
               aes(x = x, xend = x, y = y_tick_top, yend = y_tick_bot),
               color = "grey20", linewidth = 0.3) +
  geom_text(data = count_ticks, aes(x = x, y = y_label, label = val),
            size = FONT_UNI * 0.78 + 0.5, fontface = "bold", color = "grey20") +
  annotate("text", x = X_BAR_L + (count_ticks_max / 2) * BAR_SCALE,
           y = BAR_YMAX + ROW_H * 5.4,
           label = "Protein count", size = FONT_UNI - 0.8, fontface = "bold",
           color = "grey20", hjust = 0.5) +
  scale_y_reverse() +
  coord_cartesian(xlim = c(0.0, X_BAR_RIGHT + 0.3),
                  ylim = c(BAR_YMAX + ROW_H * 6.6, -ROW_H * 0.05),
                  expand = FALSE) +
  labs(title = cfg$title,
       subtitle = sprintf("%d proteins | GO Slim | %d pathways", n_total, n_pw)) +
  theme_void() +
  theme(plot.margin = margin(6, -30, 38, -12, "mm"),
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE, hjust = 0,
                                  margin = margin(l = 31.5, unit = "mm")),
        plot.subtitle = element_text(face = "italic", size = FIG_SUBTITLE_SIZE,
                                     hjust = 0, color = "grey40",
                                     margin = margin(l = 31.5, unit = "mm")),
        plot.title.position = "panel")

# 10. INSET LEGEND (cfg$inset_legend == TRUE; bottom-right inside bar plot)
inset_quad_df <- tibble(
  quadrant  = factor(QUAD_ORDER, levels = QUAD_ORDER),
  y         = rev(seq_along(QUAD_ORDER)),
  bar_color = unname(QUAD_COLORS[QUAD_ORDER]),
  bg_color  = unname(QUAD_BG[QUAD_ORDER])
)
bar_xmin   <- 0.10
bar_xmax   <- 0.10 + 0.95 * 0.8    # half of original bar width
bar_half_h <- 0.26 * 1.6           # 1.6x bar height (original)
label_x    <- bar_xmax + 0.18
label_xmax <- label_x + max(nchar(QUAD_ORDER)) * 0.18
inset_legend <- ggplot(inset_quad_df) +
  geom_rect(aes(xmin = 0, xmax = label_xmax, ymin = y - 0.5, ymax = y + 0.5),
            fill = inset_quad_df$bg_color, color = "black", linewidth = 0.55) +
  geom_rect(aes(xmin = bar_xmin, xmax = bar_xmax,
                ymin = y - bar_half_h, ymax = y + bar_half_h),
            fill = inset_quad_df$bar_color, color = "black", linewidth = 0.18) +
  geom_text(aes(x = label_x, y = y, label = quadrant),
            hjust = 0, size = 2.8, fontface = "bold", color = "grey15") +
  scale_x_continuous(limits = c(0, label_xmax), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, length(QUAD_ORDER) + 0.5), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm"))

INSET_BOUNDS <- cfg$inset_bounds %||%
                list(left = 0.63, right = 0.96, bottom = 0.09, top = 0.30)

pB_standalone <- if (isTRUE(cfg$inset_legend)) {
  pB + inset_element(inset_legend,
                     left   = INSET_BOUNDS$left,   right = INSET_BOUNDS$right,
                     bottom = INSET_BOUNDS$bottom, top   = INSET_BOUNDS$top,
                     align_to = "panel")
} else pB

# 11. SAVE
ggsave(file.path(RPT_PNG, "MAIN_panel_B_pattern_heatmap.png"), pB_standalone,
       width = PW_OUT, height = PH_OUT, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_B_pattern_heatmap.pdf"), pB_standalone,
       width = PW_OUT, height = PH_OUT, units = "mm", device = pdf_device)

# 11. DATA EXPORTS
sig_df |>
  transmute(gene, quadrant = as.character(quadrant), sig_cat, pathway,
            !!lfc_x := round(.data[[lfc_x]], 4),
            !!lfc_y := round(.data[[lfc_y]], 4)) |>
  write_csv(file.path(DAT, "panel_B_heatmap", "pattern_classification.csv"))
flow_df |> write_csv(file.path(DAT, "panel_B_heatmap", "sankey_links.csv"))
bar_data |>
  select(pathway, quadrant, n_seg, xmin, xmax) |>
  write_csv(file.path(DAT, "panel_B_heatmap", "bar_data.csv"))

# (Legend keys now embedded directly in the panel plot above — no separate PNG needed)

# Export for composite
pB_title    <- cfg$title
pB_subtitle <- sprintf("%d proteins | GO Slim | %d pathways", n_total, n_pw)
pB_legend   <- NULL
pB          <- strip_for_composite(pB)
# Tighten for composite: reduce bottom whitespace (standalone keeps full padding)
pB <- pB +
  coord_cartesian(xlim = c(0.0, X_BAR_RIGHT + 0.3),
                  ylim = c(BAR_YMAX + ROW_H * 5.6, -ROW_H * 0.05),
                  expand = FALSE) +
  theme(plot.margin = margin(-1, -28, 4, -14, "mm"))

if (isTRUE(cfg$inset_legend)) {
  pB <- pB + inset_element(inset_legend,
                            left   = INSET_BOUNDS$left,   right = INSET_BOUNDS$right,
                            bottom = INSET_BOUNDS$bottom, top   = INSET_BOUNDS$top,
                            align_to = "panel")
}

message(sprintf("%s Panel B (pattern heatmap) done", cfg$fig_id))
