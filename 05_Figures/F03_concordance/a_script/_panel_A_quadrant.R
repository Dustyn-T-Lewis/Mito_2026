#!/usr/bin/env Rscript
# F03 Panel A: YvO-style ORA composite — central quadrant scatter (Transplant_Mito
# CTLvMITO vs Rescue_Mito+Phe PHEvPHE_MITO) flanked by 4 corner ORA half-bar
# panels (top-5 enriched pathways per quadrant). Larger protein labels + points.
# Exports: pA (composite), n_total_A, n_sig_A, r_spear_A, rho_lo_A, rho_hi_A, pct_conc_A.

source(here::here("04_Figures", "shared", "print_scale_apply_380mm.R"))
source(here::here("04_Figures", "shared", "pathway_utils.R"))
library(fgsea)
library(ggrepel)
library(patchwork)

BASE    <- "05_Figures/F03_concordance"
RPT_PNG <- file.path(BASE, "b_reports", "main", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf", "panels")
DAT     <- file.path(BASE, "c_data")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_A"), recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

COMP_RED  <- unname(DIR_COLORS["Up"])
COMP_BLUE <- unname(DIR_COLORS["Down"])
N_SHOW    <- 5L

# Subcellular-compartment lookup (frozen cache; built by build_localization_lookup.R).
# Drives the point OUTLINE colour on the scatter (Mitochondrial / Nuclear /
# Cytosolic / Other) so localization is read at a glance instead of gene names.
LOC_LOOKUP <- readr::read_csv("04_Figures/shared/protein_localization_rat.csv",
                              show_col_types = FALSE)

# Data
dep_df <- readr::read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

imp_path <- "02_Imputation/c_data/02_mar_mnar_classification.csv"
imputation_df <- if (file.exists(imp_path)) {
  readr::read_csv(imp_path, show_col_types = FALSE) |>
    dplyr::transmute(gene, imputed = classification != "Complete") |>
    dplyr::distinct(gene, .keep_all = TRUE)
} else tibble::tibble(gene = character(), imputed = logical())

scatter_df <- dep_df |>
  dplyr::transmute(
    gene,
    logFC_x = logFC_CTLvMITO,      logFC_y = logFC_PHEvPHE_MITO,
    pi_Tr   = pi_score_CTLvMITO,   pi_Re   = pi_score_PHEvPHE_MITO,
    pi_Int  = pi_score_Interaction) |>
  dplyr::filter(!is.na(logFC_x), !is.na(logFC_y)) |>
  dplyr::distinct(gene, .keep_all = TRUE) |>
  dplyr::left_join(imputation_df, by = "gene", relationship = "many-to-one") |>
  dplyr::left_join(LOC_LOOKUP, by = "gene", relationship = "many-to-one") |>
  dplyr::mutate(
    imputed   = tidyr::replace_na(imputed, FALSE),
    localization = factor(tidyr::replace_na(localization, "Other"),
                          levels = names(LOC_COLORS)),
    is_mito   = localization == "Mitochondrial",
    sig_class = classify_proteins_f2(pi_Tr, pi_Re, pi_Int),
    is_sig    = sig_class != "NS",
    dist      = sqrt(logFC_x^2 + logFC_y^2),
    quadrant  = dplyr::case_when(
      logFC_x > 0 & logFC_y > 0 ~ "Concordant Up",
      logFC_x < 0 & logFC_y < 0 ~ "Concordant Down",
      logFC_x > 0 & logFC_y < 0 ~ "Discordant (Tr Up / Re Down)",
      TRUE                       ~ "Discordant (Tr Down / Re Up)"))

universe  <- scatter_df$gene
n_total_A <- nrow(scatter_df)
n_sig_A   <- sum(scatter_df$is_sig)

ct        <- cor.test(scatter_df$logFC_x, scatter_df$logFC_y, method = "spearman")
r_spear_A <- as.numeric(ct$estimate)
ci        <- fisher_z_ci(r_spear_A, n_total_A)
rho_lo_A  <- ci[1]; rho_hi_A <- ci[2]
pct_conc_A <- 100 * sum(sign(scatter_df$logFC_x) == sign(scatter_df$logFC_y)) / n_total_A

message(sprintf("  F03 panel A: %d proteins | %d DEPs | rho = %.3f [%.3f, %.3f] | %.0f%% concordant",
                n_total_A, n_sig_A, r_spear_A, rho_lo_A, rho_hi_A, pct_conc_A))

# Per-quadrant ORA (threshold-free, all proteins in the quadrant)
# Harmonized backbone (Hallmark + Reactome + MitoCarta + GO Slim) so the wings
# surface OXPHOS / TCA / mitoribosome / FAO as MitoCarta terms instead of the
# vague oversized GO:BP umbrella bins that previously won on set size alone.
pw_collection <- build_harmonized_collection(min_size = 10, max_size = 350,
                                             include_goslim = TRUE)

run_set_ora <- function(quad_name) {
  genes <- scatter_df$gene[scatter_df$quadrant == quad_name]
  if (length(genes) < 5) return(tibble::tibble())
  res <- run_ora_deduplicated(genes = genes, universe = universe,
                              pathways = pw_collection, jaccard_cutoff = 0.5,
                              min_size = 10, max_size = 350, padj_cutoff = 1)
  if (nrow(res) == 0) return(tibble::tibble())
  res |>
    dplyr::mutate(pathway_label = clean_pathway_name(pathway),
                  display_label = clean_display_label(pathway),
                  neg_log10_padj = -log10(padj), significant = padj < 0.05) |>
    dplyr::arrange(dplyr::desc(neg_log10_padj)) |>
    dplyr::slice_head(n = N_SHOW)
}
ora_q1 <- run_set_ora("Concordant Up")
ora_q2 <- run_set_ora("Discordant (Tr Down / Re Up)")
ora_q3 <- run_set_ora("Concordant Down")
ora_q4 <- run_set_ora("Discordant (Tr Up / Re Down)")
all_quad_ora <- dplyr::bind_rows(
  dplyr::mutate(ora_q1, set = "Concordant Up"),
  dplyr::mutate(ora_q2, set = "Discordant (Tr Down / Re Up)"),
  dplyr::mutate(ora_q3, set = "Concordant Down"),
  dplyr::mutate(ora_q4, set = "Discordant (Tr Up / Re Down)"))
if (nrow(all_quad_ora) > 0)
  readr::write_csv(all_quad_ora, file.path(DAT, "panel_A", "ora_quadrant.csv"))

# Central scatter (LARGER points + labels)
lim_x <- max(abs(scatter_df$logFC_x)) * 1.08
lim_y <- max(abs(scatter_df$logFC_y)) * 1.08
lim_x <- lim_y <- max(lim_x, lim_y)   # symmetric -> square panel (ratio 1, 45° diagonal)
ns_df  <- dplyr::filter(scatter_df, sig_class == "NS")
sig_df <- dplyr::filter(scatter_df, sig_class != "NS")
# draw mito points last (on top) so the mitochondrial ring group reads clearly
sig_plot <- dplyr::arrange(sig_df, localization == "Mitochondrial")
sig_alpha <- ifelse(sig_plot$sig_class == "Interaction", 0.65, 0.9)

q_tab <- scatter_df |>
  dplyr::count(quadrant, name = "n_total") |>
  dplyr::left_join(scatter_df |> dplyr::filter(is_sig) |> dplyr::count(quadrant, name = "n_sig"),
                   by = "quadrant") |>
  dplyr::mutate(n_sig = tidyr::replace_na(n_sig, 0L))
qn <- function(q, col) { v <- q_tab[[col]][q_tab$quadrant == q]; if (length(v) == 0) 0L else v }

# No gene labels: significant proteins are encoded by FILL = significance class
# and OUTLINE = subcellular compartment (ring), so the mitochondrial group is
# read at a glance without label clutter.
txt_quad <- scale_text(BASE_QUADRANT, 178) * 1.05

p_scatter <- ggplot(mapping = aes(x = logFC_x, y = logFC_y)) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_point(data = ns_df, color = "grey80", fill = "grey85", shape = 21,
             size = 0.7, alpha = 0.30, stroke = 0.12) +
  # significant proteins: fill = significance class; mitochondrial proteins get a
  # BOLD dark outline (others a thin grey outline) so the mito group reads at a glance.
  geom_point(data = sig_plot, aes(fill = sig_class), shape = 21, size = 2.1,
             alpha = sig_alpha,
             color = ifelse(sig_plot$is_mito, "grey10", "grey75"),
             stroke = ifelse(sig_plot$is_mito, 1.2, 0.3)) +
  scale_fill_manual(values = SIG_COLORS_F2, name = "Significance") +
  annotate("label", x = lim_x, y = lim_y,
           label = sprintf("Concordant Up\n%d/%d", qn("Concordant Up", "n_sig"), qn("Concordant Up", "n_total")),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold", color = COMP_RED,
           fill = alpha("white", 0.92), label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  annotate("label", x = -lim_x, y = -lim_y,
           label = sprintf("%d/%d\nConcordant Down", qn("Concordant Down", "n_sig"), qn("Concordant Down", "n_total")),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold", color = COMP_RED,
           fill = alpha("white", 0.92), label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  annotate("label", x = -lim_x, y = lim_y,
           label = sprintf("Discordant (Tr↓ Re↑)\n%d/%d", qn("Discordant (Tr Down / Re Up)", "n_sig"), qn("Discordant (Tr Down / Re Up)", "n_total")),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold", color = COMP_BLUE,
           fill = alpha("white", 0.92), label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  annotate("label", x = lim_x, y = -lim_y,
           label = sprintf("%d/%d\nDiscordant (Tr↑ Re↓)", qn("Discordant (Tr Up / Re Down)", "n_sig"), qn("Discordant (Tr Up / Re Down)", "n_total")),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold", color = COMP_BLUE,
           fill = alpha("white", 0.92), label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  # axis ticks + titles drawn INSIDE along the zero axes (YvO style) so the
  # panel borders are free for the ORA wings to butt against.
  geom_text(data = data.frame(v = setdiff(seq(-floor(lim_x), floor(lim_x)), 0)),
            aes(x = v, y = 0, label = v), inherit.aes = FALSE,
            vjust = 1.4, size = 1.4 * PRINT_SCALE, color = "grey45", fontface = "bold") +
  geom_text(data = data.frame(v = setdiff(seq(-floor(lim_y), floor(lim_y)), 0)),
            aes(x = 0, y = v, label = v), inherit.aes = FALSE,
            hjust = -0.4, size = 1.4 * PRINT_SCALE, color = "grey45", fontface = "bold") +
  annotate("text", x = 0, y = -lim_y * 0.97, vjust = 0, hjust = 0.5,
           label = "log₂FC  Transplant_Mito (Mito − Ctl)",
           size = 1.5 * PRINT_SCALE, color = "grey25", fontface = "bold") +
  annotate("text", x = -lim_x * 0.97, y = 0, angle = 90, vjust = 1, hjust = 0.5,
           label = "log₂FC  Rescue_Mito+Phe (PHE_Mito − PHE)",
           size = 1.5 * PRINT_SCALE, color = "grey25", fontface = "bold") +
  coord_fixed(ratio = lim_y / lim_x, xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y), expand = FALSE) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        legend.position = "none", plot.margin = margin(2, 0, 0, 0, "mm"))

# Half-bar ORA wing builder
make_half_bars <- function(df, fill_color, side, ylim) {
  bar_h <- 0.42
  # keep 0-enrichment (padj = 1 -> -log10 = 0) bars OFF the panel
  if (!is.null(df) && nrow(df) > 0) df <- dplyr::filter(df, neg_log10_padj > 0)
  n_bars <- if (is.null(df) || nrow(df) == 0) 0L else min(nrow(df), 5L)
  if (n_bars == 0)
    return(ggplot() + theme_void() + scale_y_continuous(limits = ylim, expand = c(0, 0)))
  y_pos <- if (ylim[1] >= 0) rev(seq(0.3, 2.3, length.out = 5))[seq_len(n_bars)]
           else seq(-0.3, -2.5, length.out = 5)[seq_len(n_bars)]
  bars <- df |>
    dplyr::arrange(dplyr::desc(neg_log10_padj)) |> dplyr::slice_head(n = 5) |>
    dplyr::mutate(
      y = y_pos,
      bar_fill = ifelse(significant, scales::alpha(fill_color, 0.85), scales::alpha(fill_color, 0.30)),
      display_name = stringr::str_wrap(stringr::str_trunc(display_label, 26), width = 15),
      star = sig_stars(padj))
  x_max <- max(bars$neg_log10_padj)
  # significant wings -> full-length bars (scale to actual max); non-significant
  # wings -> axis doubled so the weaker bars read shorter / de-emphasised.
  x_display_max <- if (any(bars$significant)) x_max * 1.10 else x_max * 2.0
  is_upper <- ylim[1] >= 0
  brk_fn <- function(limits) { b <- scales::pretty_breaks(n = 3)(limits); b[b != 0] }
  bars <- bars |>
    dplyr::mutate(
      label_inside = neg_log10_padj >= x_max * 0.50,
      label_x = ifelse(label_inside, neg_log10_padj * 0.5, neg_log10_padj + x_max * 0.03),
      label_hjust = ifelse(label_inside, 0.5, 0),
      label_color = ifelse(label_inside, ifelse(significant, "white", "grey15"), "grey20"),
      text_size = scale_text(BASE_PATHWAY, 190) * 0.72)   # smaller, wrapped
  star_x_mult <- if (side == "left") 0.12 else 0.035
  p <- ggplot(bars, aes(y = y)) +
    geom_rect(aes(xmin = 0, xmax = neg_log10_padj, ymin = y - bar_h / 2, ymax = y + bar_h / 2),
              fill = bars$bar_fill, color = "black", linewidth = 0.3) +
    geom_text(aes(x = label_x, y = y, label = display_name), hjust = bars$label_hjust,
              size = bars$text_size, fontface = "bold", color = bars$label_color, lineheight = 0.85) +
    geom_text(aes(x = neg_log10_padj + x_max * star_x_mult, label = star), hjust = 0, vjust = 0.5,
              size = 2.6 * PRINT_SCALE, fontface = "bold", color = "black") +
    labs(x = if (!is_upper) expression(-log[10](p[adj])) else NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = FIG_AXIS_TEXT, face = "bold"),
          axis.title.x = if (!is_upper) element_text(size = FIG_AXIS_TEXT, face = "bold") else element_blank(),
          axis.line.x = element_line(color = "grey50", linewidth = 0.3),
          axis.ticks.x = element_line(color = "grey50", linewidth = 0.3),
          plot.margin = if (is_upper && side == "left") margin(4, 0, 0, 3, "mm")
                        else if (is_upper) margin(4, 3, 0, 0, "mm")
                        else if (side == "left") margin(2, 0, 0, 3, "mm") else margin(2, 3, 0, 0, "mm"))
  if (side == "left")
    p + scale_x_reverse(limits = c(x_display_max, 0), breaks = brk_fn, expand = expansion(mult = c(0, 0))) +
        scale_y_continuous(limits = ylim, expand = c(0, 0)) + coord_cartesian(clip = "off")
  else
    p + scale_x_continuous(limits = c(0, x_display_max), breaks = brk_fn, expand = expansion(mult = c(0, 0))) +
        scale_y_continuous(limits = ylim, expand = c(0, 0)) + coord_cartesian(clip = "off")
}

p_ul <- make_half_bars(ora_q2, COMP_BLUE, "left",  c(0, 2.8))    # Discordant Tr↓Re↑
p_ll <- make_half_bars(ora_q3, COMP_RED,  "left",  c(-2.8, 0))   # Concordant Down
p_ur <- make_half_bars(ora_q1, COMP_RED,  "right", c(0, 2.8))    # Concordant Up
p_lr <- make_half_bars(ora_q4, COMP_BLUE, "right", c(-2.8, 0))   # Discordant Tr↑Re↓

# Significance key (mito = bold outline, explained in the panel subtitle) --
key_lvls    <- c("Sig Both", "Interaction", "Sig Intervention only", "Sig Rescue only")
key_display <- c("Sig Both", "Interaction", "Sig Int.", "Sig Res.")
key_df <- tibble::tibble(category = factor(key_lvls, levels = key_lvls), display = key_display,
                         fill_col = unname(SIG_COLORS_F2[key_lvls]), x = c(1.25, 1.80, 2.45, 3.05), y = 0)
p_key <- ggplot(key_df, aes(x = x, y = y)) +
  geom_point(aes(fill = category), shape = 21, size = 2.8 * PRINT_SCALE, color = "grey50",
             stroke = 0.6, alpha = 0.85, show.legend = FALSE) +
  geom_text(aes(label = display), nudge_x = 0.06, hjust = 0, size = 2.2 * PRINT_SCALE,
            fontface = "bold", color = "grey25") +
  scale_fill_manual(values = setNames(key_df$fill_col, key_df$category)) +
  scale_x_continuous(limits = c(0.2, 4.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.15, 0.15), expand = c(0, 0)) +
  coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(-22, 0, 0, 0, "mm"))

# Assemble composite
design <- c(area(1, 1), area(1, 2, 2, 2), area(1, 3), area(2, 1), area(2, 3), area(3, 1, 3, 3))
# scatter width share ~= its rendered height so the square fills the cell (no
# horizontal centering gap) -> ORA wings butt against the scatter borders.
composite <- p_ul + p_scatter + p_ur + p_ll + p_lr + p_key +
  plot_layout(design = design, widths = c(80, 80, 80) / 240, heights = c(85, 85, 8) / 178)

COMP_W <- 200; COMP_H <- 130
ggsave(file.path(RPT_PNG, "MAIN_panel_A_ORA_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_A_ORA_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)

readr::write_csv(
  scatter_df |> dplyr::transmute(gene, quadrant, sig_class = as.character(sig_class), is_sig, imputed,
    localization = as.character(localization),
    logFC_CTLvMITO = round(logFC_x, 4), logFC_PHEvPHE_MITO = round(logFC_y, 4),
    dist_from_origin = round(dist, 4)) |>
    dplyr::arrange(dplyr::desc(dist_from_origin)),
  file.path(DAT, "panel_A", "quadrant_proteins.csv"))

pA <- composite
message("F03 Panel A (ORA composite) done")
