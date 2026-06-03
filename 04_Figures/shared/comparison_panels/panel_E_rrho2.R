# Shared RRHO2 Panel E — used by F04 (concordance) and F05 (reversal)

source(here::here("04_Figures", "shared", "style.R"))
source(here::here("04_Figures", "shared", "print_scale_apply_380mm.R"))
source(here::here("04_Figures", "shared", "pathway_utils.R"))
library(tidyverse)
library(ggrepel)
library(patchwork)
library(msigdbr)
library(fgsea)
library(RRHO2)

PE_W <- 89

RPT_PNG <- cfg$rpt_png
RPT_PDF <- cfg$rpt_pdf
DAT     <- cfg$dat
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_E"), recursive = TRUE, showWarnings = FALSE)

# Supplementary output dirs (F05 only)
if (!is.null(cfg$supp)) {
  dir.create(cfg$supp$rpt_png, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$supp$rpt_pdf, recursive = TRUE, showWarnings = FALSE)
}

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

rr_df <- dep_df |>
  transmute(gene,
            t_1 = .data[[cfg$t_col_1]],
            t_2 = .data[[cfg$t_col_2]]) |>
  filter(!is.na(t_1) & !is.na(t_2)) |>
  distinct(gene, .keep_all = TRUE)

n_shared <- nrow(rr_df)

# stratified hypergeometric test (Cahill et al. 2018)
list1 <- data.frame(gene = rr_df$gene, score = rr_df$t_1, stringsAsFactors = FALSE)
list2 <- data.frame(gene = rr_df$gene, score = rr_df$t_2, stringsAsFactors = FALSE)

rrho_obj <- RRHO2_initialize(
  list1, list2,
  labels          = cfg$rrho_labels,
  log10.ind       = TRUE,
  multipleTesting = "none",
  boundary        = 0.02,
  method          = "hyper",
  stepsize        = 20
)

hmat <- rrho_obj$hypermat
nr <- nrow(hmat); nc <- ncol(hmat)
message(sprintf("  RRHO2 matrix: %d x %d", nr, nc))

na_rows <- which(apply(hmat, 1, function(r) all(is.na(r))))
na_cols <- which(apply(hmat, 2, function(c) all(is.na(c))))

if (length(na_rows) > 0 && length(na_cols) > 0) {
  row_before <- 1:(min(na_rows) - 1)
  row_after  <- (max(na_rows) + 1):nr
  col_before <- 1:(min(na_cols) - 1)
  col_after  <- (max(na_cols) + 1):nc
} else {
  mid <- floor(nr / 2)
  row_before <- 1:mid
  row_after  <- (mid + 1):nr
  col_before <- 1:mid
  col_after  <- (mid + 1):nc
}

# Both lists sorted descending (Up first):
# [row_before, col_before] = top of both lists = UU
# [row_after,  col_after]  = bottom of both    = DD
# [row_before, col_after]  = top list1, bottom list2 = UD
# [row_after,  col_before] = bottom list1, top list2 = DU
max_UU <- max(hmat[row_before, col_before], na.rm = TRUE)
max_DD <- max(hmat[row_after,  col_after],  na.rm = TRUE)
max_UD <- max(hmat[row_before, col_after],  na.rm = TRUE)
max_DU <- max(hmat[row_after,  col_before], na.rm = TRUE)

# Quadrant density: aggregate enrichment signal per quadrant block so the four
# quadrants can be compared directly ("more concordant-down than -up"). mean is
# the fair comparison (block pixel-counts differ); sum is the total signal mass.
.quad_density <- function(block) {
  v <- block[is.finite(block)]
  c(mean = if (length(v)) mean(v) else 0,
    sum  = sum(v),
    n_px = length(v))
}
dens_UU <- .quad_density(hmat[row_before, col_before])
dens_DD <- .quad_density(hmat[row_after,  col_after])
dens_UD <- .quad_density(hmat[row_before, col_after])
dens_DU <- .quad_density(hmat[row_after,  col_before])

hotspot_genes <- list(
  UU = rrho_obj$genelist_uu$gene_list_overlap_uu,
  DD = rrho_obj$genelist_dd$gene_list_overlap_dd,
  UD = rrho_obj$genelist_ud$gene_list_overlap_ud,
  DU = rrho_obj$genelist_du$gene_list_overlap_du
)

n_UU <- length(hotspot_genes$UU)
n_DD <- length(hotspot_genes$DD)
n_UD <- length(hotspot_genes$UD)
n_DU <- length(hotspot_genes$DU)
message(sprintf("  Hotspot genes: UU=%d, DD=%d, UD=%d, DU=%d", n_UU, n_DD, n_UD, n_DU))

# no multiplier — prevents corner overlap in composite
txt_quad <- scale_text(BASE_QUADRANT, 146)

# jet colormap (Cahill et al. 2018 canonical appearance)
JET_COLORS <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                "yellow", "#FF7F00", "red", "#7F0000")

hmat_df <- expand.grid(row = 1:nr, col = 1:nc) |>
  mutate(neg_log10_p = as.vector(hmat))

max_val <- max(hmat_df$neg_log10_p, na.rm = TRUE)

# Corner anchors
ann_x_left  <- min(row_before)
ann_x_right <- max(row_after)
ann_y_bot   <- min(col_before)
ann_y_top   <- max(col_after)

LABEL_FILL    <- scales::alpha("white", 0.85)
LABEL_PADDING <- unit(0.8 * PRINT_SCALE, "mm")

ql <- cfg$quadrant_labels

pE_heat <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_p)) +
  geom_raster() +
  scale_fill_gradientn(
    colors   = JET_COLORS,
    limits   = c(0, max_val),
    na.value = "white",
    name     = expression(-log[10](P)),
    guide    = guide_colorbar(
      barwidth = unit(15 * PRINT_SCALE, "mm"), barheight = unit(2 * PRINT_SCALE, "mm"),
      title.position = "left", title.vjust = 0.5,
      title.theme = element_text(size = FIG_LEGEND_TITLE, face = "bold", color = "grey15"))
  ) +
  annotate("label", x = ann_x_left, y = ann_y_bot,
           label = sprintf("%s\nμ=%.2f | n=%d", ql$UU, dens_UU["mean"], n_UU),
           color = "grey15", fill = LABEL_FILL, linewidth = 0,
           label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
           fontface = "bold", size = txt_quad, lineheight = 0.9,
           hjust = 0, vjust = 0) +
  annotate("label", x = ann_x_right, y = ann_y_top,
           label = sprintf("%s\nμ=%.2f | n=%d", ql$DD, dens_DD["mean"], n_DD),
           color = "grey15", fill = LABEL_FILL, linewidth = 0,
           label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
           fontface = "bold", size = txt_quad, lineheight = 0.9,
           hjust = 1, vjust = 1) +
  annotate("label", x = ann_x_left, y = ann_y_top,
           label = sprintf("%s\nμ=%.2f | n=%d", ql$UD, dens_UD["mean"], n_UD),
           color = "grey15", fill = LABEL_FILL, linewidth = 0,
           label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
           fontface = "bold", size = txt_quad, lineheight = 0.9,
           hjust = 0, vjust = 1) +
  annotate("label", x = ann_x_right, y = ann_y_bot,
           label = sprintf("%s\nμ=%.2f | n=%d", ql$DU, dens_DU["mean"], n_DU),
           color = "grey15", fill = LABEL_FILL, linewidth = 0,
           label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
           fontface = "bold", size = txt_quad, lineheight = 0.9,
           hjust = 1, vjust = 0) +
  scale_x_continuous(expand = expansion(mult = 0.015)) +
  scale_y_continuous(expand = expansion(mult = 0.015)) +
  labs(
    title    = cfg$title,
    subtitle = sprintf(cfg$subtitle_fmt, n_shared),
    x = cfg$axis_label_1,
    y = cfg$axis_label_2
  ) +
  FIG_THEME +
  theme(
    axis.text        = element_blank(),
    axis.title.x     = element_text(size = FIG_AXIS_TEXT, face = "bold", margin = margin(t = 2)),
    axis.title.y     = element_text(size = FIG_AXIS_TEXT, face = "bold", margin = margin(r = 2)),
    axis.ticks       = element_blank(),
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = FIG_LEGEND_TEXT, face = "bold"),
    legend.margin    = margin(2, 24, 0, 0, "mm"),
    plot.margin = margin(0, 0, 0, 0, "mm")
  ) +
  coord_fixed(ratio = 1, clip = "off")

he <- cfg$hotspot_export_names
hotspot_export <- bind_rows(
  tibble(quadrant = he$UU, gene = hotspot_genes$UU),
  tibble(quadrant = he$DD, gene = hotspot_genes$DD),
  tibble(quadrant = he$UD, gene = hotspot_genes$UD),
  tibble(quadrant = he$DU, gene = hotspot_genes$DU)
)
write_csv(hotspot_export, file.path(DAT, "panel_E", "rrho2_hotspot_genes.csv"))

ora_min <- cfg$ora_min_size %||% 15
# Harmonized backbone (Hallmark + Reactome + MitoCarta + GO Slim) — same as the
# panel A wings and panel C NES scatter, so RRHO2 hotspot ORA surfaces MitoCarta.
pw_collection_E <- build_harmonized_collection(min_size = ora_min, max_size = 350,
                                               include_goslim = TRUE)
all_genes_E <- rr_df$gene

run_quadrant_ora <- function(gene_set, quadrant_name) {
  if (length(gene_set) < 5) return(tibble())
  res <- tryCatch(
    run_ora_deduplicated(
      genes          = gene_set,
      universe       = all_genes_E,
      pathways       = pw_collection_E,
      jaccard_cutoff = 0.5,
      min_size       = ora_min,
      max_size       = 500,
      padj_cutoff    = 0.05
    ),
    error = function(e) { message("  ORA error: ", e$message); tibble() }
  )
  if (nrow(res) > 0) {
    res |>
      mutate(quadrant = quadrant_name,
             pathway_label = clean_pathway_name(pathway),
             ID = pathway,
             p.adjust = padj,
             GeneRatio = paste0(overlap, "/", length(gene_set)),
             geneID = sapply(overlapGenes, paste, collapse = "/")) |>
      arrange(padj, size)
  } else {
    tibble()
  }
}

oq <- cfg$ora_quadrant_names
ora_UU <- run_quadrant_ora(hotspot_genes$UU, oq$UU)
ora_DD <- run_quadrant_ora(hotspot_genes$DD, oq$DD)
ora_UD <- run_quadrant_ora(hotspot_genes$UD, oq$UD)
ora_DU <- run_quadrant_ora(hotspot_genes$DU, oq$DU)

og <- cfg$ora_grouped
ora_group_1 <- bind_rows(mget(og$file_1_quads))
ora_group_2 <- bind_rows(mget(og$file_2_quads))

write_csv(ora_group_1, file.path(DAT, "panel_E", "rrho2_ora_concordant.csv"))
write_csv(ora_group_2, file.path(DAT, "panel_E", "rrho2_ora_discordant.csv"))

if (!is.null(og$note_if_empty_2) && nrow(ora_group_2) == 0) {
  write_csv(tibble(note = og$note_if_empty_2),
            file.path(DAT, "panel_E", "rrho2_ora_discordant_note.csv"))
}

ggsave(file.path(RPT_PNG, "MAIN_panel_E_rrho2.png"), pE_heat,
       width = PE_W, height = PE_W, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_E_rrho2.pdf"), pE_heat,
       width = PE_W, height = PE_W, units = "mm", device = pdf_device)

MAX_PER_QUAD <- 12

if (!is.null(cfg$supp)) {
  txt_ora <- scale_text(BASE_STAT, PE_W)
  ora_all <- bind_rows(ora_group_1, ora_group_2)

  if (nrow(ora_all) > 0) {
    bar_df <- ora_all |>
      mutate(
        neg_log10_padj = -log10(p.adjust),
        pathway_label  = str_trunc(clean_pathway_name(pathway), 40),
        quadrant       = factor(quadrant, levels = cfg$supp$ora_quad_order)
      ) |>
      filter(!is.na(quadrant)) |>
      group_by(quadrant, .drop = FALSE) |>
      arrange(desc(neg_log10_padj)) |>
      slice_head(n = MAX_PER_QUAD) |>
      ungroup() |>
      filter(!is.na(neg_log10_padj)) |>
      arrange(quadrant, neg_log10_padj) |>
      mutate(uid = fct_inorder(paste0(pathway_label, "___", quadrant)))

    n_shown <- nrow(bar_df)
    n_total <- nrow(ora_all)

    pE_ora <- ggplot(bar_df, aes(x = neg_log10_padj, y = uid, fill = quadrant)) +
      geom_col(width = 0.75) +
      geom_text(aes(label = overlap), hjust = -0.3, size = txt_ora * 0.7,
                color = "grey30") +
      scale_y_discrete(labels = function(x) str_remove(x, "___.*$")) +
      scale_fill_manual(values = cfg$ora_colors, guide = "none") +
      scale_x_continuous(expand = expansion(mult = c(0, 0))) +
      facet_grid(quadrant ~ ., scales = "free_y", space = "free_y",
                 labeller = labeller(quadrant = cfg$supp$ora_quad_short)) +
      labs(title = cfg$supp$ora_bar_title,
           subtitle = if (n_shown < n_total)
             sprintf("Top %d per quadrant (%d terms total)", MAX_PER_QUAD, n_total)
           else sprintf("%d terms total", n_total),
           x = expression(-log[10](p[adj])),
           y = NULL) +
      FIG_THEME +
      theme(
        plot.title       = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(size = 8, hjust = 0.5),
        strip.text.y     = element_text(size = 7, face = "bold", angle = 0),
        strip.background = element_rect(fill = "grey95", color = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.text.y  = element_text(size = 7),
        plot.margin  = margin(2, 4, 2, 2, "mm")
      )

    pE_combined <- (pE_heat | pE_ora) + plot_layout(widths = c(1, 1.3))

    ggsave(file.path(cfg$supp$rpt_png, "SUPP_panel_E_rrho2_ora.png"), pE_combined,
           width = 400, height = 220, units = "mm", dpi = 300)
    ggsave(file.path(cfg$supp$rpt_pdf, "SUPP_panel_E_rrho2_ora.pdf"), pE_combined,
           width = 400, height = 220, units = "mm", device = pdf_device)
  }
}

sq <- cfg$summary_quadrant_names
quad_meta <- list(
  list(name = oq$UU, slug = sq$UU_slug, data = ora_UU, n_hot = n_UU),
  list(name = oq$DD, slug = sq$DD_slug, data = ora_DD, n_hot = n_DD),
  list(name = oq$UD, slug = sq$UD_slug, data = ora_UD, n_hot = n_UD),
  list(name = oq$DU, slug = sq$DU_slug, data = ora_DU, n_hot = n_DU)
)

for (qm in quad_meta) {
  if (nrow(qm$data) == 0) next
  q_df <- qm$data |>
    mutate(neg_log10_padj = -log10(padj),
           pathway_label  = clean_pathway_name(pathway)) |>
    arrange(desc(neg_log10_padj)) |>
    slice_head(n = MAX_PER_QUAD)
  message(sprintf("  %s: %d ORA pathways (%d hotspot genes)", qm$slug, nrow(q_df), qm$n_hot))
}

rrho2_meta <- tibble(
  quadrant = c(sq$UU, sq$DD, sq$UD, sq$DU),
  max_neg_log10_pvalue = round(c(max_UU, max_DD, max_UD, max_DU), 2),
  mean_neg_log10p = round(c(dens_UU["mean"], dens_DD["mean"], dens_UD["mean"], dens_DU["mean"]), 3),
  sum_neg_log10p  = round(c(dens_UU["sum"],  dens_DD["sum"],  dens_UD["sum"],  dens_DU["sum"]), 1),
  n_pixels = c(dens_UU["n_px"], dens_DD["n_px"], dens_UD["n_px"], dens_DU["n_px"]),
  n_hotspot_genes = c(n_UU, n_DD, n_UD, n_DU),
  n_ora_pathways = c(nrow(ora_UU), nrow(ora_DD), nrow(ora_UD), nrow(ora_DU)),
  matrix_rows = nr, matrix_cols = nc, n_shared_genes = n_shared
)
write_csv(rrho2_meta, file.path(DAT, "panel_E", "rrho2_summary.csv"))

# Focused quadrant-density table, ranked by mean signal (which quadrant is "larger").
quad_density_tbl <- rrho2_meta |>
  select(quadrant, mean_neg_log10p, sum_neg_log10p, n_pixels, n_hotspot_genes) |>
  arrange(desc(mean_neg_log10p)) |>
  mutate(density_rank = row_number())
write_csv(quad_density_tbl, file.path(DAT, "panel_E", "rrho2_quadrant_density.csv"))
message(sprintf("  Quadrant density (mean -log10P): %s",
                paste(sprintf("%s=%.2f", quad_density_tbl$quadrant, quad_density_tbl$mean_neg_log10p),
                      collapse = " | ")))

# Exposed for the composite stat-snapshot.
# Diagonal comparison (UU vs DD): for F03 = Concordant Up vs Down.
dens_concordant_cmp <- if (dens_DD["mean"] >= dens_UU["mean"])
  sprintf("%s denser (μ=%.2f vs %.2f)", sq$DD, dens_DD["mean"], dens_UU["mean"]) else
  sprintf("%s denser (μ=%.2f vs %.2f)", sq$UU, dens_UU["mean"], dens_DD["mean"])
# Off-diagonal vs diagonal (for F04 reversal = reversed vs exacerbated).
dens_reversed_cmp <- {
  rev_mean <- mean(c(dens_UD["mean"], dens_DU["mean"]))
  exa_mean <- mean(c(dens_UU["mean"], dens_DD["mean"]))
  if (rev_mean >= exa_mean)
    sprintf("reversed denser (μ=%.2f vs %.2f)", rev_mean, exa_mean) else
    sprintf("exacerbated denser (μ=%.2f vs %.2f)", exa_mean, rev_mean)
}

pE_legend_grob <- cowplot::get_plot_component(pE_heat, "guide-box-bottom", return_all = FALSE)
if (!is.null(pE_legend_grob)) {
  pE_legend_plot <- cowplot::ggdraw(pE_legend_grob)
  ggsave(file.path(RPT_PNG, "MAIN_panel_E_legend.png"), pE_legend_plot,
         width = 60, height = 14, units = "mm", dpi = 300)
}

pE_title    <- cfg$title
pE_subtitle <- sprintf(cfg$subtitle_fmt, n_shared)
pE_legend   <- NULL
pE_heat     <- pE_heat +
  labs(title = NULL, subtitle = NULL, tag = NULL) +
  coord_fixed(ratio = 1, clip = "off")

message(sprintf("%s Panel E done", cfg$fig_id))
