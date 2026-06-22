# PILOT 5 — k-means on per-protein 4-D logFC vector with gap-statistic k.
# Sourced by 01_main_panels.R; reads shared state (CORE, SEED, GAP_B, DEFAULT_C,
# MAIN_*, SUPP_*, FIG_*, pdf_dev) and writes results$pilot_logfc + sheet_specs.
message("pilot_logfc")
comb_wide <- load_combined_wide()
lf_cols <- paste0("logFC_", CORE)
stopifnot(all(lf_cols %in% names(comb_wide)))
lf <- comb_wide |>
  select(gene, all_of(lf_cols)) |>
  filter(
    !is.na(.data$gene), nzchar(.data$gene),
    if_all(all_of(lf_cols), ~ !is.na(.x))
  )
lf_mat <- as.matrix(lf[, lf_cols])
rownames(lf_mat) <- lf$gene
if (nrow(lf_mat) < 12L) {
  stop(sprintf("pilot_logfc: too few genes (%d) for clustering", nrow(lf_mat)))
}
# k selection: gap statistic firstSEmax (Tibshirani, Walther & Hastie 2001
# doi:10.1111/1467-9868.00293) — smallest k whose gap is within 1 SE of the max.
gap_pick <- pick_c_gap(lf_mat,
  k_range = 2:10, B = GAP_B, seed = SEED,
  default_c = DEFAULT_C
)
pilot_c_lf <- gap_pick$c
message(sprintf("pilot_logfc: %s", gap_pick$basis))
set.seed(SEED)
km <- kmeans(lf_mat, centers = pilot_c_lf, nstart = 50, iter.max = 100)

# Quadrant-style label from (logFC_Disease, logFC_Rescue) centroid signs.
disease_idx <- which(lf_cols == "logFC_CTLvPHE")
rescue_idx <- which(lf_cols == "logFC_PHEvPHE_MITO")
quadrant_label <- function(dis, res) {
  if (abs(dis) < 0.1 && abs(res) < 0.1) {
    return("Neutral")
  }
  if (dis > 0 && res < 0) {
    return("Reversed Down")
  } # disease up, rescue brings down
  if (dis < 0 && res > 0) {
    return("Reversed Up")
  }
  if (dis > 0 && res > 0) {
    return("Concordant Up")
  }
  if (dis < 0 && res < 0) {
    return("Concordant Down")
  }
  "Other"
}
centroids <- as_tibble(km$centers, rownames = "cluster") |>
  mutate(
    label = mapply(
      quadrant_label, km$centers[, disease_idx],
      km$centers[, rescue_idx]
    ),
    n = km$size
  )

pal <- cluster_palette(pilot_c_lf)
rows <- lapply(seq_len(pilot_c_lf), function(cl) {
  g_in <- rownames(lf_mat)[km$cluster == cl]
  color <- pal[as.character(cl)]
  # trajectory = mean-logFC bar per contrast (kind = "barlogfc")
  z_cl <- lf_mat[g_in, , drop = FALSE] # already in logFC space
  colnames(z_cl) <- gsub("^logFC_", "", colnames(z_cl))
  univ <- unique(rownames(lf_mat))
  hall <- run_hallmark_ora(g_in, universe = univ)
  mito <- run_mitocarta_ora(g_in, universe = univ)
  cent <- centroids[centroids$cluster == as.character(cl), ]
  hdr <- sprintf(
    "Cluster %d  |  n = %d  |  %s",
    cl, cent$n, cent$label
  )
  build_cluster_row(
    traj_plot = build_trajectory_panel(z_cl,
      x_levels = colnames(z_cl),
      x_lab = "contrast (centroid logFC)",
      color = color, kind = "barlogfc"
    ),
    ora_plot = build_ora_bar_panel(hall,
      color = color, max_n = 6,
      db_name = "Hallmark"
    ),
    ora_plot2 = build_ora_bar_panel(mito,
      color = color, max_n = 6,
      db_name = "MitoCarta"
    ),
    header_text = hdr, color = color
  )
})

fig <- stack_cluster_rows(rows,
  title = sprintf(
    "F06 pilot_logfc — k-means on per-protein 4-D logFC vector (c = %d; %s)",
    pilot_c_lf, gap_pick$basis
  ),
  subtitle = "rows = clusters; cluster labels derived from (Disease, Rescue) sign quadrant; middle = Hallmark top-6, right = MitoCarta top-6"
)
h_mm <- 32 + 32 * pilot_c_lf
ggsave(file.path(MAIN_PDF, "MAIN_F06_pilot_logfc.pdf"), fig,
  width = FIG_W, height = h_mm, units = "mm", device = pdf_dev, limitsize = FALSE
)
ggsave(file.path(MAIN_PNG, "MAIN_F06_pilot_logfc.png"), fig,
  width = FIG_W, height = h_mm, units = "mm", dpi = 300, limitsize = FALSE
)

# supp: gap-statistic diagnostic — gap(k) ± 1 SE bars with chosen c highlighted
gap_tbl <- gap_pick$gap_tbl
sup <- ggplot(gap_tbl, aes(.data$k, .data$gap)) +
  geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim),
    color = "grey55", width = 0.2, linewidth = 0.3
  ) +
  geom_line(color = "grey50", linewidth = 0.4) +
  geom_point(size = 1.4, color = "grey30") +
  geom_vline(
    xintercept = pilot_c_lf, color = "#D6604D",
    linetype = "dashed", linewidth = 0.4
  ) +
  geom_point(
    data = gap_tbl[gap_tbl$k == pilot_c_lf, ],
    color = "#D6604D", size = 2.8
  ) +
  scale_x_continuous(breaks = 2:10) +
  labs(
    title = "F06 pilot_logfc cluster-selection diagnostic",
    subtitle = sprintf("gap statistic (Tibshirani 2001); %s", gap_pick$basis),
    x = "number of clusters (k)", y = "gap(k) ± 1 SE"
  ) +
  FIG_THEME
ggsave(file.path(SUPP_PDF, "MAIN_F06_pilot_logfc_selection.pdf"), sup,
  width = 120, height = 80, units = "mm", device = pdf_dev
)
ggsave(file.path(SUPP_PNG, "MAIN_F06_pilot_logfc_selection.png"), sup,
  width = 120, height = 80, units = "mm", dpi = 300
)

memb_lf <- tibble(gene = rownames(lf_mat), cluster = as.integer(km$cluster))
univ_lf <- unique(rownames(lf_mat))
ora_lf <- bind_rows(lapply(seq_len(pilot_c_lf), function(cl) {
  g <- memb_lf$gene[memb_lf$cluster == cl]
  o <- run_hallmark_ora(g, universe = univ_lf)
  if (is.null(o) || nrow(o) == 0) {
    return(NULL)
  }
  mutate(o, cluster = cl)
}))
ora_lf_mito <- bind_rows(lapply(seq_len(pilot_c_lf), function(cl) {
  g <- memb_lf$gene[memb_lf$cluster == cl]
  o <- run_mitocarta_ora(g, universe = univ_lf)
  if (is.null(o) || nrow(o) == 0) {
    return(NULL)
  }
  mutate(o, cluster = cl)
}))

results$pilot_logfc <- list(
  key = "pilot_logfc",
  sheets = list(
    pilot_logfc_membership = memb_lf,
    pilot_logfc_centroids  = centroids,
    pilot_logfc_ora        = if (is.null(ora_lf)) tibble() else ora_lf,
    pilot_logfc_ora_mito   = if (is.null(ora_lf_mito)) tibble() else ora_lf_mito
  )
)

if (!is.null(results$pilot_logfc)) {
  lg <- results$pilot_logfc$sheets
  sheet_specs <- c(
    sheet_specs,
    list(list(
      name = "pilot_logfc_membership", df = lg$pilot_logfc_membership,
      role = "Gene -> k-means cluster on per-protein 4-D logFC vector",
      contents = "gene, cluster (1..c)"
    )),
    list(list(
      name = "pilot_logfc_centroids", df = lg$pilot_logfc_centroids,
      role = "k-means centroids in logFC space + quadrant label",
      contents = "cluster, logFC_<contrast>... mean centroid, quadrant label, n"
    )),
    list(list(
      name = "pilot_logfc_ora", df = lg$pilot_logfc_ora,
      role = "Hallmark ORA per logFC cluster",
      contents = "cluster, pathway, padj, overlap, size, odds_ratio"
    )),
    list(list(
      name = "pilot_logfc_ora_mito", df = lg$pilot_logfc_ora_mito,
      role = "MitoCarta ORA per logFC cluster",
      contents = "cluster, pathway, padj, overlap, size, odds_ratio"
    ))
  )
  sheet_specs[[1]]$df <- bind_rows(
    sheet_specs[[1]]$df,
    tibble(
      Pilot = "pilot_logfc", Method = "k-means on per-protein 4-D logFC vector",
      Gate = "all genes with non-NA logFC across core contrasts",
      N_genes = nrow(lg$pilot_logfc_membership),
      Fuzzifier_m = NA_real_, Cluster_c = pilot_c_lf
    )
  )
}
