# PILOT 6 — RRHO2 Disease <-> Rescue threshold-free concordance map.
# Sourced by 01_main_panels.R; reads shared state (group_mat, ALL_GENES, SEED,
# MAIN_*, FIG_*, pdf_dev) and writes results$pilot_rrho2 + sheet_specs.
if (requireNamespace("RRHO2", quietly = TRUE)) {
  message("pilot_rrho2")
  cw <- load_combined_wide()
  # Rank by signed limma moderated t (Smyth 2004 doi:10.2202/1544-6115.1027):
  # variance-stabilized at small n, matches fgsea upstream so RRHO2 and fgsea
  # agree on which gene "leads" each contrast. Replaces Cahill 2018 (PMID
  # 29942049) signed -log10(P) — the rank, not the metric, drives the
  # hypergeometric tail.
  rank_vec <- function(tstat) {
    tstat[!is.finite(tstat)] <- 0
    tstat
  }
  d_rank <- rank_vec(cw$t_CTLvPHE)
  r_rank <- rank_vec(cw$t_PHEvPHE_MITO)
  keep <- !is.na(d_rank) & !is.na(r_rank) & !is.na(cw$gene) & nzchar(cw$gene)
  rrho_df <- slice_max(
    tibble(
      gene = cw$gene[keep],
      score_d = d_rank[keep],
      score_r = r_rank[keep],
      abs_d = abs(d_rank[keep])
    ),
    abs_d,
    n = 1, by = gene, with_ties = FALSE
  )
  list1 <- data.frame(gene = rrho_df$gene, score = rrho_df$score_d)
  list2 <- data.frame(gene = rrho_df$gene, score = rrho_df$score_r)
  set.seed(SEED)
  rr <- RRHO2::RRHO2_initialize(list1, list2,
    labels = c("Disease (signed mod-t)", "Rescue (signed mod-t)"),
    log10.ind = TRUE,
    stepsize = ceiling(sqrt(nrow(list1))),
    boundary = 0.025
  ) # per spec; trims 2.5% from each tail

  # Heatmap (separate file; RRHO2_heatmap uses base graphics + its own device).
  # Axis: Disease (x) and Rescue (y) run from low rank (most UP, top-left) to
  # high rank (most DOWN, bottom-right). UU=top-left, DD=bottom-right,
  # UD=top-right (reversed: D up / R down), DU=bottom-left (reversed: D down / R up).
  # RRHO2_heatmap leaves the device active on the color-bar panel, so plain
  # text() would land on the legend strip. Reserving outer margins before the
  # call and using mtext(outer=TRUE) puts the four corner labels on the
  # device's outer frame instead.
  .add_rrho2_quadrant_labels <- function() {
    mtext("UU",
      side = 3, line = 0.4, adj = 0.04, outer = TRUE,
      font = 2, cex = 1.3, col = "grey15"
    )
    mtext("UD  (reversed)",
      side = 3, line = 0.4, adj = 0.78, outer = TRUE,
      font = 2, cex = 1.3, col = "#B2182B"
    )
    mtext("DU  (reversed)",
      side = 1, line = 0.4, adj = 0.04, outer = TRUE,
      font = 2, cex = 1.3, col = "#D6604D"
    )
    mtext("DD",
      side = 1, line = 0.4, adj = 0.78, outer = TRUE,
      font = 2, cex = 1.3, col = "grey15"
    )
  }
  hm_pdf <- file.path(MAIN_PDF, "MAIN_F06_pilot_rrho2_heatmap.pdf")
  hm_png <- file.path(MAIN_PNG, "MAIN_F06_pilot_rrho2_heatmap.png")
  pdf(hm_pdf, width = 7, height = 7)
  par(oma = c(2, 1, 2, 1))
  RRHO2::RRHO2_heatmap(rr)
  .add_rrho2_quadrant_labels()
  dev.off()
  png(hm_png, width = 1800, height = 1800, res = 300)
  par(oma = c(2, 1, 2, 1))
  RRHO2::RRHO2_heatmap(rr)
  .add_rrho2_quadrant_labels()
  dev.off()

  # Quadrant gene lists. rr$genelist_* has gene_list1_*, gene_list2_*,
  # gene_list_overlap_*. When the peak-overlap intersection is < 5 (sparse
  # concordance is common), fall back to genes ranking in the top RRHO2_PCT
  # of both lists in the quadrant's implied direction — gives biologically
  # meaningful sets for trajectory + ORA even with tiny peak overlap.
  RRHO2_PCT <- 0.20
  n_top <- ceiling(RRHO2_PCT * nrow(rrho_df))
  d_asc <- order(rrho_df$score_d)
  d_desc <- order(rrho_df$score_d, decreasing = TRUE)
  r_asc <- order(rrho_df$score_r)
  r_desc <- order(rrho_df$score_r, decreasing = TRUE)
  top_d_up <- rrho_df$gene[d_desc[seq_len(n_top)]]
  top_d_down <- rrho_df$gene[d_asc[seq_len(n_top)]]
  top_r_up <- rrho_df$gene[r_desc[seq_len(n_top)]]
  top_r_down <- rrho_df$gene[r_asc[seq_len(n_top)]]
  uu_peak <- rr$genelist_uu$gene_list_overlap_uu
  dd_peak <- rr$genelist_dd$gene_list_overlap_dd
  ud_peak <- rr$genelist_ud$gene_list_overlap_ud
  du_peak <- rr$genelist_du$gene_list_overlap_du
  pick_genes <- function(peak, pct_fallback) {
    if (length(peak) >= 5) {
      list(genes = as.character(peak), source = "rrho2_peak")
    } else {
      cands <- if (is.null(pct_fallback)) character(0) else as.character(pct_fallback)
      list(genes = cands, source = "pct_fallback")
    }
  }
  quad_res <- list(
    UU = pick_genes(uu_peak, intersect(top_d_up, top_r_up)),
    DD = pick_genes(dd_peak, intersect(top_d_down, top_r_down)),
    UD = pick_genes(ud_peak, intersect(top_d_up, top_r_down)),
    DU = pick_genes(du_peak, intersect(top_d_down, top_r_up))
  )
  quad_lists <- lapply(quad_res, `[[`, "genes")
  quad_source <- vapply(quad_res, `[[`, character(1), "source")
  message(sprintf(
    "RRHO2 quadrant sizes: UU=%d DD=%d UD=%d DU=%d",
    length(quad_lists$UU), length(quad_lists$DD),
    length(quad_lists$UD), length(quad_lists$DU)
  ))
  pal_q <- c(UU = "#2E7D32", DD = "#1565C0", UD = "#B2182B", DU = "#D6604D")
  quad_role <- c(
    UU = "Concordant Up", DD = "Concordant Down",
    UD = "Reversed (Disease Up / Rescue Down)",
    DU = "Reversed (Disease Down / Rescue Up)"
  )
  rows <- lapply(names(quad_lists), function(q) {
    g_in <- quad_lists[[q]]
    if (length(g_in) < 5) {
      return(NULL)
    }
    g_in <- intersect(g_in, rownames(group_mat))
    if (length(g_in) < 5) {
      return(NULL)
    }
    z_cl <- standardise_genes(group_mat[g_in, , drop = FALSE])
    color <- pal_q[q]
    hall <- run_hallmark_ora(g_in, universe = ALL_GENES)
    mito <- run_mitocarta_ora(g_in, universe = ALL_GENES)
    hdr <- sprintf("Quadrant %s  |  n = %d  |  %s", q, length(g_in), quad_role[q])
    build_cluster_row(
      traj_plot = build_trajectory_panel(z_cl,
        x_levels = H9C2_GROUP_LEVELS,
        x_lab = "condition (group mean)",
        color = color, kind = "line"
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
  rows <- rows[!vapply(rows, is.null, logical(1))]

  fallback_quads <- names(quad_source)[quad_source == "pct_fallback"]
  fallback_note <- if (length(fallback_quads)) {
    sprintf(
      " | %s from top-20%% percentile (RRHO2 peak overlap <5 genes)",
      paste(fallback_quads, collapse = "/")
    )
  } else {
    ""
  }
  fig <- stack_cluster_rows(rows,
    title = "F06 pilot_rrho2 — Disease<->Rescue quadrants (per-quadrant ORA)",
    subtitle = paste0(
      "heatmap saved separately; rows = RRHO2 quadrants; ",
      "Hallmark + MitoCarta top-6", fallback_note
    )
  )
  h_mm <- 32 + 32 * length(rows)
  ggsave(file.path(MAIN_PDF, "MAIN_F06_pilot_rrho2.pdf"), fig,
    width = FIG_W, height = h_mm, units = "mm", device = pdf_dev, limitsize = FALSE
  )
  ggsave(file.path(MAIN_PNG, "MAIN_F06_pilot_rrho2.png"), fig,
    width = FIG_W, height = h_mm, units = "mm", dpi = 300, limitsize = FALSE
  )

  # Workbook payloads
  genelist_long <- bind_rows(lapply(names(quad_lists), function(q) {
    tibble(
      quadrant = q, role = quad_role[q], gene = quad_lists[[q]],
      source = quad_source[q]
    )
  }))
  ora_rr <- bind_rows(lapply(names(quad_lists), function(q) {
    o <- run_hallmark_ora(quad_lists[[q]], universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) {
      return(NULL)
    }
    mutate(o, quadrant = q, role = quad_role[q])
  }))
  ora_rr_mito <- bind_rows(lapply(names(quad_lists), function(q) {
    o <- run_mitocarta_ora(quad_lists[[q]], universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) {
      return(NULL)
    }
    mutate(o, quadrant = q, role = quad_role[q])
  }))

  results$pilot_rrho2 <- list(
    key = "pilot_rrho2",
    sheets = list(
      pilot_rrho2_genelists = genelist_long,
      pilot_rrho2_ora       = if (is.null(ora_rr)) tibble() else ora_rr,
      pilot_rrho2_ora_mito  = if (is.null(ora_rr_mito)) tibble() else ora_rr_mito
    )
  )
} else {
  warning("RRHO2 package missing — skipping pilot_rrho2")
}

if (!is.null(results$pilot_rrho2)) {
  rr <- results$pilot_rrho2$sheets
  sheet_specs <- c(
    sheet_specs,
    list(list(
      name = "pilot_rrho2_genelists", df = rr$pilot_rrho2_genelists,
      role = "Per-quadrant gene lists from RRHO2 (UU/DD/UD/DU)",
      contents = "quadrant (UU=concordant up, DD=concordant down, UD/DU=reversed), role, gene, source (rrho2_peak=from gene_list_overlap_*, pct_fallback=top-20% rank intersection)"
    )),
    list(list(
      name = "pilot_rrho2_ora", df = rr$pilot_rrho2_ora,
      role = "Hallmark ORA per RRHO2 quadrant",
      contents = "quadrant, role, pathway, padj, overlap, size, odds_ratio"
    )),
    list(list(
      name = "pilot_rrho2_ora_mito", df = rr$pilot_rrho2_ora_mito,
      role = "MitoCarta ORA per RRHO2 quadrant",
      contents = "quadrant, role, pathway, padj, overlap, size, odds_ratio"
    ))
  )
  sheet_specs[[1]]$df <- bind_rows(
    sheet_specs[[1]]$df,
    tibble(
      Pilot = "pilot_rrho2",
      Method = "RRHO2 threshold-free Disease<->Rescue map",
      Gate = "all genes with non-NA P.Value+logFC in Disease and Rescue",
      N_genes = nrow(rrho_df), # input universe, not long-format output rows
      Fuzzifier_m = NA_real_, Cluster_c = 4L
    )
  )
}
