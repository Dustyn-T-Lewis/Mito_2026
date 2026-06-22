# PILOT 4 — WGCNA modules from the pre-built signed-Pearson network.
# Sourced by 01_main_panels.R; reads shared state (group_mat, gene_mat, meta,
# ALL_GENES, MAIN_*, SUPP_*, FIG_*, pdf_dev) and appends to `results` and
# `sheet_specs` in the caller's environment.
WGCNA_RDS <- here::here("04_Figures", "F05_modules", "c_data", "wgcna_network.rds")
if (file.exists(WGCNA_RDS)) {
  message("pilot_wgcna")
  w <- load_wgcna_modules(WGCNA_RDS)
  mods <- w$modules |> filter(.data$module != "grey")
  MEs  <- w$MEs

  # Build the per-sample (Group) indicator vectors over the 4 contrasts the
  # eigengenes are correlated against. Group ordering must match MEs rownames.
  me_meta <- tibble(Col_ID = rownames(MEs)) |>
    left_join(as_tibble(meta) |> select(Col_ID, Group), by = "Col_ID")
  contrast_pairs <- list(
    Disease    = c("Ctl", "PHE"),
    Transplant = c("Ctl", "Mito"),
    Rescue     = c("PHE", "PHE_Mito"))
  me_corr <- compute_me_contrast_correlations(MEs, me_meta, contrast_pairs)

  # Interaction contrast: (PHE_Mito - Mito) - (PHE - Ctl) = +Ctl - Mito - PHE + PHE_Mito
  # aligned to (Ctl, Mito, PHE, PHE_Mito) = (+1, -1, -1, +1). Cannot use
  # compute_me_contrast_correlations (binary indicator only), so per-module cor.test.
  int_order <- H9C2_GROUP_LEVELS
  int_vec_map <- c(Ctl = 1, Mito = -1, PHE = -1, PHE_Mito = 1)
  int_indic <- int_vec_map[me_meta$Group]
  int_keep  <- !is.na(int_indic)
  if (sum(int_keep) >= 4) {
    int_rows <- lapply(colnames(MEs), function(mod) {
      ct <- suppressWarnings(cor.test(MEs[int_keep, mod], int_indic[int_keep],
                                      method = "pearson"))
      tibble(module = mod, contrast = "Interaction",
             r = unname(ct$estimate), p = ct$p.value)
    })
    me_corr <- bind_rows(me_corr, bind_rows(int_rows))
  }

  # r_floor 0.25 (vs helper default 0.40): at n = 24 the strict 0.40 floor
  # collapses every module into "Other"; 0.25 surfaces the documented 4 Reversal
  # + 1 Concordant Down set. Treat sign_pattern as a ranking aid, not an
  # inferential claim — Student p for r at n=24 needs |r| ~ 0.40 for p < 0.05.
  signs <- me_corr |>
    pivot_wider(id_cols = module, names_from = contrast, values_from = r,
                names_prefix = "r_") |>
    mutate(sign_pattern = classify_module_sign_pattern(r_Disease, r_Rescue,
                                                       r_floor = 0.25))
  mod_order <- signs |>
    arrange(.data$sign_pattern, desc(abs(.data$r_Rescue))) |>
    pull(.data$module)

  pal <- setNames(viridis::turbo(length(mod_order) + 2)[seq_along(mod_order) + 1],
                  mod_order)

  rows <- lapply(mod_order, function(mod) {
    g_in <- mods$gene[mods$module == mod]
    g_in <- intersect(g_in, rownames(gene_mat))
    if (length(g_in) < 5) return(NULL)
    z_cl <- standardise_genes(group_mat[g_in, , drop = FALSE])
    color <- pal[mod]
    sp   <- signs$sign_pattern[match(mod, signs$module)]
    rD   <- signs$r_Disease[match(mod, signs$module)]
    rR   <- signs$r_Rescue [match(mod, signs$module)]
    rT   <- signs$r_Transplant[match(mod, signs$module)]
    hdr  <- sprintf("Module %s  |  n = %d  |  %s  |  r(D)=%.2f r(R)=%.2f r(T)=%.2f",
                    mod, length(g_in), sp, rD, rR, rT)
    hall <- run_hallmark_ora (g_in, universe = ALL_GENES)
    mito <- run_mitocarta_ora(g_in, universe = ALL_GENES)
    build_cluster_row(
      traj_plot = build_trajectory_panel(z_cl,
                                         x_levels = H9C2_GROUP_LEVELS,
                                         x_lab = "condition (group mean)",
                                         color = color, kind = "line"),
      ora_plot  = build_ora_bar_panel(hall, color = color, max_n = 6,
                                      db_name = "Hallmark"),
      ora_plot2 = build_ora_bar_panel(mito, color = color, max_n = 6,
                                      db_name = "MitoCarta"),
      header_text = hdr, color = color)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]

  fig <- stack_cluster_rows(rows,
    title    = "F06 pilot_wgcna — modules from F05 WGCNA artifact",
    subtitle = "rows = modules; ordered by Disease<->Rescue sign pattern (reversal first); middle = Hallmark top-6, right = MitoCarta top-6")
  h_mm <- 32 + 32 * length(rows)
  ggsave(file.path(MAIN_PDF, "MAIN_F06_pilot_wgcna.pdf"), fig,
         width = FIG_W, height = h_mm, units = "mm", device = pdf_dev, limitsize = FALSE)
  ggsave(file.path(MAIN_PNG, "MAIN_F06_pilot_wgcna.png"), fig,
         width = FIG_W, height = h_mm, units = "mm", dpi = 300, limitsize = FALSE)

  # supp ME-trait heatmap
  me_hm <- ggplot(me_corr, aes(.data$contrast, .data$module, fill = .data$r)) +
    geom_tile(color = "white", linewidth = 0.1) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 1.6, color = "grey15") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
    labs(title = "F06 pilot_wgcna — module eigengene × contrast indicator",
         x = NULL, y = NULL) + FIG_THEME +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = FIG_AXIS_TEXT))
  ggsave(file.path(SUPP_PDF, "MAIN_F06_pilot_wgcna_me_traits.pdf"), me_hm,
         width = 100, height = 6 + 4 * length(unique(me_corr$module)),
         units = "mm", device = pdf_dev, limitsize = FALSE)
  ggsave(file.path(SUPP_PNG, "MAIN_F06_pilot_wgcna_me_traits.png"), me_hm,
         width = 100, height = 6 + 4 * length(unique(me_corr$module)),
         units = "mm", dpi = 300, limitsize = FALSE)

  # ORA across all modules for the workbook — Hallmark + MitoCarta separately
  ora_w <- bind_rows(lapply(mod_order, function(mod) {
    g <- mods$gene[mods$module == mod]
    o <- run_hallmark_ora(g, universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) return(NULL)
    mutate(o, module = mod)
  }))
  ora_mito_w <- bind_rows(lapply(mod_order, function(mod) {
    g <- mods$gene[mods$module == mod]
    o <- run_mitocarta_ora(g, universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) return(NULL)
    mutate(o, module = mod)
  }))

  results$pilot_wgcna <- list(
    key = "pilot_wgcna",
    sheets = list(
      pilot_wgcna_membership = mods,
      pilot_wgcna_me_traits  = signs,
      pilot_wgcna_ora        = if (is.null(ora_w))      tibble() else ora_w,
      pilot_wgcna_ora_mito   = if (is.null(ora_mito_w)) tibble() else ora_mito_w))
} else {
  warning("WGCNA artifact not found at ", WGCNA_RDS, " — skipping pilot_wgcna")
}

if (!is.null(results$pilot_wgcna)) {
  wg <- results$pilot_wgcna$sheets
  sheet_specs <- c(sheet_specs,
    list(list(name = "pilot_wgcna_membership", df = wg$pilot_wgcna_membership,
              role = "Gene -> WGCNA module assignment (grey excluded)",
              contents = "gene, module (color label from F05 build)")),
    list(list(name = "pilot_wgcna_me_traits", df = wg$pilot_wgcna_me_traits,
              role = "Module eigengene Pearson r vs contrast indicator vectors",
              contents = "module, r_Disease, r_Transplant, r_Rescue, r_Interaction, sign_pattern")),
    list(list(name = "pilot_wgcna_ora", df = wg$pilot_wgcna_ora,
              role = "Hallmark ORA per module",
              contents = "module, pathway, padj, overlap, size, odds_ratio")),
    list(list(name = "pilot_wgcna_ora_mito", df = wg$pilot_wgcna_ora_mito,
              role = "MitoCarta ORA per module",
              contents = "module, pathway, padj, overlap, size, odds_ratio")))
  overview2 <- tibble(
    Pilot = "pilot_wgcna", Method = "WGCNA modules (existing artifact)",
    Gate = "all genes (network gate)",
    N_genes = nrow(wg$pilot_wgcna_membership),
    Fuzzifier_m = NA_real_,
    Cluster_c = length(unique(wg$pilot_wgcna_membership$module)))
  sheet_specs[[1]]$df <- bind_rows(sheet_specs[[1]]$df, overview2)
}
