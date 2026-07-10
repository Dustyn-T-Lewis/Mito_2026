# Native contrast-pair comparison panels for F03 (concordance) and F04 (reversal):
# quadrant ORA, fry barcode, NES scatter, RRHO2. Each builder takes
# a wide DE table and an old-name contrast pair and returns list(plot, table).

classify_quadrant <- function(comb, ctr_x, ctr_y, pi_thresh = H9C2_PI_THRESH) {
  lfc_x <- comb[[paste0("logFC_", ctr_x)]]
  lfc_y <- comb[[paste0("logFC_", ctr_y)]]
  sig_x <- comb[[paste0("pi_score_", ctr_x)]] < pi_thresh
  sig_y <- comb[[paste0("pi_score_", ctr_y)]] < pi_thresh
  sig_x[is.na(sig_x)] <- FALSE
  sig_y[is.na(sig_y)] <- FALSE
  quadrant <- dplyr::case_when(
    !(sig_x | sig_y) ~ "ns",
    lfc_x >= 0 & lfc_y >= 0 ~ "x_up_y_up",
    lfc_x >= 0 & lfc_y < 0 ~ "x_up_y_dn",
    lfc_x < 0 & lfc_y >= 0 ~ "x_dn_y_up",
    TRUE ~ "x_dn_y_dn"
  )
  tibble::tibble(
    uniprot_id = comb$uniprot_id, gene = comb$gene,
    lfc_x = lfc_x, lfc_y = lfc_y, sig_x = sig_x, sig_y = sig_y, quadrant = quadrant
  )
}

# RRHO2::RRHO2_initialize splits hypermat on the sign-based boundary
# (count of positive scores at each step index) plus an asymmetric NA strip
# (lenStrip = round(len * boundary)), not on the index midpoint. Re-derive
# the same boundary here so quadrant corners line up with RRHO2's own.
rrho2_quadrant_table <- function(obj, l1, l2) {
  hypermat <- obj$hypermat
  n1 <- nrow(l1)
  n2 <- nrow(l2)
  stepsize <- ceiling(min(sqrt(c(n1, n2))))
  score1 <- sort(l1$score, decreasing = TRUE)
  score2 <- sort(l2$score, decreasing = TRUE)
  step1 <- seq(1, n1, stepsize)
  step2 <- seq(1, n2, stepsize)
  len1 <- length(step1)
  len2 <- length(step2)
  boundary1 <- sum(score1[step1] > 0)
  boundary2 <- sum(score2[step2] > 0)
  strip1 <- round(len1 * 0.1)
  strip2 <- round(len2 * 0.1)

  safe_seq <- function(from, to) if (from > to) integer(0) else seq(from, to)
  corner_max <- function(rows, cols) {
    if (length(rows) == 0 || length(cols) == 0) {
      return(NA_real_)
    }
    block <- hypermat[rows, cols]
    if (all(is.na(block))) {
      return(NA_real_) # all-NA corner (empty quadrant): no cell to report
    }
    max(block, na.rm = TRUE)
  }

  up_rows <- safe_seq(1, boundary1)
  up_cols <- safe_seq(1, boundary2)
  dn_rows <- strip1 + safe_seq(boundary1 + 1, len1)
  dn_cols <- strip2 + safe_seq(boundary2 + 1, len2)

  tibble::tibble(
    quadrant = c("uu", "ud", "du", "dd"),
    max_neglog10p = c(
      corner_max(up_rows, up_cols), corner_max(up_rows, dn_cols),
      corner_max(dn_rows, up_cols), corner_max(dn_rows, dn_cols)
    ),
    n_shared = c(
      length(obj$genelist_uu$gene_list_overlap_uu),
      length(obj$genelist_ud$gene_list_overlap_ud),
      length(obj$genelist_du$gene_list_overlap_du),
      length(obj$genelist_dd$gene_list_overlap_dd)
    )
  )
}

rrho2_render_plot <- function(hypermat, labels, n_total) {
  n_row <- nrow(hypermat)
  n_col <- ncol(hypermat)
  hyper_long <- tibble::tibble(
    row = rep(seq_len(n_row), times = n_col),
    col = rep(seq_len(n_col), each = n_row),
    value = as.vector(hypermat)
  )

  ggplot2::ggplot(hyper_long, ggplot2::aes(col, row, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(
      colours = c("#2166AC", "white", "#B2182B"), na.value = "grey90",
      name = expression(-log[10](italic(p)))
    ) +
    ggplot2::scale_y_reverse(expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::annotate(
      "text",
      x = c(n_col * 0.05, n_col * 0.95, n_col * 0.05, n_col * 0.95),
      y = c(n_row * 0.05, n_row * 0.05, n_row * 0.95, n_row * 0.95),
      label = c("Exacerbated", "Reversed", "Reversed", "Exacerbated"),
      fontface = "bold", size = 2, color = "grey20"
    ) +
    ggplot2::labs(
      x = bquote(.(labels[2]) ~ "rank" ~ (Up %->% Down)),
      y = bquote(.(labels[1]) ~ "rank" ~ (Up %->% Down)),
      subtitle = sprintf(
        paste0(
          "%d shared proteins · warm off-diagonal = transplant reverses ",
          "disease · No MTC (Cahill 2018)"
        ),
        n_total
      )
    ) +
    ggplot2::coord_equal() +
    FIG_THEME +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
}

build_rrho2 <- function(comb, ctr_x, ctr_y, labels = c(ctr_x, ctr_y)) {
  col_x <- paste0("t_", ctr_x)
  col_y <- paste0("t_", ctr_y)
  d <- comb[!is.na(comb[[col_x]]) & !is.na(comb[[col_y]]), ]
  l1 <- data.frame(Genes = d$uniprot_id, score = d[[col_x]])
  l2 <- data.frame(Genes = d$uniprot_id, score = d[[col_y]])

  obj <- RRHO2::RRHO2_initialize(l1, l2, labels = labels, log10.ind = TRUE)
  quad_table <- rrho2_quadrant_table(obj, l1, l2)
  plot <- rrho2_render_plot(obj$hypermat, labels, nrow(l1))

  list(plot = plot, table = quad_table)
}

# Quadrant background tints for the NES scatter: reversal (ref_slope = -1)
# warms the off-diagonal (blue, disease reversed by transplant) and the
# diagonal (red, exacerbated); concordance (ref_slope = +1) warms only the
# diagonal (green, concordant).
nes_quadrant_tints <- function(ref_slope) {
  diag_fill <- if (ref_slope < 0) "#D6604D" else "#4DAF4A"
  quadrants <- tibble::tibble(
    xmin = c(0, -Inf, -Inf, 0),
    xmax = c(Inf, 0, 0, Inf),
    ymin = c(0, 0, -Inf, -Inf),
    ymax = c(Inf, Inf, 0, 0),
    fill = c(diag_fill, "#4393C3", diag_fill, "#4393C3")
  )
  if (ref_slope < 0) {
    return(quadrants)
  }
  quadrants[c(1, 3), ]
}

build_nes_scatter <- function(fgsea_cache, ctr_x, ctr_y, ref_slope = -1) {
  wide <- fgsea_cache |>
    dplyr::filter(contrast %in% c(ctr_x, ctr_y), database %in% CANONICAL_DBS) |>
    dplyr::select(pathway, database, contrast, NES, padj) |>
    tidyr::pivot_wider(names_from = contrast, values_from = c(NES, padj)) |>
    dplyr::rename(
      NES_x = !!paste0("NES_", ctr_x), NES_y = !!paste0("NES_", ctr_y),
      padj_x = !!paste0("padj_", ctr_x), padj_y = !!paste0("padj_", ctr_y)
    ) |>
    tidyr::drop_na(NES_x, NES_y)

  attr(wide, "spearman_rho") <- stats::cor(
    wide$NES_x, wide$NES_y,
    method = "spearman"
  )

  min_padj <- pmin(wide$padj_x, wide$padj_y, na.rm = TRUE)
  sig <- !is.na(min_padj) & min_padj < 0.05
  off_diagonal <- (wide$NES_x * wide$NES_y) < 0
  pct_reversed <- if (any(sig)) 100 * mean(off_diagonal[sig]) else NA_real_

  callouts <- wide |>
    dplyr::mutate(
      diff = abs(NES_x - NES_y),
      label = clean_display_label(pathway)
    ) |>
    dplyr::slice_max(diff, n = 8, with_ties = FALSE)

  plot <- ggplot2::ggplot(wide, ggplot2::aes(NES_x, NES_y)) +
    ggplot2::geom_rect(
      data = nes_quadrant_tints(ref_slope),
      ggplot2::aes(
        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill
      ),
      inherit.aes = FALSE, alpha = 0.15
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_abline(
      slope = ref_slope, intercept = 0, linetype = "dashed", color = "grey40"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = database, size = -log10(pmin(padj_x, padj_y))),
      alpha = 0.7
    ) +
    ggplot2::scale_color_manual(values = ORA_DB_COLORS) +
    (if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(
        data = callouts, ggplot2::aes(label = label),
        size = 2, max.overlaps = Inf, segment.size = 0.2
      )
    } else {
      ggplot2::geom_text(
        data = callouts, ggplot2::aes(label = label),
        size = 2
      )
    }) +
    ggplot2::labs(
      x = sprintf("%s NES", ctr_x), y = sprintf("%s NES", ctr_y),
      size = expression(-log[10](italic(padj))),
      subtitle = sprintf(
        "ρ = %.2f · %.0f%% reversed",
        attr(wide, "spearman_rho"), pct_reversed
      )
    ) +
    FIG_THEME

  list(plot = plot, table = wide)
}

# GSEA-style running enrichment score: a cumulative walk that steps up on a
# set member (weighted by |stat|) and down otherwise, ending at zero once
# every member has been visited. running_es() and .es_barcode() ported
# verbatim from the pilot fry-mountain script (git 8dcf07f,
# 04_Figures/test/pilot_reversal/a_script/panels/e_fry_mountain.R); the
# flanking ORA bars from that pilot are dropped here as redundant with
# Panel A's quadrant ORA.
running_es <- function(t_vals, in_set) {
  n <- length(t_vals)
  n_hit <- sum(in_set)
  if (n_hit == 0) {
    return(rep(0, n))
  }
  hit_w <- ifelse(in_set, abs(t_vals), 0)
  cumsum(ifelse(in_set, hit_w / sum(hit_w), -1 / (n - n_hit)))
}

.es_barcode <- function(t_rank, in_col, es_col, fry_row, title, color, n_all) {
  marks <- t_rank[t_rank[[in_col]], ]
  line_color <- if (isTRUE(fry_row$fdr < 0.05)) {
    color
  } else {
    scales::alpha(color, 0.4)
  }
  lab <- sprintf("fry %s, FDR %.2g", fry_row$dir, fry_row$fdr)
  es <- ggplot2::ggplot(t_rank, ggplot2::aes(rank, .data[[es_col]])) +
    ggplot2::geom_area(fill = scales::alpha(line_color, 0.15), color = NA) +
    ggplot2::geom_line(color = line_color, linewidth = 0.5) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3
    ) +
    ggplot2::annotate("text",
      x = Inf, y = Inf, label = lab, hjust = 1.05, vjust = 1.5,
      size = 1.9, fontface = "bold", color = "grey20"
    ) +
    ggplot2::labs(title = title, x = NULL, y = "ES") +
    ggplot2::scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    FIG_THEME +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 1, 0, 0, "mm")
    )
  bc <- ggplot2::ggplot(
    marks, ggplot2::aes(rank, xend = rank, y = 0, yend = 1)
  ) +
    ggplot2::geom_segment(color = line_color, linewidth = 0.3, alpha = 0.7) +
    ggplot2::scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    FIG_THEME +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "grey97"),
      plot.margin = ggplot2::margin(0, 1, 0, 2, "mm")
    )
  list(es = es, bc = bc)
}

# Design-column algebra for the six canonical DE contrasts (matches the
# add_contrasts() calls in 03_DEP/*/01_run_dep*.R); looked up by ctr_y so
# fry tests the same contrast the reported DE table was fit on.
FRY_CONTRAST_FORMULAS <- c(
  CTLvPHE       = "PHE - Ctl",
  CTLvMITO      = "Mito - Ctl",
  PHEvPHE_MITO  = "PHE_Mito - PHE",
  CTLvPHE_MITO  = "PHE_Mito - Ctl",
  Interaction   = "(PHE_Mito - PHE) - (Mito - Ctl)",
  MITOvPHE_MITO = "PHE_Mito - Mito"
)

build_fry_barcode <- function(comb, expr, design, rep_block, ctr_x, ctr_y,
                              pw_collection = NULL) {
  sig_x <- comb[[paste0("sig_pi_", ctr_x)]]
  up_ids <- comb$uniprot_id[sig_x == 1]
  dn_ids <- comb$uniprot_id[sig_x == -1]

  set_rho <- limma::duplicateCorrelation(
    expr, design,
    block = rep_block
  )$consensus
  cm <- limma::makeContrasts(
    contrasts = FRY_CONTRAST_FORMULAS[[ctr_y]], levels = design
  )

  idx_up <- which(rownames(expr) %in% up_ids)
  idx_dn <- which(rownames(expr) %in% dn_ids)
  set.seed(42)
  fr <- limma::fry(
    expr,
    index = list(up = idx_up, down = idx_dn),
    design = design, contrast = cm[, 1],
    block = rep_block, correlation = set_rho
  )
  fry_up <- tibble::tibble(
    dir = fr["up", "Direction"], fdr = fr["up", "FDR"], n = length(idx_up)
  )
  fry_dn <- tibble::tibble(
    dir = fr["down", "Direction"], fdr = fr["down", "FDR"], n = length(idx_dn)
  )

  t_rank <- comb |>
    dplyr::mutate(t_rank_stat = .data[[paste0("t_", ctr_y)]]) |>
    dplyr::filter(!is.na(t_rank_stat)) |>
    dplyr::transmute(
      uniprot_id, t_rank_stat,
      in_up = uniprot_id %in% up_ids, in_down = uniprot_id %in% dn_ids
    ) |>
    dplyr::arrange(dplyr::desc(t_rank_stat)) |>
    dplyr::mutate(
      rank = dplyr::row_number(),
      es_up = running_es(t_rank_stat, in_up),
      es_down = running_es(t_rank_stat, in_down)
    )
  n_all <- nrow(t_rank)
  color <- unname(CONTRAST_COLORS[ctr_x])

  up_track <- .es_barcode(
    t_rank, "in_up", "es_up", fry_up, "disease-up set", color, n_all
  )
  dn_track <- .es_barcode(
    t_rank, "in_down", "es_down", fry_dn, "disease-down set", color, n_all
  )
  dn_track$bc <- dn_track$bc +
    ggplot2::labs(x = "rank") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(face = "bold", size = 5)
    )

  plot <- patchwork::wrap_plots(
    up_track$es, up_track$bc, dn_track$es, dn_track$bc,
    ncol = 1, heights = c(2, 0.3, 2, 0.3)
  ) +
    patchwork::plot_annotation(
      title = sprintf(
        "fry: %s set across the %s axis",
        contrast_brief(ctr_x), contrast_brief(ctr_y)
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = FIG_TITLE_SIZE, face = "bold")
      )
    )

  table <- tibble::tibble(
    set = c("up", "down"),
    dir = c(fry_up$dir, fry_dn$dir),
    fdr = c(fry_up$fdr, fry_dn$fdr),
    n = c(fry_up$n, fry_dn$n),
    es_end = c(dplyr::last(t_rank$es_up), dplyr::last(t_rank$es_down))
  )

  list(plot = plot, table = table)
}

# Panel A quadrant tint + label: mirrors nes_quadrant_tints()'s reversal
# (ref_slope < 0: diagonal = exacerbated, off-diagonal = reversed) versus
# concordance (ref_slope >= 0: diagonal = concordant, off-diagonal untinted)
# scheme, but keyed by the discrete quadrant name instead of continuous NES.
.ora_quadrant_style <- function(quadrant, ref_slope) {
  diagonal <- quadrant %in% c("x_up_y_up", "x_dn_y_dn")
  if (ref_slope < 0) {
    label <- ifelse(diagonal, "Exacerbated", "Reversed")
    fill <- ifelse(diagonal, "#D6604D", "#4393C3")
  } else {
    label <- ifelse(diagonal, "Concordant", "Discordant")
    fill <- ifelse(diagonal, "#4DAF4A", "white")
  }
  list(label = label, tint = scales::alpha(fill, 0.15))
}

# Horizontal ORA bar for one quadrant: top hits by padj, bars colored by
# database (ORA_DB_COLORS), on a background tinted by the quadrant's
# reversal/concordance meaning. Adapted from the pilot reversal mountain's
# .ora_bars() (04_Figures/test/pilot_reversal/a_script/panels/e_fry_mountain.R,
# git f50b3ee), simplified to a single DB-colored fill instead of a fixed
# up/down color and given a quadrant background tint.
.ora_bars <- function(ora_df, title, tint, n_genes, max_hits = 5) {
  full_title <- sprintf("%s (n = %d)", title, n_genes)
  tinted <- ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = tint, color = NA),
    plot.title = ggplot2::element_text(size = FIG_SUBTITLE_SIZE, face = "bold")
  )

  if (is.null(ora_df) || nrow(ora_df) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0.5, y = 0.5, label = "no enrichment", size = 2.4,
          color = "grey50"
        ) +
        ggplot2::labs(title = full_title) +
        ggplot2::theme_void() +
        tinted
    )
  }

  bars <- ora_df |>
    dplyr::slice_min(padj, n = max_hits, with_ties = FALSE) |>
    dplyr::mutate(
      nlp = -log10(pmax(padj, 1e-20)),
      y = rev(dplyr::row_number()),
      label = clean_display_label(pathway)
    )
  xmax <- max(bars$nlp, na.rm = TRUE)

  ggplot2::ggplot(bars, ggplot2::aes(y = y)) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = 0, xmax = nlp, ymin = y - 0.42, ymax = y + 0.42, fill = database
      ),
      color = "black", linewidth = 0.25
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = xmax * 0.03, label = label),
      hjust = 0, size = 1.7, fontface = "bold", color = "black",
      lineheight = 0.85
    ) +
    ggplot2::scale_fill_manual(values = ORA_DB_COLORS, guide = "none") +
    ggplot2::labs(
      title = full_title, x = expression(-log[10](italic(padj))), y = NULL
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, xmax * 1.18), breaks = scales::pretty_breaks(3),
      expand = ggplot2::expansion(0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.4, nrow(bars) + 0.6), expand = c(0, 0)
    ) +
    FIG_THEME +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    ) +
    tinted
}

# Panel A: per-quadrant over-representation. classify_quadrant() splits the
# contrast pair into its four sign quadrants; each non-"ns" quadrant's
# real-symbol genes are tested against the harmonized pathway collection with
# run_fora_by_db(), then rendered as a 2x2 .ora_bars() grid mirroring the
# NES-scatter layout (top-left = x_dn_y_up, top-right = x_up_y_up,
# bottom-left = x_dn_y_dn, bottom-right = x_up_y_dn).
build_quadrant_ora <- function(comb, ctr_x, ctr_y, universe = NULL,
                               pw_collection = NULL, pw_by_db = NULL,
                               ref_slope = -1) {
  quad <- classify_quadrant(comb, ctr_x, ctr_y)
  quad <- quad[is_real_symbol(quad$gene), ]

  if (is.null(universe)) {
    universe <- quad$gene
  }
  if (is.null(pw_collection)) {
    pw_collection <- build_harmonized_collection()
  }
  if (is.null(pw_by_db)) {
    pw_names <- names(pw_collection)
    pw_by_db <- split(pw_names, classify_database(pw_names))
  }

  quadrants <- c("x_dn_y_up", "x_up_y_up", "x_dn_y_dn", "x_up_y_dn")
  ora_by_quadrant <- sapply(quadrants, function(q) {
    genes <- quad$gene[quad$quadrant == q]
    if (length(genes) == 0) {
      return(NULL)
    }
    res <- run_fora_by_db(genes, universe, pw_collection, pw_by_db)
    if (is.null(res) || nrow(res) == 0) {
      return(NULL)
    }
    res$quadrant <- q
    res[order(res$padj), ]
  }, simplify = FALSE)

  n_genes <- vapply(quadrants, function(q) sum(quad$quadrant == q), integer(1))

  cells <- lapply(quadrants, function(q) {
    style <- .ora_quadrant_style(q, ref_slope)
    .ora_bars(ora_by_quadrant[[q]], style$label, style$tint, n_genes[[q]])
  })

  plot <- patchwork::wrap_plots(
    cells[[1]], cells[[2]], cells[[3]], cells[[4]],
    ncol = 2
  ) +
    patchwork::plot_annotation(
      title = sprintf(
        "Quadrant ORA: %s x %s", contrast_brief(ctr_x), contrast_brief(ctr_y)
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = FIG_TITLE_SIZE, face = "bold")
      )
    )

  table <- dplyr::bind_rows(ora_by_quadrant)

  list(plot = plot, table = table)
}
