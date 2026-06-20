# 04_Figures_v2/functions/07_cluster_row_layout.R
# Shared helpers for the F06 pilot framework:
#   * standardise_genes           — per-gene z-score (finite rows only); used by
#                                   pilots 1-4 (c-means + WGCNA) and pilot 6 (RRHO2)
#   * filter_sig_in_any_contrast  — per-protein significance gate
#   * load_wgcna_modules          — slim wrapper around the F05 WGCNA artifact
#   * compute_me_contrast_correlations / classify_module_sign_pattern
#   * cluster_palette             — viridis::turbo-derived cluster colors
#   * build_trajectory_panel / build_ora_bar_panel / build_cluster_row
#   * stack_cluster_rows          — patchwork vertical stack
#   * run_hallmark_ora            — Hallmark-only ORA (fora) with no dedup

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(patchwork); library(viridis)
})

# Per-gene z-score across conditions; drops rows that are all-NA or all-constant
# after scaling (produces non-finite row sums). Used by pilots 1-4 and 6.
standardise_genes <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.finite(rowSums(z)), , drop = FALSE]
}

filter_sig_in_any_contrast <- function(comb_long, col, threshold,
                                       contrasts, op = c("lt", "le")) {
  op <- match.arg(op)
  stopifnot(col %in% names(comb_long), "gene" %in% names(comb_long))
  d <- comb_long[comb_long$contrast %in% contrasts, , drop = FALSE]
  v <- d[[col]]
  hit <- if (op == "lt") !is.na(v) & v < threshold
         else            !is.na(v) & v <= threshold
  unique(d$gene[hit & !is.na(d$gene) & nzchar(d$gene)])
}

load_wgcna_modules <- function(rds_path) {
  obj <- readRDS(rds_path)
  # The build_wgcna_network.R artifact carries either a flat list of
  # (gene, module, MEs, ...) or a nested 'net' element; handle both.
  # Fallback 1: explicit module_assignments tibble (ideal).
  # Fallback 2: net$colors named by uniprot IDs + ann$gene for symbols.
  # Fallback 3: module_colors positional vector + ann$gene for symbols.
  if (!is.null(obj$module_assignments)) {
    mods <- obj$module_assignments
  } else if (!is.null(obj$net) && !is.null(obj$net$colors) &&
             length(names(obj$net$colors)) > 0) {
    # Real artifact: net$colors is named by uniprot IDs (numeric cluster IDs);
    # module_colors (when present) holds the parallel color-name strings;
    # ann$gene gives gene symbols parallel to both.
    gene_syms <- if (!is.null(obj$ann) && !is.null(obj$ann$gene)) {
      obj$ann$gene
    } else {
      names(obj$net$colors)
    }
    color_names <- if (!is.null(obj$module_colors) &&
                       length(obj$module_colors) == length(gene_syms)) {
      obj$module_colors
    } else {
      as.character(unname(obj$net$colors))
    }
    mods <- tibble(gene   = gene_syms,
                   module = color_names)
  } else if (!is.null(obj$module_colors)) {
    # Positional fallback: module_colors parallel to ann$gene or gene_order.
    gene_syms <- if (!is.null(obj$ann) && !is.null(obj$ann$gene)) {
      obj$ann$gene
    } else {
      obj$gene_order %||% seq_along(obj$module_colors)
    }
    mods <- tibble(gene   = gene_syms,
                   module = obj$module_colors)
  } else {
    stop("Unrecognised WGCNA artifact shape: cannot find module assignments")
  }

  MEs <- if (!is.null(obj$MEs)) as.matrix(obj$MEs)
         else if (!is.null(obj$net) && !is.null(obj$net$MEs)) as.matrix(obj$net$MEs)
         else stop("Unrecognised WGCNA artifact shape: cannot find module eigengenes")

  # Strip "ME" prefix from eigengene column names so they match module color labels.
  colnames(MEs) <- sub("^ME", "", colnames(MEs))

  # Module names are colors; the color lookup is the identity unless the
  # artifact already provides a translation.
  color_lookup <- setNames(unique(mods$module), unique(mods$module))
  list(modules = tibble::as_tibble(mods), MEs = MEs, color_lookup = color_lookup)
}

# Pearson r + Student p for each ME column against each binary contrast
# indicator. `meta` must carry a Group column matching MEs rownames; contrasts
# are interpreted as named lists list("Disease" = c("Ctl","PHE"), ...) where
# the second level is the +1 group and the first is the -1 group.
compute_me_contrast_correlations <- function(MEs, meta, contrasts) {
  stopifnot(nrow(MEs) == nrow(meta))
  stopifnot("compute_me_contrast_correlations: meta$Group is all NA" = any(!is.na(meta$Group)))
  out <- list()
  for (cn in names(contrasts)) {
    pair <- contrasts[[cn]]
    indic <- ifelse(meta$Group == pair[2],  1,
              ifelse(meta$Group == pair[1], -1, NA_real_))
    keep <- !is.na(indic)
    if (sum(keep) < 4) next
    for (mod in colnames(MEs)) {
      ct <- suppressWarnings(cor.test(MEs[keep, mod], indic[keep], method = "pearson"))
      out[[length(out) + 1]] <- tibble(module = mod, contrast = cn,
                                       r = unname(ct$estimate),
                                       p = ct$p.value)
    }
  }
  bind_rows(out)
}

classify_module_sign_pattern <- function(disease_r, rescue_r, r_floor = 0.4) {
  cls <- dplyr::case_when(
    abs(disease_r) >= r_floor & abs(rescue_r) >= r_floor &
      sign(disease_r) != sign(rescue_r)              ~ "Reversal",
    disease_r >=  r_floor & rescue_r >=  r_floor     ~ "Concordant up",
    disease_r <= -r_floor & rescue_r <= -r_floor     ~ "Concordant down",
    TRUE                                             ~ "Other"
  )
  factor(cls, levels = c("Reversal", "Concordant up", "Concordant down", "Other"))
}

cluster_palette <- function(n) {
  stopifnot(n >= 1)
  pal <- viridis::turbo(n + 2)[seq_len(n + 1)[-1]]   # trim extremes
  setNames(pal, as.character(seq_len(n)))
}

build_trajectory_panel <- function(z_mat, x_levels, x_lab,
                                   color, kind = c("line", "barlogfc")) {
  kind <- match.arg(kind)
  z_df <- as.data.frame(z_mat) |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(cols = colnames(z_mat), names_to = "x", values_to = "expr") |>
    dplyr::mutate(x = factor(.data$x, levels = x_levels))
  if (kind == "line") {
    p <- ggplot(z_df, aes(x, expr, group = gene)) +
      geom_line(color = color, alpha = 0.30, linewidth = 0.25) +
      stat_summary(aes(group = 1), fun = mean, geom = "line",
                   color = color, linewidth = 0.9)
  } else {
    means <- z_df |> summarise(mean_expr = mean(.data$expr), .by = x)
    p <- ggplot(means, aes(x, mean_expr)) +
      geom_col(fill = color, color = "grey20", linewidth = 0.25, width = 0.7) +
      geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3)
  }
  p +
    labs(x = x_lab, y = if (kind == "line") "z" else "mean logFC") +
    FIG_THEME +
    theme(plot.margin = margin(2, 2, 2, 2),
          axis.text.x = element_text(angle = 30, hjust = 1, size = FIG_AXIS_TEXT))
}

build_ora_bar_panel <- function(ora_df, color, max_n = 6,
                                 db_name = "Hallmark") {
  empty <- function() ggplot() +
    annotate("text", 0, 0, label = sprintf("no %s hits", db_name),
             size = 2.4, color = "grey40") +
    theme_void()
  if (is.null(ora_df) || nrow(ora_df) == 0) return(empty())
  d <- ora_df |>
    filter(.data$padj < 0.05) |>
    arrange(.data$padj) |>
    head(max_n) |>
    mutate(label = clean_display_label(.data$pathway),
           neglog10 = -log10(.data$padj))
  if (nrow(d) == 0) return(empty())
  ggplot(d, aes(reorder(.data$label, .data$neglog10), .data$neglog10)) +
    geom_col(fill = color, color = "grey20", linewidth = 0.2, width = 0.78) +
    coord_flip() +
    labs(x = NULL, y = sprintf("%s | -log10 padj", db_name)) +
    FIG_THEME +
    theme(plot.margin = margin(2, 4, 2, 2),
          axis.text.y = element_text(size = FIG_AXIS_TEXT - 0.5),
          axis.title.x = element_text(size = FIG_AXIS_TEXT - 0.5,
                                      face = "bold", color = "grey25"))
}

# Hallmark + MitoCarta side-by-side (Reimand 2019 PMID 30664679: pair a
# curated high-level collection with a focused domain-specific resource).
# When ora_plot2 = NULL the row degrades to the original single-ORA layout.
build_cluster_row <- function(traj_plot, ora_plot, header_text, color,
                              ora_plot2 = NULL, widths = NULL) {
  header_color <- ifelse(is_light_color(color), "grey15", color)
  header <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = header_text,
             hjust = 0, vjust = 0.5, size = 2.3, color = header_color, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_void()
  if (is.null(ora_plot2)) {
    if (is.null(widths)) widths <- c(1, 1.4)
    body <- traj_plot + ora_plot + patchwork::plot_layout(widths = widths)
  } else {
    if (is.null(widths)) widths <- c(1, 1, 1)
    body <- traj_plot + ora_plot + ora_plot2 + patchwork::plot_layout(widths = widths)
  }
  header / body + patchwork::plot_layout(heights = c(0.10, 1))
}

# MitoCarta-only ORA: same Fisher 2x2 odds-ratio computation as run_hallmark_ora,
# but restricted to the MitoCarta sublist of rat_gene_sets.rds with the broad
# compartment aggregates (MITO_DROP_SETS) removed so only real MitoPathways
# survive. Pairs with Hallmark per Reimand 2019 PMID 30664679 §"two complementary
# sources"; MitoCarta 3.0 (Rath PMID 33174596) resolves OXPHOS by complex,
# mito-translation, mitophagy etc. — the granularity a mito-transplant figure needs.
run_mitocarta_ora <- function(genes, universe,
                              rat_gene_sets_path = here::here(
                                "04_Figures", "shared", "rat_gene_sets.rds")) {
  gs   <- readRDS(rat_gene_sets_path)
  mito <- gs$MitoCarta
  if (length(mito) == 0) return(NULL)
  mito <- mito[!names(mito) %in% MITO_DROP_SETS]
  if (length(mito) == 0) return(NULL)
  res <- fgsea::fora(pathways = mito, genes = genes, universe = universe,
                     minSize = 5, maxSize = 500)
  res <- as.data.frame(res)
  if (nrow(res) == 0) return(NULL)
  K <- length(intersect(genes, universe))
  N <- length(universe)
  res$odds_ratio <- vapply(seq_len(nrow(res)), function(i) {
    a <- res$overlap[i]; b <- K - a
    c <- res$size[i] - a; d <- N - K - c
    if (b <= 0 || c <= 0) Inf else (a * d) / (b * c)
  }, numeric(1))
  tibble::as_tibble(res)
}

stack_cluster_rows <- function(rows_list, title, subtitle) {
  if (length(rows_list) == 0)
    return(ggplot() + annotate("text", 0, 0, label = "no clusters") + theme_void())
  patchwork::wrap_plots(rows_list, ncol = 1) +
    patchwork::plot_annotation(title = title, subtitle = subtitle,
      theme = theme(plot.title    = element_text(face = "bold", size = 8),
                    plot.subtitle = element_text(face = "italic", size = 5,
                                                 color = "grey30")))
}

# Hallmark-only ORA: fora over the Hallmark sublist of rat_gene_sets.rds. No
# dedup needed at the Hallmark level (50 sets, low redundancy).
run_hallmark_ora <- function(genes, universe,
                             rat_gene_sets_path = here::here("04_Figures", "shared",
                                                             "rat_gene_sets.rds")) {
  gs   <- readRDS(rat_gene_sets_path)
  hall <- gs$Hallmark
  if (length(hall) == 0) return(NULL)
  res <- fgsea::fora(pathways = hall, genes = genes, universe = universe,
                     minSize = 5, maxSize = 500)
  res <- as.data.frame(res)
  if (nrow(res) == 0) return(NULL)
  K <- length(intersect(genes, universe))
  N <- length(universe)
  res$odds_ratio <- vapply(seq_len(nrow(res)), function(i) {
    a <- res$overlap[i]; b <- K - a
    c <- res$size[i] - a; d <- N - K - c
    if (b <= 0 || c <= 0) Inf else (a * d) / (b * c)
  }, numeric(1))
  tibble::as_tibble(res)
}

# %||% is defined in 06_supplementary_workbook.R (null/length-0/NA-aware);
# that helper is sourced transitively in F06 before this file is used.
