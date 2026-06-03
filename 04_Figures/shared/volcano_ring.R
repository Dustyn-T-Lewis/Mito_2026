#!/usr/bin/env Rscript
# volcano_ring.R — Circular volcano-in-ring composite plot utility
# Standard Cartesian ggplot with ggforce::geom_arc_bar(); NO coord_polar().

library(tidyverse)
library(ggforce)
library(patchwork)
library(scales)

if (!exists("FIG_THEME")) source(here::here("04_Figures", "shared", "style.R"))

clean_ring_label <- function(name) {
  # Vectorize over vector input (called from mutate() with a column).
  if (length(name) > 1) {
    return(vapply(name, clean_ring_label, character(1), USE.NAMES = FALSE))
  }
  # MitoCarta-specific shortening: take the leaf of the MitoPathway hierarchy
  # and apply mito-specific capitalizations. Prevents "MitoCarta Mito.
  # Metabolism Carb. Metab. TCA" path-clutter on F03_mito rings.
  if (startsWith(name, "MITOCARTA_")) {
    n <- sub("^MITOCARTA_", "", name)
    parts <- strsplit(n, "__")[[1]]
    leaf <- tail(parts, 1)
    # If the leaf is generic (e.g. "ASSEMBLY_FACTORS"), prefix with parent context.
    # Skip when the leaf is C[IV]+_SUBUNITS / C[IV]+_ASSEMBLY_FACTORS (e.g.
    # CI_SUBUNITS, CV_ASSEMBLY_FACTORS) — the explicit replacements below
    # ("Cv Subunits" -> "Complex V") already encode the complex identity.
    if (length(parts) >= 2 &&
        grepl("SUBUNITS$|FACTORS$|ASSOCIATED$|COMPLEX$|FAMILY$",
              leaf, perl = TRUE) &&
        nchar(leaf) <= 18 &&
        !grepl("^C[IV]+_(SUBUNITS|ASSEMBLY_FACTORS)$", leaf)) {
      leaf <- paste(parts[length(parts) - 1], leaf, sep = " ")
    }
    leaf <- gsub("_", " ", leaf)
    leaf <- tools::toTitleCase(tolower(leaf))
    leaf <- gsub("\\bOxphos\\b",          "OXPHOS",      leaf)
    leaf <- gsub("\\bTca\\b",             "TCA",         leaf)
    leaf <- gsub("\\bAa\\b",              "AA",          leaf)
    leaf <- gsub("\\bFa\\b",              "FA",          leaf)
    leaf <- gsub("\\bRos\\b",             "ROS",         leaf)
    leaf <- gsub("\\bImm\\b",             "IMM",         leaf)
    leaf <- gsub("\\bOmm\\b",             "OMM",         leaf)
    leaf <- gsub("\\bIms\\b",             "IMS",         leaf)
    leaf <- gsub("\\bCi Subunits\\b",     "Complex I",   leaf)
    leaf <- gsub("\\bCii Subunits\\b",    "Complex II",  leaf)
    leaf <- gsub("\\bCiii Subunits\\b",   "Complex III", leaf)
    leaf <- gsub("\\bCiv Subunits\\b",    "Complex IV",  leaf)
    leaf <- gsub("\\bCv Subunits\\b",     "Complex V",   leaf)
    leaf <- gsub("\\bCi Assembly Factors\\b",   "CI assembly",   leaf)
    leaf <- gsub("\\bCiii Assembly Factors\\b", "CIII assembly", leaf)
    leaf <- gsub("\\bCiv Assembly Factors\\b",  "CIV assembly",  leaf)
    leaf <- gsub("\\bCv Assembly Factors\\b",   "CV assembly",   leaf)
    leaf <- gsub("\\bSlc25a Family\\b",   "SLC25A family",       leaf)
    leaf <- gsub("\\bMitochondrial\\b",   "Mito.",               leaf)
    leaf <- gsub("\\bMitochondrion\\b",   "Mito.",               leaf)
    leaf <- gsub("\\bAmino Acid\\b",      "AA",                  leaf)
    leaf <- gsub("\\bFatty Acid\\b",      "FA",                  leaf)
    # Collapse adjacent duplicate words: "OXPHOS OXPHOS Subunits" -> "OXPHOS Subunits",
    # "Complex V Complex V" -> "Complex V" (handles parent-context + expansion overlap).
    repeat {
      new_leaf <- gsub("\\b([A-Za-z][A-Za-z0-9.]*)\\b\\s+\\1\\b",
                       "\\1", leaf, perl = TRUE, ignore.case = TRUE)
      if (identical(new_leaf, leaf)) break
      leaf <- new_leaf
    }
    leaf <- gsub("\\s+", " ", leaf)
    leaf <- trimws(leaf)
    return(stringr::str_wrap(leaf, width = 15))
  }

  name |>
    clean_pathway_name() |>
    # Scientific acronyms (mostly added 2026-05-15 to fix Mito F03 artifacts)
    str_replace_all("\\bTca\\b", "TCA") |>
    str_replace_all("\\bMapk(\\d?)\\b", "MAPK\\1") |>
    str_replace_all("\\bPtk(\\d?)\\b", "PTK\\1") |>
    str_replace_all("\\bRhoa\\b", "RhoA") |>
    str_replace_all("\\bRhob\\b", "RhoB") |>
    str_replace_all("\\bRhoc\\b", "RhoC") |>
    str_replace_all("\\bGtpase\\b", "GTPase") |>
    str_replace_all("\\bGtp\\b", "GTP") |>
    str_replace_all("\\bGdp\\b", "GDP") |>
    str_replace_all("\\bRnas\\b", "RNAs") |>
    str_replace_all("\\bMrna\\b", "mRNA") |>
    str_replace_all("\\bPirnas?\\b", "piRNAs") |>
    str_replace_all("\\bDgc\\b", "DGC") |>
    str_replace_all("\\bMpc\\b", "MPC") |>
    str_replace_all("\\bAtp\\b", "ATP") |>
    str_replace_all("\\bNadh\\b", "NADH") |>
    str_replace_all("\\bFadh\\b", "FADH") |>
    str_replace_all("\\bRos\\b", "ROS") |>
    str_replace_all("\\bDna\\b", "DNA") |>
    str_replace_all("\\bRna\\b", "RNA") |>
    str_replace_all("\\bIfn\\b", "IFN") |>
    str_replace_all("\\bIl-?(\\d+)\\b", "IL\\1") |>
    # Phrase-level cleanups
    str_replace("Trna ", "tRNA ") |>
    str_replace("Pi3k", "PI3K") |>
    str_replace("Akt", "AKT") |>
    str_replace("Mtor", "mTOR") |>
    str_replace("Unfolded Protein Response", "UPR") |>
    str_replace("Fatty Acid", "FA") |>
    str_replace("Amino Acid", "AA") |>
    str_replace("Generation Of", "Gen. of") |>
    str_replace(" And ", " & ") |>
    str_replace("Precursor Metabolites & Energy", "Precursor Metabolites") |>
    str_replace("Mitochondrion", "Mito.") |>
    str_replace("Mitochondrial", "Mito.") |>
    str_replace("Ubiquinone", "UQ") |>
    str_replace("Organization", "Org.") |>
    str_replace("Cytoskeleton", "Cytoskel.") |>
    str_replace("Microtubule", "MT") |>
    str_replace("Respiratory", "Resp.") |>
    str_replace("Electron Transport", "ETC") |>
    str_replace("Synthesis Coupled", "Synth.-Coupled") |>
    str_replace("Ubiquitin Dependent", "Ub-Dep.") |>
    str_replace("Proteasome Mediated", "Proteasome-Med.") |>
    str_replace("Proteasomal", "Proteas.") |>
    str_replace("Phosphorylation", "Phosph.") |>
    str_replace("Modification", "Mod.") |>
    str_replace("Intracellular", "Intracell.") |>
    str_replace("Regulation Of", "Reg.") |>
    str_replace("Signaling Pathway", "Signaling") |>
    str_replace("Biosynthetic Process", "Biosynthesis") |>
    str_replace("Catabolic Process", "Catabolism") |>
    str_replace("Metabolic Process", "Metabolism") |>
    str_replace("Based Process", "Process") |>
    str_replace("Response To", "Resp. to") |>
    str_replace("Extracellular Matrix", "ECM") |>
    str_replace("Epithelial Mesenchymal Transition", "EMT") |>
    str_replace("Establishment Or Maintenance Of", "Maintenance of") |>
    # Redundancies created by upstream MSigDB hierarchies
    str_replace("Citric Acid Cycle TCA Cycle", "TCA Cycle") |>
    str_replace("Aerobic Respiration & Resp. ETC", "Aerobic Resp. + ETC") |>
    str_replace("Reg\\. Endogenous Retroelements By Piwi Interacting RNAs piRNAs",
                "piRNA-Mediated Retroelement Reg.") |>
    str_replace("Formation Of The Dystrophin Glycoprotein Complex DGC",
                "Dystrophin Glycoprotein Complex") |>
    str_wrap(width = 15) |>
    # manual overrides post-wrap
    str_replace(fixed("Protein\nLocalization To\nPlasma Membrane"),
                "Protein Localiz.\nto Plasma\nMem.") |>
    str_replace("(?s).*Maintenance.*Cell.*Polarity.*", "Maintenance\nof Polarity") |>
    str_replace("^Heme Metabolism$",    "Heme\nMetabolism") |>
    str_replace("^tRNA Metabolism$",    "tRNA\nMetabolism") |>
    str_replace("^Mitotic Spindle$",    "Mitotic\nSpindle") |>
    str_replace("^MYC Targets V1$",     "MYC Targets\nV1") |>
    str_replace("^MYC Targets\nV1$",    "MYC Targets\nV1") |>
    str_replace("^UV Resp\\. to Dn$",   "UV Response\nDn") |>
    str_replace("^UV Response Dn$",     "UV Response\nDn") |>
    str_trim()
}

prepare_ring_data <- function(go_df,
                              contrast,
                              n_terms      = 12,
                              gap_degrees  = 3,
                              start_offset = 0,
                              databases    = c("Hallmark", "GO Slim")) {

  ring <- go_df |>
    filter(contrast == !!contrast,
           database %in% databases,
           padj < 0.05) |>
    arrange(padj) |>
    slice_head(n = n_terms)

  n <- nrow(ring)
  if (n == 0) {
    warning("prepare_ring_data: no significant terms for '", contrast, "'")
    return(tibble())
  }

  arc_width_deg <- (360 - n * gap_degrees) / n

  ring |>
    mutate(
      term_idx    = row_number(),
      start_deg   = start_offset + (term_idx - 1) * (arc_width_deg + gap_degrees),
      end_deg     = start_deg + arc_width_deg,
      mid_deg     = (start_deg + end_deg) / 2,
      start_rad   = start_deg * pi / 180,
      end_rad     = end_deg   * pi / 180,
      mid_rad     = mid_deg   * pi / 180,
      arc_r1_var  = 5.6,
      clean_label = clean_ring_label(pathway),
      gene_list   = str_split(leadingEdge, ";")
    )
}

build_tick_data <- function(ring_data,
                            de_df,
                            contrast,
                            tick_r0 = 4.4,
                            tick_r1 = 4.8) {

  if (nrow(ring_data) == 0) return(tibble())

  logfc_col <- paste0("logFC_", contrast)

  gene_lfc <- de_df |>
    dplyr::select(gene, lfc = all_of(logfc_col)) |>
    filter(!is.na(lfc)) |>
    distinct(gene, .keep_all = TRUE)

  pad_rad <- 0.5 * pi / 180

  map_dfr(seq_len(nrow(ring_data)), function(i) {
    row <- ring_data[i, ]
    genes_in_arc <- intersect(row$gene_list[[1]], gene_lfc$gene)
    n_genes <- length(genes_in_arc)
    if (n_genes == 0) return(tibble())

    arc_start <- row$start_rad + pad_rad
    arc_end   <- row$end_rad   - pad_rad
    if (arc_end <= arc_start) arc_end <- arc_start + pad_rad

    tick_angles <- seq(arc_start, arc_end, length.out = n_genes)

    matched_tmp <- gene_lfc |>
      filter(gene %in% genes_in_arc)
    matched <- matched_tmp[match(genes_in_arc, matched_tmp$gene), ] |>
      filter(!is.na(gene))

    n_final <- nrow(matched)
    if (n_final == 0) return(tibble())
    tick_angles <- tick_angles[seq_len(n_final)]

    tibble(
      gene      = matched$gene,
      logFC     = matched$lfc,
      direction = ifelse(matched$lfc > 0, "Up", "Down"),
      angle_rad = tick_angles,
      x0        = tick_r0 * sin(tick_angles),
      y0        = tick_r0 * cos(tick_angles),
      x1        = tick_r1 * sin(tick_angles),
      y1        = tick_r1 * cos(tick_angles),
      term_idx  = row$term_idx,
      pathway   = row$pathway
    )
  })
}

build_volcano_layers <- function(de_df,
                                 contrast,
                                 volcano_radius = 3.5,
                                 fc_thresh      = log2(1.5),
                                 p_thresh       = 0.05,
                                 up_color       = DIR_COLORS["Up"],
                                 down_color     = DIR_COLORS["Down"],
                                 ns_color       = DIR_COLORS["NS"],
                                 point_size     = 0.6,
                                 point_alpha    = 0.5,
                                 count_label_size    = 2.8,
                                 count_label_padding = 3,
                                 count_border_width  = 0.4,
                                 count_y_mult        = 1.0,
                                 count_x_mult        = 0.5) {

  logfc_col <- paste0("logFC_", contrast)
  pval_col  <- paste0("P.Value_", contrast)
  pi_col    <- paste0("pi_score_", contrast)

  vdf <- de_df |>
    transmute(
      gene       = gene,
      logFC      = .data[[logfc_col]],
      pvalue     = .data[[pval_col]],
      pi_score   = .data[[pi_col]],
      neg_log10p = -log10(pvalue)
    ) |>
    filter(!is.na(logFC), !is.na(pvalue), is.finite(neg_log10p)) |>
    mutate(
      direction = case_when(
        pi_score < H9C2_PI_THRESH & abs(logFC) > fc_thresh & logFC > 0 ~ "Up",
        pi_score < H9C2_PI_THRESH & abs(logFC) > fc_thresh & logFC < 0 ~ "Down",
        TRUE                                                            ~ "NS"
      )
    )

  n_up   <- sum(vdf$direction == "Up")
  n_down <- sum(vdf$direction == "Down")

  x_data_max <- max(abs(vdf$logFC), na.rm = TRUE)
  y_data_max <- max(vdf$neg_log10p, na.rm = TRUE)

  margin <- 0.92
  vr <- volcano_radius * margin

  vdf <- vdf |>
    mutate(
      x_plot = logFC / x_data_max * vr,
      y_plot = (neg_log10p / y_data_max) * 2 * vr - vr
    )

  vdf_ns  <- vdf  |> filter(direction == "NS")
  vdf_sig <- vdf  |> filter(direction != "NS")

  layers <- list(
    ns_points = geom_point(
      data = vdf_ns, aes(x = x_plot, y = y_plot),
      color = ns_color, size = point_size * 0.8, alpha = point_alpha * 0.35,
      inherit.aes = FALSE
    ),
    sig_points = geom_point(
      data = vdf_sig, aes(x = x_plot, y = y_plot, color = direction),
      size = point_size * 1.4, alpha = point_alpha * 1.4, stroke = 0.3,
      inherit.aes = FALSE
    ),
    color_scale = scale_color_manual(
      values = c(Up = unname(up_color), Down = unname(down_color)),
      guide  = "none"
    ),
    x_axis_line = annotate(
      "segment", x = -vr * 0.42, xend = vr * 0.42, y = -vr, yend = -vr,
      linewidth = 0.3, linetype = "dashed", color = "grey50",
      arrow = arrow(ends = "both", length = unit(1.2, "mm"), type = "closed")
    ),
    x_axis_up = annotate(
      "text", x = vr * 0.45, y = -vr,
      label = "up", size = count_label_size * 0.9,
      color = unname(up_color), fontface = "bold.italic", hjust = 0
    ),
    x_axis_down = annotate(
      "text", x = -vr * 0.45, y = -vr,
      label = "down", size = count_label_size * 0.9,
      color = unname(down_color), fontface = "bold.italic", hjust = 1
    ),
    x_axis_label = annotate(
      "text", x = 0, y = -vr - 0.35,
      label = "log2 FC", size = count_label_size * 0.9,
      color = "grey40", fontface = "bold.italic"
    ),
    y_axis_line = annotate(
      "segment", x = 0, xend = 0, y = -vr, yend = vr * 0.96,
      linewidth = 0.3, linetype = "dashed", color = "grey50",
      arrow = arrow(ends = "last", length = unit(1.2, "mm"), type = "closed")
    ),
    y_axis_label = annotate(
      "text", x = 0, y = vr * 1.04,
      label = expression(bold(-log[10])~bolditalic(p)), size = count_label_size * 0.9,
      color = "grey40"
    ),
    n_up_box = annotate(
      "label", x = vr * count_x_mult, y = vr * count_y_mult,
      label = n_up, size = count_label_size * 1.25,
      color = "black", fill = alpha(up_color, 0.9), fontface = "bold",
      label.padding = unit(count_label_padding, "pt"), label.r = unit(2, "pt"),
      linewidth = count_border_width
    ),
    n_up_text = annotate(
      "text", x = vr * count_x_mult, y = vr * count_y_mult,
      label = n_up, size = count_label_size * 1.25,
      color = "white", fontface = "bold"
    ),
    n_down_box = annotate(
      "label", x = -vr * count_x_mult, y = vr * count_y_mult,
      label = n_down, size = count_label_size * 1.25,
      color = "black", fill = alpha(down_color, 0.9), fontface = "bold",
      label.padding = unit(count_label_padding, "pt"), label.r = unit(2, "pt"),
      linewidth = count_border_width
    ),
    n_down_text = annotate(
      "text", x = -vr * count_x_mult, y = vr * count_y_mult,
      label = n_down, size = count_label_size * 1.25,
      color = "white", fontface = "bold"
    )
  )

  attr(layers, "x_data_max") <- x_data_max
  attr(layers, "y_data_max") <- y_data_max
  attr(layers, "n_up")       <- n_up
  attr(layers, "n_down")     <- n_down

  layers
}

build_ring_layers <- function(ring_data,
                              tick_data,
                              tick_r0    = 4.4,
                              tick_r1    = 4.8,
                              arc_r0     = 4.8,
                              arc_r1     = 5.6,
                              up_color   = DIR_COLORS["Up"],
                              down_color = DIR_COLORS["Down"],
                              mito_pattern = NULL) {
  # mito_pattern: optional regex (perl-style). When provided, arcs whose
  # `pathway` name matches the regex OR starts with "MITOCARTA_" get a second
  # bold-dark outline layer drawn on top of the standard NES-filled arc.
  # Visual cue: "this term is mitochondrial biology".

  if (nrow(ring_data) == 0) return(list())

  layers <- list()

  layers$tick_bg <- geom_arc_bar(
    data = ring_data,
    aes(x0 = 0, y0 = 0, r0 = tick_r0, r = tick_r1,
        start = start_rad, end = end_rad),
    fill = "grey93", color = "grey78", linewidth = 0.15,
    inherit.aes = FALSE
  )

  if (nrow(tick_data) > 0) {
    layers$ticks <- geom_segment(
      data = tick_data,
      aes(x = x0, y = y0, xend = x1, yend = y1, color = direction),
      linewidth = 0.15, alpha = 0.7,
      inherit.aes = FALSE
    )
  }

  layers$enrich_arcs <- geom_arc_bar(
    data = ring_data,
    aes(x0 = 0, y0 = 0, r0 = arc_r0, r = arc_r1_var,
        start = start_rad, end = end_rad, fill = NES),
    color = "grey40", linewidth = 0.2,
    inherit.aes = FALSE
  )

  layers$fill_scale <- scale_fill_gradientn(
    colours = c("#08306B", "#4393C3", "white", "#D6604D", "#67000D"),
    values  = scales::rescale(c(-3, -1.5, 0, 1.5, 3)),
    limits  = c(-3, 3),
    oob     = scales::squish,
    name    = "NES"
  )

  # Optional second layer: bold dark outline on mito-relevant arcs.
  if (!is.null(mito_pattern) && "pathway" %in% names(ring_data)) {
    is_mito <- grepl(mito_pattern, ring_data$pathway, perl = TRUE) |
               startsWith(ring_data$pathway, "MITOCARTA_")
    mito_arcs <- ring_data[is_mito, , drop = FALSE]
    if (nrow(mito_arcs) > 0) {
      layers$mito_outline <- geom_arc_bar(
        data = mito_arcs,
        aes(x0 = 0, y0 = 0, r0 = arc_r0, r = arc_r1_var,
            start = start_rad, end = end_rad),
        fill = NA, color = "grey10", linewidth = 0.8,
        inherit.aes = FALSE
      )
    }
  }

  layers
}

build_label_layer <- function(ring_data,
                              label_r    = 7.0,
                              label_size = 3.0,
                              label_gap  = NULL,
                              label_padding = 2,
                              min_angle_gap = 18,
                              nudge_outward = 0.8,
                              up_color   = DIR_COLORS["Up"],
                              down_color = DIR_COLORS["Down"]) {

  if (nrow(ring_data) == 0) return(list())

  lbl_df <- ring_data |>
    mutate(
      label_r_term = if (!is.null(label_gap)) arc_r1_var + label_gap else label_r
    )

  # nudge close labels outward
  if (nrow(lbl_df) >= 2) {
    lbl_df <- lbl_df  |> arrange(mid_deg)
    for (i in 2:nrow(lbl_df)) {
      if (abs(lbl_df$mid_deg[i] - lbl_df$mid_deg[i - 1]) < min_angle_gap) {
        lbl_df$label_r_term[i] <- lbl_df$label_r_term[i] + nudge_outward
      }
    }
    wrap_gap <- 360 - lbl_df$mid_deg[nrow(lbl_df)] + lbl_df$mid_deg[1]
    if (wrap_gap < min_angle_gap) {
      lbl_df$label_r_term[1] <- lbl_df$label_r_term[1] + 0.8
    }
  }

  lbl_df <- lbl_df |>
    mutate(
      lbl_x      = label_r_term * sin(mid_rad),
      lbl_y      = label_r_term * cos(mid_rad),
      lead_x     = (arc_r1_var + 0.1) * sin(mid_rad),
      lead_y     = (arc_r1_var + 0.1) * cos(mid_rad),
      lead_end_x = (label_r_term - 0.25) * sin(mid_rad),
      lead_end_y = (label_r_term - 0.25) * cos(mid_rad),
      side_x     = sin(mid_rad),
      lbl_hjust  = case_when(
        side_x >  0.15 ~ 0,
        side_x < -0.15 ~ 1,
        TRUE           ~ 0.5
      ),
      nudge_x = case_when(
        side_x >  0.15 ~  0.6,
        side_x < -0.15 ~ -0.6,
        TRUE           ~  0
      ),
      legend_label = clean_label
    )

  attr(lbl_df, "max_label_r") <- max(lbl_df$label_r_term)

  up_lbl   <- lbl_df  |> filter(NES > 0)
  down_lbl <- lbl_df  |> filter(NES <= 0)

  layers <- list()

  layers$leaders <- geom_segment(
    data = lbl_df,
    aes(x = lead_x, y = lead_y, xend = lead_end_x, yend = lead_end_y),
    linewidth = 0.5, color = "grey35",
    inherit.aes = FALSE
  )

  if (nrow(up_lbl) > 0) {
    layers$up_labels <- geom_label(
      data = up_lbl,
      aes(x = lbl_x + nudge_x, y = lbl_y, label = legend_label),
      hjust = 0.5, vjust = 0.5,
      fill = unname(up_color), color = "white",
      fontface = "bold", size = label_size,
      label.padding = unit(label_padding, "pt"), label.r = unit(1.5, "pt"),
      lineheight = 0.85,
      inherit.aes = FALSE
    )
  }

  if (nrow(down_lbl) > 0) {
    layers$down_labels <- geom_label(
      data = down_lbl,
      aes(x = lbl_x + nudge_x, y = lbl_y, label = legend_label),
      hjust = 0.5, vjust = 0.5,
      fill = unname(down_color), color = "white",
      fontface = "bold", size = label_size,
      label.padding = unit(label_padding, "pt"), label.r = unit(1.5, "pt"),
      lineheight = 0.85,
      inherit.aes = FALSE
    )
  }

  attr(layers, "max_label_r") <- attr(lbl_df, "max_label_r")
  layers
}

# min_size excludes small gene sets prone to tissue-irrelevant GO artifacts
# (Reimand et al. 2019 Nat Protocols S3.4); n_each = NULL passes all sig terms.
select_ring_terms <- function(go_df, contrast_name, n_each = NULL,
                              databases = c("Hallmark", "GO Slim"),
                              min_size = 15) {
  sig <- go_df |>
    filter(contrast == contrast_name, database %in% databases,
           padj < 0.05, size >= min_size) |>
    arrange(padj)

  up   <- sig  |> filter(NES > 0)
  down <- sig  |> filter(NES < 0)

  if (!is.null(n_each)) {
    up   <- up    |> slice_head(n = n_each)
    down <- down  |> slice_head(n = n_each)
  }

  bind_rows(up, down)
}

center_ring_angles <- function(ring, n_up) {
  n <- nrow(ring)
  if (n < 2 || n_up < 1) return(ring)
  up_mid <- (ring$start_deg[1] + ring$end_deg[min(n_up, n)]) / 2
  offset <- 90 - up_mid   # center Up block at 90° (right side)
  ring$start_deg <- ring$start_deg + offset
  ring$end_deg   <- ring$end_deg   + offset
  ring$mid_deg   <- ring$mid_deg   + offset
  ring$start_rad <- ring$start_deg * pi / 180
  ring$end_rad   <- ring$end_deg   * pi / 180
  ring$mid_rad   <- ring$mid_deg   * pi / 180
  ring
}

# YvO-style ring builder. All arcs are equal width (computed so n_total arcs +
# gaps cover the circumference). Up cluster centered at 3 o'clock (90°),
# Down cluster centered at 9 o'clock (270°). Whichever direction has more
# arcs naturally "encroaches" past 12 and 6 o'clock into the other side.
# Most-significant arc at the TOP of each cluster (toward 12 o'clock), least
# at the BOTTOM (toward 6 o'clock). Single-direction case (n_up=0 or n_down=0)
# spans the full circle starting at 12 o'clock going clockwise.
# Drop-in replacement for build_ring_with_gaps; called by Mito _build_volcano_panel.R.
build_ring_180_split <- function(top_terms, contrast_name, go_df,
                                 databases  = c("Hallmark", "Reactome"),
                                 gap_intra  = 3,    # gap between arcs within a cluster
                                 gap_split  = 8,    # extra gap at 12 and 6 o'clock between Up/Down
                                 min_height = 0.05,
                                 max_height = 1.6) {
  if (nrow(top_terms) == 0) return(tibble())

  # top_terms inherits all columns from the fgsea cache (incl. leadingEdge).
  # If a caller passes a trimmed top_terms without leadingEdge, fall back to
  # looking it up from go_df.
  ring <- top_terms
  if (!"leadingEdge" %in% names(ring)) {
    le_lookup <- go_df |>
      dplyr::filter(contrast == !!contrast_name,
                    pathway %in% ring$pathway) |>
      dplyr::distinct(pathway, .keep_all = TRUE) |>
      dplyr::select(pathway, leadingEdge)
    ring <- dplyr::left_join(ring, le_lookup, by = "pathway")
  }
  ring <- ring |>
    dplyr::mutate(
      clean_label = clean_ring_label(pathway),
      gene_list   = stringr::str_split(leadingEdge, ";")
    )

  up   <- ring |> dplyr::filter(NES >  0) |> dplyr::arrange(padj)
  down <- ring |> dplyr::filter(NES <= 0) |> dplyr::arrange(padj)
  n_up   <- nrow(up)
  n_down <- nrow(down)

  scale_arc_height <- function(df) {
    if (nrow(df) == 0) return(numeric(0))
    if (nrow(df) == 1) return(5.6 + (min_height + max_height) / 2)
    neg_lp <- -log10(pmax(df$padj, .Machine$double.xmin))
    rng    <- range(neg_lp)
    if (diff(rng) <= 0) return(rep(5.6 + (min_height + max_height) / 2, nrow(df)))
    scaled <- (neg_lp - rng[1]) / diff(rng)
    5.6 + min_height + (max_height - min_height) * sqrt(scaled)
  }
  if (n_up   > 0) up$arc_r1_var   <- scale_arc_height(up)
  if (n_down > 0) down$arc_r1_var <- scale_arc_height(down)

  # Compute arc width: all arcs (Up + Down) share the same width, sized so the
  # complete set + gaps covers the circumference.
  if (n_up > 0 && n_down > 0) {
    # Bilateral case: 2 split gaps (at 12 and 6 o'clock) + intra gaps
    total_gap <- 2 * gap_split + (n_up - 1) * gap_intra + (n_down - 1) * gap_intra
    arc_w     <- (360 - total_gap) / (n_up + n_down)

    # Up cluster centered at 90°
    up_span    <- n_up * arc_w + (n_up - 1) * gap_intra
    up_start   <- 90 - up_span / 2
    up$start_deg <- up_start + (seq_len(n_up) - 1) * (arc_w + gap_intra)
    up$end_deg   <- up$start_deg + arc_w

    # Down cluster centered at 270°. Reverse so most-sig sits near top-left.
    down_span    <- n_down * arc_w + (n_down - 1) * gap_intra
    down_start   <- 270 - down_span / 2
    down_offset  <- (n_down - seq_len(n_down)) * (arc_w + gap_intra)
    down$start_deg <- down_start + down_offset
    down$end_deg   <- down$start_deg + arc_w
  } else {
    # Single-direction case.
    only   <- if (n_up > 0) up else down
    n_only <- nrow(only)
    if (n_only == 1) {
      # Compact 30° arc at the appropriate side (don't sprawl across the circle
      # for a lone term — it would suggest a much stronger signal than exists).
      cluster_center <- if (n_up > 0) 90 else 270
      arc_w <- 30
      only$start_deg <- cluster_center - arc_w / 2
      only$end_deg   <- only$start_deg + arc_w
    } else {
      # 2+ same-direction arcs fill the whole circle, most-sig at 12 o'clock.
      arc_w <- (360 - n_only * gap_intra) / n_only
      only$start_deg <- (seq_len(n_only) - 1) * (arc_w + gap_intra)
      only$end_deg   <- only$start_deg + arc_w
    }
    if (n_up > 0) up <- only else down <- only
  }

  ring <- dplyr::bind_rows(up, down)
  ring$mid_deg   <- (ring$start_deg + ring$end_deg) / 2
  ring$start_rad <- ring$start_deg * pi / 180
  ring$end_rad   <- ring$end_deg   * pi / 180
  ring$mid_rad   <- ring$mid_deg   * pi / 180
  ring$term_idx  <- seq_len(nrow(ring))
  ring
}


build_ring_with_gaps <- function(top_terms, contrast_name, go_df,
                                 n_each = NULL,
                                 databases = c("Hallmark", "GO Slim")) {
  real_rows <- go_df |>
    filter(contrast == contrast_name, pathway %in% top_terms$pathway)
  padj_lookup <- real_rows |>
    dplyr::select(pathway, padj) |>
    distinct(pathway, .keep_all = TRUE)
  go_subset <- real_rows |>
    mutate(padj = match(pathway, top_terms$pathway) * 1e-10)

  ring <- prepare_ring_data(
    go_df = go_subset, contrast = contrast_name,
    n_terms = nrow(top_terms), gap_degrees = 3, start_offset = 0,
    databases = databases
  )

  ring$padj <- padj_lookup$padj[match(ring$pathway, padj_lookup$pathway)]

  n <- nrow(ring)
  n_up <- sum(ring$NES > 0)

  if (n >= 2) {
    gap_normal <- 3; gap_split <- 8
    gaps <- rep(gap_normal, n)
    if (n_up > 0 && n_up < n) gaps[n_up] <- gap_split
    gaps[n] <- gap_split
    arc_budget <- 360 - sum(gaps)
    arc_widths <- rep(arc_budget / n, n)
    min_height <- 0.05; max_height <- 1.6
    neg_lp <- -log10(pmax(ring$padj, .Machine$double.xmin))
    scaled <- (neg_lp - min(neg_lp)) / (max(neg_lp) - min(neg_lp))
    ring$arc_r1_var <- 4.8 + min_height +
      (max_height - min_height) * sqrt(scaled)
    cum_offset <- 0
    for (i in seq_len(n)) {
      if (i > 1) cum_offset <- cum_offset + arc_widths[i - 1] + gaps[i - 1]
      ring$start_deg[i] <- cum_offset
      ring$end_deg[i]   <- ring$start_deg[i] + arc_widths[i]
      ring$mid_deg[i]   <- (ring$start_deg[i] + ring$end_deg[i]) / 2
      ring$start_rad[i] <- ring$start_deg[i] * pi / 180
      ring$end_rad[i]   <- ring$end_deg[i]   * pi / 180
      ring$mid_rad[i]   <- ring$mid_deg[i]   * pi / 180
    }
    ring <- center_ring_angles(ring, n_up)
  }
  ring
}

make_volcano_ring <- function(de_df,
                              go_df,
                              contrast,
                              title              = NULL,
                              contrast_title     = NULL,
                              contrast_subtitle  = NULL,
                              title_size         = 22,   # scaled from F01's 12pt @ 215mm to F03's 380mm canvas
                              subtitle_size      = NULL,
                              n_terms            = 12,
                              gap_degrees        = 3,
                              start_offset       = 0,
                              databases          = c("Hallmark", "GO Slim"),
                              volcano_radius     = 3.5,
                              tick_r0            = 4.4,
                              tick_r1            = 4.8,
                              arc_r0             = 4.8,
                              arc_r1             = 5.6,
                              label_r            = 7.0,
                              label_gap          = NULL,
                              fc_thresh          = log2(1.5),
                              p_thresh           = 0.05,
                              up_color           = DIR_COLORS["Up"],
                              down_color         = DIR_COLORS["Down"],
                              ns_color           = DIR_COLORS["NS"],
                              point_size         = 0.6,
                              point_alpha        = 0.5,
                              label_size         = 3.0,
                              count_label_size    = 2.8,
                              count_label_padding = 3,
                              count_border_width  = 0.4,
                              count_y_mult        = 1.0,
                              count_x_mult        = 0.5,
                              label_padding       = 2,
                              min_angle_gap       = 18,
                              nudge_outward       = 0.8,
                              ring_data_override = NULL,
                              bg_color           = NULL,
                              bg_alpha           = 0.12,
                              mito_pattern       = NULL,
                              show_legend        = TRUE) {

  if (is.null(subtitle_size)) subtitle_size <- title_size * 0.65

  if (!is.null(ring_data_override)) {
    ring_data <- ring_data_override
  } else {
    ring_data <- prepare_ring_data(
      go_df = go_df, contrast = contrast, n_terms = n_terms,
      gap_degrees = gap_degrees, start_offset = start_offset,
      databases = databases
    )
  }

  tick_data <- build_tick_data(
    ring_data = ring_data, de_df = de_df, contrast = contrast,
    tick_r0 = tick_r0, tick_r1 = tick_r1
  )

  volcano_layers <- build_volcano_layers(
    de_df = de_df, contrast = contrast, volcano_radius = volcano_radius,
    fc_thresh = fc_thresh, p_thresh = p_thresh,
    up_color = up_color, down_color = down_color, ns_color = ns_color,
    point_size = point_size, point_alpha = point_alpha,
    count_label_size = count_label_size,
    count_label_padding = count_label_padding,
    count_border_width = count_border_width,
    count_y_mult = count_y_mult,
    count_x_mult = count_x_mult
  )

  ring_layers <- build_ring_layers(
    ring_data = ring_data, tick_data = tick_data,
    tick_r0 = tick_r0, tick_r1 = tick_r1,
    arc_r0 = arc_r0, arc_r1 = arc_r1,
    up_color = up_color, down_color = down_color,
    mito_pattern = mito_pattern
  )

  label_layers <- build_label_layer(
    ring_data = ring_data, label_r = label_r, label_size = label_size,
    label_gap = label_gap, label_padding = label_padding,
    min_angle_gap = min_angle_gap, nudge_outward = nudge_outward,
    up_color = up_color, down_color = down_color
  )

  max_label_r <- attr(label_layers, "max_label_r")
  if (is.null(max_label_r)) max_label_r <- label_r

  # disc fill behind ring (contrast color, up to tick ring)
  bg_layer <- if (!is.null(bg_color)) {
    bg_circle <- data.frame(
      x = tick_r0 * cos(seq(0, 2 * pi, length.out = 200)),
      y = tick_r0 * sin(seq(0, 2 * pi, length.out = 200))
    )
    geom_polygon(data = bg_circle, aes(x = x, y = y),
                 fill = bg_color, alpha = bg_alpha, color = NA,
                 inherit.aes = FALSE)
  }

  legend_pos <- if (show_legend) "right" else "none"

  title_lab    <- contrast_title %||% title
  subtitle_lab <- contrast_subtitle

  p <- ggplot() +
    bg_layer +
    ring_layers$tick_bg +
    volcano_layers +
    ring_layers$ticks +
    ring_layers$enrich_arcs +
    ring_layers$fill_scale +
    ring_layers$mito_outline +
    label_layers +
    labs(title = title_lab, subtitle = subtitle_lab) +
    coord_fixed(
      xlim = c(-(max_label_r + 0.8), max_label_r + 0.8),
      # Top/bottom: extra headroom so multi-line labels at 12 and 6 o'clock
      # don't clip into the composite title strip / panel border.
      ylim = c(-(max_label_r + 1.2), max_label_r + 1.2),
      clip = "off"
    ) +
    theme_void() +
    theme(plot.title    = element_text(face = "bold", size = title_size,
                                       hjust = 0.5, margin = margin(b = 0, unit = "mm")),
          plot.subtitle = element_text(face = "bold.italic", size = subtitle_size,
                                       color = "grey30", hjust = 0.5,
                                       margin = margin(b = 0.5, unit = "mm")),
          plot.tag      = element_text(face = "bold", size = 26),  # scaled from F01's 15pt @ 215mm to F03's 380mm canvas
          plot.tag.position = c(0.02, 0.99),
          plot.margin   = margin(1, 1, 1, 1, "mm"),
          legend.position = legend_pos,
          legend.title = element_text(size = 7, face = "bold", color = "grey30"),
          legend.text  = element_text(size = 6, color = "grey40"),
          legend.key.width  = unit(2, "mm"),
          legend.key.height = unit(12, "mm"),
          legend.margin = margin(l = 0, r = 0)) +
    guides(color = "none",
           fill = guide_colorbar(direction = "vertical"))

  attr(p, "ring_data")  <- ring_data
  attr(p, "tick_data")  <- tick_data
  attr(p, "n_up")       <- attr(volcano_layers, "n_up")
  attr(p, "n_down")     <- attr(volcano_layers, "n_down")

  p
}

build_nes_legend_bar <- function(text_size = 10, title_size = 11,
                                 bar_margin = margin(-9, 120, 0, 120, "mm")) {
  nes_data <- data.frame(NES = seq(-3, 3, length.out = 200), y = 1)
  ggplot(nes_data, aes(x = .data$NES, y = .data$y, fill = .data$NES)) +
    geom_raster(interpolate = TRUE) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(
      colours = c("#08306B", "#4393C3", "white", "#D6604D", "#67000D"),
      values  = scales::rescale(c(-3, -1.5, 0, 1.5, 3)),
      limits  = c(-3, 3),
      guide   = "none"
    ) +
    scale_x_continuous(
      breaks = c(-3, -1.5, 0, 1.5, 3),
      labels = c("-3", "-1.5", "0", "1.5", "3"),
      expand = c(0, 0)
    ) +
    labs(x = "NES") +
    theme_void() +
    theme(
      axis.text.x  = element_text(size = text_size, face = "bold", color = "grey25"),
      axis.title.x = element_text(size = title_size, face = "bold",
                                   color = "grey25", margin = margin(t = 1)),
      plot.margin  = bar_margin
    )
}
