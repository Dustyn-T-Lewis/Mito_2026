# Pathway enrichment engine: EnrichmentMap dedup, the harmonized gene-set
# collection, database classification, and the display-label dictionary.
# Read by the pathway-bar, ring, and cluster figures.

# EnrichmentMap combined similarity coefficient (Merico 2010 PMID 21085593;
# Reimand 2019 PMID 30664679 Box 1):
#   jaccard    = |A n B| / |A u B|
#   overlap    = |A n B| / min(|A|, |B|)         (Szymkiewicz-Simpson)
#   similarity = 0.5 * overlap + 0.5 * jaccard
# Two sets are redundant when similarity >= cutoff (default 0.375, the
# EnrichmentMap default). Greedy in padj order, so the most significant member
# of a redundant cluster is the one kept. Overlaps are computed on the full
# gene sets, so redundancy reflects database structure (contrast-independent).
deduplicate_enrichment_flat <- function(results, pathways, cutoff = 0.375) {
  if (nrow(results) == 0) {
    return(results)
  }
  results <- results[order(results$padj), ]
  kept_sets <- list()
  keep_mask <- logical(nrow(results))
  for (i in seq_len(nrow(results))) {
    pw_genes <- pathways[[results$pathway[i]]]
    if (is.null(pw_genes)) {
      keep_mask[i] <- TRUE
      next
    }
    redundant <- FALSE
    for (kept in kept_sets) {
      inter <- length(intersect(pw_genes, kept))
      if (inter == 0L) next
      jaccard <- inter / (length(pw_genes) + length(kept) - inter)
      overlap <- inter / min(length(pw_genes), length(kept))
      if (0.5 * overlap + 0.5 * jaccard >= cutoff) {
        redundant <- TRUE
        break
      }
    }
    if (!redundant) {
      keep_mask[i] <- TRUE
      kept_sets[[length(kept_sets) + 1]] <- pw_genes
    }
  }
  results[keep_mask, ]
}

# Audit trail for the supplementary table: flags each pathway as the retained
# representative or redundant, and for redundant ones names the retained pathway
# it overlaps most and the overlap+jaccard coefficient to it. Reuses
# deduplicate_enrichment for the retained set, so the flags cannot drift from the
# ring. overlap+jaccard is the EnrichmentMap combined coefficient (see above).
dedup_report <- function(results, pathways, cutoff = 0.375, cross_db = TRUE) {
  results$dedup_status <- "kept"
  results$merged_into <- NA_character_
  results$overlap_jaccard <- NA_real_
  if (nrow(results) == 0) {
    return(results)
  }
  kept <- deduplicate_enrichment(results, pathways, cutoff, cross_db)$pathway
  redundant <- setdiff(results$pathway, kept)
  results$dedup_status[results$pathway %in% redundant] <- "redundant"
  coef <- function(a, b) {
    genes_a <- pathways[[a]]
    genes_b <- pathways[[b]]
    inter <- length(intersect(genes_a, genes_b))
    if (inter == 0L) {
      return(0)
    }
    0.5 * inter / min(length(genes_a), length(genes_b)) +
      0.5 * inter / (length(genes_a) + length(genes_b) - inter)
  }
  for (p in redundant) {
    if (is.null(pathways[[p]])) next
    sims <- vapply(kept, coef, numeric(1), a = p)
    best <- which.max(sims)
    results$merged_into[results$pathway == p] <- kept[best]
    results$overlap_jaccard[results$pathway == p] <- sims[best]
  }
  results[order(results$padj), , drop = FALSE]
}

# Within-database collapse first, then a cross-database pass so the same biology
# in different DBs (REACTOME_TCA_CYCLE vs KEGG_CITRATE_CYCLE) reduces to one entry.
deduplicate_enrichment <- function(results, pathways, cutoff = 0.375, cross_db = TRUE) {
  if (nrow(results) == 0) {
    return(results)
  }
  if (!"database" %in% names(results)) {
    return(deduplicate_enrichment_flat(results, pathways, cutoff))
  }
  within <- lapply(unique(results$database), function(db) {
    deduplicate_enrichment_flat(results[results$database == db, ], pathways, cutoff)
  })
  survivors <- do.call(rbind, within)
  survivors <- survivors[order(survivors$padj), ]
  if (cross_db && nrow(survivors) > 1) {
    survivors <- deduplicate_enrichment_flat(survivors, pathways, cutoff)
  }
  survivors
}

# Canonical 5-DB pathway lens — single source of truth for every enrichment
# figure (rings, pathway bars, cluster ORA). MITO_DROP_SETS are MitoCarta
# localization/aggregate sets (not pathways) excluded everywhere.
CANONICAL_DBS <- c("Hallmark", "Reactome", "KEGG", "MitoCarta", "GO Slim")
MITO_DROP_SETS <- c(
  "MITOCARTA_ALL", "MITOCARTA_IMM", "MITOCARTA_IMS",
  "MITOCARTA_MATRIX", "MITOCARTA_OMM"
)

# Pathogen / disease sets are irrelevant to an in-vitro mito-transplant study and
# leak in as high-NES noise (Reactome "Dengue Virus Host Interactions" overlaps
# the translation machinery), so they are dropped from ORA collections.
DISEASE_VIRAL_RE <- paste0(
  "DISEASE|CANCER|TUMOR|CARCINOMA|LEUKEMIA|LYMPHOMA|MELANOMA|GLIOMA|",
  "HEPATITIS|HIV|INFECTION|INFECTIOUS|VIRAL|VIRUS|INFLUENZA|SARS|HCMV|",
  "MEASLES|BACTERIAL|LISTERIA|LEISHMANIA|PARASIT"
)

# Single gate for "is this pathway an eligible, significant hit for this
# contrast" — every enrichment figure (F01 bars, F02 rings, F03 module ORA)
# calls this instead of re-deriving the filter, so they can't drift apart the
# way the F01 bars and F02 rings once did over the disease/viral exclusion.
significant_pathways <- function(fgsea_df, ctr, padj_cut = 0.05, min_size = 10, max_size = 500) {
  fgsea_df |>
    dplyr::filter(
      contrast == ctr, database %in% CANONICAL_DBS,
      !is.na(padj), padj < padj_cut,
      size >= min_size, size <= max_size,
      !pathway %in% MITO_DROP_SETS,
      !grepl(DISEASE_VIRAL_RE, pathway, ignore.case = TRUE)
    ) |>
    dplyr::arrange(padj)
}

classify_database <- function(pathway_names) {
  dplyr::case_when(
    grepl("^HALLMARK_", pathway_names) ~ "Hallmark",
    grepl("^REACTOME_", pathway_names) ~ "Reactome",
    grepl("^KEGG_MEDICUS_", pathway_names) ~ "KEGG",
    grepl("^KEGG_", pathway_names) ~ "KEGG",
    grepl("^GOSLIM_", pathway_names) ~ "GO Slim",
    grepl("^GOBP_", pathway_names) ~ "GO:BP",
    grepl("^GOCC_", pathway_names) ~ "GO:CC",
    grepl("^GOMF_", pathway_names) ~ "GO:MF",
    grepl("^MITOCARTA_", pathway_names) ~ "MitoCarta",
    TRUE ~ "Other"
  )
}

# Per-database over-representation: fora within each DB so BH is applied within a
# database, not across it (valid even for a single-pathway DB, where BH just
# returns the raw p-value). Returns the raw fora columns (pval, padj, overlap,
# size, overlapGenes) plus the database; callers add their own leading-edge and
# column shaping. A DB with no eligible sets is skipped.
run_fora_by_db <- function(genes, universe, pw_collection, pw_by_db,
                           min_size = 10, max_size = 500) {
  g <- intersect(genes, universe)
  dplyr::bind_rows(lapply(names(pw_by_db), function(db) {
    dp <- pw_collection[pw_by_db[[db]]]
    if (length(dp) < 1) {
      return(NULL)
    }
    r <- as.data.frame(fgsea::fora(
      pathways = dp, genes = g, universe = universe,
      minSize = min_size, maxSize = max_size
    ))
    r$database <- db
    r
  }))
}

# Flat gene-set collection over CANONICAL_DBS for ORA. Reads the unified
# rat_gene_sets.rds so it enriches against the same backbone as the rings.
build_harmonized_collection <- function(
  cache = here::here("04_Figures", "shared", "c_data", "rat_gene_sets.rds"),
  min_size = 10, max_size = 500
) {
  gs <- readRDS(cache)
  pw <- do.call(c, unname(gs[CANONICAL_DBS]))
  pw <- pw[!names(pw) %in% MITO_DROP_SETS]
  pw <- pw[!grepl(DISEASE_VIRAL_RE, names(pw), ignore.case = TRUE)]
  pw <- pw[!duplicated(names(pw))]
  pw <- lapply(pw, unique)
  sizes <- vapply(pw, length, integer(1))
  pw <- pw[sizes >= min_size & sizes <= max_size]
  message(sprintf(
    "Harmonized collection: %d sets (%s), size %d-%d",
    length(pw), paste(CANONICAL_DBS, collapse = "+"), min_size, max_size
  ))
  pw
}

# Display label for any pathway via the enrichVolcano cleaner, the same logic the volcano
# rings use. Its wrapping newlines are flattened so the pathway-bar labels stay one row.
clean_display_label <- function(pathway) {
  if (!length(pathway)) {
    return(character(0))
  }
  gsub("\n", " ", enrichVolcano::ev_clean_label(pathway))
}
