# EnrichmentMap combined similarity coefficient (Merico 2010 PMID 21085593;
# Reimand 2019 PMID 30664679 Box 1):
#   jaccard    = |A n B| / |A u B|
#   overlap    = |A n B| / min(|A|, |B|)         (Szymkiewicz-Simpson)
#   similarity = 0.5 * overlap + 0.5 * jaccard
# Two sets are redundant when similarity >= cutoff (default 0.375, the
# EnrichmentMap default). Greedy in padj order, so the most significant member
# of a redundant cluster is the one kept. Overlaps are computed on the full
# gene sets, so redundancy reflects database structure (contrast-independent).

# Fisher 2x2 odds-ratio per fora row given the foreground size K and universe N.
# a = overlap (hits in pathway), b = K - a, c = size - a, d = N - K - c.
fora_odds_ratio <- function(overlap, size, K, N) {
  vapply(seq_along(overlap), function(i) {
    a <- overlap[i]; b <- K - a; c <- size[i] - a; d <- N - K - c
    if (b <= 0 || c <= 0) Inf else (a * d) / (b * c)
  }, numeric(1))
}

deduplicate_enrichment_flat <- function(results, pathways, cutoff = 0.375) {
  if (nrow(results) == 0) return(results)
  results <- results[order(results$padj), ]
  kept_sets <- list()
  keep_mask <- logical(nrow(results))
  for (i in seq_len(nrow(results))) {
    pw_genes <- pathways[[results$pathway[i]]]
    if (is.null(pw_genes)) { keep_mask[i] <- TRUE; next }
    redundant <- FALSE
    for (kept in kept_sets) {
      inter <- length(intersect(pw_genes, kept))
      if (inter == 0L) next
      jaccard <- inter / (length(pw_genes) + length(kept) - inter)
      overlap <- inter / min(length(pw_genes), length(kept))
      if (0.5 * overlap + 0.5 * jaccard >= cutoff) { redundant <- TRUE; break }
    }
    if (!redundant) {
      keep_mask[i] <- TRUE
      kept_sets[[length(kept_sets) + 1]] <- pw_genes
    }
  }
  results[keep_mask, ]
}

# Within-database collapse first, then a cross-database pass so the same biology
# in different DBs (REACTOME_TCA_CYCLE vs KEGG_CITRATE_CYCLE) reduces to one entry.
deduplicate_enrichment <- function(results, pathways, cutoff = 0.375, cross_db = TRUE) {
  if (nrow(results) == 0) return(results)
  if (!"database" %in% names(results))
    return(deduplicate_enrichment_flat(results, pathways, cutoff))
  within <- lapply(unique(results$database), function(db)
    deduplicate_enrichment_flat(results[results$database == db, ], pathways, cutoff))
  survivors <- do.call(rbind, within)
  survivors <- survivors[order(survivors$padj), ]
  if (cross_db && nrow(survivors) > 1)
    survivors <- deduplicate_enrichment_flat(survivors, pathways, cutoff)
  survivors
}
