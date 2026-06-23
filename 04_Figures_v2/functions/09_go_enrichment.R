# GO over-representation for WGCNA modules — clusterProfiler enrichGO + simplify.
# Rat (org.Rn.eg.db), gene-symbol keyed, against the detected-proteome universe.
# Semantic-similarity simplify() collapses redundant GO terms (the canonical
# clusterProfiler path, vs the suite's Jaccard dedup used for MSigDB sets).

# enrichGO + simplify for one gene set in one ontology. Returns a tidy tibble
# (term, id, padj, n_overlap, gene_ratio) ordered by padj, or an empty tibble.
run_go_simplified <- function(genes, universe, ontology = c("BP", "CC", "MF"),
                              padj_cutoff = 0.05, qval_cutoff = 0.10,
                              min_size = 10, max_size = 500,
                              simplify_cutoff = 0.7) {
  ontology <- match.arg(ontology)
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  universe <- unique(universe[!is.na(universe) & nzchar(universe)])
  if (length(genes) < min_size) {
    return(.empty_go())
  }

  ego <- clusterProfiler::enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Rn.eg.db::org.Rn.eg.db,
    keyType       = "SYMBOL",
    ont           = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff  = padj_cutoff,
    qvalueCutoff  = qval_cutoff,
    minGSSize     = min_size,
    maxGSSize     = max_size,
    readable      = FALSE
  )
  if (is.null(ego) || nrow(ego) == 0) {
    return(.empty_go())
  }

  ego <- clusterProfiler::simplify(
    ego,
    cutoff = simplify_cutoff, by = "p.adjust", select_fun = min
  )
  df <- as.data.frame(ego)
  if (nrow(df) == 0) {
    return(.empty_go())
  }

  tibble::tibble(
    ontology   = ontology,
    id         = df$ID,
    term       = df$Description,
    padj       = df$p.adjust,
    n_overlap  = df$Count,
    gene_ratio = vapply(df$GeneRatio, .parse_ratio, numeric(1))
  ) |>
    dplyr::arrange(padj)
}

# All three ontologies for one module, top-n each (post-simplify).
run_go_module <- function(genes, universe, top_n = 3, ...) {
  dplyr::bind_rows(lapply(c("BP", "CC", "MF"), function(ont) {
    run_go_simplified(genes, universe, ontology = ont, ...) |>
      head(top_n)
  }))
}

.empty_go <- function() {
  tibble::tibble(
    ontology = character(), id = character(), term = character(),
    padj = numeric(), n_overlap = integer(), gene_ratio = numeric()
  )
}

# "12/345" -> 12/345
.parse_ratio <- function(x) {
  parts <- as.numeric(strsplit(x, "/", fixed = TRUE)[[1]])
  if (length(parts) == 2 && parts[2] > 0) parts[1] / parts[2] else NA_real_
}
