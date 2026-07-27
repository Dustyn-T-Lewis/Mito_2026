# Builder helpers for the gene-set and fGSEA-cache scripts under shared/a_script.
#   run_fgsea_cache()        fGSEA (multilevel, moderated-t ranks) for one database
#   build_goslim_gene_sets() GO Slim Generic BP sets via GO.db + a species orgdb
# Dedup lives in shared_enrichment_ora.R; these two are only the
# generation helpers that engine does not carry.

# eps is fgsea's smallest resolvable p-value. The default 0 makes fgseaMultilevel refine
# until it converges, which is right for the pathway cache but pathological on WGCNA
# modules: a 1000-protein set returns p ~1e-200 and the refinement runs for minutes.
# Callers that only need NES should pass a floor.
run_fgsea_cache <- function(dep_wide, gene_sets, db_name, contrasts,
                            min_size = 10, max_size = 500, eps = 0) {
  one <- function(contrast) {
    stats <- dep_wide[[paste0("t_", contrast)]]
    names(stats) <- dep_wide$gene
    stats <- stats[!is.na(stats) & !is.na(names(stats)) & names(stats) != ""]
    if (anyDuplicated(names(stats))) stats <- tapply(stats, names(stats), mean)
    stats <- sort(stats)
    set.seed(42) # per-call so fgsea's multilevel p-values are independent of loop order
    res <- fgsea::fgseaMultilevel(
      pathways = gene_sets, stats = stats,
      minSize = min_size, maxSize = max_size, eps = eps
    )
    if (nrow(res) == 0) {
      return(NULL)
    }
    res |>
      dplyr::mutate(
        database = db_name, contrast = contrast,
        leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";")
      ) |>
      dplyr::select(pathway, pval, padj, log2err, ES, NES, size, leadingEdge, database, contrast)
  }
  dplyr::bind_rows(lapply(contrasts, one))
}

build_goslim_gene_sets <- function(species = "Rattus norvegicus",
                                   orgdb = NULL,
                                   min_size = 10, max_size = 500) {
  # orgdb auto-resolves from species when NULL.
  if (is.null(orgdb)) {
    orgdb <- switch(species,
      "Rattus norvegicus" = "org.Rn.eg.db",
      "Homo sapiens"      = "org.Hs.eg.db",
      "Mus musculus"      = "org.Mm.eg.db",
      stop("No default orgdb for species '", species, "'. Pass orgdb explicitly.")
    )
  }
  requireNamespace("GO.db", quietly = TRUE)
  requireNamespace(orgdb, quietly = TRUE)
  requireNamespace("AnnotationDbi", quietly = TRUE)
  org_pkg <- getExportedValue(orgdb, orgdb)

  # GO Slim Generic BP terms: the Biological Process subset of the GO Consortium
  # Generic GO-Slim (goslim_generic.obo). 62 IDs, expanded over GOBPOFFSPRING below.
  bp_slim <- c(
    "GO:0000278", "GO:0000910", "GO:0002181", "GO:0002376", "GO:0003012",
    "GO:0003013", "GO:0003014", "GO:0003016", "GO:0005975", "GO:0006091",
    "GO:0006260", "GO:0006281", "GO:0006310", "GO:0006325", "GO:0006351",
    "GO:0006355", "GO:0006399", "GO:0006457", "GO:0006520", "GO:0006629",
    "GO:0006766", "GO:0006886", "GO:0006913", "GO:0006914", "GO:0006954",
    "GO:0007005", "GO:0007010", "GO:0007018", "GO:0007031", "GO:0007059",
    "GO:0007126", "GO:0007155", "GO:0007163", "GO:0007586", "GO:0009100",
    "GO:0012501", "GO:0016071", "GO:0016192", "GO:0023052", "GO:0030154",
    "GO:0030163", "GO:0030198", "GO:0032200", "GO:0034330", "GO:0042060",
    "GO:0042180", "GO:0042254", "GO:0044782", "GO:0048856", "GO:0048870",
    "GO:0050877", "GO:0051604", "GO:0055085", "GO:0055086", "GO:0061024",
    "GO:0065003", "GO:0071941", "GO:0072659", "GO:0098542", "GO:0098754",
    "GO:0140014", "GO:1901135"
  )

  offspring <- as.list(GO.db::GOBPOFFSPRING)

  suppressMessages({
    go_genes <- AnnotationDbi::select(
      org_pkg,
      keys = AnnotationDbi::keys(org_pkg, keytype = "GO"),
      keytype = "GO",
      columns = c("SYMBOL", "ONTOLOGY")
    )
  })
  go_bp_genes <- go_genes[!is.na(go_genes$ONTOLOGY) & go_genes$ONTOLOGY == "BP", ]
  go_to_symbols <- split(go_bp_genes$SYMBOL, go_bp_genes$GO)

  goslim_sets <- list()
  slim_names <- vapply(bp_slim, function(id) {
    tryCatch(AnnotationDbi::Term(GO.db::GOTERM[[id]]), error = function(e) NA_character_)
  }, character(1))

  for (i in seq_along(bp_slim)) {
    go_id <- bp_slim[i]
    go_term <- slim_names[i]
    if (is.na(go_term)) next

    all_terms <- go_id
    desc <- offspring[[go_id]]
    if (!is.null(desc)) all_terms <- c(all_terms, desc)

    genes <- unique(unlist(go_to_symbols[intersect(all_terms, names(go_to_symbols))],
      use.names = FALSE
    ))
    genes <- genes[!is.na(genes)]

    if (length(genes) >= min_size && length(genes) <= max_size) {
      set_name <- paste0("GOSLIM_", toupper(gsub(" ", "_", go_term)))
      goslim_sets[[set_name]] <- genes
    }
  }

  message(sprintf(
    "GO Slim: %d/%d terms passed size filter (%d-%d)",
    length(goslim_sets), length(bp_slim), min_size, max_size
  ))
  goslim_sets
}
