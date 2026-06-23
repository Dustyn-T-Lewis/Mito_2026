# Gene-set collection + fgsea/fora runners. Depends on 03a_dedup_engine.R for
# fora_odds_ratio() and deduplicate_enrichment().

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

# Flat gene-set collection over CANONICAL_DBS for ORA. Reads the unified
# rat_gene_sets.rds so it enriches against the same backbone as the rings.
build_harmonized_collection <- function(
  cache = here::here("04_Figures", "shared", "rat_gene_sets.rds"),
  min_size = 10, max_size = 350
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

# One fGSEA pass for a single database: rank each contrast by the limma moderated t,
# run fgseaMultilevel, return cache-schema rows.
run_fgsea_cache <- function(dep_wide, gene_sets, db_name, contrasts,
                            min_size = 10, max_size = 500) {
  one <- function(contrast) {
    stats <- dep_wide[[paste0("t_", contrast)]]
    names(stats) <- dep_wide$gene
    stats <- stats[!is.na(stats) & !is.na(names(stats)) & names(stats) != ""]
    if (anyDuplicated(names(stats))) stats <- tapply(stats, names(stats), mean)
    stats <- sort(stats)
    set.seed(42) # per-call so fgsea's multilevel p-values are independent of loop order
    res <- fgsea::fgseaMultilevel(
      pathways = gene_sets, stats = stats,
      minSize = min_size, maxSize = max_size, eps = 0
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

# GO Slim is OFF by default in build_pathway_collection() for rat workflows —
# collapsePathways + Jaccard dedup already handle GO:BP redundancy.
build_goslim_gene_sets <- function(species = "Rattus norvegicus",
                                   orgdb = NULL,
                                   min_size = 10, max_size = 500) {
  # orgdb auto-resolves from species when NULL.
  if (is.null(orgdb)) {
    orgdb <- switch(species,
      "Rattus norvegicus" = "org.Rn.eg.db",
      "Homo sapiens" = "org.Hs.eg.db",
      "Mus musculus" = "org.Mm.eg.db",
      stop(
        "No default orgdb for species '", species,
        "'. Pass orgdb explicitly."
      )
    )
  }
  requireNamespace("GO.db", quietly = TRUE)
  requireNamespace(orgdb, quietly = TRUE)
  requireNamespace("AnnotationDbi", quietly = TRUE)
  org_pkg <- getExportedValue(orgdb, orgdb)

  # 62 GO Slim Generic BP terms
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
    tryCatch(AnnotationDbi::Term(GO.db::GOTERM[[id]]),
      error = function(e) NA_character_
    )
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

run_ora_deduplicated <- function(genes, universe, pathways,
                                 cutoff = 0.375,
                                 min_size = 10, max_size = 500,
                                 padj_cutoff = 0.05) {
  requireNamespace("fgsea", quietly = TRUE) # fora is exact; no seed needed
  genes <- intersect(genes, universe)

  pw_by_db <- split(names(pathways), classify_database(names(pathways)))
  db_results <- list()
  for (db in names(pw_by_db)) {
    db_pw <- pathways[pw_by_db[[db]]]
    if (length(db_pw) < 2) next
    db_res <- fgsea::fora(
      pathways = db_pw,
      genes    = genes,
      universe = universe,
      minSize  = min_size,
      maxSize  = max_size
    )
    db_res <- as.data.frame(db_res)
    db_res$database <- db
    db_results[[db]] <- db_res
  }
  res <- do.call(rbind, db_results)

  res$odds_ratio <- fora_odds_ratio(res$overlap, res$size,
    K = length(genes), N = length(universe)
  )
  res <- tibble::as_tibble(res)

  sig <- res[!is.na(res$padj) & res$padj < padj_cutoff, ]
  sig_dedup <- deduplicate_enrichment(sig, pathways, cutoff = cutoff)

  n_removed <- nrow(sig) - nrow(sig_dedup)
  pct <- if (nrow(sig) > 0) round(100 * n_removed / nrow(sig), 1) else 0
  message(sprintf(
    "ORA dedup: %d sig -> %d kept (removed %d, %.1f%%)",
    nrow(sig), nrow(sig_dedup), n_removed, pct
  ))

  sig_dedup
}
