# 04_Figures_v2/functions/03_pathway_enrichment_dedup_ora.R
# Pathway enrichment engine for the figure suite: Jaccard/overlap dedup with a
# Fisher 2x2 odds ratio, the gene-set collection + fgsea/fora runners, and the
# display-name and functional-category dictionaries. A snapshot of the validated
# 04_Figures engine; read by the pathway-bar, ring, cluster, and composite figures.

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

# Pathway display labels + functional category dictionaries. Depends on
# clean_pathway_name() from 01_style_palettes_theme.R.

# Single shared display-name dictionary (was duplicated as DISPLAY_LABELS in each
# _panel_A_quadrant.R and display_overrides in each 01_main_panels.R). Keys are
# clean_pathway_name() output (Title Case, "_"->" ", DB prefix stripped,
# MitoCarta reduced to its most-specific last segment).
PATHWAY_DISPLAY_OVERRIDES <- c(
  "Oxidative Phosphorylation"                              = "OXPHOS",
  "Oxphos Subunits"                                        = "OXPHOS Subunits",
  "Cv Subunits"                                            = "Complex V Subunits",
  "Ci Subunits"                                            = "Complex I Subunits",
  "Civ Subunits"                                           = "Complex IV Subunits",
  "Aerobic Respiration And Respiratory Electron Transport" = "Aerobic Resp. + ETC",
  "Respiratory Electron Transport"                         = "Respiratory ETC",
  "Complex I Biogenesis"                                   = "Complex I Biogenesis",
  "Citric Acid Cycle Tca Cycle"                            = "TCA Cycle",
  "Tca Cycle"                                              = "TCA Cycle",
  "Mitochondrial Translation"                              = "Mito Translation",
  "Mitochondrial Ribosome"                                 = "Mitoribosome",
  "Mitochondrial Protein Degradation"                      = "Mito Protein Deg.",
  "Mitochondrial Calcium Ion Transport"                    = "Mito Ca²⁺ Transport",
  "Pink1 Prkn Mediated Mitophagy"                          = "PINK1/PRKN Mitophagy",
  "Fatty Acid Beta Oxidation"                              = "FA β-Oxidation",
  "Fatty Acid Oxidation"                                   = "FA Oxidation",
  "Fatty Acid Metabolism"                                  = "FA Metabolism",
  "Branched Chain Amino Acid Metabolism"                   = "BCAA Metabolism",
  "Cardiolipin Biosynthesis"                               = "Cardiolipin Biosynth.",
  "Mitochondrion Organization"                             = "Mito Organization",
  "Mitochondrial Organization"                             = "Mito Org.",
  "Mitochondrial Transport"                                = "Mito Transport",
  "Mitochondrial Protein Import"                           = "Mito Protein Import",
  "Small Molecule Transport"                               = "Sm. Molecule Transport",
  "Slc25a Family"                                          = "SLC25A Carriers",
  "Proteases"                                              = "Mito Proteases",
  "Ros And Glutathione Metabolism"                         = "ROS/Glutathione Metab.",
  "Xenobiotic Metabolism"                                  = "Xenobiotic Metab.",
  "Lysine Metabolism"                                      = "Lysine Metab.",
  "Amino Acid Metabolic Process"                           = "Amino Acid Metab.",
  "Amino Acid Metabolism"                                  = "Amino Acid Metab.",
  "Extracellular Matrix Organization"                      = "ECM Organization",
  "Epithelial Mesenchymal Transition"                      = "EMT",
  "Cytoplasmic Translation"                                = "Cytoplasmic Transl.",
  "Glycolysis Gluconeogenesis"                             = "Glycolysis/Gluconeo.",
  "Generation Of Precursor Metabolites And Energy"         = "Precursor Metab. & Energy",
  "Organophosphate Metabolic Process"                      = "Organophosphate Metab.",
  "Organic Acid Metabolic Process"                         = "Organic Acid Metab.",
  "Carboxylic Acid Metabolic Process"                      = "Carboxylic Acid Metab.",
  "Monocarboxylic Acid Metabolic Process"                  = "Monocarbox. Acid Metab.",
  "Small Molecule Catabolic Process"                       = "Sm. Molecule Catab.",
  "Small Molecule Metabolic Process"                       = "Sm. Molecule Metab.",
  "Ketone Metabolic Process"                               = "Ketone Metabolism",
  "Protein Localization To Plasma Membrane"                = "Plasma Membr. Protein Loc.",
  "Striated Muscle Contraction"                            = "Striated Muscle Contr."
)

# Human-readable pathway label for any DB. MitoCarta names are a
# "MITOCARTA_A__B__TERM" hierarchy -> keep the most specific last segment, drop
# the prefix; everything else goes through clean_pathway_name(). Then apply the
# shared override dictionary (extra = figure-specific extension, merged on top).
clean_display_label <- function(pathway, extra = NULL) {
  is_mito <- grepl("^MITOCARTA_", pathway)
  lab <- ifelse(
    is_mito,
    clean_pathway_name(gsub("_", " ", sub(".*__", "", sub("^MITOCARTA_", "", pathway)))),
    clean_pathway_name(pathway))
  dict <- PATHWAY_DISPLAY_OVERRIDES
  if (!is.null(extra)) dict[names(extra)] <- extra
  hit <- match(lab, names(dict))
  ifelse(!is.na(hit), unname(dict[hit]), lab)
}

# MSigDB pathway ID -> 15 consolidated categories (keyword rules)
CONSOLIDATED_PATHWAY_ORDER <- c(
  "Muscle & Contractile", "Cytoskeleton & Motility", "ECM & Adhesion",
  "Lipid Metabolism", "Carbohydrate & Energy Metabolism",
  "Amino Acid & Cofactor Metabolism",
  "Mitochondria & Energy", "Protein Homeostasis",
  "Transport", "Translation & Ribosome", "Transcription & Chromatin",
  "Immune & Inflammation", "DNA & Cell Cycle", "Circulatory System",
  "Development", "Other"
)

CONSOLIDATED_COLORS <- c(
  "Muscle & Contractile"              = "#E57373",
  "Cytoskeleton & Motility"           = "#FFB74D",
  "ECM & Adhesion"                    = "#FFF176",
  "Lipid Metabolism"                  = "#AED581",
  "Carbohydrate & Energy Metabolism"  = "#81C784",
  "Amino Acid & Cofactor Metabolism"  = "#66BB6A",
  "Mitochondria & Energy"             = "#4DB6AC",
  "Protein Homeostasis"               = "#4FC3F7",
  "Transport"                         = "#7986CB",
  "Translation & Ribosome"            = "#BA68C8",
  "Transcription & Chromatin"         = "#AB47BC",
  "Immune & Inflammation"             = "#A1887F",
  "DNA & Cell Cycle"                  = "#90A4AE",
  "Circulatory System"                = "#CE93D8",
  "Development"                       = "#B0BEC5",
  "Other"                             = "#D0D0D0"
)

classify_pathway_func <- function(ids) {
  rules <- list(
    "Muscle & Contractile"              = "MYOGEN|MYOFIBRIL|SARCOMERE|MUSCLE_|CONTRACTILE|ACTOMYOSIN|MYOSIN|I_BAND",
    "Cytoskeleton & Motility"           = "CYTOSKELET|ACTIN_BIND|STRUCTURAL_MOLECULE|MOTIL|SUPRAMOLECUL",
    "ECM & Adhesion"                    = "EXTRACELLULAR_MATRIX|COLLAGEN|BASEMENT_MEMBRANE|ADHESION|APICAL_JUNCTION|EMT|ENCAPSULATING",
    "Lipid Metabolism"                  = "FATTY_ACID|LIPID|ADIPOGEN|STEROID|SPHINGOLIPID|PHOSPHOLIPID|KETONE",
    "Carbohydrate & Energy Metabolism"  = "GLYCOLY|GLUCONEO|CARBOHYDRATE|PENTOSE|PRECURSOR_METABOL",
    "Amino Acid & Cofactor Metabolism"  = "AMINO_ACID|VITAMIN|COFACTOR|NITROGEN|DETOXIF|DIGEST|XENOBIOT",
    "Mitochondria & Energy"             = "MITOCHOND|OXIDATIVE_PHOSPH|ELECTRON_TRANSFER|RESPIRATORY|OXIDOREDUCT",
    "Protein Homeostasis"               = "PROTEASOM|UBIQUITIN|AUTOPHAGY|MTORC1|PROTEIN_FOLD",
    "Transport"                         = "TRANSPORT(?!.*ELECTRON)|VESICLE|ENDOCYT|SECRETI",
    "Translation & Ribosome"            = "TRANSLAT|RIBOSOM|TRNA|MYC_TARGET",
    "Transcription & Chromatin"         = "TRANSCRIPT|SPLICEOSOM|E2F_TARGET|CHROMATIN|MRNA_PROC",
    "Immune & Inflammation"             = "IMMUN|INFLAMMA|INTERFERON|IL2|IL6|TNFA|NF.KB|COMPLEMENT",
    "DNA & Cell Cycle"                  = "DNA_REPAIR|CELL_CYCLE|MITOTIC|P53_PATHWAY",
    "Circulatory System"                = "ANGIOGEN|BLOOD_VESSEL|HYPOXIA",
    "Development"                       = "UV_RESPONSE|GROWTH_FACTOR|WNT|HEDGEHOG|NOTCH|TGF_BETA|KRAS"
  )
  vapply(toupper(ids), function(id) {
    matches <- character(0)
    for (cat in names(rules)) {
      if (grepl(rules[[cat]], id, perl = TRUE)) matches <- c(matches, cat)
    }
    if (length(matches) > 1) {
      warning("classify_pathway_func: '", id, "' matches multiple categories [",
              paste(matches, collapse = ", "), "]; using first match: ", matches[1])
    }
    if (length(matches) >= 1) return(matches[1])
    if (grepl("METABOL", id, perl = TRUE)) return("Amino Acid & Cofactor Metabolism")
    "Other"
  }, character(1), USE.NAMES = FALSE)
}
