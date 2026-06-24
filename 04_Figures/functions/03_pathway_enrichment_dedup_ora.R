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

# Display-label dictionary for the enrichment figures, keyed on
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
    clean_pathway_name(pathway)
  )
  dict <- PATHWAY_DISPLAY_OVERRIDES
  if (!is.null(extra)) dict[names(extra)] <- extra
  hit <- match(lab, names(dict))
  ifelse(!is.na(hit), unname(dict[hit]), lab)
}
