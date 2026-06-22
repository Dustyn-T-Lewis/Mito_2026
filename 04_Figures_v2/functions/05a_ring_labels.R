# Library + style imports for the volcano-ring suite are loaded once by the
# stub 05_volcano_ring_plot_builder.R; do not re-source them here.

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

