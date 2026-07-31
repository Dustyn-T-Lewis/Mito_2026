# Mito 2026 — install all R packages used across the pipeline.
# Run once after cloning: Rscript setup.R

cran <- c(
  "here", "pacman", "jsonlite",
  "readxl", "readr", "openxlsx",
  "dplyr", "tidyr", "tibble", "stringr", "purrr", "forcats",
  "ggplot2", "ggrepel", "patchwork", "ggtext", "ggnewscale", "ggfittext",
  "scales", "shadowtext", "ggplotify", "eulerr", "vegan",
  "missForest", "imp4p", "imputeLCMD", "msigdbr", "WGCNA",
  "ragg", "testthat"
)

bioc <- c(
  "limma", "fgsea", "AnnotationDbi", "GO.db", "MsCoreUtils", "impute",
  "org.Rn.eg.db"
)

github <- c(
  proteoDA      = "ByrumLab/proteoDA",
  enrichVolcano = "Dustyn-T-Lewis/enrichVolcano"
)

installed <- rownames(installed.packages())

need_cran <- setdiff(cran, installed)
if (length(need_cran)) install.packages(need_cran)

need_bioc <- setdiff(bioc, installed)
if (length(need_bioc)) {
  if (!"BiocManager" %in% installed) install.packages("BiocManager")
  BiocManager::install(need_bioc, ask = FALSE, update = FALSE)
}

for (pkg in names(github)) {
  if (!pkg %in% installed) {
    if (!"remotes" %in% installed) install.packages("remotes")
    remotes::install_github(github[[pkg]])
  }
}

# version record for anyone reproducing the pipeline
declared <- sort(c(cran, bioc, names(github)))
present <- rownames(installed.packages())
versions <- vapply(
  declared,
  function(pkg) {
    if (!pkg %in% present) {
      return("not installed")
    }
    version <- as.character(utils::packageVersion(pkg))
    sha <- utils::packageDescription(pkg)$RemoteSha
    if (is.null(sha)) version else paste0(version, " (", substr(sha, 1, 7), ")")
  },
  character(1)
)

writeLines(
  c(
    R.version.string,
    paste("Platform:", R.version$platform),
    paste("Recorded:", format(Sys.Date())),
    "",
    sprintf("%-16s %s", declared, versions)
  ),
  here::here("package_versions.txt")
)
