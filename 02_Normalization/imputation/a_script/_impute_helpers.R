# Shared I/O for the three imputation arms. Sourced, never run directly.
# Keeps the method-specific call in each script; the read/write boilerplate lives here once.

load_normalized_matrix <- function() {
  dal <- readRDS(here::here("02_Normalization", "c_data", "DAList_normalized.rds"))
  list(dal = dal, mat = as.matrix(dal$data))
}

save_imputed_dalist <- function(dal, imp, provenance, tag) {
  dimnames(imp) <- dimnames(dal$data)
  stopifnot(sum(is.na(imp)) == 0, identical(dim(imp), dim(dal$data)))
  dal$data <- imp
  dal$imputation <- provenance
  out_dir <- here::here("02_Normalization", "imputation", "c_data")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(dal, file.path(out_dir, sprintf("DAList_imputed_%s.rds", tag)))
}
