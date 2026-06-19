# Workbook assembly for the 90_stitch scripts. Each figure's c_data ends as a
# single supplementary .xlsx: build_workbook bundles the named sheet specs, sweeps
# in any remaining result CSVs from the same c_data, then deletes the loose CSVs
# (keeping .rds caches) so the workbook is the one tabular deliverable.

suppressPackageStartupMessages({
  library(openxlsx)
  library(readr)
})

add_sheet <- function(wb, name, data) {
  addWorksheet(wb, name)
  writeData(wb, name, data)
  addStyle(wb, name, createStyle(textDecoration = "bold"),
           rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
  cat(sprintf("    + %s: %d x %d\n", name, nrow(data), ncol(data)))
}

# Excel sheet names: <=31 chars, no []:*?/\, unique.
sheet_name <- function(base, taken) {
  nm <- substr(gsub("[\\[\\]:*?/\\\\]", "_", base), 1, 31)
  i <- 2L
  while (nm %in% taken) { nm <- paste0(substr(nm, 1, 29), "_", i); i <- i + 1L }
  nm
}

build_workbook <- function(out_file, sheet_specs) {
  wb <- createWorkbook()
  taken <- character(0)
  bundled <- character(0)
  for (spec in sheet_specs) {
    df <- if (!is.null(spec$df)) spec$df else if (file.exists(spec$path)) as.data.frame(read_csv(spec$path, show_col_types = FALSE)) else NULL
    if (is.null(df)) next
    add_sheet(wb, spec$name, df)
    taken <- c(taken, spec$name)
    if (!is.null(spec$path)) bundled <- c(bundled, normalizePath(spec$path))
  }
  # Sweep remaining result CSVs in this c_data into the workbook (filename-based
  # sheet names), skipping any already bundled by path or name.
  root  <- normalizePath(dirname(out_file))
  loose <- list.files(root, pattern = "[.]csv$", recursive = TRUE, full.names = TRUE)
  for (p in loose) {
    base <- sub("[.]csv$", "", basename(p))
    if (normalizePath(p) %in% bundled || base %in% taken) next
    add_sheet(wb, nm <- sheet_name(base, taken),
              as.data.frame(read_csv(p, show_col_types = FALSE)))
    taken <- c(taken, nm)
  }
  saveWorkbook(wb, out_file, overwrite = TRUE)
  for (p in loose) if (file.exists(p)) unlink(p)
  for (d in list.dirs(root, recursive = FALSE)) if (!length(list.files(d, all.files = TRUE, no.. = TRUE))) unlink(d, recursive = TRUE)
  cat(sprintf("  Saved: %s (%.0f KB)\n\n", out_file, file.size(out_file) / 1e3))
}
