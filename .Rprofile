# Project-level R startup. Routes the default graphics device to a tempfile in
# non-interactive runs so Rscript invocations stop leaving Rplots.pdf in cwd.
if (!interactive()) {
  options(device = function(...) {
    grDevices::pdf(tempfile(fileext = ".pdf"), ...)
  })
}
