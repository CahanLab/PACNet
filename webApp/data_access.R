# Reference data access for the PACNet web app.
#
# Native R (runApp / runGitHub): files are read from ../data relative to the
# app directory (option pacnet.data_dir overrides).
# webR (Shinylive, in the browser): files are downloaded from the data/
# directory next to the exported site, once per session, then read.

is_webr <- function() identical(R.version$os, "emscripten")

pacnet_data_url <- function() {
  url <- getOption("pacnet.data_url")
  if (!is.null(url)) return(url)
  # The webR worker script lives at <site>/shinylive/webr/webr-worker.js
  href <- webr::eval_js("self.location.href")
  paste0(sub("/shinylive/.*$", "/", href), "data/")
}

pacnet_load <- function(name) {
  if (is_webr()) {
    path <- file.path(tempdir(), "pacnet_data", name)
    if (!file.exists(path)) {
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      url <- paste0(pacnet_data_url(), name)
      ok <- tryCatch(download.file(url, path, mode = "wb", quiet = TRUE) == 0,
                     error = function(e) FALSE)
      if (!ok) {
        unlink(path)
        stop("Could not download reference data from ", url, call. = FALSE)
      }
    }
  } else {
    path <- file.path(getOption("pacnet.data_dir", "../data"), name)
    if (!file.exists(path)) {
      stop("Reference data file not found: ", normalizePath(path, mustWork = FALSE), call. = FALSE)
    }
  }
  readRDS(path)
}
