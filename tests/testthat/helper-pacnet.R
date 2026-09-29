# Shared setup for the PACNet web-app tests. Run from the repo root with:
#   Rscript -e 'testthat::test_dir("tests/testthat")'
suppressPackageStartupMessages({
  library(randomForest)
  library(ggplot2)
})
repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
app_file <- function(f) file.path(repo_root, "webApp", f)
example_file <- function(f) file.path(repo_root, "example_data", f)
golden <- function(name) readRDS(file.path(repo_root, "tests", "golden", paste0(name, ".rds")))
options(pacnet.data_dir = file.path(repo_root, "data"))

source(app_file("cnlib.R"))
