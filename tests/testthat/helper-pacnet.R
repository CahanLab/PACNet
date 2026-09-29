# Shared setup for the PACNet web-app tests. Run from the repo root with:
#   Rscript -e 'testthat::test_dir("tests/testthat")'
# Packages are attached by webApp/global.R, exactly as in the app.
repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
app_file <- function(f) file.path(repo_root, "webApp", f)
example_file <- function(f) file.path(repo_root, "example_data", f)
golden <- function(name) readRDS(file.path(repo_root, "tests", "golden", paste0(name, ".rds")))
options(pacnet.data_dir = file.path(repo_root, "data"))

local({
  owd <- setwd(file.path(repo_root, "webApp")); on.exit(setwd(owd))
  sys.source("global.R", envir = globalenv())
})

example_query <- function() {
  st <- read.csv(example_file("example_sample_metadata_table.csv"), stringsAsFactors = FALSE)
  rownames(st) <- st$sample_name
  list(sampTab = st,
       exp = normalize_query(read_expression_csv(example_file("example_counts_matrix.csv")), st))
}
real_cols <- function(m, cols) m[, cols, drop = FALSE]
