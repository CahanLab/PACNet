load_ui <- function() {
  local({
    owd <- setwd(file.path(repo_root, "webApp")); on.exit(setwd(owd))
    source("ui.R", local = TRUE)$value
  })
}
ui_html <- function() paste(as.character(htmltools::renderTags(load_ui())$html), collapse = "\n")
count <- function(pattern, x) lengths(regmatches(x, gregexpr(pattern, x, fixed = TRUE)))

test_that("the page loads plotly.js and crosstalk up front", {
  # Under Shinylive, plotly widgets that arrive in the same update rendered
  # before dynamically injected plotly.js had loaded ("Plotly is not defined").
  deps <- vapply(htmltools::findDependencies(load_ui()), function(d) d$name, character(1))
  expect_true(all(c("plotly-main", "crosstalk", "typedarray") %in% deps))
})

test_that("the navbar offers Analyze, Help and About, in that order", {
  h <- ui_html()
  pos <- vapply(c('data-value="Analyze"', 'data-value="Help"', 'data-value="About"'),
                function(p) regexpr(p, h, fixed = TRUE)[[1]], numeric(1))
  expect_true(all(pos > 0))
  expect_true(all(diff(pos) > 0))
})

test_that("links point to CahanLab and the lab website, with the current copyright", {
  h <- ui_html()
  expect_false(grepl("pcahan1", h, fixed = TRUE))
  expect_false(grepl("HTML5 UP", h, fixed = TRUE))
  for (url in c("https://cahanlab.org/", "https://github.com/CahanLab/PACNet",
                "https://github.com/CahanLab/CellNet", "https://doi.org/10.1016/j.stemcr.2023.06.008",
                "https://doi.org/10.5281/zenodo.18857326")) {
    expect_match(h, url, fixed = TRUE, label = url)
  }
  expect_match(h, "© 2019–2026", fixed = TRUE)
})

test_that("the help page explains the inputs and links the example data", {
  h <- ui_html()
  expect_match(h, "example_sample_metadata_table.csv", fixed = TRUE)
  expect_match(h, "example_counts_matrix.csv", fixed = TRUE)
  expect_match(h, "description1", fixed = TRUE)
  expect_false(grepl("Study name", h, fixed = TRUE))
})

test_that("the static inputs carry tooltips", {
  h <- ui_html()
  expect_gte(count("<bslib-tooltip", h), 2)  # action choice + submit
})

test_that("the old HTML5 UP template assets are gone", {
  www <- list.files(file.path(repo_root, "webApp", "www"), recursive = TRUE)
  expect_setequal(www, "images/PACNet_icon.jpg")
})
