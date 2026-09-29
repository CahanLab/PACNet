test_that("the page loads plotly.js and crosstalk up front", {
  # Under Shinylive, plotly widgets that arrive in the same update rendered
  # before dynamically injected plotly.js had loaded ("Plotly is not defined").
  ui <- local({
    owd <- setwd(file.path(repo_root, "webApp")); on.exit(setwd(owd))
    source("ui.R", local = TRUE)$value
  })
  deps <- vapply(htmltools::findDependencies(ui), function(d) d$name, character(1))
  expect_true(all(c("plotly-main", "crosstalk", "typedarray") %in% deps))
})
