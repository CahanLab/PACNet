test_that("pacnet_load reads a data file under native R", {
  ig <- pacnet_load("liver_studies_iGenes.rds")
  expect_type(ig, "character")
  expect_gt(length(ig), 1000)
})

test_that("pacnet_load names the missing file", {
  expect_error(pacnet_load("no_such_file.rds"), "no_such_file.rds")
})

test_that("pacnet_data_url honours the override option", {
  withr::local_options(pacnet.data_url = "https://example.org/x/data/")
  expect_identical(pacnet_data_url(), "https://example.org/x/data/")
})

test_that("native R is not detected as webR", {
  expect_false(is_webr())
})
