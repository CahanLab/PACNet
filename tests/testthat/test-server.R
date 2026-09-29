server_fn <- local({
  source(app_file("server.R"), local = TRUE)
  server
})

upload <- function(f) list(name = basename(f), datapath = example_file(f), size = file.size(example_file(f)), type = "text/csv")

test_that("CellNet analysis runs end to end for a panel tissue", {
  shiny::testServer(server_fn, {
    session$setInputs(actionChoice = "CellNet analysis", species = "Human", tissueType = "liver",
                      sampTabUploadCSV = upload("example_sample_metadata_table.csv"),
                      expMatUpload = upload("example_counts_matrix.csv"))
    session$setInputs(submit = 1)
    expect_match(output$iGenes, "Retraining not necessary")
    expect_match(output$engineeredDescrip, "compare to publicly available data")
    # Rendering the plots catches masking problems (e.g. ggplot2::margin)
    expect_no_error(output$classHm)
    expect_no_error(output$engineeredRef)
    expect_no_error(output$NIS)
  })
})

test_that("CellNet analysis runs for an other-tissue choice", {
  shiny::testServer(server_fn, {
    session$setInputs(actionChoice = "CellNet analysis", species = "Human", tissueType = "kidney",
                      sampTabUploadCSV = upload("example_sample_metadata_table.csv"),
                      expMatUpload = upload("example_counts_matrix.csv"))
    session$setInputs(submit = 1)
    expect_match(output$engineeredDescrip, "unavailable for selected cell/tissue type")
  })
})

test_that("browse mode renders a reference panel description", {
  shiny::testServer(server_fn, {
    session$setInputs(actionChoice = "Browse reference data", species = "Human", tissueTypeRef = "lung")
    session$setInputs(submit = 1)
    expect_match(output$engineeredDescrip, "Reference datasets for engineered lung studies")
    expect_no_error(output$engineeredRef)
  })
})

test_that("a missing reference file is reported to the user instead of hanging", {
  tmp <- withr::local_tempdir()
  for (f in setdiff(list.files(file.path(repo_root, "data")), "lung_engineeredRef_ranks.rds")) {
    file.symlink(file.path(repo_root, "data", f), file.path(tmp, f))
  }
  withr::local_options(pacnet.data_dir = tmp)
  shown <- NULL
  local_mocked_bindings(shinyalert = function(title, text = "", ...) shown <<- text, .package = "shinyalert")
  shiny::testServer(server_fn, {
    session$setInputs(actionChoice = "Browse reference data", species = "Human", tissueTypeRef = "lung")
    expect_no_error(session$setInputs(submit = 1))
  })
  expect_match(shown, "lung_engineeredRef_ranks.rds", fixed = TRUE)
})

# ---- user input errors: one clear alert, spinner stopped --------------------
upload_path <- function(p) list(name = basename(p), datapath = p, size = file.size(p), type = "text/csv")

run_submit <- function(sampTab_path, exp_path, tissue = "liver") {
  seen <- list(alerts = character(0), hidden = 0)
  local_mocked_bindings(shinyalert = function(title, text = "", ...) seen$alerts <<- c(seen$alerts, text),
                        .package = "shinyalert", .env = parent.frame())
  local_mocked_bindings(hide_spinner = function(...) seen$hidden <<- seen$hidden + 1,
                        .package = "shinybusy", .env = parent.frame())
  shiny::testServer(server_fn, {
    session$setInputs(actionChoice = "CellNet analysis", species = "Human", tissueType = tissue,
                      sampTabUploadCSV = upload_path(sampTab_path), expMatUpload = upload_path(exp_path))
    session$setInputs(submit = 1)
  })
  seen
}

test_that("a query sharing no genes with training gets one alert and the spinner stops", {
  x <- read.csv(example_file("example_counts_matrix.csv"), check.names = FALSE)
  x[[1]] <- paste0("NOT_A_GENE_", seq_len(nrow(x)))
  f <- withr::local_tempfile(fileext = ".csv"); write.csv(x, f, row.names = FALSE)
  seen <- run_submit(example_file("example_sample_metadata_table.csv"), f)
  expect_length(seen$alerts, 1)
  expect_match(seen$alerts, "does not share any genes")
  expect_gte(seen$hidden, 1)
})

test_that("a metadata table without description1 stops the spinner", {
  st <- read.csv(example_file("example_sample_metadata_table.csv"))
  st$description1 <- NULL
  f <- withr::local_tempfile(fileext = ".csv"); write.csv(st, f, row.names = FALSE)
  seen <- run_submit(f, example_file("example_counts_matrix.csv"))
  expect_match(seen$alerts, "description1")
  expect_gte(seen$hidden, 1)
})

test_that("a metadata sample missing from the expression matrix is named in the alert", {
  st <- read.csv(example_file("example_sample_metadata_table.csv"))
  st <- rbind(st, transform(st[1, ], sample_name = "EXTRA_SAMPLE"))
  f <- withr::local_tempfile(fileext = ".csv"); write.csv(st, f, row.names = FALSE)
  seen <- run_submit(f, example_file("example_counts_matrix.csv"))
  expect_length(seen$alerts, 1)
  expect_match(seen$alerts, "EXTRA_SAMPLE")
  expect_gte(seen$hidden, 1)
})

test_that("a classification failure is reported once and stops the analysis", {
  real_classify <- get("classify", envir = globalenv())
  assign("classify", function(...) stop("boom"), envir = globalenv())
  withr::defer(assign("classify", real_classify, envir = globalenv()))
  seen <- run_submit(example_file("example_sample_metadata_table.csv"), example_file("example_counts_matrix.csv"))
  expect_length(seen$alerts, 1)
  expect_match(seen$alerts, "Something went wrong with classification of your samples")
  expect_gte(seen$hidden, 1)
})

test_that("server-rendered inputs carry tooltips", {
  shiny::testServer(server_fn, {
    session$setInputs(actionChoice = "CellNet analysis", species = "Human")
    opts <- output$appOptions$html
    expect_gte(lengths(regmatches(opts, gregexpr("<bslib-tooltip", opts, fixed = TRUE))), 3)
    expect_match(opts, "sample_name", fixed = TRUE)
    expect_match(output$tissueType$html, "<bslib-tooltip", fixed = TRUE)
    session$setInputs(actionChoice = "Browse reference data")
    expect_match(output$tissueTypeRef$html, "<bslib-tooltip", fixed = TRUE)
  })
})
