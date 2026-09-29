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
