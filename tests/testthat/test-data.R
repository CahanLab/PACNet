panel <- c("heart", "hspc", "intestine_colon", "liver", "lung", "neuron", "skeletal_muscle")
other <- c("esc", "endothelial_cell", "fibroblast", "kidney", "b_cell", "t_cell", "monocyte_macrophage")
data_dir <- file.path(repo_root, "data")
rd <- function(name) readRDS(file.path(data_dir, name))

test_that("all expected data files exist", {
  expected <- c(
    as.vector(outer(panel, c("_broadClassifier100.rds", "_studies_iGenes.rds",
                             "_engineeredRef_sampTab_all.rds", "_engineeredRef_ranks.rds",
                             "_top_TF_display_order.rds", "_nis.rds"), paste0)),
    paste0(other, "_nis.rds"),
    "none_of_the_above_broadClassifier100.rds", "none_of_the_above_studies_iGenes.rds",
    "Hs_stTrain.rds", "Hs_expTrain_sig4.rds")
  expect_setequal(list.files(data_dir, pattern = "\\.rds$"), expected)
})

test_that("reference panels are integer ranks aligned to their sample tables", {
  for (t in panel) {
    rk <- rd(paste0(t, "_engineeredRef_ranks.rds"))
    st <- rd(paste0(t, "_engineeredRef_sampTab_all.rds"))
    expect_type(rk, "integer")
    expect_identical(colnames(rk), as.character(st$sra_id))
  }
})

test_that("NIS objects are internally consistent", {
  for (t in c(panel, other)) {
    n <- rd(paste0(t, "_nis.rds"))
    expect_named(n, c("tfTargets", "netGenes", "tVals"))
    expect_true(length(n$tfTargets) > 0)
    expect_identical(names(n$tVals$mean), n$netGenes)
    expect_identical(names(n$tVals$sd), n$netGenes)
  }
})

test_that("training matrix is rounded to 4 significant digits", {
  ex <- rd("Hs_expTrain_sig4.rds")
  expect_equal(dim(ex), c(34934L, 1003L))
  expect_identical(ex[1:50, 1:5], signif(ex[1:50, 1:5], 4))
})
