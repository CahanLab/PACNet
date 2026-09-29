panel <- c("heart", "hspc", "intestine_colon", "liver", "lung", "neuron", "skeletal_muscle")
other <- c("esc", "endothelial_cell", "fibroblast", "kidney", "b_cell", "t_cell", "monocyte_macrophage")
golden_exists <- function(name) file.exists(file.path(repo_root, "tests", "golden", paste0(name, ".rds")))
eq <- example_query()

test_that("duplicate gene symbols keep the row with the highest median", {
  f <- withr::local_tempfile(fileext = ".csv")
  writeLines(c(",s1,s2,s3", "A,1,1,1", "B,5,5,5", "A,9,9,9", "C,2,2,2"), f)
  x <- read_expression_csv(f)
  expect_setequal(rownames(x), c("A", "B", "C"))
  expect_equal(unlist(x["A", ]), c(s1 = 9, s2 = 9, s3 = 9))
})

test_that("normalize_query orders columns by the metadata sample order", {
  st <- eq$sampTab[rev(seq_len(nrow(eq$sampTab))), ]
  x <- normalize_query(read_expression_csv(example_file("example_counts_matrix.csv")), st)
  expect_identical(colnames(x), st$sample_name)
})

test_that("other tissues use the none_of_the_above classifier", {
  expect_identical(load_classifier("kidney"), load_classifier("none_of_the_above"))
  expect_false(identical(load_classifier("liver")$iGenes, load_classifier("kidney")$iGenes))
})

test_that("query classification matches the original app", {
  for (src in c(panel, "none_of_the_above")) {
    name <- paste0("classify_query_", src)
    if (!golden_exists(name)) next
    clf <- load_classifier(src)
    expect_false(needs_retrain(clf$broadReturn, eq$exp))
    got <- classify(clf$broadReturn, eq$exp, nrand = 3)
    expect_equal(real_cols(got, colnames(eq$exp)), golden(name), label = src)
  }
})

test_that("engineered reference classification matches the original app", {
  for (t in panel) {
    ref <- load_engineered_ref(t)
    got <- classify(load_classifier(t)$broadReturn, ref$expDat, nrand = 2)
    expect_equal(real_cols(got, colnames(ref$expDat)), golden(paste0("classify_ref_", t)), label = t)
  }
})

test_that("NIS matches the original app for every tissue", {
  for (t in c(panel, other)) {
    iGenes <- load_classifier(t)$iGenes
    expect_equal(run_nis(eq$exp, eq$sampTab, iGenes, t), golden(paste0("nis_", t)), label = t)
  }
})

test_that("NIS matches the original app when network genes are missing", {
  g <- golden("nis_liver_subset")
  expect_equal(run_nis(eq$exp, eq$sampTab, g$iGenes, "liver"), g$nis)
})

test_that("retraining stops with a clear message when no genes are shared", {
  x <- eq$exp
  rownames(x) <- paste0("NOT_A_GENE_", seq_len(nrow(x)))
  expect_error(retrain_classifier(x, load_classifier("liver")$iGenes),
               "does not share any genes")
})

test_that("retraining on rounded training data agrees with the original", {
  skip_if_not(identical(Sys.getenv("PACNET_SLOW"), "1"), "set PACNET_SLOW=1 to run")
  g <- golden("retrain_liver")
  res <- retrain_classifier(eq$exp, g$iGenes)
  q <- real_cols(classify(res$broadReturn, eq$exp, nrand = 2), colnames(eq$exp))
  ref <- load_engineered_ref("liver")
  r <- real_cols(classify(res$broadReturn, ref$expDat, nrand = 2), colnames(ref$expDat))
  top <- function(m) rownames(m)[apply(m, 2, which.max)]
  # Rounding may flip a genuine near-tie (e.g. rand 0.410 vs liver 0.396);
  # calls with a clear winner in the original must not change.
  margin <- apply(g$query, 2, function(x) -diff(sort(x, decreasing = TRUE)[1:2]))
  clear <- margin >= 0.05
  expect_gt(sum(clear), 0)
  expect_identical(top(q)[clear], top(g$query)[clear])
  expect_lt(max(abs(q - g$query)), 0.15)
  expect_lt(max(abs(r - g$ref)), 0.15)
  expect_identical(dim(res$classMatrix)[1], 15L)
  expect_true("Rand" %in% rownames(res$classMatrix))
})
