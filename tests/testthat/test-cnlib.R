test_that("cnlib has no parallel or snow code paths", {
  src <- readLines(app_file("cnlib.R"))
  expect_false(any(grepl("detectCores|snow::", src)))
})

test_that("quick-pairs scoring gives the same pairs whether or not it is sliced", {
  set.seed(1)
  genes <- paste0("g", 1:30)
  expDat <- matrix(rexp(30 * 40), nrow = 30, dimnames = list(genes, paste0("s", 1:40)))
  labels <- rep(c("a", "b"), each = 20)
  cg <- list(a = genes[1:15], b = genes[10:30])
  one_slice <- ptGetTop(expDat, labels, cg, topX = 10, sliceSize = 1e6, quickPairs = TRUE)
  sliced <- ptGetTop(expDat, labels, cg, topX = 10, sliceSize = 7, quickPairs = TRUE)
  expect_identical(sliced, one_slice)
})

test_that("every global referenced by cnlib resolves", {
  env <- new.env()
  sys.source(app_file("cnlib.R"), envir = env)
  globals <- unique(unlist(lapply(ls(env), function(f)
    codetools::findGlobals(get(f, env), merge = FALSE)$functions)))
  unresolved <- globals[!vapply(globals, function(g) exists(g, envir = env), logical(1))]
  expect_equal(unresolved, character(0))
})

test_that("sliced pair scoring reproduces sc_testPattern on all pairs at once", {
  set.seed(2)
  genes <- paste0("g", 1:25)
  expDat <- matrix(rexp(25 * 60), nrow = 25, dimnames = list(genes, paste0("s", 1:60)))
  pattern <- sc_sampR_to_pattern(rep(c("a", "b", "c"), each = 20))[["b"]]
  pairTab <- makePairTab(genes)
  whole <- sc_testPattern(pattern, ptSmall(expDat, pairTab))
  for (size in c(1, 13, 100, 1e6)) {
    expect_identical(sc_testPattern_sliced(pattern, expDat, pairTab, size), whole, label = paste("slice", size))
  }
})
