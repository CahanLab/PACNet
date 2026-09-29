# PACNet Serverless Web App Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run the PACNet Shiny web app entirely in the browser (Shinylive/webR) from GitHub Pages, with slimmed reference data, no CellNet/cancerCellNet package dependency, and regression tests proving results match the original app.

**Architecture:** Analysis code moves out of `server.R` into Shiny-free `pipeline.R`; the cancerCellNet functions it needs are vendored into `cnlib.R` with two webR patches; all reference data is read through `pacnet_load()` (local file under R, HTTP fetch under webR). A conversion tool produces slim `.rds` files plus golden outputs from the original code; a GitHub Action exports the app with shinylive and publishes it with the data to Pages.

**Tech Stack:** R (≥ 4.4), shiny + shinyjs/shinyWidgets/shinyalert/shinybusy/shinycssloaders/shinythemes, plotly, ggplot2, pheatmap, randomForest, testthat, shinylive (R package 0.5.0), GitHub Actions, GitHub Pages, Playwright (for browser verification).

**Spec:** `docs/superpowers/specs/2026-09-29-shinylive-port-design.md`

## Global Constraints

- Numerical results for classification and NIS must be identical to the original app on the original data (compare with `expect_equal`, default tolerance).
- Retraining uses `Hs_expTrain_sig4.rds` (`signif(x, 4)`); retrained query top calls must equal the original-data retrain; max |score difference| < 0.15.
- Tissue groups: panel `T` = heart, hspc, intestine_colon, liver, lung, neuron, skeletal_muscle; other `O` = esc, endothelial_cell, fibroblast, kidney, b_cell, t_cell, monocyte_macrophage; plus `none_of_the_above`. Groups `O` and `none_of_the_above` use the `none_of_the_above` classifier/iGenes; group `O` NIS comes from the none_of_the_above GRN with `subnet = O`; NIS is skipped for `none_of_the_above`.
- All data files: `saveRDS(..., compress = "xz")`, in `data/` at repo root (never inside `webApp/`).
- No `library(CellNet)` / `library(cancerCellNet)` anywhere in `webApp/`.
- Vendored code changes are marked with `# PACNet-web:` comments.
- Original code baseline for golden outputs: repo commit `748297c`; cancerCellNet commit `dc29a8e7982461416a12fd3372bb844d35bd0ee0`.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Work on branch `shinylive-port`. Never `git add -A`: the working tree has unrelated deleted `.DS_Store` files and a `.playwright-mcp/` tool directory that must not be committed.

**Environment for this workspace (no system R on this Mac):**

```bash
export SCRATCH=/private/tmp/claude-501/-Users-pcahan-Dropbox--Personal--Code-archive-PACNet/0361f8ee-8262-4ba2-8f0d-683309e783b1/scratchpad
export PATH="$SCRATCH/renv/bin:$PATH"   # conda R 4.5 with all packages + shinylive
cd "/Users/pcahan/Dropbox (Personal)/Code/archive/PACNet"
```

`Rscript` below always means that R. A checkout of cancerCellNet is at `$SCRATCH/spike/cancerCellNet`, and the original S3 `.rda` files are already downloaded in `$SCRATCH/spike/data/`.

## Review Focus

1. Expression CSV with duplicate gene symbols: the row with the highest median is kept (server behavior today). Test in Task 4 (`read_expression_csv`).
2. Metadata rows in a different order than expression columns: `normalize_query` must reorder columns to `sample_name` order. Test in Task 4.
3. Query sharing zero genes with the training data: retraining must stop with a clear message rather than crash deep inside randomForest. Test in Task 4 (`retrain_classifier` with no shared genes).
4. Other-tissue choices (e.g. kidney) must use the none_of_the_above classifier and their own NIS sub-network. Tests in Task 4 (`load_classifier("kidney")`, NIS golden for all 14 tissues).
5. Missing/unreachable data file: error must name the file or URL. Test in Task 3 (native path); browser path checked in Task 6.

---

### Task 1: Vendor cancerCellNet functions into `webApp/cnlib.R`

**Files:**
- Create: `tools/vendor_cnlib.R`
- Create: `webApp/cnlib.R` (generated, then patched)
- Create: `tests/testthat/helper-pacnet.R`
- Test: `tests/testthat/test-cnlib.R`

**Interfaces:**
- Produces: `webApp/cnlib.R` defining (unchanged signatures from cancerCellNet) `addRandToSampTab, broadClass_predict, broadClass_train, ccn_classAssess, cn_clPerf, cn_eval, downSampleW, findBestPairs, findClassyGenes, getClassGenes, gnrAll, logRank, makeClassifier, makePairTab, plot_class_PRs, ptGetTop, ptSmall, query_transform, randomize, rf_classPredict, sc_compAlpha, sc_compMu, sc_cov, sc_fano, sc_filterGenes, sc_sampR_to_pattern, sc_statTab, sc_testPattern, splitCommon_proportion, trans_prop, weighted_down`.
- Produces: test helpers `repo_root`, `app_file(f)`, `golden(name)`, `example_file(f)`.

- [ ] **Step 1: Write the vendoring tool**

`tools/vendor_cnlib.R`:

```r
#!/usr/bin/env Rscript
# Writes webApp/cnlib.R: the verbatim source of every cancerCellNet function
# the PACNet web app calls, directly or indirectly. Web-specific patches are
# applied by hand afterwards and marked with "# PACNet-web:" comments.
#
# Usage (from repo root): Rscript tools/vendor_cnlib.R <cancerCellNet checkout> [out]

args <- commandArgs(trailingOnly = TRUE)
ccn_dir <- args[1]
out <- if (length(args) > 1) args[2] else "webApp/cnlib.R"

fns <- sort(c(
  "addRandToSampTab", "broadClass_predict", "broadClass_train", "ccn_classAssess",
  "cn_clPerf", "cn_eval", "downSampleW", "findBestPairs", "findClassyGenes",
  "getClassGenes", "gnrAll", "logRank", "makeClassifier", "makePairTab",
  "plot_class_PRs", "ptGetTop", "ptSmall", "query_transform", "randomize",
  "rf_classPredict", "sc_compAlpha", "sc_compMu", "sc_cov", "sc_fano",
  "sc_filterGenes", "sc_sampR_to_pattern", "sc_statTab", "sc_testPattern",
  "splitCommon_proportion", "trans_prop", "weighted_down"))

sha <- system2("git", c("-C", shQuote(ccn_dir), "rev-parse", "HEAD"), stdout = TRUE)

# Collate order matches R package default (alphabetical); later definitions win.
chunks <- list()
for (f in sort(list.files(file.path(ccn_dir, "R"), pattern = "\\.R$", full.names = TRUE))) {
  exprs <- parse(f, keep.source = TRUE)
  refs <- attr(exprs, "srcref")
  for (i in seq_along(exprs)) {
    e <- exprs[[i]]
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]])) {
      name <- as.character(e[[2]])
      if (name %in% fns) chunks[[name]] <- paste(as.character(refs[[i]]), collapse = "\n")
    }
  }
}
missing <- setdiff(fns, names(chunks))
if (length(missing)) stop("Not found in ", ccn_dir, ": ", paste(missing, collapse = ", "))

header <- c(
  "# Vendored from https://github.com/pcahan1/cancerCellNet (MIT License),",
  paste0("# commit ", sha, ", by tools/vendor_cnlib.R."),
  "# Only the functions the PACNet web app calls are included. Local changes",
  "# are marked with \"# PACNet-web:\" comments.",
  "")
writeLines(c(header, unlist(lapply(fns, function(n) c(chunks[[n]], "")))), out)
cat("wrote", out, "with", length(fns), "functions\n")
```

- [ ] **Step 2: Generate and commit the verbatim copy (reviewable baseline)**

Run: `Rscript tools/vendor_cnlib.R "$SCRATCH/spike/cancerCellNet"`
Expected: `wrote webApp/cnlib.R with 31 functions`

```bash
git add tools/vendor_cnlib.R webApp/cnlib.R
git commit -m "Vendor cancerCellNet functions used by the web app (verbatim)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 3: Write test helper and the failing cnlib tests**

`tests/testthat/helper-pacnet.R`:

```r
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
```

`tests/testthat/test-cnlib.R`:

```r
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
```

- [ ] **Step 4: Run tests to verify they fail**

Run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: FAIL in "cnlib has no parallel or snow code paths" (verbatim copy still contains `detectCores` and `snow::`). The sliced test may already pass (upstream quick path is a single slice regardless); the globals test may fail if `snow::` counts, which is fine.

- [ ] **Step 5: Patch `findClassyGenes` in `webApp/cnlib.R`**

Replace this block inside `findClassyGenes`:

```r
  ncores = parallel::detectCores(logical = TRUE) # detect the number of cores in the system
  mcCores = 1
  if(ncores * coreProportion >1){
    mcCores = round(ncores * coreProportion)
  }

  cat(ncores, "logical cores in total", "  --> ", mcCores, "cores running to find classification genes...","\n")
  xdiff = gnrAll(expDat, grps, sliceSize, coreProportion)

  if(mcCores == 1) {
    cgenes = lapply(X = xdiff, FUN = getClassGenes, topX = topX)

  }
  else {
    cl = snow::makeCluster(mcCores, type="SOCK")
    cgenes = snow::parLapply(cl = cl, x = xdiff, fun = getClassGenes, topX = topX)
    snow::stopCluster(cl)
  }
```

with:

```r
  # PACNet-web: serial only. parallel::detectCores() is NA under webR and the
  # browser cannot start snow clusters; upstream defaults to serial anyway.
  xdiff = gnrAll(expDat, grps, sliceSize, coreProportion)
  cgenes = lapply(X = xdiff, FUN = getClassGenes, topX = topX)
```

- [ ] **Step 6: Patch `gnrAll`**

Replace the whole `gnrAll` definition with:

```r
gnrAll<-function(expDat, cellLabels, sliceSize, coreProportion) {
  # PACNet-web: serial only (see findClassyGenes).
  nGenes = nrow(expDat)
  cat("nGenes = ",nGenes,"\n")
  str = 1
  stp = min(c(sliceSize, nGenes)) # detect what is smaller the slice size or nGenes

  myPatternG = sc_sampR_to_pattern(as.character(cellLabels))
  statList = list()
  grps = unique(cellLabels)

  for(grp in grps){
    statList[[grp]] = data.frame()
  }

  while(str <= nGenes){
    if(stp>nGenes){
      stp = nGenes
    }

    tempExpDat = expDat[str:stp, ]
    tmpAns = lapply(X = myPatternG, FUN = sc_testPattern, expDat=tempExpDat)
    cat(str,"-", stp,"\n")

    for(gi in seq(length(myPatternG))){
      grp = grps[[gi]]
      statList[[grp]] = rbind( statList[[grp]],  tmpAns[[grp]])
    }

    str = stp+1
    stp = str + sliceSize - 1
  }

  cat("Done testing\n")
  return(statList)
}
```

- [ ] **Step 7: Patch `ptGetTop`**

Replace the whole `ptGetTop` definition (keep any comments/roxygen above it) with:

```r
ptGetTop<-function(expDat, cell_labels, cgenes_list=NA, topX=50, sliceSize = 5e3, quickPairs = FALSE, coreProportion = 1/4){
  # PACNet-web: serial only (see findClassyGenes).
  if(!quickPairs){
    ans = list()
    genes = rownames(expDat)

    # make a data frame of pairs of genes that will be sliced later
    pairTab = makePairTab(genes)

    if(topX > nrow(pairTab)) {
      stop(paste0("The data set has ", nrow(pairTab), " total combination of gene pairs. Please select a smaller topX."))
    }

    cat("setup ans and make pattern\n")
    grps = unique(cell_labels)
    myPatternG = sc_sampR_to_pattern(as.character(cell_labels))
    statList = list()
    for(grp in grps){
      statList[[grp]] = data.frame()
    }

    cat("make pairDat on slice and test\n")
    nPairs = nrow(pairTab)
    cat("nPairs = ",nPairs,"\n")
    str = 1
    stp = min(c(sliceSize, nPairs)) # detect what is smaller the slice size or npairs

    while(str <= nPairs){
      if(stp>nPairs){
        stp = nPairs
      }
      tmpTab = pairTab[str:stp,]
      tmpPdat = ptSmall(expDat, tmpTab)

      tmpAns = lapply(X = myPatternG, FUN = sc_testPattern, expDat=tmpPdat)
      cat(str,"-", stp,"\n")

      for(gi in seq(length(myPatternG))){
        grp = grps[[gi]]
        statList[[grp]] = rbind( statList[[grp]],  tmpAns[[grp]])
      }

      str = stp+1
      stp = str + sliceSize - 1
    }

    cat("compile results\n")
    for(grp in grps){
      tmpAns = findBestPairs(statList[[grp]], topX)
      ans[[grp]] = tmpAns
    }
    return(ans)
  }
  else {

    myPatternG = sc_sampR_to_pattern(as.character(cell_labels))
    ans = list()
    for(cct in names(cgenes_list)){
      genes = cgenes_list[[cct]]
      pairTab = makePairTab(genes)

      nPairs = nrow(pairTab)
      cat("nPairs = ", nPairs," for ", cct, "\n")

      # PACNet-web: score the pairs in slices. Scoring ~45k pairs at once needs
      # ~235 MB per class and exhausted browser (wasm) memory. findBestPairs()
      # ranks by cval only, so the selected pairs are unchanged.
      parts = list()
      for(str in seq(1, nPairs, by = sliceSize)){
        stp = min(str + sliceSize - 1, nPairs)
        tmpPdat = ptSmall(expDat, pairTab[str:stp, , drop = FALSE])
        parts[[length(parts) + 1]] = sc_testPattern(myPatternG[[cct]], expDat=tmpPdat)
      }

      ans[[cct]] = findBestPairs(do.call(rbind, parts), topX)
    }

    return(ans)
  }
}
```

- [ ] **Step 8: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: all 3 tests PASS. If the globals test lists a package function (e.g. `ggplot`, `randomForest`), add that package's `library()` to the helper only if the app already attaches it (it does for randomForest and ggplot2); anything else is a real missing function and must be vendored (add to `fns` in the tool, regenerate, re-apply Steps 5-7).

- [ ] **Step 9: Commit**

```bash
git add webApp/cnlib.R tests/testthat/helper-pacnet.R tests/testthat/test-cnlib.R
git commit -m "cnlib: serial-only code paths and sliced quick-pairs scoring for webR

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Conversion tool, slim data, and golden outputs

**Files:**
- Create: `tools/make_slim_data.R`
- Create: `data/*.rds` (53 files, generated)
- Create: `tests/golden/*.rds` (generated)
- Test: `tests/testthat/test-data.R`

**Interfaces:**
- Produces data files (names exact): for each panel `T`: `T_broadClassifier100.rds`, `T_studies_iGenes.rds`, `T_engineeredRef_sampTab_all.rds` (rownames = `sra_id`), `T_engineeredRef_ranks.rds` (integer matrix, columns in sample-table order), `T_top_TF_display_order.rds`, `T_nis.rds`; for each other `O`: `O_nis.rds`; `none_of_the_above_broadClassifier100.rds`, `none_of_the_above_studies_iGenes.rds`, `Hs_stTrain.rds`, `Hs_expTrain_sig4.rds`.
- NIS object: `list(tfTargets = <named list TF -> character targets>, netGenes = <character>, tVals = list(mean = <named numeric over netGenes>, sd = <named numeric over netGenes>))`.
- Produces golden files: `classify_query_<T|none_of_the_above>`, `classify_ref_<T>`, `nis_<T|O>`, `nis_liver_subset` (list(iGenes, nis)), `retrain_liver` (list(iGenes, query, ref)). Classification matrices contain real sample columns only.

- [ ] **Step 1: Write the failing data test**

`tests/testthat/test-data.R`:

```r
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
    expect_identical(colnames(rk), st$sra_id)
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
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript -e 'testthat::test_dir("tests/testthat", filter = "data")'`
Expected: FAIL (`data/` does not exist).

- [ ] **Step 3: Write the conversion tool**

`tools/make_slim_data.R`:

```r
#!/usr/bin/env Rscript
# Builds the slim reference data in data/ from the original PACNet web-app
# objects on S3, and records golden outputs of the ORIGINAL pipeline in
# tests/golden/ for the regression tests.
#
# Usage (from repo root):
#   git show 748297c:webApp/pacnet_utils.R > /tmp/pacnet_utils_orig.R
#   Rscript tools/make_slim_data.R <cancerCellNet checkout> /tmp/pacnet_utils_orig.R [download_dir]
#
# The cancerCellNet checkout supplies the original (unpatched) functions.

suppressPackageStartupMessages({ library(randomForest); library(igraph) })
options(timeout = 900)

args <- commandArgs(trailingOnly = TRUE)
ccn_dir <- args[1]
orig_utils <- args[2]
dl_dir <- if (length(args) >= 3) args[3] else file.path(tempdir(), "pacnet_orig")
s3 <- "https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/"
dir.create(dl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)
dir.create(file.path("tests", "golden"), recursive = TRUE, showWarnings = FALSE)

panel <- c("heart", "hspc", "intestine_colon", "liver", "lung", "neuron", "skeletal_muscle")
other <- c("esc", "endothelial_cell", "fibroblast", "kidney", "b_cell", "t_cell", "monocyte_macrophage")
none <- "none_of_the_above"

orig <- function(name) {
  f <- file.path(dl_dir, paste0(name, ".rda"))
  if (!file.exists(f)) download.file(paste0(s3, name, ".rda"), f, mode = "wb", quiet = TRUE)
  e <- new.env(); load(f, e); get(ls(e)[1], e)
}
put <- function(obj, name) saveRDS(obj, file.path("data", paste0(name, ".rds")), compress = "xz")
gold <- function(obj, name) saveRDS(obj, file.path("tests", "golden", paste0(name, ".rds")), compress = "xz")
step <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

# ---- original code ---------------------------------------------------------
for (f in list.files(file.path(ccn_dir, "R"), full.names = TRUE)) source(f)
source(orig_utils)  # original pacnet_nis(), subsetGRNall(), subsetTrainNormParam()

# ---- example query, read exactly as the original server.R did --------------
st <- read.csv("example_data/example_sample_metadata_table.csv", stringsAsFactors = FALSE, header = TRUE)
rownames(st) <- st$sample_name
tempDat <- read.csv("example_data/example_counts_matrix.csv", stringsAsFactors = FALSE,
                    header = TRUE, row.names = NULL, check.names = FALSE)
gene_names <- tempDat[, 1]
tempDat <- tempDat[, 2:ncol(tempDat)]
dup_counts <- table(gene_names)
if (any(dup_counts > 1)) {
  dup_genes <- names(which(dup_counts > 1))
  dup_ind <- which(gene_names %in% dup_genes)
  to_keep <- rep(0, length(dup_genes))
  for (i in seq(length(dup_genes))) {
    x2 <- tempDat[which(gene_names == dup_genes[i]), ]
    to_keep[i] <- as.numeric(names(which.max(apply(x2, 1, median))))
  }
  x_new <- tempDat[setdiff(1:nrow(tempDat), dup_ind), ]
  x_new <- rbind(x_new, tempDat[to_keep, ])
  rownames(x_new) <- c(gene_names[setdiff(1:nrow(tempDat), dup_ind)], dup_genes)
  tempDat <- x_new
} else {
  rownames(tempDat) <- gene_names
}
q <- tempDat[, match(st[, "sample_name"], colnames(tempDat))]
q <- apply(q, 2, downSampleW, 1e5)
q <- log(1 + q)
qr <- as.data.frame(logRank(q, base = 0))

# ---- classifiers + query classification golden -----------------------------
for (src in c(panel, none)) {
  step("classifier", src)
  b <- orig(paste0(src, "_broadClassifier100"))
  put(b, paste0(src, "_broadClassifier100"))
  put(orig(paste0(src, "_studies_iGenes")), paste0(src, "_studies_iGenes"))
  if (all(b$cnProc$cgenes %in% rownames(q))) {
    gold(broadClass_predict(b$cnProc, q, nrand = 3)[, colnames(q)], paste0("classify_query_", src))
  } else {
    step("  example query lacks classifier genes for", src, "- no query golden")
  }
}

# ---- engineered reference panels: integer ranks ----------------------------
for (t in panel) {
  step("engineered ref", t)
  st_ref <- orig(paste0(t, "_engineeredRef_sampTab_all"))
  rownames(st_ref) <- st_ref$sra_id
  r <- as.matrix(orig(paste0(t, "_engineeredRef_normalized_expDat_all")))[, rownames(st_ref)]
  rk <- matrix(0L, nrow(r), ncol(r), dimnames = dimnames(r))
  for (j in seq_len(ncol(r))) rk[, j] <- as.integer(rank(r[, j], ties.method = "min"))
  put(rk, paste0(t, "_engineeredRef_ranks"))
  put(st_ref, paste0(t, "_engineeredRef_sampTab_all"))
  put(orig(paste0(t, "_top_TF_display_order")), paste0(t, "_top_TF_display_order"))
  b <- orig(paste0(t, "_broadClassifier100"))
  gold(broadClass_predict(b$cnProc, r, nrand = 2)[, colnames(r)], paste0("classify_ref_", t))
}

# ---- NIS: slim objects + golden from the original server logic -------------
slim_nis <- function(g, tn, t) {
  ng <- names(g$ctGRNs$geneLists[[t]])
  list(tfTargets = g$ctGRNs$tfTargets[[t]], netGenes = ng,
       tVals = list(mean = tn$tVals[[t]]$mean[ng], sd = tn$tVals[[t]]$sd[ng]))
}
orig_nis <- function(g, tn, iGenes, t) {  # server.R lines 676-804 at 748297c
  if (!all(V(g$overallGRN$graph)$name %in% iGenes)) {
    g <- subsetGRNall(g, iGenes)
    tn <- subsetTrainNormParam(tn, g, iGenes)
  }
  pacnet_nis(expDat = qr, stQuery = st, iGenes = iGenes, grnAll = g, trainNorm = tn,
             subnet = t, ctt = t, colname_sid = "sample_name", relaWeight = 0)
}
nis_groups <- c(setNames(as.list(panel), panel), list(none_of_the_above = other))
for (src in names(nis_groups)) {
  step("GRN", src)
  g <- orig(paste0(src, "_grnAll"))
  tn <- orig(paste0(src, "_trainNormParam"))
  ig <- orig(paste0(src, "_studies_iGenes"))
  for (t in nis_groups[[src]]) {
    step("  NIS", t)
    put(slim_nis(g, tn, t), paste0(t, "_nis"))
    gold(orig_nis(g, tn, ig, t), paste0("nis_", t))
  }
  if (src == "liver") {
    ig2 <- setdiff(ig, head(names(g$ctGRNs$geneLists$liver), 25))
    gold(list(iGenes = ig2, nis = orig_nis(g, tn, ig2, "liver")), "nis_liver_subset")
  }
  rm(g, tn); invisible(gc())
}

# ---- training data + original-data retrain golden --------------------------
step("training data")
expTrain <- orig("Hs_expTrain_Jun-20-2017")
stTrain <- orig("Hs_stTrain_Jun-20-2017")
put(signif(expTrain, 4), "Hs_expTrain_sig4")
put(stTrain, "Hs_stTrain")

step("retrain (original data)")
b <- orig("liver_broadClassifier100")
ig_r <- setdiff(intersect(orig("liver_studies_iGenes"), rownames(q)), head(b$cnProc$cgenes, 50))
set.seed(39)
stList <- splitCommon_proportion(sampTab = stTrain, proportion = 0.66, dLevel = "description1")
br <- broadClass_train(stTrain = stList$trainingSet,
                       expTrain = expTrain[, stList$trainingSet$sra_id][ig_r, ],
                       colName_cat = "description1", colName_samp = "sra_id",
                       nRand = 70, nTopGenes = 100, nTopGenePairs = 100, nTrees = 2000,
                       stratify = TRUE, sampsize = 25, quickPairs = TRUE)
ref <- readRDS(file.path("data", "liver_engineeredRef_ranks.rds"))
gold(list(iGenes = ig_r,
          query = broadClass_predict(br$cnProc, q, nrand = 2)[, colnames(q)],
          ref = broadClass_predict(br$cnProc, ref, nrand = 2)[, colnames(ref)]),
     "retrain_liver")
step("done")
```

- [ ] **Step 4: Run the tool**

```bash
git show 748297c:webApp/pacnet_utils.R > "$SCRATCH/pacnet_utils_orig.R"
Rscript tools/make_slim_data.R "$SCRATCH/spike/cancerCellNet" "$SCRATCH/pacnet_utils_orig.R" "$SCRATCH/spike/data" 2>&1 | grep -E "^[0-9]{2}:|Error"
```

Expected: step lines through `done`, no `Error`. Takes ~10 min (retrain).
Then: `ls data | wc -l` → `53`; `du -sh data` → about 97M; `ls tests/golden`.

- [ ] **Step 5: Run the data tests**

Run: `Rscript -e 'testthat::test_dir("tests/testthat", filter = "data")'`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add tools/make_slim_data.R data tests/golden tests/testthat/test-data.R
git commit -m "Slim reference data and golden outputs from the original pipeline

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Data access (`pacnet_load`)

**Files:**
- Create: `webApp/data_access.R`
- Modify: `tests/testthat/helper-pacnet.R` (source the new file)
- Test: `tests/testthat/test-data-access.R`

**Interfaces:**
- Produces: `pacnet_load(name)` → object stored in `data/<name>` (`name` includes `.rds`). Options: `pacnet.data_dir` (native; default `"../data"`), `pacnet.data_url` (webR override). Errors: native missing file → message contains the file path; webR download failure → message contains the URL.
- Produces: `is_webr()` → logical; `pacnet_data_url()` → character ending in `data/`.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-data-access.R`:

```r
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
```

Append to `tests/testthat/helper-pacnet.R`:

```r
source(app_file("data_access.R"))
```

- [ ] **Step 2: Run to verify failure**

Run: `Rscript -e 'testthat::test_dir("tests/testthat", filter = "data-access")'`
Expected: FAIL (helper cannot source `data_access.R`).

- [ ] **Step 3: Implement**

`webApp/data_access.R`:

```r
# Reference data access for the PACNet web app.
#
# Native R (runApp / runGitHub): files are read from ../data relative to the
# app directory (option pacnet.data_dir overrides).
# webR (Shinylive, in the browser): files are downloaded from the data/
# directory next to the exported site, once per session, then read.

is_webr <- function() identical(R.version$os, "emscripten")

pacnet_data_url <- function() {
  url <- getOption("pacnet.data_url")
  if (!is.null(url)) return(url)
  # The webR worker script lives at <site>/shinylive/webr/webr-worker.js
  href <- webr::eval_js("self.location.href")
  paste0(sub("/shinylive/.*$", "/", href), "data/")
}

pacnet_load <- function(name) {
  if (is_webr()) {
    path <- file.path(tempdir(), "pacnet_data", name)
    if (!file.exists(path)) {
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      url <- paste0(pacnet_data_url(), name)
      ok <- tryCatch(download.file(url, path, mode = "wb", quiet = TRUE) == 0,
                     error = function(e) FALSE)
      if (!ok) {
        unlink(path)
        stop("Could not download reference data from ", url, call. = FALSE)
      }
    }
  } else {
    path <- file.path(getOption("pacnet.data_dir", "../data"), name)
    if (!file.exists(path)) {
      stop("Reference data file not found: ", normalizePath(path, mustWork = FALSE), call. = FALSE)
    }
  }
  readRDS(path)
}
```

- [ ] **Step 4: Run to verify pass**

Run: `Rscript -e 'testthat::test_dir("tests/testthat", filter = "data-access")'`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add webApp/data_access.R tests/testthat/helper-pacnet.R tests/testthat/test-data-access.R
git commit -m "Add pacnet_load(): local files under R, site-relative download under webR

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Shiny-free analysis pipeline and slim `pacnet_nis`

**Files:**
- Create: `webApp/pipeline.R`
- Modify: `webApp/pacnet_utils.R` (replace whole file)
- Modify: `tests/testthat/helper-pacnet.R`
- Test: `tests/testthat/test-pipeline.R`

**Interfaces:**
- Consumes: `pacnet_load(name)`; cnlib functions.
- Produces (all in `pipeline.R`):
  - `PANEL_TISSUES`, `OTHER_TISSUES` (character vectors).
  - `read_expression_csv(path)` → data.frame, genes × samples, rownames = unique gene symbols.
  - `normalize_query(expDat, sampTab)` → numeric matrix, columns in `sampTab$sample_name` order.
  - `load_classifier(tissue)` → `list(broadReturn, iGenes)`.
  - `needs_retrain(broadReturn, queryExp)` → logical.
  - `retrain_classifier(queryExp, iGenes)` → `list(broadReturn, iGenes, classMatrix, grps, assessment)`.
  - `classify(broadReturn, expDat, nrand)` → matrix (classes × samples incl. `rand_*`).
  - `load_engineered_ref(tissue)` → `list(expDat, sampTab, topTFs)`.
  - `run_nis(queryExp, sampTab, iGenes, tissue)` → matrix TFs × samples.
- Produces (in `pacnet_utils.R`): `pacnet_nis(expDat, stQuery, iGenes, nis, colname_sid = "sra_id", relaWeight = 0)`.

- [ ] **Step 1: Write the failing tests**

Replace the `source(...)` lines at the end of `tests/testthat/helper-pacnet.R` with:

```r
for (f in c("cnlib.R", "data_access.R", "pacnet_utils.R", "pipeline.R")) source(app_file(f))

example_query <- function() {
  st <- read.csv(example_file("example_sample_metadata_table.csv"), stringsAsFactors = FALSE)
  rownames(st) <- st$sample_name
  list(sampTab = st,
       exp = normalize_query(read_expression_csv(example_file("example_counts_matrix.csv")), st))
}
real_cols <- function(m, cols) m[, cols, drop = FALSE]
```

`tests/testthat/test-pipeline.R`:

```r
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
  expect_identical(top(q), top(g$query))
  expect_lt(max(abs(q - g$query)), 0.15)
  expect_lt(max(abs(r - g$ref)), 0.15)
  expect_identical(dim(res$classMatrix)[1], 15L)
  expect_true("Rand" %in% rownames(res$classMatrix))
})
```

- [ ] **Step 2: Run to verify failure**

Run: `Rscript -e 'testthat::test_dir("tests/testthat", filter = "pipeline")'`
Expected: FAIL (helper cannot source `pipeline.R`).

- [ ] **Step 3: Replace `webApp/pacnet_utils.R`**

```r
# Utilities for PACNet webApp only

#' network influence score
#'
#' Computes network influence score (NIS).
#' @param expDat ranked query expression matrix, output of running logRank()
#' @param stQuery query sample metadata table
#' @param iGenes character vector of gene names found in both training data and query data
#' @param nis slim network object for one cell/tissue type: list(tfTargets,
#'   netGenes, tVals = list(mean, sd)); see tools/make_slim_data.R
#' @param colname_sid colname of stQuery with unique sample labels
#' @param relaWeight whether to weight by overall expression such that TFs with higher expression in ctt are more important (1=do the weighting)
#'
#' @return numeric matrix where rows are TFs in CT GRN and columns are query samples. Negative score indicates TF should be upregulated to acquire a ctt fate. Positive score indicates TF should be downregulated to acquire a ctt fate.
pacnet_nis <- function(expDat, stQuery, iGenes, nis, colname_sid="sra_id", relaWeight=0) {
  tfTargList <- nis$tfTargets
  nTargets<-vector();
  targetScore<-vector();
  tfScore<-vector();
  totalScore<-vector();
  tfWeights<-vector();

  tfs<-names(tfTargList);
  netGenes<-intersect(nis$netGenes, iGenes)

  sids<-as.vector(stQuery[,colname_sid])

  ans<-matrix(0, nrow=length(tfs), ncol=nrow(stQuery))
  rownames(ans)<-tfs;
  colnames(ans)<-sids;

  meanVals <- nis$tVals$mean[netGenes]
  sdVals <- nis$tVals$sd[netGenes]

  # compute a matrix of zscores.
  zzzMat<-matrix(0, nrow=length(netGenes), ncol=nrow(stQuery));

  for(i in seq(length(sids))){
    sid<-sids[i];
    print(paste0("Computing zscores for ", sid))
    xvals<-as.vector(expDat[netGenes, sid])
    zzzMat[,i]<-(xvals - meanVals) / sdVals
  }

  rownames(zzzMat)<-netGenes
  colnames(zzzMat)<-stQuery[,colname_sid]

  # assign weights
  meanVect<-unlist(meanVals);
  meanVect <- meanVect / min(meanVect) # Added to make function work with rank-based tVals
  weights<-(2**meanVect)/sum(2**meanVect);

  for(sid in sids) {
    print(paste0("TF scoring for ", sid));

    for(i in seq(length(tfs))){
      tf<-tfs[i];
      # zscore of TF relative to target C/T
      tfScore[i]<-zzzMat[tf,sid];

      targs<-tfTargList[[tf]];
      targs<-intersect(targs, iGenes);

      # Zscores of TF targets, relative to C/T
      tmp<-zzzMat[targs,sid];
      targetScore[i]<-sum(tmp*weights[targs]);

      ## new one:
      totalScore[i]<-targetScore[i] + (length(targs)*tfScore[i]*weights[tf]);

      if(relaWeight!=1){ # don't weight by expression
        meanW<-mean(weights)
        totalScore[i]<- sum(tmp)*meanW + (length(targs)*tfScore[i])*meanW
      }
      nTargets[i]<-length(targs) ;
      tfWeights[i]<-weights[tf];
    }

    xxx<-data.frame(tf=tfs, tfScore=tfScore, targetScore=targetScore, nTargets=nTargets,tfWeight=tfWeights, totalScore=totalScore);
    xxx<-xxx[order(xxx$totalScore),]; # puts the worst ones at top when plotting
    xxx$tf<-factor(xxx$tf, as.vector(unique(xxx$tf)));
    ans[as.vector(xxx$tf),sid]<-as.vector(xxx$totalScore);
  }

  return(ans) # returns network influence score.
}
```

(The per-sample weights did not depend on the sample in the original, so they are computed once; the z-score loop is vectorized. Arithmetic is the same operations in the same order per element; the golden tests enforce equality.)

- [ ] **Step 4: Create `webApp/pipeline.R`**

```r
# Analysis steps of the PACNet web app, free of Shiny so they can be tested
# (tests/testthat) and reused by server.R.

PANEL_TISSUES <- c("heart", "hspc", "intestine_colon", "liver", "lung", "neuron", "skeletal_muscle")
OTHER_TISSUES <- c("esc", "endothelial_cell", "fibroblast", "kidney", "b_cell", "t_cell", "monocyte_macrophage")

classifier_source <- function(tissue) {
  if (tissue %in% PANEL_TISSUES) tissue else "none_of_the_above"
}

# Reads an uploaded expression CSV (first column = gene symbols). Duplicated
# symbols keep the row with the highest median across samples.
read_expression_csv <- function(path) {
  tempDat <- read.csv(path, stringsAsFactors = FALSE, header = TRUE,
                      row.names = NULL, check.names = FALSE)
  gene_names <- tempDat[, 1]
  tempDat <- tempDat[, 2:ncol(tempDat), drop = FALSE]
  dup_counts <- table(gene_names)
  if (any(dup_counts > 1)) {
    dup_genes <- names(which(dup_counts > 1))
    dup_ind <- which(gene_names %in% dup_genes)
    to_keep <- rep(0, length(dup_genes))
    for (i in seq(length(dup_genes))) {
      x2 <- tempDat[which(gene_names == dup_genes[i]), , drop = FALSE]
      to_keep[i] <- as.numeric(names(which.max(apply(x2, 1, median))))
    }
    x_new <- tempDat[setdiff(1:nrow(tempDat), dup_ind), , drop = FALSE]
    x_new <- rbind(x_new, tempDat[to_keep, , drop = FALSE])
    rownames(x_new) <- c(gene_names[setdiff(1:nrow(tempDat), dup_ind)], dup_genes)
    tempDat <- x_new
  } else {
    rownames(tempDat) <- gene_names
  }
  tempDat
}

normalize_query <- function(expDat, sampTab) {
  expDat <- expDat[, match(sampTab[, "sample_name"], colnames(expDat))]
  expDat <- apply(expDat, 2, downSampleW, 1e5)
  log(1 + expDat)
}

load_classifier <- function(tissue) {
  src <- classifier_source(tissue)
  list(broadReturn = pacnet_load(paste0(src, "_broadClassifier100.rds")),
       iGenes = pacnet_load(paste0(src, "_studies_iGenes.rds")))
}

needs_retrain <- function(broadReturn, queryExp) {
  !all(broadReturn$cnProc$cgenes %in% rownames(queryExp))
}

# Retrains the broad classifier on the genes shared with the query and builds
# the validation outputs shown on the "New Training Validation" tab.
retrain_classifier <- function(queryExp, iGenes) {
  iGenes <- intersect(iGenes, rownames(queryExp))
  if (length(iGenes) == 0) {
    stop("Your expression matrix does not share any genes with the training data.", call. = FALSE)
  }
  expTrain <- pacnet_load("Hs_expTrain_sig4.rds")
  stTrain <- pacnet_load("Hs_stTrain.rds")

  set.seed(39) # Set same seed every time for random number generator
  stList <- splitCommon_proportion(sampTab = stTrain, proportion = 0.66, dLevel = "description1")
  stTrainSubset <- stList$trainingSet
  stVal <- stList$validationSet
  stVal <- stVal[order(stVal$description1), ] # order by classification name
  expTrainSubset <- expTrain[iGenes, stTrainSubset$sra_id]
  expVal <- expTrain[iGenes, rownames(stVal)]
  rm(expTrain); invisible(gc()) # free the full matrix before training (browser memory)

  broadReturn <- broadClass_train(stTrain = stTrainSubset, expTrain = expTrainSubset,
                                  colName_cat = "description1", colName_samp = "sra_id",
                                  nRand = 70, nTopGenes = 100, nTopGenePairs = 100,
                                  nTrees = 2000, stratify = TRUE, sampsize = 25,
                                  quickPairs = TRUE)
  rm(expTrainSubset); invisible(gc())

  classMatrix <- broadClass_predict(broadReturn$cnProc, expVal, nrand = 60)
  # Rename "rand" to "Rand" and move it last, for the validation heatmap
  classMatrix <- classMatrix[rownames(classMatrix)[c(1:12, 14:15, 13)], ]
  rownames(classMatrix)[15] <- "Rand"
  stValRand <- addRandToSampTab(classMatrix, stVal, "description1", "sra_id")
  stValRand$description1[which(stValRand$description1 == "rand")] <- "Rand"
  grps <- as.vector(stValRand$description1)
  names(grps) <- rownames(stValRand)

  list(broadReturn = broadReturn, iGenes = iGenes, classMatrix = classMatrix, grps = grps,
       assessment = ccn_classAssess(classMatrix, stValRand, "description1", "sra_id"))
}

classify <- function(broadReturn, expDat, nrand) {
  broadClass_predict(cnProc = broadReturn$cnProc, expDat = expDat, nrand = nrand)
}

load_engineered_ref <- function(tissue) {
  sampTab <- pacnet_load(paste0(tissue, "_engineeredRef_sampTab_all.rds"))
  rownames(sampTab) <- sampTab$sra_id
  list(expDat = pacnet_load(paste0(tissue, "_engineeredRef_ranks.rds"))[, rownames(sampTab)],
       sampTab = sampTab,
       topTFs = pacnet_load(paste0(tissue, "_top_TF_display_order.rds")))
}

run_nis <- function(queryExp, sampTab, iGenes, tissue) {
  nis <- pacnet_load(paste0(tissue, "_nis.rds"))
  # Drop TFs and targets not measured in this study (the original app did this
  # by subsetting the full GRN object).
  keep <- names(nis$tfTargets) %in% iGenes
  nis$tfTargets <- lapply(nis$tfTargets[keep], function(t) t[t %in% iGenes])
  queryRanked <- as.data.frame(logRank(queryExp, base = 0))
  pacnet_nis(expDat = queryRanked, stQuery = sampTab, iGenes = iGenes, nis = nis,
             colname_sid = "sample_name", relaWeight = 0)
}
```

- [ ] **Step 5: Run the fast tests**

Run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: all PASS; the retrain test reports SKIP.

If "NIS matches" fails only by element order of `tfTargets` names vs golden row order, the fix is in `run_nis` (preserve `names(nis$tfTargets)` order, which the code above does) — do not loosen the test.

- [ ] **Step 6: Run the slow retrain test**

Run: `PACNET_SLOW=1 Rscript -e 'testthat::test_dir("tests/testthat", filter = "pipeline")'`
Expected: PASS (takes ~5 min).

- [ ] **Step 7: Commit**

```bash
git add webApp/pipeline.R webApp/pacnet_utils.R tests/testthat/helper-pacnet.R tests/testthat/test-pipeline.R
git commit -m "Shiny-free analysis pipeline; NIS on slim per-tissue network objects

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Rewire the Shiny app onto the pipeline

**Files:**
- Create: `webApp/global.R`
- Modify: `webApp/ui.R:1-20` (library block)
- Modify: `webApp/server.R` (header, `trainClassifier`, upload reader, analysis + browse branches)
- Delete: `webApp/deploy_agnosticCellNet_web.sh`
- Modify: `tests/testthat/helper-pacnet.R`
- Test: `tests/testthat/test-server.R`

**Interfaces:**
- Consumes: everything from Task 4; `pacnet_load`.
- Produces: `server` function (unchanged name) usable by `shiny::testServer`.

- [ ] **Step 1: Write the failing server test**

`tests/testthat/test-server.R`:

```r
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
  })
})
```

Replace the `for (f in c("cnlib.R", ...)) source(app_file(f))` line in `tests/testthat/helper-pacnet.R` with:

```r
local({
  owd <- setwd(file.path(repo_root, "webApp")); on.exit(setwd(owd))
  sys.source("global.R", envir = globalenv())
})
```

- [ ] **Step 2: Run to verify failure**

Run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: FAIL (`global.R` missing).

- [ ] **Step 3: Create `webApp/global.R`**

```r
# Loaded once before ui.R and server.R (locally and under Shinylive).
suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(shinycssloaders)
  library(shinyWidgets)
  library(shinyjs)
  library(shinyalert)
  library(shinybusy)
  library(ggplot2)
  library(gridExtra)
  library(RColorBrewer)
  library(pheatmap)
  library(plotly)
  library(dplyr)
  library(randomForest)
})
source("cnlib.R")
source("data_access.R")
source("pacnet_utils.R")
source("pipeline.R")
source("plotting.R")

options(shiny.maxRequestSize = 50*1024^2)
```

(`igraph` is no longer needed: it was only used for GRN subsetting. Before deleting it, confirm with `grep -n "igraph\|V(\|induced_subgraph\|vertex_attr" webApp/*.R` that no remaining code uses it.)

- [ ] **Step 4: Edit `webApp/ui.R`**

Replace lines 1-20 (the `suppressPackageStartupMessages({ library(CellNet) ... source("pacnet_utils.R") })` block) with nothing — `global.R` now loads packages. The file then begins with the existing commented-out `fluidPage` block and ends with `htmlTemplate(...)`.

- [ ] **Step 5: Edit `webApp/server.R` header and `trainClassifier`**

Delete lines 1-24 (the `suppressPackageStartupMessages({...})` block and `options(shiny.maxRequestSize = ...)`).

Replace the `trainClassifier <- function(...) { ... }` definition (original lines 31-101) with:

```r
  show_validation <- function(res) {
    updateProgressBar(session = session, id = "progress",
                      title = "Creating validation plots...",
                      value = 35, total = 100)

    output$iGenes <- renderText({
      paste("Number of intersecting genes: ", length(res$iGenes))
    })
    outputOptions(output, "iGenes", suspendWhenHidden = FALSE)

    # Classification Validation Heatmap
    Sys.setlocale("LC_COLLATE", "C") # So that sort() is case sensitive
    output$newTrainHm <- renderPlot({
      acn_queryClassHm(res$classMatrix, grps=res$grps, fontsize_row=10)
    }, height=400, width=900)
    outputOptions(output, "newTrainHm", suspendWhenHidden = FALSE)

    # Classifier Assessment
    output$newTrainPR <- renderPlot({
      plot_class_PRs(res$assessment)
    }, height=720, width=720)
    outputOptions(output, "newTrainPR", suspendWhenHidden = FALSE)
  }
```

- [ ] **Step 6: Replace the expression upload reader**

Replace the body of `queryExpDat <- eventReactive(input$expMatUpload, { ... })` with:

```r
  queryExpDat <- eventReactive(input$expMatUpload, {
    tryCatch(
      read_expression_csv(input$expMatUpload$datapath),
      error=function(condition) {
        shinyalert("Oops!", paste0("Something went wrong with your expression matrix upload: ",
                                   condition), type="error")
        return(NULL)
      })
  })
```

- [ ] **Step 7: Replace the CellNet-analysis data steps**

In the `if (actionChoice() == "CellNet analysis") {` branch, replace everything from `new_col_order <- match(...)` through the line `cnProc_broad <- broadReturn$cnProc` / `classMatrixQuery <- broadClass_predict(...)` (original lines 518-603) with:

```r
      queryExpDat_new <- normalize_query(queryExpDat(), querySampTab())

      updateProgressBar(session = session, id = "progress", title = "Loaded in query data...",
                        value = 10, total = 100)

      # Load TISSUE-SPECIFIC classifier and iGenes (other tissues use none_of_the_above)
      clf <- load_classifier(tissueType())
      broadReturn <- clf$broadReturn
      iGenes <- clf$iGenes

      # Check if retraining is necessary
      if (!needs_retrain(broadReturn, queryExpDat_new)) {
        updateProgressBar(session = session, id = "progress",
                          title = "No need to retrain! Starting classification...",
                          value = 35, total = 100)

        output$iGenes <- renderText({
          "Retraining not necessary for this study!"
        })
        outputOptions(output, "iGenes", suspendWhenHidden = FALSE)
      } else {
        updateProgressBar(session = session, id = "progress",
                          title = "Retraining classifier based on intersecting genes...This may take several minutes",
                          value = 20, total = 100)
        iGenes <- intersect(iGenes, rownames(queryExpDat_new))

        if (length(iGenes) == 0) {
          shinyalert("Error!", "Your expression matrix does not share any genes with the training data.
                          Please make sure the row names of your expression matrix are GENE SYMBOLS.", type="error")
        } else if (length(iGenes) < 1000) {
          shinyalert("Warning!", "Your expression matrix shares fewer than 1000 genes with the training data.
                          The most likely cause is that the row names of your expression matrix are not GENE SYMBOLS.
                          If you have made sure that your row names are properly gene symbols and would like to proceed,
                          please note that a classifier trained on fewer than 1000 genes may not be reliable.
                          Check out the New Training Validation tab for information on the classifier's performance.",
                     type="warning",
                     confirmButtonText="I understand")
        }
        showTab(inputId="analysis_tabs", target="newTrainingValidation")
        retrained <- tryCatch(
          retrain_classifier(queryExpDat_new, iGenes),
          error=function(condition) {
            shinyalert("Error", paste0("Something went wrong with retraining: ",
                                       condition), type="error")
            return(NULL)
          }
        )
        if (is.null(retrained)) {
          hide_spinner(spin_id = "overall_spinner", session=session)
          return(NULL)
        }
        broadReturn <- retrained$broadReturn
        show_validation(retrained)
        updateProgressBar(session = session, id = "progress",
                          title = "Done retraining! Starting classification...",
                          value = 40, total = 100)
      }

      # Load tissueType-specific engineered reference data
      if (tissueType() %in% PANEL_TISSUES) {
        ref <- load_engineered_ref(tissueType())
        tissueTypeRefDat <- ref$expDat
        tissueTypeRefSt <- ref$sampTab
        tissueTypeTopTFs <- ref$topTFs
      }

      classMatrixQuery <- classify(broadReturn, queryExpDat_new, nrand = 3)
```

Also replace the two remaining hard-coded panel lists in this branch (`if(tissueType() %in% c("heart","hspc","intestine_colon", "liver","lung","neuron","skeletal_muscle"))` around the engineered-reference heatmap) with `if (tissueType() %in% PANEL_TISSUES)`, and the list `c("esc","endothelial_cell",...,"monocyte_macrophage")` before the NIS plot with `OTHER_TISSUES`.

- [ ] **Step 8: Replace the GRN/NIS data steps**

Replace everything from `###### GRN Analysis ######` through the end of the `pacnet_nis(...)` call's `tryCatch` (original lines 664-810: GRN loading, the subsetting block, `logRank`, commented GRN-status code, and the NIS `tryCatch`) with:

```r
      #####
      # NIS

      updateProgressBar(session = session, id = "progress", title="Now analyzing transcriptional regulators...may take several minutes",
                        value = 85, total = 100)

      if(tissueType() != "none_of_the_above") {

        TF_scores <- tryCatch(
          run_nis(queryExpDat_new, querySampTab(), iGenes, tissueType()),
          error=function(condition) {
            shinyalert("Error", paste0("Something went wrong with NIS analysis of your samples: ",
                                       condition), type="error")
            return(NULL)
          })
```

Keep the following original code unchanged (from `output$nisDownload <- downloadHandler(` onward, including the top-TF selection for `OTHER_TISSUES`, `queryNIS(...)`, and the `else` branch for none_of_the_above).

- [ ] **Step 9: Replace the browse-mode loads**

In the `else if (actionChoice() == "Browse reference data")` branch, replace:

```r
      broadReturn <- utils_loadObject(paste0("data/", tissueTypeRef(), "_broadClassifier100.rda"))
      iGenes <- utils_loadObject(paste0("data/", tissueTypeRef(), "_studies_iGenes.rda"))

      tissueTypeRefDat <- utils_loadObject(paste0("data/", tissueTypeRef(), "_engineeredRef_normalized_expDat_all.rda"))
      tissueTypeRefSt <- utils_loadObject(paste0("data/", tissueTypeRef(), "_engineeredRef_sampTab_all.rda"))
      rownames(tissueTypeRefSt) <- tissueTypeRefSt$sra_id
      tissueTypeRefDat <- tissueTypeRefDat[,rownames(tissueTypeRefSt)]
```

with:

```r
      broadReturn <- load_classifier(tissueTypeRef())$broadReturn
      ref <- load_engineered_ref(tissueTypeRef())
      tissueTypeRefDat <- ref$expDat
      tissueTypeRefSt <- ref$sampTab
```

and delete the line `tissueTypeRefGRNStat <- utils_loadObject(paste0("data/", tissueTypeRef(),"_GRN_statusQuery.rda"))`.

- [ ] **Step 10: Check nothing still reads old data or packages**

Run: `grep -n "utils_loadObject\|\.rda\|data/\|library(CellNet\|library(cancerCellNet\|grnAll\|trainNormParam" webApp/*.R`
Expected: no output.

- [ ] **Step 11: Delete the AWS deploy script**

Run: `git rm webApp/deploy_agnosticCellNet_web.sh`

- [ ] **Step 12: Run all tests**

Run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: all PASS (retrain SKIP). If `testServer` fails inside a UI-only helper (e.g. `hideTab`, `show_spinner`, `delay`) because `MockShinySession` lacks a method, the server code is not at fault: wrap nothing in the app; instead add `local_mocked_bindings(hideTab = function(...) NULL, showTab = function(...) NULL, .package = "shiny")` (or the equivalent for the failing helper's package) at the top of the affected test.

- [ ] **Step 13: Run the app natively and click through once**

Run: `Rscript -e 'shiny::runApp("webApp", port = 7777, launch.browser = FALSE)'` in the background, then with Playwright open `http://127.0.0.1:7777`, upload the two example CSVs, choose liver, Submit, and confirm the classification heatmap, the engineered-reference heatmap and the NIS tab render. Stop the app.

- [ ] **Step 14: Commit**

```bash
git add webApp/global.R webApp/ui.R webApp/server.R tests/testthat/helper-pacnet.R tests/testthat/test-server.R
git commit -m "Web app reads slim data through the pipeline; drop CellNet packages and AWS deploy script

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Static site build and in-browser verification

**Files:**
- Create: `tools/build_site.sh`
- Modify: `.gitignore` (add `site/`)

**Interfaces:**
- Produces: `tools/build_site.sh [outdir]` → exported site with `data/` copied in; used by CI in Task 7.

- [ ] **Step 1: Write the build script**

`tools/build_site.sh`:

```bash
#!/usr/bin/env bash
# Builds the static PACNet web app: a shinylive export of webApp/ plus the
# reference data in data/. Usage (from repo root): tools/build_site.sh [outdir]
set -euo pipefail
out="${1:-site}"
rm -rf "$out"
Rscript -e "options(timeout = 900); shinylive::export('webApp', '$out')"
cp -R data "$out/data"
touch "$out/.nojekyll"
echo "Built $out ($(du -sh "$out" | cut -f1))"
```

Run: `chmod +x tools/build_site.sh && echo "site/" >> .gitignore`

- [ ] **Step 2: Build**

Run: `tools/build_site.sh site`
Expected: `Built site (...)`, roughly 250-300 MB. Check that the bundled package list contains no `igraph` or `snow`: `ls site/shinylive/webr/packages | grep -E "^(igraph|snow)$"` → no output.

- [ ] **Step 3: Serve under a sub-path, as GitHub Pages does**

```bash
mkdir -p "$SCRATCH/srv" && ln -sfn "$PWD/site" "$SCRATCH/srv/PACNet"
python3 -m http.server -d "$SCRATCH/srv" 8766   # background
```

- [ ] **Step 4: Browser check — CellNet analysis, no retrain**

With Playwright: open `http://localhost:8766/PACNet/`, wait for the Submit button (first load ~30-60 s), upload `example_data/example_sample_metadata_table.csv` and `example_data/example_counts_matrix.csv`, choose tissue `liver`, Submit. Wait for progress "Analysis Complete!". Verify: `#classHm .plotly`, `#engineeredRef .plotly` exist; NIS tab `#NIS .plotly` exists; "Download classification results" button visible; console has no new errors other than the pre-existing `terminal.js` `consoleText` error.

- [ ] **Step 5: Browser check — forced retrain, largest gene set**

Create a query missing classifier genes for the none_of_the_above classifier:

```bash
Rscript -e '
cg <- readRDS("data/none_of_the_above_broadClassifier100.rds")$cnProc$cgenes
x <- read.csv("example_data/example_counts_matrix.csv", check.names = FALSE)
x <- x[!(x[[1]] %in% head(cg, 50)), ]
write.csv(x, file.path(Sys.getenv("SCRATCH"), "retrain_query.csv"), row.names = FALSE)'
```

(If `write.csv` writes the first column header as an empty string, that matches the original file format.)

Upload it with the example metadata, tissue `kidney`, Submit. Expect the "New Training Validation" tab to appear, the retraining progress message, then within ~10 min: validation heatmap and PR plot, classification heatmap, NIS plots. Record wall time. If the page shows an out-of-memory error, stop and report — this is spec risk 3.

- [ ] **Step 6: Browser check — browse mode and data URL failure**

Reload; choose "Browse reference data", tissue `heart`, Submit; expect the engineered-reference heatmap. Then temporarily rename `site/data/lung_engineeredRef_ranks.rds`, browse `lung`, and confirm an error mentioning `Could not download reference data from http://localhost:8766/PACNet/data/lung_engineeredRef_ranks.rds` is shown (alert or Shiny error) rather than a silent hang; rename it back.

- [ ] **Step 7: Commit**

```bash
git add tools/build_site.sh .gitignore
git commit -m "Add static site build (shinylive export + reference data)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: CI workflows and README

**Files:**
- Create: `.github/workflows/pages.yml`
- Create: `.github/workflows/tests.yml`
- Modify: `README.md` ("See also" paragraph + new "Web app" section)

- [ ] **Step 1: Pages workflow**

`.github/workflows/pages.yml`:

```yaml
name: Deploy web app to GitHub Pages

on:
  push:
    branches: [main]
    paths: ["webApp/**", "data/**", "tools/build_site.sh", ".github/workflows/pages.yml"]
  workflow_dispatch:

permissions:
  contents: read
  pages: write
  id-token: write

concurrency:
  group: pages
  cancel-in-progress: true

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
        with:
          use-public-rspm: true
      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          packages: |
            any::shinylive
            any::shiny
            any::shinythemes
            any::shinycssloaders
            any::shinyWidgets
            any::shinyjs
            any::shinyalert
            any::shinybusy
            any::ggplot2
            any::gridExtra
            any::RColorBrewer
            any::pheatmap
            any::plotly
            any::dplyr
            any::randomForest
      - run: tools/build_site.sh site
      - uses: actions/upload-pages-artifact@v3
        with:
          path: site

  deploy:
    needs: build
    runs-on: ubuntu-latest
    environment:
      name: github-pages
      url: ${{ steps.deployment.outputs.page_url }}
    steps:
      - id: deployment
        uses: actions/deploy-pages@v4
```

- [ ] **Step 2: Tests workflow**

`.github/workflows/tests.yml`:

```yaml
name: Web app tests

on:
  push:
    paths: ["webApp/**", "data/**", "tests/**", ".github/workflows/tests.yml"]
  pull_request:
  workflow_dispatch:

jobs:
  test:
    runs-on: ubuntu-latest
    env:
      PACNET_SLOW: "1"
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
        with:
          use-public-rspm: true
      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          packages: |
            any::testthat
            any::withr
            any::codetools
            any::shiny
            any::shinythemes
            any::shinycssloaders
            any::shinyWidgets
            any::shinyjs
            any::shinyalert
            any::shinybusy
            any::ggplot2
            any::gridExtra
            any::RColorBrewer
            any::pheatmap
            any::plotly
            any::dplyr
            any::randomForest
      - run: Rscript -e 'testthat::test_dir("tests/testthat", stop_on_failure = TRUE)'
```

- [ ] **Step 3: Validate workflow YAML locally**

Run: `python3 -c "import yaml,sys; [yaml.safe_load(open(f)) for f in sys.argv[1:]]; print('ok')" .github/workflows/pages.yml .github/workflows/tests.yml`
Expected: `ok`

- [ ] **Step 4: README**

In `README.md`, replace the paragraph starting `PACNet is also available as a [web application](http://cahanlab.org/resources/agnosticCellNet_web/)` with:

```markdown
PACNet is also available as a [web application](https://pcahan1.github.io/PACNet/), which takes as input an expression matrix (counts, TPM, or FPKM) and sample meta-data. The application performs CellNet analysis. Additionally, this tool includes analysis of many state-of-the-art differentiation protocols so that you can benchmark your results against those commonly used methods. The web app runs entirely in your browser (nothing is uploaded to a server); the first visit downloads the R runtime and packages, which takes a minute or so. See [Running the web app locally](#webapp-local) to run it in your own R session instead.
```

Add a section before `### Data <a name="data"></a>`:

```markdown
### Running the web app locally <a name="webapp-local"></a>

With R (≥ 4.4) installed:

```r
install.packages(c("shiny", "shinythemes", "shinycssloaders", "shinyWidgets", "shinyjs",
                   "shinyalert", "shinybusy", "ggplot2", "gridExtra", "RColorBrewer",
                   "pheatmap", "plotly", "dplyr", "randomForest"))
shiny::runGitHub("PACNet", "pcahan1", subdir = "webApp")
```

The reference data the app uses is in `data/` (built from the original objects by `tools/make_slim_data.R`). Developers: `Rscript -e 'testthat::test_dir("tests/testthat")'` runs the regression tests (`PACNET_SLOW=1` includes the retraining test); `tools/build_site.sh` builds the static site that `.github/workflows/pages.yml` deploys.
```

(Use a four-backtick fence in the actual README edit if needed so the inner R block renders.)

- [ ] **Step 5: Commit**

```bash
git add .github/workflows/pages.yml .github/workflows/tests.yml README.md
git commit -m "CI: deploy web app to GitHub Pages and run regression tests; README

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 6: Hand-off note for the repo owner**

Not automatable from here: after merging to `main`, enable Settings → Pages → Build and deployment → Source: "GitHub Actions". Report this in the final summary.
