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
