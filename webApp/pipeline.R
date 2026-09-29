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
