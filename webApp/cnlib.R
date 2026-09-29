# Vendored from https://github.com/CahanLab/cancerCellNet (MIT License),
# commit dc29a8e7982461416a12fd3372bb844d35bd0ee0, by tools/vendor_cnlib.R.
# Only the functions the PACNet web app calls are included. Local changes
# are marked with "# PACNet-web:" comments.

addRandToSampTab<-function(classRes, sampTab, desc, id) {
  cNames = colnames(classRes)
  snames = as.vector(sampTab[, id])

  rnames = setdiff(cNames, snames)

  cat("number of random samples: ",length(rnames), "\n")

  stNew = data.frame(rid=rnames, rdesc=rep("rand", length(rnames)))
  stTop = sampTab[,c(id, desc)]
  colnames(stNew) = c(id, desc)

  ans = rbind(stTop, stNew)
  rownames(ans) = colnames(classRes)

  return(ans)
}

broadClass_predict<-function(cnProc, expDat, nrand = 2) {
   ccnList = cnProc
   expVal = expDat

   rf_tsp = ccnList[['classifier']]
   cgenes = ccnList[['cgenes']]
   xpairs = ccnList[['xpairs']]

   cat("Loaded in the cnProc\n")

   expValTrans = query_transform(expVal, xpairs)
   classRes_val = rf_classPredict(rf_tsp, expValTrans, numRand=nrand)

   cat("All Done\n")

   return(classRes_val)
}

broadClass_train<-function(stTrain, expTrain, colName_cat, colName_samp="row.names", nTopGenes = 20, nTopGenePairs = 50, nRand = 40, nTrees = 1000, stratify=FALSE, sampsize=40, weightedDown_total = 5e5, weightedDown_dThresh = 0.25, transprop_xFact = 1e5, quickPairs = FALSE, coreProportion = 0) {

   if (class(stTrain) != "data.frame") {
      stTrain = as.data.frame(stTrain)
   }

   if (colName_samp != "row.names") {
     rownames(stTrain) = stTrain[, colName_samp]
   }

   cat("Sample table has been prepared\n")

   expTnorm = trans_prop(weighted_down(expTrain, weightedDown_total, dThresh = weightedDown_dThresh), transprop_xFact)
   cat("Expression data has been normalized\n")

   system.time(cgenes<-findClassyGenes(expDat = expTnorm, sampTab = stTrain, dLevel = colName_cat, topX = nTopGenes, sliceSize=1000, coreProportion = coreProportion))
   cat("Finished finding classification genes\n")

   cgenesA = cgenes[['cgenes']]
   grps = cgenes[['grps']]
   cgenes_list = cgenes[['labelled_cgenes']]

   cat("There are ", length(cgenesA), " classification genes\n")

   system.time(xpairs_list<-ptGetTop(expTrain[cgenesA,], grps, cgenes_list, topX=nTopGenePairs, sliceSize=2000, quickPairs=quickPairs, coreProportion = coreProportion))
   cat("Finished finding top gene pairs\n")

   # compile the genepair list
   xpairs = c()
   for(item in xpairs_list) {
      xpairs = c(xpairs, item)
   }

   xpairs = names(xpairs)
   xpairs = unique(xpairs)

   # some of these might include selection cassettes; remove them
   xi = setdiff(1:length(xpairs), grep("selection", xpairs))
   xpairs = xpairs[xi]

   system.time(pdTrain<-query_transform(expTrain[cgenesA, ], xpairs))
   cat("Finished pair transforming the data\n")

   tspRF = makeClassifier(pdTrain[xpairs,], genes=xpairs, groups=grps, nRand = nRand, ntrees = nTrees, stratify=stratify, sampsize=sampsize)
   cnProc = list("cgenes"= cgenesA, "xpairs"=xpairs, "grps"= grps, "classifier" = tspRF)

   returnList = list("sampTab" = stTrain, "cgenes_list" = cgenes_list, "cnProc" = cnProc, "xpairs_list" = xpairs_list)

   cat("All Done\n")
   return(returnList)
}

ccn_classAssess<-function(ct_scores, stVal, classLevels="description2", dLevelSID="sample_id", resolution=0.005) {
  allPRs = list()
  evalAll = matrix(0, nrow=nrow(ct_scores),ncol=2)
  classifications = rownames(ct_scores)
  rownames(stVal) = as.vector(stVal[,dLevelSID])
  i = 1

  for(xname in classifications){
    classification = classifications[i]
    tmpPR = cn_eval(ct_scores[xname,],
                     stVal,
                     classLevels,
                     xname,threshs=seq(0,1, by=resolution), dLevelSID=dLevelSID)
    allPRs[[xname]] = as.data.frame(tmpPR)
    i = 1+i
  }

  return(allPRs)
}

cn_clPerf<-function(vect, sampTab, dLevel, classification, thresh, dLevelSID="sample_id"){
  TP = 0
  FN = 0
  FP = 0
  TN = 0
  sampIDs = names(vect)
  classes = as.vector(sampTab[sampIDs,dLevel])

  actualPos = as.vector(sampTab[sampTab[,dLevel]==classification,dLevelSID])
  actualNeg = setdiff(sampIDs, actualPos)

  calledPos = names(which(vect>thresh))
  calledNeg = names(which(vect<=thresh))

  TP = length(intersect(actualPos, calledPos))
  FP = length(intersect(actualNeg, calledPos))
  FN = length(actualPos)-TP
  TN = length(actualNeg)-FP

  return(c(TP, FN, FP, TN))
}

cn_eval<-function(vect, sampTab, dLevel, classification, threshs=seq(0,1,by=0.05), dLevelSID="sample_id"){
  ans = matrix(0,nrow=length(threshs), ncol=9)
  for(i in seq(length(threshs))){
    thresh = threshs[i]
    ans[i,1:4] = cn_clPerf(vect, sampTab, dLevel, classification, thresh, dLevelSID=dLevelSID)
  }
  ans[,5] = threshs
  colnames(ans) = c("TP", "FN", "FP", "TN", "thresh","FPR", "TPR", "Sens", "Prec")
  TPR = ans[,'TP']/(ans[,'TP']+ans[,'FN'])
  FPR = ans[,'FP']/(ans[,'TN']+ans[,'FP'])
  sens = ans[,"TP"]/(ans[,"TP"]+ans[,"FN"])
  prec = ans[,"TP"]/(ans[,"TP"]+ans[,"FP"])
  ans[,'TPR'] = TPR
  ans[,'FPR'] = FPR
  ans[,'Sens'] = sens
  ans[,'Prec'] = prec

  return(ans)
}

downSampleW<-function(vector, total=1e5, dThresh=0){

  totalSignal = sum(vector) # get the sum of the vector
  wAve = vector/totalSignal 

  resid = totalSignal-total #num to subtract from sample
  residW = wAve*resid 

  ans = vector-residW
  ans[which(ans<dThresh)] = 0

  return(ans)
}

findBestPairs <- function(xdiff, n=50,maxPer=3) {

  # error catching in case the number of pairs wanted is more than pairs generated
  if(nrow(xdiff) < n) {
    cat("there are only", nrow(xdiff), "genepairs generated.", "\n")

    ans = as.vector(xdiff$cval)
    names(ans) = rownames(xdiff)
    return(ans)

  }
  else {
    xdiff = xdiff[order(abs(xdiff$cval), decreasing=TRUE),]

    genes = unique(unlist(strsplit(rownames(xdiff), "_")))
    countList = rep(0, length(genes))
    names(countList) = genes

    i = 0
    ans_names = vector()
    ans_signs = vector()

    xdiff_index = 1
    pair_names = rownames(xdiff)

    backup_vector = c()
    backup_vector_sign = c()

    while(i < n ){
      tmpAns = pair_names[xdiff_index]
      tmpSigns = sign(as.numeric(xdiff[tmpAns, "cval"]))

      tgp = unlist(strsplit(tmpAns, "_"))

      if((countList[tgp[1]] < maxPer) & (countList[tgp[2]] < maxPer )){

        # record down the gene pair name and sign
        ans_names = append(ans_names, tmpAns)
        ans_signs = append(ans_signs, tmpSigns)

        countList[ tgp[1] ] = countList[ tgp[1] ]+ 1
        countList[ tgp[2] ] = countList[ tgp[2] ]+ 1

        i = i+1
      }

      else {
        backup_vector = c(backup_vector, tmpAns) # place into backup vector
        backup_vector_sign = c(backup_vector_sign, tmpSigns)
      }


      xdiff_index = xdiff_index + 1

      # in the case where the original list is exhausted, dig into the backup vector
      if(xdiff_index > length(pair_names)) {
        additional_pairs = backup_vector[1:(n - i)]
        additional_pairs_sign = backup_vector_sign[1:(n - i)]

        ans_names = c(ans_names, additional_pairs)
        ans_signs = c(ans_signs, additional_pairs_sign)
        i = length(ans_signs)
      }

    }

    # assign the signs and names to the return answer
    ans = ans_signs
    names(ans) = ans_names

    ans = na.omit(ans) # just in case there were NA
    return(ans)
  }

}

findClassyGenes <- function(expDat, sampTab, dLevel, topX=25, dThresh=0, alpha1=0.05, alpha2=.001, mu=2, sliceSize = 1000, coreProportion = 1/4) {
  if((dLevel %in% colnames(sampTab)) == FALSE) {
    stop("Please enter the correct column name for sampTab that indicates the categories")
  }

  if((topX * 3 > nrow(expDat)) == TRUE) {
    stop(paste0("Please enter a topX value smaller than ", as.integer(nrow(expDat) / 3)))
  }

  # remove duplicated genes
  expDat = expDat[!duplicated(rownames(expDat)), ]

  gsTrain = sc_statTab(expDat, dThresh=dThresh)
  ggenes = sc_filterGenes(gsTrain, alpha1=alpha1, alpha2=alpha2, mu=mu)
  grps = as.vector(sampTab[,dLevel])
  names(grps) = rownames(sampTab)

  # more robust way to avoid subscript out of bound
  if(length(ggenes) != nrow(expDat)) {
    expDat = expDat[ggenes,]
  }

  # PACNet-web: serial only. The core count is NA under webR and the browser
  # cannot start clusters; upstream defaults to serial (coreProportion = 0).
  xdiff = gnrAll(expDat, grps, sliceSize, coreProportion)
  cgenes = lapply(X = xdiff, FUN = getClassGenes, topX = topX)

  labelled_cgenes = cgenes
  cgenes = unique(unlist(cgenes))

  return(list(cgenes=cgenes, grps=grps, labelled_cgenes=labelled_cgenes))
}

getClassGenes<-function(diffRes, topX=25, bottom=TRUE, leastDiff = TRUE) {
  #exclude NAs
  xi = which(!is.na(diffRes$cval))
  diffRes = diffRes[xi,] # exclude the NAs. Select the rows that does not have NAs

  diffRes = diffRes[order(diffRes$cval, decreasing=TRUE),] #order based on classification value largest to lowest
  ans = rownames(diffRes[1:topX,]) # get the top 20 genes

  if(bottom){
    ans = append(ans, rownames( diffRes[nrow(diffRes) - ((topX-1):0),]))
  }

  # get the least differentially expressed genes as house holders
  if(leastDiff) {
    sameRes = diffRes[order(abs(diffRes$cval), decreasing = FALSE), ]
    ans = append(ans, rownames(sameRes[1:topX, ]))
  }

  return(ans)
}

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

logRank <- function(expTrain, base = 0) {
  expTrain = apply(expTrain, FUN = rank, MARGIN = 2)

  if (base == 0) {
    return(expTrain)
  }
  else {
    return(log(expTrain, base = base))

  }
}

makeClassifier<-function(expTrain, genes, groups, nRand=50, ntrees=2000, stratify=FALSE, sampsize=40){
  randDat = randomize(expTrain, num=nRand)
  #randDat = ModifiedRandomize(expTrain, num=nRand)

  expTrain = cbind(expTrain, randDat)

  allgenes = rownames(expTrain)
  missingGenes = setdiff(unique(genes), allgenes)
  cat("Number of missing genes ", length(missingGenes),"\n")
  ggenes = intersect(unique(genes), allgenes)

  # return random forest object
  if(!stratify){
    return(randomForest::randomForest(t(expTrain[ggenes,]), as.factor(c(groups, rep("rand", ncol(randDat)))), ntree=ntrees))
  }else{
    return(randomForest::randomForest(t(expTrain[ggenes,]), as.factor(c(groups, rep("rand", ncol(randDat)))), ntree=ntrees, strata = as.factor(c(groups, rep("rand", ncol(randDat)))), sampsize=rep(sampsize, length(c(unique(groups), "rand")))))
  }
}

makePairTab <- function(genes) {
  pTab = t(combn(genes, 2))
  colnames(pTab) = c("genes1", "genes2")
  pTab = cbind(pTab, pairName=paste(pTab[,1], "_",pTab[,2], sep=''))

  return(pTab)
}

plot_class_PRs<-function(assessed){
  ctts = names(assessed)
  df = data.frame()
  for(ctt in ctts){
    # to check if the category has a PR representation in validation
    if(max(assessed[[ctt]][,'TP']) == 0){
      next
    }

    tmp = assessed[[ctt]]
    tmp = cbind(tmp, ctype=ctt)
    df = rbind(df, tmp)
  }

  return(ggplot2::ggplot(data=df, aes(x=Sens, y=Prec)) + geom_point(size = .5, alpha=.5) +  geom_path(size=.5, alpha=.75) +
    theme_bw() + xlab("Recall") + ylab("Precision") + facet_wrap( ~ ctype, ncol=4) +
    theme(axis.text = element_text(size=5)) + ggtitle("Classifier performance"))
}

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

      # PACNet-web: score the pairs in slices (see sc_testPattern_sliced). Scoring
      # ~45k pairs at once needs ~235 MB per class and exhausted browser memory.
      tmpAns = sc_testPattern_sliced(myPatternG[[cct]], expDat, pairTab, sliceSize)

      ans[[cct]] = findBestPairs(tmpAns, topX)
    }

    return(ans)
  }
}

ptSmall <- function(expDat, pTab) {
  npairs = nrow(pTab)
  ans = matrix(0, nrow=npairs, ncol=ncol(expDat))
  genes1 = as.vector(pTab[, "genes1"])
  genes2 = as.vector(pTab[, "genes2"])

  for(i in seq(nrow(pTab))) {
    ans[i,] = as.numeric(expDat[genes1[i],]>expDat[genes2[i],])
  }
  colnames(ans) = colnames(expDat)
  rownames(ans) = as.vector(pTab[, "pairName"])

  return(ans)
}

query_transform <- function(expDat, genePairs) {

  expDat = as.matrix(expDat) # additional, in case that expDat is in dataframe which is very likely
  genes = strsplit(genePairs, "_") # split the gene pairs into seperate genes
  ans = matrix(0, nrow=length(genes), ncol=ncol(expDat))

  for(i in seq(length(genes))){
    ans[i,] = as.numeric(expDat[genes[[i]][1],]>expDat[genes[[i]][2],])
  }
  colnames(ans) = colnames(expDat)
  rownames(ans) = genePairs

  return(ans)
}

randomize<-function(expDat, num=50) {

  randDat = t(apply(expDat, 1, sample))
  randDat = apply(randDat, 2, sample)

  randDat = randDat[,sample(1:ncol(randDat), num)]

  colnames(randDat) = paste0(rep("rand_", num), 1:num)
  rownames(randDat) = rownames(expDat)

  return(randDat)
}

rf_classPredict<-function(rfObj, expQuery, numRand=50) {

  randDat = randomize(expQuery, num=numRand) # generate the random profiles
  expQuery = cbind(expQuery, randDat)

  preds = rownames(rfObj$importance)
  xpreds = t(predict(rfObj, t(expQuery[preds,]), type='prob'))
  colnames(xpreds) = colnames(expQuery)

  return(xpreds)
}

sc_compAlpha<-function(expMat, threshold=0,pseudo=FALSE){
  # identify the index of vectors greater than threshold
  lengthFunction<-function(vector, threshold){
    length(which(vector>threshold))
  }

  alphas = apply(expMat, 1, lengthFunction, threshold)

  ans = alphas/ncol(expMat)

  if (pseudo) {
    ans = (alphas+1) / (ncol(expMat) + 1)
  }

  return(ans)
}

sc_compMu<-function(expMat, threshold=0){

  afunct<-function(vector, threshold){
    mean(vector[which(vector>threshold)])
  }

  mus = unlist(apply(expMat, 1, afunct, threshold))
  mus[is.na(mus)] = 0

  return(mus)
}

sc_cov<-function(vector){
  return(sd(vector)/mean(vector))
}

sc_fano<-function(vector){
  return(var(vector)/mean(vector))
}

sc_filterGenes<-function(geneStats, alpha1=0.1, alpha2=0.01, mu=2){

  return(rownames(geneStats[(geneStats$alpha > alpha1) | (geneStats$alpha > alpha2 & geneStats$mu > mu), ]))
}

sc_sampR_to_pattern<-function(sampR){
  d_ids = unique(as.vector(sampR))
  nnnc = length(sampR)

  ans = list()
  for(d_id in d_ids){
    x = rep(0,nnnc)
    x[which(sampR==d_id)] = 1
    ans[[d_id]] = x
  }

  return(ans)
}

sc_statTab<-function(expDat, dThresh=0){
  statTab = data.frame()
  # below generates various statistical values
  muAll = sc_compMu(expDat, threshold=dThresh)
  alphaAll = sc_compAlpha(expDat,threshold=dThresh)
  meanAll = apply(expDat, 1, mean)
  covAll = apply(expDat, 1, sc_cov)
  fanoAll = apply(expDat,1, sc_fano)
  maxAll = apply(expDat, 1, max)
  sdAll = apply(expDat, 1, sd)

  statTabAll = data.frame(gene=rownames(expDat), mu=muAll, alpha=alphaAll, overall_mean=meanAll, cov=covAll, fano=fanoAll, max_val=maxAll, sd=sdAll)

  return(statTabAll)
}

sc_testPattern<-function(pattern, expDat){
  pval = vector()
  cval = vector()

  geneids = rownames(expDat)  # get the gene names of the exp matrix

  llfit = ls.print(lsfit(pattern, t(expDat)), digits=25, print=FALSE)

  pos_R_index = which(as.numeric(llfit$summary[,2]) >= 0) # R-squared value

  xxx = matrix(unlist(llfit$coef), ncol=8,byrow=TRUE)

  ccorr = xxx[pos_R_index,6]  # t-value of X coefficient
  cval =  sqrt(as.numeric(llfit$summary[pos_R_index,2])) * sign(ccorr)  # R squared

  pval = as.numeric(xxx[,8])  # p-value of X coefficient
  holm = p.adjust(pval, method='holm')


  return(data.frame(row.names=geneids[pos_R_index], pval=pval[pos_R_index], cval=cval,holm=holm[pos_R_index]))
}

# PACNet-web: sc_testPattern() over the gene pairs in pairTab, one slice of
# pairs at a time, returning exactly what sc_testPattern(pattern,
# ptSmall(expDat, pairTab)) returns while holding only one slice in memory.
# ls.print() rounds R-squared with format() across the whole batch, so the
# raw R-squared values are collected and formatted once, over all pairs.
sc_testPattern_sliced<-function(pattern, expDat, pairTab, sliceSize){
  geneids = character(0)
  rsquared = numeric(0)
  coefs = list()
  nPairs = nrow(pairTab)

  for(str in seq(1, nPairs, by = sliceSize)){
    stp = min(str + sliceSize - 1, nPairs)
    pDat = ptSmall(expDat, pairTab[str:stp, , drop = FALSE])
    fit = lsfit(pattern, t(pDat))
    llfit = ls.print(fit, digits=25, print.it=FALSE)

    # R-squared exactly as ls.print() computes it, before formatting
    lsqr = fit$qr
    totss = if (is.matrix(lsqr$qt)) colSums(lsqr$qt[-1L, ]^2) else sum(lsqr$qt[-1L]^2)
    resss = colSums(as.matrix(fit$residuals)^2, na.rm = TRUE)
    rsquared = c(rsquared, (totss - resss)/totss)

    coefs[[length(coefs) + 1]] = matrix(unlist(llfit$coef), ncol=8, byrow=TRUE)
    geneids = c(geneids, rownames(pDat))
  }

  rsq = as.numeric(format(round(rsquared, 25))) # as ls.print() formats a batch

  pos_R_index = which(rsq >= 0)
  xxx = do.call(rbind, coefs)
  ccorr = xxx[pos_R_index,6]  # t-value of X coefficient
  cval =  sqrt(rsq[pos_R_index]) * sign(ccorr)
  pval = as.numeric(xxx[,8])  # p-value of X coefficient
  holm = p.adjust(pval, method='holm')

  return(data.frame(row.names=geneids[pos_R_index], pval=pval[pos_R_index], cval=cval,holm=holm[pos_R_index]))
}

splitCommon_proportion <- function(sampTab, proportion, dLevel) {
  trainSampleTable = NULL
  for(myCat in unique(sampTab[, dLevel])) {

    tempSt = sampTab[sampTab[, dLevel] == myCat, ]
    tempSt = tempSt[sample(rownames(tempSt), as.integer((proportion)*nrow(tempSt))), ]

    if(is.null(trainSampleTable) == TRUE) {
      trainSampleTable = tempSt
    }
    else {
      trainSampleTable = rbind(trainSampleTable, tempSt)
    }
  }

  valSampleTable = sampTab[!(rownames(sampTab) %in% rownames(trainSampleTable)), ]

  return(list("trainingSet" = trainSampleTable, "validationSet" = valSampleTable))
}

trans_prop<-function(expDat, xFact=1e5){
  ans = matrix(0, nrow=nrow(expDat), ncol=ncol(expDat))

  for(i in seq(ncol(expDat))){
    ans[,i] = expDat[,i]/sum(expDat[,i])
  }

  ans = ans*xFact
  colnames(ans) = colnames(expDat)
  rownames(ans) = rownames(expDat)

  return(log(1+ans))
}

weighted_down<-function(expRaw, total=1e5, dThresh=0) {
  expCountDnW = apply(expRaw, 2, downSampleW, total=total, dThresh=dThresh)

  return(expCountDnW)
}

