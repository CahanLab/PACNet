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

  # tVals are named lists (one value per gene)
  meanVals <- unlist(nis$tVals$mean[netGenes])
  sdVals <- unlist(nis$tVals$sd[netGenes])

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
  meanVect<-meanVals;
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
