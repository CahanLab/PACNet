# Utilities for PACNet webApp only

subsetGRNall <- function(grnAll, iGenes) {
   # Subset grnTable based on iGenes
   allTargets <- grnAll$overallGRN$grnTable$TG
   newGRNTable <- grnAll$overallGRN$grnTable[which(allTargets %in% iGenes),]
   newTFsAll <- newGRNTable$TF
   newGRNTable <- newGRNTable[which(newTFsAll %in% iGenes),]
   grnAll$overallGRN$grnTable <- newGRNTable
   
   ## Subset overallGRN graph based on iGenes
   vertex_names <- V(grnAll$overallGRN$graph)$name
   graph_iGenes <- which(vertex_names %in% iGenes)
   newGraph <- induced_subgraph(graph=grnAll$overallGRN$graph, vids=graph_iGenes, impl="copy_and_delete")
   grnAll$overallGRN$graph <- newGraph
   
   # Subset specGenes based on iGenes and tissue type
   tissueTypes <- names(grnAll$specGenes$context$general)
   newGeneral <- grnAll$specGenes$context$general
   for (tissueType in tissueTypes) {
      tissueSpecGenes <- newGeneral[[tissueType]]
      tissueSpecGenes <- tissueSpecGenes[which(names(tissueSpecGenes) %in% iGenes)]
      newGeneral[[tissueType]] <- tissueSpecGenes
   }
   grnAll$specGenes$context$general <- newGeneral
   
   # Subset ctGRN geneLists, graphLists, and tfTargets  based on iGenes and tissue type
   grnAll$ctGRNs$geneLists <- newGeneral
   
   newGraphLists <- grnAll$ctGRNs$graphLists
   for (tissueType in tissueTypes) {
      tissueGRN <- newGraphLists[[tissueType]]
      iVertices <- vertex_attr(tissueGRN, name="name")
      iVertices <- iVertices[which(iVertices %in% iGenes)]
      tissueGRN <- induced_subgraph(graph=tissueGRN, vids=iVertices, impl="copy_and_delete")
      newGraphLists[[tissueType]] <- tissueGRN
   }
   grnAll$ctGRNs$graphLists <- newGraphLists
   
   newTFTargets <- grnAll$ctGRNs$tfTargets
   for (tissueType in tissueTypes) {
      tissueTFTargets <- newTFTargets[[tissueType]]
      tissueTFTargets <- tissueTFTargets[which(names(tissueTFTargets) %in% iGenes)]
      for (TF in names(tissueTFTargets)) {
         newTargets <- tissueTFTargets[[TF]]
         newTargets <- newTargets[which(newTargets %in% iGenes)]
         tissueTFTargets[[TF]] <- newTargets
      }
      newTFTargets[[tissueType]] <- tissueTFTargets
   }
   grnAll$ctGRNs$tfTargets <- newTFTargets
   
   return(grnAll)
}

subsetTrainNormParam <- function(trainNormParam, grnAll, iGenes) {
   tissueTypes <- names(grnAll$specGenes$context$general)
   newTVals <- trainNormParam$tVals
   for (tissueType in tissueTypes) {
      newIndices <- which(names(newTVals[[tissueType]][["mean"]]) %in% iGenes)
      newTVals[[tissueType]][["mean"]] <- newTVals[[tissueType]][["mean"]][newIndices]
      newTVals[[tissueType]][["sd"]] <- newTVals[[tissueType]][["sd"]][newIndices]
   }
   trainNormParam$tVals <- newTVals
   return(trainNormParam)
}


### compute zscore
zscore<-function(x,### numeric vector
 meanVal, ### mean of distribution to compute zscore of x against 
 sdVal ### standard deviation of distribution to compute zscore of x agains
){ 
  (x-meanVal)/sdVal;
  ### zscore
}


### Compute the mean zscore of given genes in each sample
cn_zscoreVect<-function(genes, ### genes
 xvals, ### named vector
 tVals, ### tvals
 ctt ### ctt
){
  ans<-vector();
  for(gene in genes){
    ans<-append(ans, zscore(xvals[gene], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]]));
  }
  ans;
  ### zscore vector
}



#' network influence score
#'
#' Computes network influence score (NIS).
#' @param expDat ranked query expression matrix, output of running logRank()
#' @param stQuery query sample metadata table
#' @param iGenes character vector of gene names found in both training data and query data 
#' @param grnAll GRN object, output of running ccn_makeGRN()
#' @param trainNormParam normalization parameters, output of running ccn_trainNorm()
#' @param subnet name of cell/tissue type subnetwork to evaluate
#' @param ctt string indicating the cell/tissue type to compare against
#' @param colname_sid colname of stQuery with unique sample labels
#' @param relaWeight whether to weight by overall expression such that TFs with higher expression in ctt are more important (1=do the weighting)
#'
#' @return numeric matrix where rows are TFs in CT GRN and columns are query samples. Negative score indicates TF should be upregulated to acquire a ctt fate. Positive score indicates TF should be downregulated to acquire a ctt fate.
pacnet_nis <- function(expDat, stQuery, iGenes, grnAll, trainNormParam, subnet, ctt, colname_sid="sra_id", relaWeight=0) {
  tfTargList<-grnAll$ctGRNs$tfTargets
  nTargets<-vector();
  targetScore<-vector();
  tfScore<-vector();
  totalScore<-vector();
  tfWeights<-vector();
  
  tfs<-names(tfTargList[[subnet]]);
  netGenes<-names(grnAll$ctGRNs$geneLists[[subnet]])
  netGenes<-intersect(netGenes, iGenes)
  
  sids<-as.vector(stQuery[,colname_sid])
  
  ans<-matrix(0, nrow=length(tfs), ncol=nrow(stQuery))
  rownames(ans)<-tfs;
  colnames(ans)<-sids;
  
  tVals <- trainNormParam[["tVals"]]
  
  # compute a matrix of zscores.
  zzzMat<-matrix(0, nrow=length(netGenes), ncol=nrow(stQuery));
  
  for(i in seq(length(sids))){
    sid<-sids[i];
    print(paste0("Computing zscores for ", sid))
    xvals<-as.vector(expDat[netGenes, sid])
    names(xvals)<-netGenes
    zzzMat[,i]<-cn_zscoreVect(netGenes, xvals, tVals, ctt)
  }
  
  rownames(zzzMat)<-netGenes
  colnames(zzzMat)<-stQuery[,colname_sid]
  
  for(sid in sids) {
    print(paste0("TF scoring for ", sid));
    xvals<-as.vector(expDat[,sid]);
    names(xvals)<-rownames(expDat);
    
    # assign weights
    
    ### # meanVect<-unlist(tVals[[ctt]][['mean']][netGenes]);
    meanVect<-unlist(tVals[[subnet]][['mean']][netGenes]);
    meanVect <- meanVect / min(meanVect) # Added to make function work with rank-based tVals
    weights<-(2**meanVect)/sum(2**meanVect);
    
    for(i in seq(length(tfs))){
      tf<-tfs[i];
      # zscore of TF relative to target C/T
      ##      tfScore[i]<-zscore(xvals[tf], tVals[[ctt]][['mean']][[tf]], tVals[[ctt]][['sd']][[tf]]);
      
      tfScore[i]<-zzzMat[tf,sid];
      
      targs<-tfTargList[[subnet]][[tf]];
      targs<-intersect(targs, iGenes);
      
      # Zscores of TF targets, relative to C/T
      ##      tmp<-cn_zscoreVect(targs, xvals, tVals, ctt );
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

