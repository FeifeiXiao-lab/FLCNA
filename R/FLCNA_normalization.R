#' @title FLCNA normalization
#' 
#' @description Normalization function used in FLCNA.
#' 
#' @param Y A p-dimensional data matrix. Each row is an observation.
#' @param bin_size The bin size used in the data the default is 100,000.
#' @param gc A p-dimensional vector with gc concent.
#' @param map A p-dimensional vector with mappability.
#'
#'
#' @return The log2Rdata used for main step for the FLCNA method.
#' 
#' @export
FLCNA_normalization <- function(Y, bin_size=100000, gc, map){
  
  RCTL <- Y
  RC_norm=matrix(data=NA,nrow = dim(RCTL)[1],ncol = dim(RCTL)[2])
  rownames(RC_norm)=rownames(RCTL)
  colnames(RC_norm)=colnames(RCTL)
  for(sub in 1:dim(RCTL)[2]){
    ### Bin size normalization only if InTarget###
    step <- 5
    RCLNormListIn <- CorrectSize(RCTL[,sub],L=bin_size,step)
    RCLNormIn <- RCLNormListIn$RCNorm
    
    ### Mappability normalization ###
    step <- 0.01
    RCMAPNormListIn <- CorrectMAP(RCLNormIn,MAPContent=map,step)
    RCMAPNormIn <- RCMAPNormListIn$RCNorm
    
    ### GC-content Normalization ###
    step <- 5
    RCGCNormListIn <- CorrectGC(RCMAPNormIn,GCContent=gc,step)
    RCGCNormIn <- RCGCNormListIn$RCNorm
    
    RC_norm[,sub]=RCGCNormIn
  }
  
  RC_norm <- RC_norm + matrix(rgamma(nrow(RC_norm)*ncol(RC_norm),0.1,1),nrow(RC_norm),ncol(RC_norm))
  
  lrr=matrix(data=NA,nrow = dim(RC_norm)[1],ncol = dim(RC_norm)[2])
  rownames(lrr)=rownames(RC_norm)
  colnames(lrr)=colnames(RC_norm)
  for(i in 1:dim(RC_norm)[1]){
    mean=mean(RC_norm[i,])
    rc=RC_norm[i,]
    lrr[i,]=log2((rc+1)/mean)
  }
  log2Rdata=lrr
  
  return(log2Rdata)
}



#' @title FLCNA normalization with reference
#' 
#' @description Normalization function used in FLCNA.
#' 
#' @param Y A p-dimensional data matrix. Each row is an observation.
#' @param bin_size The bin size used in the data the default is 100,000.
#' @param gc A p-dimensional vector with gc concent.
#' @param map A p-dimensional vector with mappability.
#' @param ref_id cells used as reference.
#'
#'
#' @return The log2Rdata used for main step for the FLCNA method.
#' 
#' @export
FLCNA_normalization_ref <- function(Y, bin_size=100000, gc, map, ref_id){
  
  RCTL <- Y
  RC_norm=matrix(data=NA,nrow = dim(RCTL)[1],ncol = dim(RCTL)[2])
  rownames(RC_norm)=rownames(RCTL)
  colnames(RC_norm)=colnames(RCTL)
  for(sub in 1:dim(RCTL)[2]){
    ### Bin size normalization only if InTarget###
    step <- 5
    RCLNormListIn <- CorrectSize(RCTL[,sub],L=bin_size,step)
    RCLNormIn <- RCLNormListIn$RCNorm
    
    ### Mappability normalization ###
    step <- 0.01
    RCMAPNormListIn <- CorrectMAP(RCLNormIn,MAPContent=map,step)
    RCMAPNormIn <- RCMAPNormListIn$RCNorm
    
    ### GC-content Normalization ###
    step <- 5
    RCGCNormListIn <- CorrectGC(RCMAPNormIn,GCContent=gc,step)
    RCGCNormIn <- RCGCNormListIn$RCNorm
    
    RC_norm[,sub]=RCGCNormIn
  }
  
  RC_norm <- RC_norm + matrix(rgamma(nrow(RC_norm)*ncol(RC_norm),0.1,1),nrow(RC_norm),ncol(RC_norm))
  RC_norm1 <- RC_norm[, -ref_id]
  ref_norm <- RC_norm[, ref_id]
  
  lrr=matrix(data=NA,nrow = dim(RC_norm1)[1],ncol = dim(RC_norm1)[2])
  rownames(lrr)=rownames(RC_norm1)
  colnames(lrr)=colnames(RC_norm1)
  for(i in 1:dim(RC_norm1)[1]){
    mean=mean(ref_norm[i,])
    rc=RC_norm1[i,]
    lrr[i,]=log2((rc+1)/mean)
  }
  log2Rdata=lrr
  
  return(log2Rdata)
}


#' @title GC content correction
#' 
#' @description Normalization according to GC content.
#' 
#' @param RC A p-dimensional data vector. The read counts for a sample.
#' @param GCContent A p-dimensional vector with gc concent.
#' @param step Step value, a constant, used in the GC correlation procedure.
#'
#'
#' @return The output for GC content correlation.
#' 
CorrectGC<-function(RC,GCContent,step){
  stepseq<-seq(min(GCContent),max(GCContent),by=step)
  #stepseq<-seq(0,100,by=step)
  MasterMedian<-median(RC,na.rm=T)
  MedianGC<-rep(0,length(stepseq)-1)
  RCNormMedian<-RC
  for (i in 1:(length(stepseq)-1)){
    if (i==1){
      ind<-which(GCContent>=stepseq[i] & GCContent<=stepseq[i+1])
    }
    if (i!=1){
      ind<-which(GCContent>stepseq[i] & GCContent<=stepseq[i+1])
    }
    if (length(ind)>0){
      m<-median(RC[ind],na.rm=T)
      if (m>0){
        MedianGC[i]<-m
        RCNormMedian[ind]<-RC[ind]*MasterMedian/m
      }
    }
  }
  RCNormList<-list()
  RCNormList$Median<-MedianGC
  RCNormList$StepGC<-stepseq[1:(length(stepseq)-1)]
  RCNormList$RCNorm<-RCNormMedian
  RCNormList
}




#' @title bin size correction
#' 
#' @description Normalization according to bin size. Since bin size is consistant in all the markers, bin size will not affect normalization results in this study. 
#' 
#' @param RC A P-dimensional data vector. The read counts for a sample.
#' @param L The bin size used in the data the default is 100,000.
#' @param step Step value, a constant, used in the bin size correlation procedure.
#'
#'
#' @return The output bin size correlation.
CorrectSize<-function(RC,L,step){
  stepseq<-seq(min(L),max(L),by=step)
  #stepseq<-seq(0,max(L),by=step)
  MasterMedian<-median(RC,na.rm=T)
  MedianL<-rep(0,length(stepseq)-1)
  RCNormMedian<-RC
  for (i in 1:(length(stepseq)-1)){
    if (i==1){
      ind<-which(L>=stepseq[i] & L<=stepseq[i+1])
    }
    if (i!=1){
      ind<-which(L>stepseq[i] & L<=stepseq[i+1])
    }
    if (length(ind)>0){
      m<-median(RC[ind],na.rm=T)
      if (m>0){
        MedianL[i]<-m
        RCNormMedian[ind]<-RC[ind]*MasterMedian/m
      }
    }
  }
  RCNormList<-list()
  RCNormList$Median<-MedianL
  RCNormList$RCNorm<-RCNormMedian
  RCNormList
}




#' @title Mapp correction
#' 
#' @description Normalization according to mappability.
#' 
#' @param RC A P-dimensional data vector. The read counts for a sample.
#' @param MAPContent A p-dimensional vector with mappability.
#' @param step Step value, a constant, used in the Mapp correlation procedure.
#'
#'
#' @return The output Mapp correlation.
CorrectMAP<-function(RC,MAPContent,step){
  stepseq<-seq(min(MAPContent),max(MAPContent),by=step)
  #stepseq<-seq(0,100,by=step)
  MasterMedian<-median(RC,na.rm=T)
  MedianMAP<-rep(0,length(stepseq)-1)
  RCNormMedian<-RC
  for (i in 1:(length(stepseq)-1)){
    if (i==1){
      ind<-which(MAPContent>=stepseq[i] & MAPContent<=stepseq[i+1])
    }
    if (i!=1){
      ind<-which(MAPContent>stepseq[i] & MAPContent<=stepseq[i+1])
    }
    
    if (length(ind)>0){
      m<-median(RC[ind],na.rm=T)
      if (m>0){
        MedianMAP[i]<-m
        RCNormMedian[ind]<-RC[ind]*MasterMedian/m
      }
    }
  }
  RCNormList<-list()
  RCNormList$Median<-MedianMAP
  RCNormList$StepMAP<-stepseq[1:(length(stepseq)-1)]
  RCNormList$RCNorm<-RCNormMedian
  RCNormList
}

