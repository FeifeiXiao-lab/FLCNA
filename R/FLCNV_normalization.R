FLCNV_normalization <- function(Y, norm.type="mean", bin_size, gc, map){
  
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
  for(i in 1:dim(RC_norm)[2]){
    mean=mean(RC_norm[,i])
    rc=RC_norm[,i]
    lrr[,i]=log2((rc+1)/mean)
  }
  log2Rdata=lrr
  
  return(log2Rdata)
}



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

################ Correzione dei RC Per il Bin Size ####################
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

################ Correzione dei RC dalla Mappability ####################
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

