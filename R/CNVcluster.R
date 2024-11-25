#' @title CNA output
#' 
#' @description This function clusters the identified change-points to make final CNA calling. 
#' 
#' @param mean.matrix The cluster mean matrix estimated from FLCNA R function.
#' @param LRR log2R data after normlization.
#' @param Clusters Cluster index for each cell indentified from FLCNA R function.
#' @param QC_ref Reference file after QC.
#' @param cutoff Cutoff value to further control the number of CNAs, the larger value of cutoff, the smaller number of CNAs.
#' @param L Repeat times in the EM algorithm, defaults to 100.
#' 
#' 
#' @return The return is the clustered CNA segments with start position, end position and copy number states.
#' \item{state}{The CNA states assigned.}
#' \item{start}{The start point of CNAs.}
#' \item{end}{The end point of CNAs.}
#' \item{chr}{Chromosome of CNAs.}
#' \item{width}{The width of CNAs.}
#'
#' @export
#' 
#' 
CNA.out_pool <- function(mean.matrix, LRR, Clusters, QC_ref, cutoff=0.80, L=100){
  
  Chr_index <- as.numeric(gsub("^.{0,3}", "", rep(QC_ref@seqnames@values, QC_ref@seqnames@lengths)))
  Chr_cp_index <- 1+which(abs(diff(Chr_index))>=1)
  
  ## Initial para to clustering CNV state in GMM
  set.seed(2024)
  init<- Para_init(mean.matrix = mean.matrix, LRR=LRR, ref=QC_ref, 
                   Clusters=Clusters, cutoff=0.8, nclusters=5)
  
  CNAdata <- NULL
  LRR_long <- NULL
  cp.index_long <- NULL
  sample_long <- NULL
  cluster_long <- NULL
  for (s in 1:ncol(LRR)){
    g=Clusters[s]
    cp.index1 <- 1+which(abs(diff(mean.matrix[g,],1))>=cutoff)
    cp.index <- unique(c(1, cp.index1, Chr_cp_index, nrow(LRR)))
    cp.index <- cp.index[order(cp.index)]
    cp.index
    
    LRR_long <- c(LRR_long, LRR[,s])
    cp.index_long <- rbind(cp.index_long, data.frame(cp.index=cp.index,
                                                     cp.index_long1=cp.index+nrow(LRR)*(s-1),
                                                     sample_index = s, 
                                                     cluster_index=g))
    sample_long <- c(sample_long, rep(s, each=nrow(LRR)))
    cluster_long <- c(cluster_long, rep(s, each=nrow(LRR)))
    
    #cp.index_long_update <- cp.index_long+nrow(LRR)*(s-1)
    #if ((length(cp.index1) < 1) | (length(cp.index1) > 10000)){ next }
    #cp.index <- unique(c(cp.index1, Chr_cp_index, QC_cp_index))
    #cp.index <- unique(c(cp.index1, Chr_cp_index))
    #cp.index <- unique(c(cp.index1))
  }    
  
  x.inv  <- try( res <- CNAcluster(Y = LRR_long, cp=cp.index_long$cp.index_long1, L, init),  silent=TRUE)
  res1 <- data.frame(start_long=res$CNA.start, end_long=res$CNA.end, state=res$CNA.state)
  cp.index_long1 <- merge(cp.index_long, res1, by.x="cp.index_long1", by.y="start_long", all.x=T)
  aa <- cp.index_long1[!is.na(cp.index_long1$end_long), ]
  
  colnames(aa)[c(2, 3, 4)] <- c("start", "sampleID", "Cluster")
  rownames(cp.index_long) <- as.character(cp.index_long$cp.index_long1)
  aa$end <- cp.index_long[as.character(aa$end_long), "cp.index"]
  aa$samplename=colnames(LRR)[aa$sampleID]
  aa$chr=rep(QC_ref@seqnames@values, QC_ref@seqnames@lengths)[aa$start]
  aa$start.coor=QC_ref@ranges@start[aa$start]
  aa$end.coor=(QC_ref@ranges@start+QC_ref@ranges@width)[aa$end]
  aa$width_bins=aa$end-aa$start+1
  CNAdata1 <- aa    

  CNAdata1 <- CNAdata1[CNAdata1$width_bins > 2,]
  CNAdata1 <- CNAdata1[CNAdata1$width_bins < 10000,]
  
  #CNAdata <- rbind(CNAdata, CNAdata1)
  #}
  return(CNAdata1)
}




#' @title  Generate initial parameters to cluster CNA states
#' 
#' @description This function clusters the identified change-points to make final CNA calling. The potential CNA segments between two neighbor candidate change-points are assigned to different copy number states according to the estimated mean matrix from FLCNA R function. We use three clusters including duplication, normal state and deletion. A Gaussisan Mixture Model based clustering strategy was applied to assign each segment to the most likely cluster/state.
#'
#' @param mean.matrix The cluster mean matrix estimated from FLCNA R function.
#' @param LRR log2R data after normlization.
#' @param Clusters Cluster index for each cell indentified from FLCNA R function.
#' @param ref Reference file.
#' @param cutoff Cutoff value to further control the number of CNAs, the larger value of cutoff, the smaller number of CNAs. The default is 0.35.
#' @param nclusters Number of CNA states.
#' 
#' 
#' @return Prior parameters used for GMM to cluster CNA states.
#' \item{priors}{Prior parameters used for GMM}
#'
#' @export
#' 
Para_init<-function (mean.matrix, LRR, ref, Clusters, cutoff = 0.8, nclusters=5){
  library("mclust")
  #Chr_index <- as.numeric(rep(ref@seqnames@values, ref@seqnames@lengths))
  Chr_index <- as.numeric(gsub("^.{0,3}", "", rep(ref@seqnames@values, ref@seqnames@lengths)))
  Chr_cp_index <- 1+which(abs(diff(Chr_index))>=1)
  
  
  CNAdata <- vector(mode = "list", length = nrow(mean.matrix))
  # g<-1
  seg_means<-c()
  for (i in 1:ncol(LRR)) {
    g=Clusters[i]
    cp.index1 <- 1 + which(abs(diff(mean.matrix[g, ], 1)) >= 
                             cutoff)
    cp.index <- unique(c(cp.index1, Chr_cp_index))
    cp.index <- sort(cp.index)
    for (s in 1:(length(cp.index)-1)){
      seg_means<-c(seg_means,mean(LRR[cp.index[s]:cp.index[s+1], g]))
    }
  }
  
  st = nclusters
  yMclust <- Mclust(seg_means,G=st,verbose = FALSE)
  init.clust<-yMclust$classification
  mus<-c()
  sds<-c()
  ps <-c()
  for (s1 in 1:st){
    mus<-c(mus,mean(seg_means[init.clust==s1]))
    sds<-c(sds,sd(seg_means[init.clust==s1]))
    ps<-c(ps,sum(init.clust==s1)/length(init.clust))
  }
  reorder<-order(mus)
  mu=mus[reorder]
  sd=sds[reorder]
  p=ps[reorder]
  priors <- list(p = p, mu = mu, sigma = sd)
  return(priors=priors)
  # return(list(priors=priors,init.clust=init.clust))
}





#' @title CNAcluster
#'
#' @description This function clusters CNAs into different states using a Gaussian Mixed Model based clustering strategy.
#' @param Y The numeric vector of the intensities of markers, which is the estimated mean vector in our study.
#' @param cp The numeric vector of the position index for the identified change-points.
#' @param L Repeat times in the EM algorithm, defaults to 100.
#' 
#' 
#' @return The return is the clustered CNA segments with the start position and end position, length of the CNA and the copy number states (duplication or deletion). It also returns a vector of final candidates of change-points.
#' \item{newcp}{The final list of change-points.}
#' \item{h}{The bandwidth used for the identification of change-points.}
#' \item{CNA.state}{Copy number state for each CNA.}
#' \item{CNA.start}{Start position of each CNA.}
#' \item{CNA.end}{End position of each CNA.}
CNAcluster <-function(Y, cp, L, init) {
  
  st = 5
  #mu = c(-2.3, -0.75, 0, 0.5, 1.0)
  #mu = c(-0.89, -0.24, 0.01, 0.25, 0.76)
  p=init$p
  mu=init$mu
  sigma=init$sigma
  #p = rep(1/st, st)
  #p = c(0.05, 0.05, 0.8, 0.05, 0.05)
  
  #sigma = rep(0.1, st)
  priors <- list(p = p, mu = mu, sigma = sigma)
  EM = gausianMixture(x=Y, cp, priors=priors, L, st=st)
  EM$state.new
  
  newcp = EM$cp.final   
  h = EM$index.final
  CNA.state <- getState(EM = EM)
  return(list(newcp = newcp, h = h, CNA.state = CNA.state$CNA.state, 
              CNA.start = CNA.state$CNA.start, CNA.end = CNA.state$CNA.end))
}



#' @title CNA states
#'
#' @description This function uses output of Gaussian Mixture Model to obtain different CNA states.
#' @param EM The output of Gaussian Mixture Model for clustering.
#' 
#' 
#' @return The return is the estimated CNA information.
#' \item{CNA.state}{Copy number state for each CNA.}
#' \item{CNA.start}{Start position of each CNA.}
#' \item{CNA.end}{End position of each CNA.}
getState <-function (EM = EM) {
  state      = EM$state.new
  cp.f       = EM$cp.final
  start.index  = which(state!=3)  
  CNA.start     = cp.f[start.index]
  CNA.start = CNA.start[!is.na(CNA.start)]
  CNA.end = cp.f[start.index+1]
  CNA.end = CNA.end[!is.na(CNA.end)]
  CNA.state = state[start.index]
  
  CNA.state[which(CNA.state == 1)] = "Del.d"
  CNA.state[which(CNA.state == 2)] = "Del.s"
  CNA.state[which(CNA.state == 4)] = "Dup.s"
  CNA.state[which(CNA.state == 5)] = "Dup.d"
  
  return (list(CNA.state = CNA.state, CNA.start = CNA.start, CNA.end = CNA.end))
}


