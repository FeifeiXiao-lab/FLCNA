#' @title Shared CNA output
#' 
#' @description This function clusters the identified change-points to make final CNA calling. The potential CNA segments between two neighbor candidate change-points are assigned to different copy number states according to the estimated mean matrix from FLCNA R function. We use three clusters including duplication, normal state and deletion. A Gaussisan Mixture Model based clustering strategy was applied to assign each segment to the most likely cluster/state.
#'
#' @param mean.matrix The cluster mean matrix estimated from FLCNA R function.
#' @param cutoff Cutoff value to further control the number of CNAs, the larger value of cutoff, the smaller number of CNAs. The default is 0.8.
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
CNA.out <- function(mean.matrix, LRR, Clusters, QC_ref, cutoff=0.80, L=100){
  
  Chr_index <- as.numeric(gsub("^.{0,3}", "", rep(QC_ref@seqnames@values, QC_ref@seqnames@lengths)))
  Chr_cp_index <- 1+which(abs(diff(Chr_index))>=1)

  #QC_cp_index <- 1+which(abs(diff(QC_ref@ranges@start,1))>median(QC_ref@ranges@width))

  CNAdata <- NULL
  for (s in 1:ncol(LRR)){
    g=Clusters[s]
    cp.index1 <- 1+which(abs(diff(mean.matrix[g,],1))>=cutoff)
    if ((length(cp.index1) < 1) | (length(cp.index1) > 10000)){ next }
    #cp.index <- unique(c(cp.index1, Chr_cp_index, QC_cp_index))
    cp.index <- unique(c(cp.index1, Chr_cp_index))
    cp.index <- cp.index[order(cp.index)]
    cp.index
    x.inv  <- try( res <- CNAcluster(Y = LRR[,s], cp=cp.index, L),  silent=TRUE)
    if ('try-error' %in% class(x.inv)) next
    if (length(x.inv$CNA.end)==0){
      next
    }
    CNAdata1 <- data.frame(sampleID=s,
                           samplename=colnames(LRR)[s],
                           Cluster=g,
                           chr=rep(QC_ref@seqnames@values, QC_ref@seqnames@lengths)[res$CNA.start],
                           start=res$CNA.start, 
                           end=res$CNA.end,
                           start.coor=QC_ref@ranges@start[res$CNA.start],
                           end.coor=(QC_ref@ranges@start+QC_ref@ranges@width)[res$CNA.end],
                           width_bins=(res$CNA.end-res$CNA.start+1),
                           state=res$CNA.state)
    CNAdata1 <- CNAdata1[CNAdata1$width_bins > 2,]
    CNAdata1 <- CNAdata1[CNAdata1$width_bins < 10000,]
    
    CNAdata <- rbind(CNAdata, CNAdata1)
  }
  return(CNAdata)
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
CNAcluster <-function(Y, cp, L) {
 
  st = 5
  mu = c(-2.3, -0.75, 0, 0.5, 1.0)
  p = rep(1/st, st)
  sigma = rep(0.1, st)
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


