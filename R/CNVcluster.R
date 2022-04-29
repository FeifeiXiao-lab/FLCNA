#' @title Shared CNA output
#' 
#' @description This function clusters the identified change-points to make final CNA calling. The potential CNA segments between two neighbor candidate change-points are assigned to different copy number states according to the estimated mean matrix from FLCNA R function. We use three clusters including duplication, normal state and deletion. A Gaussisan Mixture Model based clustering strategy was applied to assign each segment to the most likely cluster/state.
#'
#' @param mean.matrix The cluster mean matrix estimated from FLCNA R function.
#' @param cutoff Cutoff value to further control the number of CNAs, the larger value of cutoff, the smaller number of CNAs. The default is 0.35.
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
CNA.out <- function(mean.matrix, ref, cutoff=0.35, L=100){
  CNAdata <- vector(mode = "list", length = nrow(mean.matrix))
  for (g in 1:nrow(mean.matrix)){
    cp.index <- 1+which(abs(diff(mean.matrix[g,],1))>=cutoff)
    if ((length(cp.index) < 1) | (length(cp.index) > 1000)){ next }
    x.inv  <- try( res <- CNAcluster(Y = mean.matrix[g,], cp=cp.index, L),  silent=TRUE)
    if ('try-error' %in% class(x.inv)) next
    if (length(x.inv$CNA.end)==0){
      next
    }
    CNAdata[[g]] <- data.frame(state=res$CNA.state, 
                               start=ref@ranges@start[res$CNA.start], 
                               end=(ref@ranges@start+ref@ranges@width)[res$CNA.end], 
                               chr=rep(ref@seqnames@values, ref@seqnames@lengths)[res$CNA.start],
                               width_bins=(res$CNA.end-res$CNA.start+1))
  }
  return(CNAdata)
}



#' @title CNA output for each cell
#' 
#' @description This function clusters the identified change-points to make final CNA calling for each cell. The potential CNA segments between two neighbor candidate change-points are assigned to different copy number states according to the estimated mean matrix from FLCNA R function and log2R data for each cell. We use three clusters including duplication, normal state and deletion. A Gaussisan Mixture Model based clustering strategy was applied to assign each segment to the most likely cluster/state.
#'
#' @param mean.matrix The cluster mean matrix estimated from FLCNA R function.
#' @param log2R.NRC Log2R data from normalization of original read counts.
#' @param cluster.index Cluster index for all the cells.
#' @param cutoff Cutoff value to further control the number of CNAs, the larger value of cutoff, the smaller number of CNAs. The default is 0.35.
#' @param L Repeat times in the EM algorithm, defaults to 100.
#' 
#' 
#' @return The return is the clustered CNA segments by presenting the start position and end position using CNA marker index, and the copy number states.
#' \item{state}{The CNA states assigned.}
#' \item{start}{The start point for CNAs.}
#' \item{end}{The end point for CNAs.}
#' \item{width}{The width for CNAs.}
#' \item{sample}{Sample index.}
#'
#' @export
CNA.out.eachcell <- function(mean.matrix, log2R.NRC, cluster.index, cutoff=0.5, L=100){
  CNAdata <- NULL
  for (i in 1:ncol(log2R.NRC)){
    cp.index <- 1+which(abs(diff(mean.matrix[cluster.index[i],],1))>=cutoff)
    if (length(cp.index)==1){
      next
    }
    x.inv  <- try(res <- CNAcluster(Y = as.numeric(log2R.NRC[,i]), cp=cp.index, L),  silent=TRUE)
    if ('try-error' %in% class(x.inv)) next
    if (length(x.inv$CNA.end)==0){
      next
    }
    # tryCatch(res <- CNAcluster(Y = as.numeric(log2R.NRC[,i]), cp=cp.index, L), error = function(e), { next })
    CNA <- data.frame(state=res$CNA.state, start=res$CNA.start, end=res$CNA.end, width=(res$CNA.end-res$CNA.start+1), sample=i)
    CNAdata <- rbind(CNAdata, CNA)
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
 
  st = 3
  mu = c(-1, 0, 1)
  p = rep(1/st, st)
  sigma = rep(0.1, st)
  priors <- list(p = p, mu = mu, sigma = sigma)
  EM = gausianMixture(x=Y, cp, priors=priors, L, st=st)

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
  start.index  = which(state!=2)  
  CNA.start     = cp.f[start.index]
  CNA.start = CNA.start[!is.na(CNA.start)]
  CNA.end = cp.f[start.index+1]
  CNA.end = CNA.end[!is.na(CNA.end)]
  CNA.state = state[start.index]
  
  CNA.state[which(CNA.state == 1)] = "del"
  CNA.state[which(CNA.state == 3)] = "dup"

  return (list(CNA.state = CNA.state, CNA.start = CNA.start, CNA.end = CNA.end))
}


