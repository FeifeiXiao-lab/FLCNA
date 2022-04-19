#' @title Gaussian Mixture Model for CNA clustering
#' 
#' @description Gaussian Mixture Model is applied to assign each segment to the most likely cluster/state.
#'
#' @param x The vector of the estimated mean of markers.
#' @param cp The vector of the marker index of the identified change-points.
#' @param priors Given initial parameters for the EM algorithm.
#' @param L Repeat times in the EM algorithm. Defaults to 100.
#' @param st Number of assumed states in the EM algorithm.
#' 
#' 
#' @return The return is the clustered CNA segments with the start position and end position using CNA marker index, and the copy number states. It also returns a vector of final candidates of change-points.
#' \item{p.final}{Probability of falling into each state for each CNA segment after convergence.}
#' \item{mu.final}{Segment means of each state after convergence.}
#' \item{cp.final}{List of change-points after EM algorithm.}
#' \item{index.final}{The index of change-points.}
#' \item{state.new}{Assigned copy number state for each CNA.}
#' 
gausianMixture <- function(x, cp, priors, L, st) {
  
  T      <- length(x)
  cp.new <- cp[which(cp!=T)]
  index  <- names(cp.new)
  start  <- cp.new+1
  start <- start[-(length(start))]
  end    <- cp.new[-1]
  len    <- end - start + 1
  N      <- length(len)
  means  <- vector()
  sum.x.sq <- vector()
  state    <- vector()
  for (i in 1:N) {
    ##Find segment means
    means[i] <- mean(x[start[i]:end[i]])
    ##Find segment sum.x.sq
    sum.x.sq[i] <- sum(x[start[i]:end[i]]^2)
  }
  para.new <- updateEM(p.in=priors$p, mu.in=priors$mu, sigma.in=priors$sigma, means, sum.x.sq, N, len, st)
  for (iter in 1:L){
    para.new <- updateEM(p.in=para.new$p.new, mu.in=para.new$mu.new, sigma.in=para.new$sigma.new, means, sum.x.sq, N, len, st)
  }
  if (is.na(para.new$mu.new[1])){
    para.new <- updateEM(p.in=priors$p, mu.in=priors$mu, sigma.in=priors$sigma, means, sum.x.sq, N, len, st)
  }
  for (pt in 1:N) {
    state[pt] <- which.max(para.new$p[pt,])
  }
  index.change <- which((state[1:(N-1)] == state[2:N])== "TRUE")
  if (length(index.change) == 0) {
    cp.final     <- cp.new
    index.final  <- index
    state.new    <- state
  }else{
    cp.final     <- cp.new[-(index.change+1)]
    index.final  <- index[-(index.change+1)]
    state.new    <- state[-index.change]
  }
  # state.new <- state.new[c(length(state.new),1:(length(state.new)-1))]
  return(list(p.final = para.new$p.new, mu.final = para.new$mu.new, sigma.final = para.new$sigma.new, cp.final = cp.final, index.final = index.final, state.new = state.new))
}


#' @title Update parameters using EM algorithm
#' 
#' @description In the Gaussian Mixture Model, parameters will be updated based on EM algorithm.
#'
#' @param p.in Initial probability for each CNA cluster.
#' @param mu.in Initial mean value for each CNA cluster.
#' @param sigma.in Initial variance for each CNA cluster.
#' @param means Mean value vector for each segment.
#' @param sum.x.sq Sum of squared mean values for each segment.
#' @param N Number of candiate CNAs.
#' @param len Width of candiate CNAs.
#' @param st Number of assumed states in the EM algorithm.

#' 
#' @return The return is the updated parameters using EM algorithm
#' \item{p.new}{Updated probability for each CNA cluster.}
#' \item{mu.new}{Updated mean value for each CNA cluster.}
#' \item{sigma.new}{Updated variance for each CNA cluster.}
#' 
updateEM <- function (p.in, mu.in, sigma.in, means, sum.x.sq, N, len, st) {
    ##Calcute the prob of each segment belong to each state
    p <- dens <- matrix(NA, N, st)
    for (i in 1:N) {
      a <- rep(NA, st)
      for (j in 1:st) {
        dens[i,j] <- dnorm(means[i], mu.in[j], sqrt(sigma.in[j]/len[i]), log=TRUE)
        a[j]      <- log(p.in[j]) + dens[i,j]
      }
      max.a     <- max(a)
      for (k in 1:st) {
        p[i,k] <- (exp(a[k]-max.a))/sum(exp(a-max.a))
      }
    }
    ##Update p, mu, sigma
    p <- na.omit(p)
    p.new <- mu.new <- sigma.new <- rep(NA, st)
    for (k in 1:st) {
      p.new[k]  <- sum(p[,k])/N
      mu.new[k] <- (sum(means*len*p[,k]))/(sum(len*p[,k]))
      sigma.new[k] <- sum(p[,k]*(sum.x.sq-2*len*mu.new[k]*means+len*mu.new[k]^2))/(sum(p[,k]*len))
    }
    return(list(p.new = p.new, mu.new = mu.new, sigma.new = sigma.new, p = p))
  }
