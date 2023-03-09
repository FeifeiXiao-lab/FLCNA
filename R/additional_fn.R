#' @title dmvnorm_log
#' @description Used in sapply to find all the densities
#'
#' @param index Row index of mu.
#' @param mu K by p matrix, each row represents one cluster mean.
#' @param y n by p data matrix.
#' @param sigma p by p covariance matrix (assume same covariance for each cluster).
#' 
#' @importFrom mvtnorm dmvnorm
dmvnorm_log <- function(index, mu, sigma, y) {
  
  return(mvtnorm::dmvnorm(y, mu[index,], sigma, log=TRUE))
}


#' @title count.mu 
#' @description Computing the number of unique cluster means for each dimension, which was used in computing BIC or GIC.
#'
#' @param mu.j Mean vector.
#' @param eps.diff Lower bound of mean difference.
count.mu <- function(mu.j, eps.diff) {
  temp.dist <- as.matrix(dist(mu.j, method = 'manhattan'))
  ct <- length(mu.j[abs(mu.j)>eps.diff])
        ## initial counts (nonzero elements)
  ## --- if exists same means
  temp.dist[upper.tri(temp.dist, diag = T)] <- NA
  if(any(temp.dist < eps.diff, na.rm = T)) {
    temp.index <- which(temp.dist < eps.diff, arr.ind = TRUE)
    temp1 <- mu.j
    ## --- truncated distance so means are not exactly same, make them equal
    for(i in 1:dim(temp.index)[1]){
      temp1[temp.index[i,]] <- min(temp1[temp.index[i,]])
    }
    ct <- length(unique(temp1[abs(temp1)>eps.diff]))
  }
  return(ct)
}




# ------------------------------------------------------------------------------------------ #
dmvnorm_log_sapply<-function(seq.max,
                             y,
                             mu,
                             sigma.iter,
                             batch.size=50){
  nlength<-ncol(y)
  
  times<-nlength%/%batch.size
  redu <-nlength%%batch.size
  
  tmp.normal.density = sapply(c(1:seq.max), dmvnorm_log, 
                              y     = y[,1:batch.size], 
                              mu    = mu[, 1:batch.size], 
                              sigma = diag(sigma.iter[1:batch.size]))
  if (times>1){
    for (i in 2:times) {
      tmp.normal.density1 = sapply(c(1:seq.max), dmvnorm_log, 
                                   y     = y[,((i-1)*batch.size+1):(i*batch.size)], 
                                   mu    = mu[, ((i-1)*batch.size+1):(i*batch.size)], 
                                   sigma = diag(sigma.iter[((i-1)*batch.size+1):(i*batch.size)]))
      tmp.normal.density<-tmp.normal.density+tmp.normal.density1
    }
  }
  
  if (redu>0){
    # calculate dmvnorm_log for the remaining sequence
    tmp.normal.density1 = sapply(c(1:seq.max), dmvnorm_log, 
                                 y     = y[,(times*batch.size+1):nlength], 
                                 mu    = mu[, (times*batch.size+1):nlength], 
                                 sigma = diag(sigma.iter[(times*batch.size+1):nlength]))
    tmp.normal.density<-tmp.normal.density+tmp.normal.density1
  }
  
  
  return(tmp.normal.density)
  
}



# ------------------------------------------------------------------------------------------ #
dmvnorm_sapply<-function(y1,
                         mean1,
                         sigma1,
                         batch.size=50){
  
  if (is.matrix(y1)){nlength<-ncol(y1)}
  
  if(is.vector(y1)){
    nlength<-length(y1)
    y1<-matrix(y1,nrow = 1)
    mean1<-matrix(mean1,nrow = 1)
  }
  # dim(y1)
  # dim(mean1)
  times<-nlength%/%batch.size
  redu <-nlength%%batch.size
  
  tmp.normal.density = dmvnorm(x= y1[,1:batch.size], 
                               mean = mean1[1:batch.size], 
                               sigma = sigma1[1:batch.size,1:batch.size], 
                               log = TRUE)
  if (times>1){
    for (i in 2:times) {
      tmp.normal.density1 = dmvnorm(x= y1[,((i-1)*batch.size+1):(i*batch.size)], 
                                    mean = mean1[((i-1)*batch.size+1):(i*batch.size)], 
                                    sigma = sigma1[((i-1)*batch.size+1):(i*batch.size),((i-1)*batch.size+1):(i*batch.size)], 
                                    log = TRUE)
      tmp.normal.density<-tmp.normal.density+tmp.normal.density1
    }
  }
  
  if (redu>0){
    # calculate dmvnorm for the remaining sequence
    tmp.normal.density1 = dmvnorm(x= y1[,(times*batch.size+1):nlength], 
                                  mean = mean1[(times*batch.size+1):nlength], 
                                  sigma = sigma1[(times*batch.size+1):nlength,(times*batch.size+1):nlength], 
                                  log = TRUE)
    tmp.normal.density<-tmp.normal.density+tmp.normal.density1
  }
  
  return(tmp.normal.density)
  
}

