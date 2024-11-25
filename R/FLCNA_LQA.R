Rcpp::cppFunction(depends = "RcppArmadillo", code = '
         Rcpp::List flowCalcCpp(const arma::mat &Am, const arma::mat &Cm) {
         arma::mat B = inv(Am) * Cm;
         return Rcpp::List::create( Rcpp::Named("Imp") = B);
      }')




#' @keywords internal
## ------------ dimension-wise, based on local quadratic approximation
FLCNV_LQA <- function(k, K1, index.max, y, mu.t.all, mu.no.penal, sigma.all, alpha, lambda, iter.num = 20, eps.LQA = 1e-5, eps.diff = 1e-5){
  
  ## mu.t.all : K by p mean matrix from previous EM-step
  ## mu.no.penal : K by p mean matrix of unpenalized estimates (lambda=0)
  ## y: n by p data matrix
  ## sigma.all: p by p diagnal covariance matrix
  ## alpha: n by K posterior probability matrix
  ## iter.num: max iterations in local quadratic approximation (LQA)
  ## eps.LQA: LQA stop criterion
  ## eps.diff: lower bound of mean difference
  ## lambda: tuning parameter
  index.k <- which(index.max == k)
  if (length(index.k)<=1){
    out <- mu.t.all[k,]
  } else {
    mu.t.k <- mu.t.all[k,]
    mu.no.k <- mu.no.penal[k,]

    y.k <- y[index.k,]
    alpha.k <- as.matrix(alpha[index.k, k])
  
    n.k <- dim(y.k)[1]
    P <- ncol(y)
    # y.k.long <- as.matrix(as.vector(y.k))
  
    # A <- diag(c(rep(alpha.k, P)/(2*rep(sigma.all, n.k))))
    # X <- diag(P)[rep(1:P, each = n.k),]
    
    AA <- data.frame(value=c(rep(alpha.k, P)/(2*rep(sigma.all, n.k))), group=rep(1:P, each=n.k))
    mat1 <-  diag(aggregate(AA$value, by=list(Category=AA$group), FUN=sum)$x)
    
    BB <- data.frame(value=AA$value*as.vector(y.k), group=rep(1:P, each=n.k))
    vec1 <-  as.matrix(aggregate(BB$value, by=list(Category=BB$group), FUN=sum)$x)
    
    # mat1 <- t(X)%*%A%*%X
    # vec1 <- t(X)%*%A%*%y.k.long
    
    combn.mat <- matrix(0, P-1, P)
    for (s.i in 1:nrow(combn.mat)){
      combn.mat[s.i,s.i] <- -1
      combn.mat[s.i,s.i+1] <- 1
    }

    ## if distance < eps.diff, make them equal.
    dist.mu.no.k <- abs(combn.mat%*%mu.no.k)
    dist.mu.no.k[dist.mu.no.k < eps.diff] <- eps.diff
    log.tau <- c(-log(abs(combn.mat%*%mu.no.k)))
  
    mu.k.hat <- matrix(0, nrow=iter.num+1,ncol=P)
    mu.k.hat[1,] <- mu.t.k
  
    for(s in 1:iter.num){
      dist.old.tmp0 <- abs(combn.mat%*%mu.k.hat[s,])
      dist.old.tmp0[dist.old.tmp0 < eps.diff] <- eps.diff
    
      c.kk <- as.vector(sqrt(exp(log.tau - log(2*dist.old.tmp0))))
    
      # B <- t(combn.mat)%*%t(diag(c.kk))%*%diag(c.kk)%*%combn.mat
      B <- matrix(0, P, P)
      B[1,1] <- c(c.kk^2)[1]
      B[length(c.kk), length(c.kk)+1] <- -c(c.kk^2)[length(c.kk)]
      B[length(c.kk)+1, length(c.kk)] <- -c(c.kk^2)[length(c.kk)]
      B[length(c.kk)+1, length(c.kk)+1] <- c(c.kk^2)[length(c.kk)]
      
      for (i in 1:(length(c.kk)-1)){
        B[i, i+1] <-  -c(c.kk^2)[i]
        B[i+1, i] <-  -c(c.kk^2)[i]
        B[i+1, i+1] <- c(c.kk^2)[i] + c(c.kk^2)[i+1]
      }
    
      # mu.k.hat[(s+1),] <- solve(mat1+lambda*B)%*%vec1
      system.time(C <- flowCalcCpp(mat1+lambda*B, vec1))

      # D <- solve(mat1+lambda*B)
      
      mu.k.hat[(s+1),] <- C$Imp
      ## stopping criterion = relative distance (+1) needed when the vector is 0.
      ## for computing stability, let small value to be exact zero
      if(sum(abs(mu.k.hat[(s+1),] - mu.k.hat[s,]))/(sum(abs(mu.k.hat[s,])) + 1) < eps.LQA){
        out <- mu.k.hat[(s+1),]
        out[which(abs(out) < eps.diff)] <- 0
        break
      }
    
      else if(s == iter.num){
        out <- mu.k.hat[(s+1),]
        out[which(abs(out) < eps.diff)] <- 0
        warning(paste('LQA for estimating mu does not converge for K=', K1, 'lambda=', lambda, sep=''))
      }
    }
  }  
  return(out)
}






#' @keywords internal
## ------------ dimension-wise, based on local quadratic approximation
FLCNV_LQA_per_chr <- function(k, K1, index.max, y, mu.t.all, mu.no.penal, sigma.all, alpha, lambda, chromosome.id, iter.num = 20, eps.LQA = 1e-5, eps.diff = 1e-5){
  
  ## mu.t.all : K by p mean matrix from previous EM-step
  ## mu.no.penal : K by p mean matrix of unpenalized estimates (lambda=0)
  ## y: n by p data matrix
  ## sigma.all: p by p diagnal covariance matrix
  ## alpha: n by K posterior probability matrix
  ## iter.num: max iterations in local quadratic approximation (LQA)
  ## eps.LQA: LQA stop criterion
  ## eps.diff: lower bound of mean difference
  ## lambda: tuning parameter
  index.k <- which(index.max == k)
  if (length(index.k)<=1){
    out <- mu.t.all[k,]
  } else {
    mu.t.k <- mu.t.all[k,]
    mu.no.k <- mu.no.penal[k,]
    
    y.k <- y[index.k,]
    alpha.k <- as.matrix(alpha[index.k, k])
    
    n.k <- dim(y.k)[1]
    P <- ncol(y)
    # y.k.long <- as.matrix(as.vector(y.k))
    
    # A <- diag(c(rep(alpha.k, P)/(2*rep(sigma.all, n.k))))
    # X <- diag(P)[rep(1:P, each = n.k),]
    
    AA <- data.frame(value=c(rep(alpha.k, P)/(2*rep(sigma.all, n.k))), group=rep(1:P, each=n.k))
    mat1 <-  diag(aggregate(AA$value, by=list(Category=AA$group), FUN=sum)$x)
    
    BB <- data.frame(value=AA$value*as.vector(y.k), group=rep(1:P, each=n.k))
    vec1 <-  as.matrix(aggregate(BB$value, by=list(Category=BB$group), FUN=sum)$x)
    
    
    # mat1 <- t(X)%*%A%*%X
    # vec1 <- t(X)%*%A%*%y.k.long
    
    combn.mat <- matrix(0, P-1, P)
    for (s.i in 1:nrow(combn.mat)){
      combn.mat[s.i,s.i] <- -1
      combn.mat[s.i,s.i+1] <- 1
    }
    combn.mat_mod <- combn.mat
    combn.mat_mod[ which(diff(chromosome.id) != 0),] <- 0
    
    
    ## if distance < eps.diff, make them equal.
    dist.mu.no.k <- abs(combn.mat%*%mu.no.k)
    dist.mu.no.k[dist.mu.no.k < eps.diff] <- eps.diff
    log.tau <- c(-log(abs(combn.mat%*%mu.no.k)))
    
    mu.k.hat <- matrix(0, nrow=iter.num+1,ncol=P)
    mu.k.hat[1,] <- mu.t.k
    
    for(s in 1:iter.num){
      dist.old.tmp0 <- abs(combn.mat_mod%*%mu.k.hat[s,])
      dist.old.tmp0[dist.old.tmp0 < eps.diff] <- eps.diff
      
      c.kk <- as.vector(sqrt(exp(log.tau - log(2*dist.old.tmp0))))
      
      # B <- t(combn.mat)%*%t(diag(c.kk))%*%diag(c.kk)%*%combn.mat
      B <- matrix(0, P, P)
      B[1,1] <- c(c.kk^2)[1]
      B[length(c.kk), length(c.kk)+1] <- -c(c.kk^2)[length(c.kk)]
      B[length(c.kk)+1, length(c.kk)] <- -c(c.kk^2)[length(c.kk)]
      B[length(c.kk)+1, length(c.kk)+1] <- c(c.kk^2)[length(c.kk)]
      
      for (i in 1:(length(c.kk)-1)){
        B[i, i+1] <-  -c(c.kk^2)[i]
        B[i+1, i] <-  -c(c.kk^2)[i]
        B[i+1, i+1] <- c(c.kk^2)[i] + c(c.kk^2)[i+1]
      }
      
      # mu.k.hat[(s+1),] <- solve(mat1+lambda*B)%*%vec1
      system.time(C <- flowCalcCpp(mat1+lambda*B, vec1))
      
      # D <- solve(mat1+lambda*B)
      
      mu.k.hat[(s+1),] <- C$Imp
      ## stopping criterion = relative distance (+1) needed when the vector is 0.
      ## for computing stability, let small value to be exact zero
      if(sum(abs(mu.k.hat[(s+1),] - mu.k.hat[s,]))/(sum(abs(mu.k.hat[s,])) + 1) < eps.LQA){
        out <- mu.k.hat[(s+1),]
        out[which(abs(out) < eps.diff)] <- 0
        break
      }
      
      else if(s == iter.num){
        out <- mu.k.hat[(s+1),]
        out[which(abs(out) < eps.diff)] <- 0
        warning(paste('LQA for estimating mu does not converge for K=', K1, 'lambda=', lambda, sep=''))
      }
    }
  }  
  return(out)
}


