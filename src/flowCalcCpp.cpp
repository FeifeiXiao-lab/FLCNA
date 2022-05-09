// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


//' Matrix calculation in RcppArmadillo.
//'
//' @param Am matrix
//' @param Cm matrix
//' @return Matrix, that is \code{inv(Am)%*%Cm}
// [[Rcpp::export]]
arma::mat flowCalcCpp(const arma::mat Am, const arma::mat Cm) 
{
         arma::mat B = inv(Am) * Cm;
         return B;
}