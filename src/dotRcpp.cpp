// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

//' Matrix calculation in RcppArmadillo.
//'
//' @param Am matrix
//' @param Cm matrix
//' @return Matrix calculation, that is \code{inv(Am)%*%Cm}
// [[Rcpp::export(flowCalcCpp)]]
arma::mat flowCalcCpp(arma::mat Am, arma::mat Cm) 
{
         arma::mat B = inv(Am) * Cm;
         return B;
}