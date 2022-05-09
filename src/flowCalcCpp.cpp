#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' Matrix calculation in RcppArmadillo.
//'
//' @param Am matrix
//' @param Cm matrix
//' @return Matrix, that is \code{inv(Am)%*%Cm}
// [[Rcpp::export(flowCalcCpp)]]
Rcpp::List flowCalcCpp_I(const arma::mat &Am, const arma::mat &Cm) 
{
    arma::mat B = inv(Am) * Cm;
    return Rcpp::List::create( Rcpp::Named("Imp") = B);
}