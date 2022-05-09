#include <Rcpp.h>
using namespace Rcpp;

//' Matrix calculation in Rcpp.
//'
//' @param Am matrix
//' @param Cm matrix
//' @return Matrix calculation, that is \code{inv(Am)%*%Cm}
// [[Rcpp::export]]
arma::mat flowCalcCpp(const arma::mat &Am, const arma::mat &Cm) 
{
         arma::mat B = inv(Am) * Cm;
         return arma::mat B;
}