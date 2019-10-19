#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


//' @export
// [[Rcpp::export]]
arma::mat matInv(arma::mat x) {
  return x.i();
}


