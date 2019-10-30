
#include <RcppArmadillo.h>
#include "header.hpp"


using namespace Rcpp;
using namespace arma;



//' This function updates the block-specific parameters.
//' 
//' @param betar a vector containing the current values for the block-specific beta parameters
//' @param etar a vector containing the current values for the block-specific eta parameters
//' @param xr a matrix containing the block-specific design matrix
//' @param yr a matrix containing the block-specific outcome vector
//' @param dat a matrix containing the block-specific matrix inverse of t(xr)%*%xr + I
//' @param ur a vector containing the current values for the block-specific u variables
//' @param beta a vector containing the current values for the global beta parameters
//' @param rho a numeric scalar 
//' @param alpha a numeric scalar
//' @param tau a numeric scalar
//' @param n an integer giving the overall sample size for the full data
//' @param ni an integer giving the block-specific sample size
//' @export
// [[Rcpp::export]]
List ParamUpdates(arma::vec betar, arma::vec etar, arma::mat xr, arma::vec yr, arma::mat dat, arma::vec ur, arma::vec rr, arma::vec beta, double rho, double alpha, double tau, int n, int ni){
//arma::vec yi=y.subvec(ni*i,ni*i+ni-1),ui=u.subvec(ni*i,ni*i+ni-1);
//arma::vec beta_i=betai.col(i),ri = r.subvec(ni*i,ni*i+ni-1);
//arma::mat xi=x.rows(ni*i,ni*i+ni-1);
arma::mat xbetai=alpha*xr*betar+(1-alpha)*(yr-rr);
arma::vec ri=shrinkcpp1(ur/rho+yr-xbetai-.5*(2*tau-1)/(n*rho),.5*arma::ones<arma::vec>(ni)/(n*rho));
//update betai
arma::vec betai = dat*(xr.t()*(yr-ri+ur/rho)-etar/rho + beta);
//update u and eta
arma::vec ui=ur+rho*(yr-xbetai-ri);
arma::vec etai=etar+rho*(betai-beta);

return List::create(Named("betai") = betai, Named("etai") = etai, Named("ui") = ui);
}