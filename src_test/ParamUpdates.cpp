#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
//
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
List ParamUpdates(NumericVector betar, NumericVector etar, NumericMatrix xr, NumericVector yr, NumericMatrix dat, NumericVector ur, NumericVector beta, double rho, double alpha, double tau){
//arma::vec yi=y.subvec(ni*i,ni*i+ni-1),ui=u.subvec(ni*i,ni*i+ni-1);
//arma::vec beta_i=betai.col(i),ri = r.subvec(ni*i,ni*i+ni-1);
//arma::mat xi=x.rows(ni*i,ni*i+ni-1);
xbetai=alpha*xr*betar+(1-alpha)*(yr-rr);
rr=shrinkcpp1(ur/rho+yr-xbetai-.5*(2*tau-1)/(n*rho),.5*arma::ones<arma::vec>(ni)/(n*rho));
//update betai
betai = dat*(xr.t()*(yr-rr+ur/rho)-etar/rho + beta);
//update u and eta
ur=ur+rho*(yr-xbetai-rr);
etar=etar+rho*(betai-beta);

return List::create(Named("betai") = betai, Named("ur") = ur, Named("etar") = etar)

}