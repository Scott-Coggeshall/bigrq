#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]

List ParamUpdates(NumericVector betar, NumericVector etar, NumericMatrix xr, NumericVector yr, const betaavgr, const etaavgr){
//arma::vec yi=y.subvec(ni*i,ni*i+ni-1),ui=u.subvec(ni*i,ni*i+ni-1);
//arma::vec beta_i=betai.col(i),ri = r.subvec(ni*i,ni*i+ni-1);
//arma::mat xi=x.rows(ni*i,ni*i+ni-1);
xbetai=alpha*xr*betar+(1-alpha)*(yr-rr);
rr=shrinkcpp1(ur/rho+yr-xbetai-.5*(2*tau-1)/(n*rho),.5*arma::ones<arma::vec>(ni)/(n*rho));
//update betai
betar.col(i)=dat.slice(i)*(xi.t()*(yi-r.subvec(ni*i,ni*i+ni-1)+ui/rho)-etai.col(i)/rho+beta);
//update u and eta
u.subvec(ni*i,ni*i+ni-1)=u.subvec(ni*i,ni*i+ni-1)+rho*(yi-xbetai-r.subvec(ni*i,ni*i+ni-1));
etai.col(i)=etai.col(i)+rho*(betai.col(i)-beta);


}