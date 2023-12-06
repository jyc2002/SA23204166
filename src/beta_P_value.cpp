#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
//' @title Global test of beta for moderate-dimensional inferences using Rcpp
//' @name beta_P_value
//' @description Global test of beta for moderate-dimensional inferences using Rcpp
//' @param X a numeric matrix
//' @param Y a numeric vector
//' @param beta the null beta
//' @import microbenchmark
//' @return P value of beta global test in moderate-dimensional
//' @import Rcpp
//' @import RcppArmadillo
//' @examples
//' \dontrun{
//' n=1000
//' pn=900
//' beta<-runif(pn)
//' sigma<-runif(n)
//' X<-matrix(rnorm(pn*n),nrow = n,ncol = pn)
//' epsilon<-numeric(n)
//' for(j in 1:n)epsilon[j]<-rnorm(1,mean = 0,sd = sigma)
//' Y<-X%*%beta+epsilon
//' beta_P_value(Y,X,rep(0,pn))
//' beta_P_value(Y,X,beta)
//' }
//' @export
// [[Rcpp::export]]
double beta_P_value(arma::vec Y, arma::mat X,arma::vec beta) {
  int n = X.n_rows;
  int pn = X.n_cols;
  arma::vec beta_hat = solve(trans(X) * X, trans(X) * Y);
  arma::mat XTX_inv = inv(trans(X) * X);
  double sigma_hat = as_scalar(trans(Y - X * beta_hat) * (Y - X * beta_hat) / (n - pn));
  double var = 4 * sigma_hat * as_scalar(trans(beta_hat) * XTX_inv * beta_hat) -
    2 * pow(sigma_hat, 2) * sum(diagvec(XTX_inv * XTX_inv)) +
    2 * pow(sigma_hat, 2) * pow(sum(diagvec(XTX_inv)), 2) / (n - pn);
  double obs = (as_scalar(trans(beta_hat-beta) * (beta_hat-beta)) - sigma_hat*sum(diagvec(XTX_inv))) / sqrt(var);
  double P_value = 2 * (1 - R::pnorm(std::abs(obs), 0.0, 1.0, true, false));
  return P_value;
}