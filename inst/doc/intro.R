## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo = FALSE-------------------------------------------------------------
library("SA23204166")

## ----eval=FALSE---------------------------------------------------------------
#  RS_PA <- function(X,T,a) {
#    n=dim(X)[1]
#    sig=svd(X)$d
#    sig0=matrix(rep(0,n*T),nrow = T)
#    for (t in 1:T) {
#      mat=matrix(2*rbinom(n^2,1,0.5)-1,nrow = n)
#      upper_triangle <- upper.tri(mat)
#      s_m <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
#      s_m[upper_triangle] <- mat[upper_triangle]
#      s_m[lower.tri(s_m)] <- t(s_m)[lower.tri(s_m)]
#      sig0[t,]=svd(X*s_m)$d
#    }
#    sig1=apply(sig0, 2, quantile, probs = a)
#    k=which(sig < sig1)[1]-1
#    return(k)
#  }

## -----------------------------------------------------------------------------
set.seed(0)
mat=matrix(rbinom(50^2,1,0.2),ncol=50)
upper_triangle <- upper.tri(mat)
sim_data <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
sim_data[upper_triangle] <- mat[upper_triangle]
sim_data[lower.tri(sim_data)] <- t(sim_data)[lower.tri(sim_data)]

## -----------------------------------------------------------------------------
print(RS_PA(sim_data,300,0.95))

## -----------------------------------------------------------------------------
data(dolphins)
RS_PA(dolphins,300,0.95)

## -----------------------------------------------------------------------------
data(football)
RS_PA(football,300,0.95)

## -----------------------------------------------------------------------------
data(karate)
RS_PA(karate,300,0.95)

## -----------------------------------------------------------------------------
data(polbooks)
RS_PA(polbooks,300,0.95)

## ----echo = FALSE-------------------------------------------------------------
library(Rcpp)
library(RcppArmadillo)

## ----eval=FALSE---------------------------------------------------------------
#  #include <RcppArmadillo.h>
#  using namespace Rcpp;
#  // [[Rcpp::depends(RcppArmadillo)]]
#  // [[Rcpp::plugins(cpp11)]]
#  // [[Rcpp::export]]
#  double beta_P_value(arma::vec Y, arma::mat X,arma::vec beta) {
#    int n = X.n_rows;
#    int pn = X.n_cols;
#    arma::vec beta_hat = solve(trans(X) * X, trans(X) * Y);
#    arma::mat XTX_inv = inv(trans(X) * X);
#    double sigma_hat = as_scalar(trans(Y - X * beta_hat) * (Y - X * beta_hat) / (n - pn));
#    double var = 4 * sigma_hat * as_scalar(trans(beta_hat) * XTX_inv * beta_hat) -
#      2 * pow(sigma_hat, 2) * sum(diagvec(XTX_inv * XTX_inv)) +
#      2 * pow(sigma_hat, 2) * pow(sum(diagvec(XTX_inv)), 2) / (n - pn);
#    double obs = (as_scalar(trans(beta_hat-beta) * (beta_hat-beta)) - sigma_hat*sum(diagvec(XTX_inv))) / sqrt(var);
#    double P_value = 2 * (1 - R::pnorm(std::abs(obs), 0.0, 1.0, true, false));
#    return P_value;
#  }

## -----------------------------------------------------------------------------
set.seed(0)
n=1000
pn=900
beta<-runif(pn)
sigma<-runif(n)
X<-matrix(rnorm(pn*n),nrow = n,ncol = pn)
epsilon<-numeric(n)
for(j in 1:n)epsilon[j]<-rnorm(1,mean = 0,sd = sigma)
Y<-X%*%beta+epsilon

## -----------------------------------------------------------------------------
print(beta_P_value(Y,X,rep(0,pn)))

## -----------------------------------------------------------------------------
print(beta_P_value(Y,X,beta))

