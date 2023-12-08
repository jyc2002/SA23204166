#' @name RS_PA
#' @title Random Signflips parallel analysis using R
#' @description Random Signflips parallel analysis using R
#' @param X the sample matrix
#' @param T the number of trials
#' @param a the percentile
#' @return k the number of leading singular values above the 𝛼-percentile of their signflipped analogues
#' @examples
#' \dontrun{
#'     mat=matrix(rbinom(50^2,1,0.2),ncol=50)
#'     upper_triangle <- upper.tri(mat)
#'     sim_data <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
#'     sim_data[upper_triangle] <- mat[upper_triangle]
#'     sim_data[lower.tri(sim_data)] <- t(sim_data)[lower.tri(sim_data)]
#'     RS_PA(sim_data,300,0.95)
#'     data(dolphins)
#'     RS_PA(dolphins,300,0.95)
#'     data(football)
#'     RS_PA(football,300,0.95)
#'     data(karate)
#'     RS_PA(karate,300,0.95)
#'     data(polbooks)
#'     RS_PA(polbooks,300,0.95)
#' }
#' @export
RS_PA <- function(X,T,a) {
  n=dim(X)[1]
  sig=svd(X)$d
  sig0=matrix(rep(0,n*T),nrow = T)
  for (t in 1:T) {
    mat=matrix(2*rbinom(n^2,1,0.5)-1,nrow = n)
    upper_triangle <- upper.tri(mat)
    s_m <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    s_m[upper_triangle] <- mat[upper_triangle]
    s_m[lower.tri(s_m)] <- t(s_m)[lower.tri(s_m)]
    sig0[t,]=svd(X*s_m)$d
  }
  sig1=apply(sig0, 2, quantile, probs = a)
  k=which(sig < sig1)[1]-1
  return(k)
}