#' High Dimensional BIC.
#' 
#' @param rq_coefs matrix of coefficients.
#' @param tau scalar between 0 and 1 corresponding to the conditional quantile.
#' @param x design matrix used in fitting the model.
#' @param y vector containing the values of the outcome used in fitting the model.
#' 
#' @return a vector of HBIC values
hbic <- function(rq_coefs, tau, x, y ){
  
  n <- length(y)
  p <- nrow(rq_coefs)
  checkloss <- check_loss((y - crossprod(rq_coefs, x)), tau)
  df <- colSums(rq_coefs != 0)
  
  log(checkloss) + log(n)*log(p)/(2*n)*df
  
  
}