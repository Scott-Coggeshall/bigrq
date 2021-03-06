#' High Dimensional BIC.
#' 
#' @param coefs matrix of coefficients.
#' @param tau scalar between 0 and 1 corresponding to the conditional quantile.
#' @param x design matrix used in fitting the model.
#' @param y vector containing the values of the outcome used in fitting the model.
#' 
#' @export
#' @return a vector of HBIC values
hbic <- function(coefs, tau, x, y, largep ){
  
  n <- length(y)
  p <- nrow(coefs)
  checkloss <- colSums(check_loss((y - x%*%coefs), tau))
  df <- colSums(coefs != 0)
  
  if(largep){
  log(checkloss) + log(log(n))*log(p)/(2*n)*df
  } else{
    
    log(checkloss) + log(n)/(2*n)*df
  }
  
  
}