#' Check Loss 
#' 
#' \code{check_loss} computes the check loss function for a given quantile tau.
#' @param x a numeric scalar, vector, or matrix.
#' @param tau a numeric scalar between 0 and 1.
check_loss <- function(x, tau){
  
  
  x*(tau - (x < 0))
  
}