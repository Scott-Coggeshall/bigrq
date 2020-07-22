#' Woodbury Matrix Inversion Identity
#' 
#' \code{compute_woodbury} uses the Woodbury matrix inversion identity to compute
#' the inverse of (I + M) where M is the cross-product of the design matrix X
#' and I is the identity matrix.
#' 
#' @param mat a matrix.
#' 
#' @return a matrix, the result of inverting \code{mat}.
#' @export

compute_woodbury <- function(mat){
  
  n <- nrow(mat)
  p <- ncol(mat)
  
  I_n <- diag(n)
  
  I_p <- diag(p)
  
  B_inv <- solve(I_n + tcrossprod(mat))
  
  I_p - t(mat)%*%B_inv%*%mat
  
}