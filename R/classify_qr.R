#' Binary Classifications from Quantile Regression
#' 
#' \code{classify_qr} creates binary classifications from a quantile
#' regression fit by comparing predicted values to a prespecified threshold.
#' 
#' @param newx a design matrix containing the input values for predictions.
#' @param coef a vector or matrix containing the estimated regression coefficient values. 
#' If a matrix (e.g. if you want classifications for model fits produced by a 
#' sequence of tuning parameter values), then the each column should correspond
#' to a different fit.
#' @param y_transform a vectorized function to be applied to the predictions.
#' @param threshold the threshold to compare the predictions to.
#' 
#' @export
#' @return a vector or matrix of classifications.
classify_qr <- function(newx, coef, y_transform, threshold){
  
  
  predictions <- newx%*%coef
  
  predictions_transformed <- y_transform(predictions)
  
  (predictions_transformed  - threshold) > 0
  
}