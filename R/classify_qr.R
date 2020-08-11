classify_qr <- function(newx, coef, y_transform, threshold){
  
  
  predictions <- newx%*%coef
  
  predictions_transformed <- y_transform(predictions)
  
  (predictions_transformed  - threshold) > 0
  
}