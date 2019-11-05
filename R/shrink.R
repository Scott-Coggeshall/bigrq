shrink <- function(u, v){
  
  
  (1 + sign(u - v))/2*(u - v) - (1 + sign(-u - v))/2*(-u-v)
  
  
  
}
