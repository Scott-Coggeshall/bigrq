main <- function(X, y, M, intercept, maxiter){
  
   
  iter <- 1
  beta_i <- eta_i <- matrix(0, nrow = p, ncol = M)
  beta_avg <- eta_avg <- rowMeans(beta_i)
  while(iter < maxiter){
    
   iter_run <- foreach(dat_i = dat_iter) %dopar%{
     
     
     .Call()
     
     
   } 
    
    
  }
  
  
  
}