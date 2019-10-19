split_data <- function(data_set, n_chunks){
  
  indices <- rep_len(1:n_chunks, floor(nrow(data_set)/n_chunks) + (nrow(data_set) %% n_chunks))
  
  lapply(1:n_chunks, function(x) data_set[indices == x, ])
  
  
  
  
}

main <- function(dat, M, intercept, maxiter){
  
  dat_list <- split_data(data_set = dat, n_chunks = M)
  
  dat_inverses <- lapply(dat_list, function(x) solve(crossprod(x)))
  
  dat_iter <- foreach()
  iter <- 1
  beta_i <- eta_i <- matrix(0, nrow = p, ncol = M)
  beta_avg <- eta_avg <- rowMeans(beta_i)
  while(iter < maxiter){
    
   iter_run <- foreach(dat_i = dat_iter) %dopar%{
     
     
     .Call()
     
     
   } 
    
    
  }
  
  
  
}