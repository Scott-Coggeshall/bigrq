split_data <- function(data_set, n_chunks){
  
  indices <- rep_len(1:n_chunks, floor(nrow(data_set)/n_chunks) + (nrow(data_set) %% n_chunks))
  
  lapply(1:n_chunks, function(x) data_set[indices == x, ])
  
  
  
  
}


update_beta <- function(penalty, pen_deriv, lambda, gamma, beta_mat, eta_avg){
  
  M <- ncol(beta_mat)
  
  beta_vec <- rowMeans(beta_mat)
  
  if(penalty == "lasso"){
    
   update1 <- beta_vec + eta_avg/gamma - lambda/(M*gamma) 
   update2 <- -beta_vec - eta_avg/gamma - lambda/(M*gamma)
   pmax(0, update1) - pmin(0, update2)
    
  } else if(penalty == "scad" | penalty == "mcp"){
      
      
      
    } else{
      
      stop("You've supplied an unsupported penalty function.")
      
    }
  
  
  
}

#' @export
main <- function(dat, M, intercept, maxiter, gamma, lambda, tau){
  
  # splitting dat into M chunks so that these can be iterated over
  # in the foreach loop
  dat_list <- split_data(data_set = dat, n_chunks = M)
  
  outcome_list <- lapply(dat_list, function(x) x[,1])
  
  designmat_list <- lapply(dat_list, function(x) x[, -1])
  
  dat_inverses <- lapply(dat_list, function(x) solve(crossprod(x)))
  

  iter <- 1
  beta_mat <- eta_mat <- matrix(0, nrow = p, ncol = M)
  beta_avg <- eta_avg <- rowMeans(beta_i)
  
  u_list <- r_list <- lapply(1:M, function(x) rep(0, nrow(designmat_list[[x]])))
  
  while(iter < maxiter){
  
   beta_global_i <- update_beta(penalty, pen_deriv, lambda, beta_mat, rowMeans(eta_mat))
   
          
   iter_run <- foreach(beta_i = itertools::isplitCols(beta_mat), eta_i = itertools::isplitCols(eta_mat), dat_i = dat_iter, outcome_i = outcome_list, design_i = designmat_list, u_i = u_list, r_i = r_list, .packages = "QRADMM" ) %dopar%{
     
     n_i <- nrow(design_i)
     ParamUpdates(betar = beta_i, etar = eta_i, xr = design_i, yr = outcome_i, dat = dat_i, ur = u_i, rr = r_i, beta = beta_global_i, alpha = alpha, tau = tau, n = n, ni = n_i )
     
     
   } 
    
   beta_mat <- sapply(iter_run, function(x) x$betai)
   
   eta_mat <- sapply(iter_run, function(x) x$etai)
   
   iter <- iter + 1
    
  }
  
  
  
}