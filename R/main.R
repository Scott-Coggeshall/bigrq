





#' @export
main <- function(dat, M, intercept, maxiter, lambda, tau, rho, alpha, abstol = 1e-7, reltol = 1e-4){
  
  # splitting dat into M chunks so that these can be iterated over
  # in the foreach loop
  dat_list <- split_data(data_set = dat, n_chunks = M)
  
  outcome_list <- lapply(dat_list, function(x) x[,1])
  
  designmat_list <- lapply(dat_list, function(x) x[, -1])
  
  dat_inverses <- lapply(dat_list, function(x) solve(crossprod(x[, -1]) + diag(1, nrow = ncol(x[, -1]))))
  

  n <- nrow(dat)
  p <- ncol(dat) - 1

  beta_global_i <- rep(0, p)
  beta_mat <- eta_mat <- matrix(0, nrow = p, ncol = M)
  beta_avg <- eta_avg <- rowMeans(beta_mat)
  
  u_list <- lapply(1:M, function(x) rep(0, nrow(designmat_list[[x]])))
  
  r_list <- outcome_list
  

  
  iter <- 1
  keep_going <- T
  
  while(keep_going && iter < maxiter){
   beta_old <- beta_global_i
   beta_global_i <- update_beta(penalty, pen_deriv, lambda/n, rho/n, beta_mat, rowMeans(eta_mat))
  
          
   iter_run <- foreach(beta_i = itertools::isplitCols(beta_mat, chunks = M), eta_i = itertools::isplitCols(eta_mat, chunks = M), dat_i = dat_inverses, outcome_i = outcome_list, design_i = designmat_list, u_i = u_list, r_i = r_list, .packages = "QRADMM", .noexport = "dat" ) %dopar%{
     
     n_i <- nrow(design_i)
     ParamUpdates(betar = beta_i, etar = eta_i, xr = design_i, yr = outcome_i, dat = dat_i, ur = u_i, rr = r_i, beta = beta_global_i, rho = rho, alpha = alpha, tau = tau, n = n, ni = n_i )
     
     
   } 
    
   beta_mat <- sapply(iter_run, function(x) x$betai)
   
   eta_mat <- sapply(iter_run, function(x) x$etai)
   
   u_list <- lapply(iter_run, function(x) x$ui)
   
   r_list <- lapply(iter_run, function(x) x$ri)
   
   keep_going <- !check_convergence(dat, rho/n, as.vector(beta_global_i), as.vector(beta_old), unlist(r_list), unlist(u_list), abstol, reltol)
   iter <- iter + 1
  
  }
  
  print(iter)
  beta_global_i
  
  
}