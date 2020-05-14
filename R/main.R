





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
   
   if(iter > 1) keep_going <- !check_convergence(dat, rho/n, as.vector(beta_global_i), as.vector(beta_old), unlist(r_list), unlist(u_list), abstol, reltol)
   iter <- iter + 1
  
  }
  
  print(iter)
  beta_global_i
  
  
}



#' @export
r_main <- function(dat, M, intercept, maxiter = 500, miniter = 10, lambda, tau, rho, alpha, penalty = "lasso", abstol = 1e-7, reltol = 1e-4){

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
  
  resids_list <- u_list
  


  iter <- 1


  rhon <- rho/n
  lambdan <- lambda/n
  
  keep_going <- T

  while( iter < maxiter & keep_going){
    beta_old <- beta_global_i
    beta_global_i <- update_beta(penalty, pen_deriv, lambda/n, rho/n, beta_mat, eta_mat)


    
    for(i in 1:M){

      xbeta <- alpha*designmat_list[[i]]%*%beta_mat[,i] + (1 - alpha)*(outcome_list[[i]] - r_list[[i]])

      r <- shrink(u_list[[i]]/rhon + outcome_list[[i]] - xbeta - .5*(2*tau - 1)/(n*rhon), .5*rep(1, length(outcome_list[[i]]))/(n*rhon))

      r_list[[i]] <- as.vector(r)

      beta_mat[, i] <- dat_inverses[[i]]%*%(t(designmat_list[[i]])%*%(outcome_list[[i]] - r + u_list[[i]]/rhon) - eta_mat[, i]/rhon + beta_global_i )

      u_list[[i]] <- as.vector(u_list[[i]] + rhon*(outcome_list[[i]] - xbeta - r_list[[i]]))

      eta_mat[, i] <- eta_mat[, i] + rhon*(beta_mat[,i] - beta_global_i)
      
      resids_list[[i]] <- outcome_list[[i]] - designmat_list[[i]]%*%beta_mat[,i] - r_list[[i]]

    }



    if(iter > miniter) keep_going <- !check_convergence_standard(beta_old, beta_global_i, abstol)
    iter <- iter + 1

    }
  
    


  list(iter = iter, beta = beta_global_i, r = r_list, u = u_list, resids = resids_list)



}

r_main_parallel <- function(dat, M, intercept, maxiter = 500, miniter = 10, lambda, tau, rho, alpha, penalty = "lasso", abstol = 1e-7, reltol = 1e-4){
  
  
  indices <- rep(1:M, c(floor(n/M) + n%%M, rep(floor(n/M), M - 1)))  
  
  
  n <- nrow(dat)
  p <- ncol(dat) - 1
  
  beta_global_i <- rep(0, p)
  beta_mat <- eta_mat <- matrix(0, nrow = p, ncol = M)
  beta_avg <- eta_avg <- rowMeans(beta_mat)
  
  
  u_mat <- r_mat <- resids_mat <-  matrix(0, nrow = floor(n/M) + n%%M, ncol = M)
  
  # splitting data into M blocks
  dat_list <- split.data.frame(dat[, -1], indices)
  
  # splitting outcome into M blocks
  outcome_list <- lapply(split.data.frame(dat[,1, drop = F], indices), as.vector)
  # storing inverses
  dat_inverses <-  lapply(dat_list, function(x) solve(crossprod(x) + diag(1, nrow = p)))
  
  
  
  
  iter <- 1
  
  
  rhon <- rho/n
  lambdan <- lambda/n
  
  keep_going <- T
  
  # first iteration
  beta_old <- beta_global_i
  beta_global_i <- update_beta(penalty, pen_deriv, lambda/n, rho/n, beta_mat, eta_mat)
  
  
  
  foreach(beta_mat_i = itertools::isplitCols(beta_mat ), eta_mat_i = itertools::isplitCols(eta_mat),
         , itertools::isplitCols(dat_inverses))%dopar%{
            
            
            xbeta <- alpha*designmat_list[[i]]%*%beta_mat[,i] + (1 - alpha)*(outcome_list[[i]] - r_list[[i]])
            
            r <- shrink(u_list[[i]]/rhon + outcome_list[[i]] - xbeta - .5*(2*tau - 1)/(n*rhon), .5*rep(1, length(outcome_list[[i]]))/(n*rhon))
            
            r_list[[i]] <- as.vector(r)
            
            beta_mat[, i] <- dat_inverses[[i]]%*%(t(designmat_list[[i]])%*%(outcome_list[[i]] - r + u_list[[i]]/rhon) - eta_mat[, i]/rhon + beta_global_i )
            
            u_list[[i]] <- as.vector(u_list[[i]] + rhon*(outcome_list[[i]] - xbeta - r_list[[i]]))
            
            eta_mat[, i] <- eta_mat[, i] + rhon*(beta_mat[,i] - beta_global_i)
            
            resids_list[[i]] <- outcome_list[[i]] - designmat_list[[i]]%*%beta_mat[,i] - r_list[[i]]
            
            
            
          }
  
  
}
