

#' @export

r_main <- function(dat, M, intercept, p_linear = NULL, maxiter = 500, miniter = 10, lambda, tau, rho, alpha, penalty = "lasso", parallel = FALSE, abstol = 1e-7, reltol = 1e-4){


  # splitting dat into M chunks so that these can be iterated over
  # in the foreach loop
  dat_list <- split_data(data_set = dat, n_chunks = M)
 
  outcome_list <- lapply(dat_list, function(x) x[,1])

  designmat_list <- lapply(dat_list, function(x) x[, -1])

  dat_inverses <- lapply(dat_list, function(x) solve(crossprod(x[, -1]) + diag(1, nrow = ncol(x[, -1]))))


  n <- nrow(dat)
  
  if(is.null(p_linear)){
    
    
    p <- p_linear <- ncol(dat) - 1
    p_nonlinear <- 0
    
  } else{
    
    p <- ncol(dat) - (1 + p_linear)
    
    p_nonlinear <- p - p_linear
  
  }
  
  lambda <- c(rep(lambda, p_linear), rep(0, p_nonlinear))

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

r_main_parallel <- function(dat, M, intercept, maxiter = 500, miniter = 10, n_workers, lambda, tau, rho, alpha, penalty = "lasso", abstol = 1e-7, reltol = 1e-4){
  
  
  
  n_lambda <- length(lambda)
  print("check")
  n <- nrow(dat)
  p <- ncol(dat) - 1
  
  indices <- rep(1:M, c(floor(n/M) + n%%M, rep(floor(n/M), M - 1)))  
  
  
  beta_global_i <- matrix(0,nrow = n_lambda, ncol = p)
  beta_mat <- eta_mat  <- lapply(1:M, function(x) matrix(0, nrow = n_lambda, ncol =  p))
  
  u_list <- r_list <- resids_list <-  lapply(1:M, function(x) matrix(0, nrow = n_lambda, ncol = sum(indices == x)))
  
  # splitting data into M blocks
  dat_list <- split.data.frame(dat[, -1], indices)
  
  # splitting outcome into M blocks
  outcome_list <- lapply(split.data.frame(dat[,1, drop = F], indices), as.vector)
  # storing inverses
  if(n/M >= p){
  dat_inverses <-  foreach(dat_i = dat_list)%dopar%{
    solve(crossprod(dat_i) + diag(1, nrow = p))
  }
  
  } else{
    
  dat_inverses <- foreach(dat_i = dat_list)%dopar%{
    
    I_n <- diag(nrow(dat_i))
    I_p <- diag(ncol(dat_i))
    XXt <- tcrossprod(dat_i)
    
    B_inv <- solve(I_n + XXt)
    
   I_p - t(dat_i)%*%B_inv%*%dat_i
    
  }  
    
  }
  
  print('did inverses')
  
  
  iter <- 1
  
  
  rhon <- rho/n
  lambdan <- lambda/n
  
  keep_going <- T
  
  # first iteration
  beta_old <- beta_global_i
  
  while(keep_going){
  beta_global_i <- update_beta(penalty, pen_deriv, lambda/n, rho/n, beta_mat, eta_mat)
  
 
  
  block_update <- foreach(beta_mat_i = beta_mat, eta_mat_i = eta_mat, dat_i = dat_list, outcome_i = outcome_list,
          inverse_i = dat_inverses, u_i = u_list, r_i = r_list, resids_i = resids_list, .export = "shrink")%dopar%{
            
            # iterate over lambda vals
            for(lambda_i in seq_along(lambdan)){
              
            xbeta <- alpha*dat_i%*%beta_mat_i[lambda_i,] + (1 - alpha)*(outcome_i - r_i[lambda_i, ])
            
            r <- shrink(u_i[lambda_i, ]/rhon + outcome_i - xbeta - .5*(2*tau - 1)/(n*rhon), .5*rep(1, length(outcome_i))/(n*rhon))
            
            r_i[lambda_i, ] <- as.vector(r)
            
            beta_mat_i[lambda_i, ] <- inverse_i%*%(t(dat_i)%*%(outcome_i - r_i[lambda_i, ] + u_i[lambda_i, ]/rhon) - eta_mat_i[lambda_i, ]/rhon + beta_global_i[lambda_i, ] )
            
            u_i[lambda_i, ] <- as.vector(u_i[lambda_i, ] + rhon*(outcome_i - xbeta - r_i[lambda_i, ]))
            
            eta_mat_i[lambda_i, ] <- eta_mat_i[lambda_i, ] + rhon*(beta_mat_i[lambda_i,] - beta_global_i[lambda_i, ])
            
            resids_i[lambda_i, ] <- outcome_i - dat_i%*%beta_mat_i[lambda_i,] - r_i[lambda_i, ]
            
            }
            
            list(beta_mat_i = beta_mat_i, eta_mat_i = eta_mat_i, u_i = u_i, r_i = r_i, resids_i = resids_i )
          }
   
    beta_mat <- lapply(block_update, function(x) x[[1]])
    eta_mat <- lapply(block_update, function(x) x[[2]])
    u_list <- lapply(block_update, function(x) x[[3]])
    r_list <- lapply(block_update, function(x) x[[4]])
    resids_list <- lapply(block_update, function(x) x[[5]])
    
    print(iter)
    iter <- iter + 1
    keep_going <- (iter <= maxiter)
    
  }
  
  beta_global_i
}
