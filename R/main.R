

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

r_main_parallel <- function(dat, M, intercept, max_iter = 500, min_iter = 10, n_workers, lambda, tau, rho, alpha, penalty = "lasso", abstol = 1e-7, reltol = 1e-4){
  
  
  
  n_lambda <- length(lambda)
  n <- nrow(dat)
  p <- ncol(dat) - 1
  
  indices <- rep(1:M, c(floor(n/M) + n%%M, rep(floor(n/M), M - 1)))  
  
  
  beta_global_i <- matrix(0,nrow = p, ncol = n_lambda)
 
  # splitting data into M blocks
  dat_list <- split.data.frame(dat[, -1], indices)
  
  # splitting outcome into M blocks
  
  
  ## initializing workers
  
  cl <- parallel::makeCluster(n_workers, setup_strategy = "sequential")
  on.exit(stopCluster(cl))
  ## exporting objects shared by all workers
  
  parallel::clusterExport(cl, varlist = c("beta_global_i", "lambdan"))
  
  ## exporting chunked data to workers
  
  for(i in 1:n_workers){
    
    chunk_i <- dat_chunks[[i]]
    
    parallel::clusterExport(cl[i], "chunk_i")
    
  }
  
  ## initializing data containers on workers
  
  parallel::clusterEvalQ(cl, {
    
    p <- ncol(beta_global_i)
    
    n_lambda <- length(lambdan)
    param_list_i <- lapply(1:length(chunk_i), matrix(0, dim = c(2*p, n_lambda)))
    inverse_i <- lapply(chunk_i, function(x){
                        if(nrow(x[, -1]) < ncol(x[, -1])){
                                          
                            bigrq::compute_woodbury(x[, -1])
                                          
                          } else{
                            solve(diag(1, nrow = ncol(x[, -1])) + crossprod(x[, -1]))
                                        }
                                        
                                      })
    u_i <- lapply(seq_along(chunk_i), function(i) matrix(0, nrow = n_lambda, ncol = nrow(chunk_i[[i]])))
    r_i <- lapply(seq_along(chunk_i), function(i) matrix(chunk_i[[i]][, 1], nrow = n_lambda, ncol = length(chunk_i[[i]][, 1])), byrow = T)
    NULL
                            
  
    })

  
  ## entering while loop
  iter <- 1
  while(iter <= max_iter){
    
    param_list <- unlist(block_updates, recursive = FALSE)
    beta_global_i <- bigrq::update_beta(penalty = penalty, lambda = lambdan, rho = rhon, param_list = param_list)
    ## global beta update beta_global_i <- update_beta
    block_updates <- parallel::clusterEvalQ(cl, {
      
      for(i in seq_along(inverse_i)){
        # iterate over lambda vals
        # for(lambda_i in seq_along(lambdan)){
        #   
        #   xbeta <- alpha*chunk_i[[i]][, -1]%*%beta_mat_i[[i]][lambda_i,] + (1 - alpha)*(chunk_i[[i]][, 1] - r_i[[i]][lambda_i, ])
        #   
        #   
        #   r_i[[i]][lambda_i, ] <- shrink(u_i[[i]][lambda_i, ]/rhon + chunk_i[[i]][, 1] - xbeta - .5*(2*tau - 1)/(n*rhon), .5*rep(1, length(chunk_i[[i]][, 1]))/(n*rhon))
        #   
        #   
        #   beta_mat_i[[i]][lambda_i, ] <- inverse_i[[i]]%*%(t(chunk_i[[i]][, -1])%*%(chunk_i[[i]][, 1] - r_i[[i]][lambda_i, ] + u_i[[i]][lambda_i, ]/rhon) - eta_mat_i[[i]][lambda_i, ]/rhon + beta_global_i[lambda_i, ] )
        #   
        #   u_i[[i]][lambda_i, ] <- as.vector(u_i[[i]][lambda_i, ] + rhon*(chunk_i[[i]][, 1] - xbeta - r_i[[i]][lambda_i, ]))
        #   
        #   eta_mat_i[[i]][lambda_i, ] <- eta_mat_i[[i]][lambda_i, ] + rhon*(beta_mat_i[[i]][lambda_i,] - beta_global_i[lambda_i, ])
        #   
        #   
        #   
        # }
        ## xbeta is now a matrix
        xbeta <- alpha*chunk_i[[i]][, -1]%*%param_mat_i[1:p,,i] + (1 - alpha)*(chunk_i[[i]][, 1] - r_i[[i]])
        
        r_i[[i]] <- shrink(u_i[[i]]/rhon + chunk_i[[i]][, 1] - xbeta - .5*(2*tau - 1)/(n*rhon), .5*rep(1, length(chunk_i[[i]][, 1]))/(n*rhon))
        
        param_list_i[[i]][1:p,] <- inverse_i[[i]]%*%(t(chunk_i[[i]][, -1])%*%(chunk_i[[i]][, 1] - r_i[[i]] + u_i[[i]]/rhon) - param_list_i[[i]][(p + 1):(2*p),]/rhon + beta_global_i)
        
        u_i[[i]] <- u_i[[i]] + rhon*(chunk_i[[i]][, 1] - xbeta - r_i[[i]])
        
        param_list_i[[i]][(p+1):(2*p),] <- param_list_i[[i]][(p + 1):(2*p),] + rhon*(param_list_i[1:p,] - beta_global_i)
      }
      
      param_list_i
      
    })
    
    
    
  }
  
 beta_global_i
  
}





