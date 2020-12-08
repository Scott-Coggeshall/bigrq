



#' Fit Parallelized QRADMM
#' 
#' \code{r_main_parallel} is the model-fitting function for the parallelized QRADMM algorithm
#' implemented entirely in R.
#' 
#' @param dat a numeric matrix where the first column is the vector of responses and the rest of the columns
#' are the design matrix for the model.
#' @param M an integer specifying the number of data blocks.
#' @param intercept a logical indicating whether an intercept should be included in the model.
#' @param max_iter an integer, the maximum number of iterations to run.
#' @param min_iter an integer, the minimum number of iterations to run.
#' @param n_workers an integer, the number of worker processes to create.
#' @param lambda a numeric scalar or vector of non-zero values, the values for the tuning parameter.
#' @param tau a numeric scalar between 0 and 1, the conditional quantile.
#' @param rho a positive numeric scalar.
#' @param alpha a positive numeric scalar.
#' @param penalty a character string, the penalty to apply.
#' @param penalized a logical vector of the same length as the number of model parameters, 
#' TRUE values indicate that the corresponding parameter should be penalized.
#' @param inverses_write_path a character string, the location to save the calculated inverses to.
#' @param inverses_read_path a character string, the location to read the calculated inverses from.
#' @param abstol
#' @param reltol
#' 
#' @export
r_main_parallel <- function(dat, M, intercept, max_iter = 500, min_iter = 10, n_workers, lambda, tau, rho, alpha, penalty = "lasso", penalized , inverses_write_path = NULL, inverses_read_path = NULL, abstol = 1e-7, reltol = 1e-4){
  
  
  
  n_lambda <- length(lambda)
  n <- nrow(dat)
  p <- ncol(dat) - 1
  penalized_mat <- matrix(penalized, nrow = p, ncol = n_lambda)
  lambda_mat <- penalized_mat*rep(lambda, each = p)
  
  lambdan <- lambda_mat/n
  rhon <- rho/n
  indices <- rep(1:M, c(floor(n/M) + n%%M, rep(floor(n/M), M - 1)))  
  
  if(!is.null(inverses_read_path)){
    
    inverses <- readRDS(inverses_read_path)
    
    if(length(inverses) != M) stop("Number of cached inverses does not match number of data blocks M.")
  }
  
  
  beta_global_i <- matrix(0,nrow = p, ncol = n_lambda)
 
  # splitting data into M blocks
  dat_list <- suppressWarnings(split.data.frame(dat, indices))
  
  # chunking data
  dat_chunks <- suppressWarnings(split(dat_list, 1:n_workers))
  if(!is.null(inverses_read_path)) inverse_chunks <- suppressWarnings(split(inverses, 1:n_workers))
  ## initializing workers
  
  cl <- parallel::makeCluster(n_workers, setup_strategy = "sequential")
  on.exit(parallel::stopCluster(cl))
  ## exporting objects shared by all workers
  
  parallel::clusterExport(cl, varlist = c("beta_global_i", "lambdan", "alpha", "rhon", "shrink", "tau", "n", "compute_woodbury"), envir = environment())
  
  ## exporting chunked data to workers
  
  for(i in 1:n_workers){
    
    chunk_i <- dat_chunks[[i]]
    if(!is.null(inverses_read_path)){
    inverse_i <- inverse_chunks[[i]]
    parallel::clusterExport(cl[i], c("chunk_i", "inverse_i"), envir = environment())
    
    } else{
      
      parallel::clusterExport(cl[i], "chunk_i", envir = environment())
      
    }
  }
  if(!is.null(inverses_read_path)){
    rm(chunk_i, inverse_i)
    } else{
      
    rm(chunk_i)
        
    }
  ## initializing data containers on workers
  
  parallel::clusterEvalQ(cl, {
    
    p <- nrow(beta_global_i)
    
    n_lambda <- ncol(lambdan)
    param_list_i <- lapply(1:length(chunk_i), function(x) matrix(0, nrow = 2*p, ncol = n_lambda))
    
    if(!exists("inverse_i")){
    inverse_i <- lapply(chunk_i, function(x){
                        if(nrow(x[, -1]) < ncol(x[, -1])){
                                          
                            compute_woodbury(x[, -1])
                                          
                          } else{
                            solve(diag(1, nrow = ncol(x[, -1])) + crossprod(x[, -1]))
                                        }
                                        
                                      })
    
    }
    
    mat_i <- lapply(seq_along(inverse_i), function(i){
      
      tcrossprod(inverse_i[[i]], chunk_i[[i]][,-1])
      
    })
    
    u_i <- lapply(seq_along(chunk_i), function(i) matrix(0, nrow =  nrow(chunk_i[[i]]), ncol = n_lambda))
    r_i <- lapply(seq_along(chunk_i), function(i) matrix(chunk_i[[i]][, 1], nrow = length(chunk_i[[i]][, 1]), ncol = n_lambda))
    NULL
                            
  
    })

  
  ## entering while loop
  iter <- 1
  block_updates <- lapply(dat_chunks, function(x) lapply(1:length(x), function(y) matrix(0, nrow = 2*p, ncol = n_lambda)))
  
  
  
  while(iter <= max_iter){
  
    param_list <- unlist(block_updates, recursive = FALSE)
    beta_global_i <- update_beta(penalty = penalty, lambda = lambdan, rho = rhon, param_list = param_list)
    

    parallel::clusterExport(cl, "beta_global_i", envir = environment())
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
        
        
        xbeta <- alpha*chunk_i[[i]][, -1]%*%param_list_i[[i]][1:p,, drop = FALSE]   + (1 - alpha)*(chunk_i[[i]][, 1]- r_i[[i]])
        
        r_i[[i]] <- shrink(u_i[[i]]/rhon + chunk_i[[i]][, 1] - xbeta - .5*(2*tau - 1)/(n*rhon), .5*rep(1, length(chunk_i[[i]][, 1]))/(n*rhon))
        
       param_list_i[[i]][1:p,] <- mat_i[[i]]%*%(chunk_i[[i]][, 1] - r_i[[i]] + u_i[[i]]/rhon) - inverse_i[[i]]%*%param_list_i[[i]][(p + 1):(2*p),]/rhon + inverse_i[[i]]%*%beta_global_i
        
       u_i[[i]] <- u_i[[i]] + rhon*(chunk_i[[i]][, 1] - xbeta - r_i[[i]])
        
        param_list_i[[i]][(p+1):(2*p),] <- param_list_i[[i]][(p + 1):(2*p),] + rhon*(param_list_i[[i]][1:p,] - beta_global_i)
      }
      
      param_list_i
      
    })
    iter <- iter + 1
    print(iter)
  }

 if(!is.null(inverses_write_path)){ 
   
   inverses <- unlist(parallel::clusterEvalQ(cl, inverse_i), recursive = FALSE)
   saveRDS(inverses, file = inverses_write_path)
 
 }
 list(coef = beta_global_i, iter = iter)
  
}






#' Fit Parallelized QRADMM
#' 
#' \code{r_main_parallel_Matrix} is the model-fitting function for the parallelized QRADMM algorithm
#' implemented entirely in R.
#' 
#' @param dat a numeric matrix where the first column is the vector of responses and the rest of the columns
#' are the design matrix for the model.
#' @param M an integer specifying the number of data blocks.
#' @param intercept a logical indicating whether an intercept should be included in the model.
#' @param max_iter an integer, the maximum number of iterations to run.
#' @param min_iter an integer, the minimum number of iterations to run.
#' @param n_workers an integer, the number of worker processes to create.
#' @param lambda a numeric scalar or vector of non-zero values, the values for the tuning parameter.
#' @param tau a numeric scalar between 0 and 1, the conditional quantile.
#' @param rho a positive numeric scalar.
#' @param alpha a positive numeric scalar.
#' @param penalty a character string, the penalty to apply.
#' @param penalized a logical vector of the same length as the number of model parameters, 
#' TRUE values indicate that the corresponding parameter should be penalized.
#' @param inverses_write_path a character string, the location to save the calculated inverses to.
#' @param inverses_read_path a character string, the location to read the calculated inverses from.
#' @param abstol
#' @param reltol
#' 
#' @import Matrix
#' @importMethodsFrom Matrix %*% Arith Math
#' @importClassesFrom Matrix dgCMatrix dsCMatrix dgeMatrix
#' @export
r_main_parallel_Matrix <- function(X, y, M, intercept, max_iter = 500, min_iter = 10, n_workers, lambda, tau, rho, alpha, penalty = "lasso", penalized , inverses_write_path = NULL, inverses_read_path = NULL, abstol = 1e-7, reltol = 1e-4){
  
  
  
  n_lambda <- length(lambda)
  n <- X@Dim[1]
  p <- X@Dim[2]
  penalized_mat <- matrix(penalized, nrow = p, ncol = n_lambda)
  lambda_mat <- penalized_mat*rep(lambda, each = p)
  
  lambdan <- lambda_mat/n
  rhon <- rho/n
  indices <- rep(1:M, c(floor(n/M) + n%%M, rep(floor(n/M), M - 1)))  
  
  if(!is.null(inverses_read_path)){
    
    inverses <- readRDS(inverses_read_path)
    
    if(length(inverses) != M) stop("Number of cached inverses does not match number of data blocks M.")
  }
  
  
  beta_global_i <- Matrix::Matrix(0,nrow = p, ncol = n_lambda, sparse = TRUE)
  
  # splitting data into M blocks
  X_list <- suppressWarnings(split.data.frame(X, indices))
  y_list <- suppressWarnings(split(y, indices))
  # chunking data
  X_chunks <- suppressWarnings(split(X_list, 1:n_workers))
  y_chunks <- suppressWarnings(split(y_list, 1:n_workers))
  if(!is.null(inverses_read_path)) inverse_chunks <- suppressWarnings(split(inverses, 1:n_workers))
  ## initializing workers
  
  cl <- parallel::makeCluster(n_workers, setup_strategy = "sequential")
  on.exit(parallel::stopCluster(cl))
  ## exporting objects shared by all workers
  
  parallel::clusterExport(cl, varlist = c("beta_global_i", "lambdan", "alpha", "rhon", "shrink", "tau", "n", "compute_woodbury"), envir = environment())
  
  ## exporting chunked data to workers
  
  for(i in 1:n_workers){
    
    x_i <- X_chunks[[i]]
    y_i <- y_chunks[[i]]
    if(!is.null(inverses_read_path)){
      inverse_i <- inverse_chunks[[i]]
      parallel::clusterExport(cl[i], c("chunk_i", "inverse_i"), envir = environment())
      
    } else{
      
      parallel::clusterExport(cl[i], c("x_i", "y_i"), envir = environment())
      
    }
  }
  if(!is.null(inverses_read_path)){
    rm(chunk_i, inverse_i)
  } else{
    
    rm(x_i, y_i)
    
  }
  ## initializing data containers on workers
  
    parallel::clusterEvalQ(cl, {

    p <- beta_global_i@Dim[1]
    
    n_lambda <- ncol(lambdan)
    param_list_i <- lapply(1:length(x_i), function(x) Matrix::Matrix(0, nrow = 2*p, ncol = n_lambda, sparse = TRUE))
    
    if(!exists("inverse_i")){
      xx_i <- lapply(x_i, function(x){

           Matrix::Diagonal(p) + Matrix::crossprod(x)

        })

    }
    # 
    #  mat_i <- lapply(seq_along(inverse_i), function(i){
    #   
    #   tcrossprod(inverse_i[[i]], chunk_i[[i]][,-1])
    #   
    # })
    
    u_i <- lapply(seq_along(y_i), function(i) Matrix::Matrix(0, nrow =  length(y_i[[i]]), ncol = n_lambda))
    r_i <- lapply(seq_along(y_i), function(i) Matrix::Matrix(y_i[[i]], nrow = length(y_i[[i]]), ncol = n_lambda))
    NULL
    
  })
  
   
  
  ## entering while loop
  iter <- 1
  block_updates <- lapply(X_chunks, function(x) lapply(1:length(x), function(y) Matrix::Matrix(0, nrow = 2*p, ncol = n_lambda)))
  
  
  
  while(iter <= max_iter){
    
    param_list <- unlist(block_updates, recursive = FALSE)
    beta_global_i <- update_beta_Matrix(penalty = penalty, lambda = lambdan, rho = rhon, param_list = param_list)
    
    
    parallel::clusterExport(cl, "beta_global_i", envir = environment())
    ## global beta update beta_global_i <- update_beta
    block_updates <- parallel::clusterEvalQ(cl, {
      
      for(i in seq_along(x_i)){
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
        
        
        xbeta <- alpha*x_i[[i]]%*%param_list_i[[i]][1:p,, drop = FALSE]   + (1 - alpha)*(y_i[[i]]- r_i[[i]])

        r_i[[i]] <- shrink(u_i[[i]]/rhon + y_i[[i]] - xbeta - .5*(2*tau - 1)/(n*rhon), .5*rep(1, length(y_i[[i]]))/(n*rhon))

        #param_list_i[[i]][1:p,] <- mat_i[[i]]%*%(chunk_i[[i]][, 1] - r_i[[i]] + u_i[[i]]/rhon) - inverse_i[[i]]%*%param_list_i[[i]][(p + 1):(2*p),]/rhon + inverse_i[[i]]%*%beta_global_i

        param_list_i[[i]][1:p, ] <- Matrix::solve(xx_i[[i]], Matrix::crossprod(x_i[[i]], y_i[[i]] - r_i[[i]] + u_i[[i]]/rhon) - param_list_i[[i]][(p + 1):(2*p), ]/rhon + beta_global_i)
      # 
        u_i[[i]] <- u_i[[i]] + rhon*(y_i[[i]] - xbeta - r_i[[i]])
      # 
       param_list_i[[i]][(p+1):(2*p),] <- param_list_i[[i]][(p + 1):(2*p),] + rhon*(param_list_i[[i]][1:p,] - beta_global_i)
       }
      # 
      param_list_i
      # 
    })
    iter <- iter + 1
    print(iter)
  }
 
  if(!is.null(inverses_write_path)){ 
    
    inverses <- unlist(parallel::clusterEvalQ(cl, inverse_i), recursive = FALSE)
    saveRDS(inverses, file = inverses_write_path)
    
  }
  list(coef = beta_global_i, iter = iter)
  
}




