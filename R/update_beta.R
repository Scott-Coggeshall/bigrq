
update_beta <- function(penalty, pen_deriv = NULL, lambda, rho, param_list){
 
  # convert lists to 3d arrays for easy application of rowMeans
  M <- length(param_list)
  
  p <- dim(param_list[[1]])[1]/2
  
  beta_array <- sapply(1:M, function(x) param_list[[x]][1:p,, drop = FALSE], simplify = "array")
  eta_array <- sapply(1:M, function(x) param_list[[x]][(p + 1):(2*p), , drop = FALSE], simplify = 'array')
  
  
  beta_avg <- rowMeans(beta_array, dims = 2)
  eta_avg <- rowMeans(eta_array, dims = 2)
  if(penalty == "lasso"){
    
    shrink(beta_avg + eta_avg/rho, lambda/(M*rho))
    
  } else if(penalty == "scad"){
    
    scad_deriv <- lambda*((abs(beta_avg) <= lambda)*1 + (rho*lambda - beta_avg)/((rho - 1)*lambda)*(abs(beta_avg) <=  rho*lambda))
    
    shrink(beta_avg + eta_avg/rho, scad_deriv/(M*rho))
    } else if(penalty == "mcp"){
      mcp_deriv <- (lambda - abs(beta_avg)/rho)*(abs(beta_avg) <= rho*lambda)/rho
      
      shrink(beta_avg + eta_avg/rho, mcp_deriv/(M*rho))
  } else{ 
     stop("You've supplied an unsupported penalty function.")
    
  }
  
  
}

update_beta_Matrix <- function(penalty, pen_deriv = NULL, lambda, rho, param_list){
  
  # convert lists to 3d arrays for easy application of rowMeans
  M <- length(param_list)
  
  p <- dim(param_list[[1]])[1]/2
  
  beta_array <- sapply(1:M, function(x) as.matrix(param_list[[x]][1:p,, drop = FALSE]), simplify = "array")
  eta_array <- sapply(1:M, function(x) as.matrix(param_list[[x]][(p + 1):(2*p), , drop = FALSE]), simplify = 'array')
  
  
  beta_avg <- rowMeans(beta_array, dims = 2)
  eta_avg <- rowMeans(eta_array, dims = 2)
  if(penalty == "lasso"){
    
    shrink(beta_avg + eta_avg/rho, lambda/(M*rho))
    
  } else if(penalty == "scad"){
    
    scad_deriv <- lambda*((abs(beta_avg) <= lambda)*1 + (rho*lambda - beta_avg)/((rho - 1)*lambda)*(abs(beta_avg) <=  rho*lambda))
    
    shrink(beta_avg + eta_avg/rho, scad_deriv/(M*rho))
  } else if(penalty == "mcp"){
    mcp_deriv <- (lambda - abs(beta_avg)/rho)*(abs(beta_avg) <= rho*lambda)/rho
    
    shrink(beta_avg + eta_avg/rho, mcp_deriv/(M*rho))
  } else{ 
    stop("You've supplied an unsupported penalty function.")
    
  }
  
  
}


