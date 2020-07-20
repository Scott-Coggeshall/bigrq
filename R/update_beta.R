
update_beta <- function(penalty, pen_deriv, lambda, rho, beta_list, eta_list){
 
  # convert lists to 3d arrays for easy application of rowMeans
  M <- length(beta_list)
  
  beta_array <- sapply(1:M, function(x) beta_list[[x]], simplify = "array")
  eta_array <- sapply(1:M, function(x) eta_list[[x]], simplify = 'array')
  
  
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

