
update_beta <- function(penalty, pen_deriv, lambda, rho, beta_mat, eta_avg){
  
  M <- ncol(beta_mat)
  
  beta_vec <- rowMeans(beta_mat)
  
  if(penalty == "lasso"){
    
    shrink(beta_vec + eta_avg/rho, lambda/(M*rho))
    
  } else if(penalty == "scad" | penalty == "mcp"){
    
    
    
  } else{
    
    stop("You've supplied an unsupported penalty function.")
    
  }
  
  
  
}

