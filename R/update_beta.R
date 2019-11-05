
update_beta <- function(penalty, pen_deriv, lambda, rho, beta_mat, eta_avg){
  
  M <- ncol(beta_mat)
  
  beta_vec <- rowMeans(beta_mat)
  
  if(penalty == "lasso"){
    
    shrink(beta_vec + eta_avg/rho, lambda/(M*rho))
    
  } else if(penalty == "scad"){
    
    scad_deriv <- lambda*((abs(beta_vec) <= lambda)*1 + (rho*lambda - beta_vec)/((rho - 1)*lambda)*(abs(beta_vec) <=  rho*lambda))
    
    shrink(beta_vec + eta_avg/rho, scad_deriv/(M*rho))
    } else if(penalty == "mcp"){
      mcp_deriv <- (lambda - abs(beta_vec)/rho)*(abs(beta_vec) <= rho*lambda)/rho
      
      shrink(beta_vec + eta_avg/rho, mcp_deriv/(M*rho))
  } else{ 
     stop("You've supplied an unsupported penalty function.")
    
  }
  
  
}

