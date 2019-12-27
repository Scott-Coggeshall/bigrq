
update_beta <- function(penalty, pen_deriv, lambda, rho, beta_mat, eta_mat){
  
  M <- ncol(beta_mat)
  
  beta_avg <- rowMeans(beta_mat)
  eta_avg <- rowMeans(eta_mat)
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

