
check_convergence <- function(dat, rho, beta_vec, beta_old_vec, r_vec, u_vec, abstol, reltol){
  
  n <- nrow(dat)
  xbeta <- dat[, -1]%*%beta_vec
  xbeta_old <- dat[, -1]%*%beta_old_vec
  
  rnorm <- sqrt(sum((dat[, 1] - xbeta - r_vec)^2))
  
  snorm <- sqrt(sum((rho*(xbeta - xbeta_old))^2))
  
  compare_vals <- c(sqrt(sum((xbeta)^2)), sqrt(sum(r_vec^2)))
  
  epspri <- sqrt(n)*abstol + reltol*max(compare_vals)
  epsdual <- sqrt(n)*abstol + reltol*sqrt(sum(u_vec^2))
  
  print(paste("uvec = ", sqrt(sum(u_vec^2))))
  print(paste("rnorm = ", rnorm))
  print(paste("snorm = ", snorm))
  print(paste("epspri = ", epspri))
  print(paste("epsdual = ", epsdual))
  print(paste("beta diff = ", sqrt(sum((beta_vec-beta_old_vec)^2))))
  
  rnorm < epspri && snorm < epsdual
}

check_convergence_standard <- function(beta_old, beta_new, eps){
  
  
  sqrt(sum((beta_old - beta_new)^2)) < eps
  
  
  
  
}
