bigrq <- function(x, y, M = 1, intercept = TRUE, max_iter = 500, min_iter = 10, 
                  n_workers = floor(parallel::detectCores()/2), 
                  lambda = 1, tau=0.5, rho = 10, alpha = 10, penalty = "lasso", inverses_read_path = "file", 
                  inverses_write_path = "file", abstol = 10, reltol = 10){
  
  ### Error checking block
  if(class(try(as.matrix(x), silent = TRUE))[1] == "try-error") {
    stop("x could not be converted to matrix.")
  }
  if(class(try(as.matrix(y), silent = TRUE))[1] == "try-error") {
    stop("y could not be converted to matrix.")
  }
  if(nrow(x) != nrow(y)) {stop("Number of rows in x and y do not match.")}
  if(!is.numeric(M)) {stop("M is not an integer or numeric value.")}
  if(!is.logical(intercept)) {stop("intercept is not a logical value.")}
  if(!is.numeric(max_iter)) {stop("max_iter is not an integer or numeric value.")}
  if(!is.numeric(min_iter)) {stop("min_iter is not an integer or numeric value.")}
  if(!is.numeric(M)) {stop("M is not an integer or numeric value.")}
  if(!is.numeric(n_workers)) {stop("n_workers is not an integer or numeric value.")}
  if(n_workers > parallel::detectCores()) stop("n_workers is greater than available cores")
  if(!is.numeric(lambda)) {stop("lambda is not an integer or numeric value or vector.")}
  if(tau < 0 | tau > 1) {stop("tau must be between 0 and 1")}
  if(!penalty %in% c("lasso", "scad", "mcap")) {stop('penalty must be one of "lasso", "scad", or "mcap"')}
  
  for (i in c("abstol", "reltol", "rho", "alpha")) {
    if(!(is.numeric(get(i)[1]) & get(i)[1] >= 0)) {
      stop(paste0(i, " must be positive numeric scalar"))
    }
  }
  
  ### Intermediate variable block
  
  for (i in c("M", "intercept", "max_iter", "min_iter", "n_workers", "tau", "rho", "alpha", "penalty", "abstol", "reltol")) {
    if (length(get(i)) > 1) {
      assign(i, get(i)[1])
      warning(paste0(i, " has length greater than 1, only using first element and ignoring the rest."))
    }
  }
  
  if(is.matrix(x) == FALSE) {
    x <- as.matrix(x)
  }
  if(intercept == TRUE) {
    int_col <- matrix(data = 1, ncol = 1, nrow = nrow(x))
    x <- cbind(int_col, x)
  }
  dat <- cbind(y, x)
  
  ### Model fitting block - call r_main_parallel here
  # currently won't work due to inverse file paths not existing?
  # model <- r_main_parallel(dat = dat, M = M, intercept = intercept, max_iter = max_iter, min_iter = min_iter, 
  #                          n_workers = n_workers, lambda = lambda, tau = tau, rho = rho, alpha = alpha, penalty = penalty, 
  #                          inverses_read_path = inverses_read_path, inverses_write_path = inverses_write_path, 
  #                          abstol = 10, reltol = 10)
  
  ### Output processing block
  return(alpha)
}
