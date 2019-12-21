split_data <- function(data_set, n_chunks){
  
  n <- nrow(data_set)
  
  indices <- rep(1:n_chunks, c(floor(n/n_chunks) + n%%n_chunks, rep(floor(n/n_chunks), n_chunks - 1)))  
  lapply(1:n_chunks, function(x) data_set[indices == x, ])
  
  
  
  
}
