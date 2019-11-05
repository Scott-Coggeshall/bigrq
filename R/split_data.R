split_data <- function(data_set, n_chunks){
  
  indices <- rep_len(1:n_chunks, floor(nrow(data_set)/n_chunks) + (nrow(data_set) %% n_chunks))
  
  lapply(1:n_chunks, function(x) data_set[indices == x, ])
  
  
  
  
}