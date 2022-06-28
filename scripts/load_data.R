load_data <- function(path_to_matrix){
  expr_matrix = read.delim(path_to_matrix, sep = "\t", header = TRUE)
  
  colnames(expr_matrix) = gsub("X", "", colnames(expr_matrix))
    
  return(expr_matrix)
}
