normalize_data <- function(expr_matrix, 
                           norm_method, 
                           is_RNA1,
                           threshold_read_counts,
                           threshold_perc_samples){

  expr_matrix = expr_matrix[apply(expr_matrix,1,function(x) 
    sum(x < threshold_read_counts)) <= ncol(expr_matrix)*(1-threshold_perc_samples),]

  if(norm_method == "log2"){
    expr_matrix = log2(expr_matrix + 1)
  }
    
  else if(norm_method == "SNF"){
    expr_matrix = standardNormalization(expr_matrix)
  }
  
  data.checkDistribution(expr_matrix)  
  return(as.data.frame(expr_matrix))
}

pca_analysis <- function(expr_matrix){
  expr_matrix_PCA = prcomp(t(expr_matrix))
  fviz_eig(expr_matrix_PCA)
  fviz_pca_ind(expr_matrix_PCA,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE)   # Avoid text overlapping
}
