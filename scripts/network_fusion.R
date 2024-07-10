library(SNFtool)
library(ANF)

# list_aff_matrix is a list of each single-omic matrix, calculated with aff_matrix_calc
aff_matrix_calc <- function(dataset){
  dist_matrix = as.matrix(dist(t(dataset)))
  
  aff_matrix = affinity_matrix(dist_matrix, k = 20) # k number of nearest neighbors
  
  return(aff_matrix)
}

# K in SNF and ANF is the number of K-nearest neighbors selected by the algorithm
# t is the number of iterations for the diffusion process
# Both libraries have their associated spectral clustering functions to work with the outcome matrix
# of the associated algorithm
SNF_modelization <- function(list_aff_matrix, K = 20, t = 20, num_clusters){
  snf_model_matrix = SNF(Wall = list_aff_matrix, K, t)
  
  snf_model_clusters = spectralClustering(snf_model_matrix, K = num_clusters)
  
  names(snf_model_clusters) = colnames(snf_model_matrix)
  
  snf_model = list(snf_model_matrix, snf_model_clusters)
  
  return(snf_model)
}

ANF_modelization <- function(list_aff_matrix, K = 20, num_clusters){
  
  anf_model_matrix = ANF(Wall = list_aff_matrix, K = K, verbose = T)
  
  anf_model_clusters = spectral_clustering(anf_model_matrix, k = num_clusters)
  
  names(anf_model_clusters) = colnames(anf_model_matrix)
  
  anf_model = list(anf_model_matrix, anf_model_clusters)
  
  return(anf_model)
}
