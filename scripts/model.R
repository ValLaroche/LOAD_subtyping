aff_matrix_calc <- function(dataset){
  dist_matrix = as.matrix(dist(t(dataset)))
  
  aff_matrix = affinity_matrix(dist_matrix, k = 20)
  
  return(aff_matrix)
}

Spectrum_modelization <- function(list_datasets, spectrum_method = 3, cluster_alg = "GMM", num_clusters = 2){
  
  spec_model = Spectrum(list_datasets, method = spectrum_method, silent = FALSE,
                 showres = TRUE, diffusion = TRUE,
                 kerneltype = c("stsc"),
                 maxk = 10, NN = 10, NN2 = 10, showpca = T,
                 frac = 2, thresh = 7, fontsize = 18, dotsize = 3,
                 tunekernel = FALSE, clusteralg = cluster_alg, FASP = FALSE, FASPk = NULL, fixk = num_clusters, KNNs_p = 15)
    
  names(spec_model$assignments) = colnames(spec_model$similarity_matrix)
  
  return(spec_model)
}

SNF_modelization <- function(datasets_sim_matrix, K = 20, t = 20, num_clusters = 4){
  snf_model_matrix = SNF(datasets_sim_matrix, K, t)
  
  snf_model_clusters = spectralClustering(snf_model_matrix, num_clusters)
   # integrate_similarity_matrices(list(W1,W2,W3), KNNs_p = 10,diffusion_iters = 4, method = "TPG")
  
  snf_model_clusters = as.numeric(snf_model_clusters)
  
  names(snf_model_clusters) = colnames(snf_model_matrix)
  
  snf_model = list(snf_model_matrix, snf_model_clusters)
  
  return(snf_model)
}

ANF_modelization <- function(list_aff_matrix, K = 20, num_clusters = 2){
  
  list_ANF_model = list()
  
  for(i in seq(2:length(num_clusters))){
    anf_model_matrix = ANF(Wall = list_aff_matrix, K = K, verbose = T)
    
    anf_model_clusters = spectral_clustering(anf_model_matrix, num_clusters[i])
    
    names(anf_model_clusters) = colnames(anf_model_matrix)
    
    anf_model = list(anf_model_matrix, anf_model_clusters)
  
    list_ANF_model[[length(list_ANF_model) + 1]] = anf_model
  }
  
  return(list_ANF_model)
}

CC_modelization <- function(datasets, title_model, cluster_alg = "hc", distance = "pearson"){
  CC_model = ConsensusClusterPlus(datasets, maxK = 6, reps = 50, pItem = 0.8, pFeature = 1,
                                  clusterAlg = cluster_alg, distance = distance, 
                                  seed = 126, plot = "png", title = title_model)
  return(CC_model)
}

JIVE_modelization <- function(datasets){
  
  jive_model = jive(datasets)
  
  return(jive_model)
  
}

DIABLO_modelization <- function(datasets, Y_classes, keep_vals){
  
  for (i in seq(1:length(datasets))){
      datasets[[i]] = t(as.matrix(datasets[[i]]))
  }
  
  multi_DIABLO = block.splsda(X = datasets,
                              Y = Y_classes,
                              keepX = keep_vals)
}

ICB_modelization <- function(datasets, data_distrib, K){
  
  list_datasets = list(NULL, NULL, NULL, NULL, NULL)
  
  for (i in seq(1:length(datasets))){
    if(!is.null(datasets[[i]])){
      list_datasets[[i]] = t(as.matrix(datasets[[i]]))
    }
  }
  
  multi_iCB = iClusterBayes(dt1 = list_datasets[[1]],
                            dt2 = list_datasets[[2]],
                            dt3 = list_datasets[[3]],
                            dt4 = list_datasets[[4]],
                            dt5 = list_datasets[[5]],
                            type = data_distrib,
                            K = K)
  return(multi_iCB)
}

SNF.CC_modelization <- function(datasets, num_clusters){
  multi_SNF.CC = ExecuteSNF.CC(datasets=datasets, clusterNum=num_clusters,
                               K=20, alpha=0.5, t=20,maxK = 10, pItem = 0.8,reps=500, 
                               title = "GBM", plot = "png", finalLinkage ="average")    

  return(multi_SNF.CC)
}
