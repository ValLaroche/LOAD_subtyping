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

CC_modelization <- function(datasets, title_model, cluster_alg = "hc", distance = "euclidean"){
  CC_model = ConsensusClusterPlus(datasets, maxK = 6, reps = 50, pItem = 0.8, pFeature = 1,
                                  clusterAlg = cluster_alg, distance = distance, 
                                  innerLinkage="ward.D2", finalLinkage="ward.D2",
                                  seed = 126, plot = "pdf", title = title_model)
  return(CC_model)
}

diablo_model <- function(data_matrix, metadata, dirname){
  
  metadata = metadata[metadata$HC_clusters != "Control",]
  
  diablo_Methyls = list()
  diablo_Methyls[[1]] = mixOmics::splsda(t(data_matrix), metadata$HC_clusters, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
  diablo_Methyls[[2]] = mixOmics::splsda(t(data_matrix), metadata$KM_clusters, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
  diablo_Methyls[[3]] = mixOmics::splsda(t(data_matrix), metadata$PAM_clusters, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
  diablo_Methyls[[4]] = mixOmics::splsda(t(data_matrix), metadata$Spec_clusters, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
  # 
  # for(dMethyl in seq(1:length(diablo_Methyls))){
  #   if(dMethyl == 1){
  #     outdir = paste0("results/", dirname, "/HClust/")
  #   } else if (dMethyl == 2){
  #     outdir = paste0("results/", dirname, "/Kmeans/")
  #   } else if (dMethyl == 3){
  #     outdir = paste0("results/", dirname, "/PAM/")
  #   } else if (dMethyl == 4){
  #     outdir = paste0("results/", dirname, "/Spectrum/")
  #   }
  #   
  #   plot_mixOmics_all(model = diablo_Methyls[[dMethyl]], outdir = outdir)
  # }
  
  return(diablo_Methyls)
  
}


cv_final_model <- function(data_matrix_AD, metadata_AD, top_features, dirname){
  
  dir.create(paste0("results/", dirname))
  dir.create(paste0("results/", dirname, "/Spectrum"))
  dir.create(paste0("results/", dirname, "/HClust"))
  dir.create(paste0("results/", dirname, "/Kmeans"))
  dir.create(paste0("results/", dirname, "/PAM"))
  
  
  HC_features_all = top_features[[1]]
  KM_features_all = top_features[[2]]
  PAM_features_all = top_features[[3]]
  Spec_features_all = top_features[[4]]
  
  # Extract intersect of all CV models (Top specific features)
  HC_final_features = Reduce(intersect, list(HC_features_all[,1], HC_features_all[,2], HC_features_all[,3], HC_features_all[,4],
                                             HC_features_all[,5], HC_features_all[,6], HC_features_all[,7], HC_features_all[,8], 
                                             HC_features_all[,9], HC_features_all[,10]))
  KM_final_features = Reduce(intersect, list(KM_features_all[,1], KM_features_all[,2], KM_features_all[,3], KM_features_all[,4],
                                             KM_features_all[,5], KM_features_all[,6], KM_features_all[,7], KM_features_all[,8], 
                                             KM_features_all[,9], KM_features_all[,10]))
  PAM_final_features = Reduce(intersect, list(PAM_features_all[,1], PAM_features_all[,2], PAM_features_all[,3], PAM_features_all[,4], 
                                              PAM_features_all[,5], PAM_features_all[,6], PAM_features_all[,7], PAM_features_all[,8], 
                                              PAM_features_all[,9], PAM_features_all[,10]))
  Spec_final_features = Reduce(intersect, list(Spec_features_all[,1], Spec_features_all[,2], Spec_features_all[,3], Spec_features_all[,4], 
                                               Spec_features_all[,5], Spec_features_all[,6], Spec_features_all[,7], Spec_features_all[,8], 
                                               Spec_features_all[,9], Spec_features_all[,10]))
  
  # Build new betas matrices for each algorithm
  HC_betas = data_matrix_AD[rownames(data_matrix_AD) %in% HC_final_features,]
  KM_betas = data_matrix_AD[rownames(data_matrix_AD) %in% KM_final_features,]
  PAM_betas = data_matrix_AD[rownames(data_matrix_AD) %in% PAM_final_features,]
  Spec_betas = data_matrix_AD[rownames(data_matrix_AD) %in% Spec_final_features,]
  # 200 features per comp for consistency or lower in case less features are selected
  HC_comp = min(200, length(HC_final_features))
  KM_comp = min(200, length(KM_final_features))
  PAM_comp = min(200, length(PAM_final_features))
  Spec_comp = min(200, length(Spec_final_features))
  
  #Build final crossvalidated models
  diablo_final = list()
  diablo_final[[1]] = mixOmics::splsda(t(HC_betas), metadata_AD$HC_clusters, keepX=c(HC_comp,HC_comp,HC_comp,HC_comp,HC_comp,HC_comp), ncomp = 6)
  diablo_final[[2]] = mixOmics::splsda(t(KM_betas), metadata_AD$KM_clusters, keepX=c(KM_comp,KM_comp,KM_comp,KM_comp,KM_comp,KM_comp), ncomp = 6)
  diablo_final[[3]] = mixOmics::splsda(t(PAM_betas), metadata_AD$PAM_clusters, keepX=c(PAM_comp,PAM_comp,PAM_comp,PAM_comp,PAM_comp,PAM_comp), ncomp = 6)
  diablo_final[[4]] = mixOmics::splsda(t(Spec_betas), metadata_AD$Spec_clusters, keepX=c(Spec_comp,Spec_comp,Spec_comp,Spec_comp,Spec_comp,Spec_comp), ncomp = 6)
  #Generate figures
  for(dMethyl in seq(1:length(diablo_final))){
    if(dMethyl == 1){
      outdir = paste0("results/",dirname,"/HClust/")
    } else if (dMethyl == 2){
      outdir = paste0("results/",dirname,"/Kmeans/")
    } else if (dMethyl == 3){
      outdir = paste0("results/",dirname,"/PAM/")
    } else if (dMethyl == 4){
      outdir = paste0("results/",dirname,"/Spectrum/")
    }
    
    plot_mixOmics_all(model = diablo_final[[dMethyl]], outdir = outdir)
  }
  
  return(list(diablo_final, HC_final_features, KM_final_features, PAM_final_features, Spec_final_features))
}
