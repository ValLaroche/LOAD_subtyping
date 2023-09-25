##### Sourcing scripts, libraries and data #####
setwd("D:/valentin/main/")

sources = dir("./scripts/",full.names=TRUE)
for(s in sources){
  source(s)
}

load("eQTLDB/BDR_final_set_20-9-23-50perc.Rdata")
pheno_BDR = BDR_ok_ADCtrl
load("eQTLDB/UKBBN_final_set_20-9-23-50perc.Rdata")
pheno_UKBBN = pheno
load("eQTLDB/PITT_final_set_20-9-23-50perc.Rdata")
pheno_PITT = pheno
remove(BDR_ok_ADCtrl, pheno)

##### Preparing data #####

cpgs_BDR = rownames(betas_BDR)
cpgs_UKBBN = rownames(betas_UKBBN)
cpgs_PITT = rownames(betas_PITT)
cpgs_intersect = Reduce(intersect, list(cpgs_BDR, cpgs_PITT, cpgs_UKBBN))

output = make_dataset(data_matrix = betas_BDR, metadata = pheno_BDR, dirname = "BDR", db = "BDR", cpgs_intersect = cpgs_intersect)  
data_matrix_AD = output[[1]]
metadata_AD = output[[2]]
remove(output)

##### Make models #####

#Preview to select number of clusters
df_assignment = wrapper_models(data_matrix_AD, metadata_AD,cluster_vector = c(), dirname = "BDR")

#Optimal number of clusters
df_assignment = wrapper_models(data_matrix_AD, metadata_AD,cluster_vector = c(3,3,3,3), dirname = "BDR")

  
##### SPLSDA output of the identified subtypes

  metadata_AD$HC_clusters = df_assignment
  metadata_AD$KM_clusters = df_assignment
  metadata_AD$PAM_clusters = df_assignment
  metadata_AD$Spec_clusters = df_assignment
  
  
  diablo_Methyls = list()
  diablo_Methyls[[1]] = mixOmics::splsda(t(data_matrix_AD), metadata_AD$HC_clusters, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
  diablo_Methyls[[2]] = mixOmics::splsda(t(data_matrix_AD), metadata_AD$KM_clusters, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
  diablo_Methyls[[3]] = mixOmics::splsda(t(data_matrix_AD), metadata_AD$PAM_clusters, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
  diablo_Methyls[[4]] = mixOmics::splsda(t(data_matrix_AD), metadata_AD$Spec_clusters, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
  
  for(dMethyl in seq(1:length(diablo_Methyls))){
    if(dMethyl == 1){
      outdir = paste0("results/", dirname, "/HClust/")
    } else if (dMethyl == 2){
      outdir = paste0("results/", dirname, "/Kmeans/")
    } else if (dMethyl == 3){
      outdir = paste0("results/", dirname, "/PAM/")
    } else if (dMethyl == 4){
      outdir = paste0("results/", dirname, "/Spectrum/")
    }
    
    plot_mixOmics_all(model = diablo_Methyls[[dMethyl]], outdir = outdir)
  }

##### Random permutations for specificity #####

  df_shuffle = random_permutations(data_matrix_AD, metadata_AD, diablo_Methyls)

  colnames(df_shuffle) = c("HC_accuracy","HC_intersect",
                           "KM_accuracy","KM_intersect",
                           "PAM_accuracy","PAM_intersect",
                           "Spec_accuracy","Spec_intersect")
  
  df_shuffle = df_shuffle[-1,]

  output_permutation(df_shuffle, dirname)
    
##### 10-fold validation #####
  #10-fold CV study of models 
  top_features = cv_diablo(data_matrix_AD, metadata_AD)
  
  #Selecting shared features across CV
  final_model_features = cv_final_model(data_matrix_AD, metadata_AD, top_features, dirname = "final")
  diablo_final = final_model_features[[1]]
  HC_final_features = final_model_features[[2]]
  KM_final_features = final_model_features[[3]]
  PAM_final_features = final_model_features[[4]]
  Spec_final_features = final_model_features[[5]]
  
##### Cross-method clustering #####
  
  #Building cross-method labeling based on 4 outcomes
  output = build_crossmethod(metadata_AD, df_assignment)
  metadata_AD = output[[1]]
  df_mat_assignment_order = output[[2]]
  groups = output[[3]]
  snames_mat = output[[4]]
  
  gc()
  ##### BDR cross-method #####
  #Based on the 3-4 shared outcomes, build new cross-method clusters
  table(metadata_AD$all_clust)[order(table(metadata_AD$all_clust), decreasing = T)]
  table(metadata_AD$groups)[order(table(metadata_AD$groups), decreasing = T)]
  groups$Clusters = NA
  groups[groups$groups %in% c("Group 2","Group 7"),]$Clusters = "Cluster 1"
  groups[groups$groups %in% c("Group 1","Group 9"),]$Clusters = "Cluster 2"
  groups[groups$groups %in% c("Group 3","Group 4", "Group 6"),]$Clusters = "Cluster 3"
  groups[groups$groups %in% c("Group 5","Group 11", "Group 14"),]$Clusters = "Cluster 4"
  groups[groups$groups %in% c("Group 12","Group 13","0","Group 10","Group 8"),]$Clusters= "Excluded"
  groups$groups = NULL
  groups$Clusters = factor(groups$Clusters, levels = c("Cluster 1", "Cluster 2","Cluster 3","Cluster 4","Excluded"))
  color_list_clust = list(Clusters = c("Cluster 1" = "chocolate2",
                                       "Cluster 2" = "blue2",
                                       "Cluster 3" = "brown2",
                                       "Cluster 4" = "chartreuse2",
                                       "Excluded" = "darkgrey"))
  color_list_mat = colorRampPalette(c("white", "firebrick"))(5)
  
  diag(df_mat_assignment_order) = NA
  df_mat_assignment_order = matrix(as.numeric(df_mat_assignment_order),    # Convert to numeric matrix
                                   ncol = ncol(df_mat_assignment_order))
  rownames(df_mat_assignment_order) = snames_mat
  colnames(df_mat_assignment_order) = snames_mat
  
  pdf(paste0("./results/BDR/Crossmethod.pdf"), width  = 10, height = 10)
  pheatmap(df_mat_assignment_order[order(groups$Clusters),order(groups$Clusters)],
           col = color_list_mat, labels_row = rep("", 313), labels_col = rep("", 313),
           annotation_col = groups, cluster_rows = F, cluster_cols = F, annotation_colors = color_list_clust)
  dev.off()
  
  groups$Basename = rownames(groups)
  groups = groups[order(groups$Basename),]

  ##### Crossmethod mixomics #####
  #Build the SPLSDA model to extract top features of crossmethod
  cv_res =  cv_crossmethod(data_matrix_AD, metadata_AD, groups, dirname = "BDR", dirname_final = "final")
  
  Cross_methyl = cv_res[[1]]
  Cross_final = cv_res[[2]]
  Cross_final_features = cv_res[[3]]
  cross_metadata = cv_res[[4]]
##### Prediction to UKBBN #####
  
  top_cpg = list()
  res_pred_list = list()
  
  #For all models, make a prediction of the outcome on the UKBBN cohort
  #Extract the top features and outcomes
  project_clust = projection_to_rep(model = diablo_Methyls[[1]], data_matrix_AD, rep_betas = msetEPIC_betas, rep_metadata = pheno, nclusters = 3, discovery = "BDR", replication = "UKBBN", outdir = "BDR/HClust/")
  top_cpg[[1]] = project_clust[[1]]
  res_pred_list[[1]] = project_clust[[2]]
  project_clust = projection_to_rep(model = diablo_Methyls[[2]], data_matrix_AD, rep_betas = msetEPIC_betas, rep_metadata = pheno, nclusters = 3, discovery = "BDR", replication = "UKBBN", outdir = "BDR/KMeans/")
  top_cpg[[2]] = project_clust[[1]]
  res_pred_list[[2]] = project_clust[[2]]
  project_clust = projection_to_rep(model = diablo_Methyls[[3]], data_matrix_AD, rep_betas = msetEPIC_betas, rep_metadata = pheno, nclusters = 3, discovery = "BDR", replication = "UKBBN", outdir = "BDR/PAM/")
  top_cpg[[3]] = project_clust[[1]]
  res_pred_list[[3]] = project_clust[[2]]
  project_clust = projection_to_rep(model = diablo_Methyls[[4]], data_matrix_AD, rep_betas = msetEPIC_betas, rep_metadata = pheno, nclusters = 3, discovery = "BDR", replication = "UKBBN", outdir = "BDR/Spectrum/")
  top_cpg[[4]] = project_clust[[1]]
  res_pred_list[[4]] = project_clust[[2]]
  project_clust = projection_to_rep(model = Cross_methyl, data_matrix_AD, rep_betas = msetEPIC_betas, rep_metadata = pheno, nclusters = 4, discovery = "BDR", replication = "UKBBN", outdir = "BDR/cross/")
  top_cpg[[5]] = project_clust[[1]]
  res_pred_list[[5]] = project_clust[[2]]
  
  #For all CV models, make a prediction of the outcome on the UKBBN cohort
  #Extract the top features and outcomes
  filt_betas = data_matrix_AD[rownames(data_matrix_AD) %in% HC_final_features,]
  project_clust = projection_to_rep(diablo_final[[1]], data_matrix_AD = filt_betas, rep_betas = msetEPIC_betas, rep_metadata = pheno, nclusters = 3, discovery = "BDR", replication = "UKBBN", outdir = "final/HClust/")
  top_cpg[[6]] = project_clust[[1]]
  res_pred_list[[6]] = project_clust[[2]]
  
  filt_betas = data_matrix_AD[rownames(data_matrix_AD) %in% KM_final_features,]
  project_clust = projection_to_rep(diablo_final[[2]], data_matrix_AD = filt_betas, rep_betas = msetEPIC_betas, rep_metadata = pheno, nclusters = 3, discovery = "BDR", replication = "UKBBN", outdir = "final/KMeans/")
  top_cpg[[7]] = project_clust[[1]]
  res_pred_list[[7]] = project_clust[[2]]
  
  filt_betas = data_matrix_AD[rownames(data_matrix_AD) %in% PAM_final_features,]
  project_clust = projection_to_rep(diablo_final[[3]], data_matrix_AD = filt_betas, rep_betas = msetEPIC_betas, rep_metadata = pheno, nclusters = 3, discovery = "BDR", replication = "UKBBN", outdir = "final/PAM/")
  top_cpg[[8]] = project_clust[[1]]
  res_pred_list[[8]] = project_clust[[2]]
  
  filt_betas = data_matrix_AD[rownames(data_matrix_AD) %in% Spec_final_features,]
  project_clust = projection_to_rep(diablo_final[[4]], data_matrix_AD = filt_betas, rep_betas = msetEPIC_betas, rep_metadata = pheno, nclusters = 3, discovery = "BDR", replication = "UKBBN", outdir = "final/Spectrum/")
  top_cpg[[9]] = project_clust[[1]]
  res_pred_list[[9]] = project_clust[[2]]
  
  filt_betas = data_matrix_AD[rownames(data_matrix_AD) %in% Cross_final_features,]
  project_clust = projection_to_rep(model = Cross_final, data_matrix_AD = filt_betas, rep_betas = msetEPIC_betas, rep_metadata = pheno, nclusters = 4, discovery = "BDR", replication = "UKBBN", outdir = "final/cross/")
  top_cpg[[10]] = project_clust[[1]]
  res_pred_list[[10]] = project_clust[[2]]

  UKBBN_topcpg = top_cpg
  UKBBN_reslist = res_pred_list
  
##### Prediction to NIH #####
  
  top_cpg = list()
  res_pred_list = list()
  #For all models, make a prediction of the outcome on the PITT cohort
  #Extract the top features and outcomes
  
  project_clust = projection_to_rep(model = diablo_Methyls[[1]], data_matrix_AD, rep_betas = NIH_betas, rep_metadata = NIH_pheno, nclusters = 3, discovery = "BDR", replication = "NIH", outdir = "BDR/HClust/")
  top_cpg[[1]] = project_clust[[1]]
  res_pred_list[[1]] = project_clust[[2]]
  project_clust = projection_to_rep(model = diablo_Methyls[[2]], data_matrix_AD, rep_betas = NIH_betas, rep_metadata = NIH_pheno, nclusters = 3, discovery = "BDR", replication = "NIH", outdir = "BDR/KMeans/")
  top_cpg[[2]] = project_clust[[1]]
  res_pred_list[[2]] = project_clust[[2]]
  project_clust = projection_to_rep(model = diablo_Methyls[[3]], data_matrix_AD, rep_betas = NIH_betas, rep_metadata = NIH_pheno, nclusters = 3, discovery = "BDR", replication = "NIH", outdir = "BDR/PAM/")
  top_cpg[[3]] = project_clust[[1]]
  res_pred_list[[3]] = project_clust[[2]]
  project_clust = projection_to_rep(model = diablo_Methyls[[4]], data_matrix_AD, rep_betas = NIH_betas, rep_metadata = NIH_pheno, nclusters = 3, discovery = "BDR", replication = "NIH", outdir = "BDR/Spectrum/")
  top_cpg[[4]] = project_clust[[1]]
  res_pred_list[[4]] = project_clust[[2]]
  project_clust = projection_to_rep(model = Cross_methyl, data_matrix_AD, rep_betas = NIH_betas, rep_metadata = NIH_pheno, nclusters = 4, discovery = "BDR", replication = "NIH", outdir = "BDR/cross/")
  top_cpg[[5]] = project_clust[[1]]
  res_pred_list[[5]] = project_clust[[2]]
  #For all CV models, make a prediction of the outcome on the PITT cohort
  #Extract the top features and outcomes
  filt_betas = data_matrix_AD[rownames(data_matrix_AD) %in% HC_final_features,]
  project_clust = projection_to_rep(diablo_final[[1]], data_matrix_AD = filt_betas, rep_betas = NIH_betas, rep_metadata = NIH_pheno, nclusters = 3, discovery = "BDR", replication = "NIH", outdir = "final/HClust/")
  top_cpg[[6]] = project_clust[[1]]
  res_pred_list[[6]] = project_clust[[2]]
  
  filt_betas = data_matrix_AD[rownames(data_matrix_AD) %in% KM_final_features,]
  project_clust = projection_to_rep(diablo_final[[2]], data_matrix_AD = filt_betas, rep_betas = NIH_betas, rep_metadata = NIH_pheno, nclusters = 3, discovery = "BDR", replication = "NIH", outdir = "final/KMeans/")
  top_cpg[[7]] = project_clust[[1]]
  res_pred_list[[7]] = project_clust[[2]]
  
  filt_betas = data_matrix_AD[rownames(data_matrix_AD) %in% PAM_final_features,]
  project_clust = projection_to_rep(diablo_final[[3]], data_matrix_AD = filt_betas, rep_betas = NIH_betas, rep_metadata = NIH_pheno, nclusters = 3, discovery = "BDR", replication = "NIH", outdir = "final/PAM/")
  top_cpg[[8]] = project_clust[[1]]
  res_pred_list[[8]] = project_clust[[2]]
  
  filt_betas = data_matrix_AD[rownames(data_matrix_AD) %in% Spec_final_features,]
  project_clust = projection_to_rep(diablo_final[[4]], data_matrix_AD = filt_betas, rep_betas = NIH_betas, rep_metadata = NIH_pheno, nclusters = 3, discovery = "BDR", replication = "NIH", outdir = "final/Spectrum/")
  top_cpg[[9]] = project_clust[[1]]
  res_pred_list[[9]] = project_clust[[2]]
  
  filt_betas = data_matrix_AD[rownames(data_matrix_AD) %in% Cross_final_features,]
  project_clust = projection_to_rep(model = Cross_final, data_matrix_AD = filt_betas, rep_betas = NIH_betas, rep_metadata = NIH_pheno, nclusters = 4, discovery = "BDR", replication = "NIH", outdir = "final/cross/")
  top_cpg[[10]] = project_clust[[1]]
  res_pred_list[[10]] = project_clust[[2]]
  
  NIH_topcpg = top_cpg
  NIH_reslist = res_pred_list
  
  ##### all CPG stats #####
  
  BDR_metadata_AD$cross_clusters = NA
  BDR_metadata_AD[BDR_metadata_AD$Basename %in% cross_metadata$Basename,]$cross_clusters = cross_metadata$Cluster
  #Generate the EWAS for all models in BDR
  df_res_lm_BDR = cpg_stats_BDR(metadata_AD = BDR_metadata_AD, metadata_ctrl = BDR_Ctrl, mSet_betas = mSet_betas, top_cpg)

  UKBBN_AD = pheno[pheno$AD == 1,]
  UKBBN_AD$HC_clusters = UKBBN_reslist[[1]]$class
  UKBBN_AD$KM_clusters = UKBBN_reslist[[2]]$class
  UKBBN_AD$PAM_clusters = UKBBN_reslist[[3]]$class
  UKBBN_AD$Spec_clusters = UKBBN_reslist[[4]]$class
  UKBBN_AD$cross_clusters = UKBBN_reslist[[5]]$class
  UKBBN_AD$CVHC_clusters = UKBBN_reslist[[6]]$class
  UKBBN_AD$CVKM_clusters = UKBBN_reslist[[7]]$class
  UKBBN_AD$CVPAM_clusters = UKBBN_reslist[[8]]$class
  UKBBN_AD$CVSpec_clusters = UKBBN_reslist[[9]]$class
  UKBBN_AD$CVcross_clusters = UKBBN_reslist[[10]]$class
  UKBBN_Ctrl = pheno[pheno$AD == 0,]
  
  #Generate the EWAS for all models in UKBBN
  df_res_lm_UKBBN = cpg_stats_UKBBN(metadata_AD = UKBBN_AD, metadata_ctrl = UKBBN_Ctrl, mSet_betas = msetEPIC_betas, UKBBN_topcpg)
  
  NIH_AD = NIH_pheno[NIH_pheno$Phenotype %in% c("AD+P","AD-P"),]
  NIH_AD$HC_clusters = NIH_reslist[[1]]$class
  NIH_AD$KM_clusters = NIH_reslist[[2]]$class
  NIH_AD$PAM_clusters = NIH_reslist[[3]]$class
  NIH_AD$Spec_clusters = NIH_reslist[[4]]$class
  NIH_AD$cross_clusters = NIH_reslist[[5]]$class
  NIH_AD$CVHC_clusters = NIH_reslist[[6]]$class
  NIH_AD$CVKM_clusters = NIH_reslist[[7]]$class
  NIH_AD$CVPAM_clusters = NIH_reslist[[8]]$class
  NIH_AD$CVSpec_clusters = NIH_reslist[[9]]$class
  NIH_AD$CVcross_clusters = NIH_reslist[[10]]$class
  NIH_Ctrl = NIH_pheno[NIH_pheno$Phenotype == "C",]
  
  #Generate the EWAS for all models in PITT
  df_res_lm_NIH = cpg_stats_UKBBN(metadata_AD = NIH_AD, metadata_ctrl = NIH_Ctrl, mSet_betas = NIH_betas, NIH_topcpg)
  
  length(unique(df_res_lm_BDR$cpg))
  length(unique(df_res_lm_UKBBN$cpg))
  length(unique(df_res_lm_NIH$cpg))

  ###########
  #Check cpg consistency
  length(intersect(df_res_lm_NIH$cpg, intersect(df_res_lm_BDR$cpg, df_res_lm_UKBBN$cpg)))
  intersect_cohorts = intersect(df_res_lm_NIH$cpg, intersect(df_res_lm_BDR$cpg, df_res_lm_UKBBN$cpg))
  
  df_res_lm_BDR = df_res_lm_BDR[df_res_lm_BDR$cpg %in% intersect_cohorts,]
  df_res_lm_UKBBN = df_res_lm_UKBBN[df_res_lm_UKBBN$cpg %in% intersect_cohorts,]
  df_res_lm_NIH = df_res_lm_NIH[df_res_lm_NIH$cpg %in% intersect_cohorts,]
  
  #Complete EWAS table
  df_res_final = merge(df_res_lm_BDR, df_res_lm_UKBBN, by = c("model","cluster","cpg"))
  df_res_final = merge(df_res_final, df_res_lm_NIH, by = c("model","cluster","cpg"))
  colnames(df_res_final) = c("model","cluster","cpg","BDR_pvalue","BDR_rsquared","BDR_effectsize","BDR_cohend",
                             "UKBBN_pvalue","UKBBN_rsquared","UKBBN_effectsize","UKBBN_cohend",
                             "NIH_pvalue","NIH_rsquared","NIH_effectsize","NIH_cohend")
  length(unique(df_res_final$cpg))
  

  df_res_final$BDR_pvalue = as.numeric(df_res_final$BDR_pvalue)
  df_res_final$UKBBN_pvalue = as.numeric(df_res_final$UKBBN_pvalue)
  df_res_final$NIH_pvalue = as.numeric(df_res_final$NIH_pvalue)
  
  df_res_final$BDR_rsquared = as.numeric(df_res_final$BDR_rsquared)
  df_res_final$UKBBN_rsquared = as.numeric(df_res_final$UKBBN_rsquared)
  df_res_final$NIH_rsquared = as.numeric(df_res_final$NIH_rsquared)
  
  df_res_final$BDR_effectsize = as.numeric(df_res_final$BDR_effectsize)
  df_res_final$UKBBN_effectsize = as.numeric(df_res_final$UKBBN_effectsize)
  df_res_final$NIH_effectsize = as.numeric(df_res_final$NIH_effectsize)
  
  df_res_final$BDR_cohend = as.numeric(df_res_final$BDR_cohend)
  df_res_final$UKBBN_cohend = as.numeric(df_res_final$UKBBN_cohend)
  df_res_final$NIH_cohend = as.numeric(df_res_final$NIH_cohend)
  
  length(df_res_final$BDR_pvalue[df_res_final$BDR_pvalue < 0.05])
  
  #### Hypergeo test #### 
  
  df_res_effsize = df_res_final
   # Are we still doing this since not reliable ?
  ##### UKBBN #####
  dir.create("./results/jaccard_UKBBN")
  
  hyper_UKBBN = list()
  hyper_idx = 1
  
  for(mod in unique(df_res_lm_BDR$model)){
    list_cpgs_resBDR = list()
    list_cpgs_resUKBBN = list()
    for(nclust in c("1", "3")){
      
      print(paste0(mod, " - ", nclust, "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"))
      
      tmp_both_sig = df_res_effsize[df_res_effsize$model == mod & df_res_effsize$cluster == nclust,]
      max_value = nrow(tmp_both_sig)
      
      tmp_both_sig$eff_sizeok = FALSE
      
      tmp_both_sig[tmp_both_sig$BDR_effectsize > 0 & tmp_both_sig$UKBBN_effectsize > 0,]$eff_sizeok = TRUE
      tmp_both_sig[tmp_both_sig$BDR_effectsize < 0 & tmp_both_sig$UKBBN_effectsize < 0,]$eff_sizeok = TRUE
      tmp_both_sig = tmp_both_sig[tmp_both_sig$eff_sizeok == TRUE,]
      
      tmp_res_BDR = tmp_both_sig[tmp_both_sig$BDR_pvalue < 0.05 & tmp_both_sig$BDR_cohend >= 0.2,]
      tmp_res_UKBBN = tmp_both_sig[tmp_both_sig$UKBBN_pvalue < 0.05 & tmp_both_sig$UKBBN_cohend >= 0.2,]
      
      list_cpgs_resBDR[[paste0("BDR;",mod,";",nclust)]] = tmp_res_BDR$cpg
      list_cpgs_resUKBBN[[paste0("UKBBN;",mod,";",nclust)]] = tmp_res_UKBBN$cpg
      
      
    }
    
    GOM_res = newGOM(list_cpgs_resBDR, list_cpgs_resUKBBN,
                     max_value)
    GOM_matrix = getMatrix(GOM_res, name="pval")
    GOM_melt = melt(GOM_matrix)
    
    GOM_inter = getNestedList(GOM_res, name="intersection")
    
    
    GOM_melt$jacc = NA
    
    
    GOM_index = 1
    for(UKBBNi in seq(1:length(list_cpgs_resUKBBN))){
      for(BDRi in seq(1:length(list_cpgs_resBDR))){
        
        BDRcpg = list_cpgs_resBDR[[BDRi]]
        UKBBNcpg = list_cpgs_resUKBBN[[UKBBNi]]
        
        jaccI = length(intersect(BDRcpg, UKBBNcpg)) / length(union(BDRcpg, UKBBNcpg))
        
        GOM_melt[GOM_index,]$jacc = jaccI
        GOM_index = GOM_index + 1
      }
    }
    
    GOM_melt$value =  as.numeric(format(GOM_melt$value, scientific = TRUE, digits = 2))
    GOM_melt[GOM_melt$value > 0.05,]$value = "N.S."
    colnames(GOM_melt) = c("BDR", "UKBBN", "value", "Jaccard_index")
    
    GOM_plot = ggplot(GOM_melt, aes(x = BDR, y = UKBBN, fill = Jaccard_index))+
      ggtitle("Subtype specifity - BDR/UKBBN cohort",
              subtitle = mod) +
      geom_tile() +
      geom_text(aes(label = value)) +
      scale_fill_gradient2(low = "white", high = "chocolate2")+
      theme_minimal() + xlab("") + ylab("") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    pdf(paste0("results/jaccard_UKBBN/",mod,"_Jaccard.pdf"), width = 7, height = 7)
    plot(GOM_plot)
    dev.off()

    hyper_UKBBN[[hyper_idx]] = list(GOM_res, GOM_inter, GOM_melt)
    hyper_idx = hyper_idx + 1 
  }
  
  ##### NIH #####
  dir.create("./results/jaccard_NIH")
  
  hyper_NIH = list()
  hyper_idx = 1
  
  for(mod in unique(df_res_lm_BDR$model)){
    list_cpgs_resBDR = list()
    list_cpgs_resNIH = list()
    for(nclust in c("1","2","3")){
      
      print(paste0(mod, " - ", nclust, "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"))
      
      tmp_both_sig = df_res_effsize[df_res_effsize$model == mod & df_res_effsize$cluster == nclust,]
      max_value = nrow(tmp_both_sig)
      
      tmp_both_sig$eff_sizeok = FALSE
      
      tmp_both_sig[tmp_both_sig$BDR_effectsize > 0 & tmp_both_sig$NIH_effectsize > 0,]$eff_sizeok = TRUE
      tmp_both_sig[tmp_both_sig$BDR_effectsize < 0 & tmp_both_sig$NIH_effectsize < 0,]$eff_sizeok = TRUE
      tmp_both_sig = tmp_both_sig[tmp_both_sig$eff_sizeok == TRUE,]
      
      tmp_res_BDR = tmp_both_sig[tmp_both_sig$BDR_pvalue < 0.001 & tmp_both_sig$BDR_cohend >= 0.2,]
      tmp_res_NIH = tmp_both_sig[tmp_both_sig$NIH_pvalue < 0.001 & tmp_both_sig$NIH_cohend >= 0.2,]
      
      print(nrow(tmp_res_NIH))
      GO_res = newGeneOverlap(tmp_res_BDR$cpg,
                              tmp_res_NIH$cpg,
                              genome.size = max_value)
      GO_res = testGeneOverlap(GO_res)
      print(GO_res)
      
      phyper(length(intersect(tmp_res_BDR$cpg, tmp_res_NIH$cpg))-1, length(tmp_res_BDR$cpg), 
             max_value - length(tmp_res_BDR$cpg), length(tmp_res_NIH$cpg), lower.tail = F, log.p = F)
      
      list_cpgs_resBDR[[paste0("BDR;",mod,";",nclust)]] = tmp_res_BDR$cpg
      list_cpgs_resNIH[[paste0("PITT;",mod,";",nclust)]] = tmp_res_NIH$cpg
      
      
    }
    
    GOM_res = newGOM(list_cpgs_resBDR, list_cpgs_resNIH,
                     max_value)
    GOM_matrix = getMatrix(GOM_res, name="pval")
    GOM_melt = melt(GOM_matrix)
    
    GOM_inter = getNestedList(GOM_res, name="intersection")
    
    GOM_melt$jacc = NA
    
    
    GOM_index = 1
    for(NIHi in seq(1:length(list_cpgs_resNIH))){
      for(BDRi in seq(1:length(list_cpgs_resBDR))){
        
        BDRcpg = list_cpgs_resBDR[[BDRi]]
        NIHcpg = list_cpgs_resNIH[[NIHi]]
        
        jaccI = length(intersect(BDRcpg, NIHcpg)) / length(union(BDRcpg, NIHcpg))
        
        GOM_melt[GOM_index,]$jacc = jaccI
        GOM_index = GOM_index + 1
      }
    }
    
    GOM_melt$value =  as.numeric(format(GOM_melt$value, scientific = TRUE, digits = 2))
    GOM_melt[GOM_melt$value > 0.05,]$value = "N.S."
    colnames(GOM_melt) = c("BDR", "PITT", "value", "Jaccard_index")
    
    GOM_plot = ggplot(GOM_melt, aes(x = BDR, y = PITT, fill = Jaccard_index))+
      ggtitle("Subtype specifity - BDR/PITTADR cohort",
              subtitle = mod) +
      geom_tile() +
      geom_text(aes(label = value)) +
      scale_fill_gradient2(low = "white", high = "chocolate2")+
      theme_minimal() + xlab("") + ylab("") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    pdf(paste0("results/jaccard_NIH/",mod,"_Jaccard.pdf"), width = 7, height = 7)
    plot(GOM_plot)
    dev.off()
    
    hyper_NIH[[hyper_idx]] = list(GOM_res, GOM_inter, GOM_melt)
    hyper_idx = hyper_idx + 1 
    
  }

  table(BDR_metadata_AD$cross_clusters)
  table(NIH_AD$CVcross_clusters)
  table(UKBBN_AD$CVcross_clusters)
  
  sum(table(BDR_metadata_AD$cross_clusters))
  sum(table(NIH_AD$CVcross_clusters))
  sum(table(UKBBN_AD$CVcross_clusters))
  
  ##### Within cohort check ##### 
  
  #The goal of this entire step is to find out how many cpgs are shared between outcomes in a single cohort
  # and how many are shared across cohorts within a single outcome
  # We want to maximise the latter, and minimize the first (Outcome specifity)
  
  res_method = list()
  # Within a cohort, for each model, we compute the significant CPGs for each outcome
  # We then compare those 1-to-1 to verify the overlap
  for(method in unique(df_res_lm_BDR$model)){ 
    BDR_EWAS_spec = df_res_lm_BDR[df_res_lm_BDR$model == method,]
    UKBBN_EWAS_spec = df_res_lm_UKBBN[df_res_lm_UKBBN$model == method,]
    PITT_EWAS_spec = df_res_lm_NIH[df_res_lm_NIH$model == method,]

    listx = list()
    listy = list()
    
    for(nclustx in unique(BDR_EWAS_spec$cluster)){
      for(nclusty in unique(BDR_EWAS_spec$cluster)){
        #BDR
        tmp_EWAS_specx = BDR_EWAS_spec[BDR_EWAS_spec$cluster == nclustx,]
        tmp_EWAS_specy = BDR_EWAS_spec[BDR_EWAS_spec$cluster == nclusty,]
        max_value = nrow(tmp_EWAS_specx)
        # We select only the cpgs having the same direction to be included in the overlap
        tmp_EWAS_specx$eff_sizeok = FALSE
        tmp_EWAS_specx$effy = tmp_EWAS_specy$effect
        tmp_EWAS_specx[tmp_EWAS_specx$effect < 0 & tmp_EWAS_specx$effy < 0,]$eff_sizeok = TRUE
        tmp_EWAS_specx[tmp_EWAS_specx$effect > 0 & tmp_EWAS_specx$effy > 0,]$eff_sizeok = TRUE
        
        tmp_EWAS_specy$eff_sizeok = FALSE
        tmp_EWAS_specy$effx = tmp_EWAS_specx$effect
        tmp_EWAS_specy[tmp_EWAS_specy$effect < 0 & tmp_EWAS_specy$effx < 0,]$eff_sizeok = TRUE
        tmp_EWAS_specy[tmp_EWAS_specy$effect > 0 & tmp_EWAS_specy$effx > 0,]$eff_sizeok = TRUE      
        
        tmp_EWAS_specx = tmp_EWAS_specx[tmp_EWAS_specx$eff_sizeok == TRUE,]
        tmp_EWAS_specy = tmp_EWAS_specy[tmp_EWAS_specy$eff_sizeok == TRUE,]
        
        #Significance filter
        tmp_EWAS_specx = tmp_EWAS_specx[tmp_EWAS_specx$pvalue < 0.05 & tmp_EWAS_specx$cohend >= 0.2,]
        tmp_EWAS_specy = tmp_EWAS_specy[tmp_EWAS_specy$pvalue < 0.05 & tmp_EWAS_specy$cohend >= 0.2,]
        
        listx[[paste0("BDR;",nclustx)]] = tmp_EWAS_specx$cpg
        listy[[paste0("BDR;",nclusty)]] = tmp_EWAS_specy$cpg
        
        #UKBBN
        tmp_EWAS_specx = UKBBN_EWAS_spec[UKBBN_EWAS_spec$cluster == nclustx,]
        tmp_EWAS_specy = UKBBN_EWAS_spec[UKBBN_EWAS_spec$cluster == nclusty,]
        max_value = nrow(tmp_EWAS_specx)
        # We select only the cpgs having the same direction to be included in the overlap
        tmp_EWAS_specx$eff_sizeok = FALSE
        tmp_EWAS_specx$effy = tmp_EWAS_specy$effect
        tmp_EWAS_specx[tmp_EWAS_specx$effect < 0 & tmp_EWAS_specx$effy < 0,]$eff_sizeok = TRUE
        tmp_EWAS_specx[tmp_EWAS_specx$effect > 0 & tmp_EWAS_specx$effy > 0,]$eff_sizeok = TRUE
        
        tmp_EWAS_specy$eff_sizeok = FALSE
        tmp_EWAS_specy$effx = tmp_EWAS_specx$effect
        tmp_EWAS_specy[tmp_EWAS_specy$effect < 0 & tmp_EWAS_specy$effx < 0,]$eff_sizeok = TRUE
        tmp_EWAS_specy[tmp_EWAS_specy$effect > 0 & tmp_EWAS_specy$effx > 0,]$eff_sizeok = TRUE      
        
        tmp_EWAS_specx = tmp_EWAS_specx[tmp_EWAS_specx$eff_sizeok == TRUE,]
        tmp_EWAS_specy = tmp_EWAS_specy[tmp_EWAS_specy$eff_sizeok == TRUE,]
        #Significance filter
        tmp_EWAS_specx = tmp_EWAS_specx[tmp_EWAS_specx$pvalue < 0.05 & tmp_EWAS_specx$cohend >= 0.2,]
        tmp_EWAS_specy = tmp_EWAS_specy[tmp_EWAS_specy$pvalue < 0.05 & tmp_EWAS_specy$cohend >= 0.2,]
        
        listx[[paste0("UKBBN;",nclustx)]] = tmp_EWAS_specx$cpg
        listy[[paste0("UKBBN;",nclusty)]] = tmp_EWAS_specy$cpg
        
        #PITT
        tmp_EWAS_specx = PITT_EWAS_spec[PITT_EWAS_spec$cluster == nclustx,]
        tmp_EWAS_specy = PITT_EWAS_spec[PITT_EWAS_spec$cluster == nclusty,]
        max_value = nrow(tmp_EWAS_specx)
        # We select only the cpgs having the same direction to be included in the overlap
        tmp_EWAS_specx$eff_sizeok = FALSE
        tmp_EWAS_specx$effy = tmp_EWAS_specy$effect
        tmp_EWAS_specx[tmp_EWAS_specx$effect < 0 & tmp_EWAS_specx$effy < 0,]$eff_sizeok = TRUE
        tmp_EWAS_specx[tmp_EWAS_specx$effect > 0 & tmp_EWAS_specx$effy > 0,]$eff_sizeok = TRUE
        
        tmp_EWAS_specy$eff_sizeok = FALSE
        tmp_EWAS_specy$effx = tmp_EWAS_specx$effect
        tmp_EWAS_specy[tmp_EWAS_specy$effect < 0 & tmp_EWAS_specy$effx < 0,]$eff_sizeok = TRUE
        tmp_EWAS_specy[tmp_EWAS_specy$effect > 0 & tmp_EWAS_specy$effx > 0,]$eff_sizeok = TRUE      
        
        tmp_EWAS_specx = tmp_EWAS_specx[tmp_EWAS_specx$eff_sizeok == TRUE,]
        tmp_EWAS_specy = tmp_EWAS_specy[tmp_EWAS_specy$eff_sizeok == TRUE,]
        #Significance filter
        tmp_EWAS_specx = tmp_EWAS_specx[tmp_EWAS_specx$pvalue < 0.05 & tmp_EWAS_specx$cohend >= 0.2,]
        tmp_EWAS_specy = tmp_EWAS_specy[tmp_EWAS_specy$pvalue < 0.05 & tmp_EWAS_specy$cohend >= 0.2,]
    
        listx[[paste0("PITT;",nclustx)]] = tmp_EWAS_specx$cpg
        listy[[paste0("PITT;",nclusty)]] = tmp_EWAS_specy$cpg
      }
    }
    
    GOM_res = newGOM(listx, listy,
                     max_value)
    GOM_matrix = getMatrix(GOM_res, name="pval")
    GOM_melt = melt(GOM_matrix)
    
    GOM_inter = getNestedList(GOM_res, name="intersection")
    #Within a cohort, we look into the overlap of each significant EWAS results with another
    intracohort = list("BDR" = list("s1;s1" = GOM_inter$`BDR;1`$`BDR;1`,
                                    "s2;s2" = GOM_inter$`BDR;2`$`BDR;2`,
                                    "s1;s2" = GOM_inter$`BDR;1`$`BDR;2`,
                                    "s1;s3" = GOM_inter$`BDR;1`$`BDR;3`,
                                    "s2;s3" = GOM_inter$`BDR;2`$`BDR;3`,
                                    "s3;s3" = GOM_inter$`BDR;3`$`BDR;3`),
                       "UKBBN" = list("s1;s1" = GOM_inter$`UKBBN;1`$`UKBBN;1`,
                                    "s2;s2" = GOM_inter$`UKBBN;2`$`UKBBN;2`,
                                    "s1;s2" = GOM_inter$`UKBBN;1`$`UKBBN;2`,
                                    "s1;s3" = GOM_inter$`UKBBN;1`$`UKBBN;3`,
                                    "s2;s3" = GOM_inter$`UKBBN;2`$`UKBBN;3`,
                                    "s3;s3" = GOM_inter$`UKBBN;3`$`UKBBN;3`),
                       "PITT" = list("s1;s1" = GOM_inter$`PITT;1`$`PITT;1`,
                                    "s2;s2" = GOM_inter$`PITT;2`$`PITT;2`,
                                    "s1;s2" = GOM_inter$`PITT;1`$`PITT;2`,
                                    "s1;s3" = GOM_inter$`PITT;1`$`PITT;3`,
                                    "s2;s3" = GOM_inter$`PITT;2`$`PITT;3`,
                                    "s3;s3" = GOM_inter$`PITT;3`$`PITT;3`))
    
    #########
    
    # Within a cohort, for each model, we compute the significant CPGs for each cohort
    # We then compare the results for each outcome to verify the overlap across cohorts
    
    #Here we select only the significant cpgs we detected earlier within each cohort to work with
    BDR_s1 = BDR_EWAS_spec[BDR_EWAS_spec$cpg %in% GOM_inter$`BDR;1`$`BDR;1` & BDR_EWAS_spec$cluster == "1",]
    BDR_s2 = BDR_EWAS_spec[BDR_EWAS_spec$cpg %in% GOM_inter$`BDR;2`$`BDR;2` & BDR_EWAS_spec$cluster == "2",]
    BDR_s3 = BDR_EWAS_spec[BDR_EWAS_spec$cpg %in% GOM_inter$`BDR;3`$`BDR;3` & BDR_EWAS_spec$cluster == "3",]
    
    UKBBN_s1 = UKBBN_EWAS_spec[UKBBN_EWAS_spec$cpg %in% GOM_inter$`UKBBN;1`$`UKBBN;1` & UKBBN_EWAS_spec$cluster == "1",]
    UKBBN_s2 = UKBBN_EWAS_spec[UKBBN_EWAS_spec$cpg %in% GOM_inter$`UKBBN;2`$`UKBBN;2` & UKBBN_EWAS_spec$cluster == "2",]
    UKBBN_s3 = UKBBN_EWAS_spec[UKBBN_EWAS_spec$cpg %in% GOM_inter$`UKBBN;3`$`UKBBN;3` & UKBBN_EWAS_spec$cluster == "3",]
    
    PITT_s1 = PITT_EWAS_spec[PITT_EWAS_spec$cpg %in% GOM_inter$`PITT;1`$`PITT;1` & PITT_EWAS_spec$cluster == "1",]
    PITT_s2 = PITT_EWAS_spec[PITT_EWAS_spec$cpg %in% GOM_inter$`PITT;2`$`PITT;2` & PITT_EWAS_spec$cluster == "2",]
    PITT_s3 = PITT_EWAS_spec[PITT_EWAS_spec$cpg %in% GOM_inter$`PITT;3`$`PITT;3` & PITT_EWAS_spec$cluster == "3",]
    
    list_s1 = list(BDR_s1, UKBBN_s1, PITT_s1)
    list_s2 = list(BDR_s2, UKBBN_s2, PITT_s2)
    list_s3 = list(BDR_s3, UKBBN_s3, PITT_s3)
    
    list_s1_allinter = list()
    list_s2_allinter = list()
    list_s3_allinter = list()
    for(cohortx in c(1,2,3)){
      for(cohorty in c(1,2,3)){
        #We select two cohorts to work with
        if(cohortx == 1){
          ncohortx = "BDR"
        } else if (cohortx == 2){
          ncohortx = "UKBBN"
        } else if (cohortx == 3){
          ncohortx = "PITT"
        }
        
        if(cohorty == 1){
          ncohorty = "BDR"
        } else if (cohorty == 2){
          ncohorty = "UKBBN"
        } else if (cohorty == 3){
          ncohorty = "PITT"
        }
        #Outcome subtype 1
        tmp_cohort1 = list_s1[[cohortx]]
        tmp_cohort2 = list_s1[[cohorty]]
        tmp_cohort1 = tmp_cohort1[tmp_cohort1$cpg %in% tmp_cohort2$cpg,]
        tmp_cohort1 = tmp_cohort1[order(tmp_cohort1$cpg),]
        tmp_cohort2 = tmp_cohort2[tmp_cohort2$cpg %in% tmp_cohort1$cpg,]
        tmp_cohort2 = tmp_cohort2[order(tmp_cohort2$cpg),]
        #Effect size direction filter, some groups don't have any cpg so we clear for that
        tmp_cohort1$effect_ok = FALSE
        tryCatch({
        tmp_cohort1[tmp_cohort1$effect > 0 & tmp_cohort2$effect > 0,]$effect_ok = TRUE
        }, error=function(e){print("no")})
        tryCatch({
        tmp_cohort1[tmp_cohort1$effect < 0 & tmp_cohort2$effect < 0,]$effect_ok = TRUE
        }, error=function(e){print("no")})
        tmp_cohort1= tmp_cohort1[tmp_cohort1$effect_ok == TRUE,]  
        
        list_s1_allinter[[paste0(ncohortx, ";", ncohorty)]] = tmp_cohort1$cpg
        #Outcome subtype 2
        tmp_cohort1 = list_s2[[cohortx]]
        tmp_cohort2 = list_s2[[cohorty]]
        tmp_cohort1 = tmp_cohort1[tmp_cohort1$cpg %in% tmp_cohort2$cpg,]
        tmp_cohort1 = tmp_cohort1[order(tmp_cohort1$cpg),]
        tmp_cohort2 = tmp_cohort2[tmp_cohort2$cpg %in% tmp_cohort1$cpg,]
        tmp_cohort2 = tmp_cohort2[order(tmp_cohort2$cpg),]
        #Effect size direction filter, some groups don't have any cpg so we clear for that
        tmp_cohort1$effect_ok = FALSE
        tryCatch({
          tmp_cohort1[tmp_cohort1$effect > 0 & tmp_cohort2$effect > 0,]$effect_ok = TRUE
        }, error=function(e){print("no")})
        tryCatch({
        tmp_cohort1[tmp_cohort1$effect < 0 & tmp_cohort2$effect < 0,]$effect_ok = TRUE
        }, error=function(e){print("no")})
        tmp_cohort1= tmp_cohort1[tmp_cohort1$effect_ok == TRUE,]  
        
        list_s2_allinter[[paste0(ncohortx, ";", ncohorty)]] = tmp_cohort1$cpg
        
        #Outcome subtype 3
        tmp_cohort1 = list_s3[[cohortx]]
        tmp_cohort2 = list_s3[[cohorty]]
        tmp_cohort1 = tmp_cohort1[tmp_cohort1$cpg %in% tmp_cohort2$cpg,]
        tmp_cohort1 = tmp_cohort1[order(tmp_cohort1$cpg),]
        tmp_cohort2 = tmp_cohort2[tmp_cohort2$cpg %in% tmp_cohort1$cpg,]
        tmp_cohort2 = tmp_cohort2[order(tmp_cohort2$cpg),]
        #Effect size direction filter, some groups don't have any cpg so we clear for that
        tmp_cohort1$effect_ok = FALSE
        tryCatch({
          tmp_cohort1[tmp_cohort1$effect > 0 & tmp_cohort2$effect > 0,]$effect_ok = TRUE
        }, error=function(e){print("no")})
        tryCatch({
          tmp_cohort1[tmp_cohort1$effect < 0 & tmp_cohort2$effect < 0,]$effect_ok = TRUE
        }, error=function(e){print("no")})
        tmp_cohort1= tmp_cohort1[tmp_cohort1$effect_ok == TRUE,]  
        
        list_s3_allinter[[paste0(ncohortx, ";", ncohorty)]] = tmp_cohort1$cpg
        
      }
    }
    #Outcome subtype 4, if crossmethod
    if(length(unique(BDR_EWAS_spec$cluster)) == 4){
      BDR_s4 = BDR_EWAS_spec[BDR_EWAS_spec$cpg %in% GOM_inter$`BDR;4`$`BDR;4` & BDR_EWAS_spec$cluster == "4",]
      UKBBN_s4 = UKBBN_EWAS_spec[UKBBN_EWAS_spec$cpg %in% GOM_inter$`UKBBN;4`$`UKBBN;4` & UKBBN_EWAS_spec$cluster == "4",]  
      PITT_s4 = PITT_EWAS_spec[PITT_EWAS_spec$cpg %in% GOM_inter$`PITT;4`$`PITT;4` & PITT_EWAS_spec$cluster == "4",]
      list_s4 = list(BDR_s4, UKBBN_s4, PITT_s4)
      
      list_s4_allinter = list()
      for(cohortx in c(1,2,3)){
        for(cohorty in c(1,2,3)){
          tmp_cohort1 = list_s4[[cohortx]]
          tmp_cohort2 = list_s4[[cohorty]]
          tmp_cohort1 = tmp_cohort1[tmp_cohort1$cpg %in% tmp_cohort2$cpg,]
          tmp_cohort1 = tmp_cohort1[order(tmp_cohort1$cpg),]
          tmp_cohort2 = tmp_cohort2[tmp_cohort2$cpg %in% tmp_cohort1$cpg,]
          tmp_cohort2 = tmp_cohort2[order(tmp_cohort2$cpg),]
          
          tmp_cohort1$effect_ok = FALSE
          tryCatch({
            tmp_cohort1[tmp_cohort1$effect > 0 & tmp_cohort2$effect > 0,]$effect_ok = TRUE
          }, error=function(e){print("no")})
          tryCatch({
            tmp_cohort1[tmp_cohort1$effect < 0 & tmp_cohort2$effect < 0,]$effect_ok = TRUE
          }, error=function(e){print("no")})
          tmp_cohort1= tmp_cohort1[tmp_cohort1$effect_ok == TRUE,]  
          
          if(cohortx == 1){
            ncohortx = "BDR"
          } else if (cohortx == 2){
            ncohortx = "UKBBN"
          } else if (cohortx == 3){
            ncohortx = "PITT"
          }
          
          if(cohorty == 1){
            ncohorty = "BDR"
          } else if (cohorty == 2){
            ncohorty = "UKBBN"
          } else if (cohorty == 3){
            ncohorty = "PITT"
          }
          
          list_s4_allinter[[paste0(ncohortx, ";",ncohorty)]] = tmp_cohort1$cpg
          
        }
      }
    }
    
    
    if(length(unique(BDR_EWAS_spec$cluster)) == 4){
      res_method[[method]] = list("Intracohort" = intracohort, 
                                  "IntercohortS1" = list_s1_allinter,
                                  "IntercohortS2" = list_s2_allinter,
                                  "IntercohortS3" = list_s3_allinter,
                                  "IntercohortS4" = list_s4_allinter)
    } else {
      res_method[[method]] = list("Intracohort" = intracohort, 
                                  "IntercohortS1" = list_s1_allinter,
                                  "IntercohortS2" = list_s2_allinter,
                                  "IntercohortS3" = list_s3_allinter)
    }
  }
  
  
  
  #########
  #Compiling the results of the previous steps in a table
  res_intracohort = data.frame(matrix(nrow = 0, ncol = 11))
  for(method in names(res_method)){
    for(cohort in names(res_method$PAM$Intracohort)){
      
      s1s1 = length(res_method[[method]]$Intracohort[[cohort]]$`s1;s1`)
      s2s2 = length(res_method[[method]]$Intracohort[[cohort]]$`s2;s2`)
      s3s3 = length(res_method[[method]]$Intracohort[[cohort]]$`s3;s3`)
      s1s2 = length(res_method[[method]]$Intracohort[[cohort]]$`s1;s2`)
      s1s3 = length(res_method[[method]]$Intracohort[[cohort]]$`s1;s3`)
      s2s3 = length(res_method[[method]]$Intracohort[[cohort]]$`s2;s3`)
      
      # Within a cohort, we compute the ratio of shared cpgs between two clusters
      # divided by the largest list of the pair
      if(s1s1 > s2s2){
        ratios1s2 = s1s2/ s2s2
      } else if(s1s1 < s2s2){
        ratios1s2 = s1s2 / s1s1
      }
      
      if(s1s1 > s3s3){
        ratios1s3 = s1s3 / s3s3 
      } else if(s1s1 < s3s3){
        ratios1s3 = s1s3 / s1s1
      }
      
      if(s2s2 > s3s3){
        ratios2s3 = s2s3 / s3s3  
      } else if(s2s2 < s3s3){
        ratios2s3 = s2s3 / s2s2
      }
      newline = c(method,cohort,s1s1,s2s2,s3s3,s1s2,ratios1s2,s1s3,ratios1s3,s2s3,ratios2s3)
      res_intracohort = rbind(res_intracohort, newline)
    }
  }
    
  colnames(res_intracohort) = c("method","cohort","s1;s1","s2;s2","s3;s3","s1;s2","ratio;s1;s2","s1;s3","ratio;s1;s3","s2;s3","ratio;s2;s3")
  
  #######
  
  #We do the same in between cohorts 
  res_intercohort = data.frame(matrix(nrow = 0, ncol = 11))
  for(method in names(res_method)){
      
      BDRs1 = length(res_method[[method]]$IntercohortS1$`BDR;BDR`)
      UKBBNs1 = length(res_method[[method]]$IntercohortS1$`UKBBN;UKBBN`)
      PITTs1 = length(res_method[[method]]$IntercohortS1$`PITT;PITT`)
      BU1 = length(res_method[[method]]$IntercohortS1$`BDR;UKBBN`)
      BP1 = length(res_method[[method]]$IntercohortS1$`BDR;PITT`)
      UP1 = length(res_method[[method]]$IntercohortS1$`UKBBN;PITT`)
      
      BDRs2 = length(res_method[[method]]$IntercohortS2$`BDR;BDR`)
      UKBBNs2 = length(res_method[[method]]$IntercohortS2$`UKBBN;UKBBN`)
      PITTs2 = length(res_method[[method]]$IntercohortS2$`PITT;PITT`)
      BU2 = length(res_method[[method]]$IntercohortS2$`BDR;UKBBN`)
      BP2 = length(res_method[[method]]$IntercohortS2$`BDR;PITT`)
      UP2 = length(res_method[[method]]$IntercohortS2$`UKBBN;PITT`)
      
      BDRs3 = length(res_method[[method]]$IntercohortS3$`BDR;BDR`)
      UKBBNs3 = length(res_method[[method]]$IntercohortS3$`UKBBN;UKBBN`)
      PITTs3 = length(res_method[[method]]$IntercohortS3$`PITT;PITT`)
      BU3 = length(res_method[[method]]$IntercohortS3$`BDR;UKBBN`)
      BP3 = length(res_method[[method]]$IntercohortS3$`BDR;PITT`)
      UP3 = length(res_method[[method]]$IntercohortS3$`UKBBN;PITT`)
      
      # How much is shared vs how much is non-shared
      if(BDRs1 > UKBBNs1){
        ratioBU1 = BU1 / UKBBNs1
      } else if(BDRs1 < UKBBNs1){
        ratioBU1 = BU1 / BDRs1
      }
      
      if(BDRs1 > PITTs1){
        ratioBP1 = BP1 / PITTs1 
      } else if(BDRs1 < PITTs1){
        ratioBP1 = BP1 / BDRs1
      }
      
      if(UKBBNs1 > PITTs1){
        ratioUP1 = UP1 / PITTs1  
      } else if(UKBBNs1 < PITTs1){
        ratioUP1 = UP1 / UKBBNs1
      }
      
      if(BDRs2 > UKBBNs2){
        ratioBU2 = BU2 / UKBBNs2
      } else if(BDRs2 < UKBBNs2){
        ratioBU2 = BU2 / BDRs2
      }
      
      if(BDRs2 > PITTs2){
        ratioBP2 = BP2 / PITTs2 
      } else if(BDRs2 < PITTs2){
        ratioBP2 = BP2 / BDRs2
      }
      
      if(UKBBNs2 > PITTs2){
        ratioUP2 = UP2 / PITTs2  
      } else if(UKBBNs2 < PITTs2){
        ratioUP2 = UP2 / UKBBNs2
      }
      
      if(BDRs3 > UKBBNs3){
        ratioBU3 = BU3 / UKBBNs3
      } else if(BDRs3 < UKBBNs3){
        ratioBU3 = BU3 / BDRs3
      }
      
      if(BDRs3 > PITTs3){
        ratioBP3 = BP3 / PITTs3 
      } else if(BDRs3 < PITTs3){
        ratioBP3 = BP3 / BDRs3
      }
      
      if(UKBBNs3 > PITTs3){
        ratioUP3 = UP3 / PITTs3  
      } else if(UKBBNs3 < PITTs3){
        ratioUP3 = UP3 / UKBBNs3
      }
      
      
      
      newline1 = c(method,1,BDRs1,UKBBNs1,PITTs1,BU1,ratioBU1,BP1,ratioBP1,UP1,ratioUP1)
      newline2 = c(method,2,BDRs2,UKBBNs2,PITTs2,BU2,ratioBU2,BP2,ratioBP2,UP2,ratioUP2)
      newline3 = c(method,3,BDRs3,UKBBNs3,PITTs3,BU3,ratioBU3,BP3,ratioBP3,UP3,ratioUP3)
      res_intercohort = rbind(res_intercohort, newline1)
      res_intercohort = rbind(res_intercohort, newline2)
      res_intercohort = rbind(res_intercohort, newline3)
  
  }
  colnames(res_intercohort) = c("method","subtype","BDR","UKBBN","PITT","BDR/UKBBN","ratio;BDR/UKBBN","BDR/PITT","ratio;BDR/PITT","UKBBN/PITT","ratio;UKBBN/PITT")
  
  
  ##### PCA #####
  
  #Merging betas sets of all cohorts
  CV_spec_BDR = mSet_betas_AD[rownames(mSet_betas) %in% Spec_final_features, !is.na(BDR_metadata_AD$Spec_clusters)]
  CV_spec_UKBBN = msetEPIC_betas[rownames(msetEPIC_betas) %in% Spec_final_features, colnames(msetEPIC_betas) %in% UKBBN_AD$Row.names]
  CV_spec_NIH = NIH_betas[rownames(NIH_betas) %in% Spec_final_features, colnames(NIH_betas) %in% NIH_AD$Basename]
  CV_spec_all = cbind(CV_spec_BDR, CV_spec_UKBBN, CV_spec_NIH)

  #Generating crosscohort PCA
  labels_pc = c(BDR_metadata_AD$Spec_clusters[!is.na(BDR_metadata_AD$Spec_clusters)], UKBBN_AD$CVSpec_clusters, NIH_AD$CVSpec_clusters)
  cohort_pc = c(rep("BDR", ncol(CV_spec_BDR)), rep("UKBBN", ncol(CV_spec_UKBBN)), rep("NIH", ncol(CV_spec_NIH)))
  col_pc = labels_pc
  col_pc[col_pc == 1] = "steelblue2"
  col_pc[col_pc == 2] = "darkorange1"
  col_pc[col_pc == 3] = "azure3"
  pc_all = prcomp(t(CV_spec_all))
  
  pdf(paste0("results/CV_spectrum_map.pdf"), width = 7.5, height = 7.5)
  fviz_pca_ind(pc_all, geom = "point", alpha.ind = 0) +
    geom_point(aes(shape = factor(cohort_pc), colour = as.character(labels_pc))) +
    scale_color_manual(values = c("1" = "steelblue2",
                          "2" = "darkorange1",
                          "3" = "azure3")) +
    scale_shape_manual(values = c(19, 2, 3)) +
    labs(shape = "Cohort", color = "Subtype") + ggtitle("Cohorts map of individuals")
  dev.off()
  
  ##### Median samples #####
  
  #Creating median samples for each outcome for all cohorts (CV crossmethod)
  cor_median = data.frame(matrix(NA, nrow = 0, ncol = nrow(data_matrix_AD)))
  allrownames = c()
  for(subtype in unique(UKBBN_AD$CVSpec_clusters)){
    subtype_samples = BDR_metadata_AD[BDR_metadata_AD$Spec_clusters == subtype,]$Basename
    BDR_CVCM_median = data_matrix_AD[,colnames(data_matrix_AD) %in% subtype_samples]
    BDR_CVCM_median = rowMedians(BDR_CVCM_median)
    
    subtype_samples = UKBBN_AD[UKBBN_AD$CVSpec_clusters == subtype,]$Basename
    UKBBN_CVCM_median = msetEPIC_betas[rownames(msetEPIC_betas) %in% cpgs_intersect,
                                       colnames(msetEPIC_betas) %in% subtype_samples]
    UKBBN_CVCM_median = rowMedians(UKBBN_CVCM_median)
    
    subtype_samples = as.character(NIH_AD[NIH_AD$CVSpec_clusters == subtype,]$Basename)
    NIH_CVCM_median = NIH_betas[rownames(NIH_betas) %in% cpgs_intersect,
                                colnames(NIH_betas) %in% subtype_samples]
    NIH_CVCM_median = rowMedians(NIH_CVCM_median)
    
    cor_median = rbind(cor_median, BDR_CVCM_median) 
    cor_median = rbind(cor_median, UKBBN_CVCM_median)
    cor_median = rbind(cor_median, NIH_CVCM_median)
    
    allrownames = c(allrownames, 
                    paste0("BDRS", subtype),
                    paste0("UKBBNS", subtype),
                    paste0("NIHS", subtype))
  }
  
  rownames(cor_median) = allrownames
  colnames(cor_median) = rownames(data_matrix_AD)
  
  pdf("./results/corr_median_sample.pdf", width = 7, height = 7)
  corrplot(cor(t(cor_median)))
  dev.off()
  
  ##### RNAseq analysis #####
  
  # load("./data/UKBBN_RNAseq_reads.Rdata")
  # UKBBN_counts = counts
  # UKBBN_pheno = pheno
  # load("./data/PITTS_RNAseq_counts.Rdata")
  # PITT_counts = counts
  # PITT_pheno = pheno
  #  
  load("./RNA_study.Rdata")
  
  # Preparing RNAseq data to match previous steps and labeling
  UKBBN_TINs = read.csv("./data/UKBBNtin_score.csv")
  UKBBN_TINs = UKBBN_TINs[UKBBN_TINs$ID %in% UKBBN_pheno$BBNId,]
  colnames(UKBBN_TINs)[1] = "BBNId"
  UKBBN_celltypes = read.csv("./data/UKBBNCSH_cibersort_results.csv")
  UKBBN_celltypes = UKBBN_celltypes[UKBBN_celltypes$X %in% UKBBN_pheno$BBNId,]
  colnames(UKBBN_celltypes)[1] = "BBNId"
  UKBBN_pheno = UKBBN_pheno[order(UKBBN_pheno$BBNId),]
  UKBBN_celltypes = UKBBN_celltypes[order(UKBBN_celltypes$BBNId),]
  UKBBN_TINs = UKBBN_TINs[order(UKBBN_TINs$BBNId),]
  
  
  UKBBN_pheno = cbind(UKBBN_pheno, UKBBN_celltypes[,-1], UKBBN_TINs[,-1])
  colnames(UKBBN_pheno)[69] = "TIN"
  
  ctrl_filt = UKBBN_Ctrl[!UKBBN_Ctrl$Row.names %in% UKBBN_Ctrl_filt$Row.names,]$Row.names
  UKBBN_pheno = UKBBN_pheno[!UKBBN_pheno$X %in% ctrl_filt,]
  
  colnames(PITT_pheno)[1] = "Individual_ID"
  
  PITT_pheno = merge(PITT_pheno, NIH_pheno, by = "Individual_ID")
  
  UKBBN_pheno$subtype= UKBBN_pheno$CVSpec_clusters
  
  PITT_pheno$CVSpec_clusters = NA
  
  PITT_pheno$CVSpec_clusters = NIH_AD$CVSpec_clusters[match(PITT_pheno$Individual_ID, NIH_AD$Individual_ID)]
  PITT_pheno[is.na(PITT_pheno$CVSpec_clusters),]$CVSpec_clusters = "Control"
  table(PITT_pheno$CVSpec_clusters)
  table(PITT_pheno$subtype)
  PITT_pheno$subtype= PITT_pheno$CVSpec_clusters
  PITT_pheno$Sex.y = factor(PITT_pheno$Sex.y, levels = c("M","F"))
  
  #Running subtype X vs controls
  UKBBN_edgeR_s1 = run_edgeR(UKBBN_counts, UKBBN_pheno, subtype = "1", cohort = "UKBBN")
  UKBBN_edgeR_s1$subtype = "1"
  UKBBN_edgeR_s1$cohort = "UKBBN"
  UKBBN_edgeR_s1$ENSG = rownames(UKBBN_edgeR_s1)
  UKBBN_edgeR_s3 = run_edgeR(UKBBN_counts, UKBBN_pheno, subtype = "3", cohort = "UKBBN")
  UKBBN_edgeR_s3$subtype = "3"
  UKBBN_edgeR_s3$cohort = "UKBBN"
  UKBBN_edgeR_s3$ENSG = rownames(UKBBN_edgeR_s3)
  PITT_edgeR_s1 = run_edgeR(PITT_counts, PITT_pheno, subtype = "1", cohort = "PITT")
  PITT_edgeR_s1$subtype = "1"
  PITT_edgeR_s1$cohort = "PITT"
  PITT_edgeR_s1$ENSG = rownames(PITT_edgeR_s1)
  PITT_edgeR_s3 = run_edgeR(PITT_counts, PITT_pheno, subtype = "3", cohort = "PITT")
  PITT_edgeR_s3$subtype = "3"
  PITT_edgeR_s3$cohort = "PITT"
  PITT_edgeR_s3$ENSG = rownames(PITT_edgeR_s3)
  ##### RNAseq PITT vs UKBBN #####
  
  UKBBN_edgeR = rbind(UKBBN_edgeR_s1, UKBBN_edgeR_s3)
  PITT_edgeR = rbind(PITT_edgeR_s1, PITT_edgeR_s3)
  
  UKBBN_edgeR_s1 = UKBBN_edgeR[UKBBN_edgeR$subtype == "1",]
  UKBBN_edgeR_s3 = UKBBN_edgeR[UKBBN_edgeR$subtype == "3",]
  PITT_edgeR_s3 = PITT_edgeR[PITT_edgeR$subtype == "3",]
  PITT_edgeR_s1 = PITT_edgeR[PITT_edgeR$subtype == "1",]
  
  # Take shared ENSGs between cohorts and subtypes
  UKBBN_edgeR_red = UKBBN_edgeR[UKBBN_edgeR$ENSG %in% Reduce(intersect, list(PITT_edgeR_s1$ENSG,PITT_edgeR_s3$ENSG,UKBBN_edgeR_s1$ENSG,UKBBN_edgeR_s3$ENSG)),]
  PITT_edgeR_red = PITT_edgeR[PITT_edgeR$ENSG %in% Reduce(intersect, list(PITT_edgeR_s1$ENSG,PITT_edgeR_s3$ENSG,UKBBN_edgeR_s1$ENSG,UKBBN_edgeR_s3$ENSG)),]
  
  all_RNA_merge = merge(UKBBN_edgeR_red, PITT_edgeR_red, by = c("subtype","ENSG"))
  n_RNA = nrow(all_RNA_merge)
  all_RNA_merge$fc_ok = FALSE

  # FC direction check
  all_RNA_merge[all_RNA_merge$logFC.x > 0 & all_RNA_merge$logFC.y > 0,]$fc_ok = TRUE
  all_RNA_merge[all_RNA_merge$logFC.x < 0 & all_RNA_merge$logFC.y < 0,]$fc_ok = TRUE

  all_RNA_merge = all_RNA_merge[all_RNA_merge$fc_ok == TRUE,]

  all_RNA_merge_s1 = all_RNA_merge[all_RNA_merge$subtype == "1",]
  all_RNA_merge_s3 = all_RNA_merge[all_RNA_merge$subtype == "3",]
  
  PITT_edgeR_s1 = PITT_edgeR_s1[PITT_edgeR_s1$ENSG %in% all_RNA_merge_s1$ENSG,]
  UKBBN_edgeR_s1 = UKBBN_edgeR_s1[UKBBN_edgeR_s1$ENSG %in% all_RNA_merge_s1$ENSG,]
  PITT_edgeR_s3 = PITT_edgeR_s3[PITT_edgeR_s3$ENSG %in% all_RNA_merge_s3$ENSG,]
  UKBBN_edgeR_s3 = UKBBN_edgeR_s3[UKBBN_edgeR_s3$ENSG %in% all_RNA_merge_s3$ENSG,]
  
  
  PITT_edgeR_s1$BHp_value = p.adjust(PITT_edgeR_s1$PValue, method = "BH")
  PITT_edgeR_s3$BHp_value = p.adjust(PITT_edgeR_s3$PValue, method = "BH")
  UKBBN_edgeR_s1$BHp_value = p.adjust(UKBBN_edgeR_s1$PValue, method = "BH")
  UKBBN_edgeR_s3$BHp_value = p.adjust(UKBBN_edgeR_s3$PValue, method = "BH")
  
  # Significance of each subtype/cohort
  RNA_subtype1_PITTS = PITT_edgeR_s1[PITT_edgeR_s1$PValue < 0.05,]
  # RNA_subtype1_PITTS = RNA_subtype1_PITTS[abs(RNA_subtype1_PITTS$logFC) >= 0.2,]
  RNA_subtype1_UKBBN = UKBBN_edgeR_s1[UKBBN_edgeR_s1$PValue < 0.05,]
  # RNA_subtype1_UKBBN = RNA_subtype1_UKBBN[abs(RNA_subtype1_UKBBN$logFC) >= 0.2,]
  RNA_subtype3_PITTS = PITT_edgeR_s3[PITT_edgeR_s3$PValue < 0.05,]
  # RNA_subtype3_PITTS = RNA_subtype3_PITTS[abs(RNA_subtype3_PITTS$logFC) >= 0.2,]
  RNA_subtype3_UKBBN = UKBBN_edgeR_s3[UKBBN_edgeR_s3$PValue < 0.05,]
  # RNA_subtype3_UKBBN = RNA_subtype3_UKBBN[abs(RNA_subtype3_UKBBN$logFC) >= 0.2,]
  
  
  list_RNA_PITT = list(RNA_subtype1_PITTS$ENSG, RNA_subtype3_PITTS$ENSG, RNA_subtype1_UKBBN$ENSG, RNA_subtype3_UKBBN$ENSG)
  list_RNA_UKBBN = list(RNA_subtype1_UKBBN$ENSG, RNA_subtype3_UKBBN$ENSG, RNA_subtype1_PITTS$ENSG, RNA_subtype3_PITTS$ENSG)
  
  GOM_res = newGOM(list_RNA_PITT, list_RNA_UKBBN, genome.size = n_RNA)
  GOM_matrix = getMatrix(GOM_res, name="pval")
  GOM_melt = melt(GOM_matrix)
  
  GOM_melt[GOM_melt$Var1 == "1",]$Var1 = "PITT1"
  GOM_melt[GOM_melt$Var1 == "2",]$Var1 = "PITT3"
  GOM_melt[GOM_melt$Var1 == "3",]$Var1 = "UK1"
  GOM_melt[GOM_melt$Var1 == "4",]$Var1 = "UK3"
  GOM_melt[GOM_melt$Var2 == "1",]$Var2 = "UK1"
  GOM_melt[GOM_melt$Var2 == "2",]$Var2 = "UK3"
  GOM_melt[GOM_melt$Var2 == "3",]$Var2 = "PITT1"
  GOM_melt[GOM_melt$Var2 == "4",]$Var2 = "PITT3"
  
  GOM_inter = getNestedList(GOM_res, name="intersection")
  
  
  
  GOM_melt$jacc = NA
  GOM_index = 1
  for(UKBBNi in seq(1:length(list_RNA_UKBBN))){
    for(BDRi in seq(1:length(list_RNA_PITT))){
      
      BDRcpg = list_RNA_PITT[[BDRi]]
      UKBBNcpg = list_RNA_UKBBN[[UKBBNi]]
      
      jaccI = length(intersect(BDRcpg, UKBBNcpg)) / length(union(BDRcpg, UKBBNcpg))
      
      GOM_melt[GOM_index,]$jacc = jaccI
      GOM_index = GOM_index + 1
    }
  }
  
  GOM_melt$value =  as.numeric(format(GOM_melt$value, scientific = TRUE, digits = 2))
  GOM_melt[GOM_melt$value > 0.05,]$value = "N.S."
  colnames(GOM_melt) = c("PITTADR", "UKBBN", "value", "Jaccard_index")
  
  GOM_plot = ggplot(GOM_melt, aes(x = PITTADR, y = UKBBN, fill = Jaccard_index))+
    ggtitle("Subtype specifity - PITT/UKBBN cohort") +
    geom_tile() +
    geom_text(aes(label = value)) +
    scale_fill_gradient2(low = "white", high = "chocolate2")+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  GOM_melt
  # pdf(paste0("results/Jaccard_RNAseq.pdf"), width = 7, height = 7)
  plot(GOM_plot)
  # dev.off()
  
  sub13 = all_RNA_merge[all_RNA_merge$ENSG %in% GOM_res@go.nested.list[[1]][[2]]@intersection,]
  sub11 = all_RNA_merge[all_RNA_merge$ENSG %in% GOM_res@go.nested.list[[1]][[1]]@intersection,]
  sub11 = sub11[sub11$subtype == "1",]
  sub11 = sub11[sub11$ENSG %in% sub13$ENSG,]
  sub33 = all_RNA_merge[all_RNA_merge$ENSG %in% GOM_res@go.nested.list[[2]][[2]]@intersection,]
  sub33 = sub33[sub33$subtype == "3",]
  sub33 = sub33[sub33$ENSG %in% sub13$ENSG,]
  
  
  plot(sub11$logFC.x, sub11$logFC.y)
  plot(-log10(sub11$PValue.x), -log10(sub11$PValue.y))
  plot(sub33$logFC.x, sub33$logFC.y)
  plot(-log10(sub33$PValue.x), -log10(sub33$PValue.y))
  plot(sub13$logFC.x, sub13$logFC.y)
  plot(-log10(sub13$PValue.x), -log10(sub13$PValue.y))
  
  
  biomart_names = read.table("./data/togenenames.txt", sep = "\t")
  colnames(biomart_names) = biomart_names[1,]
  biomart_names = biomart_names[-1,]
  
  ##### RNAseq PITT / UKBBN #####
  
  #Volcano plots of single subtype vs control
  UKBBN_vplot_s1 = UKBBN_edgeR[UKBBN_edgeR$subtype == "1",]
  UKBBN_vplot_s3 = UKBBN_edgeR[UKBBN_edgeR$subtype == "3",]
  PITT_vplot_s3 = PITT_edgeR[PITT_edgeR$subtype == "3",]
  PITT_vplot_s1 = PITT_edgeR[PITT_edgeR$subtype == "1",]
  
  rownames(UKBBN_vplot_s1) = UKBBN_vplot_s1$ENSG
  rownames(UKBBN_vplot_s3) = UKBBN_vplot_s3$ENSG
  rownames(PITT_vplot_s1) = PITT_vplot_s1$ENSG
  rownames(PITT_vplot_s3) = PITT_vplot_s3$ENSG
  
  vplot_edgeR(UKBBN_vplot_s1, stat_value = "PValue", stat_threshold = 2, FC_threshold = 0.5)
  vplot_edgeR(UKBBN_vplot_s3, stat_value = "PValue", stat_threshold = 2, FC_threshold = 0.5)
  vplot_edgeR(PITT_vplot_s1, stat_value = "PValue", stat_threshold = 2, FC_threshold = 0.5)
  vplot_edgeR(PITT_vplot_s3, stat_value = "PValue", stat_threshold = 2, FC_threshold = 0.5)
  
  #Associate ENSGs to gene names
  UKBBN_vplot_s1 = biomart_gene_assoc(UKBBN_vplot_s1, biomart_names, Pval_threshold = 2, FC_threshold = 0.5)
  UKBBN_vplot_s3 = biomart_gene_assoc(UKBBN_vplot_s3, biomart_names, Pval_threshold = 2, FC_threshold = 0.5)
  PITT_vplot_s3 = biomart_gene_assoc(PITT_vplot_s3, biomart_names, Pval_threshold = 2, FC_threshold = 0.5)
  PITT_vplot_s1 = biomart_gene_assoc(PITT_vplot_s1, biomart_names, Pval_threshold = 2, FC_threshold = 0.5)
  
  ##### Coloc #####
  
  ADGWAS_coloc = read.table("C:/Users/p70077107/Desktop/AD_kunkle_formatted.txt", sep = "\t")
  
  mQTL = read.table("C:/Users/p70077107/Desktop/ROSMAP.mQTL.2021.csv", sep = ",")
  
  MAFs = read.table("C:/Users/p70077107/Desktop/allchrMAF.txt")
  
  
  # Generating coloc mQTL list
  EWAS_to_mQTL = df_res_final[df_res_final$model == "CV_Crossmethod",]
  
  # Selecting significant CPGs across cohorts 
  EWAS_to_mQTL = EWAS_to_mQTL[EWAS_to_mQTL$BDR_pvalue < 0.05 &
                                EWAS_to_mQTL$UKBBN_pvalue < 0.05 &
                                EWAS_to_mQTL$NIH_pvalue < 0.05,]
  EWAS_to_mQTL = EWAS_to_mQTL[EWAS_to_mQTL$cluster == c("1", "3"),]
  EWAS_to_mQTL_s1 = EWAS_to_mQTL[EWAS_to_mQTL$cluster == "1",]  
  EWAS_to_mQTL_s3 = EWAS_to_mQTL[EWAS_to_mQTL$cluster == "3",]  
  EWAS_to_mQTL_s1$effok = FALSE
  # Effect size correction
  EWAS_to_mQTL_s1[EWAS_to_mQTL_s1$BDR_effectsize > 0 &
                    EWAS_to_mQTL_s1$UKBBN_effectsize > 0 &
                    EWAS_to_mQTL_s1$NIH_effectsize > 0,]$effok = TRUE
  EWAS_to_mQTL_s1[EWAS_to_mQTL_s1$BDR_effectsize < 0 &
                    EWAS_to_mQTL_s1$UKBBN_effectsize < 0 &
                    EWAS_to_mQTL_s1$NIH_effectsize < 0,]$effok = TRUE
  EWAS_to_mQTL_s1 = EWAS_to_mQTL_s1[EWAS_to_mQTL_s1$effok == TRUE,]
  
  
  EWAS_to_mQTL_s3$effok = FALSE
  EWAS_to_mQTL_s3[EWAS_to_mQTL_s3$BDR_effectsize > 0 &
                    EWAS_to_mQTL_s3$UKBBN_effectsize > 0 &
                    EWAS_to_mQTL_s3$NIH_effectsize > 0,]$effok = TRUE
  EWAS_to_mQTL_s3[EWAS_to_mQTL_s3$BDR_effectsize < 0 &
                    EWAS_to_mQTL_s3$UKBBN_effectsize < 0 &
                    EWAS_to_mQTL_s3$NIH_effectsize < 0,]$effok = TRUE
  EWAS_to_mQTL_s3 = EWAS_to_mQTL_s3[EWAS_to_mQTL_s3$effok == TRUE,]
  
  # Matching with mQTL database
  colnames(mQTL) = mQTL[1,]
  mQTL = mQTL[-1,]  
  mQTL_s1 = mQTL[mQTL$CpG %in% EWAS_to_mQTL_s1$cpg,]
  mQTL_s3 = mQTL[mQTL$CpG %in% EWAS_to_mQTL_s3$cpg,]
  
  # Adding MAFs to database
  colnames(MAFs) = c("rsID", "MAF")
  MAFs = MAFs[MAFs$rsID %in% c(mQTL_s1$rsID, mQTL_s3$rsID),]  
  
  # Coloc list with subtype 1
  mQTL_s1 = mQTL_s1[mQTL_s1$rsID %in% MAFs$rsID,]  
  mQTL_s1 = mQTL_s1[order(mQTL_s1$rsID),]
  MAFs = MAFs[order(MAFs$rsID),]
  
  mQTL_s1$MAF = MAFs$MAF
  
  mQTL_coloc_s1 = list(beta = as.numeric(mQTL_s1$beta), 
                       varbeta = as.numeric(mQTL_s1$se)^2,
                       type = "quant",
                       snp = mQTL_s1$rsID, 
                       MAF = as.numeric(mQTL_s1$MAF),
                       N = 543)
  
  #Filtering the ADGWAS and list to coloc
  colnames(ADGWAS_coloc) = ADGWAS_coloc[1,]
  ADGWAS_coloc = ADGWAS_coloc[-1,]
  ADGWAS_coloc = ADGWAS_coloc[!is.na(ADGWAS_coloc$RSid),]
  ADGWAS_coloc_data = list(beta = as.numeric(ADGWAS_coloc$beta),
                           varbeta = as.numeric(ADGWAS_coloc$se)^2,
                           type = "quant",
                           snp = ADGWAS_coloc$RSid,
                           MAF = as.numeric(ADGWAS_coloc$maf),
                           N = as.numeric(max(ADGWAS_coloc$n_sample)))
  check_dataset(mQTL_coloc_s1)
  check_dataset(ADGWAS_coloc_data)
  
  coloc_s1 = coloc.abf(mQTL_coloc_s1, ADGWAS_coloc_data)  
  
  ##### Metadata analysis #####
  
  UKBBN_newpheno = read.csv("./data/pheno_UKBBN.txt", sep = "\t")
  PITT_newpheno = read.csv("./data/pheno_nih.txt", sep = "\t")
  
  UKBBN_newpheno = UKBBN_newpheno[UKBBN_newpheno$BBNId %in% pheno$BBNId,]
  pheno = pheno[order(pheno$BBNId),]
  UKBBN_newpheno = UKBBN_newpheno[order(UKBBN_newpheno$BBNId),]
  UKBBN_newpheno$Basename = pheno$Row.names
  pheno = pheno[order(pheno$Row.names),]
  UKBBN_newpheno = UKBBN_newpheno[order(UKBBN_newpheno$Basename),]
  
  PITT_newpheno = PITT_newpheno[PITT_newpheno$SENTRIX_ID %in% NIH_pheno$Basename,]
  PITT_newpheno = PITT_newpheno[order(PITT_newpheno$SENTRIX_ID),]
  
  all_cohort_metadata = data.frame(matrix(NA, nrow = 843, ncol = 0))
  BDR_ok_ADCtrl = BDR_ok_ADCtrl[order(BDR_ok_ADCtrl$Basename),]
  all_cohort_metadata$Cohort = c(rep("BDR", nrow(BDR_ok_ADCtrl)),
                                 rep("UKBBN", nrow(pheno)),
                                 rep("PITT", nrow(NIH_pheno)))
  all_cohort_metadata$Basename = c(BDR_ok_ADCtrl$Basename,
                                   pheno$Row.names,
                                   as.character(NIH_pheno$Basename))
  all_cohort_metadata$Age = c(BDR_ok_ADCtrl$Age,
                              pheno$Age,
                              as.numeric(as.character(NIH_pheno$Age)))
  all_cohort_metadata$Gender = c(as.character(BDR_ok_ADCtrl$Gender),
                                 as.character(pheno$Gender),
                                 as.character(NIH_pheno$Sex))
  all_cohort_metadata[all_cohort_metadata$Gender == "M",]$Gender = "male"
  all_cohort_metadata[all_cohort_metadata$Gender == "F",]$Gender = "female"
  
  mSet_ecc = as.data.frame(mSet_ecc)
  
  all_cohort_metadata$NeuNpos = c(mSet_ecc$NeuN_pos,
                                  pheno$NeuN_pos,
                                  NIH_pheno$NeuN_pos)
  
  all_cohort_metadata$NeuNneg = c(mSet_ecc$NeuN_neg,
                                  pheno$NeuN_neg,
                                  NIH_pheno$NeuN_neg)
  
  all_cohort_metadata$Institute = c(BDR_ok_ADCtrl$Institute,
                                    as.character(pheno$Brain.Bank),
                                    as.character(NIH_pheno$Institute))
  
  all_cohort_metadata$plate = c(BDR_ok_ADCtrl$Plate,
                                pheno$plate,
                                as.character(NIH_pheno$Plate))
  
  all_cohort_metadata$Pathology = c(BDR_ok_ADCtrl$diag,
                                    pheno$AD,
                                    as.character(NIH_pheno$Phenotype))
  all_cohort_metadata[all_cohort_metadata$Pathology == "C",]$Pathology = "Control"
  all_cohort_metadata[all_cohort_metadata$Pathology == "0",]$Pathology = "Control"
  all_cohort_metadata[all_cohort_metadata$Pathology == "1",]$Pathology = "AD"
  
  
  BDR_ok_ADCtrl$subtype = NA
  BDR_ok_ADCtrl$subtype = BDR_metadata_AD$Spec_clusters[match(BDR_ok_ADCtrl$Basename, BDR_metadata_AD$Basename)]
  BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "Control",]$subtype = "Control"
  # BDR_ok_ADCtrl[is.na(BDR_ok_ADCtrl$subtype),]$subtype = "Undefined"
  table(BDR_ok_ADCtrl$subtype)
  
  UKBBN_newpheno$subtype = NA
  UKBBN_newpheno$subtype = UKBBN_AD$CVSpec_clusters[match(UKBBN_newpheno$Basename, UKBBN_AD$Basename)]
  UKBBN_newpheno[UKBBN_newpheno$AD == "0",]$subtype = "Control"
  table(UKBBN_newpheno$subtype)
  
  NIH_pheno$subtype = NA
  NIH_pheno$subtype = NIH_AD$CVSpec_clusters[match(NIH_pheno$Basename, NIH_AD$Basename)]
  NIH_pheno[NIH_pheno$Phenotype == "C",]$subtype = "Control"
  table(NIH_pheno$subtype)
  
  all_cohort_metadata$subtype = c(BDR_ok_ADCtrl$subtype,
                                  pheno$subtype,
                                  NIH_pheno$subtype)
  table(all_cohort_metadata$subtype)
  
  all_cohort_metadata$clinical_assessments = c(BDR_ok_ADCtrl$pathology,
                                               paste(UKBBN_newpheno$Pathology,
                                                     UKBBN_newpheno$Pathology2,
                                                     UKBBN_newpheno$Pathology3,
                                                     sep = " / "),
                                               paste(PITT_newpheno$NP.2000.Dx1,
                                                     PITT_newpheno$NP.2000.Dx2.1,
                                                     PITT_newpheno$NP.2000.Dx3,
                                                     PITT_newpheno$NP.2000.Dx4,
                                                     sep = " / "))
  
  all_cohort_metadata$BraakTangle = c(BDR_ok_ADCtrl$BraakTangle_numeric,
                                      UKBBN_newpheno$Braak_Tangle,
                                      PITT_newpheno$Braak.NFT.stage)
  all_cohort_metadata$BraakLB = c(BDR_ok_ADCtrl$Braak_LB,
                                  UKBBN_newpheno$BraakLB,
                                  PITT_newpheno$Lewy.Bodies)
  all_cohort_metadata[all_cohort_metadata$BraakLB == "",]$BraakLB = NA
  all_cohort_metadata$Cerad = c(BDR_ok_ADCtrl$Cerad,
                                UKBBN_newpheno$Cerad,
                                PITT_newpheno$CERAD.score..NPNEUR.)
  all_cohort_metadata[!is.na(all_cohort_metadata$Cerad) & all_cohort_metadata$Cerad == "",]$Cerad = NA
  all_cohort_metadata$MMSE = c(BDR_ok_ADCtrl$MMSE,
                               rep(NA, nrow(pheno)),
                               PITT_newpheno$Last.MMSE.score)
  all_cohort_metadata$CDR = c(BDR_ok_ADCtrl$CDR,
                              rep(NA, nrow(pheno)),
                              PITT_newpheno$Highest.CDR.score)
  all_cohort_metadata$Thal = c(BDR_ok_ADCtrl$Thal_number,
                               UKBBN_newpheno$Thal,
                               PITT_newpheno$Thal)
  all_cohort_metadata[!is.na(all_cohort_metadata$Thal) & all_cohort_metadata$Thal == "",]$Thal = NA
  all_cohort_metadata$APOE = c(BDR_ok_ADCtrl$APOE,
                               UKBBN_newpheno$APOE,
                               PITT_newpheno$ApoE.1)
  
  
  