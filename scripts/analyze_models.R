plot_model <- function(model, model_type, label_samples, outdir, numCC_clusters = 0){
  
  if(model_type == "Spectrum"){
    model_clusters = model$assignments
    model_matrix = model$similarity_matrix
  } else if(model_type == "CC"){
    model_clusters = model[[numCC_clusters]]$consensusClass
    model_matrix = model[[numCC_clusters]]$consensusMatrix
  }
  
  sil_model = silhouette(model_clusters, model_matrix)
  pdf(paste0(outdir, "silhouette_plot.pdf"), width = 7, height = 7)
  plot(sil_model,col = sil_model[,1][order(sil_model[,1])])
  dev.off()
  
  plot_clusters = model_heatmap(model_matrix = model_matrix,
                                model_clusters = model_clusters,
                                label = label_samples)
  
  pdf(paste0(outdir, "clusters_heatmap.pdf"), width = 7, height = 7)
  plot(plot_clusters)
  dev.off()
  
}

models_tables <- function(df_assignment, df_NMI_models, outdir){
  tables_res = list()
  
  for(i in c(1:length(df_assignment))){
    for(j in c(2:length(df_assignment))){
      if(i != j && i < j){
        
        table_res = tableGrob(table(df_assignment[[i]], df_assignment[[j]]), 
                              theme = ttheme_minimal(base_size = 7, padding = unit(c(4,2),"mm"),
                                rowhead=list(fg_params=list(col="black", fontface=2L,
                                                            hjust = 1.5))))
        
        table_res = gtable_add_grob(table_res,
                             grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                             t = 2, b = nrow(table_res), l = 1, r = ncol(table_res))
        table_res = gtable_add_grob(table_res,
                             grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                             t = 1, b = nrow(table_res), l = 2, r = ncol(table_res))
        
        h = grobHeight(table_res)
        w = grobWidth(table_res)
        
        
        table1_label = textGrob(names(df_assignment)[i], x=unit(0.5, "npc") - w*1, 
                                y=unit(0.5,"npc"), vjust=1, gp=gpar(fontsize=8, fontface = "bold"), rot = 90)
        table2_label = textGrob(names(df_assignment)[j], x=unit(0.5, "npc") + w*0.1, 
                                y=unit(0.5,"npc") + 1.05*h, vjust=0, gp=gpar(fontsize=8, fontface = "bold"))
        
        final_res = gTree(children = gList(table_res, table1_label, table2_label))
        tables_res[[length(tables_res)+1]] = final_res
      }
    }
  }
  
  plot_clusters = ggplot(data = df_NMI_models, aes(x = model_1, y = model_2, fill = nmi)) +
    ggtitle("NMI between models") +
    geom_tile() +
    geom_text(aes(label = round(nmi, 2), size = 13)) +
    scale_y_discrete(limits = rev(levels(factor(colnames(df_assignment), levels = colnames(df_assignment))))) +
    scale_x_discrete(limits = levels(factor(colnames(df_assignment), levels = colnames(df_assignment)))) +
    scale_fill_gradient2(low="#ffffb2", mid="#fd8d3c", high="firebrick", #colors in the scale
                         midpoint = 0.5, limits=c(0, 1)) +
    theme_bw() + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          panel.border = element_blank(),
          legend.position = "none")
  
  tables_res[[length(tables_res)+1]] = plot_clusters
  
  if(length(tables_res) == 7) {
    final_plot = arrangeGrob(
      grobs = tables_res,
      layout_matrix = rbind(c(7,7,7,1,4),
                            c(7,7,7,2,5),
                            c(7,7,7,3,6)),
      top = textGrob("Comparative results of clustering methods", gp = gpar(fontsize = 14))
    )
    ggsave(paste0(outdir, "Comparative_clustering.pdf"), final_plot, units = "px", width  = 3000, height = 1500)
  } else if(length(tables_res) == 11){
    final_plot = arrangeGrob(
      grobs = tables_res,
      layout_matrix = rbind(c(11,11,11,11,1,2),
                            c(11,11,11,11,3,4),
                            c(11,11,11,11,5,6),
                            c(11,11,11,11,7,8),
                            c(11,11,11,11,9,10)),
      
      top = textGrob("Comparative results of clustering methods", gp = gpar(fontsize = 16))
    )
    ggsave(paste0(outdir, "Comparative_clustering.pdf"), final_plot, units = "px", width  = 3000, height = 1500)
  }
}

calcNMI_all <- function(df_assignment){
  df_NMI_models = data.frame(matrix(ncol = 3))
  colnames(df_NMI_models) = c("model_1","model_2","nmi")
  
  if(ncol(df_assignment >= 2)){
    row_index = 1
    for(i in c(1:ncol(df_assignment))){
      for(j in c(i:ncol(df_assignment))){  
        df_NMI_models[row_index,] = c(colnames(df_assignment)[i],
                               colnames(df_assignment)[j],
                               calNMI(df_assignment[,i], df_assignment[,j]))
        row_index = row_index + 1
      }
    }
  }
  
  df_NMI_models$nmi = as.numeric(df_NMI_models$nmi)
  return(df_NMI_models)
}

plot_mixOmics_all <- function(model, outdir){
  
  pdf(paste0(outdir, "diablo_indiv_plotbackground12.pdf"), width = 7, height = 7)
  plotIndiv(model,ind.names = FALSE, comp = c(1,2), 
            legend=TRUE, ellipse = TRUE,
            title = 'Methyl') 
  dev.off()
  
  pdf(paste0(outdir, "diablo_indiv_plotbackground13.pdf"), width = 7, height = 7)
  plotIndiv(model,ind.names = FALSE, comp = c(1,3), 
            legend=TRUE, ellipse = TRUE,
            title = 'Methyl')
  dev.off()
  
  pdf(paste0(outdir, "diablo_indiv_plotbackground14.pdf"), width = 7, height = 7)
  plotIndiv(model,ind.names = FALSE, comp = c(1,4), 
            legend=TRUE, ellipse = TRUE,
            title = 'Methyl') 
  dev.off()
  
  pdf(paste0(outdir, "diablo_var_plot.pdf"), width = 7, height = 7)
  plotVar(model)
  dev.off()
  
  pdf(paste0(outdir, "diablo_loading_comp1.pdf"), width = 7, height = 7)
  plotLoadings(model, comp = 1, contrib = "max")
  dev.off()
  
  pdf(paste0(outdir, "diablo_loading_comp2.pdf"), width = 7, height = 7)
  plotLoadings(model, comp = 2, contrib = "max")
  dev.off()
  
  pdf(paste0(outdir, "diablo_loading_comp3.pdf"), width = 7, height = 7)
  plotLoadings(model, comp = 3, contrib = "max")
  dev.off()
  
  pdf(paste0(outdir, "diablo_loading_comp4.pdf"), width = 7, height = 7)
  plotLoadings(model, comp = 4, contrib = "max")
  dev.off()

  rocroc = auroc(model, roc.comp = 6)
  pdf(paste0(outdir, "diablo_roc_curves.pdf"), width = 7, height = 7)
  plot(rocroc$graph.Comp6 + theme_minimal())
  dev.off()
  
  
  colorscim = c("#388ECC","#F68B33","#C2C2C2","#009E73","#CC79A7","#F0E442")
  length(colorscim) = length(levels(model$Y))
  # set the styling of the legend to be homogeneous with previous plots
  legend=list(legend = levels(model$Y), # set of classes
              col = colorscim, # set of colours
              title = "Predicted subtypes", # legend title
              cex = 0.7) # legend size
  pdf(paste0(outdir, "diablo_heatmap.pdf"), width = 20, height = 17)
  # generate the CIM, using the legend and colouring rows by each sample's class
  
  cim(model, row.sideColors = color.mixo(model$Y), 
             legend = legend, mapping = "X")
  
  dev.off()
}

random_permutations <- function(data_matrix_AD, metadata_AD, diablo_Methyls){
  #Get original 1200 features stats
  HC_roc = auroc(diablo_Methyls[[1]], plot = FALSE, print = FALSE)
  KM_roc = auroc(diablo_Methyls[[2]], plot = FALSE, print = FALSE)
  PAM_roc = auroc(diablo_Methyls[[3]], plot = FALSE, print = FALSE)
  Spec_roc = auroc(diablo_Methyls[[4]], plot = FALSE, print = FALSE)
  final_aucs = c(mean(HC_roc$Comp6[,1]),
                 mean(KM_roc$Comp6[,1]),
                 mean(PAM_roc$Comp6[,1]),
                 mean(Spec_roc$Comp6[,1]))
  
  #Get top features
  HC_var = plotVar(diablo_Methyls[[1]], comp.select = c(1:6), plot = FALSE)
  HC_features = HC_var$names
  
  KM_var = plotVar(diablo_Methyls[[2]], comp.select = c(1:6), plot = FALSE)
  KM_features = KM_var$names
  
  PAM_var = plotVar(diablo_Methyls[[3]], comp.select = c(1:6), plot = FALSE)
  PAM_features = PAM_var$names
  
  Spec_var = plotVar(diablo_Methyls[[4]], comp.select = c(1:6), plot = FALSE)
  Spec_features = Spec_var$names
  
  df_shuffle = data.frame(matrix(nrow = 1, ncol = 12))
  colnames(df_shuffle) = c("HC_auc","HC_ratio", "HC_trueprop",
                           "KM_auc","KM_ratio", "KM_trueprop",
                           "PAM_auc","PAM_ratio", "PAM_trueprop",
                           "Spec_auc","Spec_ratio", "Spec_trueprop")
  #Generate random permutations of labeling from 100% random to 1% random
for(i in c(rep(seq(from = 1, to = 0.1, by = -0.1), each = 10), 
                             rep(seq(from = 0.09, to = 0.01, by = -0.01),each = 10))){
      #Randomise outcomes
      to_random = sample(seq(1:nrow(metadata_AD)), size = i*nrow(metadata_AD))
      shuffle_clust = sample(metadata_AD[to_random,]$HC_clusters)
      HC_shuffle = metadata_AD$HC_clusters
      HC_shuffle[to_random] = shuffle_clust
      df_random = data.frame(metadata_AD$HC_clusters, HC_shuffle)
      HC_prop = 1 - length(which(df_random[,1] == df_random[,2])) / nrow(metadata_AD)
      #Randomise outcomes
      to_random = sample(seq(1:nrow(metadata_AD)), size = i*nrow(metadata_AD))
      shuffle_clust = sample(metadata_AD[to_random,]$KM_clusters)
      KM_shuffle = metadata_AD$KM_clusters
      KM_shuffle[to_random] = shuffle_clust
      df_random = data.frame(metadata_AD$KM_clusters, KM_shuffle)
      KM_prop = 1 - length(which(df_random[,1] == df_random[,2])) / nrow(metadata_AD)
      #Randomise outcomes
      to_random = sample(seq(1:nrow(metadata_AD)), size = i*nrow(metadata_AD))
      shuffle_clust = sample(metadata_AD[to_random,]$PAM_clusters)
      PAM_shuffle = metadata_AD$PAM_clusters
      PAM_shuffle[to_random] = shuffle_clust
      df_random = data.frame(metadata_AD$PAM_clusters, PAM_shuffle)
      PAM_prop = 1 - length(which(df_random[,1] == df_random[,2])) / nrow(metadata_AD)
      #Randomise outcomes
      to_random = sample(seq(1:nrow(metadata_AD)), size = i*nrow(metadata_AD))
      shuffle_clust = sample(metadata_AD[to_random,]$Spec_clusters)
      Spec_shuffle = metadata_AD$Spec_clusters
      Spec_shuffle[to_random] = shuffle_clust
      df_random = data.frame(metadata_AD$Spec_clusters, Spec_shuffle)
      Spec_prop = 1 - length(which(df_random[,1] == df_random[,2])) / nrow(metadata_AD)
      #Make new models with random outcome
      diablo_shuffles = list()
      diablo_shuffles[[1]] = mixOmics::splsda(t(data_matrix_AD), HC_shuffle, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
      diablo_shuffles[[2]] = mixOmics::splsda(t(data_matrix_AD), KM_shuffle, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
      diablo_shuffles[[3]] = mixOmics::splsda(t(data_matrix_AD), PAM_shuffle, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
      diablo_shuffles[[4]] = mixOmics::splsda(t(data_matrix_AD), Spec_shuffle, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
      #Calculate stats
      HC_roc = auroc(diablo_shuffles[[1]], plot = FALSE, print = FALSE)
      KM_roc = auroc(diablo_shuffles[[2]], plot = FALSE, print = FALSE)
      PAM_roc = auroc(diablo_shuffles[[3]], plot = FALSE, print = FALSE)
      Spec_roc = auroc(diablo_shuffles[[4]], plot = FALSE, print = FALSE)
      
      #Extract top features of random models
      HC_shuffle_var = plotVar(diablo_shuffles[[1]], comp.select = c(1:6), plot = FALSE)
      HC_shuffle_features = HC_shuffle_var$names
      
      KM_shuffle_var = plotVar(diablo_shuffles[[2]], comp.select = c(1:6), plot = FALSE)
      KM_shuffle_features = KM_shuffle_var$names
      
      PAM_shuffle_var = plotVar(diablo_shuffles[[3]], comp.select = c(1:6), plot = FALSE)
      PAM_shuffle_features = PAM_shuffle_var$names
      
      Spec_shuffle_var = plotVar(diablo_shuffles[[4]], comp.select = c(1:6), plot = FALSE)
      Spec_shuffle_features = Spec_shuffle_var$names
      
      #Intersect of features between original model and randomized model
      HC_inter = intersect(HC_features, HC_shuffle_features)
      KM_inter = intersect(KM_features, KM_shuffle_features)
      PAM_inter = intersect(PAM_features, PAM_shuffle_features)
      Spec_inter = intersect(Spec_features, Spec_shuffle_features)
      
      shuffle_res = c(mean(HC_roc$Comp4[,1]), length(HC_inter)/length(HC_features), HC_prop,
                      mean(KM_roc$Comp4[,1]), length(KM_inter)/length(KM_features), KM_prop,
                      mean(PAM_roc$Comp4[,1]), length(PAM_inter)/length(PAM_features), PAM_prop,
                      mean(Spec_roc$Comp4[,1]), length(Spec_inter)/length(Spec_features), Spec_prop)
      df_shuffle = rbind(df_shuffle,shuffle_res)
  }
  return(df_shuffle) 
}

output_permutation <- function(df_shuffle, dirname){
  
  # Output models of accuracy vs intersect of features, useless?
  ggplot(data = df_shuffle, aes(x = HC_accuracy, y = HC_intersect)) +
    geom_point() +
    xlim(0,1) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/",dirname,"/HClust/HC_acc.jpeg"))
  ggplot(data = df_shuffle, aes(x = KM_accuracy, y = KM_intersect)) +
    geom_point() +
    xlim(0,1) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/",dirname,"/Kmeans/KM_acc.jpeg"))
  ggplot(data = df_shuffle, aes(x = PAM_accuracy, y = PAM_intersect)) +
    geom_point() +
    xlim(0,1) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/",dirname,"/PAM/PAM_acc.jpeg"))
  ggplot(data = df_shuffle, aes(x = Spec_accuracy, y = Spec_intersect)) +
    geom_point() +
    xlim(0,1) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/",dirname,"/Spectrum/Spec_acc.jpeg"))
  
  # Output models of proportion of randomized outcome vs prop of shared features
  df_shuffle$prop = rep(c(seq(from = 1, to = 0.1, by = -0.1), seq(from = 0.09, to = 0.01, by = -0.01)), each = 10)
  
  ggplot(data = df_shuffle, aes(x = prop, y = HC_intersect)) +
    geom_point() +
    xlim(1,0) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/", dirname, "/HClust/HC_prop.jpeg"))
  ggplot(data = df_shuffle, aes(x = prop, y = KM_intersect)) +
    geom_point() +
    xlim(1,0) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/", dirname, "/Kmeans/KM_prop.jpeg"))
  ggplot(data = df_shuffle, aes(x = prop, y = PAM_intersect)) +
    geom_point() +
    xlim(1,0) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/", dirname, "/PAM/PAM_prop.jpeg"))
  ggplot(data = df_shuffle, aes(x = prop, y = Spec_intersect)) +
    geom_point() +
    xlim(1,0) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/", dirname, "/Spectrum/Spec_prop.jpeg"))
}

cv_diablo <- function(data_matrix_AD, metadata_AD){
  #Create 10folds
  flds <- createFolds(seq(1:nrow(metadata_AD)), k = 10, list = TRUE, returnTrain = FALSE)
  
  #Create tables to hold all the top features up until 1200
  HC_features_all = data.frame(matrix(NA, nrow = 1200, ncol = 10))
  KM_features_all = data.frame(matrix(NA, nrow = 1200, ncol = 10))
  PAM_features_all = data.frame(matrix(NA, nrow = 1200, ncol = 10))
  Spec_features_all = data.frame(matrix(NA, nrow = 1200, ncol = 10))
  
  for(i in seq(1:10)){
    #Select 90% samples
    sample_set = flds[seq(1:10)[-i]]
    sample_set = unlist(sample_set)
    sub_betas_AD = data_matrix_AD[,sample_set]
    sub_clusters = metadata_AD[sample_set,]
    #Build CV models
    HC_val = mixOmics::splsda(t(sub_betas_AD), sub_clusters$HC_clusters, keepX=c(200,200,200,200,200,200), ncomp = 6)
    KM_val = mixOmics::splsda(t(sub_betas_AD), sub_clusters$KM_clusters, keepX=c(200,200,200,200,200,200), ncomp = 6)
    PAM_val = mixOmics::splsda(t(sub_betas_AD), sub_clusters$PAM_clusters, keepX=c(200,200,200,200,200,200), ncomp = 6)
    Spec_val = mixOmics::splsda(t(sub_betas_AD), sub_clusters$Spec_clusters, keepX=c(200,200,200,200,200,200), ncomp = 6)
    
    #Select all top features from the CV
    HC_val_var = plotVar(HC_val, comp.select = c(1:6))
    HC_val_features = HC_val_var$names
    
    KM_val_var = plotVar(KM_val, comp.select = c(1:6))
    KM_val_features = KM_val_var$names
    
    PAM_val_var = plotVar(PAM_val, comp.select = c(1:6))
    PAM_val_features = PAM_val_var$names
    
    Spec_val_var = plotVar(Spec_val, comp.select = c(1:6))
    Spec_val_features = Spec_val_var$names
    
    #Make sure that length is consistent with table
    length(HC_val_features) = 1200
    length(KM_val_features) = 1200
    length(PAM_val_features) = 1200
    length(Spec_val_features) = 1200
    
    HC_features_all[,i] = HC_val_features
    KM_features_all[,i] = KM_val_features
    PAM_features_all[,i] = PAM_val_features
    Spec_features_all[,i] = Spec_val_features
    
  }
  
  return(list(HC_features_all,
              KM_features_all,
              PAM_features_all,
              Spec_features_all))
}

build_crossmethod <- function(metadata_AD, df_assignment){
  metadata_AD$all_clust = NA
  n_profile = 0
  for(i in seq(1:nrow(metadata_AD))){#Generate a single 4-outcome profile for each sample
    tmp_profile = paste(c(metadata_AD[i,]$HC_clusters,
                          metadata_AD[i,]$KM_clusters,
                          metadata_AD[i,]$PAM_clusters,
                          metadata_AD[i,]$Spec_clusters),
                        collapse = "-")
    metadata_AD[i,]$all_clust = tmp_profile
  }
  #Checking results
  table(metadata_AD$all_clust)[order(table(metadata_AD$all_clust), decreasing = T)]
  
  df_mat_assignment = data.frame(matrix(0, nrow = nrow(df_assignment), ncol = nrow(df_assignment)))
  rownames(df_mat_assignment) = rownames(df_assignment)
  colnames(df_mat_assignment) = rownames(df_assignment)
  
  for(i in seq(1,nrow(df_assignment)-1)){# Compute how much outcomes from models are shared
    for(j in seq(i+1,nrow(df_assignment))){
      
      which(df_assignment[i,] == df_assignment[j,])
      
      sum(df_assignment[i,] == df_assignment[j,])
      
      df_mat_assignment[i,j] = sum(df_assignment[i,] == df_assignment[j,])
      df_mat_assignment[j,i] = sum(df_assignment[i,] == df_assignment[j,])
    }
  }
  
  idx = 1
  txt_idx = paste("Group", idx)
  for(i in seq(1,nrow(df_assignment))){#Assign a distinct group number for each 4-outcome
    if(length(which(df_mat_assignment[i,] == 4)) > 0 & df_mat_assignment[i,i] == 0){
      for(j in c(i, which(df_mat_assignment[i,] == 4))){
        df_mat_assignment[j,j] = txt_idx
      }
      
      idx = idx + 1
      txt_idx = paste("Group", idx)
      
    }
  }
  #Assign groups to each sample in diagonal for double checking
  metadata_AD$groups = diag(as.matrix(df_mat_assignment))
  
  table(diag(as.matrix(df_mat_assignment)))[order(table(diag(as.matrix(df_mat_assignment))), decreasing = T)]
  df_mat_assignment_order = as.matrix(df_mat_assignment[order(diag(as.matrix(df_mat_assignment))),order(diag(as.matrix(df_mat_assignment)))])
  #Build a table of cross-method clusters
  snames_mat = rownames(df_mat_assignment_order)
  groups = as.data.frame(diag(as.matrix(df_mat_assignment_order)))
  colnames(groups) = "groups"
  
  return(list(metadata_AD, df_mat_assignment_order, groups, snames_mat))
  
}

cv_crossmethod <- function(data_matrix_AD, metadata_AD, groups, dirname, dirname_final){
  # Remove excluded samples
  to_mixo = groups[groups$Clusters != "Excluded",]
  
  to_mixo$Nclusters = gsub("Cluster ", "", to_mixo$Clusters)
  to_mixo$Nclusters = as.numeric(to_mixo$Nclusters)
  
  #Select all betas for selected samples
  cross_betas = data_matrix_AD[,colnames(data_matrix_AD) %in% to_mixo$Basename]
  
  # 1200 features selected for crossmethod outcome
  Cross_methyl = mixOmics::splsda(t(cross_betas), to_mixo$Nclusters, keepX=c(200,200,200,200,200,200), ncomp = 6)

  outdir = paste0("results/", dirname, "/cross/")
  dir.create(outdir)
  plot_mixOmics_all(model = Cross_methyl, outdir = outdir)
  
  #10-fold Cross validating as before for single method outcome
  flds <- createFolds(seq(1:nrow(to_mixo)), k = 10, list = TRUE, returnTrain = FALSE)
  
  Cross_features_all = data.frame(matrix(NA, nrow = 1200, ncol = 10))
  
  for(i in seq(1:10)){
    
    sample_set = flds[seq(1:10)[-i]]
    sample_set = unlist(sample_set)
    sub_cross_betas = cross_betas[,sample_set]
    sub_clusters = to_mixo[sample_set,]
    Cross_val = mixOmics::splsda(t(sub_cross_betas), sub_clusters$Nclusters, keepX=c(200,200,200,200,200,200), ncomp = 6)
    
    Cross_val_var = plotVar(Cross_val, comp.select = c(1:6), plot = FALSE)
    Cross_val_features = Cross_val_var$names
    
    length(Cross_val_features) = 1200
    
    Cross_features_all[,i] = Cross_val_features
    
    
  }
  # Selecting shared features
  Cross_final_features = Reduce(intersect, list(Cross_features_all[,1], Cross_features_all[,2], Cross_features_all[,3], Cross_features_all[,4],
                                                Cross_features_all[,5], Cross_features_all[,6], Cross_features_all[,7], Cross_features_all[,8], 
                                                Cross_features_all[,9], Cross_features_all[,10]))
  
  Cross_mSet = cross_betas[rownames(cross_betas) %in% Cross_final_features,]
  # Build CV crossmethod model
  Cross_comp = min(200, length(Cross_final_features))
  Cross_final = mixOmics::splsda(t(Cross_mSet), to_mixo$Nclusters, keepX=c(Cross_comp,Cross_comp,Cross_comp,
                                                                           Cross_comp,Cross_comp,Cross_comp), ncomp = 6)
  
  outdir = paste0("results/", dirname_final, "/cross/")
  dir.create(outdir)
  plot_mixOmics_all(model = Cross_final, outdir = outdir)
  
  dd = plotIndiv(Cross_final,ind.names = FALSE, comp = c(1,2), 
                 legend=TRUE, ellipse = TRUE,
                 title = 'Methyl', plot = FALSE) 
  df_cor_clusters = data.frame(dd$df$group, dd$df$x, dd$df$y)
  dd = plotIndiv(Cross_final,ind.names = FALSE, comp = c(3,4), 
                 legend=TRUE, ellipse = TRUE,
                 title = 'Methyl', plot = FALSE) 
  df_cor_clusters = cbind(df_cor_clusters, data.frame(dd$df$x, dd$df$y))
  dd = plotIndiv(Cross_final,ind.names = FALSE, comp = c(5,6), 
                 legend=TRUE, ellipse = TRUE,
                 title = 'Methyl', plot = FALSE) 
  df_cor_clusters = cbind(df_cor_clusters, data.frame(dd$df$x, dd$df$y))
  
  colnames(df_cor_clusters) = c("Cluster", "Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Comp6")
  rownames(df_cor_clusters) = rownames(to_mixo)
  
  df_cor_clusters$Cluster = as.numeric(df_cor_clusters$Cluster)
  
  cross_metadata = metadata_AD[metadata_AD$Basename %in% to_mixo$Basename,]
  
  cross_metadata = cbind(cross_metadata, df_cor_clusters)
  
  return(list(Cross_methyl, Cross_final, Cross_final_features, cross_metadata))
}

projection_to_rep <- function(model, data_matrix_AD, rep_betas, rep_metadata, nclusters, discovery, replication, outdir){
  
  # Select the correct samples
  if(replication == "UKBBN"){
    pred_cohort = rep_metadata[rep_metadata$AD == 1,]
    pred_betas = rep_betas[rownames(rep_betas) %in% rownames(data_matrix_AD), colnames(rep_betas) %in% pred_cohort$Row.names]
  } else if(replication == "BDR"){
    pred_cohort = rep_metadata
    pred_betas = rep_betas[rownames(rep_betas) %in% rownames(data_matrix_AD),]
  } else if(replication == "NIH"){
    rep_metadata$Phenotype = as.character(rep_metadata$Phenotype)
    pred_cohort = rep_metadata[rep_metadata$Phenotype %in% c("AD+P","AD-P"),]
    pred_betas = rep_betas[rownames(rep_betas) %in% rownames(data_matrix_AD), colnames(rep_betas) %in% pred_cohort$Basename]
  } 
  
  diablo_pred = predict(model, newdata = t(pred_betas)) 
  
  # Clusters are attributed according to the most probable outcome
  if(nclusters == 3){
  res_pred = data.frame("1" = rowSums(diablo_pred[["predict"]][,1,]),
                        "2" = rowSums(diablo_pred[["predict"]][,2,]),
                        "3" = rowSums(diablo_pred[["predict"]][,3,]))
  res_pred$class = apply(res_pred ,1, which.max)
  } else if(nclusters == 4){
    res_pred = data.frame("1" = rowSums(diablo_pred[["predict"]][,1,]),
                          "2" = rowSums(diablo_pred[["predict"]][,2,]),
                          "3" = rowSums(diablo_pred[["predict"]][,3,]),
                          "4" = rowSums(diablo_pred[["predict"]][,4,]))
    res_pred$class = apply(res_pred ,1, which.max)
  }
  
  # Exporting figure and result
  pdf(paste0("./results/", outdir, "/", replication, "projection_indiv.pdf"), width = 7, height = 7)
  plotIndiv(model, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE)
  points(diablo_pred$variates[, 1], diablo_pred$variates[, 2], pch = 19, cex = 1.2, col = res_pred$class)
  dev.off()
  tmp = plotLoadings(model, comp = 1, contrib = "max", plot = FALSE)
  tmp = rbind(tmp, plotLoadings(model, comp = 2, contrib = "max", plot = FALSE))
  tmp = rbind(tmp, plotLoadings(model, comp = 3, contrib = "max", plot = FALSE))
  tmp = rbind(tmp, plotLoadings(model, comp = 4, contrib = "max", plot = FALSE))
  tmp = rbind(tmp, plotLoadings(model, comp = 5, contrib = "max", plot = FALSE))
  tmp = rbind(tmp, plotLoadings(model, comp = 6, contrib = "max", plot = FALSE))
  
  return(list(tmp, res_pred))
}

cpg_stats_BDR <- function(metadata_AD, metadata_ctrl, mSet_betas, top_cpg){
  
  # Very questionable loop ?
  all_cpgs = data.frame(matrix(NA, nrow = 1200, ncol = 40))
  col_idx = 1
  max_val = 0
  for(model_n in seq(1:10)){
    print(model_n)
    for(clust_n in unique(top_cpg[[model_n]]$GroupContrib)[order(unique(top_cpg[[model_n]]$GroupContrib))]){
      print(clust_n)
      colname = paste0("model", model_n, ";clust", clust_n)
      tmp_vec = rownames(top_cpg[[model_n]][top_cpg[[model_n]]$GroupContrib == clust_n,])
      
      if(length(tmp_vec) > max_val){
        max_val = length(tmp_vec)
      }
      
      length(tmp_vec) = 1200
      
      all_cpgs[,col_idx] = tmp_vec
      colnames(all_cpgs)[col_idx] = colname
      
      col_idx = col_idx + 1
      
    }
  }
  
  all_cpgs = all_cpgs[1:max_val, 1:col_idx-1]
  allclust = data.frame(metadata_AD$Basename,
                            metadata_AD$HC_clusters,
                            metadata_AD$KM_clusters,
                            metadata_AD$PAM_clusters,
                            metadata_AD$Spec_clusters,
                            metadata_AD$cross_clusters)
  colnames(allclust) = c("Basename","HC_clusters","KM_clusters","PAM_clusters","Spec_clusters","cross_clusters")
  Ctrl_clust = data.frame("Basename" = metadata_ctrl$Basename,
                              "Outcome" = rep(99, nrow(metadata_ctrl)))
  

  df_res_lm = data.frame(matrix(nrow = 1, ncol = 7))
  colnames(df_res_lm) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
  df_res_lm_final = data.frame(matrix(nrow = 1, ncol = 7))
  colnames(df_res_lm_final) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
  #For each model select the appropriate samples and model/outcome and run the EWAS
  for(n_model in seq(1:ncol(all_cpgs))){
    output = colnames(all_cpgs)[n_model]
    output = unlist(strsplit(output, split = ";"))
    output = gsub("model", "", output)
    output = as.numeric(gsub("clust", "", output))
    
    print(output)
    if(output[1] %in% c(1,6)){
      list_samples = allclust[allclust$HC_clusters == output[2],]$Basename
    } else if(output[1] %in% c(2,7)){
      list_samples = allclust[allclust$KM_clusters == output[2],]$Basename
    } else if(output[1] %in% c(3,8)){
      list_samples = allclust[allclust$PAM_clusters == output[2],]$Basename
    } else if(output[1] %in% c(4,9)){
      list_samples = allclust[allclust$Spec_clusters == output[2],]$Basename
    } else if(output[1] %in% c(5,10)){
      list_samples = allclust[allclust$cross_clusters == output[2],]$Basename
      list_samples = list_samples[!is.na(list_samples)]
    }
    
    df_tolm = data.frame("Basename" = list_samples,
                         "Outcome" = output[2])
    df_tolm = rbind(df_tolm, Ctrl_clust)
    df_tolm$Outcome = as.numeric(df_tolm$Outcome)
    print(df_tolm)
    
    betas_tolm = as.data.frame(t(mSet_betas[,colnames(mSet_betas) %in% df_tolm$Basename]))
    betas_tolm = cbind(df_tolm$Outcome, betas_tolm)
    colnames(betas_tolm)[1] = "Outcome"
    
    cores=detectCores()
    cl <- makeCluster(cores[1]-2) #not to overload your computer
    registerDoParallel(cl)
    clusterExport(cl, "cohensD")
    
    df_res_lm <- foreach(ncpg=c(2:ncol(betas_tolm)), .combine=rbind) %dopar% {
      input_betas = betas_tolm[,c(1,ncpg)]
      
      cpg_name = colnames(input_betas)[2]
      colnames(input_betas)[2] = "predictor"
      
      res_lm = lm(formula = Outcome ~ predictor, data = input_betas)
      sum_lm = summary(res_lm)
      cohd = cohensD(formula = predictor ~ Outcome, data = input_betas)
      
      new_line = c(output[1],output[2],cpg_name, sum_lm$coefficients[2,4], sum_lm$r.squared, res_lm$effects[2], cohd)
      
      new_line #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
    }
    stopCluster(cl)
    colnames(df_res_lm) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
    df_res_lm_final = rbind(df_res_lm_final, df_res_lm)
    
  }
  
  df_res_lm_final = as.data.frame(df_res_lm_final)
  
  df_res_lm_final = df_res_lm_final[-1,]
  colnames(df_res_lm_final) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
  df_res_lm_final[df_res_lm_final$model == "1",]$model = "HClust"
  df_res_lm_final[df_res_lm_final$model == "2",]$model = "KMeans"
  df_res_lm_final[df_res_lm_final$model == "3",]$model = "PAM"
  df_res_lm_final[df_res_lm_final$model == "4",]$model = "Spectrum"
  df_res_lm_final[df_res_lm_final$model == "5",]$model = "Crossmethod"
  df_res_lm_final[df_res_lm_final$model == "6",]$model = "CV_HClust"
  df_res_lm_final[df_res_lm_final$model == "7",]$model = "CV_KMeans"
  df_res_lm_final[df_res_lm_final$model == "8",]$model = "CV_PAM"
  df_res_lm_final[df_res_lm_final$model == "9",]$model = "CV_Spectrum"
  df_res_lm_final[df_res_lm_final$model == "10",]$model = "CV_Crossmethod"
  
  return(df_res_lm_final)
}

cpg_stats_UKBBN <- function(metadata_AD, metadata_ctrl, mSet_betas, top_cpg){
  all_cpgs = data.frame(matrix(NA, nrow = 1200, ncol = 40))
  col_idx = 1
  max_val = 0
  for(model_n in seq(1:10)){
    print(model_n)
    for(clust_n in unique(top_cpg[[model_n]]$GroupContrib)[order(unique(top_cpg[[model_n]]$GroupContrib))]){
      print(clust_n)
      colname = paste0("model", model_n, ";clust", clust_n)
      tmp_vec = rownames(top_cpg[[model_n]][top_cpg[[model_n]]$GroupContrib == clust_n,])
      
      if(length(tmp_vec) > max_val){
        max_val = length(tmp_vec)
      }
      
      length(tmp_vec) = 1200
      
      all_cpgs[,col_idx] = tmp_vec
      colnames(all_cpgs)[col_idx] = colname
      
      col_idx = col_idx + 1
      
    }
  }
  
  all_cpgs = all_cpgs[1:max_val, 1:col_idx-1]
  
  allclust = data.frame(metadata_AD$Basename,
                        metadata_AD$HC_clusters,
                        metadata_AD$KM_clusters,
                        metadata_AD$PAM_clusters,
                        metadata_AD$Spec_clusters,
                        metadata_AD$cross_clusters,
                        metadata_AD$CVHC_clusters,
                        metadata_AD$CVKM_clusters,
                        metadata_AD$CVPAM_clusters,
                        metadata_AD$CVSpec_clusters,
                        metadata_AD$CVcross_clusters)
  colnames(allclust) = c("Basename","HC_clusters","KM_clusters","PAM_clusters","Spec_clusters","cross_clusters",
                         "CVHC_clusters","CVKM_clusters","CVPAM_clusters","CVSpec_clusters","CVcross_clusters")
  Ctrl_clust = data.frame("Basename" = metadata_ctrl$Basename,
                          "Outcome" = rep(99, nrow(metadata_ctrl)))
  
  
  df_res_lm = data.frame(matrix(nrow = 1, ncol = 7))
  colnames(df_res_lm) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
  df_res_lm_final = data.frame(matrix(nrow = 1, ncol = 7))
  colnames(df_res_lm_final) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
  
  for(n_model in seq(1:ncol(all_cpgs))){
    output = colnames(all_cpgs)[n_model]
    output = unlist(strsplit(output, split = ";"))
    output = gsub("model", "", output)
    output = as.numeric(gsub("clust", "", output))
    
    print(output)
    if(output[1] == 1){
      list_samples = allclust[allclust$HC_clusters == output[2],]$Basename
    } else if(output[1] == 2){
      list_samples = allclust[allclust$KM_clusters == output[2],]$Basename
    } else if(output[1] == 3){
      list_samples = allclust[allclust$PAM_clusters == output[2],]$Basename
    } else if(output[1] == 4){
      list_samples = allclust[allclust$Spec_clusters == output[2],]$Basename
    } else if(output[1] == 5){
      list_samples = allclust[allclust$cross_clusters == output[2],]$Basename
    } else if(output[1] == 6){
      list_samples = allclust[allclust$CVHC_clusters == output[2],]$Basename
    } else if(output[1] == 7){
      list_samples = allclust[allclust$CVKM_clusters == output[2],]$Basename
    } else if(output[1] == 8){
      list_samples = allclust[allclust$CVPAM_clusters == output[2],]$Basename
    } else if(output[1] == 9){
      list_samples = allclust[allclust$CVSpec_clusters == output[2],]$Basename
    } else if(output[1] == 10){
      list_samples = allclust[allclust$CVcross_clusters == output[2],]$Basename
    }
    
    if(length(list_samples) > 0){
      df_tolm = data.frame("Basename" = list_samples,
                           "Outcome" = output[2])
      df_tolm = rbind(df_tolm, Ctrl_clust)
      df_tolm$Outcome = as.numeric(df_tolm$Outcome)
      print(df_tolm)
      
      betas_tolm = as.data.frame(t(mSet_betas[,colnames(mSet_betas) %in% df_tolm$Basename]))
      betas_tolm = cbind(df_tolm$Outcome, betas_tolm)
      colnames(betas_tolm)[1] = "Outcome"
      
      cores=detectCores()
      cl <- makeCluster(cores[1]-2) #not to overload your computer
      registerDoParallel(cl)
      clusterExport(cl, "cohensD")
      
      
      df_res_lm <- foreach(ncpg=c(2:ncol(betas_tolm)), .combine=rbind) %dopar% {
        input_betas = betas_tolm[,c(1,ncpg)]
        
        cpg_name = colnames(input_betas)[2]
        colnames(input_betas)[2] = "predictor"
        
        res_lm = lm(formula = Outcome ~ predictor, data = input_betas)
        sum_lm = summary(res_lm)
        cohd = cohensD(formula = predictor ~ Outcome, data = input_betas)
        
        new_line = c(output[1],output[2],cpg_name, sum_lm$coefficients[2,4], sum_lm$r.squared, res_lm$effects[2], cohd)
        
        new_line #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
      }
      stopCluster(cl)
      colnames(df_res_lm) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
      df_res_lm_final = rbind(df_res_lm_final, df_res_lm)
    }
    
  }
  
  df_res_lm_final = as.data.frame(df_res_lm_final)
  df_res_lm_final = df_res_lm_final[-1,]
  colnames(df_res_lm_final) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
  df_res_lm_final[df_res_lm_final$model == "1",]$model = "HClust"
  df_res_lm_final[df_res_lm_final$model == "2",]$model = "KMeans"
  df_res_lm_final[df_res_lm_final$model == "3",]$model = "PAM"
  df_res_lm_final[df_res_lm_final$model == "4",]$model = "Spectrum"
  df_res_lm_final[df_res_lm_final$model == "5",]$model = "Crossmethod"
  df_res_lm_final[df_res_lm_final$model == "6",]$model = "CV_HClust"
  df_res_lm_final[df_res_lm_final$model == "7",]$model = "CV_KMeans"
  df_res_lm_final[df_res_lm_final$model == "8",]$model = "CV_PAM"
  df_res_lm_final[df_res_lm_final$model == "9",]$model = "CV_Spectrum"
  df_res_lm_final[df_res_lm_final$model == "10",]$model = "CV_Crossmethod"
  
  return(df_res_lm_final)
}

cpg_stats_BDR_s1vs3 <- function(metadata_AD, mSet_betas, top_cpg){
  all_cpgs = data.frame(matrix(NA, nrow = 1200, ncol = 40))
  col_idx = 1
  max_val = 0
  for(model_n in seq(1:10)){
    print(model_n)
    for(clust_n in unique(top_cpg[[model_n]]$GroupContrib)[order(unique(top_cpg[[model_n]]$GroupContrib))]){
      print(clust_n)
      colname = paste0("model", model_n, ";clust", clust_n)
      tmp_vec = rownames(top_cpg[[model_n]][top_cpg[[model_n]]$GroupContrib == clust_n,])
      
      if(length(tmp_vec) > max_val){
        max_val = length(tmp_vec)
      }
      
      length(tmp_vec) = 1200
      
      all_cpgs[,col_idx] = tmp_vec
      colnames(all_cpgs)[col_idx] = colname
      
      col_idx = col_idx + 1
      
    }
  }
  
  all_cpgs = all_cpgs[1:max_val, 1:col_idx-1]
  allclust = data.frame(metadata_AD$Basename,
                        metadata_AD$HC_clusters,
                        metadata_AD$KM_clusters,
                        metadata_AD$PAM_clusters,
                        metadata_AD$Spec_clusters,
                        metadata_AD$cross_clusters)
  colnames(allclust) = c("Basename","HC_clusters","KM_clusters","PAM_clusters","Spec_clusters","cross_clusters")
  
  df_res_lm = data.frame(matrix(nrow = 1, ncol = 7))
  colnames(df_res_lm) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
  df_res_lm_final = data.frame(matrix(nrow = 1, ncol = 7))
  colnames(df_res_lm_final) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
  all_cpgs = all_cpgs[,grep(";clust1",colnames(all_cpgs))]
  for(n_model in seq(1:ncol(all_cpgs))){
    output = colnames(all_cpgs)[n_model]
    output = unlist(strsplit(output, split = ";"))
    output = gsub("model", "", output)
    output = as.numeric(gsub("clust", "", output))
    
    print(output)
    if(output[1] %in% c(1,6)){
      list_samples = allclust[allclust$HC_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,2)]
    } else if(output[1] %in% c(2,7)){
      list_samples = allclust[allclust$KM_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,3)]
    } else if(output[1] %in% c(3,8)){
      list_samples = allclust[allclust$PAM_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,4)]
    } else if(output[1] %in% c(4,9)){
      list_samples = allclust[allclust$Spec_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,5)]
    } else if(output[1] %in% c(5,10)){
      list_samples = allclust[allclust$cross_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,6)]
    }
    
    df_tolm = data.frame("Basename" = list_samples[,1],
                         "Outcome" = list_samples[,2])
    df_tolm$Outcome = as.numeric(df_tolm$Outcome)
    print(df_tolm)
  
    betas_tolm = as.data.frame(t(mSet_betas[,colnames(mSet_betas) %in% df_tolm$Basename]))
    betas_tolm = cbind(df_tolm$Outcome, betas_tolm)
    colnames(betas_tolm)[1] = "Outcome"
    
    cores=detectCores()
    cl <- makeCluster(cores[1]-2) #not to overload your computer
    registerDoParallel(cl)
    clusterExport(cl, "cohensD")
    
    df_res_lm <- foreach(ncpg=c(2:ncol(betas_tolm)), .combine=rbind) %dopar% {
      input_betas = betas_tolm[,c(1,ncpg)]
      
      cpg_name = colnames(input_betas)[2]
      colnames(input_betas)[2] = "predictor"
      
      res_lm = lm(formula = Outcome ~ predictor, data = input_betas)
      sum_lm = summary(res_lm)
      cohd = cohensD(formula = predictor ~ Outcome, data = input_betas)
      
      new_line = c(output[1],output[2],cpg_name, sum_lm$coefficients[2,4], sum_lm$r.squared, res_lm$effects[2], cohd)
      
      new_line #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
    }
    stopCluster(cl)
    colnames(df_res_lm) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
    df_res_lm_final = rbind(df_res_lm_final, df_res_lm)
    
  }
  
  df_res_lm_final = as.data.frame(df_res_lm_final)
  
  df_res_lm_final = df_res_lm_final[-1,]
  colnames(df_res_lm_final) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
  df_res_lm_final[df_res_lm_final$model == "1",]$model = "HClust"
  df_res_lm_final[df_res_lm_final$model == "2",]$model = "KMeans"
  df_res_lm_final[df_res_lm_final$model == "3",]$model = "PAM"
  df_res_lm_final[df_res_lm_final$model == "4",]$model = "Spectrum"
  df_res_lm_final[df_res_lm_final$model == "5",]$model = "Crossmethod"
  df_res_lm_final[df_res_lm_final$model == "6",]$model = "CV_HClust"
  df_res_lm_final[df_res_lm_final$model == "7",]$model = "CV_KMeans"
  df_res_lm_final[df_res_lm_final$model == "8",]$model = "CV_PAM"
  df_res_lm_final[df_res_lm_final$model == "9",]$model = "CV_Spectrum"
  df_res_lm_final[df_res_lm_final$model == "10",]$model = "CV_Crossmethod"
  
  return(df_res_lm_final)
}



cpg_stats_UKBBN_s1vs3 <- function(metadata_AD, betas_AD, top){
  all_cpgs = data.frame(matrix(NA, nrow = 1200, ncol = 40))
  col_idx = 1
  max_val = 0
  for(model_n in seq(1:10)){
    print(model_n)
    for(clust_n in unique(top[[model_n]]$GroupContrib)[order(unique(top[[model_n]]$GroupContrib))]){
      print(clust_n)
      colname = paste0("model", model_n, ";clust", clust_n)
      tmp_vec = rownames(top[[model_n]][top[[model_n]]$GroupContrib == clust_n,])
      
      if(length(tmp_vec) > max_val){
        max_val = length(tmp_vec)
      }
      
      length(tmp_vec) = 1200
      
      all_cpgs[,col_idx] = tmp_vec
      colnames(all_cpgs)[col_idx] = colname
      
      col_idx = col_idx + 1
      
    }
  }
  
  all_cpgs = all_cpgs[1:max_val, 1:col_idx-1]
  
  allclust = data.frame(metadata_AD$Basename,
                        metadata_AD$HC_clusters,
                        metadata_AD$KM_clusters,
                        metadata_AD$PAM_clusters,
                        metadata_AD$Spec_clusters,
                        metadata_AD$cross_clusters,
                        metadata_AD$CVHC_clusters,
                        metadata_AD$CVKM_clusters,
                        metadata_AD$CVPAM_clusters,
                        metadata_AD$CVSpec_clusters,
                        metadata_AD$CVcross_clusters)
  colnames(allclust) = c("Basename","HC_clusters","KM_clusters","PAM_clusters","Spec_clusters","cross_clusters",
                         "CVHC_clusters","CVKM_clusters","CVPAM_clusters","CVSpec_clusters","CVcross_clusters")
  
  
  df_res_lm = data.frame(matrix(nrow = 1, ncol = 7))
  colnames(df_res_lm) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
  df_res_lm_final = data.frame(matrix(nrow = 1, ncol = 7))
  colnames(df_res_lm_final) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
  all_cpgs = all_cpgs[,grep(";clust1",colnames(all_cpgs))]
  for(n_model in seq(1:ncol(all_cpgs))){
    output = colnames(all_cpgs)[n_model]
    output = unlist(strsplit(output, split = ";"))
    output = gsub("model", "", output)
    output = as.numeric(gsub("clust", "", output))
    
    print(output)
    if(output[1] == 1){
      list_samples = allclust[allclust$HC_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,2)]
    } else if(output[1] == 2){
      list_samples = allclust[allclust$KM_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,3)]
    } else if(output[1] == 3){
      list_samples = allclust[allclust$PAM_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,4)]
    } else if(output[1] == 4){
      list_samples = allclust[allclust$Spec_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,5)]
    } else if(output[1] == 5){
      list_samples = allclust[allclust$cross_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,6)]
    } else if(output[1] == 6){
      list_samples = allclust[allclust$CVHC_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,7)]
    } else if(output[1] == 7){
      list_samples = allclust[allclust$CVKM_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,8)]
    } else if(output[1] == 8){
      list_samples = allclust[allclust$CVPAM_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,9)]
    } else if(output[1] == 9){
      list_samples = allclust[allclust$CVSpec_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,10)]
    } else if(output[1] == 10){
      list_samples = allclust[allclust$CVcross_clusters %in% c(1,3),]
      list_samples = list_samples[,c(1,11)]
    }
    
    if(length(list_samples) > 0){
      df_tolm = data.frame("Basename" = list_samples[,1],
                           "Outcome" = list_samples[,2])
      df_tolm$Outcome = as.numeric(df_tolm$Outcome)
      print(df_tolm)
      
      betas_tolm = as.data.frame(t(betas_AD[,colnames(betas_AD) %in% df_tolm$Basename]))
      betas_tolm = cbind(df_tolm$Outcome, betas_tolm)
      colnames(betas_tolm)[1] = "Outcome"
      
      cores=detectCores()
      cl <- makeCluster(cores[1]-2) #not to overload your computer
      registerDoParallel(cl)
      clusterExport(cl, "cohensD")
      
      
      df_res_lm <- foreach(ncpg=c(2:ncol(betas_tolm)), .combine=rbind) %dopar% {
        input_betas = betas_tolm[,c(1,ncpg)]
        
        cpg_name = colnames(input_betas)[2]
        colnames(input_betas)[2] = "predictor"
        
        res_lm = lm(formula = Outcome ~ predictor, data = input_betas)
        sum_lm = summary(res_lm)
        cohd = cohensD(formula = predictor ~ Outcome, data = input_betas)
        
        new_line = c(output[1],output[2],cpg_name, sum_lm$coefficients[2,4], sum_lm$r.squared, res_lm$effects[2], cohd)
        
        new_line #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
      }
      stopCluster(cl)
      colnames(df_res_lm) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
      df_res_lm_final = rbind(df_res_lm_final, df_res_lm)
    }
    
  }
  
  df_res_lm_final = as.data.frame(df_res_lm_final)
  df_res_lm_final = df_res_lm_final[-1,]
  colnames(df_res_lm_final) = c("model","cluster","cpg","pvalue","rsquared","effect","cohend")
  df_res_lm_final[df_res_lm_final$model == "1",]$model = "HClust"
  df_res_lm_final[df_res_lm_final$model == "2",]$model = "KMeans"
  df_res_lm_final[df_res_lm_final$model == "3",]$model = "PAM"
  df_res_lm_final[df_res_lm_final$model == "4",]$model = "Spectrum"
  df_res_lm_final[df_res_lm_final$model == "5",]$model = "Crossmethod"
  df_res_lm_final[df_res_lm_final$model == "6",]$model = "CV_HClust"
  df_res_lm_final[df_res_lm_final$model == "7",]$model = "CV_KMeans"
  df_res_lm_final[df_res_lm_final$model == "8",]$model = "CV_PAM"
  df_res_lm_final[df_res_lm_final$model == "9",]$model = "CV_Spectrum"
  df_res_lm_final[df_res_lm_final$model == "10",]$model = "CV_Crossmethod"
  
  return(df_res_lm_final)
}

